"""
Shared pipeline: Welch FFT → CV-weighted smoothing → correction → IIR fit → C DSP.

This module is the single source of truth for the EQ design and audio processing
pipeline.  All CLI entry points (eqgen, audition, wire, export) delegate here.

Exports:
  run_pipeline()   — WAV inputs → (freqs, gains_db, sample_rate)
  design_eq()      — EQ curve → quantized biquad coefficients
  process_track()  — audio file → C enhancer → WAV
  curve_to_json()  — serialise (freqs, gains_db) for JSON output
"""

import json
import os
import random
import struct
import subprocess
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from eqgen.eq_fit import BiquadCoeffs, cascade_response_db, fit_eq_curve
from eqgen.quantize import BiquadQ28, q28_to_float, quantize_biquads_q28
from eqgen import enhancer_ffi

ROOT = Path(__file__).resolve().parent.parent

# ═══════════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════════

BOTTOM_F = 20.0       # Hz — lowest frequency analyzed
TOP_F = 14000.0       # Hz — highest frequency analyzed
WELCH_FFT_SIZE = 16384
WELCH_OVERLAP = 0.5

# Evaluation grid for the output curve — matches the IIR fitter's internal
# log-spaced grid so no re-interpolation is needed downstream.
N_EVAL = 256


# ═══════════════════════════════════════════════════════════════════════════════
# dB conversion helpers
# ═══════════════════════════════════════════════════════════════════════════════

def db_to_ratio(db: float) -> float:
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio: float) -> float:
    return 20.0 * np.log10(np.maximum(ratio, 1e-20))


# ═══════════════════════════════════════════════════════════════════════════════
# WAV I/O
# ═══════════════════════════════════════════════════════════════════════════════

def read_wav(path: str) -> Tuple[np.ndarray, float]:
    """Read a mono or stereo WAV, returning mono float64 samples in [-1, 1]."""
    from scipy.io import wavfile

    rate, data = wavfile.read(path)

    if data.dtype == np.int16:
        data = data.astype(np.float64) / 32768.0
    elif data.dtype == np.int32:
        data = data.astype(np.float64) / 2147483648.0
    elif data.dtype == np.uint8:
        data = (data.astype(np.float64) - 128.0) / 128.0
    elif data.dtype == np.float32:
        data = data.astype(np.float64)
    else:
        raise ValueError(f"Unsupported WAV dtype: {data.dtype}")

    if data.ndim == 2:
        data = data.mean(axis=1)
    elif data.ndim != 1:
        raise ValueError(f"Unexpected WAV shape: {data.shape}")

    return np.ascontiguousarray(data, dtype=np.float64), float(rate)


# ═══════════════════════════════════════════════════════════════════════════════
# Welch's method FFT
# ═══════════════════════════════════════════════════════════════════════════════

def welch_stats(samples: np.ndarray, sample_rate: float) -> List[dict]:
    """Run Welch's method and return per-bin mean magnitude, power, and CV."""
    fft_size = WELCH_FFT_SIZE
    noverlap = int(fft_size * WELCH_OVERLAP)

    step = fft_size - noverlap
    win = np.hanning(fft_size)
    bin_count = fft_size // 2 + 1

    raw_sum = np.zeros(bin_count)
    raw_sum_sq = np.zeros(bin_count)
    norm_sum = np.zeros(bin_count)
    norm_sum_sq = np.zeros(bin_count)
    window_count = 0

    for start in range(0, len(samples) - fft_size + 1, step):
        segment = samples[start:start + fft_size] * win
        spectrum = np.fft.rfft(segment)
        mags_window = np.abs(spectrum) / fft_size
        raw_sum += mags_window
        raw_sum_sq += mags_window * mags_window

        mean_mag = np.mean(mags_window)
        if mean_mag > 0:
            nm = mags_window / mean_mag
            norm_sum += nm
            norm_sum_sq += nm * nm
        else:
            norm_sum += 1.0
            norm_sum_sq += 1.0
        window_count += 1

    n = float(window_count)
    out_freqs = np.fft.rfftfreq(fft_size, 1.0 / sample_rate)

    results = []
    for i in range(bin_count):
        f = out_freqs[i]
        if f < BOTTOM_F or f > TOP_F:
            continue
        mean_val = max(raw_sum[i] / n, 0.0)
        mean_sq = max(raw_sum_sq[i] / n, 0.0)
        norm_mean = norm_sum[i] / n
        norm_var = (norm_sum_sq[i] / n) - (norm_mean * norm_mean)
        cv = (np.sqrt(norm_var) / norm_mean) if norm_mean > 0 else float('inf')
        results.append({
            "freq": f,
            "count": window_count,
            "mean": mean_val,
            "mean_sq": mean_sq,
            "cv": cv,
        })
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# CV-weighted kernel smoother
# ═══════════════════════════════════════════════════════════════════════════════

def _find_viable_range(freqs: np.ndarray, cv: np.ndarray,
                       cv_threshold: float = 2.0,
                       margin_oct: float = 0.25) -> Tuple[float, float]:
    """Find the frequency range where CV stays below *cv_threshold*.

    Works outward from the middle: the last frequency where a running
    maximum of CV (over ~⅓ octave) exceeds the threshold defines the
    boundary.  A *margin_oct* pad is subtracted to stay safely inside.
    If no shelf is found, returns the full range.
    """
    log_f = np.log2(freqs)
    n = len(freqs)

    # Running max over ~⅓ octave to catch sustained CV shelves, not single-bin spikes
    bin_width = (log_f[-1] - log_f[0]) / (n - 1) if n > 1 else 0.0
    radius = max(1, int(np.ceil((1.0 / 3.0) / bin_width))) if bin_width > 0 else 1
    kernel = np.ones(2 * radius + 1)
    cv_smooth = np.convolve(np.maximum(cv, 0), kernel, mode='same') / len(kernel)

    mid = n // 2
    margin_bins = int(np.ceil(margin_oct / bin_width)) if bin_width > 0 else 0

    # Low boundary: walk left from middle, find first bin exceeding threshold
    lo = 0
    for i in range(mid, 0, -1):
        if cv_smooth[i] > cv_threshold:
            lo = min(i + 1 + margin_bins, n - 1)
            break

    # High boundary: walk right from middle
    hi = n - 1
    for i in range(mid, n):
        if cv_smooth[i] > cv_threshold:
            hi = max(i - 1 - margin_bins, 0)
            break

    return freqs[lo], freqs[hi]


def _smooth_kernel(bin_freqs: np.ndarray, bin_values: np.ndarray,
                   bin_cv: np.ndarray, eval_freqs: np.ndarray,
                   bandwidth_oct = 0.5) -> np.ndarray:
    """Evaluate a CV-weighted kernel smoother on a log-spaced grid.

    Each output point is a weighted average of input bins within
    ±bandwidth_oct octaves.  *bandwidth_oct* may be a scalar or a
    per-evaluation-point array.  Distance weight (tricube) decays to
    zero at the bandwidth edge.  Confidence weight is 1/CV, so low-CV
    bins dominate their local region while high-CV bins are averaged out.

    Unlike a cubic spline, this cannot oscillate below zero or collapse
    to a flat line — the output is always a convex combination of input
    values.
    """
    bandwidth_oct = np.asarray(bandwidth_oct, dtype=float)
    scalar_bw = bandwidth_oct.ndim == 0

    n_out = len(eval_freqs)
    log_freqs = np.log2(bin_freqs)
    log_eval = np.log2(eval_freqs)
    conf = 1.0 / np.clip(bin_cv, 0.01, None)

    result = np.empty(n_out)
    # Process one point at a time since bandwidth may vary per point
    for i in range(n_out):
        d = np.abs(log_eval[i] - log_freqs)
        bw = float(bandwidth_oct) if scalar_bw else float(bandwidth_oct[i])
        if bw <= 0.0:
            # Zero bandwidth — nearest-bin interpolation
            j = np.argmin(d)
            result[i] = bin_values[j]
            continue
        d_norm = d / bw
        w_dist = np.where(d_norm < 1.0, (1.0 - d_norm ** 3) ** 3, 0.0)
        w = w_dist * conf
        w_sum = w.sum()
        if w_sum > 0:
            result[i] = np.average(bin_values, weights=w)
        else:
            result[i] = np.average(bin_values, weights=conf)
    return np.maximum(result, 0.0)


# ═══════════════════════════════════════════════════════════════════════════════
# Core pipeline: WAV inputs → EQ curve
# ═══════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    measurement_paths: List[str],
    target_path: str,
    noise_path: Optional[str] = None,
    bass_enhancer_cutoff: Optional[float] = None,
    h2: float = 1.0,
    h3: float = 1.0,
    high_rolloffs: Optional[List[Tuple[float, float]]] = None,
    low_rolloffs: Optional[List[Tuple[float, float]]] = None,
    sample_rate_override: Optional[float] = None,
    smooth_exponent: float = 1.0,
    detailed: bool = False,
):
    """Run the full EQ correction pipeline.

    Without detailed: returns (freqs_hz, gains_db, sample_rate).
    With detailed=True: returns a dict with all intermediate data for
    visualization.

    Steps:
      1. Read & validate WAV files
      2. Welch FFT on each input
      3. Pool per-bin stats across measurement files
      4. Merge noise-floor CV + inflate CV in noise-dominated bins
      5. Build CV-weighted smoothing splines (measurement, target, noise)
      6. Evaluate on uniform log-spaced grid + spectral subtraction
      7. Apply rolloffs to target
      8. Correction = target / measurement
      9. Bass enhancer preprocessing (if cutoff provided)
    """
    if high_rolloffs is None:
        high_rolloffs = []
    if low_rolloffs is None:
        low_rolloffs = []

    # 1. Read & validate
    wavs = []
    for path in measurement_paths:
        samples, rate = read_wav(path)
        wavs.append((samples, rate, path))

    target_samples, target_rate = read_wav(target_path)
    sample_rate = sample_rate_override or wavs[0][1]
    for s, r, p in wavs[1:]:
        assert r == sample_rate, f"Sample rate mismatch in {p}: {r} vs {sample_rate}"
    assert target_rate == sample_rate, \
        f"Target sample rate mismatch: {target_rate} vs {sample_rate}"

    # 2. Welch FFT
    all_stats = [welch_stats(s, sample_rate) for s, _, _ in wavs]
    bin_count = len(all_stats[0])
    for s in all_stats[1:]:
        assert len(s) == bin_count, "Frequency ranges differ"

    # 3. Pool stats across measurement files
    pooled = []
    for i in range(bin_count):
        total_n = sum(s[i]["count"] for s in all_stats)
        total_sum = sum(s[i]["mean"] * s[i]["count"] for s in all_stats)
        total_sum_sq = sum(s[i]["mean_sq"] * s[i]["count"] for s in all_stats)
        pooled_cv = sum(s[i]["cv"] * s[i]["count"] for s in all_stats) / total_n
        pooled.append({
            "freq": all_stats[0][i]["freq"],
            "count": total_n,
            "mean": total_sum / total_n,
            "mean_sq": total_sum_sq / total_n,
            "cv": pooled_cv,
        })

    # 4. Noise merge: max-CV + noise-floor inflation
    noise_stats = None
    if noise_path:
        noise_samples, noise_rate = read_wav(noise_path)
        assert noise_rate == sample_rate
        noise_stats = welch_stats(noise_samples, sample_rate)

        # 4a. max() merge — high-CV intermittent noise dominates
        for i in range(min(len(pooled), len(noise_stats))):
            pooled[i]["cv"] = max(pooled[i]["cv"], noise_stats[i]["cv"])

        # 4b. Noise-floor inflation — stationary noise with poor SNR
        for i in range(len(pooled)):
            P_meas = max(pooled[i]["mean_sq"], 1e-20)
            P_noise = max(noise_stats[i]["mean_sq"], 1e-20)
            cv_n = noise_stats[i]["cv"]
            # Inflate noise power by 1σ — accounts for peak, not mean, noise
            effective_noise = P_noise * (1.0 + cv_n)
            noise_ratio = min(0.99, effective_noise / P_meas)
            pooled[i]["cv"] = pooled[i]["cv"] / np.sqrt(max(0.01, 1.0 - noise_ratio))

    # 5. Extract arrays for kernel smoothing
    meas_freqs = np.array([s["freq"] for s in pooled])
    meas_raw = np.array([s["mean"] for s in pooled])
    meas_cv = np.array([s["cv"] for s in pooled])

    target_stats = welch_stats(target_samples, sample_rate)
    targ_freqs = np.array([s["freq"] for s in target_stats])
    targ_raw = np.array([s["mean"] for s in target_stats])
    targ_cv = np.array([s["cv"] for s in target_stats])

    noise_freqs = None
    noise_raw = None
    noise_cv = None
    if noise_stats is not None:
        noise_freqs = np.array([s["freq"] for s in noise_stats])
        noise_raw = np.array([s["mean"] for s in noise_stats])
        noise_cv = np.array([s["cv"] for s in noise_stats])

    # 6. CV-weighted kernel smoothing on uniform log-spaced grid.
    # Bandwidth scales with local CV: at the Rayleigh floor (0.52) it
    # is tiny (essentially no smoothing); at CV=2+ it widens significantly.
    MIN_CV = 0.52   # Rayleigh CV — the noise floor for stationary signals
    BASE_BW = 0.08  # octaves at MIN_CV (~3 FFT bins at 400 Hz)

    eval_freqs = np.logspace(np.log10(BOTTOM_F), np.log10(TOP_F), N_EVAL)

    # Estimate local CV at each eval point to scale bandwidth
    cv_at_eval = _smooth_kernel(meas_freqs, meas_cv, meas_cv, eval_freqs,
                                bandwidth_oct=0.3)
    bw_meas = np.maximum(BASE_BW * (cv_at_eval / MIN_CV) ** smooth_exponent, 0.01)
    meas_vals = _smooth_kernel(meas_freqs, meas_raw, meas_cv, eval_freqs,
                               bandwidth_oct=bw_meas)

    # Spectral subtraction
    if noise_freqs is not None:
        noise_vals = _smooth_kernel(noise_freqs, noise_raw, noise_cv, eval_freqs)
        meas_vals = np.maximum(meas_vals - noise_vals, meas_vals * 0.01)

    # Target uses its own CV-scaled bandwidth
    cv_targ = _smooth_kernel(targ_freqs, targ_cv, targ_cv, eval_freqs,
                             bandwidth_oct=0.3)
    bw_targ = np.maximum(BASE_BW * (cv_targ / MIN_CV) ** smooth_exponent, 0.01)
    target_vals = _smooth_kernel(targ_freqs, targ_raw, targ_cv, eval_freqs,
                                 bandwidth_oct=bw_targ)

    # 7. Rolloffs
    for freq, db in high_rolloffs:
        atten = db_to_ratio(db)
        mask = eval_freqs > freq
        octaves = np.log2(np.maximum(eval_freqs[mask] / freq, 1.0))
        target_vals[mask] *= atten ** octaves
    for freq, db in low_rolloffs:
        atten = db_to_ratio(db)
        mask = eval_freqs < freq
        octaves = np.log2(np.maximum(freq / eval_freqs[mask], 1.0))
        target_vals[mask] *= atten ** octaves

    # 8. Correction
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = np.where(meas_vals > 1e-20, target_vals / meas_vals, 1e6)

    # 9. Bass enhancer preprocessing
    if bass_enhancer_cutoff is not None:
        fc = bass_enhancer_cutoff
        from eqgen.model import model_gain_needed

        def meas_at(f: float) -> float:
            """Measurement magnitude (kernel-smoothed) at frequency f."""
            if f <= 0 or f > TOP_F:
                return 0.0
            return float(_smooth_kernel(meas_freqs, meas_raw, meas_cv,
                                         np.array([f]), bandwidth_oct=0.08).item())

        # Compute flat correction value at fc/2 before modifying corr
        idx_fc2 = int(np.searchsorted(eval_freqs, fc / 2.0))
        m_fc2 = meas_at(fc / 2.0)
        G_flat = model_gain_needed(fc / 2.0, corr[idx_fc2] * m_fc2, fc, h2, h3, meas_at)
        flat_val = min(G_flat, corr[idx_fc2])

        for i, f in enumerate(eval_freqs):
            G = model_gain_needed(f, corr[i] * meas_at(f), fc, h2, h3, meas_at)
            if G > 1e-12:
                corr[i] = min(G, corr[i])

        # Below fc/2 the enhancer handles bass perception via harmonics.
        # Hold correction flat at the model gain computed at fc/2.
        corr[:idx_fc2 + 1] = flat_val

    # 10. Clamp correction flat outside the CV-defined viable range.
    # CV shelves (sudden sustained increases) indicate frequencies where
    # the measurement is too unreliable to correct — hold the correction
    # constant beyond those boundaries.
    lo_f, hi_f = _find_viable_range(meas_freqs, meas_cv)
    lo_idx = int(np.searchsorted(eval_freqs, lo_f))
    hi_idx = int(np.searchsorted(eval_freqs, hi_f))
    if lo_idx > 0:
        corr[:lo_idx] = corr[lo_idx]
    if hi_idx < len(corr) - 1:
        corr[hi_idx + 1:] = corr[hi_idx]

    freqs = eval_freqs
    gains_db = ratio_to_db(corr)

    # Normalize: re-center correction so midrange mean = 0 dB.
    # Raw FFT magnitudes are uncalibrated — measurement and target WAVs
    # may have different recording levels, creating a spurious DC offset
    # in the correction curve.  Subtracting the midrange mean gives the
    # IIR fitter a balanced mix of positive and negative peaks around
    # 0 dB instead of an all-positive curve with a large offset.
    mid_mask = (freqs >= 500.0) & (freqs <= 2000.0)
    if np.any(mid_mask):
        mid_mean = float(np.mean(gains_db[mid_mask]))
        gains_db = gains_db - mid_mean

    if detailed:
        # Raw measurement (Welch bins)
        raw_resp = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                    for s in pooled]
        # Raw target (Welch bins)
        raw_target = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                      for s in target_stats]
        # CV per bin (after merge + inflation)
        cv_data = [{"freq": p["freq"], "cv": p["cv"]} for p in pooled]
        # Noise floor CV
        noise_data = []
        if noise_stats is not None:
            noise_data = [{"freq": s["freq"], "cv": s["cv"]} for s in noise_stats]
        # Subtracted measurement (eval grid)
        meas_sub = [{"freq": float(f), "db": ratio_to_db(v)}
                    for f, v in zip(eval_freqs, meas_vals)]
        # Target after rolloffs (eval grid)
        targ_sub = [{"freq": float(f), "db": ratio_to_db(v)}
                    for f, v in zip(eval_freqs, target_vals)]

        return {
            "freqs": freqs.tolist(),
            "gains_db": gains_db.tolist(),
            "sample_rate": sample_rate,
            "raw_measurement": raw_resp,
            "raw_target": raw_target,
            "noise_cv": cv_data,
            "noise_floor": noise_data,
            "meas_subtracted": meas_sub,
            "target_resampled": targ_sub,
        }

    return freqs, gains_db, sample_rate


# ═══════════════════════════════════════════════════════════════════════════════
# Serialisation
# ═══════════════════════════════════════════════════════════════════════════════

def curve_to_json(freqs: np.ndarray, gains_db: np.ndarray) -> list:
    """Convert (freqs, gains_db) to JSON-serialisable list of {freq, gain_db}."""
    return [
        {"freq": round(float(f), 2), "gain_db": round(float(g), 3)}
        for f, g in zip(freqs, gains_db)
    ]


# ═══════════════════════════════════════════════════════════════════════════════
# IIR design: EQ curve → quantized biquads
# ═══════════════════════════════════════════════════════════════════════════════

def design_eq(
    freqs: np.ndarray,
    target_db: np.ndarray,
    fs: float,
    max_bands: int = 40,
) -> Tuple[List[int], List[dict], np.ndarray, np.ndarray, np.ndarray]:
    """Design a quantized IIR EQ to match the target curve.

    Returns (coeffs_flat_q28, bands, freqs, target_db, fitted_db).
    coeffs_flat_q28 is a flat list of int32 Q4.28 coefficients:
    [b0, b1, b2, a1, a2, b0, b1, b2, a1, a2, ...].
    """
    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=freqs[0])

    bq_q28 = quantize_biquads_q28(fit.biquads)
    coeffs = []
    for bq in bq_q28:
        coeffs.extend([bq.b0, bq.b1, bq.b2, bq.a1, bq.a2])

    q28_floats = [BiquadCoeffs(b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                               b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                               a2=q28_to_float(b.a2)) for b in bq_q28]
    fitted_db = cascade_response_db(q28_floats, freqs, fs)

    return coeffs, fit.bands, freqs, target_db, fitted_db


# ═══════════════════════════════════════════════════════════════════════════════
# Default EQ coefficients (flat 3-HP cascade)
# ═══════════════════════════════════════════════════════════════════════════════

def build_default_eq_coeffs(fs: float = 44100.0) -> List[int]:
    """Build a flat 3-HP-cascade default EQ as Q4.28 int list."""
    coeffs = []
    SQRT2 = np.sqrt(2.0)
    for fc in [25, 35, 45]:
        omega = np.tan(np.pi * fc / fs)
        c = 1.0 + SQRT2 * omega + omega * omega
        b0 = int(np.round((1.0 / c) * (1 << 28)))
        b1 = int(np.round((-2.0 / c) * (1 << 28)))
        b2 = int(np.round((1.0 / c) * (1 << 28)))
        a1 = int(np.round((2.0 * (omega * omega - 1.0) / c) * (1 << 28)))
        a2 = int(np.round(((1.0 - SQRT2 * omega + omega * omega) / c) * (1 << 28)))
        coeffs.extend([b0, b1, b2, a1, a2])
    return coeffs


# ═══════════════════════════════════════════════════════════════════════════════
# Audio processing: file → C enhancer → WAV
# ═══════════════════════════════════════════════════════════════════════════════

def process_track(
    input_path: str,
    output_path: str,
    coeffs_q28: List[int],
    n_biquads: int,
    cutoff_hz: float = 60.0,
    h2: float = 0.5,
    h3: float = 1.0,
    start_sec: float = 30.0,
    duration_sec: float = 40.0,
    release_secs: float = 0.2,
) -> bool:
    """Decode audio via ffmpeg, process through C enhancer, write WAV.

    Returns True on success.
    """
    cmd = ["ffmpeg", "-y", "-v", "error",
           "-ss", str(start_sec), "-t", str(duration_sec),
           "-i", input_path,
           "-f", "s16le", "-acodec", "pcm_s16le",
           "-ar", "44100", "-ac", "2", "pipe:1"]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0 or len(result.stdout) < 4:
        print(f"  ffmpeg failed")
        return False

    pcm = result.stdout
    n_frames = len(pcm) // 4
    print(f"  44100Hz stereo, {n_frames} frames ({n_frames/44100:.1f}s)")

    enh = enhancer_ffi.create_enhancer(
        cutoff_hz=cutoff_hz, h2_amp=h2, h3_amp=h3,
        release_secs=release_secs,
        fs=44100.0, coeffs_q28=coeffs_q28)

    out_data = bytearray(len(pcm))
    peak_in = 0
    peak_out = 0
    for i in range(0, len(pcm), 4):
        l = struct.unpack_from('<h', pcm, i)[0]
        r = struct.unpack_from('<h', pcm, i+2)[0]
        l_out, r_out = enhancer_ffi.process_stereo_frame(enh, l, r)
        struct.pack_into('<hh', out_data, i, l_out, r_out)
        if abs(l_out) > peak_out: peak_out = abs(l_out)
        if abs(r_out) > peak_out: peak_out = abs(r_out)
        if abs(l) > peak_in: peak_in = abs(l)
        if abs(r) > peak_in: peak_in = abs(r)

    enhancer_ffi.destroy_enhancer(enh)

    datasize = len(out_data)
    wav = bytearray(44 + datasize)
    wav[0:4] = b'RIFF'
    struct.pack_into('<I', wav, 4, 36 + datasize)
    wav[8:16] = b'WAVEfmt '
    struct.pack_into('<I', wav, 16, 16)
    struct.pack_into('<H', wav, 20, 1)
    struct.pack_into('<H', wav, 22, 2)
    struct.pack_into('<I', wav, 24, 44100)
    struct.pack_into('<I', wav, 28, 44100 * 4)
    struct.pack_into('<H', wav, 32, 4)
    struct.pack_into('<H', wav, 34, 16)
    wav[36:40] = b'data'
    struct.pack_into('<I', wav, 40, datasize)
    wav[44:] = out_data

    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    with open(output_path, 'wb') as f:
        f.write(wav)

    gain_db = 20.0 * np.log10(max(peak_out, 1) / max(peak_in, 1))
    print(f"  Peak: in={peak_in} out={peak_out} ({gain_db:+.1f} dB) → {output_path}")
    if peak_out >= 32767:
        print("  ⚠️  CLIPPING")
    return True
