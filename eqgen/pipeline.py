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
from eqgen.presets import MAX_IIR_BANDS
from eqgen import enhancer_ffi
from eqgen.dsp import first_order_lp_mag
from scipy.interpolate import interp1d

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
N_EVAL = 512


# ═══════════════════════════════════════════════════════════════════════════════
# dB conversion helpers
# ═══════════════════════════════════════════════════════════════════════════════

def db_to_ratio(db: float) -> float:
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio: float) -> float:
    return 20.0 * np.log10(np.maximum(ratio, 1e-20))


def pre_gain_from_max_gain(max_gain_db: float) -> float:
    """Attenuation factor that keeps biquad internal states below overflow.

    Input is attenuated by this factor before the biquad cascade;
    biquads fit the full (unshifted) correction curve and boost
    back up to unity at peak-gain frequencies.
    """
    return 1.0 / max(1.0, db_to_ratio(max_gain_db))


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
    sample_rate_override: Optional[float] = None,
    smooth_exponent: float = 1.0,
    detailed: bool = False,
    n_eval: int = 512,
    house_curve: Optional[list] = None,
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
      5. CV-weighted kernel smoothing on measurement + spectral subtraction
      6. Target = linear fit of raw target (dB vs log-freq),
         shifted to measurement midrange level
      7. Correction = target / measurement
      8. Bass enhancer preprocessing (if cutoff provided)
    """

    # 1. Read & validate
    wavs = []
    for path in measurement_paths:
        samples, rate = read_wav(path)
        wavs.append((samples, rate, path))

    sample_rate = sample_rate_override or wavs[0][1]
    for s, r, p in wavs[1:]:
        assert r == sample_rate, f"Sample rate mismatch in {p}: {r} vs {sample_rate}"

    if target_path:
        target_samples, target_rate = read_wav(target_path)
        assert target_rate == sample_rate, \
            f"Target sample rate mismatch: {target_rate} vs {sample_rate}"
    else:
        target_samples = None

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

    # 5. CV-weighted kernel smoothing on measurement only.
    # Bandwidth scales with local CV: at the Rayleigh floor (0.52) it
    meas_freqs = np.array([s["freq"] for s in pooled])
    meas_raw = np.array([s["mean"] for s in pooled])
    meas_cv = np.array([s["cv"] for s in pooled])

    if target_path:
        target_stats = welch_stats(target_samples, sample_rate)
        targ_freqs = np.array([s["freq"] for s in target_stats])
        targ_raw_db = ratio_to_db(np.array([s["mean"] for s in target_stats]))
    else:
        # No target WAV — use a flat 0 dB line spanning the measurement range
        target_stats = None
        targ_freqs = np.array([BOTTOM_F, TOP_F])
        targ_raw_db = np.array([0.0, 0.0])

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

    eval_freqs = np.logspace(np.log10(BOTTOM_F), np.log10(TOP_F), n_eval)

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

    # 6. Target from WAV (linear fit in dB vs log-freq) or flat 0 dB line.
    #    House curve is an additive adjustment applied on top.
    eval_log10 = np.log10(eval_freqs)
    mid = (eval_freqs >= 500) & (eval_freqs <= 2000)

    log10_f_targ = np.log10(targ_freqs)
    slope, intercept = np.polyfit(log10_f_targ, targ_raw_db, 1)
    targ_db_eval = slope * eval_log10 + intercept

    if house_curve is not None:
        curve_freqs = np.array([p[0] for p in house_curve])
        curve_db = np.array([p[1] for p in house_curve])
        targ_db_eval = targ_db_eval + np.interp(eval_freqs, curve_freqs, curve_db)

    meas_mid_db = float(np.mean(ratio_to_db(meas_vals[mid]))) if mid.any() else 0.0
    targ_mid_db = float(np.mean(targ_db_eval[mid])) if mid.any() else 0.0
    targ_db_eval += meas_mid_db - targ_mid_db
    target_vals = 10.0 ** (targ_db_eval / 20.0)

    # 7. Correction
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = np.where(meas_vals > 1e-20, target_vals / meas_vals, 1e6)

    # 8. Bass enhancer preprocessing (skip if cutoff ≤ 0 — bypass)
    if bass_enhancer_cutoff is not None and bass_enhancer_cutoff > 0:
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

    # 9. Clamp correction flat outside the CV-defined viable range.
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

    # Pre-gain: the max positive gain of the *normalized* correction curve.
    # By applying this as a uniform gain before the EQ biquads, the
    # fitted EQ curve only needs cuts (negative gains), avoiding
    # internal clipping from large boosts in the biquad cascade.
    # Must be computed AFTER normalization so recording-level offsets
    # don't inflate the pre-gain.
    max_gain_db = max(0.0, float(np.max(gains_db)))

    if detailed:
        # Raw measurement (Welch bins)
        raw_resp = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                    for s in pooled]
        # Raw target (Welch bins from WAV, or empty if no target WAV)
        if target_stats is not None:
            raw_target = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                          for s in target_stats]
        else:
            raw_target = []
        # CV per bin (after merge + inflation)
        cv_data = [{"freq": p["freq"], "cv": p["cv"]} for p in pooled]
        # Noise floor CV
        noise_data = []
        if noise_stats is not None:
            noise_data = [{"freq": s["freq"], "cv": s["cv"]} for s in noise_stats]
        # Subtracted measurement (eval grid)
        meas_sub = [{"freq": float(f), "db": ratio_to_db(v)}
                    for f, v in zip(eval_freqs, meas_vals)]
        # Target = linear fit (eval grid, dB)
        targ_sub = [{"freq": float(f), "db": float(d)}
                    for f, d in zip(eval_freqs, targ_db_eval)]

        return {
            "freqs": freqs.tolist(),
            "gains_db": gains_db.tolist(),
            "sample_rate": sample_rate,
            "max_gain_db": max_gain_db,
            "raw_measurement": raw_resp,
            "raw_target": raw_target,
            "noise_cv": cv_data,
            "noise_floor": noise_data,
            "meas_subtracted": meas_sub,
            "target_resampled": targ_sub,
        }

    return freqs, gains_db, sample_rate, max_gain_db


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
# FFT + IIR hybrid: FFT handles broad shape, IIR fits narrow residual
# ═══════════════════════════════════════════════════════════════════════════════

FFT_N = 256  # window size for overlap-add FFT EQ
FFT_RESIDUAL_THRESHOLD_DB = 1.0  # IIR only corrects where |residual| > this


def compute_fft_residual(
    freqs: np.ndarray,
    target_db: np.ndarray,
    fs: float,
    fft_n: int = FFT_N,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute the FFT approximation of a target curve and the residual.

    In the FFT + IIR hybrid architecture, the FFT overlap-add EQ applies
    per-bin scalar gains at ~187.5 Hz resolution (256-point @ 48 kHz).
    The IIR biquads then surgically correct the narrow features the FFT
    can't resolve — the *residual*.

    The residual is self-normalizing: it oscillates around 0 dB because
    the FFT handles the broad shape.  No pre-gain offset is needed —
    biquads handle both boosts and cuts, and a uniform post-gain after
    the cascade compensates for any overall level shift.

    This function:
      1. Samples the target curve at FFT bin center frequencies
      2. Cubic-spline interpolates back to the eval grid as the FFT
         approximation (Hann window blends ~3 bins smoothly)
      3. Subtracts to get the residual
      4. Thresholds residual at ±1 dB

    Returns:
        fft_bin_freqs:  bin center frequencies (Hz)
        fft_gains_linear: per-bin scalar gains (linear)
        fft_approx_db: FFT approximation on the eval grid (dB)
        residual_db:   target - FFT approx, thresholded (dB)
    """
    # 1. Bin center frequencies
    fft_bin_freqs = np.arange(fft_n // 2 + 1, dtype=float) * fs / fft_n

    # 2. Target curve sampled at bin centers, converted to linear gains.
    #    The DC bin (0 Hz) is pinned to 1.0 — there is no musical content at DC,
    #    and np.interp extrapolates an arbitrary value from the first eval point
    #    (~20 Hz).  A large DC gain would distort all bass frequencies through
    #    the Hann window's spectral leakage (main lobe ~750 Hz wide).
    target_at_bins = np.interp(fft_bin_freqs, freqs, target_db)
    target_at_bins[0] = 0.0  # 0 dB → linear gain = 1.0
    gains_linear = 10.0 ** (target_at_bins / 20.0)

    # Ramp DC → first few bins to prevent time-domain aliasing in the
    # FFT overlap-add.  A step from 0 dB (DC) to the first bin's gain
    # creates an impulse response far longer than the 256-point window,
    # causing wrap-around distortion on bass tones.
    # The IIR picks up whatever the FFT ramp leaves behind.
    ramp_end = 4  # ramp bins 0..3 linearly in dB toward bin 4
    if ramp_end < len(gains_linear):
        ramp_target_db = float(20.0 * np.log10(max(gains_linear[ramp_end], 1e-12)))
        for i in range(1, ramp_end):
            t = i / ramp_end
            gains_linear[i] = 10.0 ** (t * ramp_target_db / 20.0)

    fft_gains_db = 20.0 * np.log10(np.maximum(gains_linear, 1e-12))

    # 3. Cubic spline interpolates the piecewise-constant bin gains back to
    #    the eval grid, approximating the Hann window's ~3-bin blending.
    fft_interp = interp1d(fft_bin_freqs, fft_gains_db, kind='cubic',
                          bounds_error=False, fill_value='extrapolate')
    fft_approx_db = fft_interp(freqs)

    # 4. Residual = what the FFT missed (narrow resonances, sharp features)
    residual_db = target_db - fft_approx_db

    # 5. Threshold: IIR only corrects where residual exceeds ±1 dB.
    #    Below that, the FFT is close enough — don't waste biquads on noise.
    residual_db[np.abs(residual_db) < FFT_RESIDUAL_THRESHOLD_DB] = 0.0

    return fft_bin_freqs, gains_linear, fft_approx_db, residual_db


# ═══════════════════════════════════════════════════════════════════════════════
# IIR design: EQ curve → quantized biquads
# ═══════════════════════════════════════════════════════════════════════════════

def design_eq(
    freqs: np.ndarray,
    target_db: np.ndarray,
    fs: float,
    max_bands: int = MAX_IIR_BANDS,
    fft_n: int = 0,
    min_peaking_freq: float = 0.0,
) -> Tuple[List[float], List[dict], np.ndarray, np.ndarray, np.ndarray]:
    """Design an IIR EQ to match a target curve.

    When fft_n > 0 (FFT hybrid mode), the IIR biquads are fitted to
    the *residual* after subtracting the FFT's broad correction.
    The FFT handles the overall shape at coarse bin resolution;
    the IIR surgically corrects narrow peaks and dips.

    Returns (coeffs_flat, bands, freqs, fit_target_db, fitted_db).
    fit_target_db is the curve actually fitted (full target or residual).
    coeffs_flat is a flat list of float coefficients:
    [b0, b1, b2, a1, a2, b0, b1, b2, a1, a2, ...].
    """
    if fft_n > 0:
        # Residual is self-normalizing — centered around 0 dB by definition.
        # FFT handles the broad shape, IIR surgically corrects narrow features.
        _, _, _, residual_db = compute_fft_residual(freqs, target_db, fs, fft_n)
        shifted_db = residual_db
        fit_target = residual_db
    else:
        shifted_db = target_db
        fit_target = target_db

    fit = fit_eq_curve(freqs, shifted_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=min_peaking_freq if min_peaking_freq > 0 else freqs[0])

    coeffs = []
    for bc in fit.biquads:
        coeffs.extend([bc.b0, bc.b1, bc.b2, bc.a1, bc.a2])

    fitted_db = cascade_response_db(fit.biquads, freqs, fs)

    return coeffs, fit.bands, freqs, fit_target, fitted_db


# ═══════════════════════════════════════════════════════════════════════════════
# Default EQ coefficients (flat 3-HP cascade)
# ═══════════════════════════════════════════════════════════════════════════════

def build_default_eq_coeffs(fs: float = 44100.0) -> List[float]:
    """Build a flat 3-HP-cascade default EQ as float list."""
    coeffs = []
    SQRT2 = np.sqrt(2.0)
    for fc in [25, 35, 45]:
        omega = np.tan(np.pi * fc / fs)
        c = 1.0 + SQRT2 * omega + omega * omega
        b0 = 1.0 / c
        b1 = -2.0 / c
        b2 = 1.0 / c
        a1 = (2.0 * (omega * omega - 1.0)) / c
        a2 = (1.0 - SQRT2 * omega + omega * omega) / c
        coeffs.extend([b0, b1, b2, a1, a2])
    return coeffs


# ═══════════════════════════════════════════════════════════════════════════════
# Smart volume simulation: model the correction curve at different volume levels
# ═══════════════════════════════════════════════════════════════════════════════

# Constants matching eq_coeffs.h / firmware
SV_LOUDNESS_FC = 200.0    # Hz — one-pole shelf corner
SV_SHELF_MAX_DB = 8.0     # dB — max shelf boost at vol→0
SV_BLEED_MAX = 0.25       # linear — max LP bleed at vol→0


def compute_smart_volume_curves(
    freqs: np.ndarray,
    correction_db: np.ndarray,
    h2: float,
    h3: float,
    pre_gain_db_loud: float = 0.0,
) -> dict:
    """Compute effective correction curves at different smart-volume levels.

    The correction_db is already normalized (midrange mean = 0 dB).
    The effective curve shows correction + shelf only — pre-gain is reported
    as absolute text, not applied to the curve.

    Returns dict with keys: vol_levels, curves, shelf_curves, params.
    """
    levels = [
        ("Quiet (vol→0)", 0.0),
        ("25% vol", 0.25),
        ("50% vol", 0.5),
        ("75% vol", 0.75),
        ("Loud (vol→127)", 1.0),
    ]

    max_shelf_linear = 10.0 ** (SV_SHELF_MAX_DB / 20.0)

    curves = []
    shelf_curves = []

    for label, t in levels:
        boost_db = SV_SHELF_MAX_DB * (1.0 - t)

        # Absolute pre-gain: pg_loud / [1/max + t*(1 - 1/max)]
        pg_fraction = 1.0 / max_shelf_linear + t * (1.0 - 1.0 / max_shelf_linear)
        pg_db_abs = pre_gain_db_loud + 20.0 * np.log10(max(pg_fraction, 1e-12))

        # Shelf contribution: one-pole low shelf at 200 Hz
        if boost_db > 0.01:
            gain_linear = 10.0 ** (boost_db / 20.0)
            shelf_gain = np.array([
                1.0 + (gain_linear - 1.0) * first_order_lp_mag(f, SV_LOUDNESS_FC)
                for f in freqs
            ])
            shelf_db = 20.0 * np.log10(np.maximum(shelf_gain, 1e-12))
        else:
            shelf_db = np.zeros_like(freqs)

        # Effective: correction + shelf + absolute pre_gain.
        # This shifts the entire curve so you can see that at t=0,
        # bass stays at the loud pre-gain level (−4.8 dB) because
        # shelf (+8) + pre_gain (−12.8) = −4.8 dB — same as loud.
        effective_db = correction_db + shelf_db + pg_db_abs

        curves.append({
            "label": label,
            "t": round(t, 2),
            "pre_gain_db": round(pg_db_abs, 1),
            "shelf_max_db": round(boost_db, 1),
            "freqs": freqs.tolist(),
            "effective_db": effective_db.tolist(),
        })
        shelf_curves.append({
            "label": label,
            "t": round(t, 2),
            "freqs": freqs.tolist(),
            "shelf_db": shelf_db.tolist(),
        })

    return {
        "vol_levels": [(lbl, t) for lbl, t in levels],
        "curves": curves,
        "shelf_curves": shelf_curves,
        "params": {
            "fc": float(SV_LOUDNESS_FC),
            "shelf_max_db": float(SV_SHELF_MAX_DB),
            "bleed_max": float(SV_BLEED_MAX),
            "h2_loud": float(h2),
            "h2_quiet": float(h2 * 0.5),
            "h3_loud": float(h3),
            "h3_quiet": float(h3 * 0.5),
            "pre_gain_db_loud": round(pre_gain_db_loud, 1),
        },
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Audio processing: file → C enhancer → WAV
# ═══════════════════════════════════════════════════════════════════════════════

def process_track(
    input_path: str,
    output_path: str,
    coeffs: List[float],
    n_biquads: int,
    cutoff_hz: float = 60.0,
    h2: float = 0.5,
    h3: float = 1.0,
    pre_gain: float = 1.0,
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
        pre_gain=pre_gain,
        fs=44100.0, coeffs=coeffs)

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
