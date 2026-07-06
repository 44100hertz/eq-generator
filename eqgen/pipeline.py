"""
Shared pipeline: Welch FFT → adaptive EQ points → correction → IIR fit → C DSP.

This module is the single source of truth for the EQ design and audio processing
pipeline.  All CLI entry points (eqgen, audition, wire, export) delegate here.

Exports:
  run_pipeline()   — WAV inputs → (freqs, gains_db, sample_rate)
  design_eq()      — EQ curve → quantized biquad coefficients
  process_track()  — audio file → C enhancer → WAV
  Graph            — frequency-response curve with interpolation and transforms
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


# ═══════════════════════════════════════════════════════════════════════════════
# dB conversion helpers
# ═══════════════════════════════════════════════════════════════════════════════

def db_to_ratio(db: float) -> float:
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio: float) -> float:
    return 20.0 * np.log10(max(ratio, 1e-20))


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
    """Run Welch's method and return per-bin mean and per-bin CV."""
    fft_size = WELCH_FFT_SIZE
    noverlap = int(fft_size * WELCH_OVERLAP)

    step = fft_size - noverlap
    win = np.hanning(fft_size)
    bin_count = fft_size // 2 + 1

    raw_sum = np.zeros(bin_count)
    norm_sum = np.zeros(bin_count)
    norm_sum_sq = np.zeros(bin_count)
    window_count = 0

    for start in range(0, len(samples) - fft_size + 1, step):
        segment = samples[start:start + fft_size] * win
        spectrum = np.fft.rfft(segment)
        mags_window = np.abs(spectrum) / fft_size
        raw_sum += mags_window

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
        norm_mean = norm_sum[i] / n
        norm_var = (norm_sum_sq[i] / n) - (norm_mean * norm_mean)
        cv = (np.sqrt(norm_var) / norm_mean) if norm_mean > 0 else float('inf')
        results.append({
            "freq": f,
            "count": window_count,
            "mean": mean_val,
            "noise": cv,
        })
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Graph: frequency-response curve
# ═══════════════════════════════════════════════════════════════════════════════

class Graph:
    """A frequency-response curve: list of (x=freq, y=value) points."""

    def __init__(self, points):
        if isinstance(points, list) and len(points) > 0:
            if isinstance(points[0], dict):
                self.points = [
                    {"x": p.get("freq", p.get("x")), "y": p.get("mean", p.get("y"))}
                    for p in points
                ]
            else:
                self.points = [{"x": p[0], "y": p[1]} for p in points]
        else:
            self.points = []

    def point(self, x: float) -> float:
        """Linear interpolation to find y at given x."""
        if not self.points:
            return 0.0
        if x <= self.points[0]["x"]:
            return self.points[0]["y"]
        if x >= self.points[-1]["x"]:
            return self.points[-1]["y"]
        lo, hi = 0, len(self.points) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.points[mid]["x"] <= x:
                lo = mid
            else:
                hi = mid
        p1 = self.points[lo]
        p2 = self.points[hi]
        slope = (p2["y"] - p1["y"]) / (p2["x"] - p1["x"])
        return p1["y"] + (x - p1["x"]) * slope

    def emit(self, x_values: np.ndarray) -> "Graph":
        """Return Graph with y values interpolated to exactly x_values."""
        result = []
        for i, x in enumerate(x_values):
            left = x_values[i - 1] if i > 0 else x
            right = x_values[i + 1] if i + 1 < len(x_values) else x
            x1 = (x + left) / 2.0
            x2 = (x + right) / 2.0
            total = 0.0
            for j in range(20):
                xj = x1 + (x2 - x1) * j / 20.0
                yj = self.point(xj)
                total += yj * yj
            result.append({"x": x, "y": np.sqrt(total / 20.0)})
        return Graph(result)

    def rolloff(self, freq: float, db_per_octave: float, direction: str) -> "Graph":
        """Apply a shelving rolloff above ('high') or below ('low') freq."""
        atten = db_to_ratio(db_per_octave)
        new_points = []
        for p in self.points:
            ratio = p["x"] / freq if direction == "high" else freq / p["x"]
            octaves = np.log2(max(ratio, 1.0))
            new_points.append({"x": p["x"], "y": p["y"] * (atten ** octaves)})
        return Graph(new_points)

    def bass_enhancer_preprocess(self, fc: float, measurement: "Graph",
                                  h2: float = 1.0, h3: float = 1.0) -> "Graph":
        """Preprocess EQ curve for the harmonic bass enhancer.

        For each frequency f, computes the EQ gain G such that the
        perceived output A(f, G) matches the target curve.
        Delegates to model.model_gain_needed for the core solve.
        """
        from eqgen.model import model_gain_needed

        def solve(f: float, comp: float) -> float:
            if f > 2.0 * fc:
                return comp
            m_f = measurement.point(f)
            G = model_gain_needed(f, comp * m_f, fc, h2, h3, measurement.point)
            if G < 1e-12:
                return 0.01
            return min(G, comp)

        new_points = [
            {"x": p["x"], "y": solve(p["x"], p["y"])}
            for p in self.points
        ]

        # Below fc/2 the enhancer handles bass perception via harmonics.
        # Hold correction flat at the model gain computed AT fc/2.
        flat_val = solve(fc / 2.0, self.point(fc / 2.0))
        for p in new_points:
            if p["x"] <= fc / 2.0:
                p["y"] = flat_val

        return Graph(new_points)

    def max_value(self) -> float:
        return max(p["y"] for p in self.points) if self.points else 0.0

    def to_json(self) -> list:
        """Serialize as list of {freq, gain_db} for JSON output."""
        return [
            {"freq": round(p["x"], 2), "gain_db": round(ratio_to_db(p["y"]), 3)}
            for p in self.points
        ]


# ═══════════════════════════════════════════════════════════════════════════════
# Adaptive EQ point generation
# ═══════════════════════════════════════════════════════════════════════════════

def _erb_hz(freq: float) -> float:
    """Equivalent Rectangular Bandwidth at freq (Hz)."""
    return 24.7 * (4.37 * freq / 1000.0 + 1.0)


def eq_points_adaptive(pooled: List[dict], max_noise: float,
                        min_erb_fraction: float = 0.11,
                        max_span_octaves: float = 2.0) -> np.ndarray:
    """Generate EQ points greedily from pooled Welch FFT stats."""
    bin_hz = float(pooled[1]["freq"] - pooled[0]["freq"]) if len(pooled) > 1 else 1.0
    max_ratio = 2.0 ** max_span_octaves

    points = []
    n = sx = sx2 = sv = 0.0
    span_start = 0
    bins_in_span = 0

    def finalize(end_idx: int):
        nonlocal n, sx, sx2, sv, span_start, bins_in_span
        start = span_start
        f_lo = BOTTOM_F if start == 0 else pooled[start]["freq"] - bin_hz / 2.0
        f_hi = TOP_F if end_idx >= len(pooled) - 1 else pooled[end_idx]["freq"] + bin_hz / 2.0
        points.append((f_lo + f_hi) / 2.0)
        span_start = end_idx + 1
        n = sx = sx2 = sv = 0.0
        bins_in_span = 0

    for i, b in enumerate(pooled):
        c = float(b["count"])
        m = b["mean"]
        var_i = (b["noise"] * m) ** 2

        f_lo = BOTTOM_F if span_start == 0 else pooled[span_start]["freq"] - bin_hz / 2.0
        f_hi = pooled[i]["freq"] + bin_hz / 2.0
        if bins_in_span > 0 and f_hi / f_lo > max_ratio:
            finalize(i - 1)

        n += c
        sx += c * m
        sx2 += c * m * m
        sv += c * var_i
        bins_in_span += 1

        span_hz = f_hi - f_lo
        center = (f_lo + f_hi) / 2.0
        min_hz = _erb_hz(center) * min_erb_fraction
        if span_hz >= min_hz:
            wm = sx / n
            if wm > 1e-20:
                within_var = sv / n
                noise = np.sqrt(within_var) / wm
                if noise <= max_noise:
                    finalize(i)

    if bins_in_span > 0:
        finalize(len(pooled) - 1)

    return np.array(points)


# ═══════════════════════════════════════════════════════════════════════════════
# Core pipeline: WAV inputs → EQ curve
# ═══════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    measurement_paths: List[str],
    target_path: str,
    noise_path: Optional[str] = None,
    max_noise: float = 0.65,
    bass_enhancer_cutoff: Optional[float] = None,
    h2: float = 1.0,
    h3: float = 1.0,
    high_rolloffs: Optional[List[Tuple[float, float]]] = None,
    low_rolloffs: Optional[List[Tuple[float, float]]] = None,
    sample_rate_override: Optional[float] = None,
    detailed: bool = False,
):
    """Run the full EQ correction pipeline.

    Without detailed: returns (freqs_hz, gains_db, sample_rate).
    With detailed=True: returns a dict with all intermediate graphs for visualization.

    Steps:
      1. Read & validate WAV files
      2. Welch FFT on each measurement
      3. Pool per-bin stats across measurements
      4. Merge noise floor CV (if noise file provided)
      5. Adaptive EQ points (noise-weighted frequency grid)
      6. Resample measurement + spectral subtraction
      7. Resample target + rolloffs
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

    # 3. Pool stats
    pooled = []
    for i in range(bin_count):
        total_n = sum(s[i]["count"] for s in all_stats)
        total_sum = sum(s[i]["mean"] * s[i]["count"] for s in all_stats)
        pooled_noise = sum(s[i]["noise"] * s[i]["count"] for s in all_stats) / total_n
        pooled.append({
            "freq": all_stats[0][i]["freq"],
            "count": total_n,
            "mean": total_sum / total_n,
            "noise": pooled_noise,
        })

    # 4. Build graph
    measurement_full = Graph(pooled)

    # 5. Noise merge
    noise_full_res = None
    if noise_path:
        noise_samples, noise_rate = read_wav(noise_path)
        assert noise_rate == sample_rate
        noise_stats = welch_stats(noise_samples, sample_rate)
        for i in range(min(len(pooled), len(noise_stats))):
            pooled[i]["noise"] = max(pooled[i]["noise"], noise_stats[i]["noise"])
        noise_full_res = Graph(noise_stats)

    # 6. Adaptive EQ points
    eq_pts = eq_points_adaptive(pooled, max_noise)

    # 7. Resample measurement + spectral subtraction
    measurement_mean = measurement_full.emit(eq_pts)
    if noise_full_res is not None:
        noise_mean = noise_full_res.emit(eq_pts)
        for i in range(len(measurement_mean.points)):
            m = measurement_mean.points[i]["y"]
            n = noise_mean.points[i]["y"]
            measurement_mean.points[i]["y"] = max(m - n, m * 0.01)

    # 8. Target
    target_stats = welch_stats(target_samples, sample_rate)
    target_full = Graph(target_stats)
    target = target_full.emit(eq_pts)
    for freq, db in high_rolloffs:
        target = target.rolloff(freq, db, "high")
    for freq, db in low_rolloffs:
        target = target.rolloff(freq, db, "low")

    # 9. Correction
    comp_points = []
    for i, x in enumerate(eq_pts):
        t = target.points[i]["y"]
        m = measurement_mean.points[i]["y"]
        comp_points.append({"x": x, "y": t / m if m > 1e-20 else 1e6})
    comp = Graph(comp_points)

    # 10. Bass enhancer preprocessing
    if bass_enhancer_cutoff is not None:
        comp = comp.bass_enhancer_preprocess(
            bass_enhancer_cutoff, measurement_mean, h2=h2, h3=h3)

    freqs = np.array([p["x"] for p in comp.points])
    gains_db = np.array([ratio_to_db(p["y"]) for p in comp.points])

    if detailed:
        # Gather all intermediate graphs for visualization
        raw_resp = []
        for p in measurement_full.points:
            raw_resp.append({"freq": p["x"], "db": ratio_to_db(p["y"])})
        raw_target = []
        for p in target_full.points:
            raw_target.append({"freq": p["x"], "db": ratio_to_db(p["y"])})
        noise_data = []
        if noise_full_res is not None:
            for p in noise_full_res.points:
                noise_data.append({"freq": p["x"], "cv": p["y"]})
        # Also get per-bin CV from pooled for the raw measurement noise
        pooled_noise = [{"freq": p["freq"], "cv": p["noise"]} for p in pooled]
        meas_sub = []
        for p in measurement_mean.points:
            meas_sub.append({"freq": p["x"], "db": ratio_to_db(p["y"])})
        targ_sub = []
        for p in target.points:
            targ_sub.append({"freq": p["x"], "db": ratio_to_db(p["y"])})
        return {
            "freqs": freqs.tolist(),
            "gains_db": gains_db.tolist(),
            "sample_rate": sample_rate,
            "raw_measurement": raw_resp,
            "raw_target": raw_target,
            "noise_cv": pooled_noise,
            "noise_floor": noise_data,
            "meas_subtracted": meas_sub,
            "target_resampled": targ_sub,
            "eq_points_freqs": eq_pts.tolist(),
        }

    return freqs, gains_db, sample_rate


def pipeline_to_graph(freqs: np.ndarray, gains_db: np.ndarray) -> Graph:
    """Convert pipeline output arrays back to a Graph for serialization."""
    return Graph([{"x": f, "y": db_to_ratio(g)} for f, g in zip(freqs, gains_db)])


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
