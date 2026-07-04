#!/usr/bin/env python3
"""
Speaker EQ correction suite for JamesDSP.

Replaces the Rust autoeq binary (src/main.rs) with a pure Python
implementation using scipy for signal processing and the existing
python/ modules for the psychoacoustic bass enhancer model.

Usage:
    python eqgen.py -m measurements/technics/standing/measurement2.wav \
                     -t measurements/technics/standing/target.wav \
                     --noise measurements/technics/standing/noise2.wav \
                     -o sanitycheck \
                     --bass-enhancer-cutoff 50.0

    python eqgen.py -m m1.wav m2.wav -t target.wav -o autoeq --max-noise 0.10
"""

import argparse
import sys
import os
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import welch, windows

# Add python/ to path for our modules
sys.path.insert(0, str(Path(__file__).parent / "python"))

from model import (
    butterworth_lp_mag,
    butterworth_hp_mag,
    ear_sensitivity,
    model_perceived_amplitude,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Constants (from src/main.rs)
# ═══════════════════════════════════════════════════════════════════════════════

BOTTOM_F = 20.0       # Hz — lowest frequency analyzed
TOP_F = 14000.0       # Hz — highest frequency analyzed
BASE_RESOLUTION = 8000.0  # base EQ point density before CV adaptation

WELCH_FFT_SIZE = 16384
WELCH_OVERLAP = 0.5



# ═══════════════════════════════════════════════════════════════════════════════
# WAV I/O
# ═══════════════════════════════════════════════════════════════════════════════

def read_wav(path):
    """Read a mono or stereo WAV, returning mono float64 samples in [-1, 1].

    Handles 8/16/24/32-bit integer and 32-bit float formats via scipy.
    """
    rate, data = wavfile.read(path)

    # Convert to float64 in [-1, 1]
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

    # Convert to mono
    if data.ndim == 1:
        pass  # already mono
    elif data.ndim == 2:
        data = data.mean(axis=1)
    else:
        raise ValueError(f"Unexpected WAV shape: {data.shape}")

    # Ensure contiguous float64 array
    samples = np.ascontiguousarray(data, dtype=np.float64)
    return samples, rate


# ═══════════════════════════════════════════════════════════════════════════════
# Welch's method FFT
# ═══════════════════════════════════════════════════════════════════════════════

def welch_stats(samples, sample_rate):
    """Run Welch's method and return per-bin mean (raw) and per-bin CV.

    Normalization removes common-mode gain variation between windows
    (mic distance drift) so CV only reflects frequency-specific variance.

    Returns list of (freq, count, mean, cv) dicts.
    """
    fft_size = WELCH_FFT_SIZE
    noverlap = int(fft_size * WELCH_OVERLAP)
    freqs, pxx = welch(
        samples, sample_rate,
        window=windows.hann(fft_size),
        nperseg=fft_size,
        noverlap=noverlap,
        return_onesided=True,
        scaling="spectrum",
    )
    # pxx is power spectral density; convert to linear magnitude
    mags = np.sqrt(pxx)
    n_windows = int(np.floor((len(samples) - noverlap) / (fft_size - noverlap)))

    # For CV: we need per-window normalized magnitudes. Re-run manually.
    step = fft_size - noverlap
    win = windows.hann(fft_size)
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
# dB conversions
# ═══════════════════════════════════════════════════════════════════════════════

def db_to_ratio(db):
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio):
    return 20.0 * np.log10(max(ratio, 1e-20))


# ═══════════════════════════════════════════════════════════════════════════════
# Graph: our basic data structure (replaced Rust Graph struct)
# ═══════════════════════════════════════════════════════════════════════════════

class Graph:
    """A frequency-response curve: list of (x=freq, y=value) points."""

    def __init__(self, points):
        # points: list of {"freq": f, "value": v} or {"x": x, "y": y}
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

    def point(self, x):
        """Linear interpolation to find y at given x."""
        if not self.points:
            return 0.0
        if x <= self.points[0]["x"]:
            return self.points[0]["y"]
        if x >= self.points[-1]["x"]:
            return self.points[-1]["y"]
        # Binary search
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

    def emit(self, x_values):
        """Return Graph with y values interpolated to exactly x_values."""
        result = []
        for i, x in enumerate(x_values):
            left = x_values[i - 1] if i > 0 else x
            right = x_values[i + 1] if i + 1 < len(x_values) else x
            x1 = (x + left) / 2.0
            x2 = (x + right) / 2.0
            # RMS average over 20 samples in the bin
            total = 0.0
            for j in range(20):
                xj = x1 + (x2 - x1) * j / 20.0
                yj = self.point(xj)
                total += yj * yj
            result.append({"x": x, "y": np.sqrt(total / 20.0)})
        return Graph(result)

    def rolloff(self, freq, db_per_octave, direction):
        """Apply a shelving rolloff above (High) or below (Low) freq."""
        atten = db_to_ratio(db_per_octave)
        new_points = []
        for p in self.points:
            ratio = p["x"] / freq if direction == "high" else freq / p["x"]
            octaves = np.log2(max(ratio, 1.0))
            new_points.append({"x": p["x"], "y": p["y"] * (atten ** octaves)})
        return Graph(new_points)

    def bass_enhancer_preprocess(self, fc, measurement, ramp_db=-60.0, h2=1.0, h3=1.0):
        """Preprocess EQ curve for the linearised harmonic bass enhancer.

        For each frequency f, computes the EQ gain G such that the
        perceived output A(f, G) matches the target curve.

        h2, h3 must match the enhancer plugin settings.

        Port of the Rust bass_enhancer_preprocess method.
        """
        h2sq = h2 * h2
        h3sq = h3 * h3

        # Equal-loudness: at bass frequencies the ear is ~5× more
        # sensitive to harmonics (100-150 Hz) than fundamentals (50 Hz).
        # Weight the harmonic terms so the model doesn't over-boost.
        EAR_W = 5.0
        EAR_W2 = EAR_W * EAR_W

        def solve(f, comp):
            if f > 2.0 * fc:
                return comp

            m_f = measurement.point(f)
            m_2f = measurement.point(2.0 * f)
            m_3f = measurement.point(3.0 * f)

            hp_f = butterworth_hp_mag(f, fc)
            hp_2f = butterworth_hp_mag(2.0 * f, fc)
            hp_3f = butterworth_hp_mag(3.0 * f, fc)

            a = (hp_f * m_f) ** 2
            b = h2sq * butterworth_lp_mag(f, fc) ** 2 * hp_2f ** 2 * m_2f ** 2
            c = h3sq * butterworth_lp_mag(f, fc / 2.0) ** 2 * hp_3f ** 2 * m_3f ** 2

            denom = a + EAR_W2 * b + EAR_W2 * c
            if denom < 1e-12:
                return 0.01
            target = comp * m_f
            return min(target / np.sqrt(denom), comp)

        new_points = [
            {"x": p["x"], "y": solve(p["x"], p["y"])}
            for p in self.points
        ]

        # Linear ramp below fc/2: gently attenuate sub-bass EQ gains
        # rather than hard-truncating at fc/3.
        ramp_start = max(20.0, fc / 3.0)
        ramp_end = fc / 2.0
        if ramp_end > ramp_start:
            ramp_ratio = db_to_ratio(ramp_db)
            for p in new_points:
                f = p["x"]
                if f <= ramp_start:
                    p["y"] *= ramp_ratio
                elif f < ramp_end:
                    t = (f - ramp_start) / (ramp_end - ramp_start)
                    atten_db = ramp_db * (1.0 - t)
                    p["y"] *= db_to_ratio(atten_db)

        return Graph(new_points)

    def max_value(self):
        return max(p["y"] for p in self.points) if self.points else 0.0

    def to_jamesdsp_eq(self, boost=0.0):
        """Format as JamesDSP desktop CSV: freq<TAB>gain."""
        dbs = [ratio_to_db(p["y"]) for p in self.points]
        max_db = max(dbs) if dbs else 0.0
        lines = []
        for i, p in enumerate(self.points):
            val = (dbs[i] - max_db + boost)
            lines.append(f"{p['x']:.3f}\t{min(val, 0.0):.3f}")
        return "\n".join(lines) + "\n"

    def to_jamesdsp_eq_android(self, boost=0.0):
        """Format as JamesDSP Android GraphicEQ string."""
        dbs = [ratio_to_db(p["y"]) for p in self.points]
        max_db = max(dbs) if dbs else 0.0
        parts = []
        for i, p in enumerate(self.points):
            val = (dbs[i] - max_db + boost)
            parts.append(f"{p['x']:.3f} {min(val, 0.0):.3f};")
        return "GraphicEQ: " + " ".join(parts)


# ═══════════════════════════════════════════════════════════════════════════════
# EQ point generation
# ═══════════════════════════════════════════════════════════════════════════════

def _erb_hz(freq):
    """Equivalent Rectangular Bandwidth at freq (Hz)."""
    return 24.7 * (4.37 * freq / 1000.0 + 1.0)


def eq_points_adaptive(pooled, max_noise, min_erb_fraction=0.11, max_span_octaves=2.0):
    """Generate EQ points greedily from pooled Welch FFT stats.

    Single forward pass. Accumulates Welch bins into spans. A span is
    finalized when (a) it is at least min_erb_fraction of an ERB wide,
    and (b) its pooled time-domain noise drops below max_noise.

    Guards:
      min_erb_fraction  – minimum span width as fraction of ERB
                           (prevents psychoacoustically redundant points)
      max_span_octaves  – force-finalize if span exceeds this width

    Clean regions → many narrow spans at ERB resolution.
    Noisy regions → few wide spans at max_span_octaves resolution.

    Higher max_noise → spans finalize sooner → more EQ points.
    """
    bin_hz = float(pooled[1]["freq"] - pooled[0]["freq"]) if len(pooled) > 1 else 1.0
    max_ratio = 2.0 ** max_span_octaves

    points = []

    # Running sums for the current span
    n   = 0.0   # Σ count_i
    sx  = 0.0   # Σ count_i · mean_i
    sx2 = 0.0   # Σ count_i · mean_i²
    sv  = 0.0   # Σ count_i · (noise_i · mean_i)²

    span_start = 0
    bins_in_span = 0

    def finalize(end_idx):
        """Commit span [span_start, end_idx] as an EQ point, reset acc."""
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

        # ── max-span guard: if adding this bin would make the span
        #    too wide, force-finalize the current span first ──
        f_lo = BOTTOM_F if span_start == 0 else pooled[span_start]["freq"] - bin_hz / 2.0
        f_hi = pooled[i]["freq"] + bin_hz / 2.0
        if bins_in_span > 0 and f_hi / f_lo > max_ratio:
            finalize(i - 1)
            # fall through: start new span with this bin

        # ── accept this bin into the span ──
        n   += c
        sx  += c * m
        sx2 += c * m * m
        sv  += c * var_i
        bins_in_span += 1

        # ── noise check: if span is wide enough AND clean enough,
        #    finalize now (data is trustworthy at this resolution) ──
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

    # ── finalize trailing span ──
    if bins_in_span > 0:
        finalize(len(pooled) - 1)

    return np.array(points)


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def parse_rolloff(spec):
    """Parse a "FREQ,DB_OCTAVE" rolloff spec."""
    parts = spec.split(",")
    freq = float(parts[0])
    db = float(parts[1])
    return freq, db


def main():
    parser = argparse.ArgumentParser(
        description="Speaker EQ correction suite for JamesDSP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python eqgen.py -m meas1.wav meas2.wav -t target.wav -o autoeq
    python eqgen.py -m meas.wav -t target.wav --noise noise.wav -b 3.0 \\
        --bass-enhancer-cutoff 50.0 \\
        --low-rolloff 220,2 --low-rolloff 50,-2
""",
    )

    parser.add_argument(
        "-m", "--measurement", nargs="+", required=True,
        help="One or more measurement WAV files (headphones, brown noise)",
    )
    parser.add_argument(
        "-t", "--target", required=True,
        help="Target WAV file (reference system recording)",
    )
    parser.add_argument(
        "-n", "--noise", default=None,
        help="Background noise WAV file for spectral subtraction",
    )
    parser.add_argument(
        "-o", "--output", default="autoeq",
        help="Output prefix (generates <prefix>_desktop.csv and <prefix>_mobile.csv)",
    )
    parser.add_argument(
        "--high-rolloff", action="append", default=[], metavar="FREQ,DB_OCTAVE",
        help="High-frequency rolloff (e.g. 4000,-6). Repeatable.",
    )
    parser.add_argument(
        "--low-rolloff", action="append", default=[], metavar="FREQ,DB_OCTAVE",
        help="Low-frequency rolloff (e.g. 100,-12). Repeatable.",
    )
    parser.add_argument(
        "-b", "--boost", type=float, default=0.0,
        help="Boost in dB to apply to the final EQ [default: 0.0]",
    )
    parser.add_argument(
        "--max-noise", type=float, default=0.65,
        help="Maximum per-span noise (coeff. of variation). Higher tolerates more uncertainty → fewer EQ points. Lower insists on cleaner data → more EQ points. [default: 0.65]",
    )
    parser.add_argument(
        "--bass-enhancer-cutoff", type=float, default=None,
        help="Cutoff frequency for harmonic bass enhancer compensation.",
    )
    parser.add_argument(
        "--h2", type=float, default=1.0,
        help="2nd harmonic amplitude (must match plugin) [default: 1.0]",
    )
    parser.add_argument(
        "--h3", type=float, default=1.0,
        help="3rd harmonic amplitude (must match plugin) [default: 1.0]",
    )

    args = parser.parse_args()

    # ── 1. Read & validate all WAV files ──────────────────────────────
    print("Reading WAV files...", file=sys.stderr)
    wavs = []
    for path in args.measurement:
        samples, rate = read_wav(path)
        wavs.append((samples, rate, path))
        print(f"  {path}: {len(samples)} samples @ {rate} Hz", file=sys.stderr)

    target_samples, target_rate = read_wav(args.target)
    print(f"  {args.target}: {len(target_samples)} samples @ {target_rate} Hz", file=sys.stderr)

    sample_rate = wavs[0][1]
    for s, r, p in wavs[1:]:
        assert r == sample_rate, f"Sample rate mismatch in {p}: {r} vs {sample_rate}"
    assert target_rate == sample_rate, \
        f"Target sample rate mismatch: {target_rate} vs {sample_rate}"

    # ── 2. Welch FFT each measurement → per-bin stats ─────────────────
    print("Computing FFT (Welch's method)...", file=sys.stderr)
    all_stats = [welch_stats(s, sample_rate) for s, _, _ in wavs]

    bin_count = len(all_stats[0])
    for s in all_stats[1:]:
        assert len(s) == bin_count, "Frequency ranges differ — sample rates must match"

    # ── 3. Pool per-bin stats across measurements ─────────────────────
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

    print(f"  {pooled[0]['count']} windows across {len(wavs)} measurement(s)", file=sys.stderr)

    # ── 4. Build full-resolution graphs ───────────────────────────────
    measurement_full = Graph(pooled)

    # ── 5. Noise: spectral subtraction + merge noise CV into pooled ───
    noise_full_res = None
    if args.noise:
        noise_samples, noise_rate = read_wav(args.noise)
        assert noise_rate == sample_rate
        noise_stats = welch_stats(noise_samples, sample_rate)

        # Merge noise CV into pooled stats so pruning sees worst-case CV
        for i in range(min(len(pooled), len(noise_stats))):
            pooled[i]["noise"] = max(pooled[i]["noise"], noise_stats[i]["noise"])

        noise_full_res = Graph(noise_stats)

    # ── 6. Generate adaptive EQ points ────────────────────────────────
    eq_pts = eq_points_adaptive(pooled, args.max_noise)
    print(f"  {len(eq_pts)} EQ points (max-noise={args.max_noise})", file=sys.stderr)

    # ── 7. Resample measurement + spectral subtraction ─────────────────
    measurement_mean = measurement_full.emit(eq_pts)

    if noise_full_res is not None:
        noise_mean = noise_full_res.emit(eq_pts)
        for i in range(len(measurement_mean.points)):
            m = measurement_mean.points[i]["y"]
            n = noise_mean.points[i]["y"]
            measurement_mean.points[i]["y"] = max(m - n, m * 0.01)

    # ── 8. Target: Welch → full-res → resample → rolloffs ────────────
    target_stats = welch_stats(target_samples, sample_rate)
    target_full = Graph(target_stats)
    target = target_full.emit(eq_pts)

    for spec in args.high_rolloff:
        freq, db = parse_rolloff(spec)
        target = target.rolloff(freq, db, "high")
    for spec in args.low_rolloff:
        freq, db = parse_rolloff(spec)
        target = target.rolloff(freq, db, "low")

    # ── 9. Compute correction = target / measurement ──────────────────
    comp_points = []
    for i, x in enumerate(eq_pts):
        t = target.points[i]["y"]
        m = measurement_mean.points[i]["y"]
        comp_points.append({"x": x, "y": t / m if m > 1e-20 else 1e6})
    comp = Graph(comp_points)

    # ── 10. Bass enhancer preprocessing ────────────────────────────────
    if args.bass_enhancer_cutoff is not None:
        comp = comp.bass_enhancer_preprocess(
            args.bass_enhancer_cutoff, measurement_mean,
            h2=args.h2, h3=args.h3,
        )

    # ── 11. Write outputs ──────────────────────────────────────────────
    desktop_path = f"{args.output}_desktop.csv"
    with open(desktop_path, "w") as f:
        f.write(comp.to_jamesdsp_eq(args.boost))
    print(f"Wrote {desktop_path}", file=sys.stderr)

    mobile_path = f"{args.output}_mobile.csv"
    with open(mobile_path, "w") as f:
        f.write(comp.to_jamesdsp_eq_android(args.boost))
    print(f"Wrote {mobile_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
