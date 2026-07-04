#!/usr/bin/env python3
"""
3D plot: EQ correction vs frequency vs bass enhancer cutoff.

Sweeps the cutoff frequency from 50 Hz to 150 Hz, re-running the
full correction pipeline at each step, and plots the resulting
EQ gain surface.
"""

import sys
from pathlib import Path
import numpy as np
from scipy.io import wavfile
from scipy.signal import welch, windows

sys.path.insert(0, str(Path(__file__).parent.parent))

from model import butterworth_lp_mag, butterworth_hp_mag

# ═══════════════════════════════════════════════════════════════════════
# Constants (mirror eqgen.py)
# ═══════════════════════════════════════════════════════════════════════

BOTTOM_F = 20.0
TOP_F = 14000.0
BASE_RESOLUTION = 8000.0
H2_AMP = 0.33
H3_AMP = 0.33
WELCH_FFT_SIZE = 16384
WELCH_OVERLAP = 0.5


# ═══════════════════════════════════════════════════════════════════════
# Helpers (mirrored from eqgen.py to keep the script self-contained)
# ═══════════════════════════════════════════════════════════════════════

def db_to_ratio(db):
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio):
    return 20.0 * np.log10(max(ratio, 1e-20))


def read_wav(path):
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
    return np.ascontiguousarray(data, dtype=np.float64), rate


def welch_stats(samples, sample_rate):
    fft_size = WELCH_FFT_SIZE
    noverlap = int(fft_size * WELCH_OVERLAP)
    step = fft_size - noverlap
    win = windows.hann(fft_size)
    bin_count = fft_size // 2 + 1

    raw_sum = np.zeros(bin_count)
    norm_sum = np.zeros(bin_count)
    norm_sum_sq = np.zeros(bin_count)
    window_count = 0

    for start in range(0, len(samples) - fft_size + 1, step):
        segment = samples[start : start + fft_size] * win
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
        cv = (np.sqrt(norm_var) / norm_mean) if norm_mean > 0 else float("inf")
        results.append(
            {"freq": f, "count": window_count, "mean": mean_val, "cv": cv}
        )
    return results


class Graph:
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

    def point(self, x):
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
        p1, p2 = self.points[lo], self.points[hi]
        slope = (p2["y"] - p1["y"]) / (p2["x"] - p1["x"])
        return p1["y"] + (x - p1["x"]) * slope

    def emit(self, x_values):
        result = []
        for i, x in enumerate(x_values):
            left = x_values[i - 1] if i > 0 else x
            right = x_values[i + 1] if i + 1 < len(x_values) else x
            x1, x2 = (x + left) / 2.0, (x + right) / 2.0
            total = 0.0
            for j in range(20):
                xj = x1 + (x2 - x1) * j / 20.0
                yj = self.point(xj)
                total += yj * yj
            result.append({"x": x, "y": np.sqrt(total / 20.0)})
        return Graph(result)

    def rolloff(self, freq, db_per_octave, direction):
        atten = db_to_ratio(db_per_octave)
        new_points = []
        for p in self.points:
            ratio = p["x"] / freq if direction == "high" else freq / p["x"]
            octaves = np.log2(max(ratio, 1.0))
            new_points.append({"x": p["x"], "y": p["y"] * (atten ** octaves)})
        return Graph(new_points)


def eq_points_adaptive(cv_graph, target_cv):
    equal_loudness = Graph(
        [
            {"x": 20.0, "y": 109.0},
            {"x": 80.0, "y": 82.0},
            {"x": 400.0, "y": 62.0},
            {"x": 1000.0, "y": 60.0},
            {"x": 1500.0, "y": 64.0},
            {"x": 2500.0, "y": 57.0},
            {"x": 4000.0, "y": 57.0},
            {"x": 8500.0, "y": 73.0},
            {"x": 15000.0, "y": 72.0},
            {"x": 19000.0, "y": 68.0},
            {"x": 30000.0, "y": 130.0},
        ]
    )

    def base_density(freq):
        return BASE_RESOLUTION / equal_loudness.point(freq)

    freq = BOTTOM_F
    out = []
    while freq < TOP_F:
        out.append(freq)
        base_step = freq / base_density(freq) * 2.0 / np.e
        cv = max(cv_graph.point(freq), 0.001)
        freq += base_step * (cv / target_cv)
    out.append(TOP_F)
    return np.array(out)


def bass_enhancer_preprocess(comp_graph, measurement, fc, ramp_db=-12.0):
    """Standalone version of the preprocess, returning a list of (freq, gain_dB)."""
    h2sq = H2_AMP**2
    h3sq = H3_AMP**2

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
        b = h2sq * butterworth_lp_mag(f, fc) ** 2 * hp_2f**2 * m_2f**2
        c = h3sq * butterworth_lp_mag(f, fc / 2.0) ** 2 * hp_3f**2 * m_3f**2

        denom = a + b + c
        if denom < 1e-12:
            return 0.01
        target = comp * m_f
        return min(target / np.sqrt(denom), comp)

    points = []
    for p in comp_graph.points:
        points.append({"x": p["x"], "y": solve(p["x"], p["y"])})

    # Ramp (mirrors eqgen.py Graph.bass_enhancer_preprocess)
    ramp_start = max(20.0, fc / 3.0)
    ramp_end = fc / 2.0
    if ramp_end > ramp_start:
        ramp_ratio = db_to_ratio(ramp_db)
        for p in points:
            f = p["x"]
            if f <= ramp_start:
                p["y"] *= ramp_ratio
            elif f < ramp_end:
                t = (f - ramp_start) / (ramp_end - ramp_start)
                atten_db = ramp_db * (1.0 - t)
                p["y"] *= db_to_ratio(atten_db)

    # Normalize and format as dB
    dbs = [ratio_to_db(p["y"]) for p in points]
    max_db = max(dbs) if dbs else 0.0
    return [(p["x"], dbs[i] - max_db) for i, p in enumerate(points)]


# ═══════════════════════════════════════════════════════════════════════
# Main sweep
# ═══════════════════════════════════════════════════════════════════════

def main():
    import matplotlib.pyplot as plt
    from matplotlib import cm

    ROOT = Path(__file__).parent.parent.parent
    meas_path = str(ROOT / "measurements/technics/standing/measurement2.wav")
    target_path = str(ROOT / "measurements/technics/standing/target.wav")
    noise_path = str(ROOT / "measurements/technics/standing/noise2.wav")

    print("Reading WAVs...", file=sys.stderr)
    meas_samples, rate = read_wav(meas_path)
    target_samples, _ = read_wav(target_path)
    noise_samples, _ = read_wav(noise_path)

    # ── One-time measurements ─────────────────────────────────────────
    print("Welch stats...", file=sys.stderr)
    stats = welch_stats(meas_samples, rate)
    measurement_full = Graph(stats)
    cv_full = Graph([{"freq": p["freq"], "y": p["cv"]} for p in stats])

    # Noise CV merge
    noise_stats = welch_stats(noise_samples, rate)
    for i in range(min(len(cv_full.points), len(noise_stats))):
        cv_full.points[i]["y"] = max(cv_full.points[i]["y"], noise_stats[i]["cv"])

    noise_full = Graph(noise_stats)

    # EQ points
    eq_pts = eq_points_adaptive(cv_full, target_cv=0.12)
    print(f"{len(eq_pts)} EQ points", file=sys.stderr)

    # Resample measurement + noise subtract
    measurement_mean = measurement_full.emit(eq_pts)
    noise_mean = noise_full.emit(eq_pts)
    for i in range(len(measurement_mean.points)):
        m = measurement_mean.points[i]["y"]
        n = noise_mean.points[i]["y"]
        measurement_mean.points[i]["y"] = max(m - n, m * 0.01)

    # Target
    target_stats = welch_stats(target_samples, rate)
    target = Graph(target_stats).emit(eq_pts)

    # ── Sweep cutoff frequencies ──────────────────────────────────────
    cutoffs = np.arange(50, 155, 5)  # 50, 55, ..., 150
    all_curves = {}  # cutoff → [(freq, dB), ...]

    for fc in cutoffs:
        print(f"Cutoff {fc:.0f} Hz...", file=sys.stderr)

        # Compute bare correction
        comp_points = []
        for i in range(len(eq_pts)):
            t = target.points[i]["y"]
            m = measurement_mean.points[i]["y"]
            comp_points.append({"x": eq_pts[i], "y": t / m if m > 1e-20 else 1e6})
        comp = Graph(comp_points)

        # Apply enhancer preprocess
        curve = bass_enhancer_preprocess(comp, measurement_mean, float(fc))
        all_curves[float(fc)] = curve

    # ── Build 3D surface data ─────────────────────────────────────────
    # Use a common frequency grid for all curves (the EQ points themselves)
    freqs = eq_pts
    Z = np.zeros((len(cutoffs), len(freqs)))
    for i, fc in enumerate(cutoffs):
        curve_dict = dict(all_curves[float(fc)])
        for j, f in enumerate(freqs):
            Z[i, j] = curve_dict.get(f, np.nan)

    X, Y = np.meshgrid(freqs, cutoffs)

    # ── Plot ───────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(
        X, Y, Z,
        cmap=cm.viridis,
        alpha=0.92,
        linewidth=0,
        antialiased=True,
        vmin=-50,
        vmax=0,
    )

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Cutoff (Hz)")
    ax.set_zlabel("EQ Gain (dB, normalized)")
    ax.set_title("EQ Correction vs Bass Enhancer Cutoff\n(Technics standing, sanity check data)")

    ax.set_xlim(20, 500)
    ax.set_ylim(50, 150)
    ax.set_zlim(-50, 0)

    ax.view_init(elev=25, azim=-55)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0f}"))

    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=15, pad=0.08)
    cbar.set_label("dB (rel. max)")

    outpath = str(ROOT / "python/archive/cutoff_sweep.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nSaved {outpath}", file=sys.stderr)

    # Also print text summary for the bass region
    print("\n── Bass region (20–200 Hz) at selected cutoffs ──")
    print(f"{'Freq':>8s}", end="")
    for fc in [50, 70, 90, 110, 130, 150]:
        print(f"  {'fc=' + str(fc):>10s}", end="")
    print()
    for j, f in enumerate(freqs):
        if f > 200:
            break
        print(f"{f:8.1f}", end="")
        for i, fc in enumerate(cutoffs):
            if fc in [50, 70, 90, 110, 130, 150]:
                print(f"  {Z[i, j]:+9.2f}", end="")
        print()


if __name__ == "__main__":
    main()
