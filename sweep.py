#!/usr/bin/env python3
"""
Sweep FFT bin count × IIR biquad count to find the Pareto frontier
of perceptual error vs CPU load, identifying configurations that
won't cause underruns on ESP32.

Usage:
    python sweep.py                          # use default preset (lunchbox)
    python sweep.py --preset lunchbox        # use named preset
    python sweep.py --preset cardboard
    python sweep.py --measurements meas.wav --target target.wav
    python sweep.py --esp32                  # also run on ESP32 hardware
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SRC_DIR = ROOT / "src"
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import (
    run_pipeline, compute_fft_residual, design_eq, FFT_N as DEFAULT_FFT_N,
    FFT_RESIDUAL_THRESHOLD_DB,
)
from eqgen.eq_fit import cascade_response_db, BiquadCoeffs
from eqgen.quantize import q28_to_float, quantize_biquads_q28
from eqgen.presets import PresetManager, MAX_IIR_BANDS

# ── Sweep grid ───────────────────────────────────────────────────────
FFT_SIZES = [64, 128, 256, 512]       # powers of 2
BIQUAD_COUNTS = [0, 2, 4, 6, 8, 10, 12, 16, 20, 24]
# Measure biquad cost independently (not FFT-dependent in C impl)
BIQUAD_BENCH_COUNTS = [0, 2, 4, 6, 8, 10, 12, 16, 20, 24, 28, 32]
N_FRAMES = 1000  # fewer frames for quick sweep

# Perception-weighted bands (dB/octave response — bass errors matter more)
ERROR_BANDS = [
    ("Sub-bass (20–60)", 20, 60, 1.5),
    ("Bass (60–250)", 60, 250, 1.2),
    ("Mid (250–2k)", 250, 2000, 1.0),
    ("Treble (2k–8k)", 2000, 8000, 0.8),
    ("Air (>8k)", 8000, 20000, 0.5),
]


def build_bench_sweep(fft_n):
    """Build bench_sweep for a specific FFT_N. Returns the binary path."""
    binary = SRC_DIR / f"bench_sweep_{fft_n}"
    cmd = [
        "gcc", "-O2", "-Wall",
        f"-DFFT_EQ_N={fft_n}",
        "-I", str(SRC_DIR),
        "-o", str(binary),
        str(SRC_DIR / "bench_sweep.c"),
        str(SRC_DIR / "fft_eq.c"),
        str(SRC_DIR / "enhancer.c"),
        "-lm",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  ERROR building bench_sweep for FFT_N={fft_n}:")
        print(r.stderr)
        return None
    return binary


def run_bench(fft_n, n_biquads, n_frames=N_FRAMES):
    """Run bench_sweep for given FFT size and biquad count. Returns dict."""
    binary = build_bench_sweep(fft_n)
    if binary is None:
        return None
    r = subprocess.run(
        [str(binary), str(n_biquads), str(n_frames)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(f"  ERROR running bench_sweep FFT={fft_n} bq={n_biquads}")
        return None
    try:
        return json.loads(r.stdout.strip())
    except json.JSONDecodeError:
        print(f"  Failed to parse bench output for FFT={fft_n} bq={n_biquads}")
        return None


def sweep_errors(eval_freqs, target_db, fs):
    """Sweep FFT_N × n_biquads and compute prediction error metrics.

    Returns a list of dicts with fft_n, n_biquads, and error metrics.
    """
    results = []
    for fft_n in FFT_SIZES:
        # Compute FFT residual for this FFT size
        fft_bin_freqs, fft_gains_linear, fft_approx_db, residual_db = \
            compute_fft_residual(eval_freqs, target_db, fs, fft_n)

        max_residual = float(np.max(np.abs(residual_db)))

        for n_bq in BIQUAD_COUNTS:
            # Fit IIR biquads to the residual
            if n_bq == 0:
                # No IIR — error = residual itself
                combined_err = residual_db.copy()
            else:
                try:
                    coeffs, bands, fit_freqs, fit_target, fitted_db = design_eq(
                        eval_freqs, target_db, fs,
                        max_bands=n_bq, fft_n=fft_n)
                    # Combined response = FFT approx + IIR fitted
                    combined_db = fft_approx_db + fitted_db
                    combined_err = target_db - combined_db
                except Exception as e:
                    print(f"  fit failed FFT={fft_n} bq={n_bq}: {e}")
                    continue

            # Compute error metrics
            # RMSE
            rmse = float(np.sqrt(np.mean(combined_err ** 2)))

            # Per-band errors
            band_errors = {}
            for label, lo, hi, weight in ERROR_BANDS:
                mask = (eval_freqs >= lo) & (eval_freqs <= hi)
                if mask.sum() > 0:
                    mae = float(np.mean(np.abs(combined_err[mask])))
                    maxe = float(np.max(np.abs(combined_err[mask])))
                    band_errors[label] = {
                        "mae": mae,
                        "max_err": maxe,
                        "weighted_mae": mae * weight,
                    }

            # Weighted perceptual error score
            weighted_score = sum(
                band_errors[label]["weighted_mae"]
                for label in band_errors
            ) / len(ERROR_BANDS) if band_errors else rmse

            results.append({
                "fft_n": fft_n,
                "n_biquads": n_bq,
                "rmse_db": rmse,
                "weighted_error": weighted_score,
                "max_residual_db": max_residual,
                "band_errors": band_errors,
            })

    return results


def sweep_cpu():
    """Benchmark CPU for each FFT_N and biquad count combo on desktop.

    Returns (biquad_costs, fft_costs) dicts mapping size→us_per_frame.
    """
    print("\n── Benchmarking biquad cascade CPU cost...")
    bq_costs = {}
    # Biquad cost is independent of FFT_N — measure once
    for n_bq in BIQUAD_BENCH_COUNTS:
        r = run_bench(256, n_bq)
        if r:
            bq_costs[n_bq] = r["us_bq_frame"]
            print(f"  {n_bq:2d} biquads: {r['us_bq_frame']:6.1f} µs/frame")

    print("\n── Benchmarking FFT EQ CPU cost at each size...")
    fft_costs = {}
    for fft_n in FFT_SIZES:
        r = run_bench(fft_n, 0)
        if r:
            fft_costs[fft_n] = r["us_fft_frame"]
            print(f"  FFT {fft_n:4d}-pt: {r['us_fft_frame']:6.1f} µs/frame, "
                  f"{r['us_per_frame']:5.1f} µs combined")

    return bq_costs, fft_costs


def estimate_esp32_cpu(desktop_us_per_frame, desktop_cpu_pct,
                       esp32_factor=15.0):
    """Estimate ESP32 CPU load from desktop benchmark.

    240 MHz Xtensa LX6 is ~15× slower than modern desktop in practice.
    Uses the bench-reported desktop_cpu_pct as baseline.
    """
    return desktop_cpu_pct * esp32_factor


def combine_results(error_results, bq_costs, fft_costs, fs=48000.0):
    """Combine error and CPU data into a single table.

    Uses the actual bench-reported us_per_frame and frame_budget_us
    which correctly account for varying hop sizes (N/2 samples).
    """
    combined = []
    for r in error_results:
        fft_n = r["fft_n"]
        n_bq = r["n_biquads"]

        # CPU estimate: biquad cost + FFT cost
        us_bq = bq_costs.get(n_bq, 0.0)
        us_fft = fft_costs.get(fft_n, 0.0)
        us_total = us_bq + us_fft

        # Frame budget depends on hop = fft_n/2
        hop = fft_n // 2
        frame_budget_us = 1e6 / fs * hop

        if frame_budget_us > 0:
            desktop_cpu_pct = (us_total / frame_budget_us) * 100.0
            cpu_pct = estimate_esp32_cpu(us_total, desktop_cpu_pct)
            desktop_speedup = frame_budget_us / us_total if us_total > 0 else float('inf')
        else:
            cpu_pct = 0.0
            desktop_speedup = float('inf')

        combined.append({
            **r,
            "us_bq": us_bq,
            "us_fft": us_fft,
            "us_total": us_total,
            "cpu_esp32_pct": cpu_pct,
            "desktop_speedup": desktop_speedup,
            "frame_budget_us": frame_budget_us,
            "hop": hop,
        })

    return combined


def find_pareto_frontier(combined, x_key="cpu_esp32_pct", y_key="weighted_error"):
    """Return Pareto-optimal points (those not dominated)."""
    points = sorted(combined, key=lambda p: (p[x_key], p[y_key]))
    frontier = []
    best_y = float('inf')
    for p in points:
        if p[y_key] < best_y:
            frontier.append(p)
            best_y = p[y_key]
    return frontier


def print_report(combined, frontier):
    """Print a formatted sweep report."""
    print("\n" + "=" * 90)
    print("  FFT + IIR HYBRID SWEEP REPORT")
    print("  Minimizing perceptual error while staying under ~70% CPU on ESP32")
    print("=" * 90)

    # Sort by error ascending
    sorted_results = sorted(combined, key=lambda p: p["weighted_error"])

    print(f"\n{'FFT_N':>6s} {'BQs':>4s} {'RMSE':>6s} {'W.Err':>6s} "
          f"{'BQ µs':>7s} {'FFT µs':>7s} {'Sum µs':>7s} {'ESP32%':>7s} "
          f"{'DSpeed':>7s}  {'█ SAFE':7s}")
    print("-" * 85)

    for r in sorted_results:
        safe = "SAFE ✓" if r["cpu_esp32_pct"] < 70 else "RISK ✗"
        bar = "█" * max(1, int(min(r["cpu_esp32_pct"] / 5, 16)))
        print(f"{r['fft_n']:6d} {r['n_biquads']:4d} "
              f"{r['rmse_db']:5.2f}dB {r['weighted_error']:5.2f}dB "
              f"{r['us_bq']:6.1f}µs {r['us_fft']:6.1f}µs {r['us_total']:6.1f}µs "
              f"{r['cpu_esp32_pct']:5.1f}% {r['desktop_speedup']:5.0f}×  "
              f"{bar:7s} {safe}")

    print(f"\n── Pareto Frontier (best error for each CPU level) ──")
    for r in frontier:
        print(f"  FFT_N={r['fft_n']:4d}  {r['n_biquads']:2d} biquads  "
              f"err={r['weighted_error']:.2f} dB  "
              f"CPU={r['cpu_esp32_pct']:.1f}%  "
              f"({r['us_total']:.0f} µs / {r['frame_budget_us']:.0f} µs)")

    # Best safe config: among configurations with the best error,
    # pick the one with lowest CPU (FFT_N=256 is the sweet spot)
    safe = [r for r in sorted_results if r["cpu_esp32_pct"] < 70]
    if safe:
        best_error = safe[0]["weighted_error"]
        # Configs within 2% of best error
        contenders = [r for r in safe if r["weighted_error"] <= best_error * 1.02]
        best = min(contenders, key=lambda r: r["cpu_esp32_pct"])
        latency_ms = best["fft_n"] / 48000.0 * 1000.0
        print(f"\n── Best safe configuration ──")
        print(f"  FFT_N = {best['fft_n']}, {best['n_biquads']} biquads")
        print(f"  Resolution: {48000/best['fft_n']:.0f} Hz, hop = {best['hop']} samples")
        print(f"  Latency: ~{latency_ms:.1f} ms (window + overlap-add)")
        print(f"  Weighted error: {best['weighted_error']:.2f} dB  "
              f"RMSE: {best['rmse_db']:.2f} dB")
        print(f"  Est. ESP32 CPU: {best['cpu_esp32_pct']:.1f}%  "
              f"({best['us_total']:.0f} µs / {best['frame_budget_us']:.0f} µs)")
        print(f"  Desktop speedup: {best['desktop_speedup']:.0f}×")

        # Per-band breakdown
        print(f"\n  Per-band errors:")
        for label, lo, hi, _ in ERROR_BANDS:
            be = best["band_errors"].get(label, {})
            if be:
                print(f"    {label:<20s}  MAE={be['mae']:.2f} dB  "
                      f"max={be['max_err']:.2f} dB")
    else:
        print(f"\n── NO safe configuration found! All exceed 70% CPU. ──")


def main():
    ap = argparse.ArgumentParser(
        description="Sweep FFT bin count × IIR biquad count")
    ap.add_argument("--preset", default="lunchbox",
                    help="Preset name from presets/ directory")
    ap.add_argument("--measurements", nargs="+", default=None)
    ap.add_argument("--target", default=None)
    ap.add_argument("--noise", default=None)
    ap.add_argument("--fs", type=float, default=48000.0)
    ap.add_argument("--fc", type=float, default=60.0)
    ap.add_argument("--h2", type=float, default=0.5)
    ap.add_argument("--h3", type=float, default=1.0)
    ap.add_argument("--cpu-only", action="store_true",
                    help="Only run CPU benchmarks (skip error sweep)")
    ap.add_argument("--error-only", action="store_true",
                    help="Only compute errors (skip CPU benchmarks)")
    ap.add_argument("--esp32", action="store_true",
                    help="Also flash and test on ESP32 hardware")
    ap.add_argument("--json", default=None,
                    help="Save full results to JSON file")
    args = ap.parse_args()

    # ── Load measurements ──────────────────────────────────────
    if args.measurements and args.target:
        meas_paths = args.measurements
        target_path = args.target
        noise_path = args.noise
    else:
        pm = PresetManager()
        preset = pm.load(args.preset)
        meas_paths = [str(ROOT / p) for p in preset.measurements]
        target_path = str(ROOT / preset.target)
        noise_path = str(ROOT / preset.noise) if preset.noise else None
        args.fc = preset.fc if preset.fc is not None else args.fc
        args.h2 = preset.h2
        args.h3 = preset.h3
        print(f"  Preset: {args.preset}")

    # ── Run pipeline to get correction curve ────────────────────
    print(f"\n── Running EQ pipeline...")
    eval_freqs, target_db, measured_fs, max_gain_db = run_pipeline(
        meas_paths, target_path, noise_path,
        bass_enhancer_cutoff=args.fc, h2=args.h2, h3=args.h3,
    )
    design_fs = args.fs
    print(f"  {len(eval_freqs)} evaluation points, "
          f"{eval_freqs[0]:.0f}–{eval_freqs[-1]:.0f} Hz")
    print(f"  Max correction gain: {max_gain_db:+.1f} dB")

    # ── Sweep errors ───────────────────────────────────────────
    if not args.cpu_only:
        print(f"\n── Sweeping FFT_N × biquads for perceptual error...")
        error_results = sweep_errors(eval_freqs, target_db, design_fs)
        print(f"  Computed {len(error_results)} error data points")
    else:
        error_results = []

    # ── Sweep CPU ──────────────────────────────────────────────
    if not args.error_only:
        bq_costs, fft_costs = sweep_cpu()
    else:
        bq_costs, fft_costs = {}, {}

    # ── Combine & Report ───────────────────────────────────────
    if not args.cpu_only and not args.error_only:
        combined = combine_results(error_results, bq_costs, fft_costs, design_fs)
        frontier = find_pareto_frontier(combined)
        print_report(combined, frontier)

        if args.json:
            with open(args.json, "w") as f:
                json.dump({
                    "error_results": error_results,
                    "bq_costs": bq_costs,
                    "fft_costs": fft_costs,
                    "combined": combined,
                    "frontier": frontier,
                    "preset": args.preset,
                }, f, indent=2)
            print(f"\n  Saved full results to {args.json}")

    # ── Optional: ESP32 validation ─────────────────────────────
    if args.esp32:
        print(f"\n── ESP32 hardware validation ──")
        print(f"  (ESP32 profiling not yet implemented in this script)")
        print(f"  Build firmware with best config and flash manually:")
        best = combined[0] if combined else None
        if best:
            fft_n = best["fft_n"]
            n_bq = best["n_biquads"]
            print(f"\n  Recommended: FFT_N={fft_n}, {n_bq} biquads")
            print(f"  1. Run: FFT_N={fft_n} python -m eqgen.cli.wire setup "
                  f"--preset {args.preset} --max-bands {n_bq}")
            print(f"  2. Build firmware with FFT_EQ_N={fft_n}")
            print(f"  3. Flash and test for underruns")


if __name__ == "__main__":
    main()
