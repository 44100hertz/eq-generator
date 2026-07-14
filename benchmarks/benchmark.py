#!/usr/bin/env python3
"""Benchmark harness for IIR biquad fitters.

Measures perceptual (ISO 226-weighted) error across all presets.
Run:  python benchmarks/benchmark.py [--preset X] [--bands N]
"""

import sys
import time
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import run_pipeline
from eqgen.presets import PresetManager, Preset, MAX_IIR_BANDS, get_house_curve
from eqgen.eq_fit import FitResult, BiquadCoeffs, cascade_response_db
from eqgen.model import ear_sensitivity

# ── Psychoacoustic bands for reporting ────────────────────────────────
PERCEPTUAL_BANDS = [
    ("Sub-bass (20-40)",      20,   40),
    ("Enhancer (40-120)",     40,  120),
    ("Upper bass (120-250)",  120,  250),
    ("Low mid (250-500)",     250,  500),
    ("Mid (500-2k)",         500, 2000),
    ("Upper mid (2k-8k)",   2000, 8000),
    ("Treble (8k-14k)",     8000, 14000),
    ("FULL (20-14k)",          20, 14000),
]


def perceptual_rms_error(freqs: np.ndarray, target_db: np.ndarray,
                          fitted_db: np.ndarray) -> float:
    """Perceptually-weighted RMS error in dB.

    Weights each frequency bin by ear_sensitivity(f) so errors in regions
    the ear is most sensitive to (1-5 kHz) count more than bass/treble.
    Returns RMS error in dB.
    """
    weights = np.array([ear_sensitivity(f) for f in freqs])
    err = fitted_db - target_db
    weighted_err = err * weights
    return float(np.sqrt(np.mean(weighted_err ** 2)))


def unweighted_rms_error(freqs: np.ndarray, target_db: np.ndarray,
                          fitted_db: np.ndarray) -> float:
    """Standard unweighted RMS error in dB."""
    err = fitted_db - target_db
    return float(np.sqrt(np.mean(err ** 2)))


def max_abs_error(freqs: np.ndarray, target_db: np.ndarray,
                   fitted_db: np.ndarray) -> float:
    """Maximum absolute error in dB."""
    return float(np.max(np.abs(fitted_db - target_db)))


# ── Fitter protocol ───────────────────────────────────────────────────
# Each fitter module must expose:
#   fit_bands(freqs, target_db, fs, max_bands) -> (list_of_BiquadCoeffs, metadata_dict)
# where metadata_dict has at least {"name": "fitter-name"}

DEFAULT_FITTER = "eqgen.eq_fit"


def load_presets() -> List[Preset]:
    """Load all available presets."""
    pm = PresetManager()
    names = pm.list_presets()
    presets = []
    for name in names:
        try:
            p = pm.load(name)
            presets.append(p)
        except Exception as e:
            print(f"  SKIP {name}: {e}", file=sys.stderr)
    return presets


def run_benchmark(fitter_module: str, presets: Optional[List[str]] = None,
                   band_counts: Optional[List[int]] = None,
                   ) -> dict:
    """Run the benchmark across all (or selected) presets.

    Returns dict: {preset_name: {n_bands: {perceptual_rms, unweighted_rms, max_err, time_s, n_bands_used}}}
    """
    import importlib
    mod = importlib.import_module(fitter_module)

    all_presets = load_presets()
    if presets:
        all_presets = [p for p in all_presets if p.name in presets]

    if not band_counts:
        band_counts = [8, 16, 24]

    results = {}

    for preset in all_presets:
        print(f"\n  ── {preset.name} ──")
        sys.stdout.flush()

        try:
            meas_paths = preset.resolve_measurements()
            targ_path = preset.resolve_target()
            noise_path = preset.resolve_noise()

            hc = get_house_curve(preset.house_curve) if preset.house_curve else None

            freqs, gains_db, fs, _max_gain, _efficacy = run_pipeline(
                meas_paths, targ_path, noise_path,
                bass_enhancer_cutoff=preset.fc,
                h2=preset.h2, h3=preset.h3,
                smooth_exponent=preset.smooth_exponent,
                house_curve=hc,
            )

            # Clamp gains to safe range
            target_db = np.clip(gains_db, -24.0, 24.0)

        except Exception as e:
            print(f"    SKIP: pipeline failed: {e}", file=sys.stderr)
            continue

        preset_results = {}
        for n_bands in band_counts:
            if n_bands > MAX_IIR_BANDS:
                continue

            t0 = time.perf_counter()
            biquads, meta = mod.fit_bands(freqs, target_db, fs, max_bands=n_bands)
            elapsed = time.perf_counter() - t0

            q28_floats = biquads
            fitted_db = cascade_response_db(q28_floats, freqs, fs)

            p_rms = perceptual_rms_error(freqs, target_db, fitted_db)
            u_rms = unweighted_rms_error(freqs, target_db, fitted_db)
            mx = max_abs_error(freqs, target_db, fitted_db)

            # Per-band breakdown
            band_errors = {}
            for bl, lo, hi in PERCEPTUAL_BANDS:
                m = (freqs >= lo) & (freqs <= hi)
                if m.any():
                    w = np.array([ear_sensitivity(f) for f in freqs[m]])
                    err = fitted_db[m] - target_db[m]
                    band_errors[bl] = float(np.sqrt(np.mean((err * w) ** 2)))

            preset_results[n_bands] = {
                "perceptual_rms": p_rms,
                "unweighted_rms": u_rms,
                "max_err": mx,
                "time_s": elapsed,
                "n_bands_used": len(biquads),
                "band_errors": band_errors,
            }

            print(f"    {n_bands:2d} bands → pRMS={p_rms:.3f}  uRMS={u_rms:.3f}  "
                  f"max={mx:.3f}  ({elapsed:.2f}s, used {len(biquads)})")

        results[preset.name] = preset_results

    return results


def summarize(results: dict):
    """Print summary table comparing fitters."""
    print("\n" + "=" * 90)
    print(f"{' SUMMARY ':~^90}")
    print("=" * 90)

    for preset_name, preset_data in sorted(results.items()):
        print(f"\n  {preset_name}:")
        print(f"    {'Bands':>5s}  {'pRMS':>8s}  {'uRMS':>8s}  {'Max':>8s}  {'Time':>8s}  {'Used':>5s}")
        print(f"    {'':->5s}  {'':->8s}  {'':->8s}  {'':->8s}  {'':->8s}  {'':->5s}")
        for n_b, r in sorted(preset_data.items()):
            print(f"    {n_b:5d}  {r['perceptual_rms']:8.3f}  {r['unweighted_rms']:8.3f}  "
                  f"{r['max_err']:8.3f}  {r['time_s']:8.2f}  {r['n_bands_used']:5d}")

    # Aggregate: mean perceptual RMS across presets at each band count
    band_counts = sorted(set(n for pd in results.values() for n in pd))
    print("\n" + "=" * 90)
    print("  AGGREGATE (mean perceptual RMS across all presets)")
    print("=" * 90)
    print(f"  {'Bands':>5s}  {'Mean pRMS':>12s}  {'Mean uRMS':>12s}  {'Mean Time':>10s}")
    print(f"  {'':->5s}  {'':->12s}  {'':->12s}  {'':->10s}")
    for n_b in band_counts:
        vals_p = []
        vals_u = []
        vals_t = []
        for pd in results.values():
            if n_b in pd:
                vals_p.append(pd[n_b]['perceptual_rms'])
                vals_u.append(pd[n_b]['unweighted_rms'])
                vals_t.append(pd[n_b]['time_s'])
        if vals_p:
            print(f"  {n_b:5d}  {np.mean(vals_p):12.3f}  {np.mean(vals_u):12.3f}  "
                  f"{np.mean(vals_t):10.2f}")


# ── CLI ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Benchmark IIR biquad fitters")
    parser.add_argument("--preset", nargs="*", default=None,
                        help="Preset(s) to benchmark (default: all)")
    parser.add_argument("--fitter", default=DEFAULT_FITTER,
                        help=f"Fitter module (default: {DEFAULT_FITTER})")
    parser.add_argument("--bands", nargs="*", type=int, default=[8, 16, 24],
                        help="Band counts to test")
    args = parser.parse_args()

    results = run_benchmark(args.fitter, presets=args.preset, band_counts=args.bands)
    summarize(results)
