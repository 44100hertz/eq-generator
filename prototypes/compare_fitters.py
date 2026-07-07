#!/usr/bin/env python3
"""Profile the IIR fitter against the actual JamesDSP EQ curve from eqgen.py.

Uses eqgen.run_pipeline() — the same Welch + adaptive-point + bass-enhancer
pipeline that produces the CSV — as the IIR fitter's target.  No longer
re-derives a separate target from raw FFT data.
"""
import sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import run_pipeline
from eqgen.eq_fit import fit_eq_curve, cascade_response_db, BiquadCoeffs
from eqgen.quantize import BiquadQ28, quantize_biquads_q28, q28_to_float

# ── Psychoacoustic bands ───────────────────────────────────────────────
BANDS = [
    ("Sub-bass",             20,   40),
    ("Bass (enhancer)",      40,  120),
    ("Upper bass",           120,  250),
    ("Low mid",              250,  500),
    ("Mid (critical)",       500, 2000),
    ("Upper mid",           2000, 8000),
    ("Treble",              8000, 16000),
    ("FULL",                  20, 16000),
]


def spot_freqs():
    return [40, 60, 80, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 16000]


# ── Run profile ────────────────────────────────────────────────────────


def profile(freqs, target_db, fs, label, band_counts):
    f_min = freqs[0]
    f_max = freqs[-1]

    print(f"\n{'='*80}")
    print(f"  {label}")
    print(f"  {len(freqs)} adaptive EQ points  {f_min:.0f}–{f_max:.0f} Hz  fs={fs:.0f}")
    print(f"{'='*80}")

    # Header
    hdr = f"  {'Bands':>5s}  {'Time':>6s}" + "".join(
        f"  {bl:>8s}" for bl, _, _ in BANDS
    )
    sep = f"  {'':->5s}  {'':->6s}" + "".join(f"  {'':->8s}" for _ in BANDS)
    print(hdr)
    print(sep)

    for n in band_counts:
        t0 = time.perf_counter()
        fit = fit_eq_curve(
            freqs, target_db, fs,
            max_bands=n, min_freq=f_min, max_freq=f_max,
            min_peaking_freq=f_min,
        )
        elapsed = time.perf_counter() - t0

        # Quantized
        bq_q28 = quantize_biquads_q28(fit.biquads)
        q28f = [
            BiquadCoeffs(
                b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                a2=q28_to_float(b.a2),
            )
            for b in bq_q28
        ]
        fitted_db = cascade_response_db(q28f, freqs, fs)
        float_db = cascade_response_db(fit.biquads, freqs, fs)
        err = fitted_db - target_db
        q_delta = np.max(np.abs(fitted_db - float_db))

        parts = []
        for bl, lo, hi in BANDS:
            m = (freqs >= lo) & (freqs <= hi)
            parts.append(
                f"  {np.max(np.abs(err[m])):7.2f}dB" if m.any() else f"  {'—':>7s} "
            )

        print(
            f"  {len(fit.bands):5d}  {elapsed:5.1f}s"
            + "".join(parts)
            + f"            {q_delta:7.3f}dB"
        )

    # Detailed breakdown at highest band count
    n_max = band_counts[-1]
    fit = fit_eq_curve(
        freqs, target_db, fs,
        max_bands=n_max, min_freq=f_min, max_freq=f_max,
        min_peaking_freq=f_min,
    )
    bq_q28 = quantize_biquads_q28(fit.biquads)
    q28f = [
        BiquadCoeffs(
            b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
            b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
            a2=q28_to_float(b.a2),
        )
        for b in bq_q28
    ]
    fitted_db = cascade_response_db(q28f, freqs, fs)
    err = fitted_db - target_db

    print(f"\n  ── RMS error at {n_max} bands (Q4.28) ──")
    for bl, lo, hi in BANDS:
        m = (freqs >= lo) & (freqs <= hi)
        if m.any():
            rms = np.sqrt(np.mean(err[m] ** 2))
            mx = np.max(np.abs(err[m]))
            print(f"    {bl:<20s}  RMS {rms:5.1f} dB   max {mx:5.1f} dB")

    print(f"\n  ── Spot check (target → fitted → error) ──")
    print(f"    {'Freq':>6s}  {'Target':>8s}  {'Fitted':>8s}  {'Error':>8s}")
    for f_chk in spot_freqs():
        idx = np.argmin(np.abs(freqs - f_chk))
        print(
            f"    {freqs[idx]:6.0f}  {target_db[idx]:+7.1f} dB"
            f"  {fitted_db[idx]:+7.1f} dB  {err[idx]:+7.2f} dB"
        )

    bass_n = sum(1 for b in fit.bands if b["f0"] <= 250)
    mid_n  = sum(1 for b in fit.bands if 250 < b["f0"] <= 2000)
    high_n = sum(1 for b in fit.bands if b["f0"] > 2000)
    print(f"\n  ── {len(fit.bands)} bands ──")
    print(f"    Bass (≤250 Hz):  {bass_n}")
    print(f"    Mid  (250-2k):   {mid_n}")
    print(f"    High (>2 kHz):   {high_n}")

    print(f"\n  ── Target stats ──")
    print(f"    Range:  {np.min(target_db):+.1f} .. {np.max(target_db):+.1f} dB")
    print(f"    Mean:   {np.mean(target_db):+.1f} dB")
    n_clip = np.sum(np.abs(target_db) > 24.0)
    print(f"    Bins > ±24 dB:  {n_clip}/{len(target_db)}")


# ── Main ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    COUNTS = [4, 8, 16, 32, 48]

    scenarios = [
        # (label,                   measurements,              target,                     noise,                    fc,  h2,  h3)
        ("technics — NO enhancer",
         ["measurements/technics/standing/measurement2.wav"],
         "measurements/technics/standing/target.wav",
         None,
         None, 0.0, 0.0),

        ("technics — enhancer fc=60 h2=0.5 h3=1.0",
         ["measurements/technics/standing/measurement2.wav"],
         "measurements/technics/standing/target.wav",
         None,
         60.0, 0.5, 1.0),

        ("cardboard — enhancer fc=50 h2=0.5 h3=1.0",
         ["measurements/cardboard/response.wav"],
         "measurements/cardboard/target.wav",
         "measurements/cardboard/noise.wav",
         50.0, 0.5, 1.0),
    ]

    for label, meas, targ, noise, fc, h2, h3 in scenarios:
        try:
            freqs, gains_db, fs, _max_gain = run_pipeline(
                meas, targ, noise,
                max_noise=0.65,
                bass_enhancer_cutoff=fc,
                h2=h2, h3=h3,
            )
        except Exception as e:
            print(f"  SKIP {label}: {e}")
            continue
        profile(freqs, gains_db, fs, label, COUNTS)
