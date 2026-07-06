#!/usr/bin/env python3
"""Test optimized fitter with bass enhancer model — focus on what matters.
   Clamp model target to ±24 dB (peaking filter range) and report
   error in the critical enhancer region vs full band."""
import sys, struct
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.model import preprocess_eq_curve
from eqgen.eq_fit import fit_eq_curve, cascade_response_db, BiquadCoeffs
from eqgen.quantize import BiquadQ28, quantize_biquads_q28, q28_to_float
from eqgen.io import load_measurement


mp = ROOT/"measurements/technics/standing/measurement2.wav"
tp = ROOT/"measurements/technics/standing/target.wav"
speaker_fn, fs = load_measurement(mp, tp)
nyquist = min(16000.0, fs/2 * 0.95)

FC = 60.0; H2, H3 = 0.5, 1.0

print("="*70)
print(f"  technics/standing — enhancer ON: fc={FC} h2={H2} h3={H3}")
print("="*70)

freqs = np.logspace(np.log10(FC/2.0), np.log10(nyquist), 300)
target_curve = {float(f): 0.5 for f in freqs}

# ── Model: raw (no clamp) and clamped to ±24 dB ───────────────────────
eq_curve = preprocess_eq_curve(target_curve, speaker_fn, FC, H2, H3)
eq_freqs = np.array(sorted(eq_curve.keys()))
target_raw_db = 20.0 * np.log10(np.array([max(eq_curve[f], 1e-12) for f in eq_freqs]))
target_db = np.clip(target_raw_db, -24.0, 24.0)

# Show the model target: raw vs clamped
print(f"\n  Model target curve (raw → clamped to ±24 dB):")
print(f"  {'Freq':>8s}  {'Raw':>8s}  {'Clamped':>8s}  {'Lost':>8s}")
for fc in [20, 30, 40, 50, 60, 80, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 20000]:
    idx = np.argmin(np.abs(eq_freqs - fc))
    raw = target_raw_db[idx]
    clp = target_db[idx]
    lost = raw - clp
    print(f"  {eq_freqs[idx]:8.0f}  {raw:+7.1f} dB  {clp:+7.1f} dB  {lost:+7.1f} dB")

# ── Fit both clamped and raw ──────────────────────────────────────────
for label, tgt in [("CLAMPED ±24 dB", target_db), ("RAW (unclamped)", target_raw_db)]:
    fit = fit_eq_curve(
        eq_freqs, tgt, fs, max_bands=40,
        min_freq=FC/2.0, max_freq=nyquist, min_peaking_freq=max(40.0, FC/2.0),
        gain_range=(-60.0, 60.0), q_range=(0.3, 6.0), stop_db=0.3,
    )
    bq_q28 = quantize_biquads_q28(fit.biquads)
    q28_floats = [BiquadCoeffs(b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                               b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                               a2=q28_to_float(b.a2)) for b in bq_q28]
    fitted_db = cascade_response_db(q28_floats, eq_freqs, fs)
    err = fitted_db - tgt  # error vs whichever target we used

    bands_def = [
        ("Enhancer region (40-120 Hz)", 40, 120),
        ("Upper bass (120-250 Hz)",      120, 250),
        ("Low mid (250-500 Hz)",         250, 500),
        ("Mid (500-2k Hz)",              500, 2000),
        ("Upper mid (2k-8k Hz)",         2000, 8000),
        ("Treble (8k-20k Hz)",           8000, 22050),
        ("FULL BAND (20-20k Hz)",         20, 22050),
    ]
    print(f"\n  ═══ {label} — {len(fit.bands)} bands ═══")
    print(f"  {'Band':>28s}  {'Max err':>8s}  {'RMS err':>8s}")
    print(f"  {'':->28s}  {'':->8s}  {'':->8s}")
    for bl, lo, hi in bands_def:
        m = (eq_freqs >= lo) & (eq_freqs <= hi)
        if not m.any(): continue
        print(f"  {bl:>28s}:  {np.max(np.abs(err[m])):+7.2f} dB  "
              f"{np.sqrt(np.mean(err[m]**2)):+7.2f} dB")

    # Spot check
    print(f"\n  Spot check:")
    for fc in [20, 40, 60, 80, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 20000]:
        idx = np.argmin(np.abs(eq_freqs - fc))
        print(f"    {eq_freqs[idx]:5.0f} Hz: target={tgt[idx]:+6.1f}  "
              f"fit={fitted_db[idx]:+6.1f}  err={err[idx]:+6.2f} dB")

    # How many bands in the critical enhancer region?
    enhancer_bands = [b for b in fit.bands if 40 <= b['f0'] <= 200]
    print(f"\n  Bands in enhancer region (40-200 Hz): {len(enhancer_bands)}/{len(fit.bands)}")
    print(f"  Bands above 200 Hz: {len(fit.bands)-len(enhancer_bands)}")
