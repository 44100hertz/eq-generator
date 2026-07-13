#!/usr/bin/env python3
"""Test the IIR fitter with real measurement data — enhancer region focus.
   Uses the pipeline to compute a correction curve from real measurements,
   then fits IIR biquads and reports error by frequency band."""
import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import run_pipeline
from eqgen.eq_fit import fit_eq_curve, cascade_response_db, BiquadCoeffs
from eqgen.quantize import BiquadQ28, quantize_biquads_q28, q28_to_float

mp = str(ROOT/"measurements/technics/standing/measurement2.wav")
tp = str(ROOT/"measurements/technics/standing/target.wav")

FC = 60.0; H2, H3 = 0.5, 1.0

print("="*70)
print(f"  technics/standing — IIR fit from measurement pipeline")
print("="*70)

# Run pipeline: measurement → correction curve (no model preprocessing)
freqs, target_db, fs, _max_gain, _efficacy = run_pipeline([mp], tp)

nyquist = min(16000.0, fs/2 * 0.95)

# Clamp to ±24 dB (peaking filter range)
target_db = np.clip(target_db, -24.0, 24.0)

print(f"\n  Correction curve (clamped to ±24 dB):")
for fc in [20, 30, 40, 50, 60, 80, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 20000]:
    idx = np.argmin(np.abs(freqs - fc))
    print(f"    {freqs[idx]:8.0f}  {target_db[idx]:+7.1f} dB")

# Fit IIR biquads
fit = fit_eq_curve(
    freqs, target_db, fs, max_bands=40,
    min_freq=freqs[0], max_freq=freqs[-1], min_peaking_freq=max(40.0, freqs[0]),
    gain_range=(-60.0, 60.0), q_range=(0.3, 6.0), stop_db=0.3,
)

bq_q28 = quantize_biquads_q28(fit.biquads)
q28_floats = [BiquadCoeffs(b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                           b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                           a2=q28_to_float(b.a2)) for b in bq_q28]
fitted_db = cascade_response_db(q28_floats, freqs, fs)
err = fitted_db - target_db

bands_def = [
    ("Enhancer region (40-120 Hz)", 40, 120),
    ("Upper bass (120-250 Hz)",      120, 250),
    ("Low mid (250-500 Hz)",         250, 500),
    ("Mid (500-2k Hz)",              500, 2000),
    ("Upper mid (2k-8k Hz)",         2000, 8000),
    ("Treble (8k-20k Hz)",           8000, 22050),
    ("FULL BAND (20-20k Hz)",         20, 22050),
]
print(f"\n  ═══ Real measurement fit — {len(fit.bands)} bands ═══")
print(f"  {'Band':>28s}  {'Max err':>8s}  {'RMS err':>8s}")
print(f"  {'':->28s}  {'':->8s}  {'':->8s}")
for bl, lo, hi in bands_def:
    m = (freqs >= lo) & (freqs <= hi)
    if not m.any(): continue
    print(f"  {bl:>28s}:  {np.max(np.abs(err[m])):+7.2f} dB  "
          f"{np.sqrt(np.mean(err[m]**2)):+7.2f} dB")

# Spot check
print(f"\n  Spot check:")
for fc in [20, 40, 60, 80, 100, 150, 200, 500, 1000, 2000, 5000, 10000, 20000]:
    idx = np.argmin(np.abs(freqs - fc))
    print(f"    {freqs[idx]:5.0f} Hz: target={target_db[idx]:+6.1f}  "
          f"fit={fitted_db[idx]:+6.1f}  err={err[idx]:+6.2f} dB")

enhancer_bands = [b for b in fit.bands if 40 <= b['f0'] <= 200]
print(f"\n  Bands in enhancer region (40-200 Hz): {len(enhancer_bands)}/{len(fit.bands)}")
print(f"  Bands above 200 Hz: {len(fit.bands)-len(enhancer_bands)}")
