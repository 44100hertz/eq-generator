#!/usr/bin/env python3
"""Quick test: 40 peaking filters, full-band, how bad is the error?"""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "python"))

from eq_fit import fit_eq_curve, cascade_response_db, BiquadCoeffs
from quantize import BiquadQ28, quantize_biquads_q28, q28_to_float

# ── Synthetic speaker: -12 dB @ 50 Hz, flat above 100 Hz ──────────────
def small_speaker(f):
    if f <= 50:   return 0.25
    if f >= 100:  return 1.0
    return 0.25 + 0.75 * (f - 50) / 50

fs = 44100.0
freqs = np.logspace(np.log10(20), np.log10(20000), 300)

# Target curve: inverse speaker (flat perceived output)
speaker_linear = np.array([small_speaker(f) for f in freqs])
target_db = 20.0 * np.log10(1.0 / np.maximum(speaker_linear, 1e-12))

print("=" * 70)
print("  40 PEAKING FILTERS — FULL-BAND FIT (20 Hz – 20 kHz)")
print("  Target: flat output → inverse of speaker (-12 dB @ 50 Hz)")
print("=" * 70)

# ── Fit ───────────────────────────────────────────────────────────────
fit = fit_eq_curve(
    freqs, target_db, fs, max_bands=40,
    min_freq=20.0, max_freq=20000.0, min_peaking_freq=40.0,
)

# ── Q4.28 quantization ────────────────────────────────────────────────
bq_q28 = quantize_biquads_q28(fit.biquads)
q28_floats = [
    BiquadCoeffs(
        b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
        b2=q28_to_float(b.b2), a1=q28_to_float(b.a1), a2=q28_to_float(b.a2))
    for b in bq_q28
]
fitted_db = cascade_response_db(q28_floats, freqs, fs)

# Float-precision (before quantization) for comparison
fitted_float_db = cascade_response_db(fit.biquads, freqs, fs)

# ── Error analysis ────────────────────────────────────────────────────
err_float = fitted_float_db - target_db
err_q28   = fitted_db - target_db

def band_stats(label, lo, hi):
    m = (freqs >= lo) & (freqs <= hi)
    if not m.any():
        return (0, 0, 0, 0)
    max_abs = np.max(np.abs(err_q28[m]))
    rms = np.sqrt(np.mean(err_q28[m] ** 2))
    max_float = np.max(np.abs(err_float[m]))
    return max_abs, rms, max_float

bands = [
    ("Sub-bass  (20-40 Hz)",    20, 40),
    ("Bass      (40-100 Hz)",   40, 100),
    ("Upper bass(100-250 Hz)",  100, 250),
    ("Low mid   (250-500 Hz)",  250, 500),
    ("Mid       (500-2k Hz)",   500, 2000),
    ("Upper mid (2k-8k Hz)",    2000, 8000),
    ("Treble    (8k-20k Hz)",   8000, 20000),
]

print(f"\n{'Octave band':>25s}  {'Max err (float)':>12s}  {'Max err (Q4.28)':>12s}  {'RMS err (Q4.28)':>12s}")
print(f"{'':->25s}  {'':->12s}  {'':->12s}  {'':->12s}")
for label, lo, hi in bands:
    max_abs, rms, max_f = band_stats(label, lo, hi)
    print(f"  {label:>23s}:  {max_f:+11.3f} dB  {max_abs:+11.3f} dB  {rms:+11.3f} dB")

# ── Overall ───────────────────────────────────────────────────────────
overall_max  = np.max(np.abs(err_q28))
overall_rms  = np.sqrt(np.mean(err_q28 ** 2))
overall_maxf = np.max(np.abs(err_float))
print(f"\n  {'OVERALL':>23s}:  {overall_maxf:+11.3f} dB  {overall_max:+11.3f} dB  {overall_rms:+11.3f} dB")

# ── Spot checks ───────────────────────────────────────────────────────
print(f"\n  Spot check (freq → target → fitted → error):")
for f_check in [20, 40, 60, 100, 250, 500, 1000, 2000, 5000, 10000, 20000]:
    idx = np.argmin(np.abs(freqs - f_check))
    print(f"    {freqs[idx]:5.0f} Hz: target={target_db[idx]:+6.1f}  "
          f"fit={fitted_db[idx]:+6.1f}  err={err_q28[idx]:+6.2f} dB")

# ── Band distribution ─────────────────────────────────────────────────
print(f"\n  Band distribution ({len(fit.bands)} total):")
for i, b in enumerate(fit.bands):
    print(f"    [{i:2d}] f0={b['f0']:6.1f} Hz  gain={b['gain_db']:+5.1f} dB  Q={b['Q']:.2f}")

# ── Quantization impact stats ─────────────────────────────────────────
q_err = np.abs(fitted_db - fitted_float_db)
q_max = np.max(q_err)
q_bins = [0.01, 0.05, 0.10, 0.50, 1.00]
print(f"\n  Q4.28 quantization impact: max {q_max:.3f} dB")
for thresh in q_bins:
    frac = np.sum(q_err <= thresh) / len(q_err) * 100
    print(f"    ≤ {thresh:.2f} dB: {frac:.0f}% of freq bins")
