#!/usr/bin/env python3
"""Sine sweep sanity test: C enhancer → measure output across frequencies.
   Checks whether the output level is approximately flat across frequency."""
import sys, struct
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.eq_fit import fit_eq_curve, cascade_response_db, BiquadCoeffs
from eqgen.analysis import goertzel_magnitude
from eqgen.io import load_measurement
from eqgen.pipeline import run_pipeline
from eqgen import enhancer_ffi as effi

# ── Design EQ from measurement data ────────────────────────────────────
FC = 60.0; H2, H3 = 0.5, 1.0; FS = 44100.0; MAX_BANDS = 40

mp = str(ROOT/"measurements/technics/standing/measurement2.wav")
tp = str(ROOT/"measurements/technics/standing/target.wav")

# Run pipeline: measurement → correction curve
freqs, gains_db, fs, _max_gain, _efficacy = run_pipeline([mp], tp)

# Fit IIR biquads to the correction curve
fit = fit_eq_curve(freqs, gains_db, FS, max_bands=MAX_BANDS,
                   min_freq=freqs[0], max_freq=freqs[-1],
                   min_peaking_freq=40.0,
                   gain_range=(-60.0, 60.0), q_range=(0.3, 6.0))

coeffs = [v for bc in fit.biquads for v in [bc.b0, bc.b1, bc.b2, bc.a1, bc.a2]]

print("="*70)
print(f"  SINE SWEEP — C enhancer, technics/standing, fc={FC} h2={H2} h3={H3}")
print(f"  {len(fit.bands)} bands")
print("="*70)

# ── Build enhancer ────────────────────────────────────────────────────
enh = effi.create_enhancer(
    cutoff_hz=FC, h2_amp=H2, h3_amp=H3,
    release_secs=0.2, fs=FS, coeffs=coeffs,
)

# ── Sweep ─────────────────────────────────────────────────────────────
AMP_IN = 0.001  # -60 dBFS — very low to survive EQ boost without clipping
DURATION = 1.0
N_SAMPLES = int(DURATION * FS)
SS_START = N_SAMPLES // 4  # skip startup transient

test_freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120,
               150, 200, 300, 500, 1000, 2000]

print(f"\n  Input: {AMP_IN:.1f} ({20*np.log10(AMP_IN):+.0f} dBFS), "
      f"{DURATION}s per tone, steady-state from {SS_START/FS:.0f}ms")
print(f"\n  {'Freq':>6s}  {'RMS':>8s}  {'dBFS':>7s}  "
      f"{'Fund':>7s}  {'H2':>7s}  {'H3':>7s}")
print(f"  {'':->6s}  {'':->8s}  {'':->7s}  "
      f"{'':->7s}  {'':->7s}  {'':->7s}")

results = []
for f_test in test_freqs:
    t = np.arange(N_SAMPLES) / FS
    sine = AMP_IN * np.sin(2.0 * np.pi * f_test * t)

    # → int16 stereo
    pcm = bytearray(N_SAMPLES * 4)
    for i in range(N_SAMPLES):
        v = int(np.clip(sine[i] * 32767, -32768, 32767))
        struct.pack_into('<hh', pcm, i*4, v, v)

    effi.reset_enhancer(enh)
    for i in range(0, len(pcm), 4):
        l = struct.unpack_from('<h', pcm, i)[0]
        r = struct.unpack_from('<h', pcm, i+2)[0]
        l_out, r_out = effi.process_stereo_frame(enh, l, r)
        struct.pack_into('<hh', pcm, i, int(l_out), int(r_out))

    # Decode steady-state
    out_float = np.array([struct.unpack_from('<h', pcm, i*4)[0]/32768.0
                          for i in range(SS_START, N_SAMPLES)])

    rms_out = np.sqrt(np.mean(out_float**2))
    fund_amp = goertzel_magnitude(out_float, f_test, FS)
    h2_amp = goertzel_magnitude(out_float, 2*f_test, FS) if 2*f_test < FS/2 else 0
    h3_amp = goertzel_magnitude(out_float, 3*f_test, FS) if 3*f_test < FS/2 else 0

    results.append((f_test, rms_out, fund_amp, h2_amp, h3_amp))

    print(f"  {f_test:6.0f}  {rms_out:8.4f}  {20*np.log10(max(rms_out,1e-6)):+6.1f}  "
          f"{fund_amp:7.4f}  {h2_amp:7.4f}  {h3_amp:7.4f}")

effi.destroy_enhancer(enh)

# ── Flatness summary ───────────────────────────────────────────────────
print(f"\n\n  ── Flatness check ──")
rmss = np.array([r[1] for r in results])
ref = np.mean(rmss)
print(f"  Mean output RMS: {ref:.4f} ({20*np.log10(ref):+.1f} dBFS)")
spread = 20*np.log10(np.max(rmss)/max(np.min(rmss), 1e-6))
print(f"  Min/Max spread:  {spread:.1f} dB")

bass_rms  = [r[1] for r in results if r[0] <= 120]
treble_rms = [r[1] for r in results if r[0] > 120]
if bass_rms:
    print(f"  Bass (20-120 Hz): {20*np.log10(max(bass_rms)/max(min(bass_rms), 1e-6)):.1f} dB spread")
if treble_rms:
    print(f"  Treble (>120 Hz): {20*np.log10(max(treble_rms)/max(min(treble_rms), 1e-6)):.1f} dB spread")
