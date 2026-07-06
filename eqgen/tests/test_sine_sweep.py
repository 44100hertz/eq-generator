#!/usr/bin/env python3
"""Sine sweep sanity test: EQ → C enhancer → measure perceived output.
   Checks whether the output level is approximately flat across frequency."""
import sys, struct
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.model import preprocess_eq_curve
from eqgen.eq_fit import fit_eq_curve
from eqgen.quantize import BiquadQ28, quantize_biquads_q28
from eqgen.analysis import goertzel_magnitude
from eqgen.io import load_measurement
from eqgen import enhancer_ffi as effi

# ── Design EQ ──────────────────────────────────────────────────────────
FC = 60.0; H2, H3 = 0.5, 1.0; FS = 44100.0; MAX_BANDS = 40

mp = ROOT/"measurements/technics/standing/measurement2.wav"
tp = ROOT/"measurements/technics/standing/target.wav"
speaker_fn, _ = load_measurement(mp, tp)

f_min = FC/2.0  # below fc/2 the model ramp handles sub-bass
f_max = min(16000.0, FS/2*0.95)
freqs = np.logspace(np.log10(f_min), np.log10(f_max), 200)

target_curve = {float(f): 0.5 for f in freqs}
eq_curve = preprocess_eq_curve(target_curve, speaker_fn, FC, H2, H3)
eq_freqs = np.array(sorted(eq_curve.keys()))
target_db = 20.0 * np.log10(np.array([max(eq_curve[f], 1e-12) for f in eq_freqs]))

fit = fit_eq_curve(eq_freqs, target_db, FS, max_bands=MAX_BANDS,
                             min_freq=f_min, max_freq=f_max,
                             min_peaking_freq=40.0,
                             gain_range=(-24.0, 24.0), q_range=(0.3, 6.0))

bq_q28 = quantize_biquads_q28(fit.biquads)
coeffs = [v for bq in bq_q28 for v in [bq.b0, bq.b1, bq.b2, bq.a1, bq.a2]]

print("="*70)
print(f"  SINE SWEEP — technics/standing, fc={FC} h2={H2} h3={H3}")
print(f"  {len(fit.bands)} bands, optimized fitter")
print("="*70)

# ── Build enhancer ────────────────────────────────────────────────────
enh = effi.create_enhancer(
    cutoff_hz=FC, h2_amp=H2, h3_amp=H3,
    release_secs=0.2, fs=FS, coeffs_q28=coeffs,
)

# ── Sweep ─────────────────────────────────────────────────────────────
AMP_IN = 0.001  # -60 dBFS — very low to survive EQ boost + enhancer without clipping
DURATION = 1.0
N_SAMPLES = int(DURATION * FS)
SS_START = N_SAMPLES // 4  # skip startup transient

test_freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120, 150, 200, 300, 500, 1000, 2000]

print(f"\n  Input: {AMP_IN:.1f} ({20*np.log10(AMP_IN):+.0f} dBFS), "
      f"{DURATION}s per tone, steady-state from {SS_START/FS:.0f}ms")
print(f"\n  {'Freq':>6s}  {'RMS':>8s}  {'dBFS':>7s}  "
      f"{'Fund':>7s}  {'H2':>7s}  {'H3':>7s}  "
      f"{'Model EQ':>9s}  {'Δ flat':>7s}")
print(f"  {'':->6s}  {'':->8s}  {'':->7s}  "
      f"{'':->7s}  {'':->7s}  {'':->7s}  "
      f"{'':->9s}  {'':->7s}")

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
        struct.pack_into('<hh', pcm, i, l_out, r_out)

    # Decode steady-state
    out_float = np.array([struct.unpack_from('<h', pcm, i*4)[0]/32768.0
                          for i in range(SS_START, N_SAMPLES)])

    rms_out = np.sqrt(np.mean(out_float**2))
    fund_amp = goertzel_magnitude(out_float, f_test, FS)
    h2_amp = goertzel_magnitude(out_float, 2*f_test, FS) if 2*f_test < FS/2 else 0
    h3_amp = goertzel_magnitude(out_float, 3*f_test, FS) if 3*f_test < FS/2 else 0

    # Model EQ gain at this frequency
    idx = np.argmin(np.abs(eq_freqs - f_test))
    eq_gain_db = target_db[idx]

    results.append((f_test, rms_out, fund_amp, h2_amp, h3_amp, eq_gain_db))

    print(f"  {f_test:6.0f}  {rms_out:8.4f}  {20*np.log10(max(rms_out,1e-6)):+6.1f}  "
          f"{fund_amp:7.4f}  {h2_amp:7.4f}  {h3_amp:7.4f}  "
          f"{eq_gain_db:+8.1f}dB  ", end="")

effi.destroy_enhancer(enh)

# ── Flatness summary ───────────────────────────────────────────────────
print(f"\n\n  ── Flatness check ──")
rmss = np.array([r[1] for r in results])
ref = np.mean(rmss)
print(f"  Mean output RMS: {ref:.4f} ({20*np.log10(ref):+.1f} dBFS)")
print(f"  Min/Max ratio:   {np.min(rmss)/np.max(rmss):.3f} ({20*np.log10(np.min(rmss)/np.max(rmss)):+.1f} dB spread)")

# Per-region flatness
bass_rms  = [r[1] for r in results if r[0] <= 120]
treble_rms = [r[1] for r in results if r[0] > 120]
if bass_rms:
    print(f"  Bass (20-120 Hz) spread: {20*np.log10(max(bass_rms)/max(min(bass_rms),1e-6)):.1f} dB")
if treble_rms:
    print(f"  Treble (>120 Hz) spread: {20*np.log10(max(treble_rms)/max(min(treble_rms),1e-6)):.1f} dB")
