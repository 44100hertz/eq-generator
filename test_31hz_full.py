#!/usr/bin/env python3
"""Full end-to-end 31Hz test — generates WAV, runs through C enhancer, FFT analysis.

The user says there's extreme distortion at 31 Hz and their speakers are big enough
to reproduce it. This test runs the ACTUAL pipeline (not just a simulation).

Steps:
  1. Generate a 31 Hz sine WAV at various amplitudes
  2. Process through process_track() → C enhancer on real float samples
  3. FFT analysis to measure harmonics vs input level
"""

import sys, struct, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from eqgen.presets import PresetManager
from eqgen.pipeline import run_pipeline, design_eq
from eqgen.dsp import pre_gain_from_max_gain
from eqgen.sweep import run_sine_sweep, build_vol_lut, compute_smart_volume
from eqgen.analysis import magnitude_spectrum, find_peaks

# ── Load preset ───────────────────────────────────────────────────
pm = PresetManager()
p = pm.load("technics-standing")

print(f"Preset: {p.name}")
print(f"  fc={p.fc}Hz  max_bands={p.max_bands}  release={p.release}s")
print(f"  speaker_level={p.speaker_level}dB  overboost={p.overboost_db}dB")

# ── Run EQ pipeline ───────────────────────────────────────────────
result = run_pipeline(
    p.resolve_measurements(), p.resolve_target(), p.resolve_noise(),
    bass_enhancer_cutoff=p.fc,
    smooth_exponent=p.smooth_exponent,
    detailed=True,
)
freqs = np.array(result["freqs"])
gains_db = np.array(result["gains_db"])
fs = result["sample_rate"]
max_gain_db = result["max_gain_db"]
efficacy = result["efficacy"]

FC = p.fc
H2 = efficacy["h2_amp"]
H3 = efficacy["h3_amp"]
OVERBOOST = p.overboost_db
SPEAKER_LEVEL = p.speaker_level

print(f"Efficacy: h2={H2:.4f}  h3={H3:.4f}")

# ── Design IIR EQ ─────────────────────────────────────────────────
coeffs_flat, bands, _, _, fitted_db = design_eq(freqs, gains_db, fs, max_bands=p.max_bands)
pre_gain = pre_gain_from_max_gain(max_gain_db)
eq_gain_31hz_db = float(np.interp(31.0, freqs, gains_db))
eq_gain_31hz_lin = 10.0 ** (eq_gain_31hz_db / 20.0)

print(f"\nEQ curve: max_gain={max_gain_db:+.1f}dB  pre_gain={pre_gain:.4f} ({20*np.log10(pre_gain):+.1f}dB)")
print(f"Correction @ 31Hz: {eq_gain_31hz_db:+.2f} dB  ({eq_gain_31hz_lin:.2f}×)")

# ── Volume LUT ────────────────────────────────────────────────────
VOL = 127
vol_gain = build_vol_lut(VOL, SPEAKER_LEVEL, OVERBOOST)
sv = compute_smart_volume(VOL, pre_gain)
print(f"\nVolume LUT (vol={VOL}): {vol_gain:.2f}× ({20*np.log10(vol_gain):+.1f} dB)")
print(f"Smart volume:  pre_gain={sv['pre_gain']:.4f}  shelf={sv['shelf_db']:.1f}dB  boost={sv['boost']:.3f}")

# ══════════════════════════════════════════════════════════════════════
#  FULL PIPELINE: run_sine_sweep with DETAILED FFT output
# ══════════════════════════════════════════════════════════════════════
print(f"\n{'='*75}")
print(f"  FULL PIPELINE: 31Hz tone through actual C enhancer + tanh clamp")
print(f"{'='*75}")

for amp_db in [-60, -40, -20, -12, -6, -3, 0]:
    amp_lin = 10.0 ** (amp_db / 20.0)
    results = run_sine_sweep(
        freqs_hz=[31.0],
        eq_coeffs=coeffs_flat,
        fc=FC, h2=H2, h3=H3,
        amplitude=amp_lin, duration_sec=2.0,
        vol_gain=vol_gain,
        pre_gain=sv['pre_gain'],
        loudness_boost=sv['boost'],
    )
    r = results[31.0]
    fund = r["fundamental"]
    h2 = r["h2"]
    h3 = r["h3"]
    rms_out = r["rms"]

    rms_in = amp_lin / np.sqrt(2)
    gain = 20 * np.log10(max(rms_out / max(rms_in, 1e-20), 1e-12))
    thd = 100.0 * np.sqrt(h2**2 + h3**2) / max(fund, 1e-20)
    h2_rel = 20 * np.log10(max(h2, 1e-12) / max(fund, 1e-12))
    h3_rel = 20 * np.log10(max(h3, 1e-12) / max(fund, 1e-12))

    print(f"\n  In: {amp_lin:6.4f} ({amp_db:+3.0f} dBFS)")
    print(f"  Out RMS: {rms_out:.4f} ({20*np.log10(max(rms_out,1e-12)):+.1f} dBFS)  gain={gain:+.1f} dB")
    print(f"  Fund(31Hz)={fund:.4f}  H2(62Hz)={h2:.4f} ({h2_rel:+.1f} dB)  H3(93Hz)={h3:.4f} ({h3_rel:+.1f} dB)")
    print(f"  THD: {thd:.1f}%")

    # Add more harmonics if we can capture them
    if fund > 0.01 and amp_db >= -12:
        # Check if the tanh is the culprit: what is the peak signal before tanh?
        # After pre_gain, vol_gain, and EQ, the peak signal at 31 Hz should be:
        eq_peak = amp_lin * vol_gain * sv['pre_gain'] * eq_gain_31hz_lin
        print(f"  → Pre-tanh peak estimate: {eq_peak:.4f} (tanh({eq_peak:.3f}) = {np.tanh(eq_peak):.4f})")
        print(f"  → Tanh saturation: {'YES ⚠️' if eq_peak > 0.5 else 'minimal'}")

print(f"\n{'='*75}")
print(f"  ROOT CAUSE ANALYSIS")
print(f"{'='*75}")

# Show the correction curve shape around the problematic region
print(f"\nCorrection curve shape (30-50 Hz band):")
print(f"  Freq    Gain(dB)  Change")
prev = None
for i in range(len(freqs)):
    f = freqs[i]
    if 30 <= f <= 50:
        g = gains_db[i]
        delta = ""
        if prev is not None:
            d = g - prev
            if abs(d) > 3:
                delta = f"  {'🚀' if d > 0 else '🚀'} {d:+.1f}dB"
        print(f"  {f:6.1f}  {g:+8.2f}  {delta}")
        prev = g

# Show the IIR fit error in this region
print(f"\nIIR fit accuracy at 31Hz:")
fit_31hz = float(np.interp(31.0, freqs, fitted_db))
target_31hz = float(np.interp(31.0, freqs, gains_db))
print(f"  Target correction @ 31Hz: {target_31hz:+.2f} dB")
print(f"  IIR fit @ 31Hz:          {fit_31hz:+.2f} dB")
print(f"  Fit error:               {fit_31hz - target_31hz:+.2f} dB")

# Check the EQ response at the IIR output
print(f"\nEQ biquad cascade frequency response at key frequencies:")
for f_test in [20, 25, 31, 35, 43, 50, 70]:
    from eqgen.eq_fit import cascade_response_db
    resp = cascade_response_db(bands, np.array([f_test]), fs)
    print(f"  {f_test:3.0f} Hz: IIR gain={resp[0]:+.2f} dB, target={float(np.interp(f_test, freqs, gains_db)):+.2f} dB")
