#!/usr/bin/env python3
"""What happens to a 31Hz pure tone through the FULL technics-standing DSP pipeline.

Includes: BT volume LUT (overboost), pre-gain, 32-band EQ, LP/HP split,
crossfader, Chebyshev harmonics, loudness shelf, tanh clamp.

Matches exactly what runs on the ESP32 firmware and desktop PipeWire filter."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from eqgen.presets import PresetManager
from eqgen.pipeline import run_pipeline, design_eq, pre_gain_from_max_gain
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

FC = p.fc  # 70
H2 = efficacy["h2_amp"]
H3 = efficacy["h3_amp"]
OVERBOOST = p.overboost_db  # 12
SPEAKER_LEVEL = p.speaker_level  # 60

# ── Design IIR EQ ─────────────────────────────────────────────────
coeffs_flat, bands, _, _, _ = design_eq(freqs, gains_db, fs, max_bands=p.max_bands)
pre_gain = pre_gain_from_max_gain(max_gain_db)
eq_gain_31hz_db = float(np.interp(31.0, freqs, gains_db))
eq_gain_31hz_lin = 10.0 ** (eq_gain_31hz_db / 20.0)

print(f"\nEQ curve: max_gain={max_gain_db:+.1f}dB  pre_gain={pre_gain:.4f} ({20*np.log10(pre_gain):+.1f}dB)")
print(f"Correction @ 31Hz: {eq_gain_31hz_db:+.2f} dB  ({eq_gain_31hz_lin:.2f}×)")
print(f"Harmonic efficacy: h2={H2:.4f}  h3={H3:.4f}")

# ── Volume LUT ────────────────────────────────────────────────────
VOL = 127  # max volume (what filter.c uses; firmware varies)
vol_gain = build_vol_lut(VOL, SPEAKER_LEVEL, OVERBOOST)
sv = compute_smart_volume(VOL, pre_gain)
print(f"\nVolume LUT (vol={VOL}): {vol_gain:.2f}× ({20*np.log10(vol_gain):+.1f} dB)")
print(f"Smart volume:  pre_gain={sv['pre_gain']:.4f}  shelf={sv['shelf_db']:.1f}dB  boost={sv['boost']:.3f}")

# ── Net gain at 31Hz ─────────────────────────────────────────────
net_gain = vol_gain * sv['pre_gain'] * eq_gain_31hz_lin
print(f"\nTotal chain gain @ 31Hz:  vol({vol_gain:.2f}) × pre({sv['pre_gain']:.4f}) × EQ({eq_gain_31hz_lin:.2f})")
print(f"                         = {net_gain:.4f}  ({20*np.log10(net_gain):+.2f} dB)")

# ── Run tone through full pipeline ────────────────────────────────
print(f"\n{'='*70}")
print(f"  FULL PIPELINE: 31Hz tone, vol={VOL}, overboost={OVERBOOST}dB")
print(f"{'='*70}")

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

    rms_in = amp_lin / np.sqrt(2)
    rms_out = r["rms"]
    gain = 20 * np.log10(max(rms_out / max(rms_in, 1e-20), 1e-12))

    print(f"\n  In:  {amp_lin:.4f} ({amp_db:+3.0f} dBFS)  peak={amp_lin:.4f}")
    print(f"  Out: RMS={rms_out:.4f} ({20*np.log10(max(rms_out,1e-12)):+.1f} dBFS)  "
          f"gain={gain:+.1f} dB")
    print(f"  Fund(31Hz)={r['fundamental']:.4f}  "
          f"H2(62Hz)={r['h2']:.6f} ({20*np.log10(max(r['h2'],1e-12)/max(r['fundamental'],1e-12)):+.1f} dB)  "
          f"H3(93Hz)={r['h3']:.6f} ({20*np.log10(max(r['h3'],1e-12)/max(r['fundamental'],1e-12)):+.1f} dB)")
