"""
Level-linearity verification using Goertzel single-frequency analysis.

No FFT bin misalignment. Measures exact amplitude at f, 2f, 3f.
Anything else in the signal is classified as noise.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    process_original,
    process_linear,
    generate_test_sine,
)
from fix_verify import process_fixed
from goertzel_analyzer import measure_tones, format_tone_report


def analyze_variant(
    label: str,
    process_fn,
    freq: float,
    cfg: BassEnhancerConfig,
    amplitudes: list = None,
):
    """Full linearity sweep for one variant using Goertzel analysis."""
    if amplitudes is None:
        amplitudes = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09]

    duration = 0.5
    pad = int(0.3 * cfg.fs)

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  f={freq} Hz, cutoff={cfg.cutoff} Hz, h2={cfg.h2_amp}, h3={cfg.h3_amp}")
    print(f"{'='*70}")

    print(f"\n  {'In dBFS':>8s}  {'Fund dBFS':>10s}  "
          f"{'H2 rel':>8s}  {'H3 rel':>8s}  "
          f"{'Noise dB':>9s}  {'SNR dB':>8s}")
    print(f"  {'-'*56}")

    h2_rel_list = []
    h3_rel_list = []
    noise_list = []
    snr_list = []

    for amp in amplitudes:
        signal = generate_test_sine([freq], [amp], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out = process_fn(stereo, cfg)
        start = pad + len(signal) // 2
        end = pad + len(signal)
        mono = (out[0, start:end] + out[1, start:end]) / 2.0

        m = measure_tones(mono, freq, cfg.fs, harmonics=(1, 2, 3))
        tones = m["tones"]

        fund_db = 20 * np.log10(max(tones.get(1, 0), 1e-12))
        h2_rel = 20 * np.log10(max(tones.get(2, 0), 1e-12) / max(tones.get(1, 0), 1e-12))
        h3_rel = 20 * np.log10(max(tones.get(3, 0), 1e-12) / max(tones.get(1, 0), 1e-12))
        noise_db = 20 * np.log10(max(m["noise_rms"], 1e-12))
        snr = m.get("snr_db", -np.inf)

        h2_rel_list.append(h2_rel)
        h3_rel_list.append(h3_rel)
        noise_list.append(noise_db)
        snr_list.append(snr)

        db_in = 20 * np.log10(amp)
        print(f"  {db_in:+8.2f}  {fund_db:+10.2f}  "
              f"{h2_rel:+8.2f}  {h3_rel:+8.2f}  "
              f"{noise_db:+9.2f}  {snr:+8.1f}")

    # Summary stats
    valid_h2 = [v for v in h2_rel_list if v > -200]
    valid_h3 = [v for v in h3_rel_list if v > -200]

    if valid_h2:
        h2_mean, h2_std, h2_range = np.mean(valid_h2), np.std(valid_h2), max(valid_h2) - min(valid_h2)
        linear = "✅ LINEAR" if h2_range < 1.5 else f"❌ NONLINEAR (±{h2_range:.1f} dB spread)"
        print(f"\n  H2:  mean={h2_mean:+.2f} dB  σ={h2_std:.2f} dB  range={h2_range:.2f} dB  {linear}")

    if valid_h3:
        h3_mean, h3_std, h3_range = np.mean(valid_h3), np.std(valid_h3), max(valid_h3) - min(valid_h3)
        linear = "✅ LINEAR" if h3_range < 1.5 else f"❌ NONLINEAR (±{h3_range:.1f} dB spread)"
        print(f"  H3:  mean={h3_mean:+.2f} dB  σ={h3_std:.2f} dB  range={h3_range:.2f} dB  {linear}")

    print(f"  SNR: mean={np.mean(snr_list):+.1f} dB  "
          f"min={min(snr_list):+.1f} dB  (higher = cleaner)")

    return {"h2_mean": h2_mean if valid_h2 else None,
            "h2_std": h2_std if valid_h2 else None,
            "h3_mean": h3_mean if valid_h3 else None,
            "h3_std": h3_std if valid_h3 else None,
            "snr_mean": np.mean(snr_list)}


def run():
    """Run the level-linearity verification."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, corr_limit=20.0,
    )

    print("=" * 70)
    print("  LEVEL LINEARITY VERIFICATION  (Goertzel — exact frequency)")
    print("=" * 70)

    analyze_variant("ORIGINAL (harmonic_bass_enhancer.eel)", process_original, 50.0, cfg)
    analyze_variant("LINEAR (broken — divide after)", process_linear, 50.0, cfg)
    analyze_variant("FIXED (normalize before)", process_fixed, 50.0, cfg)

    # ── Detailed tone purity check at a single level ──
    print(f"\n{'='*70}")
    print(f"  DETAILED TONE PURITY: 50 Hz, A=0.5 (−6 dBFS)")
    print(f"{'='*70}")

    for label, fn in [("Original", process_original),
                       ("Linear (broken)", process_linear),
                       ("Fixed", process_fixed)]:
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([50.0], [0.5], 0.5, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95
        out = fn(stereo, cfg)
        start = pad + len(signal) // 2
        end = pad + len(signal)
        mono = (out[0, start:end] + out[1, start:end]) / 2.0
        m = measure_tones(mono, 50.0, cfg.fs, harmonics=(1, 2, 3))
        format_tone_report(m, label, fundamental=50.0)


if __name__ == "__main__":
    run()
