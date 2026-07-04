"""
Fixed v2: separate LP at fc/2 for T₃ (like the original plugin).

The single-LP design lets too much above-cutoff content reach T₃.
A 2nd-order Butterworth LP at 60 Hz only gives −19 dB at 180 Hz,
so T₃ generates an audible 540 Hz "overtone" at −39 dB rel.

Giving T₃ its own LP at fc/2 adds 12 dB/octave more rejection
above cutoff, matching the original plugin's spectral behavior
while keeping the normalize-before-Chebyshev linearisation.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    process_original,
    design_butter_lp,
    design_butter_hp,
    biquad_tick,
    EnvFollower,
    generate_test_sine,
)
from goertzel_analyzer import goertzel_magnitude, format_tone_report, measure_tones
from fix_verify import process_fixed


def process_fixed_v2(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """FIXED v2: separate LP at fc/2 for T₃ path.

    Signal chain per channel:
      1. LP at fc       → T₂ (normalize → 2x²−1 → scale)
      2. LP at fc/2     → T₃ (normalize → 4x³−3x → scale)
      3. HP at fc       → original pass-through
      4. HP at fc       → harmonic mix
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    # T₂ LP at fc, T₃ LP at fc/2
    lp_t2_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    lp_t3_coeffs = design_butter_lp(cfg.cutoff / 2.0, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    # Filter states
    lp_t2_L = np.zeros(4); lp_t2_R = np.zeros(4)
    lp_t3_L = np.zeros(4); lp_t3_R = np.zeros(4)
    hp_L = np.zeros(4); hp_R = np.zeros(4)
    hp_harm_L = np.zeros(4); hp_harm_R = np.zeros(4)

    # Envelope followers — one per path per channel
    env_t2_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t2_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t3_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t3_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    floor = 0.0001

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        # LP filters
        l_t2, lp_t2_L = biquad_tick(sL, lp_t2_coeffs, lp_t2_L)
        r_t2, lp_t2_R = biquad_tick(sR, lp_t2_coeffs, lp_t2_R)
        l_t3, lp_t3_L = biquad_tick(sL, lp_t3_coeffs, lp_t3_L)
        r_t3, lp_t3_R = biquad_tick(sR, lp_t3_coeffs, lp_t3_R)

        # Envelope followers
        env_t2_L.tick(l_t2); env_t2_R.tick(r_t2)
        env_t3_L.tick(l_t3); env_t3_R.tick(r_t3)

        # T₂: normalize → Chebyshev → scale
        e2L = max(env_t2_L.read(), floor)
        e2R = max(env_t2_R.read(), floor)
        norm2L = l_t2 / e2L
        norm2R = r_t2 / e2R
        t2L = e2L * cfg.h2_amp * (2.0 * norm2L * norm2L - 1.0)
        t2R = e2R * cfg.h2_amp * (2.0 * norm2R * norm2R - 1.0)

        # T₃: normalize → Chebyshev → scale (using narrower LP at fc/2)
        e3L = max(env_t3_L.read(), floor)
        e3R = max(env_t3_R.read(), floor)
        norm3L = l_t3 / e3L
        norm3R = r_t3 / e3R
        t3L = e3L * cfg.h3_amp * (4.0 * norm3L * norm3L * norm3L - 3.0 * norm3L)
        t3R = e3R * cfg.h3_amp * (4.0 * norm3R * norm3R * norm3R - 3.0 * norm3R)

        harmL = t2L + t3L
        harmR = t2R + t3R

        wetL, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        wetR, hp_R = biquad_tick(sR, hp_coeffs, hp_R)
        harmL_hp, hp_harm_L = biquad_tick(harmL, hp_harm_coeffs, hp_harm_L)
        harmR_hp, hp_harm_R = biquad_tick(harmR, hp_harm_coeffs, hp_harm_R)

        out[0, i] = wetL + harmL_hp
        out[1, i] = wetR + harmR_hp

    return out


def compare_t3_rejection():
    """Compare T₃ output between v1 (single LP) and v2 (separate fc/2 LP)."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  T₃ REJECTION ABOVE CUTOFF: v1 (single LP) vs v2 (fc/2 LP)")
    print("=" * 70)

    for freq in [50, 60, 80, 100, 120, 180, 300]:
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [0.5], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out_v1 = process_fixed(stereo, cfg)
        out_v2 = process_fixed_v2(stereo, cfg)

        start = pad + len(signal) // 2
        end = pad + len(signal)
        mono_v1 = (out_v1[0, start:end] + out_v1[1, start:end]) / 2.0
        mono_v2 = (out_v2[0, start:end] + out_v2[1, start:end]) / 2.0

        # Measure 3rd harmonic content
        h3_v1 = goertzel_magnitude(mono_v1, 3 * freq, cfg.fs)
        h3_v2 = goertzel_magnitude(mono_v2, 3 * freq, cfg.fs)
        fund_v1 = goertzel_magnitude(mono_v1, freq, cfg.fs)
        fund_v2 = goertzel_magnitude(mono_v2, freq, cfg.fs)

        h3_rel_v1 = 20 * np.log10(max(h3_v1, 1e-12) / max(fund_v1, 1e-12))
        h3_rel_v2 = 20 * np.log10(max(h3_v2, 1e-12) / max(fund_v2, 1e-12))
        delta = h3_rel_v2 - h3_rel_v1

        rms_v1 = np.sqrt(np.mean(mono_v1**2))
        rms_v2 = np.sqrt(np.mean(mono_v2**2))
        rms_delta = 20 * np.log10(rms_v2 / rms_v1) if rms_v1 > 0 else 0

        status = "✅ quieter" if delta < -1 else ("⚠️ same" if abs(delta) < 2 else "❌ LOUDER")
        print(f"  {freq:3d} Hz:  H3 rel: v1={h3_rel_v1:+.1f} dB  v2={h3_rel_v2:+.1f} dB  "
              f"Δ={delta:+.1f} dB  {status}  (RMS Δ={rms_delta:+.1f} dB)")

    # ── Detailed tone purity comparison at 180 Hz ─────────────────────
    print(f"\n  ── Detailed T₃ output at 180 Hz ──")
    freq = 180.0
    pad = int(0.3 * cfg.fs)
    signal = generate_test_sine([freq], [0.5], 0.5, cfg.fs, stereo=False)
    stereo = np.zeros((2, len(signal) + 2 * pad))
    stereo[0, pad:pad + len(signal)] = signal
    stereo[1, :] = stereo[0, :] * 0.95

    for label, fn in [("v1 (single LP at fc)", process_fixed),
                       ("v2 (separate fc/2 LP)", process_fixed_v2)]:
        out = fn(stereo, cfg)
        start = pad + len(signal) // 2
        end = pad + len(signal)
        mono = (out[0, start:end] + out[1, start:end]) / 2.0
        m = measure_tones(mono, freq, cfg.fs, harmonics=(1, 2, 3))
        format_tone_report(m, label, fundamental=freq)


def run():
    """Run the v2 fix analysis."""
    compare_t3_rejection()


if __name__ == "__main__":
    run()
