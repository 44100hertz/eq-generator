"""
Isolate the T₃ (4x³ − 3x) path to find why it produces an audible overtone
at 180 Hz / 60 Hz cutoff.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    design_butter_lp,
    biquad_tick,
    EnvFollower,
    generate_test_sine,
)
from goertzel_analyzer import goertzel_magnitude


def analyze_t3_output(
    freq: float,
    amplitude: float,
    cfg: BassEnhancerConfig,
    duration: float = 1.0,
):
    """Extract and measure ONLY the T₃ contribution at the final output."""
    fs = cfg.fs
    n = int(duration * fs)
    pad = int(0.3 * fs)

    # Generate test signal
    signal = generate_test_sine([freq], [amplitude], duration, fs, stereo=False)
    stereo = np.zeros((2, len(signal) + 2 * pad))
    stereo[0, pad:pad + len(signal)] = signal
    stereo[1, :] = stereo[0, :] * 0.95

    # ── Run the full fixed plugin, but also capture T₃ contribution ──
    lp_coeffs = design_butter_lp(cfg.cutoff, fs)
    hp_coeffs = design_butter_lp(cfg.cutoff, fs)  # not used directly here
    hp_out_coeffs = design_butter_lp(cfg.cutoff, fs)
    from harmonic_bass import design_butter_hp
    hp_coeffs = design_butter_hp(cfg.cutoff, fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, fs)

    total_n = stereo.shape[1]
    lp_L = np.zeros(4); hp_L = np.zeros(4); hp_harm_L = np.zeros(4)
    env_L = EnvFollower.from_params(cfg.env_release, fs)
    floor = 0.0001

    t3_only = np.zeros(total_n)
    final_output = np.zeros(total_n)

    for i in range(total_n):
        sL = stereo[0, i]
        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        env_L.tick(l_lp)
        env_s = max(env_L.read(), floor)
        norm = l_lp / env_s

        # T₂ contribution (for reference, not used in T₃-only)
        t2 = cfg.h2_amp * (2.0 * norm * norm - 1.0)
        # T₃ contribution
        t3 = cfg.h3_amp * (4.0 * norm * norm * norm - 3.0 * norm)

        harm_total = env_s * (t2 + t3)
        harm_t3_only = env_s * t3

        wet, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        harm_hp, hp_harm = biquad_tick(harm_total, hp_harm_coeffs, hp_harm_L)

        # Also run T₃ through the same HP
        _, hp_harm_t3 = biquad_tick(harm_t3_only, hp_harm_coeffs,
                                     np.zeros(4) if i == 0 else np.zeros(4))
        # Need separate state for T₃-only path
        # Actually let me just compute final output and T₃ separately

        final_output[i] = wet + harm_hp

    # Re-run to get T₃-only output with its own HP state
    lp_L2 = np.zeros(4)
    hp_harm_t3 = np.zeros(4)
    env_L2 = EnvFollower.from_params(cfg.env_release, fs)
    t3_final = np.zeros(total_n)

    for i in range(total_n):
        sL = stereo[0, i]
        l_lp, lp_L2 = biquad_tick(sL, lp_coeffs, lp_L2)
        env_L2.tick(l_lp)
        env_s = max(env_L2.read(), floor)
        norm = l_lp / env_s
        t3 = cfg.h3_amp * (4.0 * norm * norm * norm - 3.0 * norm)
        harm_t3 = env_s * t3
        harm_hp, hp_harm_t3 = biquad_tick(harm_t3, hp_harm_coeffs, hp_harm_t3)
        t3_final[i] = harm_hp

    # Steady state
    ss = slice(pad + len(signal) // 2, pad + len(signal))

    full_ss = final_output[ss]
    t3_ss = t3_final[ss]

    # ── Measure with Goertzel ─────────────────────────────────────────
    full_rms = np.sqrt(np.mean(full_ss**2))
    t3_rms = np.sqrt(np.mean(t3_ss**2))

    # What does T₃ contribute at each frequency?
    tones_in = [freq, 2*freq, 3*freq, 4*freq, 5*freq, 6*freq]
    print(f"\n  ── T₃-only analysis: {freq} Hz @ {20*np.log10(amplitude):+.0f} dBFS ──")
    print(f"  cutoff={cfg.cutoff} Hz, h3={cfg.h3_amp}")
    print(f"")
    print(f"  Full output RMS:  {full_rms:.6f}  ({20*np.log10(full_rms):+.2f} dBFS)")
    print(f"  T₃-only RMS:      {t3_rms:.6f}  ({20*np.log10(max(t3_rms,1e-12)):+.2f} dBFS)")
    print(f"  T₃/Full ratio:    {20*np.log10(max(t3_rms,1e-12)/max(full_rms,1e-12)):+.2f} dB")
    print(f"")
    print(f"  T₃ spectral content (Goertzel at exact frequencies):")
    print(f"  {'Freq (Hz)':>10s}  {'Amplitude':>10s}  {'dBFS':>8s}  {'Note':>30s}")
    print(f"  {'-'*60}")

    t3_max = 0
    for tone in tones_in:
        amp = goertzel_magnitude(t3_ss, tone, fs)
        db = 20 * np.log10(max(amp, 1e-12))
        t3_max = max(t3_max, amp)

        note = ""
        if abs(tone - freq) < 1:
            note = "← FUNDAMENTAL (leakage!)"
        elif abs(tone - 2*freq) < 2:
            note = "← 2nd harmonic"
        elif abs(tone - 3*freq) < 3:
            note = "← 3rd harmonic (desired)"
        elif abs(tone - 4*freq) < 4:
            note = "← 4th harmonic"
        elif abs(tone - 5*freq) < 5:
            note = "← 5th harmonic"
        elif abs(tone - 6*freq) < 6:
            note = "← 6th harmonic"
        print(f"  {tone:10.0f}  {amp:10.6f}  {db:+8.2f}  {note}")

    # Relate to full signal's fundamental
    full_f = goertzel_magnitude(full_ss, freq, fs)
    if full_f > 1e-12:
        print(f"")
        print(f"  T₃ components relative to output fundamental:")
        for tone in tones_in:
            amp = goertzel_magnitude(t3_ss, tone, fs)
            db_rel = 20 * np.log10(max(amp, 1e-12) / full_f)
            if db_rel > -80:  # only show audible ones
                print(f"    {tone:.0f} Hz: {db_rel:+.1f} dB rel fundamental")

    # Also: measure what the ORIGINAL plugin's T₃ would produce
    lp_t3_coeffs = design_butter_lp(cfg.cutoff / 2.0, fs)
    lp_t3_state = np.zeros(4)
    original_t3 = np.zeros(total_n)
    for i in range(total_n):
        sL = stereo[0, i]
        l_t3, lp_t3_state = biquad_tick(sL, lp_t3_coeffs, lp_t3_state)
        original_t3[i] = cfg.h3_amp * (4.0 * l_t3 * l_t3 * l_t3 - 3.0 * l_t3)

    # HP filter the original T₃
    hp_state = np.zeros(4)
    orig_t3_final = np.zeros(total_n)
    for i in range(total_n):
        orig_t3_final[i], hp_state = biquad_tick(original_t3[i], hp_harm_coeffs, hp_state)

    orig_t3_ss = orig_t3_final[ss]
    orig_t3_rms = np.sqrt(np.mean(orig_t3_ss**2))

    print(f"")
    print(f"  COMPARISON: original plugin T₃ path (LP at fc/2={cfg.cutoff/2:.0f} Hz)")
    print(f"  Original T₃ RMS:  {orig_t3_rms:.6f}  ({20*np.log10(max(orig_t3_rms,1e-12)):+.2f} dBFS)")
    print(f"  Fixed T₃ RMS:     {t3_rms:.6f}  ({20*np.log10(max(t3_rms,1e-12)):+.2f} dBFS)")
    print(f"  Fixed/Original:   {20*np.log10(max(t3_rms,1e-12)/max(orig_t3_rms,1e-12)):+.2f} dB")

    return {"t3_rms": t3_rms, "orig_t3_rms": orig_t3_rms, "t3_tone_amps": {
        tone: goertzel_magnitude(t3_ss, tone, fs) for tone in tones_in
    }}


def run():
    """Run the T₃ path isolation analysis."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  T₃ PATH ISOLATION")
    print("=" * 70)

    # The user's problematic case
    analyze_t3_output(180.0, 0.5, cfg)

    # Normal bass operation for comparison
    analyze_t3_output(50.0, 0.5, cfg)

    # What about at the cutoff?
    analyze_t3_output(60.0, 0.5, cfg)

    # Test: fc/2 LP for T₃ as a fix?
    print(f"\n{'='*70}")
    print(f"  FIX TEST: use fc/2={cfg.cutoff/2:.0f} Hz LP for T₃ only")
    print(f"{'='*70}")

    for freq in [50, 60, 100, 180]:
        r = analyze_t3_output(freq, 0.5, cfg)
        print(f"  {freq:3d} Hz: T₃ RMS = {20*np.log10(max(r['t3_rms'],1e-12)):+.1f} dBFS")


if __name__ == "__main__":
    run()
