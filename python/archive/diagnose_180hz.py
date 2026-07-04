"""
Diagnose the 180 Hz overtone / slow fade-in bug.

User reports: cutoff=60 Hz, input=180 Hz → audible overtone that
fades in over seconds when the plugin is enabled.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    design_butter_lp,
    design_butter_hp,
    biquad_tick,
    EnvFollower,
    generate_test_sine,
)
from goertzel_analyzer import measure_tones


def trace_detailed(freq: float = 180.0, amp: float = 1.0, cutoff: float = 60.0, fs: float = 44100.0):
    """Trace every internal signal for a pure tone through the fixed plugin."""

    cfg = BassEnhancerConfig(cutoff=cutoff, h2_amp=0.13, h3_amp=0.10, fs=fs)

    duration = 0.5
    n = int(duration * fs)

    # Generate stereo signal
    t = np.arange(n) / fs
    signal = amp * np.sin(2.0 * np.pi * freq * t)
    stereo = np.zeros((2, n))
    stereo[0, :] = signal
    stereo[1, :] = signal * 0.95

    # Run the fixed plugin, recording internals
    lp_coeffs = design_butter_lp(cutoff, fs)
    hp_coeffs = design_butter_hp(cutoff, fs)

    lp_L = np.zeros(4)
    hp_L = np.zeros(4)
    hp_harm_L = np.zeros(4)
    env_L = EnvFollower.from_params(cfg.env_release, fs)
    floor = 0.0001

    # Recording arrays
    lp_out = np.zeros(n)
    env_out = np.zeros(n)
    norm_out = np.zeros(n)
    t2_out = np.zeros(n)
    t3_out = np.zeros(n)
    harm_out = np.zeros(n)
    final_out = np.zeros(n)

    for i in range(n):
        sL = stereo[0, i]

        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        lp_out[i] = l_lp

        env_L.tick(l_lp)
        env_val = env_L.read()
        env_out[i] = env_val

        env_s = max(env_val, floor)
        norm = l_lp / env_s
        norm_out[i] = norm

        t2 = 2.0 * norm * norm - 1.0
        t3 = 4.0 * norm * norm * norm - 3.0 * norm
        t2_out[i] = t2
        t3_out[i] = t3

        harm = env_s * (cfg.h2_amp * t2 + cfg.h3_amp * t3)
        harm_out[i] = harm

        wet, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        harm_hp, hp_harm_L = biquad_tick(harm, hp_coeffs, hp_harm_L)
        final_out[i] = wet + harm_hp

    # ── Analysis ────────────────────────────────────────────────────
    # Steady-state (last 25%)
    ss_start = 3 * n // 4
    ss_lp = lp_out[ss_start:]
    ss_env = env_out[ss_start:]
    ss_norm = norm_out[ss_start:]
    ss_t2 = t2_out[ss_start:]
    ss_t3 = t3_out[ss_start:]
    ss_harm = harm_out[ss_start:]
    ss_final = final_out[ss_start:]

    # LP filter gain at this frequency
    lp_gain = np.sqrt(2.0 * np.mean(ss_lp**2)) / (amp / np.sqrt(2))

    print(f"\n{'='*70}")
    print(f"  DIAGNOSIS: {freq} Hz tone, A={amp} ({20*np.log10(amp):+.1f} dBFS)")
    print(f"  cutoff={cutoff} Hz, h2={cfg.h2_amp}, h3={cfg.h3_amp}")
    print(f"{'='*70}")

    print(f"\n  ── Steady-state levels ──")
    print(f"  Input amplitude:          {amp:.4f}  ({20*np.log10(amp):+.1f} dBFS)")
    print(f"  LP filter gain:           {lp_gain:.4f}  ({20*np.log10(lp_gain):+.1f} dB)")
    print(f"  LP output amplitude:      {np.sqrt(2*np.mean(ss_lp**2)):.6f}")
    print(f"  Envelope (mean):          {np.mean(ss_env):.6f}")
    print(f"  Envelope (max):           {np.max(ss_env):.6f}")
    print(f"  Envelope (min):           {np.min(ss_env):.6f}")
    print(f"  Envelope ripple:          {20*np.log10(np.max(ss_env)/np.min(ss_env)):.2f} dB")
    print(f"  Is envelope > floor?      {'YES' if np.min(ss_env) > floor else 'NO — floor dominates!'}")

    # Check if norm is actually ~1
    norm_amp = np.sqrt(2 * np.mean(ss_norm**2))
    print(f"  Norm amplitude (should ~1): {norm_amp:.4f}")

    t2_amp = np.sqrt(2 * np.mean(ss_t2**2))
    t3_amp = np.sqrt(2 * np.mean(ss_t3**2))
    print(f"  T₂ amplitude:             {t2_amp:.4f}")
    print(f"  T₃ amplitude:             {t3_amp:.4f}")
    print(f"  Harmonic (pre-HP) RMS:    {np.sqrt(np.mean(ss_harm**2)):.6f}")
    print(f"  Final output RMS:         {np.sqrt(np.mean(ss_final**2)):.6f}")

    # Goertzel analysis of harmonic signal
    m_harm = measure_tones(ss_harm, freq, fs, harmonics=(1, 2, 3))
    print(f"\n  ── Harmonic signal spectrum (Goertzel) ──")
    for h in (1, 2, 3):
        amp_val = m_harm["tones"].get(h, 0)
        db = 20 * np.log10(max(amp_val, 1e-12))
        print(f"  H{h} ({freq*h:.0f} Hz):  {amp_val:.8f}  ({db:+.2f} dBFS)")

    # Goertzel analysis of final output
    m_final = measure_tones(ss_final, freq, fs, harmonics=(1, 2, 3))
    print(f"\n  ── Final output spectrum (Goertzel) ──")
    for h in (1, 2, 3):
        amp_val = m_final["tones"].get(h, 0)
        fund = m_final["tones"].get(1, 1e-12)
        db = 20 * np.log10(max(amp_val, 1e-12))
        db_rel = 20 * np.log10(max(amp_val, 1e-12) / max(fund, 1e-12))
        print(f"  H{h} ({freq*h:.0f} Hz):  {amp_val:.8f}  ({db:+.2f} dBFS, {db_rel:+.2f} dB rel f)")

    print(f"  Noise RMS:  {20*np.log10(max(m_final['noise_rms'], 1e-12)):+.2f} dBFS")
    print(f"  SNR:        {m_final.get('snr_db', -np.inf):+.1f} dB")

    # ── Startup transient ───────────────────────────────────────────
    print(f"\n  ── Startup transient (first 100 ms) ──")
    startup_n = int(0.1 * fs)

    # Find when envelope crosses floor
    floor_cross = np.argmax(env_out[:startup_n] > floor)
    print(f"  Envelope exceeds floor at sample {floor_cross} ({floor_cross/fs*1000:.1f} ms)")

    # Find when LP output reaches steady state (within 10%)
    ss_level = np.sqrt(2 * np.mean(ss_lp**2))
    lp_mag = np.abs(lp_out[:startup_n])
    settle_idx = np.argmax(lp_mag > 0.9 * ss_level)
    if settle_idx > 0:
        print(f"  LP reaches 90% steady-state at sample {settle_idx} ({settle_idx/fs*1000:.1f} ms)")

    # Check the norm signal during startup
    norm_startup = norm_out[:startup_n]
    print(f"  Norm max during startup:  {np.max(np.abs(norm_startup)):.4f}")
    print(f"  Norm RMS during startup:  {np.sqrt(np.mean(norm_startup**2)):.4f}")
    print(f"  Harm max during startup:  {np.max(np.abs(harm_out[:startup_n])):.6f}")

    # ── Envelope decay analysis ─────────────────────────────────────
    print(f"\n  ── Envelope dynamics ──")
    # Simulate: what if there's a sudden silence after loud bass?
    # Check how long the envelope takes to decay to floor
    decay_to_floor = np.log(floor / np.max(ss_env)) / np.log(cfg.env_release * 0.001 * fs)
    # Actually compute properly
    decay_per_sample = np.exp(-1.0 / (max(cfg.env_release, 10.0) * 0.001 * fs))
    env_peak = np.max(ss_env)
    samples_to_floor = int(np.log(floor / env_peak) / np.log(decay_per_sample))
    print(f"  Release time:             {cfg.env_release} ms")
    print(f"  Decay per sample:         {decay_per_sample:.8f}")
    print(f"  Decay per second:         {decay_per_sample**fs:.4f}")
    print(f"  Time for env to decay from {env_peak:.4f} to floor: {samples_to_floor/fs*1000:.0f} ms")

    # ── Frequency sweep around 180 Hz ───────────────────────────────
    print(f"\n  ── Frequency response of harmonic path ──")
    print(f"  {'Freq':>8s}  {'LP gain':>10s}  {'Env≈':>10s}  {'H2 out':>10s}  {'H3 out':>10s}")
    for test_freq in [30, 50, 80, 120, 180, 250, 360, 540]:
        # LP + HP magnitude
        lp_mag = 1.0 / np.sqrt(1.0 + (test_freq / cutoff)**4)
        hp_mag = (test_freq / cutoff)**2 / np.sqrt(1.0 + (test_freq / cutoff)**4)
        hp_2f = (2*test_freq/cutoff)**2 / np.sqrt(1.0 + (2*test_freq/cutoff)**4)
        hp_3f = (3*test_freq/cutoff)**2 / np.sqrt(1.0 + (3*test_freq/cutoff)**4)

        env_ss = amp * lp_mag  # steady-state envelope
        h2_out = env_ss * cfg.h2_amp * hp_2f
        h3_out = env_ss * cfg.h3_amp * hp_3f
        print(f"  {test_freq:8.0f}  {20*np.log10(lp_mag):+10.2f}  {env_ss:10.4f}  "
              f"{20*np.log10(max(h2_out,1e-12)):+10.2f}  {20*np.log10(max(h3_out,1e-12)):+10.2f}")

    return {
        "lp_gain": lp_gain,
        "env_ss": np.mean(ss_env),
        "harm_rms": np.sqrt(np.mean(ss_harm**2)),
        "final_rms": np.sqrt(np.mean(ss_final**2)),
        "tones_final": m_final["tones"],
        "floor_cross_ms": floor_cross / fs * 1000 if floor_cross > 0 else 0,
    }


def run():
    """Run the 180 Hz overtone diagnosis."""
    # Test: 180 Hz at various amplitudes
    for amp in [1.0, 0.5, 0.25]:
        trace_detailed(freq=180.0, amp=amp, cutoff=60.0)

    # Compare: what about 50 Hz (in-band)?
    print("\n" + "="*70)
    print("  COMPARISON: 50 Hz (in-band) vs 180 Hz (out-of-band)")
    print("="*70)
    trace_detailed(freq=50.0, amp=0.5, cutoff=60.0)


if __name__ == "__main__":
    run()
