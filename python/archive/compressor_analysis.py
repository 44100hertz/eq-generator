"""
Compressor transparency analysis.

The harmonic bass enhancer's envelope follower acts as a compressor:
  - Instant attack  (max operation)
  - 200ms release   (exponential decay)
  - Effectively infinite ratio (hard limit via division)

This module isolates the compressor behavior from the harmonic generation.
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass

from harmonic_bass import (
    BassEnhancerConfig,
    design_butter_lp,
    biquad_tick,
    EnvFollower,
)


@dataclass
class CompressorTrace:
    """Complete trace of the compressor chain for one signal."""
    input_signal: np.ndarray
    lp_signal: np.ndarray
    envelope: np.ndarray
    norm_signal: np.ndarray       # lp / max(env, floor)
    gain: np.ndarray              # 1 / max(env, floor) — the compression gain
    fs: float


def trace_compressor(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
    analyze_steady_state: bool = True,
) -> CompressorTrace:
    """Run the compressor chain and capture every intermediate signal.

    Only runs the LP filter + envelope follower + normalization.
    Does NOT apply Chebyshev or output HP — purely the compressor.
    """
    n = len(samples)
    lp_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    lp_state = np.zeros(4)
    env = EnvFollower.from_params(cfg.env_release, cfg.fs)
    floor = 0.0001

    lp_sig = np.zeros(n)
    env_sig = np.zeros(n)
    norm_sig = np.zeros(n)
    gain_sig = np.zeros(n)

    for i in range(n):
        lp_sig[i], lp_state = biquad_tick(samples[i], lp_coeffs, lp_state)
        env.tick(lp_sig[i])
        env_sig[i] = env.read()
        env_s = max(env_sig[i], floor)
        norm_sig[i] = lp_sig[i] / env_s
        gain_sig[i] = 1.0 / env_s

    return CompressorTrace(
        input_signal=samples,
        lp_signal=lp_sig,
        envelope=env_sig,
        norm_signal=norm_sig,
        gain=gain_sig,
        fs=cfg.fs,
    )


def analyze_compressor_transparency(
    freq: float,
    amplitude: float,
    cfg: BassEnhancerConfig,
    duration: float = 1.0,
):
    """Full compressor analysis for a single sine tone.

    Measures:
      - LP filter attenuation (should be heavy above cutoff)
      - Envelope ripple (DC offset + AC ripple from rectification)
      - Gain modulation depth (how much does 1/env vary?)
      - Norm signal purity (should be a clean sine at amplitude ~1)
    """
    fs = cfg.fs
    n = int(duration * fs)
    t = np.arange(n) / fs
    signal = amplitude * np.sin(2.0 * np.pi * freq * t)

    trace = trace_compressor(signal, cfg)

    # Skip first 200ms for steady-state (filter + envelope settling)
    settle = int(0.2 * fs)
    ss = slice(settle, n)

    lp = trace.lp_signal[ss]
    env = trace.envelope[ss]
    norm = trace.norm_signal[ss]
    gain = trace.gain[ss]

    # ── LP filter analysis ────────────────────────────────────────────
    lp_rms = np.sqrt(np.mean(lp**2))
    lp_peak = np.max(np.abs(lp))
    lp_attenuation_db = 20 * np.log10(lp_rms / (amplitude / np.sqrt(2)))

    # ── Envelope analysis ─────────────────────────────────────────────
    env_mean = np.mean(env)
    env_std = np.std(env)
    env_ripple_db = 20 * np.log10(env_std / max(env_mean, 1e-12))
    env_min = np.min(env)
    env_max = np.max(env)
    env_crest = env_max / max(env_min, 1e-12)

    # ── Gain modulation ───────────────────────────────────────────────
    gain_mean = np.mean(gain)
    gain_std = np.std(gain)
    gain_modulation_pct = 100.0 * gain_std / gain_mean
    gain_min = np.min(gain)
    gain_max = np.max(gain)

    # ── Norm signal analysis ──────────────────────────────────────────
    norm_rms = np.sqrt(np.mean(norm**2))
    norm_peak = np.max(np.abs(norm))
    # How close is norm to a perfect sine?
    # Ideal: norm ≈ sin(2πft), RMS = 1/√2 ≈ 0.707
    norm_ideal_rms = 1.0 / np.sqrt(2.0)
    norm_deviation_db = 20 * np.log10(norm_rms / norm_ideal_rms)

    # ── Spectral analysis of norm signal ──────────────────────────────
    # Use Goertzel at f, 2f, 3f to check for distortion in norm
    from goertzel_analyzer import goertzel_magnitude

    norm_f = goertzel_magnitude(norm, freq, fs)
    norm_2f = goertzel_magnitude(norm, 2 * freq, fs)
    norm_3f = goertzel_magnitude(norm, 3 * freq, fs)
    norm_total_rms = np.sqrt(np.mean(norm**2))

    # ── Print report ──────────────────────────────────────────────────
    print(f"\n  ── Compressor analysis: {freq} Hz @ {20*np.log10(amplitude):+.0f} dBFS ──")
    print(f"  Cutoff: {cfg.cutoff} Hz, release: {cfg.env_release} ms")
    print(f"")
    print(f"  LP FILTER:")
    print(f"    Attenuation:   {lp_attenuation_db:+.1f} dB  (RMS: {lp_rms:.6f})")
    print(f"    Peak output:   {lp_peak:.6f}")
    print(f"")
    print(f"  ENVELOPE FOLLOWER:")
    print(f"    Mean:          {env_mean:.6f}")
    print(f"    Std dev:       {env_std:.6f}  ({env_ripple_db:+.1f} dB ripple)")
    print(f"    Range:         {env_min:.6f} – {env_max:.6f}  (crest: {env_crest:.1f}x)")
    print(f"    Release:       {cfg.env_release} ms (decay coeff: {np.exp(-1.0/(cfg.env_release*0.001*fs)):.6f})")
    print(f"")
    print(f"  GAIN MODULATION (1/env):")
    print(f"    Mean:          {gain_mean:.1f}  ({20*np.log10(gain_mean):+.1f} dB)")
    print(f"    Modulation:    ±{gain_modulation_pct:.2f}%  ({20*np.log10((gain_mean+gain_std)/gain_mean):+.2f} dB)")
    print(f"    Range:         {gain_min:.1f} – {gain_max:.1f}")
    print(f"")
    print(f"  NORM SIGNAL (lp / env):")
    print(f"    RMS:           {norm_rms:.4f}  (ideal: {norm_ideal_rms:.4f}, deviation: {norm_deviation_db:+.2f} dB)")
    print(f"    Peak:          {norm_peak:.4f}")
    print(f"    Fund ampl:     {norm_f:.4f}  ({20*np.log10(norm_f/norm_ideal_rms/np.sqrt(2)):+.2f} dB rel ideal)")
    print(f"    2f distortion: {norm_2f:.4f}  ({20*np.log10(max(norm_2f,1e-12)/max(norm_f,1e-12)):+.2f} dBc)")
    print(f"    3f distortion: {norm_3f:.4f}  ({20*np.log10(max(norm_3f,1e-12)/max(norm_f,1e-12)):+.2f} dBc)")

    # ── Envelope startup analysis ─────────────────────────────────────
    startup_samples = int(0.5 * fs)  # first 500ms
    env_startup = trace.envelope[:startup_samples]
    gain_startup = trace.gain[:startup_samples]

    # Time to reach 90% of final value
    final_env = env_mean
    t90_idx = np.argmax(env_startup >= 0.9 * final_env) if final_env > 0 else 0
    t90_ms = t90_idx / fs * 1000

    print(f"")
    print(f"  STARTUP BEHAVIOR:")
    print(f"    Envelope at t=0:   {trace.envelope[0]:.6f}")
    print(f"    Envelope at t=10ms: {trace.envelope[int(0.01*fs)]:.6f}")
    print(f"    Envelope at t=50ms: {trace.envelope[int(0.05*fs)]:.6f}")
    print(f"    Envelope at t=200ms:{trace.envelope[int(0.2*fs)]:.6f}")
    print(f"    Time to 90% steady: {t90_ms:.1f} ms")
    print(f"    Initial gain (1/floor): {1.0/0.0001:.0f}  ({80:.0f} dB!)")

    # ── Attack test: sudden onset ─────────────────────────────────────
    # Place a burst of tone after 100ms of silence
    burst_n = int(0.5 * fs)
    burst = np.zeros(burst_n)
    silence = int(0.1 * fs)
    burst_start = silence
    burst_end = silence + int(0.2 * fs)
    burst[burst_start:burst_end] = amplitude * np.sin(
        2.0 * np.pi * freq * np.arange(burst_end - burst_start) / fs
    )

    burst_trace = trace_compressor(burst, cfg)
    attack_env = burst_trace.envelope[burst_start:burst_start + int(0.05 * fs)]
    attack_gain = burst_trace.gain[burst_start:burst_start + int(0.05 * fs)]

    print(f"")
    print(f"  ATTACK TRANSIENT (silence → tone burst):")
    print(f"    Env at onset+0ms:  {burst_trace.envelope[burst_start]:.6f}")
    print(f"    Env at onset+1ms:  {burst_trace.envelope[burst_start + int(0.001*fs)]:.6f}")
    print(f"    Env at onset+10ms: {burst_trace.envelope[burst_start + int(0.01*fs)]:.6f}")
    print(f"    Gain at onset+0ms: {burst_trace.gain[burst_start]:.1f} ({20*np.log10(burst_trace.gain[burst_start]):+.0f} dB)")
    print(f"    Gain settled:      {gain_mean:.1f} ({20*np.log10(gain_mean):+.0f} dB)")

    # ── Release test: tone → silence ──────────────────────────────────
    release_start = burst_end
    release_len = min(int(0.4 * fs), len(burst_trace.envelope) - release_start)
    release_env_decay = burst_trace.envelope[release_start:release_start + release_len]

    # Time for envelope to fall to 1/e
    peak_env = np.max(burst_trace.envelope[burst_start:release_start])
    target_env = peak_env / np.e
    decay_idx = np.argmax(release_env_decay <= target_env)
    decay_ms = decay_idx / fs * 1000 if decay_idx > 0 else 0

    print(f"")
    print(f"  RELEASE (tone → silence, {release_len/fs*1000:.0f}ms window):")
    print(f"    Peak envelope:     {peak_env:.6f}")
    t50 = min(int(0.05*fs), release_len - 1)
    t200 = min(int(0.2*fs), release_len - 1)
    print(f"    Env at release+50ms:  {release_env_decay[t50]:.4f}")
    print(f"    Env at release+200ms: {release_env_decay[t200]:.4f}")
    print(f"    Time to 1/e ({target_env:.4f}): {decay_ms:.0f} ms (target: {cfg.env_release} ms)")
    print(f"    Gain at release+50ms:  {burst_trace.gain[release_start + t50]:.1f}")
    print(f"    Gain at release+200ms: {burst_trace.gain[release_start + t200]:.1f}")

    return {
        "lp_attenuation_db": lp_attenuation_db,
        "env_ripple_db": env_ripple_db,
        "gain_modulation_pct": gain_modulation_pct,
        "norm_deviation_db": norm_deviation_db,
        "norm_2f_distortion_dbc": 20 * np.log10(max(norm_2f, 1e-12) / max(norm_f, 1e-12)),
        "norm_3f_distortion_dbc": 20 * np.log10(max(norm_3f, 1e-12) / max(norm_f, 1e-12)),
        "t90_ms": t90_ms,
        "decay_ms": decay_ms,
    }


def compressor_sweep(cfg: BassEnhancerConfig):
    """Analyze compressor behavior across frequencies (above and below cutoff)."""
    print(f"\n{'='*70}")
    print(f"  COMPRESSOR FREQUENCY SWEEP")
    print(f"  cutoff={cfg.cutoff} Hz, release={cfg.env_release} ms")
    print(f"{'='*70}")

    freqs = [30, 50, 60, 80, 120, 180, 300, 500, 1000]
    amp = 0.5

    print(f"\n  {'Freq (Hz)':>10s}  {'LP atten':>10s}  "
          f"{'Gain mean':>10s}  {'Gain mod %':>12s}  "
          f"{'Norm dev':>10s}  {'Norm 2f':>10s}  {'Norm 3f':>10s}")
    print(f"  {'-'*72}")

    for freq in freqs:
        r = analyze_compressor_transparency(freq, amp, cfg, duration=0.5)
        print(f"  {freq:10.0f}  {r['lp_attenuation_db']:+9.1f}  "
              f"{20*np.log10(1.0/(amp*np.sqrt(0.5)*10**(r['lp_attenuation_db']/20))):+9.1f}  "
              f"{r['gain_modulation_pct']:+11.2f}  "
              f"{r['norm_deviation_db']:+9.2f}  "
              f"{r['norm_2f_distortion_dbc']:+9.2f}  "
              f"{r['norm_3f_distortion_dbc']:+9.2f}")


def run():
    """Run the compressor transparency analysis."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    # ── Single tone — the user's 180 Hz case ──────────────────────────
    print("=" * 70)
    print("  COMPRESSOR TRANSPARENCY ANALYSIS")
    print("=" * 70)

    # Below cutoff: bass region (normal operation)
    analyze_compressor_transparency(50.0, 0.5, cfg)

    # At cutoff
    analyze_compressor_transparency(60.0, 0.5, cfg)

    # Above cutoff: the user's 180 Hz case
    analyze_compressor_transparency(180.0, 0.5, cfg)

    # ── Sweep ─────────────────────────────────────────────────────────
    compressor_sweep(cfg)

    # ── Release time sensitivity ──────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  RELEASE TIME SENSITIVITY (180 Hz, A=0.5)")
    print(f"{'='*70}")

    for release_ms in [50, 100, 200, 500, 1000]:
        cfg_rel = BassEnhancerConfig(
            cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
            fs=44100.0, env_release=release_ms,
        )
        r = analyze_compressor_transparency(180.0, 0.5, cfg_rel, duration=0.5)
        print(f"  release={release_ms:4d} ms:  gain_mod={r['gain_modulation_pct']:+.2f}%  "
              f"norm_dev={r['norm_deviation_db']:+.2f} dB  "
              f"norm_2f={r['norm_2f_distortion_dbc']:+.2f} dBc  "
              f"decay={r['decay_ms']:.0f} ms")


if __name__ == "__main__":
    run()
