"""
Compressor and dynamics tests.

Tests the envelope follower's behavior as a dynamics processor:
  - Compressor transparency (envelope ripple, gain modulation, attack/release)
  - LP filter startup transients and their interaction with the envelope

The harmonic bass enhancer's envelope follower acts as a compressor:
  - Instant attack  (max operation)
  - Configurable release (exponential decay)
  - Effectively infinite ratio (hard limit via division)
"""

import numpy as np
from dataclasses import dataclass

from dsp import (
    BassEnhancerConfig,
    design_butter_lp,
    design_butter_hp,
    biquad_tick,
    EnvFollower,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Compressor trace (from compressor_analysis.py)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CompressorTrace:
    """Complete trace of the compressor chain for one signal."""
    input_signal: np.ndarray
    lp_signal: np.ndarray
    envelope: np.ndarray
    norm_signal: np.ndarray       # lp / max(env, floor)
    gain: np.ndarray              # 1 / max(env, floor) — the compression gain
    fs: float


def trace_compressor(samples, cfg, analyze_steady_state=True):
    """Run the compressor chain and capture every intermediate signal."""
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


def analyze_compressor_transparency(freq, amplitude, cfg, duration=1.0):
    """Full compressor analysis for a single sine tone."""
    fs = cfg.fs
    n = int(duration * fs)
    t = np.arange(n) / fs
    signal = amplitude * np.sin(2.0 * np.pi * freq * t)

    trace = trace_compressor(signal, cfg)

    settle = int(0.2 * fs)
    ss = slice(settle, n)

    lp = trace.lp_signal[ss]
    env = trace.envelope[ss]
    norm = trace.norm_signal[ss]
    gain = trace.gain[ss]

    lp_rms = np.sqrt(np.mean(lp**2))
    lp_attenuation_db = 20 * np.log10(lp_rms / (amplitude / np.sqrt(2)))

    env_mean = np.mean(env)
    env_std = np.std(env)
    env_ripple_db = 20 * np.log10(env_std / max(env_mean, 1e-12))
    env_min = np.min(env)
    env_max = np.max(env)
    env_crest = env_max / max(env_min, 1e-12)

    gain_mean = np.mean(gain)
    gain_std = np.std(gain)
    gain_modulation_pct = 100.0 * gain_std / gain_mean
    gain_min = np.min(gain)
    gain_max = np.max(gain)

    norm_rms = np.sqrt(np.mean(norm**2))
    norm_peak = np.max(np.abs(norm))
    norm_ideal_rms = 1.0 / np.sqrt(2.0)
    norm_deviation_db = 20 * np.log10(norm_rms / norm_ideal_rms)

    from analysis import goertzel_magnitude
    norm_f = goertzel_magnitude(norm, freq, fs)
    norm_2f = goertzel_magnitude(norm, 2 * freq, fs)
    norm_3f = goertzel_magnitude(norm, 3 * freq, fs)

    print(f"\n  ── Compressor analysis: {freq} Hz @ {20*np.log10(amplitude):+.0f} dBFS ──")
    print(f"  Cutoff: {cfg.cutoff} Hz, release: {cfg.env_release} ms")
    print(f"")
    print(f"  LP FILTER:")
    print(f"    Attenuation:   {lp_attenuation_db:+.1f} dB  (RMS: {lp_rms:.6f})")
    print(f"    Peak output:   {np.max(np.abs(lp)):.6f}")
    print(f"")
    print(f"  ENVELOPE FOLLOWER:")
    print(f"    Mean:          {env_mean:.6f}")
    print(f"    Std dev:       {env_std:.6f}  ({env_ripple_db:+.1f} dB ripple)")
    print(f"    Range:         {env_min:.6f} – {env_max:.6f}  (crest: {env_crest:.1f}x)")
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

    # Startup behavior
    startup_samples = int(0.5 * fs)
    env_startup = trace.envelope[:startup_samples]
    final_env = env_mean
    t90_idx = np.argmax(env_startup >= 0.9 * final_env) if final_env > 0 else 0
    t90_ms = t90_idx / fs * 1000

    print(f"")
    print(f"  STARTUP BEHAVIOR:")
    print(f"    Envelope at t=0:    {trace.envelope[0]:.6f}")
    print(f"    Envelope at t=10ms: {trace.envelope[int(0.01*fs)]:.6f}")
    print(f"    Envelope at t=50ms: {trace.envelope[int(0.05*fs)]:.6f}")
    print(f"    Envelope at t=200ms:{trace.envelope[int(0.2*fs)]:.6f}")
    print(f"    Time to 90% steady: {t90_ms:.1f} ms")
    print(f"    Initial gain (1/floor): {1.0/0.0001:.0f}  ({80:.0f} dB!)")

    # Attack transient
    burst_n = int(0.5 * fs)
    burst = np.zeros(burst_n)
    silence = int(0.1 * fs)
    burst_start = silence
    burst_end = silence + int(0.2 * fs)
    burst[burst_start:burst_end] = amplitude * np.sin(
        2.0 * np.pi * freq * np.arange(burst_end - burst_start) / fs
    )
    burst_trace = trace_compressor(burst, cfg)

    print(f"")
    print(f"  ATTACK TRANSIENT (silence → tone burst):")
    print(f"    Env at onset+0ms:  {burst_trace.envelope[burst_start]:.6f}")
    print(f"    Env at onset+1ms:  {burst_trace.envelope[burst_start + int(0.001*fs)]:.6f}")
    print(f"    Env at onset+10ms: {burst_trace.envelope[burst_start + int(0.01*fs)]:.6f}")
    print(f"    Gain at onset+0ms: {burst_trace.gain[burst_start]:.1f} ({20*np.log10(burst_trace.gain[burst_start]):+.0f} dB)")
    print(f"    Gain settled:      {gain_mean:.1f} ({20*np.log10(gain_mean):+.0f} dB)")

    # Release
    release_start = burst_end
    release_len = min(int(0.4 * fs), len(burst_trace.envelope) - release_start)
    release_env_decay = burst_trace.envelope[release_start:release_start + release_len]
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


def compressor_sweep(cfg):
    """Analyze compressor behavior across frequencies."""
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
        gain_mean_db = 20 * np.log10(1.0 / (amp * np.sqrt(0.5) * 10 ** (r['lp_attenuation_db'] / 20)))
        print(f"  {freq:10.0f}  {r['lp_attenuation_db']:+9.1f}  "
              f"{gain_mean_db:+9.1f}  "
              f"{r['gain_modulation_pct']:+11.2f}  "
              f"{r['norm_deviation_db']:+9.2f}  "
              f"{r['norm_2f_distortion_dbc']:+9.2f}  "
              f"{r['norm_3f_distortion_dbc']:+9.2f}")


# ═══════════════════════════════════════════════════════════════════════════════
# LP filter transient analysis (from transient_analysis.py)
# ═══════════════════════════════════════════════════════════════════════════════

def trace_lp_transient(freq, amplitude, cfg, duration=0.5):
    """Trace the LP filter's startup transient and its interaction with
    the envelope follower."""
    fs = cfg.fs
    n = int(duration * fs)
    t = np.arange(n) / fs
    signal = amplitude * np.sin(2.0 * np.pi * freq * t)

    lp_coeffs = design_butter_lp(cfg.cutoff, fs)
    lp_state = np.zeros(4)
    lp_sig = np.zeros(n)
    for i in range(n):
        lp_sig[i], lp_state = biquad_tick(signal[i], lp_coeffs, lp_state)

    env = EnvFollower.from_params(cfg.env_release, fs)
    env_sig = np.zeros(n)
    gain_sig = np.zeros(n)
    floor = 0.0001
    for i in range(n):
        env.tick(lp_sig[i])
        env_sig[i] = env.read()
        gain_sig[i] = 1.0 / max(env_sig[i], floor)

    print(f"\n{'='*70}")
    print(f"  LP FILTER TRANSIENT: {freq} Hz tone, cutoff={cfg.cutoff} Hz")
    print(f"{'='*70}")

    ss_start = int(0.4 * fs)
    lp_ss_rms = np.sqrt(np.mean(lp_sig[ss_start:]**2))
    env_ss = np.mean(env_sig[ss_start:])
    gain_ss = np.mean(gain_sig[ss_start:])

    early = lp_sig[:int(0.1 * fs)]
    zero_crossings = np.where(np.diff(np.signbit(early)))[0]
    ring_freq = None
    if len(zero_crossings) >= 2:
        ring_period = np.mean(np.diff(zero_crossings)) / fs * 2
        ring_freq = 1.0 / ring_period if ring_period > 0 else 0
        print(f"  LP filter rings at ~{ring_freq:.0f} Hz during startup")

    env_peak = np.max(env_sig)
    env_overshoot = env_peak / max(env_ss, 1e-12)
    t_peak_ms = np.argmax(env_sig) / fs * 1000

    settle_threshold = env_ss * 1.1
    settle_idx = np.argmax(env_sig[int(0.02*fs):] <= settle_threshold)
    settle_ms = (settle_idx + int(0.02*fs)) / fs * 1000 if settle_idx > 0 else 0

    print(f"")
    print(f"  Steady state (after {ss_start/fs*1000:.0f}ms):")
    print(f"    LP RMS:        {lp_ss_rms:.6f}  ({20*np.log10(lp_ss_rms/(amplitude/np.sqrt(2))):+.1f} dB atten)")
    print(f"    Envelope:      {env_ss:.6f}")
    print(f"    Gain (1/env):  {gain_ss:.1f}  ({20*np.log10(gain_ss):+.1f} dB)")
    print(f"")
    print(f"  Startup transient:")
    print(f"    Envelope peak: {env_peak:.6f} at t={t_peak_ms:.1f}ms")
    print(f"    Overshoot:     {env_overshoot:.1f}× steady state")
    print(f"    Settle time:   ~{settle_ms:.0f}ms (to within 10% of steady)")
    print(f"")
    print(f"  Consequence: during the first {settle_ms:.0f}ms, the envelope is")
    print(f"  {env_overshoot:.1f}× too high → norm is {env_overshoot:.1f}× too small →")
    print(f"  Chebyshev polynomials are underdriven → harmonics are {(env_overshoot**2):.1f}× to {(env_overshoot**3):.1f}× quieter.")
    print(f"")
    if ring_freq:
        print(f"  BUT: the LP filter's ringing at ~{ring_freq:.0f} Hz (near cutoff)")
        print(f"  IS passed to the Chebyshev polynomials at full amplitude,")
        print(f"  generating harmonics at {ring_freq*2:.0f} Hz and {ring_freq*3:.0f} Hz")
        print(f"  that slowly decay as the transient ring dies out.")

    from analysis import goertzel_magnitude
    transient_slice = lp_sig[:int(0.2 * fs)]
    ss_slice = lp_sig[ss_start:]

    lp_60hz_transient = goertzel_magnitude(transient_slice, cfg.cutoff, fs)
    lp_60hz_steady = goertzel_magnitude(ss_slice, cfg.cutoff, fs)
    lp_f_transient = goertzel_magnitude(transient_slice, freq, fs)
    lp_f_steady = goertzel_magnitude(ss_slice, freq, fs)

    print(f"")
    print(f"  LP frequency content at {cfg.cutoff:.0f} Hz (cutoff frequency):")
    print(f"    During transient (0-200ms): {lp_60hz_transient:.6f}")
    print(f"    Steady state:               {lp_60hz_steady:.6f}")
    print(f"    Ratio transient/steady:      {lp_60hz_transient/max(lp_60hz_steady,1e-12):.0f}×")
    print(f"")
    print(f"  LP content at input frequency ({freq:.0f} Hz):")
    print(f"    During transient (0-200ms): {lp_f_transient:.6f}")
    print(f"    Steady state:               {lp_f_steady:.6f}")

    return {
        "ring_freq": ring_freq,
        "env_overshoot": env_overshoot,
        "settle_ms": settle_ms,
        "lp_60hz_ratio": lp_60hz_transient / max(lp_60hz_steady, 1e-12),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all compressor and dynamics tests."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    # ── Compressor transparency ───────────────────────────────────────
    print("=" * 70)
    print("  COMPRESSOR TRANSPARENCY ANALYSIS")
    print("=" * 70)

    analyze_compressor_transparency(50.0, 0.5, cfg)
    analyze_compressor_transparency(60.0, 0.5, cfg)
    analyze_compressor_transparency(180.0, 0.5, cfg)
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

    # ── LP filter transients ──────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  LP FILTER TRANSIENT ANALYSIS")
    print(f"{'='*70}")
    trace_lp_transient(180.0, 0.5, cfg)
    trace_lp_transient(50.0, 0.5, cfg)

    print(f"\n{'='*70}")
    print(f"  HIGHER CUTOFF = LESS RINGING?")
    print(f"{'='*70}")
    for fc in [60, 80, 100, 120]:
        c = BassEnhancerConfig(
            cutoff=fc, h2_amp=0.13, h3_amp=0.10,
            fs=44100.0, env_release=200.0,
        )
        r = trace_lp_transient(180.0, 0.5, c)
        print(f"  cutoff={fc:3d} Hz: overshoot={r['env_overshoot']:.1f}×  "
              f"settle={r['settle_ms']:.0f}ms  ring={r.get('ring_freq','?'):.0f} Hz  "
              f"60Hz_ratio={r['lp_60hz_ratio']:.0f}×")


if __name__ == "__main__":
    run()
