"""
Visualize the LP filter transient and its effect on the compressor/envelope.

When a tone above the cutoff starts, the Butterworth LP rings at its
natural frequency (~cutoff) for ~50-200ms. The envelope follower tracks
this transient as if it were real bass content, then slowly releases.
This creates harmonics of the TRANSIENT frequency, not the input tone.
"""

import numpy as np
from harmonic_bass import (
    BassEnhancerConfig,
    design_butter_lp,
    biquad_tick,
    EnvFollower,
)


def trace_lp_transient(
    freq: float,
    amplitude: float,
    cfg: BassEnhancerConfig,
    duration: float = 0.5,
):
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

    # Envelope
    env = EnvFollower.from_params(cfg.env_release, fs)
    env_sig = np.zeros(n)
    gain_sig = np.zeros(n)
    floor = 0.0001
    for i in range(n):
        env.tick(lp_sig[i])
        env_sig[i] = env.read()
        gain_sig[i] = 1.0 / max(env_sig[i], floor)

    # Normalized
    norm_sig = np.zeros(n)
    for i in range(n):
        e = max(env_sig[i], floor)
        norm_sig[i] = lp_sig[i] / e

    print(f"\n{'='*70}")
    print(f"  LP FILTER TRANSIENT: {freq} Hz tone, cutoff={cfg.cutoff} Hz")
    print(f"{'='*70}")

    # Steady-state values (average over last 100ms)
    ss_start = int(0.4 * fs)
    lp_ss_rms = np.sqrt(np.mean(lp_sig[ss_start:]**2))
    env_ss = np.mean(env_sig[ss_start:])
    gain_ss = np.mean(gain_sig[ss_start:])

    # Find the LP filter ringing frequency by looking at zero crossings
    # in the first 100ms (where transient dominates)
    early = lp_sig[:int(0.1 * fs)]
    zero_crossings = np.where(np.diff(np.signbit(early)))[0]
    if len(zero_crossings) >= 2:
        ring_period = np.mean(np.diff(zero_crossings)) / fs * 2  # period in seconds
        ring_freq = 1.0 / ring_period if ring_period > 0 else 0
        print(f"  LP filter rings at ~{ring_freq:.0f} Hz during startup")

    # Envelope overshoot
    env_peak = np.max(env_sig)
    env_overshoot = env_peak / max(env_ss, 1e-12)
    t_peak_ms = np.argmax(env_sig) / fs * 1000

    # Time for envelope to settle within 10% of steady state
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
    print(f"  Chebyshev polynomials are underdriven → harmonics are {(env_overshoot**2):.1f}× to {(env_overshoot**3):.1f}× quieter than steady state.")
    print(f"")
    print(f"  BUT: the LP filter's ringing at ~{ring_freq:.0f} Hz (near cutoff)")
    print(f"  IS passed to the Chebyshev polynomials at full amplitude,")
    print(f"  generating harmonics at {ring_freq*2:.0f} Hz and {ring_freq*3:.0f} Hz")
    print(f"  that slowly decay as the transient ring dies out.")
    print(f"  This is the 'overtone that fades in slowly' — it's the FILTER's")
    print(f"  natural frequency being harmonically enriched, not the input tone.")

    # Verify by Goertzel: measure 60 Hz component in LP signal during transient
    from goertzel_analyzer import goertzel_magnitude
    transient_slice = lp_sig[:int(0.2 * fs)]  # first 200ms
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
        "ring_freq": ring_freq if len(zero_crossings) >= 2 else None,
        "env_overshoot": env_overshoot,
        "settle_ms": settle_ms,
        "lp_60hz_ratio": lp_60hz_transient / max(lp_60hz_steady, 1e-12),
    }


def run():
    """Run the LP filter transient analysis."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    # The user's exact scenario
    trace_lp_transient(180.0, 0.5, cfg)

    # Also test: 50 Hz bass (normal operation)
    trace_lp_transient(50.0, 0.5, cfg)

    # Test: does higher-order LP help?
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
