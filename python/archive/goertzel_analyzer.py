"""
Precise single-frequency analysis using the Goertzel algorithm.

Avoids FFT bin-misalignment issues. Measures exact amplitude at the
fundamental, 2nd harmonic, and 3rd harmonic, then computes a "noise"
residual: any energy NOT at those three frequencies.
"""

import numpy as np
from typing import Tuple


def goertzel(samples: np.ndarray, freq: float, fs: float) -> complex:
    """Goertzel algorithm — precise single-frequency DFT bin.

    Returns the complex DFT coefficient at exactly `freq`.
    Equivalent to: sum_n samples[n] * exp(-j*2*pi*freq*n/fs)
    """
    n = len(samples)
    omega = 2.0 * np.pi * freq / fs
    coeff = 2.0 * np.cos(omega)

    s0 = 0.0
    s1 = 0.0
    for x in samples:
        s2 = s1
        s1 = s0
        s0 = x + coeff * s1 - s2

    # Final complex value
    real = s0 - s1 * np.cos(omega)
    imag = s1 * np.sin(omega)
    return complex(real, -imag)


def goertzel_magnitude(samples: np.ndarray, freq: float, fs: float) -> float:
    """Amplitude of a single frequency component (linear, not dB)."""
    c = goertzel(samples, freq, fs)
    return abs(c) * 2.0 / len(samples)


def measure_tones(
    samples: np.ndarray,
    fundamental: float,
    fs: float,
    harmonics: Tuple[int, ...] = (1, 2, 3),
    nearby_bw_hz: float = 0.5,
) -> dict:
    """Measure amplitude of fundamental and its harmonics using Goertzel.

    Also measures energy at frequencies just off the tones to detect
    bandwidth smearing (indicates intermodulation / non-tonal content).

    Returns:
        tones: dict of harmonic_number -> amplitude (linear)
        nearby: dict of harmonic_number -> amplitude at freq ± bw (for smear check)
        rms_total: total RMS of the signal
        noise_rms: sqrt(RMS² − Σ tones²)  — energy NOT at any harmonic
    """
    results = {"tones": {}, "nearby": {}, "rms_total": 0.0, "noise_rms": 0.0}

    rms_total = np.sqrt(np.mean(samples**2))
    results["rms_total"] = rms_total

    tone_power = 0.0
    for h in harmonics:
        f = fundamental * h
        amp = goertzel_magnitude(samples, f, fs)
        results["tones"][h] = amp
        tone_power += amp**2 / 2.0  # power of sine wave = A²/2

        # Check nearby frequencies for spectral smear
        amp_lo = goertzel_magnitude(samples, f - nearby_bw_hz, fs)
        amp_hi = goertzel_magnitude(samples, f + nearby_bw_hz, fs)
        results["nearby"][h] = {"lo": amp_lo, "hi": amp_hi}

    noise_power = max(0.0, rms_total**2 - tone_power)
    results["noise_rms"] = np.sqrt(noise_power)

    if tone_power > 0:
        results["snr_db"] = 10.0 * np.log10(tone_power / max(noise_power, 1e-20))
    else:
        results["snr_db"] = -np.inf

    return results


def format_tone_report(measurements: dict, label: str = "", fundamental: float = 50.0):
    """Pretty-print a tone measurement report."""
    tones = measurements["tones"]
    nearby = measurements["nearby"]

    print(f"\n  ── {label} ──")
    print(f"  Total RMS:  {20*np.log10(measurements['rms_total']):+.2f} dBFS")
    print(f"  Noise RMS:  {20*np.log10(max(measurements['noise_rms'],1e-20)):+.2f} dBFS  "
          f"(SNR: {measurements.get('snr_db', -np.inf):+.1f} dB)")

    # Find fundamental amplitude for relative dB
    fund_amp = tones.get(1, 1e-10)

    print(f"  {'Harmonic':>10s}  {'Freq (Hz)':>10s}  {'Amplitude':>10s}  "
          f"{'dBFS':>8s}  {'dB rel f':>10s}  {'Smear ±':>8s}")
    print(f"  {'-'*60}")

    for h in sorted(tones.keys()):
        amp = tones[h]
        db_fs = 20 * np.log10(amp) if amp > 0 else -np.inf
        db_rel = 20 * np.log10(amp / fund_amp) if fund_amp > 0 else -np.inf

        n = nearby.get(h, {})
        smear = max(n.get("lo", 0), n.get("hi", 0))
        smear_db = 20 * np.log10(smear / max(amp, 1e-10)) if amp > 1e-10 else 0

        print(f"  {f'H{h}':>10s}  {fundamental*h:10.1f}  {amp:10.6f}  "
              f"{db_fs:+8.2f}  {db_rel:+10.2f}  {smear_db:+8.2f}")

    # Check for tone purity
    unwanted = 0.0
    for h, amp in tones.items():
        if h not in (1, 2, 3):
            unwanted += amp**2 / 2.0

    noise_power = measurements["noise_rms"]**2
    intermod_power = max(0.0, noise_power - unwanted)
    print(f"\n  Tonal purity: {10*np.log10(max(measurements['rms_total']**2 - noise_power, 1e-20) / max(measurements['rms_total']**2, 1e-20)):+.1f} dB "
          f"(fraction of energy in harmonics f/2f/3f)")

    return measurements
