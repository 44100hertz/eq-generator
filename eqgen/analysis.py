"""
Analysis utilities for the harmonic bass enhancer.

Includes:
  - Goertzel algorithm for precise single-frequency measurement
  - FFT-based magnitude spectrum and peak finding
  - Harmonic description and intermodulation detection
  - Test signal generation
"""

import numpy as np
from typing import Tuple, List


# ─────────────────────────────────────────────────────────────────────────────
# Goertzel algorithm — precise single-frequency DFT
# ─────────────────────────────────────────────────────────────────────────────

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


# ─────────────────────────────────────────────────────────────────────────────
# FFT-based analysis
# ─────────────────────────────────────────────────────────────────────────────

def magnitude_spectrum(
    samples: np.ndarray,
    fs: float,
    fft_size: int = 16384,
    window: str = "hann",
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute magnitude spectrum (linear, not dB) via FFT with windowing.

    Returns (freqs, magnitudes).
    """
    window_fn = np.hanning(fft_size) if window == "hann" else np.ones(fft_size)
    segment = samples[:fft_size] * window_fn
    spectrum = np.fft.rfft(segment)
    mags = np.abs(spectrum) / fft_size
    freqs = np.fft.rfftfreq(fft_size, 1.0 / fs)
    return freqs, mags


def find_peaks(
    freqs: np.ndarray,
    mags: np.ndarray,
    min_height_ratio: float = 0.01,
    min_dist_hz: float = 5.0,
) -> List[Tuple[float, float]]:
    """Find spectral peaks above `min_height_ratio * max(mags)`."""
    max_mag = mags.max()
    if max_mag == 0:
        return []

    threshold = min_height_ratio * max_mag
    min_dist_bins = max(1, int(min_dist_hz / (freqs[1] - freqs[0])))

    peaks = []
    for i in range(1, len(mags) - 1):
        if mags[i] > threshold and mags[i] > mags[i - 1] and mags[i] > mags[i + 1]:
            peaks.append((freqs[i], mags[i]))

    # Deduplicate within min_dist_bins
    if len(peaks) > 1:
        filtered = [peaks[0]]
        for p in peaks[1:]:
            if p[0] - filtered[-1][0] > min_dist_hz:
                filtered.append(p)
        peaks = filtered

    return peaks


def describe_harmonics(
    fundamental: float,
    peaks: List[Tuple[float, float]],
    tolerance_hz: float = 3.0,
) -> dict:
    """Describe harmonic content relative to a known fundamental.

    Returns dict mapping harmonic number → (freq, magnitude, dB rel fundamental).
    """
    result = {}
    if not peaks:
        return result

    fundamental_peak = None
    for f, m in peaks:
        if abs(f - fundamental) <= tolerance_hz:
            fundamental_peak = m
            break

    ref_mag = fundamental_peak if fundamental_peak else peaks[0][1]

    for harmonic in range(1, 13):  # up to 12th harmonic
        target_f = fundamental * harmonic
        for f, m in peaks:
            if abs(f - target_f) <= tolerance_hz:
                db = 20.0 * np.log10(m / ref_mag) if ref_mag > 0 else -np.inf
                result[harmonic] = (f, m, db)
                break

    return result


def compute_intermodulation(
    freqs: np.ndarray,
    mags: np.ndarray,
    fundamentals: List[float],
    tolerance_hz: float = 3.0,
) -> dict:
    """Identify intermodulation products from a set of fundamentals.

    Checks for peaks at all ωᵢ ± ωⱼ, 2ωᵢ ± ωⱼ, ωᵢ ± 2ωⱼ combinations.
    """
    result = {}
    for i, f1 in enumerate(fundamentals):
        for f2 in fundamentals:  # include self (generates 2f, 3f harmonics)
            # Sum and difference
            for combo_name, target_f in [
                (f"sum_{f1:.0f}+{f2:.0f}", f1 + f2),
                (f"diff_{f1:.0f}-{f2:.0f}", abs(f1 - f2)),
                (f"2*{f1:.0f}+{f2:.0f}", 2 * f1 + f2),
                (f"2*{f1:.0f}-{f2:.0f}", abs(2 * f1 - f2)),
                (f"{f1:.0f}+2*{f2:.0f}", f1 + 2 * f2),
                (f"{f1:.0f}-2*{f2:.0f}", abs(f1 - 2 * f2)),
            ]:
                if target_f < tolerance_hz:
                    continue
                for f, m in find_peaks(freqs, mags):
                    if abs(f - target_f) <= tolerance_hz:
                        result[combo_name] = (target_f, f, m)
                        break
    return result


def generate_test_sine(
    freqs: List[float],
    amplitudes: List[float],
    duration: float,
    fs: float,
    stereo: bool = True,
) -> np.ndarray:
    """Generate a stereo test signal: sum of sine waves at given frequencies.

    Returns shape (2, N) for stereo or (N,) for mono.
    """
    n = int(duration * fs)
    t = np.arange(n) / fs
    signal = np.zeros(n)
    for f, a in zip(freqs, amplitudes):
        signal += a * np.sin(2.0 * np.pi * f * t)

    if stereo:
        # Slight phase offset to prevent perfect correlation
        out = np.zeros((2, n))
        out[0, :] = signal
        out[1, :] = signal * 0.95  # slightly lower right channel
        return out
    return signal
