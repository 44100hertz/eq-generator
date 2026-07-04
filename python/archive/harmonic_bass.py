"""
Harmonic Bass Enhancer — DSP implementation in Python.

Implements both the original EEL plugin and the new linear variant
for offline analysis, comparison, and debugging.

Key insight for the "overtone" problem:
  Chebyshev polynomials T₂ and T₃ generate perfect harmonics ONLY for
  pure sinusoids. With real multi-frequency bass content, they produce
  intermodulation products (sum & difference frequencies) that leak
  through the output HP filter as "overtones."

  T₂(x) = 2x² − 1:  for x = Σ cos(ωᵢt), gives desired 2ωᵢ harmonics
  PLUS ωᵢ±ωⱼ intermodulation terms at double amplitude.

  T₃(x) = 4x³ − 3x: additionally produces 2ωᵢ±ωⱼ and ωᵢ±2ωⱼ terms.

  The linear variant compounds this by feeding the SAME LP signal to both
  T₂ and T₃ (the original uses fc/2 for T₃, reducing T₃ intermodulation).
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List


# ─────────────────────────────────────────────────────────────────────────────
# 2nd-order Butterworth biquad filters (bilinear transform)
# Matches the EEL implementation exactly.
# ─────────────────────────────────────────────────────────────────────────────

SQRT2 = np.sqrt(2.0)  # 1.4142135623730951


def design_butter_lp(fc: float, fs: float) -> np.ndarray:
    """2nd-order Butterworth LP: returns [b0, b1, b2, a1, a2]."""
    omega = np.tan(np.pi * fc / fs)
    c = 1.0 + SQRT2 * omega + omega * omega
    b0 = omega * omega / c
    b1 = 2.0 * omega * omega / c
    b2 = omega * omega / c
    a1 = 2.0 * (omega * omega - 1.0) / c
    a2 = (1.0 - SQRT2 * omega + omega * omega) / c
    return np.array([b0, b1, b2, a1, a2])


def design_butter_hp(fc: float, fs: float) -> np.ndarray:
    """2nd-order Butterworth HP: returns [b0, b1, b2, a1, a2]."""
    omega = np.tan(np.pi * fc / fs)
    c = 1.0 + SQRT2 * omega + omega * omega
    b0 = 1.0 / c
    b1 = -2.0 / c
    b2 = 1.0 / c
    a1 = 2.0 * (omega * omega - 1.0) / c
    a2 = (1.0 - SQRT2 * omega + omega * omega) / c
    return np.array([b0, b1, b2, a1, a2])


def biquad_tick(x: float, coeffs: np.ndarray, state: np.ndarray) -> Tuple[float, np.ndarray]:
    """Process one sample through a biquad filter.

    Args:
        x: input sample
        coeffs: [b0, b1, b2, a1, a2]
        state: [x1, x2, y1, y2] — will be mutated in-place

    Returns (y, state).
    """
    b0, b1, b2, a1, a2 = coeffs
    y = b0 * x + b1 * state[0] + b2 * state[1] - a1 * state[2] - a2 * state[3]
    state[1] = state[0]
    state[0] = x
    state[3] = state[2]
    state[2] = y
    return y, state


# ─────────────────────────────────────────────────────────────────────────────
# Envelope follower — peak-hold with exponential decay
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EnvFollower:
    """Peak-hold envelope follower with instant attack (max) and release decay."""

    decay_coeff: float
    floor: float = 1e-4
    value: float = 0.0

    @classmethod
    def from_params(cls, release_ms: float, fs: float, floor: float = 1e-4) -> "EnvFollower":
        r_ms = max(release_ms, 10.0)
        decay_coeff = np.exp(-1.0 / (r_ms * 0.001 * fs))
        return cls(decay_coeff=decay_coeff, floor=floor)

    def tick(self, x: float) -> float:
        self.value = max(self.value * self.decay_coeff, abs(x))
        return self.value

    def read(self) -> float:
        return self.value


# ─────────────────────────────────────────────────────────────────────────────
# Plugin variants
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BassEnhancerConfig:
    """Configuration shared by both plugin variants."""
    cutoff: float = 60.0       # Hz
    h2_amp: float = 0.13       # 2nd harmonic amplitude
    h3_amp: float = 0.10       # 3rd harmonic amplitude
    env_attack: float = 5.0    # ms (linear variant only)
    env_release: float = 200.0 # ms (linear variant only)
    corr_limit: float = 20.0   # dB (linear variant only)
    fs: float = 44100.0        # sample rate


def process_original(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """Process stereo audio through the ORIGINAL harmonic_bass_enhancer.eel.

    Uses separate LP filters:
      - T₂ fed by LP at cutoff
      - T₃ fed by LP at cutoff/2 (narrower band = less intermodulation)
    No envelope following or level linearisation.
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    # Design filters
    lp_t2_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    lp_t3_coeffs = design_butter_lp(cfg.cutoff / 2.0, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    # Filter states per channel
    lp_t2_L = np.zeros(4); lp_t2_R = np.zeros(4)
    lp_t3_L = np.zeros(4); lp_t3_R = np.zeros(4)
    hp_L = np.zeros(4);    hp_R = np.zeros(4)

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        # LP for harmonics
        l_t2, lp_t2_L = biquad_tick(sL, lp_t2_coeffs, lp_t2_L)
        r_t2, lp_t2_R = biquad_tick(sR, lp_t2_coeffs, lp_t2_R)
        l_t3, lp_t3_L = biquad_tick(sL, lp_t3_coeffs, lp_t3_L)
        r_t3, lp_t3_R = biquad_tick(sR, lp_t3_coeffs, lp_t3_R)

        # Chebyshev harmonic generation
        harmL = cfg.h2_amp * (2.0 * l_t2 * l_t2 - 1.0) \
              + cfg.h3_amp * (4.0 * l_t3 * l_t3 * l_t3 - 3.0 * l_t3)
        harmR = cfg.h2_amp * (2.0 * r_t2 * r_t2 - 1.0) \
              + cfg.h3_amp * (4.0 * r_t3 * r_t3 * r_t3 - 3.0 * r_t3)

        # Output = HP(original + harmonics)
        wetL, hp_L = biquad_tick(sL + harmL, hp_coeffs, hp_L)
        wetR, hp_R = biquad_tick(sR + harmR, hp_coeffs, hp_R)

        out[0, i] = wetL
        out[1, i] = wetR

    return out


def process_linear(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """Process stereo audio through the LINEAR harmonic_bass_enhancer_linear.eel.

    Same LP signal fed to both T₂ and T₃ (no separate fc/2 filter).
    Peak-hold envelope followers with level linearisation.
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    # Design filters
    lp_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    # Filter states
    lp_L = np.zeros(4); lp_R = np.zeros(4)
    hp_L = np.zeros(4); hp_R = np.zeros(4)
    hp_harm_L = np.zeros(4); hp_harm_R = np.zeros(4)

    # Envelope followers
    env_L2 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R2 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_L3 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R3 = EnvFollower.from_params(cfg.env_release, cfg.fs)

    max_corr = 10.0 ** (cfg.corr_limit / 20.0)

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        # LP for bass extraction
        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        r_lp, lp_R = biquad_tick(sR, lp_coeffs, lp_R)

        # Envelope followers
        env_L2.tick(l_lp)
        env_R2.tick(r_lp)
        env_L3.tick(l_lp)
        env_R3.tick(r_lp)

        # Level linearisation
        envL2_s = max(env_L2.read(), env_follower_floor := 0.0001)
        envR2_s = max(env_R2.read(), env_follower_floor)
        envL3_s = max(env_L3.read(), env_follower_floor)
        envR3_s = max(env_R3.read(), env_follower_floor)

        h2_corrL = min(1.0 / envL2_s, max_corr)
        h3_corrL = min(1.0 / (envL3_s * envL3_s), max_corr)
        h2_corrR = min(1.0 / envR2_s, max_corr)
        h3_corrR = min(1.0 / (envR3_s * envR3_s), max_corr)

        # Harmonic generation
        harmL = cfg.h2_amp * h2_corrL * (2.0 * l_lp * l_lp - 1.0) \
              + cfg.h3_amp * h3_corrL * (4.0 * l_lp * l_lp * l_lp - 3.0 * l_lp)
        harmR = cfg.h2_amp * h2_corrR * (2.0 * r_lp * r_lp - 1.0) \
              + cfg.h3_amp * h3_corrR * (4.0 * r_lp * r_lp * r_lp - 3.0 * r_lp)

        # Output = HP(original) + HP(harmonics)
        wetL, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        harmL_hp, hp_harm_L = biquad_tick(harmL, hp_harm_coeffs, hp_harm_L)
        wetR, hp_R = biquad_tick(sR, hp_coeffs, hp_R)
        harmR_hp, hp_harm_R = biquad_tick(harmR, hp_harm_coeffs, hp_harm_R)

        out[0, i] = wetL + harmL_hp
        out[1, i] = wetR + harmR_hp

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Analysis utilities
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


def find_peaks(freqs: np.ndarray, mags: np.ndarray, min_height_ratio: float = 0.01, min_dist_hz: float = 5.0) -> List[Tuple[float, float]]:
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
                for f, m in peaks(freqs, mags):
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
