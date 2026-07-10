"""
Filter design and analysis utilities.

All DSP processing now goes through the C enhancer (src/enhancer.c) via
eqgen.enhancer_ffi.  This module retains only the design-time and
analysis helpers that Python is responsible for:

  - Butterworth biquad design (bilinear transform, float precision)
  - Sample-by-sample biquad tick (for compressor/envelope analysis)
  - Envelope follower (peak-hold, for offline trace analysis)
  - BassEnhancerConfig (parameter container for documentation / tests)

The four process_* Python variants have been removed; the production
pipeline is the C implementation exposed through enhancer_ffi.
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

# ─────────────────────────────────────────────────────────────────────────────
# 2nd-order Butterworth biquad filters (bilinear transform)
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
# Butterworth filter magnitude responses (for analysis)
# ─────────────────────────────────────────────────────────────────────────────

def first_order_lp_mag(f: float, fc: float) -> float:
    """1st-order LP magnitude at frequency f.

    Matches the first-order LP filters used in the C enhancer
    (alpha = 1 - exp(-2π·fc/fs), approximated for fc ≪ fs/2).
    """
    w = f / fc
    return 1.0 / np.sqrt(1.0 + w**2)


def butterworth_lp_mag(f: float, fc: float) -> float:
    """2nd-order Butterworth LP magnitude at frequency f."""
    w = f / fc
    return 1.0 / np.sqrt(1.0 + w**4)


def butterworth_hp_mag(f: float, fc: float) -> float:
    """2nd-order Butterworth HP magnitude at frequency f."""
    w = f / fc
    return w * w / np.sqrt(1.0 + w**4)


# ─────────────────────────────────────────────────────────────────────────────
# Configuration — parameter container for documentation / test annotation
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BassEnhancerConfig:
    """Configuration for the C harmonic bass enhancer.

    This is a documentation / test-annotation container.  The actual
    DSP runs through enhancer_ffi.create_enhancer() which passes these
    parameters to src/enhancer.c.
    """
    cutoff: float = 60.0       # Hz
    h2_amp: float = 1.0       # 2nd harmonic amplitude
    h3_amp: float = 1.0       # 3rd harmonic amplitude
    env_attack: float = 5.0    # ms
    env_release: float = 200.0 # ms
    fs: float = 44100.0        # sample rate
    limiter_release_ms: float = 50.0   # ms — output limiter release
    limiter_threshold: float = 1.0     # linear — peak threshold for harmonic attenuation
