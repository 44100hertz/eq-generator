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
from typing import List, Tuple

# ─────────────────────────────────────────────────────────────────────────────
# 2nd-order Butterworth biquad filters (bilinear transform)
# ─────────────────────────────────────────────────────────────────────────────

SQRT2 = np.sqrt(2.0)  # 1.4142135623730951


# Used only by test_chebyshev.py — prefer eq_fit.design_butter_hp for new code.
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
    """2nd-order Butterworth HP: returns [b0, b1, b2, a1, a2].

    Thin wrapper around eq_fit.design_butter_hp for backward compatibility
    with callers that expect ndarray instead of BiquadCoeffs."""
    from eqgen.eq_fit import design_butter_hp as _design_butter_hp
    bc = _design_butter_hp(fc, fs)
    return np.array([bc.b0, bc.b1, bc.b2, bc.a1, bc.a2])


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


# ═══════════════════════════════════════════════════════════════════════════════
# dB conversion helpers
# ═══════════════════════════════════════════════════════════════════════════════

def db_to_ratio(db: float) -> float:
    return 10.0 ** (db / 20.0)


def ratio_to_db(ratio: float) -> float:
    return 20.0 * np.log10(np.maximum(ratio, 1e-20))


def pre_gain_from_max_gain(max_gain_db: float) -> float:
    """Attenuation factor that keeps biquad internal states below overflow.

    Input is attenuated by this factor before the biquad cascade;
    biquads fit the full (unshifted) correction curve and boost
    back up to unity at peak-gain frequencies.
    """
    return 1.0 / max(1.0, db_to_ratio(max_gain_db))


def compute_overboost_ceiling(
    h2_amp: float, h3_amp: float,
    fc: float, coeffs: List[float],
    pre_gain: float, fs: float = 44100.0,
    push_gain: float = 0.3,
    test_freq: float = 50.0,
    nominal_amp: float = 0.5,
    n_iters: int = 10,
) -> float:
    """Find max safe overboost (dB) by sweeping input amplitude.

    Fixes push_gain low enough that harmonics engage, then binary-searches
    input amplitude to find where the bass-sum limiter begins to clip.
    Returns overboost_db = 20·log10(amp_at_boundary / nominal_amp).
    """
    from eqgen import enhancer_ffi as effi

    def _measure(amp: float) -> float:
        """Return peak output for a sine at test_freq with given amplitude."""
        duration = 0.5
        n_samples = int(duration * fs)
        t = np.arange(n_samples) / fs
        sine = amp * np.sin(2.0 * np.pi * test_freq * t)
        eq_coeffs = list(coeffs) if coeffs else []

        enh = effi.create_enhancer(
            cutoff_hz=fc, h2_amp=h2_amp, h3_amp=h3_amp,
            release_secs=0.2, fs=fs, pre_gain=pre_gain,
            push_gain=push_gain, coeffs=eq_coeffs,
        )
        out = np.zeros(len(sine), dtype=np.float64)
        for i, x in enumerate(sine):
            lf, _rf = effi.process_stereo_frame(enh, float(x), float(x))
            out[i] = float(lf)
        effi.destroy_enhancer(enh)

        settle = int(0.05 * fs)
        return float(np.max(np.abs(out[settle:])))

    lo_amp, hi_amp = nominal_amp, 8.0
    # Expand hi_amp until we find clipping
    for _ in range(4):
        if _measure(hi_amp) >= 1.0:
            break
        hi_amp *= 2.0

    safe_amp = lo_amp
    for _ in range(n_iters):
        mid_amp = np.sqrt(lo_amp * hi_amp)
        if _measure(mid_amp) >= 1.0:
            hi_amp = mid_amp
        else:
            safe_amp = mid_amp
            lo_amp = mid_amp

    return float(20.0 * np.log10(safe_amp / nominal_amp))


# ═══════════════════════════════════════════════════════════════════════════════
# Default EQ coefficients (flat 3-HP cascade)
# ═══════════════════════════════════════════════════════════════════════════════

def build_default_eq_coeffs(fs: float = 44100.0) -> List[float]:
    """Build a flat 3-HP-cascade default EQ as float list."""
    coeffs = []
    for fc in [25, 35, 45]:
        omega = np.tan(np.pi * fc / fs)
        c = 1.0 + SQRT2 * omega + omega * omega
        b0 = 1.0 / c
        b1 = -2.0 / c
        b2 = 1.0 / c
        a1 = (2.0 * (omega * omega - 1.0)) / c
        a2 = (1.0 - SQRT2 * omega + omega * omega) / c
        coeffs.extend([b0, b1, b2, a1, a2])
    return coeffs
