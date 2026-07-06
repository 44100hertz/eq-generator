"""
DSP engine for the harmonic bass enhancer — Python simulation.

Implements all plugin variants (original EEL, linear, fixed, fixed v2)
for offline analysis, comparison, and debugging.

Variant summary:
  process_original  — separate LP at fc & fc/2, no envelope (original .eel)
  process_linear    — single LP, envelope + level linearisation via division  [BROKEN]
  process_fixed     — normalize-before-Chebyshev, single LP
  process_fixed_v2  — normalize-before-Chebyshev, separate LP at fc & fc/2  [BEST]

Key insight for the "overtone" problem:
  Chebyshev polynomials T₂ and T₃ generate perfect harmonics ONLY for
  pure sinusoids. With real multi-frequency bass content, they produce
  intermodulation products (sum & difference frequencies) that leak
  through the output HP filter as "overtones."
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple


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
    """Configuration shared by all plugin variants."""
    cutoff: float = 60.0       # Hz
    h2_amp: float = 1.0       # 2nd harmonic amplitude
    h3_amp: float = 1.0       # 3rd harmonic amplitude
    env_attack: float = 5.0    # ms (linear variant only)
    env_release: float = 200.0 # ms (linear variant only)
    corr_limit: float = 20.0   # dB (linear variant only)
    fs: float = 44100.0        # sample rate
    limiter_release_ms: float = 50.0   # ms — output limiter release
    limiter_threshold: float = 1.0     # linear — peak threshold for harmonic attenuation


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
        env_follower_floor = 0.0001
        envL2_s = max(env_L2.read(), env_follower_floor)
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


def process_fixed(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """FIXED variant: normalize input BEFORE Chebyshev, scale back after.

    Core fix:
        _norm = _lp / max(env, floor)   # normalize to amplitude ≈ 1
        T₂(_norm) * env                  # clean 2nd harmonic, no DC offset
        T₃(_norm) * env                  # clean 3rd harmonic, no fund. leakage

    Only ONE envelope per channel needed (T₂ and T₃ share it).
    No correction limit needed (no divisions by tiny numbers that blow up).
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    lp_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    lp_L = np.zeros(4); lp_R = np.zeros(4)
    hp_L = np.zeros(4); hp_R = np.zeros(4)
    hp_harm_L = np.zeros(4); hp_harm_R = np.zeros(4)

    env_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    floor = 0.0001

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        r_lp, lp_R = biquad_tick(sR, lp_coeffs, lp_R)

        env_L.tick(l_lp)
        env_R.tick(r_lp)

        envL_s = max(env_L.read(), floor)
        envR_s = max(env_R.read(), floor)

        # Normalize to amplitude ~1 before Chebyshev
        normL = l_lp / envL_s
        normR = r_lp / envR_s

        # T₂(norm) = 2·norm² − 1 = cos(2ωt)  [exact for |norm| ≤ 1]
        # T₃(norm) = 4·norm³ − 3·norm = cos(3ωt)  [exact for |norm| ≤ 1]
        # Scale back by envelope to restore amplitude
        harmL = envL_s * (
            cfg.h2_amp * (2.0 * normL * normL - 1.0)
            + cfg.h3_amp * (4.0 * normL * normL * normL - 3.0 * normL)
        )
        harmR = envR_s * (
            cfg.h2_amp * (2.0 * normR * normR - 1.0)
            + cfg.h3_amp * (4.0 * normR * normR * normR - 3.0 * normR)
        )

        wetL, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        harmL_hp, hp_harm_L = biquad_tick(harmL, hp_harm_coeffs, hp_harm_L)
        wetR, hp_R = biquad_tick(sR, hp_coeffs, hp_R)
        harmR_hp, hp_harm_R = biquad_tick(harmR, hp_harm_coeffs, hp_harm_R)

        out[0, i] = wetL + harmL_hp
        out[1, i] = wetR + harmR_hp

    return out


def process_fixed_v2(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """FIXED v2: normalize-before-Chebyshev + separate LP at fc/2 for T₃.

    Matches the production model (src/main.rs bass_enhancer_preprocess).

    Signal chain per channel:
      1. LP at fc       → T₂ (normalize → 2x²−1 → scale)
      2. LP at fc/2     → T₃ (normalize → 4x³−3x → scale)
      3. HP at fc       → original pass-through
      4. HP at fc       → harmonic mix
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    # T₂ LP at fc, T₃ LP at fc/2
    lp_t2_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    lp_t3_coeffs = design_butter_lp(cfg.cutoff / 2.0, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    # Filter states
    lp_t2_L = np.zeros(4); lp_t2_R = np.zeros(4)
    lp_t3_L = np.zeros(4); lp_t3_R = np.zeros(4)
    hp_L = np.zeros(4); hp_R = np.zeros(4)
    hp_harm_L = np.zeros(4); hp_harm_R = np.zeros(4)

    # Envelope followers — one per path per channel
    env_t2_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t2_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t3_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_t3_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    floor = 0.0001

    # Output peak limiter envelope (linked stereo)
    lim_env = EnvFollower.from_params(cfg.limiter_release_ms, cfg.fs)
    lim_gain = 1.0

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        # LP filters
        l_t2, lp_t2_L = biquad_tick(sL, lp_t2_coeffs, lp_t2_L)
        r_t2, lp_t2_R = biquad_tick(sR, lp_t2_coeffs, lp_t2_R)
        l_t3, lp_t3_L = biquad_tick(sL, lp_t3_coeffs, lp_t3_L)
        r_t3, lp_t3_R = biquad_tick(sR, lp_t3_coeffs, lp_t3_R)

        # Envelope followers
        env_t2_L.tick(l_t2); env_t2_R.tick(r_t2)
        env_t3_L.tick(l_t3); env_t3_R.tick(r_t3)

        # T₂: normalize → Chebyshev → scale
        e2L = max(env_t2_L.read(), floor)
        e2R = max(env_t2_R.read(), floor)
        norm2L = l_t2 / e2L
        norm2R = r_t2 / e2R
        t2L = e2L * cfg.h2_amp * (2.0 * norm2L * norm2L - 1.0)
        t2R = e2R * cfg.h2_amp * (2.0 * norm2R * norm2R - 1.0)

        # T₃: normalize → Chebyshev → scale (using narrower LP at fc/2)
        e3L = max(env_t3_L.read(), floor)
        e3R = max(env_t3_R.read(), floor)
        norm3L = l_t3 / e3L
        norm3R = r_t3 / e3R
        t3L = e3L * cfg.h3_amp * (4.0 * norm3L * norm3L * norm3L - 3.0 * norm3L)
        t3R = e3R * cfg.h3_amp * (4.0 * norm3R * norm3R * norm3R - 3.0 * norm3R)

        harmL = t2L + t3L
        harmR = t2R + t3R

        wetL, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        wetR, hp_R = biquad_tick(sR, hp_coeffs, hp_R)
        harmL_hp, hp_harm_L = biquad_tick(harmL, hp_harm_coeffs, hp_harm_L)
        harmR_hp, hp_harm_R = biquad_tick(harmR, hp_harm_coeffs, hp_harm_R)

        # Apply limiter to harmonics only — dry signal passes through
        # untouched. The envelope rides the total output peak but gain
        # reduction touches ONLY the harmonic path.
        thresh = cfg.limiter_threshold

        totalL = wetL + harmL_hp
        totalR = wetR + harmR_hp

        peak = max(abs(totalL), abs(totalR))
        lim_env.tick(peak)
        env_peak = max(lim_env.read(), thresh)
        lim_gain = thresh / env_peak

        out[0, i] = wetL + harmL_hp * lim_gain
        out[1, i] = wetR + harmR_hp * lim_gain

    return out
