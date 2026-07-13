"""
Random-restart stochastic hill climbing IIR biquad fitter.

Optimizes peaking filter parameters (f0, gain_db, Q) to minimize
perceptually-weighted (ISO 226:2023) RMS error against a target curve.

Key optimization: caches per-band complex responses so mutations
only recompute the changed band, not the whole cascade.
"""

import numpy as np
from typing import List, Tuple

from eqgen.eq_fit import (BiquadCoeffs, design_peaking, cascade_response_db,
                           fit_eq_curve, biquad_response)
from eqgen.model import ear_sensitivity

N_EVAL = 512
GAIN_LO, GAIN_HI = -24.0, 24.0
Q_LO, Q_HI = 0.3, 6.0


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _build_eval_grid(freqs: np.ndarray, target_db: np.ndarray):
    """Log-spaced eval grid matching the greedy fitter's internal grid."""
    f_min, f_max = float(freqs[0]), float(freqs[-1])
    eval_freqs = np.logspace(np.log10(f_min), np.log10(f_max), N_EVAL)
    eval_target = np.interp(eval_freqs, freqs, target_db)
    ear_weights = np.array([ear_sensitivity(f) for f in eval_freqs])
    return eval_freqs, eval_target, ear_weights


def _bands_to_biquads(bands: List[Tuple[float, float, float]],
                       fs: float) -> List[BiquadCoeffs]:
    """Convert list of (f0, gain_db, Q) tuples to BiquadCoeffs list."""
    return [design_peaking(f0, g, q, fs) for f0, g, q in bands]


# ─────────────────────────────────────────────────────────────────────────────
# Fast cascade with per-band response caching
# ─────────────────────────────────────────────────────────────────────────────

class CascadeState:
    """State for fast cascade response evaluation.

    Stores the complex response of each band separately.
    The cascade response is the element-wise product of all band responses.
    When one band changes, only its response needs to be recomputed.
    """

    def __init__(self, bands: List[Tuple[float, float, float]],
                 eval_freqs: np.ndarray, fs: float):
        self.eval_freqs = eval_freqs
        self.fs = fs
        self.n_bands = len(bands)
        # Per-band complex response: shape (n_bands, n_freqs)
        self._band_H = np.empty((self.n_bands, len(eval_freqs)), dtype=complex)
        for i, (f0, gain, Q) in enumerate(bands):
            bc = design_peaking(f0, gain, Q, fs)
            self._band_H[i] = biquad_response(bc, eval_freqs, fs)

    @property
    def cascade_H(self) -> np.ndarray:
        """Element-wise product of all band responses (complex)."""
        return np.prod(self._band_H, axis=0)

    @property
    def cascade_db(self) -> np.ndarray:
        """Magnitude response in dB."""
        mag = np.abs(self.cascade_H)
        mag = np.maximum(mag, 1e-12)
        return 20.0 * np.log10(mag)

    def update_band(self, idx: int, f0: float, gain: float, Q: float):
        """Replace band idx and update its cached response."""
        bc = design_peaking(f0, gain, Q, self.fs)
        self._band_H[idx] = biquad_response(bc, self.eval_freqs, self.fs)

    def score(self, eval_target: np.ndarray, ear_weights: np.ndarray) -> float:
        """Perceptual RMS error."""
        err = self.cascade_db - eval_target
        return float(np.sqrt(np.mean((err * ear_weights) ** 2)))

    def residual_db(self, eval_target: np.ndarray) -> np.ndarray:
        """Residual in dB: target - cascade."""
        return eval_target - self.cascade_db


def _score_bands(bands: List[Tuple[float, float, float]],
                 eval_freqs: np.ndarray, eval_target: np.ndarray,
                 fs: float, ear_weights: np.ndarray) -> float:
    """Perceptual RMS error (convenience, not for inner loop)."""
    cs = CascadeState(bands, eval_freqs, fs)
    return cs.score(eval_target, ear_weights)


# ─────────────────────────────────────────────────────────────────────────────
# Random / greedy initialization
# ─────────────────────────────────────────────────────────────────────────────

def _random_bands(n_bands: int, eval_freqs: np.ndarray,
                  ) -> List[Tuple[float, float, float]]:
    """Random initial bands: log-spaced f0s, small random gains/Q."""
    log_lo = np.log10(eval_freqs[0])
    log_hi = np.log10(eval_freqs[-1])
    base = np.linspace(log_lo, log_hi, n_bands + 2)[1:-1]
    jitter = np.random.uniform(-0.15, 0.15, n_bands)
    f0s = 10 ** (base + jitter)
    f0s = np.clip(f0s, eval_freqs[0], eval_freqs[-1])
    f0s = np.sort(f0s)
    gains = np.random.uniform(-6, 6, n_bands)
    Qs = np.random.uniform(Q_LO, Q_HI, n_bands)
    return [(float(f0s[i]), float(gains[i]), float(Qs[i])) for i in range(n_bands)]


def _greedy_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
                  max_bands: int) -> List[Tuple[float, float, float]]:
    """Run greedy fitter, return (f0, gain, Q) tuples."""
    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=freqs[0])
    return [(b["f0"], b["gain_db"], b["Q"]) for b in fit.bands]


# ─────────────────────────────────────────────────────────────────────────────
# Hill climbing core (fast, with CascadeState caching)
# ─────────────────────────────────────────────────────────────────────────────

def _try_gain_adjust(cs: CascadeState, band_idx: int,
                     current_bands: List[Tuple[float, float, float]],
                     eval_target: np.ndarray, ear_weights: np.ndarray,
                     dg_values: List[float]) -> float:
    """Try gain adjustments on one band. Returns best new gain, or None."""
    f0, gain, Q = current_bands[band_idx]
    best_g = None
    best_s = float('inf')
    for dg in dg_values:
        new_g = np.clip(gain + dg, GAIN_LO, GAIN_HI)
        if abs(new_g - gain) < 0.01:
            continue
        cs.update_band(band_idx, f0, new_g, Q)
        s = cs.score(eval_target, ear_weights)
        if s < best_s:
            best_s = s
            best_g = new_g
    # Restore original
    cs.update_band(band_idx, f0, gain, Q)
    return best_g if best_g is not None and abs(best_g - gain) > 0.01 else None


def _try_q_adjust(cs: CascadeState, band_idx: int,
                  current_bands: List[Tuple[float, float, float]],
                  eval_target: np.ndarray, ear_weights: np.ndarray) -> float:
    """Try Q adjustments. Returns best new Q, or None."""
    f0, gain, Q = current_bands[band_idx]
    best_q = None
    best_s = float('inf')
    for new_Q in [0.5, 1.0, 2.0, 4.0, 6.0]:
        new_Q = np.clip(new_Q, Q_LO, Q_HI)
        if abs(new_Q - Q) < 0.02:
            continue
        cs.update_band(band_idx, f0, gain, new_Q)
        s = cs.score(eval_target, ear_weights)
        if s < best_s:
            best_s = s
            best_q = new_Q
    # Also try ×√2 and ÷√2
    for factor in [1.4, 0.7]:
        new_Q = np.clip(Q * factor, Q_LO, Q_HI)
        if abs(new_Q - Q) < 0.02:
            continue
        cs.update_band(band_idx, f0, gain, new_Q)
        s = cs.score(eval_target, ear_weights)
        if s < best_s:
            best_s = s
            best_q = new_Q
    cs.update_band(band_idx, f0, gain, Q)
    return best_q


def _try_f0_adjust(cs: CascadeState, band_idx: int,
                   current_bands: List[Tuple[float, float, float]],
                   eval_target: np.ndarray, ear_weights: np.ndarray,
                   freq_limits: Tuple[float, float]) -> float:
    """Try f0 shifts. Returns best new f0, or None."""
    f0, gain, Q = current_bands[band_idx]
    best_f = None
    best_s = float('inf')
    f_lo, f_hi = freq_limits
    for df_oct in [-0.15, -0.05, 0.05, 0.15]:
        new_f0 = np.clip(f0 * (2 ** df_oct), f_lo, f_hi)
        if abs(new_f0 - f0) < 0.5:
            continue
        cs.update_band(band_idx, new_f0, gain, Q)
        s = cs.score(eval_target, ear_weights)
        if s < best_s:
            best_s = s
            best_f = new_f0
    cs.update_band(band_idx, f0, gain, Q)
    return best_f


def _hill_climb(initial_bands: List[Tuple[float, float, float]],
                eval_freqs: np.ndarray, eval_target: np.ndarray,
                fs: float, ear_weights: np.ndarray,
                max_stuck: int = 6) -> Tuple[List[Tuple[float, float, float]], float]:
    """Hill climb from initial band placement using CascadeState caching.

    Each iteration:
      1. Compute residual, find frequency with worst error.
      2. Try mutations on bands nearest to that frequency.
      3. Accept first improvement.

    Returns (best_bands, best_score).
    """
    bands = list(initial_bands)
    nb = len(bands)
    if nb == 0:
        return bands, float('inf')

    cs = CascadeState(bands, eval_freqs, fs)
    best_score = cs.score(eval_target, ear_weights)
    freq_limits = (float(eval_freqs[0]), float(eval_freqs[-1]))

    stuck = 0
    iters = 0

    while stuck < max_stuck:
        iters += 1
        if iters > 200:
            break

        # Find worst residual frequency
        residual_db = cs.residual_db(eval_target)
        abs_res = np.abs(residual_db)
        idx_worst = int(np.argmax(abs_res))
        f_worst = float(eval_freqs[idx_worst])
        resid_worst = float(residual_db[idx_worst])

        # Order bands by log-distance to f_worst
        log_dists = np.abs(np.log10([b[0] for b in bands]) - np.log10(f_worst))
        order = np.argsort(log_dists)

        improved = False

        # Phase 1: try nearest bands with focused mutations
        for band_idx in order[:min(4, nb)]:
            f0, gain, Q = bands[band_idx]

            # Try gain: bias toward -residual
            gain_bias = -resid_worst * 0.4
            dgs = sorted([gain_bias, -1.5, -0.5, 0.5, 1.5],
                         key=lambda x: -abs(x))  # bigger first
            new_g = _try_gain_adjust(cs, band_idx, bands, eval_target, ear_weights, dgs)
            if new_g is not None:
                bands[band_idx] = (f0, new_g, Q)
                best_score = cs.score(eval_target, ear_weights)
                improved = True
                stuck = 0
                break

            # Try Q
            new_q = _try_q_adjust(cs, band_idx, bands, eval_target, ear_weights)
            if new_q is not None:
                bands[band_idx] = (f0, gain, new_q)
                best_score = cs.score(eval_target, ear_weights)
                improved = True
                stuck = 0
                break

            # Try f0 shift toward worst
            new_f = _try_f0_adjust(cs, band_idx, bands, eval_target, ear_weights,
                                    freq_limits)
            if new_f is not None:
                bands[band_idx] = (new_f, gain, Q)
                best_score = cs.score(eval_target, ear_weights)
                improved = True
                stuck = 0
                break

        if improved:
            continue

        # Phase 2: wider search — try all bands with larger steps
        for band_idx in order[:min(nb, 8)]:
            f0, gain, Q = bands[band_idx]

            # Wider gain range
            for dg in [-3.0, -1.0, 1.0, 3.0]:
                new_g = np.clip(gain + dg, GAIN_LO, GAIN_HI)
                if abs(new_g - gain) < 0.1:
                    continue
                cs.update_band(band_idx, f0, new_g, Q)
                s = cs.score(eval_target, ear_weights)
                if s < best_score:
                    bands[band_idx] = (f0, new_g, Q)
                    best_score = s
                    improved = True
                    stuck = 0
                    break
            if improved:
                break

        if improved:
            continue

        # Phase 3: replace weakest band with one at worst freq
        weakest_idx = int(np.argmin(np.abs([b[1] for b in bands])))
        f0_w, gain_w, Q_w = bands[weakest_idx]
        if abs(np.log10(f0_w) - np.log10(f_worst)) > 0.1:
            # Try various replacements
            found = False
            for new_g in np.clip([resid_worst * m for m in [0.6, 0.8, 0.5, 1.2]],
                                  GAIN_LO, GAIN_HI):
                new_g = float(new_g)
                for new_Q in [1.5, 2.0, 3.0]:
                    cs.update_band(weakest_idx, f_worst, new_g, new_Q)
                    s = cs.score(eval_target, ear_weights)
                    if s < best_score:
                        bands[weakest_idx] = (f_worst, new_g, new_Q)
                        best_score = s
                        improved = True
                        stuck = 0
                        found = True
                        break
                if found:
                    break
            if not found:
                # Restore original
                cs.update_band(weakest_idx, f0_w, gain_w, Q_w)

        if not improved:
            stuck += 1

    return bands, best_score


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def fit_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
              max_bands: int = 24) -> Tuple[List[BiquadCoeffs], dict]:
    """Fit peaking biquads to target curve using hill climbing.

    Uses greedy seed + random restarts with CascadeState caching
    for fast score evaluation.

    Returns:
        (list_of_BiquadCoeffs, metadata_dict)
    """
    eval_freqs, eval_target, ear_weights = _build_eval_grid(freqs, target_db)

    # Restart count scales inversely with band count
    if max_bands <= 8:
        n_restarts = 4
    elif max_bands <= 16:
        n_restarts = 2
    else:
        n_restarts = 2

    # ── Greedy seed (always included) ──
    greedy_bands = _greedy_bands(freqs, target_db, fs, max_bands)
    best_bands, best_score = _hill_climb(
        greedy_bands, eval_freqs, eval_target, fs, ear_weights)

    # ── Random restarts ──
    for ri in range(n_restarts - 1):
        rbands = _random_bands(max_bands, eval_freqs)
        bands, score = _hill_climb(
            rbands, eval_freqs, eval_target, fs, ear_weights)
        if score < best_score:
            best_score = score
            best_bands = bands

    biquads = _bands_to_biquads(best_bands, fs)
    bands_meta = [
        {"type": "peaking", "f0": b[0], "gain_db": b[1], "Q": b[2]}
        for b in best_bands
    ]

    return biquads, {
        "name": "hillclimb",
        "bands": bands_meta,
        "score": best_score,
    }
