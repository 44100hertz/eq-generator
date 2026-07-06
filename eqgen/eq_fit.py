"""
IIR EQ design: fit a target gain curve to cascaded biquad filters.

Each biquad is a parametric peaking filter, low-shelf, or high-shelf
using the standard Audio EQ Cookbook formulas (RBJ).

The fitting strategy is greedy:
  1. Fit a low-shelf to handle the broad bass boost trend.
  2. Iteratively add peaking filters at the frequency with the largest
     residual error.
  3. Optionally refine all parameters jointly via least-squares.

All coefficients are stored normalized (a0 = 1): [b0, b1, b2, a1, a2].
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional


# ─────────────────────────────────────────────────────────────────────────────
# Biquad design (float precision) — Audio EQ Cookbook
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BiquadCoeffs:
    """Normalized biquad coefficients: b0 + b1·z⁻¹ + b2·z⁻²
                                     ────────────────────────
                                      1 + a1·z⁻¹ + a2·z⁻²
    """
    b0: float
    b1: float
    b2: float
    a1: float
    a2: float

    def to_array(self) -> np.ndarray:
        return np.array([self.b0, self.b1, self.b2, self.a1, self.a2])

    @classmethod
    def from_array(cls, arr: np.ndarray) -> "BiquadCoeffs":
        return cls(arr[0], arr[1], arr[2], arr[3], arr[4])


def design_peaking(f0: float, gain_db: float, Q: float, fs: float) -> BiquadCoeffs:
    """Parametric peaking/cut filter at f0 with gain_db and Q."""
    A = 10.0 ** (gain_db / 40.0)
    omega = 2.0 * np.pi * f0 / fs
    alpha = np.sin(omega) / (2.0 * Q)

    b0 = 1.0 + alpha * A
    b1 = -2.0 * np.cos(omega)
    b2 = 1.0 - alpha * A
    a0 = 1.0 + alpha / A
    a1 = -2.0 * np.cos(omega)
    a2 = 1.0 - alpha / A

    return BiquadCoeffs(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)


def design_low_shelf(f0: float, gain_db: float, Q: float, fs: float) -> BiquadCoeffs:
    """Low-shelf filter at f0 (boost below f0). Q ≈ 0.707 for Butterworth."""
    A = 10.0 ** (gain_db / 40.0)
    omega = 2.0 * np.pi * f0 / fs
    alpha = np.sin(omega) / (2.0 * Q)
    sqrtA = np.sqrt(A)

    b0 = A * ((A + 1.0) - (A - 1.0) * np.cos(omega) + 2.0 * sqrtA * alpha)
    b1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * np.cos(omega))
    b2 = A * ((A + 1.0) - (A - 1.0) * np.cos(omega) - 2.0 * sqrtA * alpha)
    a0 = (A + 1.0) + (A - 1.0) * np.cos(omega) + 2.0 * sqrtA * alpha
    a1 = -2.0 * ((A - 1.0) + (A + 1.0) * np.cos(omega))
    a2 = (A + 1.0) + (A - 1.0) * np.cos(omega) - 2.0 * sqrtA * alpha

    return BiquadCoeffs(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)


def design_butter_hp(fc: float, fs: float) -> BiquadCoeffs:
    """2nd-order Butterworth high-pass (bilinear transform).

    Numerator sums to zero at DC (b0 + b1 + b2 = 0), which survives Q16
    quantization perfectly — no near-cancellation problem.
    """
    SQRT2 = np.sqrt(2.0)
    omega = np.tan(np.pi * fc / fs)
    c = 1.0 + SQRT2 * omega + omega * omega

    b0 = 1.0 / c
    b1 = -2.0 / c
    b2 = 1.0 / c
    a1 = 2.0 * (omega * omega - 1.0) / c
    a2 = (1.0 - SQRT2 * omega + omega * omega) / c
    return BiquadCoeffs(b0, b1, b2, a1, a2)


def design_high_shelf(f0: float, gain_db: float, Q: float, fs: float) -> BiquadCoeffs:
    """High-shelf filter at f0 (boost above f0). Q ≈ 0.707 for Butterworth."""
    A = 10.0 ** (gain_db / 40.0)
    omega = 2.0 * np.pi * f0 / fs
    alpha = np.sin(omega) / (2.0 * Q)
    sqrtA = np.sqrt(A)

    b0 = A * ((A + 1.0) + (A - 1.0) * np.cos(omega) + 2.0 * sqrtA * alpha)
    b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * np.cos(omega))
    b2 = A * ((A + 1.0) + (A - 1.0) * np.cos(omega) - 2.0 * sqrtA * alpha)
    a0 = (A + 1.0) - (A - 1.0) * np.cos(omega) + 2.0 * sqrtA * alpha
    a1 = 2.0 * ((A - 1.0) - (A + 1.0) * np.cos(omega))
    a2 = (A + 1.0) - (A - 1.0) * np.cos(omega) - 2.0 * sqrtA * alpha

    return BiquadCoeffs(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)


def design_highpass(f0: float, Q: float, fs: float) -> BiquadCoeffs:
    """2nd-order high-pass (Butterworth when Q = 0.7071)."""
    omega = 2.0 * np.pi * f0 / fs
    alpha = np.sin(omega) / (2.0 * Q)

    b0 = (1.0 + np.cos(omega)) / 2.0
    b1 = -(1.0 + np.cos(omega))
    b2 = (1.0 + np.cos(omega)) / 2.0
    a0 = 1.0 + alpha
    a1 = -2.0 * np.cos(omega)
    a2 = 1.0 - alpha

    return BiquadCoeffs(b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0)


# ─────────────────────────────────────────────────────────────────────────────
# Cascaded biquad frequency response
# ─────────────────────────────────────────────────────────────────────────────

def biquad_response(coeffs: BiquadCoeffs, freqs: np.ndarray, fs: float) -> np.ndarray:
    """Compute complex frequency response of a single biquad at given frequencies."""
    z = np.exp(-1j * 2.0 * np.pi * freqs / fs)
    H = (coeffs.b0 + coeffs.b1 * z + coeffs.b2 * z**2) / \
        (1.0 + coeffs.a1 * z + coeffs.a2 * z**2)
    return H


def cascade_response(biquads: List[BiquadCoeffs], freqs: np.ndarray, fs: float) -> np.ndarray:
    """Compute magnitude response (linear) of cascaded biquads."""
    H = np.ones(len(freqs), dtype=complex)
    for bc in biquads:
        H *= biquad_response(bc, freqs, fs)
    return np.abs(H)


def cascade_response_db(biquads: List[BiquadCoeffs], freqs: np.ndarray, fs: float) -> np.ndarray:
    """Compute magnitude response in dB of cascaded biquads."""
    mag = cascade_response(biquads, freqs, fs)
    mag = np.maximum(mag, 1e-12)
    return 20.0 * np.log10(mag)


# ─────────────────────────────────────────────────────────────────────────────
# Greedy IIR fitter
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FitResult:
    biquads: List[BiquadCoeffs] = field(default_factory=list)
    bands: List[dict] = field(default_factory=list)  # metadata per band
    residual_db: Optional[np.ndarray] = None

    @property
    def n_bands(self) -> int:
        return len(self.biquads)


def fit_eq_curve(
    freqs: np.ndarray,
    target_db: np.ndarray,
    fs: float,
    max_bands: int = 40,
    min_freq: float = 20.0,
    max_freq: float = 20000.0,
    min_peaking_freq: float = 40.0,
    gain_range: Tuple[float, float] = (-60.0, 60.0),
    q_range: Tuple[float, float] = (0.3, 32.0),
    stop_db: float = 0.3,
) -> FitResult:
    """Fit a target gain curve (dB) to cascaded peaking biquads.

    Greedy per-band (gain, Q) optimization via golden-section search:
      1. Find f0 = frequency with largest |residual|.
      2. Search multiple gain candidates around the residual.
      3. For each gain, golden-section search for optimal Q.
      4. Pick the (gain, Q) that minimizes total RMS error.
      5. Apply the band, update residual, repeat.

    Early-exits when RMS improvement < 0.01 dB and current RMS < 2 × stop_db.
    Multiple bands may land near the same frequency when the target curve
    demands composite shapes (e.g. narrow cut + broad boost).

    Args:
        freqs: frequency points (Hz), sorted ascending.
        target_db: target gain in dB at each frequency.
        fs: sample rate (Hz).
        max_bands: maximum number of biquad bands.
        min_freq, max_freq: frequency range to consider for error.
        min_peaking_freq: minimum centre frequency for peaking filters.
        gain_range: (min_gain_db, max_gain_db) per band.
        q_range: (min_Q, max_Q) search range.
        stop_db: stop when max |residual| falls below this threshold.

    Returns:
        FitResult with biquads, metadata, and final residual.
    """
    result = FitResult()

    # Evaluate on a uniform log-spaced grid so residual search isn't
    # biased by the adaptive measurement grid's point density.
    N_EVAL = 256
    eval_freqs = np.logspace(np.log10(min_freq), np.log10(max_freq), N_EVAL)
    eval_target = np.interp(eval_freqs, freqs, target_db)
    eval_mask = np.ones(N_EVAL, dtype=bool)  # all points already in [min,max]

    cascade_db = np.zeros(N_EVAL)
    residual = eval_target.copy()

    for i in range(max_bands):
        actionable = np.abs(residual)
        too_low = eval_freqs < min_peaking_freq
        actionable[too_low] = 0.0

        if np.max(actionable) < stop_db:
            break

        idx = np.argmax(actionable)
        f0 = eval_freqs[idx]
        error_db = residual[idx]

        best_gain, best_Q, best_err = _search_band_params(
            f0, error_db, cascade_db, eval_target, eval_freqs, fs, eval_mask,
            gain_range=gain_range, q_range=q_range,
        )

        current_rms = float(np.sqrt(np.mean(residual ** 2)))
        if current_rms - best_err < 0.01 and current_rms < stop_db * 2.0:
            break

        pk = design_peaking(f0, best_gain, best_Q, fs)
        result.biquads.append(pk)
        result.bands.append({
            "type": "peaking",
            "f0": f0,
            "gain_db": best_gain,
            "Q": best_Q,
        })

        pk_db = _peaking_response_db(f0, best_gain, best_Q, eval_freqs, fs)
        cascade_db = cascade_db + pk_db
        residual = eval_target - cascade_db

    # Store residual on the original freqs for API compatibility
    result.residual_db = target_db - cascade_response_db(result.biquads, freqs, fs)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: fit to dict-based EQ curve (from model.preprocess_eq_curve)
# ─────────────────────────────────────────────────────────────────────────────

def fit_from_dict(
    eq_curve: Dict[float, float],  # {freq_hz: gain_linear}
    fs: float,
    max_bands: int = 40,
    **kwargs,
) -> FitResult:
    """Fit a dict-based EQ curve (from model) to cascaded biquads."""
    freqs = np.array(sorted(eq_curve.keys()))
    gains_linear = np.array([eq_curve[f] for f in freqs])
    target_db = 20.0 * np.log10(np.maximum(gains_linear, 1e-12))
    return fit_eq_curve(freqs, target_db, fs, max_bands=max_bands, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Per-band (gain, Q) search helpers
# ─────────────────────────────────────────────────────────────────────────────

# Golden ratio
_PHI = (1.0 + np.sqrt(5.0)) / 2.0


def _peaking_response_db(f0: float, gain_db: float, Q: float,
                          freqs: np.ndarray, fs: float) -> np.ndarray:
    """Magnitude response in dB of a single peaking filter."""
    pk = design_peaking(f0, gain_db, Q, fs)
    z = np.exp(-1j * 2.0 * np.pi * freqs / fs)
    H = (pk.b0 + pk.b1 * z + pk.b2 * z**2) / (1.0 + pk.a1 * z + pk.a2 * z**2)
    mag = np.abs(H)
    mag = np.maximum(mag, 1e-12)
    return 20.0 * np.log10(mag)


def _rms_error_db(cascade_db: np.ndarray, peaking_db: np.ndarray,
                  target_db: np.ndarray, mask: np.ndarray) -> float:
    """RMS error after adding a peaking filter to the current cascade."""
    new_db = cascade_db + peaking_db
    err = new_db - target_db
    return float(np.sqrt(np.mean(err[mask] ** 2)))


def _golden_search_Q(f0: float, gain_db: float, cascade_db: np.ndarray,
                     target_db: np.ndarray, freqs: np.ndarray, fs: float,
                     mask: np.ndarray,
                     q_lo: float = 0.3, q_hi: float = 6.0,
                     tol: float = 0.05, max_iter: int = 15) -> float:
    """Golden-section search for optimal Q at given (f0, gain_db).

    Assumes the RMS error as a function of Q is unimodal.
    Returns the Q that minimizes RMS error.
    """
    # Evaluate a few points to find a rough bracket
    test_Qs = np.array([q_lo, 1.0, 2.0, 4.0, q_hi])
    test_Qs = test_Qs[(test_Qs >= q_lo) & (test_Qs <= q_hi)]
    if len(test_Qs) < 3:
        test_Qs = np.linspace(q_lo, q_hi, 5)

    test_errs = np.array([
        _rms_error_db(cascade_db,
                      _peaking_response_db(f0, gain_db, q, freqs, fs),
                      target_db, mask)
        for q in test_Qs
    ])

    best_idx = int(np.argmin(test_errs))

    if best_idx == 0:
        a, b = q_lo, test_Qs[1] if len(test_Qs) > 1 else (q_lo + q_hi) / 2.0
    elif best_idx == len(test_Qs) - 1:
        a, b = test_Qs[-2] if len(test_Qs) > 1 else (q_lo + q_hi) / 2.0, q_hi
    else:
        a, b = test_Qs[best_idx - 1], test_Qs[best_idx + 1]

    # Golden-section within [a, b]
    c = b - (b - a) / _PHI
    d = a + (b - a) / _PHI

    ec = _rms_error_db(cascade_db,
                       _peaking_response_db(f0, gain_db, c, freqs, fs),
                       target_db, mask)
    ed = _rms_error_db(cascade_db,
                       _peaking_response_db(f0, gain_db, d, freqs, fs),
                       target_db, mask)

    for _ in range(max_iter):
        if b - a < tol:
            break
        if ec < ed:
            b = d
            d = c
            ed = ec
            c = b - (b - a) / _PHI
            ec = _rms_error_db(cascade_db,
                               _peaking_response_db(f0, gain_db, c, freqs, fs),
                               target_db, mask)
        else:
            a = c
            c = d
            ec = ed
            d = a + (b - a) / _PHI
            ed = _rms_error_db(cascade_db,
                               _peaking_response_db(f0, gain_db, d, freqs, fs),
                               target_db, mask)

    return (a + b) / 2.0


def _search_band_params(f0: float, residual_at_f0: float,
                         cascade_db: np.ndarray, target_db: np.ndarray,
                         freqs: np.ndarray, fs: float, mask: np.ndarray,
                         gain_range: Tuple[float, float] = (-60.0, 60.0),
                         q_range: Tuple[float, float] = (0.3, 32.0),
                         ) -> Tuple[float, float, float]:
    """Search for optimal (gain_db, Q) for a peaking filter at f0.

    Tries gain values around residual_at_f0, and for each does golden-section
    search on Q. Returns (best_gain, best_Q, best_rms_error).
    """
    # Candidate gains: residual scaled by various factors
    multipliers = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    gain_candidates = np.clip(multipliers * residual_at_f0, gain_range[0], gain_range[1])
    gain_candidates = np.unique(np.round(gain_candidates, decimals=1))

    best_gain = gain_candidates[0]
    best_Q = 2.0
    best_err = float('inf')

    for g in gain_candidates:
        if abs(g) < 0.5 and abs(residual_at_f0) > 2.0:
            continue  # skip gains too weak relative to error

        q_opt = _golden_search_Q(f0, float(g), cascade_db, target_db,
                                  freqs, fs, mask,
                                  q_lo=q_range[0], q_hi=q_range[1])
        pk_db = _peaking_response_db(f0, float(g), q_opt, freqs, fs)
        err = _rms_error_db(cascade_db, pk_db, target_db, mask)

        if err < best_err:
            best_err = err
            best_gain = float(g)
            best_Q = q_opt

    return best_gain, best_Q, best_err



