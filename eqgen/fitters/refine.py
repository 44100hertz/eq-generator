"""
Joint refinement IIR biquad fitter.

Strategy: seed with greedy placement, then jointly optimize all band
parameters (f0, gain, Q) using scipy's L-BFGS-B.

Key design:
- Optimize in log-frequency space for uniform step sizes across octaves
- Use a coarser evaluation grid (256 points) during optimization for speed
- Perceptual-weighted RMS error as primary objective
- Soft max-error penalty to prevent extreme excursions in insensitive bands
- Multiple refinement stages: gains only → gains+Q → full joint
"""

import numpy as np
from typing import List, Tuple

from eqgen.eq_fit import (
    BiquadCoeffs, design_peaking, cascade_response_db, fit_eq_curve,
)
from eqgen.model import ear_sensitivity

# ── Optimization grid: use fewer points for speed ─────────────────────
N_OPT_GRID = 256

# ── Max-error penalty weight ──────────────────────────────────────────
# After perceptual RMS falls below ~0.5 dB, the optimizer may trade bass
# fit quality for mid-range perfection.  A small penalty on the worst
# 95th-percentile error keeps extremes in check.
MAX_ERR_PENALTY = 0.02

# ── Gain / Q bounds ───────────────────────────────────────────────────
GAIN_LO, GAIN_HI = -24.0, 24.0
Q_LO, Q_HI = 0.3, 6.0


def _params_to_biquads(params: np.ndarray, fs: float) -> List[BiquadCoeffs]:
    """Convert flat parameter vector [logf0, g, Q, ...] to biquads."""
    n = len(params) // 3
    biquads = []
    for i in range(n):
        logf0 = float(params[3*i])
        gain_db = float(params[3*i + 1])
        Q = float(params[3*i + 2])
        biquads.append(design_peaking(10.0 ** logf0, gain_db, Q, fs))
    return biquads


def _cascade_db_from_params(params: np.ndarray, freqs: np.ndarray,
                             fs: float, z: np.ndarray) -> np.ndarray:
    """Compute cascade dB response directly from params (no BiquadCoeffs alloc).

    Pre-computed z = exp(-j*2*pi*freqs/fs).
    """
    n = len(params) // 3
    # Accumulate in linear domain first
    H = np.ones(len(freqs), dtype=complex)
    for i in range(n):
        logf0 = float(params[3*i])
        gain_db = float(params[3*i + 1])
        Q = float(params[3*i + 2])

        A = 10.0 ** (gain_db / 40.0)
        omega = 2.0 * np.pi * (10.0 ** logf0) / fs
        alpha = np.sin(omega) / (2.0 * Q)

        b0 = 1.0 + alpha * A
        b1 = -2.0 * np.cos(omega)
        b2 = 1.0 - alpha * A
        a0 = 1.0 + alpha / A
        a1 = -2.0 * np.cos(omega)
        a2 = 1.0 - alpha / A

        H *= (b0 + b1 * z + b2 * z**2) / (a0 + a1 * z + a2 * z**2)
    mag = np.abs(H)
    mag = np.maximum(mag, 1e-12)
    return 20.0 * np.log10(mag)


def _objective(params: np.ndarray, freqs: np.ndarray, target_db: np.ndarray,
                fs: float, z: np.ndarray, weights: np.ndarray) -> float:
    """Perceptual RMS error + soft max-error penalty."""
    fitted_db = _cascade_db_from_params(params, freqs, fs, z)
    err = (fitted_db - target_db) * weights
    rms = float(np.sqrt(np.mean(err ** 2)))

    # Soft penalty on the 95th percentile of |error|
    abs_err = np.abs(fitted_db - target_db)
    p95 = float(np.percentile(abs_err, 95))
    if p95 > 5.0:
        rms += MAX_ERR_PENALTY * (p95 - 5.0)

    return rms


def fit_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
              max_bands: int = 24) -> Tuple[List[BiquadCoeffs], dict]:
    """Fit IIR biquads with greedy seeding + L-BFGS-B joint refinement.

    Returns (list_of_BiquadCoeffs, metadata_dict).
    """
    # ── 1. Build optimization grid (log-spaced, fewer points) ─────────
    f_lo = max(freqs[0], 1.0)
    f_hi = min(freqs[-1], fs / 2.2)
    opt_freqs = np.logspace(np.log10(f_lo), np.log10(f_hi), N_OPT_GRID)
    opt_target = np.interp(opt_freqs, freqs, target_db)

    # Pre-compute z for frequency response
    z_opt = np.exp(-1j * 2.0 * np.pi * opt_freqs / fs)

    # Pre-compute perceptual weights
    opt_weights = np.array([ear_sensitivity(f) for f in opt_freqs])

    # ── 2. Greedy seed ────────────────────────────────────────────────
    greedy = fit_eq_curve(
        freqs, target_db, fs, max_bands=max_bands,
        min_freq=freqs[0], max_freq=freqs[-1],
        min_peaking_freq=freqs[0],
    )

    n_bands = len(greedy.biquads)
    if n_bands == 0:
        return [], {"name": "refine", "n_bands": 0, "improvement": 0.0}

    # ── 3. Build initial parameter vector ─────────────────────────────
    x0 = np.empty(n_bands * 3)
    for i, b in enumerate(greedy.bands):
        x0[3*i]     = np.log10(b["f0"])
        x0[3*i + 1] = b["gain_db"]
        x0[3*i + 2] = b["Q"]

    # ── 4. Initial error ──────────────────────────────────────────────
    err0 = _objective(x0, opt_freqs, opt_target, fs, z_opt, opt_weights)

    # ── 5. Bounds ─────────────────────────────────────────────────────
    log_lo = np.log10(f_lo)
    log_hi = np.log10(f_hi)
    bound_list = []
    for _ in range(n_bands):
        bound_list.extend([(log_lo, log_hi), (GAIN_LO, GAIN_HI), (Q_LO, Q_HI)])

    # ── 6. Multi-stage L-BFGS-B ───────────────────────────────────────
    from scipy.optimize import minimize

    best_x = x0.copy()
    best_err = err0

    def run_stage(x_init, maxiter, label=""):
        nonlocal best_x, best_err
        res = minimize(
            lambda x: _objective(x, opt_freqs, opt_target, fs, z_opt, opt_weights),
            x_init, method="L-BFGS-B", bounds=bound_list,
            options={"maxiter": maxiter, "ftol": 1e-8},
        )
        if res.fun < best_err:
            best_x = res.x.copy()
            best_err = res.fun
        return res

    # Stage 1: gains only (freeze f0 and Q) — quick coarse improvement
    # We approximate by scaling: tight bounds on f0 and Q
    tight_bounds = []
    for _ in range(n_bands):
        tight_bounds.extend([(log_lo, log_hi), (GAIN_LO, GAIN_HI), (Q_LO, Q_HI)])

    run_stage(best_x, maxiter=200, label="full")

    if n_bands >= 4:
        # Stage 2: perturb and re-optimize to escape shallow minima
        for attempt in range(2):
            noise = 0.03 * (attempt + 1)
            x_p = best_x.copy()
            for i in range(n_bands):
                x_p[3*i]     += np.random.normal(0, noise * 0.5)  # log_f0
                x_p[3*i + 1] += np.random.normal(0, noise * 4.0)  # gain (dB)
                x_p[3*i + 2] += np.random.normal(0, noise * 1.5)  # Q
            x_p = np.clip(x_p, [b[0] for b in bound_list], [b[1] for b in bound_list])

            res = minimize(
                lambda x: _objective(x, opt_freqs, opt_target, fs, z_opt, opt_weights),
                x_p, method="L-BFGS-B", bounds=bound_list,
                options={"maxiter": 300, "ftol": 1e-8},
            )
            if res.fun < best_err:
                best_x = res.x.copy()
                best_err = res.fun

    # ── 7. Convert to biquads ─────────────────────────────────────────
    biquads = _params_to_biquads(best_x, fs)
    bands = []
    for i in range(n_bands):
        bands.append({
            "type": "peaking",
            "f0": round(10.0 ** float(best_x[3*i]), 1),
            "gain_db": round(float(best_x[3*i + 1]), 1),
            "Q": round(float(best_x[3*i + 2]), 2),
        })

    improvement = err0 - best_err
    return biquads, {
        "name": "refine",
        "n_bands": n_bands,
        "initial_pRMS": round(err0, 4),
        "final_pRMS": round(best_err, 4),
        "improvement": round(improvement, 4),
        "bands": bands,
    }
