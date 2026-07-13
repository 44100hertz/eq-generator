"""
Balanced joint refinement IIR biquad fitter.

Like refine.py, but uses a combined objective (perceptual + unweighted)
to prevent the optimizer from sacrificing bass/treble accuracy.

Objective: 0.5 * pRMS + 0.5 * uRMS
"""

import numpy as np
from typing import List, Tuple

from eqgen.eq_fit import (
    BiquadCoeffs, design_peaking, fit_eq_curve,
)
from eqgen.model import ear_sensitivity

N_OPT_GRID = 256
GAIN_LO, GAIN_HI = -24.0, 24.0
Q_LO, Q_HI = 0.3, 6.0


def _cascade_db_from_params(params: np.ndarray, freqs: np.ndarray,
                             fs: float, z: np.ndarray) -> np.ndarray:
    n = len(params) // 3
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
    """50/50 blend of perceptual and unweighted RMS error."""
    fitted_db = _cascade_db_from_params(params, freqs, fs, z)
    err = fitted_db - target_db
    p_rms = float(np.sqrt(np.mean((err * weights) ** 2)))
    u_rms = float(np.sqrt(np.mean(err ** 2)))
    return 0.5 * p_rms + 0.5 * u_rms


def fit_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
              max_bands: int = 24) -> Tuple[List[BiquadCoeffs], dict]:
    f_lo = max(freqs[0], 1.0)
    f_hi = min(freqs[-1], fs / 2.2)
    opt_freqs = np.logspace(np.log10(f_lo), np.log10(f_hi), N_OPT_GRID)
    opt_target = np.interp(opt_freqs, freqs, target_db)

    z_opt = np.exp(-1j * 2.0 * np.pi * opt_freqs / fs)
    opt_weights = np.array([ear_sensitivity(f) for f in opt_freqs])

    greedy = fit_eq_curve(
        freqs, target_db, fs, max_bands=max_bands,
        min_freq=freqs[0], max_freq=freqs[-1],
        min_peaking_freq=freqs[0],
    )

    n_bands = len(greedy.biquads)
    if n_bands == 0:
        return [], {"name": "blend", "n_bands": 0}

    x0 = np.empty(n_bands * 3)
    for i, b in enumerate(greedy.bands):
        x0[3*i]     = np.log10(b["f0"])
        x0[3*i + 1] = b["gain_db"]
        x0[3*i + 2] = b["Q"]

    log_lo, log_hi = np.log10(f_lo), np.log10(f_hi)
    bound_list = [(log_lo, log_hi), (GAIN_LO, GAIN_HI), (Q_LO, Q_HI)] * n_bands

    from scipy.optimize import minimize

    res = minimize(
        lambda x: _objective(x, opt_freqs, opt_target, fs, z_opt, opt_weights),
        x0, method="L-BFGS-B", bounds=bound_list,
        options={"maxiter": 500, "ftol": 1e-8},
    )
    x_best = res.x

    # Try a perturbed restart
    x_p = x_best.copy()
    for i in range(n_bands):
        x_p[3*i]     += np.random.normal(0, 0.02)
        x_p[3*i + 1] += np.random.normal(0, 0.5)
        x_p[3*i + 2] += np.random.normal(0, 0.2)
    bounds_arr = np.array(bound_list)
    x_p = np.clip(x_p, bounds_arr[:, 0], bounds_arr[:, 1])

    res2 = minimize(
        lambda x: _objective(x, opt_freqs, opt_target, fs, z_opt, opt_weights),
        x_p, method="L-BFGS-B", bounds=bound_list,
        options={"maxiter": 500, "ftol": 1e-8},
    )
    if res2.fun < res.fun:
        x_best = res2.x

    biquads = []
    bands = []
    for i in range(n_bands):
        f0 = 10.0 ** float(x_best[3*i])
        gain_db = float(x_best[3*i + 1])
        Q = float(x_best[3*i + 2])
        if abs(gain_db) >= 0.05:
            biquads.append(design_peaking(f0, gain_db, Q, fs))
            bands.append({"type": "peaking", "f0": f0, "gain_db": gain_db, "Q": Q})

    return biquads, {"name": "blend", "n_bands": n_bands, "bands": bands}
