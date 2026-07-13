"""
CMA-ES IIR biquad fitter.

Covariance Matrix Adaptation Evolution Strategy — from scratch, no external deps.
Uses multi-start (random + greedy seeding) to escape local optima.
Optimizes [log10(f0), gain, Q] × max_bands with perceptually-weighted
RMS error (ISO 226:2023).
"""

import numpy as np
from typing import List, Tuple

from eqgen.eq_fit import BiquadCoeffs, design_peaking
from eqgen.model import ear_sensitivity

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

_GAIN_MIN = -24.0
_GAIN_MAX = 24.0
_Q_MIN = 0.3
_Q_MAX = 6.0


# ─────────────────────────────────────────────────────────────────────────────
# Objective function (permutation-invariant via f0-sorting)
# ─────────────────────────────────────────────────────────────────────────────

def _make_objective(freqs, target_db, fs, max_bands, percept_weights):
    """Build a fast objective function. Returns f(x) → perceptual RMS."""
    z = np.exp(-1j * 2.0 * np.pi * freqs / fs)
    n_freqs = len(freqs)

    def obj(x: np.ndarray) -> float:
        bands = x.reshape(max_bands, 3)
        order = np.argsort(bands[:, 0])
        sb = bands[order]

        H = np.ones(n_freqs, dtype=complex)
        for i in range(max_bands):
            log10_f0, gain_db, Q_val = sb[i, 0], sb[i, 1], sb[i, 2]
            if abs(gain_db) < 0.05:
                continue
            f0 = 10.0 ** log10_f0
            pk = design_peaking(f0, gain_db, Q_val, fs)
            H *= (pk.b0 + pk.b1 * z + pk.b2 * z**2) / \
                 (1.0 + pk.a1 * z + pk.a2 * z**2)

        mag = np.abs(H)
        mag = np.maximum(mag, 1e-12)
        fitted_db = 20.0 * np.log10(mag)
        err = (fitted_db - target_db) * percept_weights
        return float(np.sqrt(np.mean(err ** 2)))

    return obj


# ─────────────────────────────────────────────────────────────────────────────
# Initialization
# ─────────────────────────────────────────────────────────────────────────────

def _greedy_x0(freqs, target_db, fs, max_bands):
    """Get greedy solution as a flat parameter vector."""
    from eqgen.eq_fit import fit_eq_curve

    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=freqs[0])

    f0_log_min = np.log10(max(freqs[0], 20.0))
    f0_log_max = np.log10(min(freqs[-1], 14000.0))

    x = np.zeros(max_bands * 3)
    for i, b in enumerate(fit.bands):
        if i >= max_bands:
            break
        x[i * 3]     = np.log10(b["f0"])
        x[i * 3 + 1] = b["gain_db"]
        x[i * 3 + 2] = b["Q"]

    for i in range(len(fit.bands), max_bands):
        x[i * 3]     = f0_log_min + (f0_log_max - f0_log_min) * i / max_bands
        x[i * 3 + 1] = 0.0
        x[i * 3 + 2] = 1.5
    return x


def _random_x0(freqs, target_db, max_bands):
    """Random initial parameter vector, guided by the target curve.

    Places bands near frequencies where target_db deviates from zero.
    """
    rng = np.random.default_rng()
    f0_log_min = np.log10(max(freqs[0], 20.0))
    f0_log_max = np.log10(min(freqs[-1], 14000.0))

    # Use the absolute target_db as a "weight" for where to place bands
    abs_target = np.abs(target_db)

    x = np.zeros(max_bands * 3)
    for i in range(max_bands):
        # Weighted random selection of frequency
        idx = rng.choice(len(freqs), p=abs_target / (abs_target.sum() + 1e-10))
        log10_f0 = np.log10(max(freqs[idx], 10.0 ** f0_log_min))

        # Jitter a bit
        log10_f0 += rng.uniform(-0.15, 0.15)
        log10_f0 = np.clip(log10_f0, f0_log_min, f0_log_max)

        # Gain: roughly target at that frequency, with some scatter
        gain = target_db[idx] * rng.uniform(0.3, 1.2) / max_bands
        gain = np.clip(gain, _GAIN_MIN, _GAIN_MAX)

        Q = rng.uniform(0.5, 4.0)

        x[i * 3]     = log10_f0
        x[i * 3 + 1] = gain
        x[i * 3 + 2] = Q

    return x


# ─────────────────────────────────────────────────────────────────────────────
# CMA-ES single run
# ─────────────────────────────────────────────────────────────────────────────

def _cmaes_run(x0, bounds_low, bounds_high, obj_fn,
               max_gens=150, sigma0=0.5, seed=None):
    """Run one CMA-ES optimization from x0. Returns (best_x, best_f, n_evals)."""
    if seed is not None:
        np.random.seed(seed)

    N = len(x0)
    lamb = max(6, 4 + int(3.0 * np.log(N)))
    mu = lamb // 2
    w_raw = np.log(mu + 0.5) - np.log(np.arange(1, mu + 1))
    weights = w_raw / w_raw.sum()
    mueff = 1.0 / (weights ** 2).sum()

    cc = (4.0 + mueff / N) / (N + 4.0 + 2.0 * mueff / N)
    cs = (mueff + 2.0) / (N + mueff + 5.0)
    c1 = 2.0 / ((N + 1.3) ** 2 + mueff)
    cmu = min(1.0 - c1, 2.0 * (mueff - 2.0 + 1.0 / mueff) / ((N + 2.0) ** 2 + mueff))
    damps = 1.0 + 2.0 * max(0.0, np.sqrt((mueff - 1.0) / (N + 1.0)) - 1.0) + cs
    chiN = np.sqrt(N) * (1.0 - 1.0 / (4.0 * N) + 1.0 / (21.0 * N * N))

    xmean = x0.copy()
    sigma = sigma0
    C = np.eye(N)
    B = np.eye(N)
    D = np.ones(N)
    pc = np.zeros(N)
    ps = np.zeros(N)

    best_x = x0.copy()
    best_f = obj_fn(x0)
    n_evals = 1
    stale = 0
    restart_count = 0

    for gen in range(1, max_gens + 1):
        # Sample + evaluate
        zs = np.random.randn(lamb, N)
        arx = np.zeros((lamb, N))
        for k in range(lamb):
            arx[k] = np.clip(xmean + sigma * np.dot(B, D * zs[k]), bounds_low, bounds_high)
        arf = np.array([obj_fn(arx[k]) for k in range(lamb)])
        n_evals += lamb

        idx = np.argsort(arf)
        if arf[idx[0]] < best_f - 1e-10:
            best_f = arf[idx[0]]
            best_x = arx[idx[0]].copy()
            stale = 0
        else:
            stale += 1

        # Weighted recombination
        xold = xmean.copy()
        xmean = np.sum(arx[idx[:mu]] * weights.reshape(-1, 1), axis=0)

        # Cumulation
        zmean = np.sum(zs[idx[:mu]] * weights.reshape(-1, 1), axis=0)
        ps = (1.0 - cs) * ps + np.sqrt(cs * (2.0 - cs) * mueff) * np.dot(B, zmean)
        hsig = float(
            np.linalg.norm(ps) / np.sqrt(1.0 - (1.0 - cs) ** (2 * gen)) / chiN
            < 1.4 + 2.0 / (N + 1.0)
        )
        pc = (1.0 - cc) * pc + hsig * np.sqrt(cc * (2.0 - cc) * mueff) * \
             np.dot(B, D * zmean)

        # Covariance update (rank-mu + rank-one)
        artmp = np.dot(B, D[:, np.newaxis] * zs[idx[:mu]].T)
        C = (1.0 - c1 - cmu) * C \
            + c1 * (np.outer(pc, pc) + (1.0 - hsig) * cc * (2.0 - cc) * C) \
            + cmu * np.dot(artmp * weights, artmp.T)

        # Step-size adaptation
        sigma = sigma * np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1.0))

        # Eigendecomposition
        try:
            D2, B = np.linalg.eigh(C)
            D2 = np.maximum(D2, 1e-30)
            D = np.sqrt(D2)
        except np.linalg.LinAlgError:
            C = np.eye(N)
            B = np.eye(N)
            D = np.ones(N)

        # Restart if stalled
        if stale > 60:
            if restart_count < 2:
                sigma = 1.0
                pc[:] = 0
                ps[:] = 0
                C = np.eye(N)
                B = np.eye(N)
                D = np.ones(N)
                stale = 0
                restart_count += 1
            else:
                break

        if sigma < 1e-8 and stale > 15:
            break

    return best_x, best_f, n_evals


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def fit_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
              max_bands: int = 24) -> Tuple[List[BiquadCoeffs], dict]:
    """CMA-ES multi-start fitter.

    Strategy:
      1. Run CMA-ES from greedy seed (refine the local optimum).
      2. Run CMA-ES from 3-5 random seeds (explore other basins).
      3. Return the best found.
    """
    N = max_bands * 3
    if N == 0:
        return [], {"name": "cma_es", "bands": []}

    percept_weights = np.array([ear_sensitivity(float(f)) for f in freqs])
    obj_fn = _make_objective(freqs, target_db, fs, max_bands, percept_weights)

    f0_log_min = np.log10(max(freqs[0], 20.0))
    f0_log_max = np.log10(min(freqs[-1], 14000.0))
    bounds_low = np.empty(N)
    bounds_high = np.empty(N)
    for i in range(max_bands):
        bounds_low[i * 3]     = f0_log_min
        bounds_low[i * 3 + 1] = _GAIN_MIN
        bounds_low[i * 3 + 2] = _Q_MIN
        bounds_high[i * 3]     = f0_log_max
        bounds_high[i * 3 + 1] = _GAIN_MAX
        bounds_high[i * 3 + 2] = _Q_MAX

    # Number of random restarts scales inversely with problem size
    n_random = max(1, min(5, 120 // N))
    gens_per = max(80, min(250, 500 // (N // 12)))

    best_x = None
    best_f = float('inf')
    total_evals = 0

    # Run greedily-seeded
    x0 = _greedy_x0(freqs, target_db, fs, max_bands)
    x, f, evals = _cmaes_run(x0, bounds_low, bounds_high, obj_fn,
                              max_gens=gens_per, sigma0=0.3, seed=42)
    total_evals += evals
    if f < best_f:
        best_f, best_x = f, x

    # Run randomly-seeded
    for r in range(n_random):
        x0 = _random_x0(freqs, target_db, max_bands)
        x, f, evals = _cmaes_run(x0, bounds_low, bounds_high, obj_fn,
                                  max_gens=gens_per, sigma0=0.5, seed=100 + r)
        total_evals += evals
        if f < best_f:
            best_f, best_x = f, x

    # Build biquads from best solution
    bands_3 = best_x.reshape(max_bands, 3)
    order = np.argsort(bands_3[:, 0])
    sorted_bands = bands_3[order]

    biquads = []
    bands = []
    for i in range(max_bands):
        log10_f0, gain_db, Q_val = sorted_bands[i]
        if abs(gain_db) < 0.05:
            continue
        f0 = float(10.0 ** log10_f0)
        pk = design_peaking(f0, float(gain_db), float(Q_val), fs)
        biquads.append(pk)
        bands.append({
            "type": "peaking",
            "f0": f0,
            "gain_db": float(gain_db),
            "Q": float(Q_val),
        })

    return biquads, {"name": "cma_es", "bands": bands}
