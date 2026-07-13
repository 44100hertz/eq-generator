#!/usr/bin/env python3
"""Differential Evolution IIR biquad fitter.

Hybrid approach: seed the population from the greedy fitter's solution,
then use DE + L-BFGS-B to jointly refine all bands.

Strategy:
- Phase 0: Run greedy to get a good starting point
- Phase 1: Populate DE individuals as perturbations of greedy solution
- Phase 2: Run DE with best/1/bin mutation for fast convergence
- Phase 3: Final L-BFGS-B polish

Key insight: DE from random start in 24+ dimensions is too slow.
Starting near the greedy optimum lets DE escape local optima the greedy
algorithm missed, without needing to explore the full search space.
"""

import numpy as np
from typing import List, Tuple
from dataclasses import dataclass, field

from eqgen.eq_fit import (
    BiquadCoeffs, design_peaking, cascade_response_db,
    fit_eq_curve,
)
from eqgen.model import ear_sensitivity

# ── DE parameters ─────────────────────────────────────────────────────
POP_MULTIPLIER = 5       # NP = POP_MULTIPLIER * D
MAX_GENERATIONS = 80
EARLY_STOP_GENS = 20     # stop if no improvement for this many gens
EARLY_STOP_TOL = 1e-6    # relative improvement threshold
POLISH_GENS = 12         # run L-BFGS-B on best every POLISH_GENS
SIGMA_INIT = 0.15        # perturbation std dev in parameter space

# jDE adaptation rates
TAU1 = 0.1               # prob of updating F_i
TAU2 = 0.1               # prob of updating CR_i
F_LOWER = 0.1
F_UPPER = 0.8

# Bounds in parameter space
GAIN_MIN, GAIN_MAX = -24.0, 24.0
Q_MIN, Q_MAX = 0.3, 6.0


# ── Fast vectorized response ──────────────────────────────────────────

def _cascade_response_db_fast(biquads: List[BiquadCoeffs],
                               z: np.ndarray) -> np.ndarray:
    """Compute magnitude response in dB for cascaded biquads.

    z = exp(-1j * 2*pi * freqs / fs) — precomputed complex exponentials.
    """
    H = np.ones(len(z), dtype=complex)
    for bc in biquads:
        num = bc.b0 + bc.b1 * z + bc.b2 * z * z
        den = 1.0 + bc.a1 * z + bc.a2 * z * z
        H *= num / den
    mag = np.abs(H)
    mag = np.maximum(mag, 1e-12)
    return 20.0 * np.log10(mag)


# ── Objective ─────────────────────────────────────────────────────────

def _params_to_biquads(params: np.ndarray, fs: float,
                        n_bands: int) -> List[BiquadCoeffs]:
    """Convert flat parameter vector to list of BiquadCoeffs (sorted by f0)."""
    band_data = []
    for i in range(n_bands):
        log_f0 = params[i * 3]
        gain = params[i * 3 + 1]
        q_val = params[i * 3 + 2]
        f0 = 10.0 ** log_f0
        band_data.append((log_f0, f0, gain, q_val))

    # Sort by f0 to break permutation symmetry
    band_data.sort(key=lambda x: x[0])

    biquads = []
    for _, f0, gain, q_val in band_data:
        biquads.append(design_peaking(f0, gain, q_val, fs))
    return biquads


def _error(params: np.ndarray, z: np.ndarray, target_db: np.ndarray,
           weights: np.ndarray, fs: float, n_bands: int) -> float:
    """Perceptually-weighted RMS error."""
    biquads = _params_to_biquads(params, fs, n_bands)
    fitted_db = _cascade_response_db_fast(biquads, z)
    err = fitted_db - target_db
    return float(np.sqrt(np.mean((err * weights) ** 2)))


# ── Greedy → parameter vector ────────────────────────────────────────

def _greedy_to_params(freqs: np.ndarray, target_db: np.ndarray, fs: float,
                       max_bands: int) -> np.ndarray:
    """Run the greedy fitter and convert result to DE parameter vector.

    Returns flat array [log10(f0), gain, Q] * max_bands, sorted by f0.
    If greedy produces fewer than max_bands, pad with flat bands.
    """
    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=freqs[0],
                       gain_range=(GAIN_MIN, GAIN_MAX),
                       q_range=(Q_MIN, Q_MAX))

    params = np.zeros(max_bands * 3)
    for i, b in enumerate(fit.bands):
        if i >= max_bands:
            break
        params[i * 3] = np.log10(b['f0'])
        params[i * 3 + 1] = b['gain_db']
        params[i * 3 + 2] = b['Q']

    # Sort by f0
    band_data = [(params[i*3], params[i*3+1], params[i*3+2]) for i in range(max_bands)]
    band_data.sort(key=lambda x: x[0])
    for i, (lf0, g, q) in enumerate(band_data):
        params[i*3] = lf0
        params[i*3+1] = g
        params[i*3+2] = q

    return params


# ── Bounds helpers ────────────────────────────────────────────────────

def _make_bounds(min_freq: float, max_freq: float, n_bands: int) -> list:
    """Create bounds list for L-BFGS-B."""
    log_f0_min = np.log10(min_freq)
    log_f0_max = np.log10(max_freq)
    bounds = []
    for _ in range(n_bands):
        bounds.append((log_f0_min, log_f0_max))
        bounds.append((GAIN_MIN, GAIN_MAX))
        bounds.append((Q_MIN, Q_MAX))
    return bounds


def _filled_population(NP: int, seed_params: np.ndarray,
                        min_freq: float, max_freq: float,
                        n_bands: int, sigma: float) -> np.ndarray:
    """Create population by perturbing seed_params with Gaussian noise.

    First individual is the seed itself; rest are seed + N(0, sigma).
    Sigma is in parameter space units (log for f0, dB for gain, linear for Q).
    """
    D = n_bands * 3
    pop = np.tile(seed_params, (NP, 1))

    log_f0_min = np.log10(min_freq)
    log_f0_max = np.log10(max_freq)

    for j in range(1, NP):
        noise = np.random.normal(0, sigma, D)
        pop[j] = pop[j] + noise

    # Clip all to bounds
    for k in range(n_bands):
        i = k * 3
        pop[:, i] = np.clip(pop[:, i], log_f0_min, log_f0_max)
        pop[:, i + 1] = np.clip(pop[:, i + 1], GAIN_MIN, GAIN_MAX)
        pop[:, i + 2] = np.clip(pop[:, i + 2], Q_MIN, Q_MAX)

    return pop


# ── jDE mutation / crossover (best/1/bin) ─────────────────────────────

def _best1_bin(pop: np.ndarray, fitness: np.ndarray,
               F_vals: np.ndarray, CR_vals: np.ndarray,
               best_idx: int, min_freq: float, max_freq: float,
               n_bands: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """DE/best/1/bin with jDE self-adaptation."""
    NP, D = pop.shape
    trial = np.empty((NP, D))
    new_F = F_vals.copy()
    new_CR = CR_vals.copy()

    log_f0_min = np.log10(min_freq)
    log_f0_max = np.log10(max_freq)
    best = pop[best_idx]

    for i in range(NP):
        if np.random.random() < TAU1:
            new_F[i] = F_LOWER + np.random.random() * (F_UPPER - F_LOWER)
        if np.random.random() < TAU2:
            new_CR[i] = np.random.random()

        Fi = new_F[i]
        CRi = new_CR[i]

        # Pick 2 distinct random individuals, both ≠ i and ≠ best_idx
        candidates = [j for j in range(NP) if j != i]
        r1, r2 = np.random.choice(candidates, 2, replace=False)

        # Mutation: v = best + F * (x_r1 - x_r2)
        donor = best + Fi * (pop[r1] - pop[r2])

        # Crossover: binomial
        j_rand = np.random.randint(D)
        trial_i = np.where(
            (np.random.random(D) <= CRi) | (np.arange(D) == j_rand),
            donor, pop[i]
        )

        # Clip to bounds
        for k in range(n_bands):
            j = k * 3
            trial_i[j] = np.clip(trial_i[j], log_f0_min, log_f0_max)
            trial_i[j + 1] = np.clip(trial_i[j + 1], GAIN_MIN, GAIN_MAX)
            trial_i[j + 2] = np.clip(trial_i[j + 2], Q_MIN, Q_MAX)

        trial[i] = trial_i

    return trial, new_F, new_CR


# ── L-BFGS-B polish ───────────────────────────────────────────────────

def _polish_best(params: np.ndarray, z: np.ndarray, target_db: np.ndarray,
                  weights: np.ndarray, fs: float, n_bands: int,
                  min_freq: float, max_freq: float) -> np.ndarray:
    """Run L-BFGS-B on the best candidate for local refinement."""
    from scipy.optimize import minimize

    bounds = _make_bounds(min_freq, max_freq, n_bands)

    def objective(x):
        return _error(x, z, target_db, weights, fs, n_bands)

    result = minimize(
        objective, params, method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': 80, 'ftol': 1e-10, 'gtol': 1e-8},
    )

    if result.fun < objective(params):
        return result.x
    return params


# ── Public API ────────────────────────────────────────────────────────

def fit_bands(freqs: np.ndarray, target_db: np.ndarray, fs: float,
              max_bands: int = 24,
              ) -> Tuple[List[BiquadCoeffs], dict]:
    """Fit target curve to cascaded biquads using differential evolution.

    Starts from greedy seed, uses DE/best/1/bin to escape local optima,
    with periodic L-BFGS-B polish.

    Returns (list_of_BiquadCoeffs, metadata_dict).
    """
    n_bands = max_bands
    D = n_bands * 3
    NP = max(30, min(POP_MULTIPLIER * D, 300))
    sigma = SIGMA_INIT

    # Frequency limits from data
    min_freq = freqs[0]
    max_freq = freqs[-1]

    # Precompute
    z = np.exp(-1j * 2.0 * np.pi * freqs / fs)
    weights = np.array([ear_sensitivity(f) for f in freqs])

    # ── Phase 0: Greedy seed ──────────────────────────────────────────
    seed_params = _greedy_to_params(freqs, target_db, fs, n_bands)
    seed_fitness = _error(seed_params, z, target_db, weights, fs, n_bands)

    # ── Phase 1: Initialize population around greedy seed ─────────────
    pop = _filled_population(NP, seed_params, min_freq, max_freq, n_bands, sigma)
    F_vals = np.full(NP, 0.5)
    CR_vals = np.full(NP, 0.9)

    # Evaluate
    fitness = np.array([
        _error(pop[j], z, target_db, weights, fs, n_bands)
        for j in range(NP)
    ])

    # Seed might be best, but perturbation might have found better
    best_idx = int(np.argmin(fitness))
    best_fitness = fitness[best_idx]
    best_params = pop[best_idx].copy()

    gen_since_improvement = 0
    total_generations = 0

    # Scale max gens with problem size
    max_gens = min(MAX_GENERATIONS, max(30, 300 // max(1, n_bands // 2)))

    for gen in range(max_gens):
        total_generations = gen + 1

        # Generate trial population
        trial, F_vals, CR_vals = _best1_bin(
            pop, fitness, F_vals, CR_vals, best_idx,
            min_freq, max_freq, n_bands,
        )

        # Evaluate trials
        trial_fitness = np.array([
            _error(trial[j], z, target_db, weights, fs, n_bands)
            for j in range(NP)
        ])

        # Selection
        improved = trial_fitness < fitness
        pop[improved] = trial[improved]
        fitness[improved] = trial_fitness[improved]

        # Track best
        current_best_idx = int(np.argmin(fitness))
        current_best_fit = fitness[current_best_idx]

        if current_best_fit < best_fitness - EARLY_STOP_TOL:
            best_fitness = current_best_fit
            best_params = pop[current_best_idx].copy()
            best_idx = current_best_idx
            gen_since_improvement = 0
        else:
            gen_since_improvement += 1

        # Periodic polish
        if gen % POLISH_GENS == (POLISH_GENS - 1) and gen > 0:
            best_params = _polish_best(
                best_params, z, target_db, weights, fs, n_bands,
                min_freq, max_freq,
            )
            new_fit = _error(best_params, z, target_db, weights, fs, n_bands)
            if new_fit < best_fitness:
                best_fitness = new_fit
                pop[best_idx] = best_params
                fitness[best_idx] = best_fitness
                gen_since_improvement = 0

        # Early stop
        if gen_since_improvement >= EARLY_STOP_GENS:
            # Final polish attempt
            best_params = _polish_best(
                best_params, z, target_db, weights, fs, n_bands,
                min_freq, max_freq,
            )
            new_fit = _error(best_params, z, target_db, weights, fs, n_bands)
            if new_fit < best_fitness - EARLY_STOP_TOL:
                best_fitness = new_fit
                gen_since_improvement = 0
                pop[best_idx] = best_params
                fitness[best_idx] = best_fitness
                continue
            break

    # ── Phase 3: Final polish ─────────────────────────────────────────
    best_params = _polish_best(
        best_params, z, target_db, weights, fs, n_bands,
        min_freq, max_freq,
    )
    final_fitness = _error(best_params, z, target_db, weights, fs, n_bands)

    # ── Convert to biquads ────────────────────────────────────────────
    biquads = _params_to_biquads(best_params, fs, n_bands)

    # Build metadata
    bands = []
    band_data = [(best_params[i*3], best_params[i*3+1], best_params[i*3+2])
                 for i in range(n_bands)]
    band_data.sort(key=lambda x: x[0])
    for lf0, g, q in band_data:
        bands.append({
            "type": "peaking",
            "f0": float(10.0 ** lf0),
            "gain_db": float(g),
            "Q": float(q),
        })

    return biquads, {
        "name": "diffevo",
        "bands": bands,
        "fitness": final_fitness,
        "greedy_fitness": seed_fitness,
        "improvement": seed_fitness - final_fitness,
        "generations": total_generations,
        "np": NP,
    }
