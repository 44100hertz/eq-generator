"""
Shared pipeline: Welch FFT \u2192 CV-weighted smoothing \u2192 correction \u2192 IIR fit \u2192 C DSP.

This module is the single source of truth for the EQ design and audio processing
pipeline.  All CLI entry points (eqgen, audition, wire, export) delegate here.

Exports:
  run_pipeline()   \u2014 WAV inputs \u2192 (freqs, gains_db, sample_rate, max_gain_db, efficacy)
  design_eq()      \u2014 EQ curve \u2192 biquad coefficients
  curve_to_json()  \u2014 serialise (freqs, gains_db) for JSON output
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from eqgen.eq_fit import BiquadCoeffs, cascade_response_db, fit_eq_curve
from eqgen.presets import MAX_IIR_BANDS
from eqgen.io import read_wav
from eqgen.dsp import ratio_to_db, pre_gain_from_max_gain
from eqgen.analysis import (
    welch_stats,
    _find_viable_range,
    _smooth_kernel,
    BOTTOM_F,
    TOP_F,
)

ROOT = Path(__file__).resolve().parent.parent

# Evaluation grid for the output curve \u2014 matches the IIR fitter's internal
# log-spaced grid so no re-interpolation is needed downstream.
N_EVAL = 512


# ═══════════════════════════════════════════════════════════════════════════════
# Core pipeline: WAV inputs \u2192 EQ curve
# ═══════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    measurement_paths: List[str],
    target_path: str,
    noise_path: Optional[str] = None,
    bass_enhancer_cutoff: Optional[float] = None,
    sample_rate_override: Optional[float] = None,
    smooth_exponent: float = 1.0,
    detailed: bool = False,
    n_eval: int = 512,
    house_curve: Optional[list] = None,
):
    """Run the full EQ correction pipeline.

    Without detailed: returns (freqs_hz, gains_db, sample_rate, max_gain_db, efficacy).
    With detailed=True: returns a dict with all intermediate data for
    visualization.

    Steps:
      1. Read & validate WAV files
      2. Welch FFT on each input
      3. Pool per-bin stats across measurement files
      4. Merge noise-floor CV + inflate CV in noise-dominated bins
      5. CV-weighted kernel smoothing on measurement + spectral subtraction
      6. Target = linear fit of raw target (dB vs log-freq),
         shifted to measurement midrange level
      7. Correction = target / measurement
      8. Compute harmonic efficacy from correction curve (if cutoff provided)
    """

    # 1. Read & validate
    wavs = []
    for path in measurement_paths:
        samples, rate = read_wav(path)
        wavs.append((samples, rate, path))

    sample_rate = sample_rate_override or wavs[0][1]
    for s, r, p in wavs[1:]:
        assert r == sample_rate, f"Sample rate mismatch in {p}: {r} vs {sample_rate}"

    if target_path:
        target_samples, target_rate = read_wav(target_path)
        assert target_rate == sample_rate, \
            f"Target sample rate mismatch: {target_rate} vs {sample_rate}"
    else:
        target_samples = None

    # 2. Welch FFT
    all_stats = [welch_stats(s, sample_rate) for s, _, _ in wavs]
    bin_count = len(all_stats[0])
    for s in all_stats[1:]:
        assert len(s) == bin_count, "Frequency ranges differ"

    # 3. Pool stats across measurement files
    pooled = []
    for i in range(bin_count):
        total_n = sum(s[i]["count"] for s in all_stats)
        total_sum = sum(s[i]["mean"] * s[i]["count"] for s in all_stats)
        total_sum_sq = sum(s[i]["mean_sq"] * s[i]["count"] for s in all_stats)
        pooled_cv = sum(s[i]["cv"] * s[i]["count"] for s in all_stats) / total_n
        pooled.append({
            "freq": all_stats[0][i]["freq"],
            "count": total_n,
            "mean": total_sum / total_n,
            "mean_sq": total_sum_sq / total_n,
            "cv": pooled_cv,
        })

    # 4. Noise merge: max-CV + noise-floor inflation
    noise_stats = None
    if noise_path:
        noise_samples, noise_rate = read_wav(noise_path)
        assert noise_rate == sample_rate
        noise_stats = welch_stats(noise_samples, sample_rate)

        # 4a. max() merge \u2014 high-CV intermittent noise dominates
        for i in range(min(len(pooled), len(noise_stats))):
            pooled[i]["cv"] = max(pooled[i]["cv"], noise_stats[i]["cv"])

        # 4b. Noise-floor inflation \u2014 stationary noise with poor SNR
        for i in range(len(pooled)):
            P_meas = max(pooled[i]["mean_sq"], 1e-20)
            P_noise = max(noise_stats[i]["mean_sq"], 1e-20)
            cv_n = noise_stats[i]["cv"]
            # Inflate noise power by 1\u03c3 \u2014 accounts for peak, not mean, noise
            effective_noise = P_noise * (1.0 + cv_n)
            noise_ratio = min(0.99, effective_noise / P_meas)
            pooled[i]["cv"] = pooled[i]["cv"] / np.sqrt(max(0.01, 1.0 - noise_ratio))

    # 5. CV-weighted kernel smoothing on measurement only.
    meas_freqs = np.array([s["freq"] for s in pooled])
    meas_raw = np.array([s["mean"] for s in pooled])
    meas_cv = np.array([s["cv"] for s in pooled])

    if target_path:
        target_stats = welch_stats(target_samples, sample_rate)
        targ_freqs = np.array([s["freq"] for s in target_stats])
        targ_raw_db = ratio_to_db(np.array([s["mean"] for s in target_stats]))
    else:
        # No target WAV \u2014 use a flat 0 dB line spanning the measurement range
        target_stats = None
        targ_freqs = np.array([BOTTOM_F, TOP_F])
        targ_raw_db = np.array([0.0, 0.0])

    noise_freqs = None
    noise_raw = None
    noise_cv = None
    if noise_stats is not None:
        noise_freqs = np.array([s["freq"] for s in noise_stats])
        noise_raw = np.array([s["mean"] for s in noise_stats])
        noise_cv = np.array([s["cv"] for s in noise_stats])

    # 6. CV-weighted kernel smoothing on uniform log-spaced grid.
    MIN_CV = 0.52   # Rayleigh CV \u2014 the noise floor for stationary signals
    BASE_BW = 0.08  # octaves at MIN_CV (~3 FFT bins at 400 Hz)

    eval_freqs = np.logspace(np.log10(BOTTOM_F), np.log10(TOP_F), n_eval)

    # Estimate local CV at each eval point to scale bandwidth
    cv_at_eval = _smooth_kernel(meas_freqs, meas_cv, meas_cv, eval_freqs,
                                bandwidth_oct=0.3)
    bw_meas = np.maximum(BASE_BW * (cv_at_eval / MIN_CV) ** smooth_exponent, 0.01)
    meas_vals = _smooth_kernel(meas_freqs, meas_raw, meas_cv, eval_freqs,
                               bandwidth_oct=bw_meas)

    # Spectral subtraction
    if noise_freqs is not None:
        noise_vals = _smooth_kernel(noise_freqs, noise_raw, noise_cv, eval_freqs)
        meas_vals = np.maximum(meas_vals - noise_vals, meas_vals * 0.01)

    # 6. Target from WAV (linear fit in dB vs log-freq) or flat 0 dB line.
    #    House curve is an additive adjustment applied on top.
    eval_log10 = np.log10(eval_freqs)
    mid = (eval_freqs >= 500) & (eval_freqs <= 2000)

    log10_f_targ = np.log10(targ_freqs)
    slope, intercept = np.polyfit(log10_f_targ, targ_raw_db, 1)
    targ_db_eval = slope * eval_log10 + intercept

    if house_curve is not None:
        curve_freqs = np.array([p[0] for p in house_curve])
        curve_db = np.array([p[1] for p in house_curve])
        targ_db_eval = targ_db_eval + np.interp(eval_freqs, curve_freqs, curve_db)

    meas_mid_db = float(np.mean(ratio_to_db(meas_vals[mid]))) if mid.any() else 0.0
    targ_mid_db = float(np.mean(targ_db_eval[mid])) if mid.any() else 0.0
    targ_db_eval += meas_mid_db - targ_mid_db
    target_vals = 10.0 ** (targ_db_eval / 20.0)

    # 7. Correction
    with np.errstate(divide='ignore', invalid='ignore'):
        corr = np.where(meas_vals > 1e-20, target_vals / meas_vals, 1e6)

    # 8. Bass enhancer efficacy: compute h2/h3 from measurement data.
    efficacy = {"h2_amp": 0.0, "h3_amp": 0.0}
    if bass_enhancer_cutoff is not None and bass_enhancer_cutoff > 0:
        from eqgen.model import compute_harmonic_efficacy
        efficacy = compute_harmonic_efficacy(eval_freqs, ratio_to_db(corr),
                                             bass_enhancer_cutoff)

    # 9. Clamp correction flat outside the CV-defined viable range.
    lo_f, hi_f = _find_viable_range(meas_freqs, meas_cv)
    lo_idx = int(np.searchsorted(eval_freqs, lo_f))
    hi_idx = int(np.searchsorted(eval_freqs, hi_f))
    if lo_idx > 0:
        corr[:lo_idx] = corr[lo_idx]
    if hi_idx < len(corr) - 1:
        corr[hi_idx + 1:] = corr[hi_idx]

    freqs = eval_freqs
    gains_db = ratio_to_db(corr)

    # Normalize: re-center correction so midrange mean = 0 dB.
    mid_mask = (freqs >= 500.0) & (freqs <= 2000.0)
    if np.any(mid_mask):
        mid_mean = float(np.mean(gains_db[mid_mask]))
        gains_db = gains_db - mid_mean

    # Pre-gain: the max positive gain of the *normalized* correction curve.
    max_gain_db = max(0.0, float(np.max(gains_db)))

    if detailed:
        raw_resp = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                    for s in pooled]
        if target_stats is not None:
            raw_target = [{"freq": float(s["freq"]), "db": ratio_to_db(s["mean"])}
                          for s in target_stats]
        else:
            raw_target = []
        cv_data = [{"freq": p["freq"], "cv": p["cv"]} for p in pooled]
        noise_data = []
        if noise_stats is not None:
            noise_data = [{"freq": s["freq"], "cv": s["cv"]} for s in noise_stats]
        meas_sub = [{"freq": float(f), "db": ratio_to_db(v)}
                    for f, v in zip(eval_freqs, meas_vals)]
        targ_sub = [{"freq": float(f), "db": float(d)}
                    for f, d in zip(eval_freqs, targ_db_eval)]

        return {
            "freqs": freqs.tolist(),
            "gains_db": gains_db.tolist(),
            "sample_rate": sample_rate,
            "max_gain_db": max_gain_db,
            "efficacy": efficacy,
            "raw_measurement": raw_resp,
            "raw_target": raw_target,
            "noise_cv": cv_data,
            "noise_floor": noise_data,
            "meas_subtracted": meas_sub,
            "target_resampled": targ_sub,
        }

    return freqs, gains_db, sample_rate, max_gain_db, efficacy


# ═══════════════════════════════════════════════════════════════════════════════
# Serialisation
# ═══════════════════════════════════════════════════════════════════════════════

def curve_to_json(freqs: np.ndarray, gains_db: np.ndarray) -> list:
    """Convert (freqs, gains_db) to JSON-serialisable list of {freq, gain_db}."""
    return [
        {"freq": round(float(f), 2), "gain_db": round(float(g), 3)}
        for f, g in zip(freqs, gains_db)
    ]


# ═══════════════════════════════════════════════════════════════════════════════
# IIR design: EQ curve \u2192 biquad coefficients
# ═══════════════════════════════════════════════════════════════════════════════

def design_eq(
    freqs: np.ndarray,
    target_db: np.ndarray,
    fs: float,
    max_bands: int = MAX_IIR_BANDS,
    min_peaking_freq: float = 0.0,
) -> Tuple[List[float], List[dict], np.ndarray, np.ndarray, np.ndarray]:
    """Design an IIR EQ to match a target curve using cascaded biquads.

    Returns (coeffs_flat, bands, freqs, fit_target_db, fitted_db).
    """
    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=min_peaking_freq if min_peaking_freq > 0 else freqs[0])

    coeffs = []
    for bc in fit.biquads:
        coeffs.extend([bc.b0, bc.b1, bc.b2, bc.a1, bc.a2])

    fitted_db = cascade_response_db(fit.biquads, freqs, fs)

    return coeffs, fit.bands, freqs, target_db, fitted_db
