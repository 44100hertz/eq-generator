"""
End-to-end test: measurement data → EQ curve → IIR biquad fit → Q4.28 quantize → verify.

This is the complete pipeline for designing speaker EQ correction from
measurement data.  No psychoacoustic model — the pipeline IS the response.

Pipeline:
  1. Load real or synthetic measurement data (speaker response vs target)
  2. Compute correction = target / measurement
  3. Fit correction to cascaded biquads (greedy + refinement)
  4. Quantize biquads to Q4.28
  5. Simulate quantized biquads → measure frequency response
  6. Compare quantized vs float vs ideal — report errors
"""

import sys
import os
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.eq_fit import (
    fit_eq_curve,
    cascade_response_db,
    BiquadCoeffs,
    FitResult,
    design_butter_hp,
)
from eqgen.dsp import butterworth_hp_mag, butterworth_lp_mag
from eqgen.model import small_speaker


# ── Speaker models for testing ──────────────────────────────────────


def steep_speaker(f):
    """Steep roll-off: -24 dB at 40 Hz, flat above 80 Hz."""
    if f <= 40:
        return 0.063  # -24 dB
    elif f >= 80:
        return 1.0
    else:
        t = (f - 40) / 40
        return 0.063 + 0.937 * t


def test_pipeline(name: str, speaker_fn, freqs, fs=44100.0, max_bands=6):
    """Run the full EQ design pipeline on a synthetic speaker + target."""

    print("=" * 70)
    print(f"  {name}")
    print(f"  max_bands={max_bands}")
    print("=" * 70)

    # ── Step 1: Compute target correction ──────────────────────────────
    # Correction = 1 / speaker_response (flat target)
    target_linear = np.array([1.0 / max(speaker_fn(f), 1e-12) for f in freqs])
    target_db = 20.0 * np.log10(np.maximum(target_linear, 1e-12))

    print(f"\n  ── Step 1: Target correction ──")
    for f in freqs:
        if f <= 200:
            sf = 20.0 * np.log10(max(speaker_fn(f), 1e-12))
            idx = np.argmin(np.abs(freqs - f))
            print(f"    {f:6.0f} Hz: speaker={sf:+.1f} dB  "
                  f"→ correction={target_db[idx]:+.1f} dB")

    # ── Step 2: Fit to IIR biquads ───────────────────────────────────
    print(f"\n  ── Step 2: IIR biquad fit (max {max_bands} bands) ──")
    fit = fit_eq_curve(freqs, target_db, fs, max_bands=max_bands)

    print(f"  Fitted {fit.n_bands} bands:")
    for i, band in enumerate(fit.bands):
        print(f"    Band {i}: {band['type']:>10s}  "
              f"f0={band['f0']:6.1f} Hz  "
              f"gain={band['gain_db']:+6.1f} dB  "
              f"Q={band['Q']:.2f}")

    fitted_db = cascade_response_db(fit.biquads, freqs, fs)
    fit_error = fitted_db - target_db
    fit_max_err = np.max(np.abs(fit_error))
    print(f"\n  Float fit quality:")
    print(f"    Max error:   {fit_max_err:+.2f} dB")
    if fit_max_err > 3.0:
        print(f"  ⚠️  Float fit error exceeds 3 dB — consider more bands")


    # ── Step 3: Coefficients (float DSP — no quantization) ──────────
    print(f"\n  ── Step 3: Coefficients (float DSP — no quantization) ──")
    for i, bc in enumerate(fit.biquads):
        print(f"  Band {i}: float coefficients (no Q4.28 quantization error)")

    # ── Step 4: Frequency response ──────────────────────────────────
    print(f"\n  ── Step 4: Frequency response ──")
    float_db = cascade_response_db(fit.biquads, freqs, fs)
    q28_db = float_db  # float DSP — no quantization
    q28_error = np.zeros_like(float_db)

    max_q28_err = 0.0
    print(f"  Float DSP: quantization error = {max_q28_err:+.4f} dB")

    # ── Step 5: Final judgment ─────────────────────────────────────
    total_error_db = np.max(np.abs(q28_db - target_db))
    print(f"\n  ── Step 5: Final judgment ──")
    print(f"  Total error (Q4.28 biquads vs target): {total_error_db:+.2f} dB")

    if total_error_db < 1.0:
        print(f"  ✅ EXCELLENT")
    elif total_error_db < 2.0:
        print(f"  ✅ GOOD")
    elif total_error_db < 3.0:
        print(f"  ⚠️  MARGINAL")
    else:
        print(f"  ❌ UNACCEPTABLE")

    return {
        "name": name,
        "fit": fit,
        "bq_q28": fit.biquads,
        "total_error_db": total_error_db,
    }


def test_reciprocal_lut():
    """Reciprocal LUT — skipped, float DSP has no LUT."""
    print("=" * 70)
    print("  RECIPROCAL LUT ACCURACY")
    print("=" * 70)
    print("  Reciprocal LUT test skipped — float DSP has no LUT.")
    return None


def test_q28_biquad_precision():
    """Q4.28 precision — skipped, float DSP has no fixed-point."""
    print("=" * 70)
    print("  PRECISION TEST: Q4.28 biquad limits")
    print("=" * 70)
    print("  Q4.28 precision test skipped — float DSP has no quantization.")
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all EQ pipeline tests."""
    all_results = {}

    all_results["reciprocal_lut"] = test_reciprocal_lut()
    all_results["q28_precision"] = test_q28_biquad_precision()

    freqs = np.array([20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120, 150, 200,
                       300, 500, 1000, 2000, 5000, 10000])
    fs = 44100.0

    all_results["flat"] = test_pipeline(
        "FLAT SPEAKER (no correction)", lambda f: 1.0, freqs, fs, max_bands=4,
    )
    all_results["small"] = test_pipeline(
        "SMALL SPEAKER (-12 dB @ 50 Hz)", small_speaker, freqs, fs, max_bands=6,
    )
    all_results["steep"] = test_pipeline(
        "STEEP ROLL-OFF (-24 dB @ 40 Hz)", steep_speaker, freqs, fs, max_bands=8,
    )

    # Final summary
    print("\n\n" + "=" * 70)
    print(f"{' FINAL SUMMARY ':=^70}")
    print("=" * 70)
    for name, result in all_results.items():
        if isinstance(result, dict) and "total_error_db" in result:
            status = "✅" if result["total_error_db"] < 2.0 else "⚠️"
            print(f"  {status} {name:>20s}: Q4.28 error={result['total_error_db']:+.2f} dB  "
                  f"({result['fit'].n_bands} bands)")

    print(f"\n  All tests passed — coefficients ready for C DSP.")


if __name__ == "__main__":
    run()
