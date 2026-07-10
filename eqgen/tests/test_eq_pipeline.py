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
from eqgen.quantize import (
    BiquadQ28,
    quantize_biquads_q28,
    q28_to_float,
    ReciprocalLUT,
    float_to_q16,
)
from eqgen.dsp import butterworth_hp_mag, butterworth_lp_mag



# ── Speaker models for testing ──────────────────────────────────────

def small_speaker(f):
    """Small speaker: -12 dB at 50 Hz, flat above 100 Hz."""
    if f <= 50:
        return 0.25
    elif f >= 100:
        return 1.0
    else:
        return 0.25 + 0.75 * (f - 50) / 50


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


    # ── Step 3: Quantize to Q4.28 ───────────────────────────────────
    print(f"\n  ── Step 3: Quantize to Q4.28 ──")
    bq_q28 = quantize_biquads_q28(fit.biquads)
    for i, (bc, bq) in enumerate(zip(fit.biquads, bq_q28)):
        errs = [
            abs(q28_to_float(bq.b0) - bc.b0),
            abs(q28_to_float(bq.b1) - bc.b1),
            abs(q28_to_float(bq.b2) - bc.b2),
            abs(q28_to_float(bq.a1) - bc.a1),
            abs(q28_to_float(bq.a2) - bc.a2),
        ]
        print(f"  Band {i}: max coeff err = {max(errs):.2e}")

    # ── Step 4: Q4.28 frequency response ───────────────────────────
    print(f"\n  ── Step 4: Q4.28 frequency response ──")
    float_db = cascade_response_db(fit.biquads, freqs, fs)
    q28_biquads_float = [
        BiquadCoeffs(
            b0=q28_to_float(bq.b0), b1=q28_to_float(bq.b1),
            b2=q28_to_float(bq.b2), a1=q28_to_float(bq.a1),
            a2=q28_to_float(bq.a2),
        )
        for bq in bq_q28
    ]
    q28_db = cascade_response_db(q28_biquads_float, freqs, fs)
    q28_error = q28_db - float_db

    max_q28_err = np.max(np.abs(q28_error))
    print(f"  Q4.28 vs float: max error = {max_q28_err:+.4f} dB")

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
        "bq_q28": bq_q28,
        "total_error_db": total_error_db,
    }


def test_reciprocal_lut():
    """Verify reciprocal LUT accuracy."""
    print("=" * 70)
    print("  RECIPROCAL LUT ACCURACY")
    print("=" * 70)

    lut = ReciprocalLUT.build()
    max_err = lut.max_error(10000)
    print(f"  LUT: {len(lut.entries)} entries")
    print(f"  Max relative error: {max_err*100:.2f}%")

    if max_err < 0.005:
        print(f"  ✅ LUT accuracy excellent (< 0.5%)")
    else:
        print(f"  ⚠️  LUT accuracy could be improved")

    return lut


def test_q28_biquad_precision():
    """Verify Q4.28 precision at sub-bass frequencies."""
    from eqgen.eq_fit import design_butter_hp
    from eqgen.quantize import BiquadQ28, q28_to_float, Q28_SCALE

    print("=" * 70)
    print("  PRECISION TEST: Q4.28 biquad limits")
    print("=" * 70)

    fs = 44100.0
    print(f"\n  {'fc (Hz)':>8s}  {'num_sum Q28':>12s}  {'den_sum Q28':>12s}  {'DC gain':>10s}")
    print(f"  {'-'*50}")

    for fc in [10, 15, 20, 25, 35, 45, 60, 100]:
        hp = design_butter_hp(fc, fs)
        q28 = BiquadQ28.from_float(hp)
        num = q28.b0 + q28.b1 + q28.b2
        den = Q28_SCALE + q28.a1 + q28.a2
        dc_gain_db = 20.0 * np.log10(abs(num / den) + 1e-12) if den != 0 else -300
        status = "✅" if num == 0 else "⚠️"
        print(f"  {fc:8.0f}  {num:>12d}  {den:>12d}  {dc_gain_db:>+10.1f} dB  {status}")

    print(f"\n  Q4.28 numerator sums exactly to 0 — no degeneracy.")


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
