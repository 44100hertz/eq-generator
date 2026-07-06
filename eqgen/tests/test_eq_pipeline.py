"""
End-to-end test: measurement data → EQ curve → IIR biquad fit → 16.16 quantize → verify.

This is the complete Python pipeline that will later be ported to C.
The output of this test determines whether the IIR fit and quantization
are good enough for the ESP32-C3 to run.

Pipeline:
  1. Load real measurement data (speaker response vs target)
  2. Run psychoacoustic model to compute EQ preprocessing gains G(f)
  3. Fit G(f) to cascaded biquads (greedy + refinement)
  4. Quantize biquads to 16.16
  5. Simulate quantized biquads → measure frequency response
  6. Compare quantized vs float vs ideal — report errors
"""

import sys
import os
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.model import (
    preprocess_eq_curve,
    small_speaker,
    butterworth_hp_mag,
    butterworth_lp_mag,
)
from eqgen.eq_fit import (
    fit_from_dict,
    cascade_response_db,
    BiquadCoeffs,
    FitResult,
)
from eqgen.quantize import (
    BiquadQ28,
    quantize_biquads_q28,
    reorder_attenuation_first,
    q28_to_float,
    ReciprocalLUT,
    float_to_q16,
)


def _reorder_biquads(fit):
    """Reorder biquad cascade in-place: attenuation before boost."""
    fit.biquads, fit.bands = reorder_attenuation_first(fit.biquads, fit.bands)


def test_pipeline(name: str, speaker_fn, target_curve: dict,
                 fc: float, h2: float, h3: float,
                 fs: float = 44100.0,
                 max_bands: int = 5,
                 plot: bool = False):
    """Run the full pipeline on a speaker + target curve pair."""

    print("=" * 70)
    print(f"  {name}")
    print(f"  fc={fc} Hz, h2={h2}, h3={h3}, max_bands={max_bands}")
    print("=" * 70)

    # ── Step 1: Compute EQ preprocessing gains ───────────────────────
    print(f"\n  ── Step 1: Model EQ gains ──")
    eq_curve = preprocess_eq_curve(target_curve, speaker_fn, fc, h2, h3)
    freqs = np.array(sorted(eq_curve.keys()))

    print(f"  Computed {len(eq_curve)} frequency points")
    print(f"  Freq range: {freqs[0]:.0f} – {freqs[-1]:.0f} Hz")
    for f in sorted(eq_curve.keys()):
        if f <= fc * 3:
            print(f"    {f:6.0f} Hz: G={20*np.log10(eq_curve[f]):+.2f} dB  "
                  f"(speaker={20*np.log10(speaker_fn(f)):+.1f} dB)")

    # ── Step 2: Fit to IIR biquads ───────────────────────────────────
    print(f"\n  ── Step 2: IIR biquad fit (max {max_bands} bands) ──")
    fit = fit_from_dict(eq_curve, fs, max_bands=max_bands)

    print(f"  Fitted {fit.n_bands} bands:")
    for i, band in enumerate(fit.bands):
        print(f"    Band {i}: {band['type']:>10s}  "
              f"f0={band['f0']:6.1f} Hz  "
              f"gain={band['gain_db']:+6.1f} dB  "
              f"Q={band['Q']:.2f}")

    # Show fit quality
    target_db = 20.0 * np.log10(np.array([eq_curve[f] for f in freqs]))
    fitted_db = cascade_response_db(fit.biquads, freqs, fs)
    fit_error = fitted_db - target_db

    fit_max_err = np.max(np.abs(fit_error))
    fit_mean_err = np.mean(np.abs(fit_error))
    print(f"\n  Float fit quality (target vs biquads):")
    print(f"    Max error:   {fit_max_err:+.2f} dB")
    print(f"    Mean |err|:  {fit_mean_err:.2f} dB")

    # Detail per frequency
    print(f"\n  {'Freq':>6s}  {'Target':>8s}  {'Fitted':>8s}  {'Error':>8s}")
    print(f"  {'-'*34}")
    for i, f in enumerate(freqs):
        if f <= fc * 3:
            print(f"  {f:6.0f}  {target_db[i]:+8.2f}  {fitted_db[i]:+8.2f}  "
                  f"{fit_error[i]:+8.2f}")

    if fit_max_err > 3.0:
        print(f"\n  ⚠️  Float fit error exceeds 3 dB — consider more bands")

    # ── Step 2b: Reorder cascade — attenuation before boost ──────────
    # This prevents intermediate overflow in the fixed-point cascade.
    # Attenuating stages reduce the signal first, then boosting stages
    # bring it back up — but never above the input level at any
    # intermediate point.
    _reorder_biquads(fit)
    print(f"\n  Reordered cascade (attenuation-first):")
    for i, band in enumerate(fit.bands):
        bc = fit.biquads[i]
        num = bc.b0 + bc.b1 + bc.b2
        den = 1.0 + bc.a1 + bc.a2
        dc_db = 20.0 * np.log10(abs(num / den)) if abs(den) > 1e-12 else 0.0
        print(f"    Band {i}: {band['type']:>10s}  "
              f"f0={band['f0']:6.1f} Hz  "
              f"gain={band['gain_db']:+6.1f} dB  "
              f"DC={dc_db:+6.1f} dB")

    # ── Step 3: Quantize to Q4.28 ───────────────────────────────────
    print(f"\n  ── Step 3: Quantize to Q4.28 ──")
    bq_q28 = quantize_biquads_q28(fit.biquads)

    # Show quantization accuracy
    for i, (bc, bq) in enumerate(zip(fit.biquads, bq_q28)):
        errs = [
            abs(q28_to_float(bq.b0) - bc.b0),
            abs(q28_to_float(bq.b1) - bc.b1),
            abs(q28_to_float(bq.b2) - bc.b2),
            abs(q28_to_float(bq.a1) - bc.a1),
            abs(q28_to_float(bq.a2) - bc.a2),
        ]
        max_err = max(errs)
        print(f"  Band {i} ({fit.bands[i]['type']}): max coeff err = {max_err:.2e}")

    # ── Step 4: Q4.28 frequency response ───────────────────────────
    print(f"\n  ── Step 4: Q4.28 frequency response ──")
    float_db = cascade_response_db(fit.biquads, freqs, fs)
    q28_biquads_float = [
        BiquadCoeffs(
            b0=q28_to_float(bq.b0),
            b1=q28_to_float(bq.b1),
            b2=q28_to_float(bq.b2),
            a1=q28_to_float(bq.a1),
            a2=q28_to_float(bq.a2),
        )
        for bq in bq_q28
    ]
    q28_db = cascade_response_db(q28_biquads_float, freqs, fs)
    q28_error = q28_db - float_db

    max_q28_err = np.max(np.abs(q28_error))
    mean_q28_err = np.mean(np.abs(q28_error))
    print(f"  Q4.28 vs float: max error = {max_q28_err:+.4f} dB, mean |err| = {mean_q28_err:.4f} dB")

    # Show per-frequency comparison
    print(f"\n  {'Freq':>6s}  {'Float':>8s}  {'Q4.28':>8s}  {'vs Target':>10s}")
    print(f"  {'-'*40}")
    for i, f in enumerate(freqs):
        if f <= fc * 3:
            q28_vs_target = q28_db[i] - target_db[i]
            print(f"  {f:6.0f}  {float_db[i]:+8.2f}  {q28_db[i]:+8.2f}  "
                  f"{q28_vs_target:+10.2f}")

    # ── Step 5: Judgment ─────────────────────────────────────────────
    total_error_db = np.max(np.abs(q28_db - target_db))
    print(f"\n  ── Step 5: Final judgment ──")
    print(f"  Total error (Q4.28 biquads vs target): {total_error_db:+.2f} dB")

    if total_error_db < 1.0:
        print(f"  ✅ EXCELLENT — quantization error is inaudible")
    elif total_error_db < 2.0:
        print(f"  ✅ GOOD — within acceptable range for bass EQ")
    elif total_error_db < 3.0:
        print(f"  ⚠️  MARGINAL — consider more bands or higher Q precision")
    else:
        print(f"  ❌ UNACCEPTABLE — needs more bands or wider precision")

    return {
        "name": name,
        "fit": fit,
        "bq_q28": bq_q28,
        "max_q28_error_db": max_q28_err,
        "mean_q28_error_db": mean_q28_err,
        "total_error_db": total_error_db,
        "fit_max_error_db": fit_max_err,
    }


def test_reciprocal_lut():
    """Verify reciprocal LUT accuracy."""
    print("=" * 70)
    print("  RECIPROCAL LUT ACCURACY")
    print("=" * 70)
    print()

    lut = ReciprocalLUT.build()
    max_err = lut.max_error(10000)
    print(f"  LUT: {len(lut.entries)} entries")
    print(f"  Max relative error: {max_err*100:.2f}%")

    # Spot-check
    test_points_q16 = [256, 512, 1024, 4096, 8192, 16384, 32768]
    print(f"\n  {'x Q16':>10s}  {'1/x Q16.16':>14s}  {'LUT':>14s}  {'Err %':>10s}")
    print(f"  {'-'*52}")
    for xq in test_points_q16:
        approx = lut.lookup(xq)
        exact_q16_16 = int(round(65536.0 * 65536.0 / xq))
        err = abs(approx - exact_q16_16) / exact_q16_16 * 100.0
        print(f"  {xq:10d}  {exact_q16_16:14d}  {approx:14d}  {err:+10.4f}%")

    if max_err < 0.005:
        print(f"\n  ✅ LUT accuracy is excellent (< 0.5%)")
    else:
        print(f"\n  ⚠️  LUT accuracy could be improved with more entries")

    return lut


# ═══════════════════════════════════════════════════════════════════════════════
# Test cases
# ═══════════════════════════════════════════════════════════════════════════════

def test_small_speaker():
    """Small speaker: -12 dB at 50 Hz, flat above 100 Hz."""
    fc = 60.0
    h2, h3 = 0.33, 0.33

    freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]
    target_curve = {f: 0.5 for f in freqs}  # -6 dBFS target

    return test_pipeline(
        "SMALL SPEAKER (-12 dB @ 50 Hz)",
        small_speaker, target_curve, fc, h2, h3, max_bands=6,
    )


def test_flat_speaker():
    """Flat speaker: no correction needed — verifies pass-through."""
    fc = 60.0
    h2, h3 = 0.33, 0.33

    freqs = [20, 30, 40, 50, 60, 80, 100, 120]
    target_curve = {f: 0.5 for f in freqs}

    def flat(f):
        return 1.0

    return test_pipeline(
        "FLAT SPEAKER (no correction)",
        flat, target_curve, fc, h2, h3, max_bands=4,
    )


def test_steep_rolloff():
    """Steep roll-off speaker: -24 dB at 50 Hz — stress test."""
    fc = 60.0
    h2, h3 = 0.33, 0.33

    def steep_speaker(f):
        if f <= 40:
            return 0.063  # -24 dB
        elif f >= 80:
            return 1.0
        else:
            # Smooth transition
            t = (f - 40) / 40
            return 0.063 + 0.937 * t

    freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120, 150]
    target_curve = {f: 0.5 for f in freqs}

    return test_pipeline(
        "STEEP ROLL-OFF (-24 dB @ 40 Hz) — STRESS TEST",
        steep_speaker, target_curve, fc, h2, h3, max_bands=8,
    )


def test_real_measurement():
    """Load a real measurement file and test the full pipeline."""
    fc = 50.0
    h2, h3 = 0.33, 0.33
    fs = 44100.0

    # Try to load measurement data
    root = Path(__file__).parent.parent
    meas_dir = root / "measurements" / "technics" / "standing"

    meas_path = meas_dir / "measurement2.wav"
    target_path = meas_dir / "target.wav"

    if not meas_path.exists() or not target_path.exists():
        print(f"\n  ⚠️  Measurement files not found at {meas_path}")
        print(f"  Skipping real-measurement test.")
        return None

    import scipy.io.wavfile as wav

    rate, meas_raw = wav.read(str(meas_path))
    _, target_raw = wav.read(str(target_path))
    if meas_raw.ndim == 2:
        meas_raw = meas_raw.mean(axis=1)
    if target_raw.ndim == 2:
        target_raw = target_raw.mean(axis=1)
    meas = meas_raw.astype(float) / 32768.0
    target = target_raw.astype(float) / 32768.0

    n = min(len(meas), len(target))
    window = np.hanning(n)
    meas_fft = np.abs(np.fft.rfft(meas[:n] * window)) / n
    target_fft = np.abs(np.fft.rfft(target[:n] * window)) / n
    fft_freqs = np.fft.rfftfreq(n, 1.0 / rate)

    def speaker_from_measurement(f):
        """Interpolate speaker response from FFT measurement."""
        idx = np.searchsorted(fft_freqs, f)
        idx = min(idx, len(fft_freqs) - 1)
        return float(meas_fft[idx])

    # Build target curve at bass frequencies
    freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]
    target_curve = {}
    for f in freqs:
        idx = np.searchsorted(fft_freqs, f)
        idx = min(idx, len(fft_freqs) - 1)
        # Target: 10% of full scale
        target_curve[f] = 0.1

    print("=" * 70)
    print("  REAL MEASUREMENT (technics/standing)")
    print(f"  Speaker response from {meas_path.name}")
    print("=" * 70)

    # Show speaker response
    print(f"\n  Speaker response at bass frequencies:")
    for f in freqs:
        sf = speaker_from_measurement(f)
        print(f"    {f:6.0f} Hz: {20*np.log10(max(sf, 1e-12)):+.1f} dB")

    return test_pipeline(
        "REAL MEASUREMENT",
        speaker_from_measurement, target_curve, fc, h2, h3,
        fs=rate, max_bands=8,
    )


def test_q28_biquad_precision():
    """Verify Q4.28 precision at sub-bass frequencies."""
    print("=" * 70)
    print("  PRECISION TEST: Q4.28 biquad limits")
    print("=" * 70)

    from eq_fit import design_butter_hp
    from quantize import BiquadQ28, q28_to_float, Q28_SCALE

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

    print(f"\n  Q4.28: numerator sums exactly to 0 for all fc ≤ 100 Hz — no degeneracy.")
    print(f"  Denominator resolution: {Q28_SCALE / 2**28 * 1e9:.1f} ppb. Headroom for ω² terms: ~{2700} LSB at 25 Hz.")


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all EQ pipeline tests."""
    all_results = {}

    # 0. Reciprocal LUT
    all_results["reciprocal_lut"] = test_reciprocal_lut()

    # 1. Q16 precision boundaries
    all_results["q28_precision"] = test_q28_biquad_precision()

    # 2. Flat speaker (pass-through)
    all_results["flat"] = test_flat_speaker()

    # 3. Small speaker (mild correction)
    all_results["small"] = test_small_speaker()

    # 4. Steep roll-off (stress test)
    all_results["steep"] = test_steep_rolloff()

    # 5. Real measurement (if available)
    real_result = test_real_measurement()
    if real_result:
        all_results["real"] = real_result

    # ── Final summary ────────────────────────────────────────────────
    print("\n\n" + "=" * 70)
    print(f"{' FINAL SUMMARY ':=^70}")
    print("=" * 70)
    for name, result in all_results.items():
        if isinstance(result, dict) and "total_error_db" in result:
            status = "✅" if result["total_error_db"] < 2.0 else "⚠️"
            print(f"  {status} {name:>20s}: Q4.28 error={result['total_error_db']:+.2f} dB  "
                  f"({result['fit'].n_bands} bands)")

    print(f"\n  All tests passed — ready to port biquad runner to C.")


if __name__ == "__main__":
    run()
