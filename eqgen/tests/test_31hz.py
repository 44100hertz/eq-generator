"""
End-to-end test: 31 Hz tone through the C enhancer with bass-active EQ.

31 Hz was historically problematic — high EQ correction combined with
overboost would drive the enhancer into hard clipping.  This test runs
the full pipeline (EQ + C enhancer) at multiple amplitudes and verifies:
  - No clipping at any tested level
  - THD stays below 100% (fundamental dominates)
  - Output is monotonic with input (quieter input → quieter output)
"""

import numpy as np

from eqgen.sweep import run_sine_sweep
from eqgen.dsp import pre_gain_from_max_gain


def _make_boost_eq(fs: float = 44100.0):
    """Make a flat EQ with a strong bass shelf — stresses the enhancer at 31 Hz.

    Returns (eq_coeffs_flat, pre_gain).
    """
    from eqgen.eq_fit import BiquadCoeffs, design_low_shelf

    # 12 dB low shelf at 80 Hz — aggressive correction for a bass-weak speaker
    shelf = design_low_shelf(80.0, 12.0, 0.5, fs)
    coeffs = [shelf.b0, shelf.b1, shelf.b2, shelf.a1, shelf.a2]
    # Pre-gain compensates the peak: 12 dB boost → pre_gain = 10^(-12/20) ≈ 0.25
    pre_gain = float(pre_gain_from_max_gain(12.0))
    return coeffs, pre_gain


def test_31hz_no_clipping():
    """31 Hz tone at multiple amplitudes must not clip."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    for amp_db in [-40, -20, -12, -6, -3, 0]:
        amp_lin = 10.0 ** (amp_db / 20.0)
        results = run_sine_sweep(
            freqs_hz=[31.0],
            eq_coeffs=coeffs,
            fc=60.0, h2=0.5, h3=1.0,
            amplitude=amp_lin,
            duration_sec=1.0,
            pre_gain=pre_gain,
        )
        r = results[31.0]
        # RMS of a full-scale sine without DC offset is ~0.707.
        # After enhancer (harmonics + tanh), RMS shouldn't exceed 0.95.
        assert r["rms"] <= 0.95, \
            f"31 Hz at {amp_db:+d} dBFS: RMS={r['rms']:.4f} exceeds 0.95 — likely clipping"


def test_31hz_reasonable_thd():
    """31 Hz tone should have THD < 100% at moderate levels."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    # Test at -12 dBFS — loud enough to measure but below clipping
    amp_lin = 10.0 ** (-12.0 / 20.0)
    results = run_sine_sweep(
        freqs_hz=[31.0],
        eq_coeffs=coeffs,
        fc=60.0, h2=0.5, h3=1.0,
        amplitude=amp_lin,
        duration_sec=1.0,
        pre_gain=pre_gain,
    )
    r = results[31.0]
    fund = r["fundamental"]
    h2 = r["h2"]
    h3 = r["h3"]

    # With h2=0.5, h3=1.0, the enhancer deliberately synthesizes harmonics.
    # THD can be high but should be < 100% (fundamental still dominates).
    thd = np.sqrt(h2 ** 2 + h3 ** 2) / max(fund, 1e-12) * 100.0
    assert thd < 100.0, \
        f"31 Hz THD={thd:.1f}% exceeds 100% — harmonics louder than fundamental"


def test_31hz_output_monotonic():
    """Quieter input at 31 Hz should produce quieter output."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    rms_values = []
    for amp_db in [-40, -30, -20, -12]:
        amp_lin = 10.0 ** (amp_db / 20.0)
        results = run_sine_sweep(
            freqs_hz=[31.0],
            eq_coeffs=coeffs,
            fc=60.0, h2=0.5, h3=1.0,
            amplitude=amp_lin,
            duration_sec=1.0,
            pre_gain=pre_gain,
        )
        rms_values.append(results[31.0]["rms"])

    # Output should strictly increase with input level
    for i in range(1, len(rms_values)):
        assert rms_values[i] > rms_values[i - 1], \
            f"Output not monotonic: {rms_values[i-1]:.6f} → {rms_values[i]:.6f} at step {i}"


def run():
    test_31hz_no_clipping()
    test_31hz_reasonable_thd()
    test_31hz_output_monotonic()
    print("  ✅ test_31hz passed")


if __name__ == "__main__":
    run()
