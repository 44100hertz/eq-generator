"""
End-to-end input signal test: problematic tones through the C enhancer with bass EQ.

Runs sine tones at multiple frequencies through the full pipeline (EQ + C enhancer)
and verifies:
  - No clipping at any tested level
  - THD stays below 100% (fundamental dominates)
  - Output is monotonic with input (quieter input → quieter output)

31 Hz was historically problematic — high EQ correction with overboost drove
the enhancer into hard clipping.  This test generalizes that check to other
bass frequencies.
"""

import numpy as np

from eqgen.sweep import run_sine_sweep
from eqgen.dsp import pre_gain_from_max_gain

# Frequencies to test — all bass tones that stress the enhancer crossover
TEST_FREQS = [31.0, 60.0, 100.0]


def _make_boost_eq(fs: float = 44100.0):
    """12 dB low shelf at 80 Hz — simulates aggressive bass correction."""
    from eqgen.eq_fit import design_low_shelf

    shelf = design_low_shelf(80.0, 12.0, 0.5, fs)
    coeffs = [shelf.b0, shelf.b1, shelf.b2, shelf.a1, shelf.a2]
    pre_gain = float(pre_gain_from_max_gain(12.0))
    return coeffs, pre_gain


def _test_no_clipping(freq: float):
    """Tone at multiple amplitudes must not clip."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    for amp_db in [-40, -20, -12, -6, -3, 0]:
        amp_lin = 10.0 ** (amp_db / 20.0)
        results = run_sine_sweep(
            freqs_hz=[freq],
            eq_coeffs=coeffs,
            fc=60.0, h2=0.5, h3=1.0,
            amplitude=amp_lin,
            duration_sec=1.0,
            pre_gain=pre_gain,
        )
        r = results[freq]
        assert r["rms"] <= 0.95, \
            f"{freq:.0f} Hz at {amp_db:+d} dBFS: RMS={r['rms']:.4f} exceeds 0.95 — likely clipping"


def _test_reasonable_thd(freq: float):
    """Tone should have THD < 100% at -12 dBFS."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    amp_lin = 10.0 ** (-12.0 / 20.0)
    results = run_sine_sweep(
        freqs_hz=[freq],
        eq_coeffs=coeffs,
        fc=60.0, h2=0.5, h3=1.0,
        amplitude=amp_lin,
        duration_sec=1.0,
        pre_gain=pre_gain,
    )
    r = results[freq]
    fund = r["fundamental"]
    h2 = r["h2"]
    h3 = r["h3"]

    thd = np.sqrt(h2 ** 2 + h3 ** 2) / max(fund, 1e-12) * 100.0
    assert thd < 100.0, \
        f"{freq:.0f} Hz THD={thd:.1f}% exceeds 100% — harmonics louder than fundamental"


def _test_output_monotonic(freq: float):
    """Quieter input should produce quieter output."""
    fs = 44100.0
    coeffs, pre_gain = _make_boost_eq(fs)

    rms_values = []
    for amp_db in [-40, -30, -20, -12]:
        amp_lin = 10.0 ** (amp_db / 20.0)
        results = run_sine_sweep(
            freqs_hz=[freq],
            eq_coeffs=coeffs,
            fc=60.0, h2=0.5, h3=1.0,
            amplitude=amp_lin,
            duration_sec=1.0,
            pre_gain=pre_gain,
        )
        rms_values.append(results[freq]["rms"])

    for i in range(1, len(rms_values)):
        assert rms_values[i] > rms_values[i - 1], \
            f"{freq:.0f} Hz output not monotonic: {rms_values[i-1]:.6f} → {rms_values[i]:.6f} at step {i}"


def run():
    for freq in TEST_FREQS:
        _test_no_clipping(freq)
        _test_reasonable_thd(freq)
        _test_output_monotonic(freq)
    print(f"  ✅ test_input_signals passed ({len(TEST_FREQS)} frequencies)")


if __name__ == "__main__":
    run()
