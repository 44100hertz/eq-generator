"""
Table-generating analysis tools for the harmonic bass enhancer model.

These were extracted from eqgen.model to keep production code clean.
They are standalone analysis scripts, not part of the EQ pipeline.

Usage:
    python -c "from eqgen.analysis_tools.loudness import run_loudness_analysis; run_loudness_analysis()"
"""

import numpy as np

from eqgen.model import (
    model_gain,
    ear_sensitivity,
    small_speaker,
    butterworth_hp_mag,
    first_order_lp_mag,
)


def run_loudness_analysis():
    """Equal-loudness weighting effect on model gain."""
    fc = 60.0
    h2, h3 = 0.33, 0.33

    print("=" * 65)
    print("  EQUAL-LOUDNESS WEIGHTING EFFECT ON MODEL GAIN")
    print("  Speaker: -12 dB @ 50 Hz, flat > 100 Hz")
    print(f"  h2={h2}, h3={h3}, cutoff={fc} Hz")
    print("=" * 65)

    print(f"\n  {'Freq':>6s}  {'s(f)':>8s}  {'s(2f)':>8s}  {'s(3f)':>8s}  "
          f"{'Unweighted':>12s}  {'Weighted':>12s}  {'Delta':>8s}")
    print(f"  {'-'*64}")

    for f in [25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]:
        G_unw = model_gain(f, fc, h2, h3, small_speaker, weighted=False)
        G_w = model_gain(f, fc, h2, h3, small_speaker, weighted=True)
        s_f = ear_sensitivity(f)
        s_2f = ear_sensitivity(2 * f)
        s_3f = ear_sensitivity(3 * f)
        print(f"  {f:5.0f} Hz  {s_f:8.3f}  {s_2f:8.3f}  {s_3f:8.3f}  "
              f"{20*np.log10(G_unw):+11.2f} dB  "
              f"{20*np.log10(G_w):+11.2f} dB  "
              f"{20*np.log10(G_w/G_unw):+7.2f} dB")

    # Ear sensitivity curve
    print(f"\n  ── Ear sensitivity curve ──")
    for f in [20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]:
        s = ear_sensitivity(f)
        bar = '█' * int(s * 40)
        print(f"  {f:5.0f} Hz:  {s:.3f}  {bar}")


def run_model_gain_analysis():
    """Analyze how h2/h3 changes affect model gain."""
    fc = 60.0
    freqs = [30, 40, 50, 60, 80, 100, 120]

    print("Speaker: -12 dB at 50 Hz, flat > 100 Hz")
    print(f"Cutoff: {fc} Hz")
    print()
    print(f"{'Freq':>6s}  {'Old G':>8s} {'Old dB':>8s}  {'New G':>8s} {'New dB':>8s}  {'Delta':>8s}")
    print(f"{'':->50s}")

    for f in freqs:
        Sf = small_speaker(f)
        S2f = small_speaker(2*f) if 2*f < 20000 else 1.0
        S3f = small_speaker(3*f) if 3*f < 20000 else 1.0

        def _gain(h2, h3):
            hp = butterworth_hp_mag
            lp = first_order_lp_mag
            fund_factor = 1.0 + 3.0 * h3 * (h3 * h3 - 1.0) * lp(f, fc / 2.0)
            a = (hp(f, fc) * Sf * fund_factor)**2
            b = h2**4 * lp(f, fc)**2 * hp(2 * f, fc)**2 * S2f**2
            c = h3**6 * lp(f, fc / 2.0)**2 * hp(3 * f, fc)**2 * S3f**2
            G = 1.0 / np.sqrt(a + b + c) if a + b + c > 1e-12 else 1.0
            return min(G, 1.0 / Sf if Sf > 1e-12 else 1.0)

        G_old = _gain(1.0, 2.0)
        G_new = _gain(0.33, 0.33)
        dB_old = 20*np.log10(G_old)
        dB_new = 20*np.log10(G_new)
        print(f"  {f:4.0f} Hz  {G_old:8.3f} {dB_old:+7.2f}  {G_new:8.3f} {dB_new:+7.2f}  {dB_new-dB_old:+7.2f}")

    # What h2/h3 give the SAME gains as the old model?
    print(f"\n{'='*50}")
    print("What h2/h3 produce the SAME gain as old model?")
    print("(Old model was tuned with h2=1.0, h3=2.0)")
    print(f"{'='*50}")
    print(f"{'Freq':>6s}  {'Old G':>8s}  ", end="")
    for ratio in [0.33, 0.5, 0.7, 1.0]:
        print(f"{'h='+str(ratio):>8s}", end="  ")
    print()
    for f in [30, 40, 50, 60, 80]:
        Sf = small_speaker(f)
        S2f = small_speaker(2*f) if 2*f < 20000 else 1.0
        S3f = small_speaker(3*f) if 3*f < 20000 else 1.0

        def _gain_h(h2, h3):
            hp = butterworth_hp_mag
            lp = first_order_lp_mag
            fund_factor = 1.0 + 3.0 * h3 * (h3 * h3 - 1.0) * lp(f, fc / 2.0)
            a = (hp(f, fc) * Sf * fund_factor)**2
            b = h2**4 * lp(f, fc)**2 * hp(2 * f, fc)**2 * S2f**2
            c = h3**6 * lp(f, fc / 2.0)**2 * hp(3 * f, fc)**2 * S3f**2
            G = 1.0 / np.sqrt(a + b + c) if a + b + c > 1e-12 else 1.0
            return min(G, 1.0 / Sf if Sf > 1e-12 else 1.0)

        G_old = _gain_h(1.0, 2.0)
        print(f"  {f:4.0f} Hz  {G_old:8.3f}  ", end="")
        for ratio in [0.33, 0.5, 0.7, 1.0]:
            G = _gain_h(ratio*2.0, ratio*2.0)
            print(f"{G:8.3f}  ", end="")
        print()


if __name__ == "__main__":
    print("── Loudness analysis ──")
    run_loudness_analysis()
    print()
    print("── Model gain analysis ──")
    run_model_gain_analysis()
