"""
Equal-loudness weighting for the bass enhancer model.

The ear is much more sensitive at 2f and 3f than at f for bass frequencies.
The model's RMS sum treats all frequencies equally, so it underestimates
the harmonics' contribution to perceived loudness → computes too much gain.

ISO 226:2023 equal-loudness contours at 60 phon (typical listening):
  Freq     SPL for 60 phon    Sensitivity rel 1 kHz
  30 Hz    88 dB              s(30)  = 10^(-28/20) = 0.040
  50 Hz    78 dB              s(50)  = 10^(-18/20) = 0.126
  100 Hz   65 dB              s(100) = 10^(-5/20)  = 0.562
  150 Hz   62 dB              s(150) = 10^(-2/20)  = 0.794
  1 kHz    60 dB              s(1k)  = 1.000

Weighted model:
  L² = G² · [s(f)²·a  +  s(2f)²·b  +  s(3f)²·c]
  G  = target / sqrt(s(f)²·a + s(2f)²·b + s(3f)²·c)

where a,b,c are the physical amplitude terms from before.
"""

import numpy as np

# ── Equal-loudness sensitivity (relative to 1 kHz) at 60 phon ────────
# From ISO 226:2023, approximate fit.
def ear_sensitivity(f: float) -> float:
    """Ear sensitivity at frequency f, relative to 1 kHz = 1.0.
    Values < 1 mean the ear is LESS sensitive (needs more SPL).
    Approximate fit to ISO 226:2023 at ~60 phon."""
    if f <= 0:
        return 0.0
    # Key data points at 60 phon:
    # 20Hz: 98dB → s = 10^((60-98)/20) = 0.013
    # 30Hz: 88dB → s = 0.040
    # 50Hz: 78dB → s = 0.126
    # 100Hz: 65dB → s = 0.562
    # 200Hz: 61dB → s = 0.891
    # 1kHz: 60dB → s = 1.000
    # Piecewise log-linear interpolation
    points = [
        (20, 0.013),
        (30, 0.040),
        (50, 0.126),
        (100, 0.562),
        (200, 0.891),
        (1000, 1.000),
        (20000, 1.000),
    ]
    for i in range(len(points) - 1):
        f1, s1 = points[i]
        f2, s2 = points[i + 1]
        if f1 <= f <= f2:
            # Log-linear: log(s) = linear in log(f)
            log_f1, log_f2 = np.log10(f1), np.log10(f2)
            log_s1, log_s2 = np.log10(s1), np.log10(s2)
            t = (np.log10(f) - log_f1) / (log_f2 - log_f1)
            log_s = log_s1 + t * (log_s2 - log_s1)
            return 10.0 ** log_s
    return 1.0


def hp(f, fc):
    w = f / fc
    return w * w / np.sqrt(1 + w**4)


def lp(f, fc):
    w = f / fc
    return 1.0 / np.sqrt(1 + w**4)


def model_gain(f, fc, h2, h3, S, weighted=False):
    """Compute model gain G for frequency f."""
    Sf = S(f)
    S2f = S(2 * f)
    S3f = S(3 * f)

    a = hp(f, fc)**2 * Sf**2
    b = h2**2 * lp(f, fc)**2 * hp(2 * f, fc)**2 * S2f**2
    c = h3**2 * lp(f, fc / 2.0)**2 * hp(3 * f, fc)**2 * S3f**2

    if weighted:
        sf = ear_sensitivity(f)
        s2f = ear_sensitivity(2 * f)
        s3f = ear_sensitivity(3 * f)
        denom = np.sqrt(sf**2 * a + s2f**2 * b + s3f**2 * c)
    else:
        denom = np.sqrt(a + b + c)

    if denom < 1e-12:
        return 1.0

    target = 1.0  # desired perceived output = 1
    G = target / denom
    comp = target / Sf if Sf > 1e-12 else 1.0
    return min(G, comp)


# ── Speaker model: -12 dB at 50 Hz, flat above 100 Hz ────────────────
def speaker(f):
    if f <= 50:
        return 0.25
    elif f >= 100:
        return 1.0
    else:
        return 0.25 + 0.75 * (f - 50) / 50


def run():
    """Run the equal-loudness weighting analysis."""
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
        G_unw = model_gain(f, fc, h2, h3, speaker, weighted=False)
        G_w = model_gain(f, fc, h2, h3, speaker, weighted=True)
        s_f = ear_sensitivity(f)
        s_2f = ear_sensitivity(2 * f)
        s_3f = ear_sensitivity(3 * f)
        print(f"  {f:5.0f} Hz  {s_f:8.3f}  {s_2f:8.3f}  {s_3f:8.3f}  "
              f"{20*np.log10(G_unw):+11.2f} dB  "
              f"{20*np.log10(G_w):+11.2f} dB  "
              f"{20*np.log10(G_w/G_unw):+7.2f} dB")

    # ── Show the ear sensitivity curve ───────────────────────────────────
    print(f"\n  ── Ear sensitivity curve ──")
    for f in [20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500, 1000]:
        s = ear_sensitivity(f)
        bar = '█' * int(s * 40)
        print(f"  {f:5.0f} Hz:  {s:.3f}  {bar}")


if __name__ == "__main__":
    run()
