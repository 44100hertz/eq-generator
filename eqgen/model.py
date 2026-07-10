"""
Frequency-domain model for the harmonic bass enhancer.

Predicts the enhancer's RMS output for a sine at frequency f with EQ gain G.
The harmonic terms use h2⁴ / h3⁶ to match the Chebyshev polynomial generator:
  T₂(x) = 2x² − 1  →  2nd harmonic amplitude ∝ h2_amp²  (power ∝ h2⁴)
  T₃(x) = 4x³ − 3x →  3rd harmonic amplitude ∝ h3_amp³  (power ∝ h3⁶)

  A(f, G) = G · √[ HP(f,fc)² · S(f)² · (1 + 3·h3·(h3²−1)·LP(f,fc/2))²
                   + h2⁴ · LP(f,fc)²  · HP(2f,fc)² · S(2f)²
                   + h3⁶ · LP(f,fc/2)² · HP(3f,fc)² · S(3f)² ]

  The first term accounts for T₃ fundamental leakage: T₃(x) = 4x³−3x
  produces a fundamental component 3h3(h3²−1)·sin(ωt) that is coherent
  with the dry fundamental.  When h3 < 1 this reduces the perceived
  fundamental level near the crossover.

  LP filters are first-order (matching the C enhancer implementation).

Where S(f) is the speaker's frequency response.

Solving A(f, G) = target(f) for G gives the EQ gain that, after the
enhancer, produces the desired output level.
"""

import numpy as np

from eqgen.dsp import first_order_lp_mag, butterworth_hp_mag


# ── Core model ────────────────────────────────────────────────────────

# Equal-loudness: at bass frequencies the ear is ~5× more sensitive
# to harmonics (100-150 Hz) than fundamentals (50 Hz). Weight the
# harmonic terms so the model doesn't over-boost.
EAR_W = 5.0
EAR_W2 = EAR_W * EAR_W

def model_perceived_amplitude(
    f: float,
    G: float,
    fc: float,
    h2: float,
    h3: float,
    speaker_response: callable = None,
) -> float:
    """RMS amplitude the model predicts at the output for a sine at f with gain G.

    If speaker_response is None, assumes flat S(f)=1.
    """
    S = speaker_response(f) if speaker_response else 1.0
    S2 = speaker_response(2*f) if speaker_response else 1.0
    S3 = speaker_response(3*f) if speaker_response else 1.0

    hp_f  = butterworth_hp_mag(f, fc)
    hp_2f = butterworth_hp_mag(2*f, fc)
    hp_3f = butterworth_hp_mag(3*f, fc)
    lp_f  = first_order_lp_mag(f, fc)
    lp_f2 = first_order_lp_mag(f, fc / 2.0)

    # T₃ fundamental leakage: 3h3(h3²−1) is negative for h3<1, so the
    # coherent combination with the dry fundamental REDUCES the level.
    fund_factor = 1.0 + 3.0 * h3 * (h3 * h3 - 1.0) * lp_f2

    a = (hp_f * S * fund_factor)**2
    b = h2**4 * lp_f**2  * hp_2f**2 * S2**2
    c = h3**6 * lp_f2**2 * hp_3f**2 * S3**2

    return G * np.sqrt(a + EAR_W2 * b + EAR_W2 * c)


def model_gain_needed(
    f: float,
    target_amplitude: float,
    fc: float,
    h2: float,
    h3: float,
    speaker_response: callable = None,
) -> float:
    """The gain G needed so that model output = target_amplitude at frequency f."""
    S = speaker_response(f) if speaker_response else 1.0

    hp_f  = butterworth_hp_mag(f, fc)
    hp_2f = butterworth_hp_mag(2*f, fc)
    hp_3f = butterworth_hp_mag(3*f, fc)
    lp_f  = first_order_lp_mag(f, fc)
    lp_f2 = first_order_lp_mag(f, fc / 2.0)

    S2 = speaker_response(2*f) if speaker_response else 1.0
    S3 = speaker_response(3*f) if speaker_response else 1.0

    fund_factor = 1.0 + 3.0 * h3 * (h3 * h3 - 1.0) * lp_f2

    a = (hp_f * S * fund_factor)**2
    b = h2**4 * lp_f**2  * hp_2f**2 * S2**2
    c = h3**6 * lp_f2**2 * hp_3f**2 * S3**2

    denom = np.sqrt(a + EAR_W2 * b + EAR_W2 * c)
    if denom < 1e-12:
        return 0.0
    return target_amplitude / denom


# ── Model gain (legacy — used by loudness analysis) ───────────────────

def model_gain(f, fc, h2, h3, S, weighted=False):
    """Compute model gain G for frequency f (with optional ear weighting)."""
    Sf = S(f)
    S2f = S(2 * f)
    S3f = S(3 * f)

    hp = butterworth_hp_mag
    lp = first_order_lp_mag

    fund_factor = 1.0 + 3.0 * h3 * (h3 * h3 - 1.0) * lp(f, fc / 2.0)

    a = (hp(f, fc) * Sf * fund_factor)**2
    b = h2**4 * lp(f, fc)**2 * hp(2 * f, fc)**2 * S2f**2
    c = h3**6 * lp(f, fc / 2.0)**2 * hp(3 * f, fc)**2 * S3f**2

    if weighted:
        sf = ear_sensitivity(f)
        s2f = ear_sensitivity(2 * f)
        s3f = ear_sensitivity(3 * f)
        denom = np.sqrt(sf**2 * a + s2f**2 * b + s3f**2 * c)
    else:
        denom = np.sqrt(a + b + c)

    if denom < 1e-12:
        return 1.0

    target = 1.0
    G = target / denom
    comp = target / Sf if Sf > 1e-12 else 1.0
    return min(G, comp)


# ── Equal-loudness sensitivity (ISO 226:2023, ~60 phon) ───────────────

def ear_sensitivity(f: float) -> float:
    """Ear sensitivity at frequency f, relative to 1 kHz = 1.0.

    Values < 1 mean the ear is LESS sensitive (needs more SPL).
    Approximate fit to ISO 226:2023 at ~60 phon.
    """
    if f <= 0:
        return 0.0
    # Key data points at 60 phon:
    # 20Hz: 98dB → s = 10^((60-98)/20) = 0.013
    # 30Hz: 88dB → s = 0.040
    # 50Hz: 78dB → s = 0.126
    # 100Hz: 65dB → s = 0.562
    # 200Hz: 61dB → s = 0.891
    # 1kHz: 60dB → s = 1.000
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
            log_f1, log_f2 = np.log10(f1), np.log10(f2)
            log_s1, log_s2 = np.log10(s1), np.log10(s2)
            t = (np.log10(f) - log_f1) / (log_f2 - log_f1)
            log_s = log_s1 + t * (log_s2 - log_s1)
            return 10.0 ** log_s
    return 1.0


# ── Example speaker models ────────────────────────────────────────────

def flat_speaker(f: float) -> float:
    """Flat speaker response (identity)."""
    return 1.0


def small_speaker(f: float) -> float:
    """Small speaker: -12 dB at 50 Hz, flat above 100 Hz."""
    if f <= 50:
        return 0.25
    elif f >= 100:
        return 1.0
    else:
        return 0.25 + 0.75 * (f - 50) / 50


# ── Full preprocessing curve (mirrors Rust bass_enhancer_preprocess) ─

def preprocess_eq_curve(
    target_curve: dict,       # {freq: target_gain}
    speaker_response: callable,
    fc: float,
    h2: float,
    h3: float,
) -> dict:
    """Compute preprocessed EQ gains matching the Rust pipeline.

    For each frequency in target_curve, computes G such that:
        A(f, G) = target_curve[f]
    where A is model_perceived_amplitude with the given speaker.

    Below fc/2 the enhancer handles bass perception via harmonics —
    EQ should stay flat there.
    """
    result = {}

    for f, target in sorted(target_curve.items()):
        S_f = speaker_response(f)
        model_factor = model_perceived_amplitude(f, 1.0, fc, h2, h3, speaker_response)
        G = target / model_factor if model_factor > 1e-12 else 0.01

        bare_correction = target / S_f if S_f > 1e-12 else 0.01
        G = min(G, bare_correction)

        result[f] = G

    # Below fc/2: hold correction flat at the model gain computed AT fc/2.
    freqs_sorted = sorted(target_curve.keys())
    target_at_fc2 = float(np.interp(fc / 2.0, freqs_sorted,
                                     [target_curve[f] for f in freqs_sorted]))
    S_fc2 = speaker_response(fc / 2.0)
    model_g = model_perceived_amplitude(fc / 2.0, 1.0, fc, h2, h3, speaker_response)
    flat_val = target_at_fc2 / model_g if model_g > 1e-12 else 0.01
    flat_val = min(flat_val, target_at_fc2 / S_fc2 if S_fc2 > 1e-12 else 0.01)
    for f in result:
        if f <= fc / 2.0:
            result[f] = flat_val

    return dict(sorted(result.items()))


# ── Analysis runs ─────────────────────────────────────────────────────

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
