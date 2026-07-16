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
#
# NOTE: No additional ear-weighting constant here.
# compute_harmonic_efficacy() builds Ear(f)/Ear(nf) directly into
# h2/h3 via:  h = Correction(nf)/Correction(f) * Ear(f)/Ear(nf).
# The crossfade in enhancer.c then uses h as perceptual efficiency.
# An extra EAR_W on the harmonic term would double-count the ear.
#
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

    return G * np.sqrt(a + b + c)


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

    denom = np.sqrt(a + b + c)
    if denom < 1e-12:
        return 0.0
    return target_amplitude / denom


# ── Equal-loudness sensitivity (ISO 226:2023, ~60 phon) ───────────────

def ear_sensitivity(f: float) -> float:
    """Ear sensitivity at frequency f, relative to 1 kHz = 1.0.

    Values < 1 mean the ear is LESS sensitive (needs more SPL).
    Based on ISO 226:2023 contours from iso226_table.csv.
    """
    from eqgen.iso226 import sensitivity as _iso226_sensitivity
    return _iso226_sensitivity(f, phon=60.0)


# ── Perceptual weighting (ear × music spectrum) ───────────────────

# Long-term average music spectrum slope: ~-4.5 dB/octave relative to 1 kHz.
# This is the approximate tilt observed across large mixed-music corpora —
# slightly steeper than pink noise (-3 dB/oct) because instrument fundamentals
# cluster at lower frequencies while harmonics roll off.
_MUSIC_TILT_DB_PER_OCTAVE = 4.5


def music_spectrum_weight(f: float) -> float:
    """Relative energy in average music at frequency f (1 kHz = 1.0).

    Models the long-term average spectrum of mixed music as a simple
    log-linear slope.  Values > 1 mean more energy than at 1 kHz.
    """
    octaves = np.log2(max(f, 1.0) / 1000.0)
    db = -_MUSIC_TILT_DB_PER_OCTAVE * octaves
    return float(10.0 ** (db / 20.0))


def perceptual_weight(f: float) -> float:
    """Combined perceptual importance of frequency f.

    ear_sensitivity(f) × music_spectrum_weight(f).

    Uses 80 phon (typical listening level) rather than 60 phon — the
    equal-loudness contours flatten significantly at higher SPL, and
    music is rarely auditioned near threshold.

    Peak is ~500 Hz where both curves cross; rolls off at extremes.
    """
    from eqgen.iso226 import sensitivity as _iso226_sensitivity
    ear = _iso226_sensitivity(f, phon=80.0)
    return ear * music_spectrum_weight(f)


# ── Harmonic efficacy: compute h2/h3 from measurement data ────

def compute_harmonic_efficacy(
    freqs: np.ndarray,
    correction_db: np.ndarray,
    fc: float,
) -> dict:
    """Compute h2_amp, h3_amp from speaker measurement data.

    For each harmonic n (2 or 3), average over the bass band (f < fc/2):

        h_amp = avg( Correction(n·f) / Correction(f) · Ear(f) / Ear(n·f) )

    Correction(f) = 10^(correction_db/20) — EQ gain at f
    Ear(f) = ear sensitivity at f

    When the speaker rolls off sharply, Correction(f) ≫ Correction(n·f)
    → denominator dominates → h is small → harmonics are perceptually less
    needed (the fundamental gets through fine at n·f frequencies).

    When the speaker has flat bass, Correction(f) ≈ Correction(n·f)
    → ratio ≈ Ear(f)/Ear(n·f) ≈ 0.35 — pure psychoacoustic rolloff.
    """
    correction_ratio = 10.0 ** (correction_db / 20.0)

    mask = freqs <= fc / 2.0
    if not mask.any():
        return {"h2_amp": 0.0, "h3_amp": 0.0}

    bass_f = freqs[mask]

    def avg_ratio(mult):
        vals = []
        for f in bass_f:
            hf = f * mult
            if hf > freqs[-1]:
                continue
            cr_f  = np.interp(f,  freqs, correction_ratio)
            cr_hf = np.interp(hf, freqs, correction_ratio)
            ear_f  = ear_sensitivity(f)
            ear_hf = ear_sensitivity(hf)
            vals.append(cr_hf / cr_f * ear_f / ear_hf)
        return float(np.mean(vals)) if vals else 0.0

    h2 = avg_ratio(2)
    h3 = avg_ratio(3)

    return {
        "h2_amp": round(float(np.clip(h2, 0.0, 0.5)), 4),
        "h3_amp": round(float(np.clip(h3, 0.0, 0.5)), 4),
    }


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



