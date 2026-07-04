"""
Verify the Rust bass_enhancer_preprocess model against actual v2 output.

The model says:
  A(f, G) = G · √[ HP(f,fc)²·S(f)²
                  + h2²·LP(f,fc)²  ·HP(2f,fc)²·S(2f)²
                  + h3²·LP(f,fc/2)²·HP(3f,fc)²·S(3f)² ]

We test: feed a sine at frequency f through the actual v2 plugin,
measure the RMS output, and compare against what the model predicts.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    design_butter_lp,
    design_butter_hp,
    generate_test_sine,
)
from fix_v2 import process_fixed_v2
from goertzel_analyzer import goertzel_magnitude


def butterworth_lp_mag(f: float, fc: float) -> float:
    w = f / fc
    return 1.0 / np.sqrt(1.0 + w**4)


def butterworth_hp_mag(f: float, fc: float) -> float:
    w = f / fc
    return w * w / np.sqrt(1.0 + w**4)


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
    lp_f  = butterworth_lp_mag(f, fc)
    lp_f2 = butterworth_lp_mag(f, fc / 2.0)

    a = (hp_f * S)**2
    b = h2**2 * lp_f**2  * hp_2f**2 * S2**2
    c = h3**2 * lp_f2**2 * hp_3f**2 * S3**2

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
    S2 = speaker_response(2*f) if speaker_response else 1.0
    S3 = speaker_response(3*f) if speaker_response else 1.0

    hp_f  = butterworth_hp_mag(f, fc)
    hp_2f = butterworth_hp_mag(2*f, fc)
    hp_3f = butterworth_hp_mag(3*f, fc)
    lp_f  = butterworth_lp_mag(f, fc)
    lp_f2 = butterworth_lp_mag(f, fc / 2.0)

    a = (hp_f * S)**2
    b = h2**2 * lp_f**2  * hp_2f**2 * S2**2
    c = h3**2 * lp_f2**2 * hp_3f**2 * S3**2

    denom = np.sqrt(a + b + c)
    if denom < 1e-12:
        return 0.0
    return target_amplitude / denom


def verify_model():
    """Compare model predictions against actual v2 plugin output."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.33, h3_amp=0.33,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  MODEL VERIFICATION: predicted vs actual v2 output")
    print(f"  h2={cfg.h2_amp}, h3={cfg.h3_amp}, cutoff={cfg.cutoff} Hz")
    print("=" * 70)

    freqs = [30, 40, 50, 60, 80, 100, 150, 200]

    # For each frequency: run v2 on a sine at amplitude=0.5,
    # measure the total RMS and the fundamental amplitude.
    # The model predicts total RMS = model_perceived_amplitude(f, 0.5, ...)

    # The model assumes the plugin preserves amplitude linearly, so:
    # output_RMS ≈ input_amplitude * model_perceived_amplitude(f, 1.0, ...)
    # But the input goes through HP too, so we need to be careful.

    # The model says: with gain G=1.0, perceived amplitude at ear is
    # sqrt(a+b+c) where each term is the contribution at each frequency.
    # The actual total RMS of the output signal should be
    # (A/√2) * sqrt(a+b+c) because each sine component has RMS = amplitude/√2
    # and they're at different frequencies (power sum).

    # So: output_RMS = (A/√2) * sqrt(a + b + c)
    # where A is the input amplitude

    print(f"\n  {'Freq':>6s}  {'Actual RMS':>10s}  {'Model RMS':>10s}  "
          f"{'Ratio':>8s}  {'Error dB':>10s}  {'F pred':>10s}  {'F actual':>10s}")
    print(f"  {'-'*72}")

    for freq in freqs:
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        A = 0.5
        signal = generate_test_sine([freq], [A], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out = process_fixed_v2(stereo, cfg)
        ss = slice(pad + len(signal)//2, pad + len(signal))
        mono = (out[0, ss] + out[1, ss]) / 2.0

        actual_rms = np.sqrt(np.mean(mono**2))

        # Model prediction
        model_factor = model_perceived_amplitude(freq, 1.0, cfg.cutoff,
                                                  cfg.h2_amp, cfg.h3_amp)
        # For the model, the output is G * sqrt(a+b+c) where each
        # term in a,b,c already has the amplitude factored in.
        # The actual per-frequency amplitude from the plugin for input A·sin(2πft):
        #   Fundamental: A · HP(f,fc)
        #   H2: h2 · A · LP(f,fc) · HP(2f,fc)
        #   H3: h3 · A · LP(f,fc/2) · HP(3f,fc)
        # These are at different frequencies, so total power is sum of powers.
        # RMS = sqrt(fund²/2 + h2²/2 + h3²/2) = (1/√2)·sqrt(fund² + h2² + h3²)
        #      = (A/√2)·sqrt(HP²·S² + h2²·LP²·HP²(2f)·S²(2f) + h3²·LP²(fc/2)·HP²(3f)·S²(3f))
        #      = (A/√2)·model_factor
        # With S(f)=1 (flat speaker):
        predicted_rms = (A / np.sqrt(2.0)) * model_factor

        ratio = actual_rms / predicted_rms if predicted_rms > 0 else 0
        error_db = 20 * np.log10(ratio) if ratio > 0 else float('-inf')

        # Also check fundamental amplitude specifically
        fund_actual = goertzel_magnitude(mono, freq, cfg.fs)
        hp_f = butterworth_hp_mag(freq, cfg.cutoff)
        fund_predicted = A * hp_f

        print(f"  {freq:6.0f}  {actual_rms:10.6f}  {predicted_rms:10.6f}  "
              f"{ratio:8.4f}  {error_db:+10.2f}  {fund_predicted:10.6f}  {fund_actual:10.6f}")

    # ── Test the gain solver ──────────────────────────────────────────
    print(f"\n  ── Gain solver test ──")
    print(f"  For each frequency, compute the G needed so that model")
    print(f"  output at f equals a target amplitude of 0.5.")
    print(f"  Then feed that G into the plugin and verify output ≈ target.")
    print(f"")
    print(f"  {'Freq':>6s}  {'G needed':>10s}  {'Actual RMS':>10s}  "
          f"{'Target':>8s}  {'Error dB':>10s}")
    print(f"  {'-'*52}")

    for freq in freqs:
        target = 0.5
        # The model gain solver returns G such that model_perceived_amplitude(f, G, ...) = target
        # But model_perceived_amplitude returns RMS-like value (already divided by √2)
        # Actually no, it returns G * sqrt(a+b+c). For A·sin input with G=1:
        # output RMS = (A/√2) * sqrt(a+b+c). For arbitrary input gain G:
        # output RMS = (G/√2) * sqrt(a+b+c) = model_perceived_amplitude(f, G, ...)/√2
        # Wait, model_perceived_amplitude returns G*sqrt(a+b+c), which is the PEAK amplitude sum.
        # The actual RMS is model_perceived_amplitude(f, G, ...)/√2.

        # So we want output RMS = target, which means:
        # G * sqrt(a+b+c) / √2 = target
        # G = target * √2 / sqrt(a+b+c)

        model_factor = model_perceived_amplitude(freq, 1.0, cfg.cutoff,
                                                  cfg.h2_amp, cfg.h3_amp)
        G = target * np.sqrt(2.0) / model_factor if model_factor > 1e-12 else 0.0

        # Now run v2 with input amplitude = G (this simulates the EQ applying gain G)
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [G], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out = process_fixed_v2(stereo, cfg)
        ss = slice(pad + len(signal)//2, pad + len(signal))
        mono = (out[0, ss] + out[1, ss]) / 2.0

        actual_rms = np.sqrt(np.mean(mono**2))
        error_db = 20 * np.log10(actual_rms / target) if target > 0 else 0

        print(f"  {freq:6.0f}  {G:10.4f}  {actual_rms:10.6f}  "
              f"{target:8.4f}  {error_db:+10.2f}")


def run():
    """Run the model verification."""
    verify_model()


if __name__ == "__main__":
    run()
