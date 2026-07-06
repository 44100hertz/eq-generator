"""
Verify the Rust bass_enhancer_preprocess model against actual DSP output.

The model says:
  A(f, G) = G · √[ HP(f,fc)²·S(f)²
                  + h2²·LP(f,fc)²  ·HP(2f,fc)²·S(2f)²
                  + h3²·LP(f,fc/2)²·HP(3f,fc)²·S(3f)² ]

We test: feed a sine at frequency f through the actual DSP (process_fixed_v2),
measure the RMS output, and compare against what the model predicts.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.dsp import BassEnhancerConfig, process_fixed_v2
from eqgen.analysis import generate_test_sine, goertzel_magnitude
from eqgen.model import (
    butterworth_lp_mag,
    butterworth_hp_mag,
    model_perceived_amplitude,
    model_gain_needed,
    preprocess_eq_curve,
    small_speaker,
)


def verify_model():
    """Compare model predictions against actual process_fixed_v2 output."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=1.0, h3_amp=1.0,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  MODEL VERIFICATION: predicted vs actual v2 output")
    print(f"  h2={cfg.h2_amp}, h3={cfg.h3_amp}, cutoff={cfg.cutoff} Hz")
    print("=" * 70)

    freqs = [30, 40, 50, 60, 80, 100, 150, 200]

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

        model_factor = model_perceived_amplitude(freq, 1.0, cfg.cutoff,
                                                  cfg.h2_amp, cfg.h3_amp)
        predicted_rms = (A / np.sqrt(2.0)) * model_factor

        ratio = actual_rms / predicted_rms if predicted_rms > 0 else 0
        error_db = 20 * np.log10(ratio) if ratio > 0 else float('-inf')

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
        model_factor = model_perceived_amplitude(freq, 1.0, cfg.cutoff,
                                                  cfg.h2_amp, cfg.h3_amp)
        G = target * np.sqrt(2.0) / model_factor if model_factor > 1e-12 else 0.0

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


def verify_full_stack():
    """End-to-end test mirroring the complete Rust pipeline.

    1. Model a speaker with heavy bass roll-off (small_speaker)
    2. Compute preprocessed EQ gains using the bass enhancer model
    3. Apply EQ → enhancer → speaker to sine tones
    4. Verify perceived output RMS matches the target curve
    """
    fc = 60.0
    h2, h3 = 1.0, 1.0
    cfg = BassEnhancerConfig(
        cutoff=fc, h2_amp=h2, h3_amp=h3,
        fs=44100.0, env_release=200.0,
    )

    # ── 1. Create a flat target curve ─────────────────────────────────
    freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120]
    target_curve = {f: 1.0 for f in freqs}

    # ── 2. Run preprocessor (mirrors Rust bass_enhancer_preprocess) ──
    eq_curve = preprocess_eq_curve(target_curve, small_speaker, fc, h2, h3)

    print("=" * 70)
    print("  FULL STACK: EQ preprocess → enhancer → speaker")
    print(f"  Speaker: -12 dB @ 50 Hz, flat > 100 Hz")
    print(f"  Target: flat, fc={fc} Hz, h2={h2}, h3={h3}")
    print("=" * 70)

    # Show the computed EQ curve
    print(f"\n  ── Computed EQ curve ──")
    print(f"  {'Freq':>6s}  {'Speaker':>8s}  {'EQ gain':>8s}  "
          f"{'Bare corr':>10s}")
    print(f"  {'-'*38}")
    for f in sorted(eq_curve.keys()):
        S_db = 20 * np.log10(small_speaker(f))
        G_db = 20 * np.log10(eq_curve[f])
        bare_db = 20 * np.log10(1.0 / max(small_speaker(f), 1e-12))
        print(f"  {f:6.0f}  {S_db:+8.1f}  {G_db:+8.2f}  {bare_db:+10.2f}")

    # ── 3. Verify: feed EQ'd sines through the full pipeline ──────────
    print(f"\n  ── Pipeline verification ──")
    print(f"  {'Freq':>6s}  {'EQ (dB)':>8s}  {'Spkr RMS':>10s}  "
          f"{'Target':>8s}  {'Error dB':>10s}  {'Note':>20s}")
    print(f"  {'-'*64}")

    errors = []
    third = fc / 3.0
    for freq in sorted(eq_curve.keys()):
        G = eq_curve[freq]
        target = target_curve.get(freq, 1.0)

        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [G], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out = process_fixed_v2(stereo, cfg)
        ss = slice(pad + len(signal) // 2, pad + len(signal))
        mono = (out[0, ss] + out[1, ss]) / 2.0

        h1 = goertzel_magnitude(mono, freq, cfg.fs)
        h2 = goertzel_magnitude(mono, 2 * freq, cfg.fs)
        h3 = goertzel_magnitude(mono, 3 * freq, cfg.fs)

        Sf  = small_speaker(freq)
        S2f = small_speaker(2 * freq)
        S3f = small_speaker(3 * freq)

        perceived_rms = np.sqrt(
            (h1 * Sf)**2 / 2 +
            (h2 * S2f)**2 / 2 +
            (h3 * S3f)**2 / 2
        )
        perceived_db = 20 * np.log10(max(perceived_rms, 1e-12))
        target_db = 20 * np.log10(target / np.sqrt(2.0))
        error = perceived_db - target_db

        # Classify the error
        bare_limit = 20 * np.log10(target / max(Sf, 1e-12))
        note = ""
        if freq < third:
            note = "sub-fc/3 (expected)"
        elif abs(20 * np.log10(G) - bare_limit) < 0.1:
            note = "clamped (expected)"
        elif abs(error) < 0.5:
            note = "✅"

        print(f"  {freq:6.0f}  {20*np.log10(G):+8.2f}  {perceived_db:+10.2f}  "
              f"{target_db:+8.2f}  {error:+10.2f}  {note:>20s}")

        # Only count errors in the valid region (>= fc/3, not clamped)
        if freq >= third and abs(20 * np.log10(G) - bare_limit) > 0.5:
            errors.append(error)

    err_arr = np.array(errors) if errors else np.array([0.0])
    print(f"\n  ── Summary (valid region: ≥ fc/3, not clamped) ──")
    print(f"  Error: mean={err_arr.mean():+.2f} dB  "
          f"std={err_arr.std():.2f} dB  "
          f"max={np.max(np.abs(err_arr)):.2f} dB")

    max_err = np.max(np.abs(err_arr))
    if max_err < 0.5:
        print(f"  ✅ Full stack verified — model, EQ preprocessor, and DSP agree")
    else:
        print(f"  ⚠️  Max error {max_err:.2f} dB exceeds 0.5 dB threshold")
    print(f"  (sub-fc/3 and clamped frequencies excluded — those are design limits)")

    return {"errors": err_arr, "eq_curve": eq_curve}


def run():
    """Run the model verification."""
    verify_model()
    verify_full_stack()


if __name__ == "__main__":
    run()
