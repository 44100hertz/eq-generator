#!/usr/bin/env python3
"""
Bass evenness test: measures perceived output flatness across bass
frequencies using the full pipeline (EQ preprocess → enhancer → speaker).

Uses the Technics sanity check measurement data to model the real speaker.
Reveals how h2/h3 values affect harmonic-to-fundamental ratio and
perceived level evenness.
"""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from model import (
    butterworth_lp_mag, butterworth_hp_mag,
    model_perceived_amplitude, preprocess_eq_curve,
)
from dsp import BassEnhancerConfig, process_fixed_v2
from analysis import generate_test_sine, goertzel_magnitude

# ── Load real speaker measurement from sanity check data ──────────────

def load_speaker_model():
    """Build an interpolated speaker response from the sanity check
    measurement + target WAVs."""
    from scipy.io import wavfile

    ROOT = Path(__file__).parent.parent.parent
    meas_path = ROOT / "measurements/technics/standing/measurement2.wav"
    target_path = ROOT / "measurements/technics/standing/target.wav"

    rate_m, meas = wavfile.read(str(meas_path))
    rate_t, target = wavfile.read(str(target_path))
    if meas.ndim == 2: meas = meas.mean(axis=1)
    if target.ndim == 2: target = target.mean(axis=1)
    meas = meas.astype(float) / 32768.0
    target = target.astype(float) / 32768.0

    n = min(len(meas), len(target))
    window = np.hanning(n)
    meas_fft = np.abs(np.fft.rfft(meas[:n] * window))
    target_fft = np.abs(np.fft.rfft(target[:n] * window))
    freqs = np.fft.rfftfreq(n, 1.0 / rate_m)

    # Speaker response = measurement / target (normalized)
    # We want S(f) such that speaker_output = input * S(f)
    def S(f):
        if f <= 0 or f >= rate_m / 2:
            return 0.0
        idx = np.searchsorted(freqs, f)
        idx = min(idx, len(freqs) - 1)
        # Smooth: average nearby bins
        lo = max(0, idx - 3)
        hi = min(len(freqs) - 1, idx + 3)
        m = np.mean(meas_fft[lo:hi+1])
        t = np.mean(target_fft[lo:hi+1])
        return m / max(t, 1e-12)

    return S, freqs, meas_fft, target_fft


def run_evenness_test(fc=60.0, h2=5.0, h3=5.0, ramp_db=-60.0):
    """Test perceived output flatness across bass frequencies."""
    S, freqs_fft, meas_fft, target_fft = load_speaker_model()
    fs = 44100.0
    cfg = BassEnhancerConfig(
        cutoff=fc, h2_amp=h2, h3_amp=h3,
        fs=fs, env_release=200,
        limiter_release_ms=200, limiter_threshold=1.0,
    )

    # Build a flat target curve and compute EQ
    bass_freqs = np.arange(25, 180, 5)
    target_curve = {f: 1.0 for f in bass_freqs}
    eq_curve = preprocess_eq_curve(target_curve, S, fc, h2, h3, ramp_db=ramp_db)

    print(f"fc={fc:.0f} Hz  h2={h2:.2f}  h3={h3:.2f}  ramp={ramp_db:.0f} dB")
    print(f"{'Freq':>6s}  {'S(f)dB':>8s}  {'EQ dB':>8s}  {'H1 ear':>8s}  "
          f"{'H2 ear':>8s}  {'H3 ear':>8s}  {'H2/H1':>8s}  {'H3/H1':>8s}  "
          f"{'Total':>8s}  {'Flat?':>6s}")
    print("-" * 88)

    results = []
    for freq in bass_freqs:
        G = eq_curve.get(freq, 0.01)

        # Run through DSP
        duration = 0.1
        pad = int(0.2 * fs)
        sig = generate_test_sine([freq], [G], duration, fs, stereo=False)
        stereo = np.zeros((2, len(sig) + 2 * pad))
        stereo[0, pad:pad + len(sig)] = sig
        stereo[1, :] = stereo[0, :]
        out = process_fixed_v2(stereo, cfg)
        ss = slice(pad + len(sig) // 2, pad + len(sig))
        mono = (out[0, ss] + out[1, ss]) / 2.0

        # Measure harmonics at ear (after speaker)
        h1 = goertzel_magnitude(mono, freq, fs)
        h2_out = goertzel_magnitude(mono, 2 * freq, fs)
        h3_out = goertzel_magnitude(mono, 3 * freq, fs)

        Sf  = S(freq)
        S2f = S(2 * freq)
        S3f = S(3 * freq)

        h1_ear = h1 * Sf
        h2_ear = h2_out * S2f
        h3_ear = h3_out * S3f

        # Total perceived RMS (fundamental + harmonics)
        total_rms = np.sqrt(h1_ear**2 + h2_ear**2 + h3_ear**2)

        Sf_db = 20 * np.log10(max(Sf, 1e-12))
        G_db = 20 * np.log10(max(G, 1e-12))
        h1_db = 20 * np.log10(max(h1_ear, 1e-12))
        h2_db = 20 * np.log10(max(h2_ear, 1e-12))
        h3_db = 20 * np.log10(max(h3_ear, 1e-12))
        total_db = 20 * np.log10(max(total_rms, 1e-12))

        h2_rel = 20 * np.log10(h2_ear / max(h1_ear, 1e-12))
        h3_rel = 20 * np.log10(h3_ear / max(h1_ear, 1e-12))

        # Flatness: total_rms relative to 1.0 (target)
        flatness_db = total_db - 20 * np.log10(1.0)

        flag = ""
        if abs(flatness_db) > 3:
            flag = "⚠️"
        elif abs(flatness_db) <= 1:
            flag = "✓"

        print(f"{freq:6.0f}  {Sf_db:+8.1f}  {G_db:+8.2f}  {h1_db:+8.2f}  "
              f"{h2_db:+8.2f}  {h3_db:+8.2f}  {h2_rel:+8.2f}  {h3_rel:+8.2f}  "
              f"{total_db:+8.2f}  {flag:>6s}")

        results.append({
            "freq": freq,
            "Sf_db": Sf_db, "G_db": G_db,
            "h1_ear_db": h1_db, "h2_ear_db": h2_db, "h3_ear_db": h3_db,
            "h2_rel": h2_rel, "h3_rel": h3_rel,
            "total_db": total_db, "flatness_db": flatness_db,
        })

    # Summary
    arr = np.array([r["total_db"] for r in results])
    valid = arr[arr > -80]
    if len(valid) > 0:
        print(f"\n  Perceived output (valid region):")
        print(f"    mean={valid.mean():+.1f} dB  std={valid.std():.1f} dB  "
              f"range={valid.max()-valid.min():.1f} dB")
        if valid.std() > 3:
            print(f"    ⚠️  UNEVEN — >3 dB standard deviation across bass")
        else:
            print(f"    ✓  Even — <3 dB variation across bass")

    # Harmonic dominance check
    h2_rel_arr = np.array([r["h2_rel"] for r in results])
    h3_rel_arr = np.array([r["h3_rel"] for r in results])
    print(f"\n  Harmonic-to-fundamental ratios:")
    print(f"    H2/H1: mean={h2_rel_arr.mean():+.1f} dB  max={h2_rel_arr.max():+.1f} dB")
    print(f"    H3/H1: mean={h3_rel_arr.mean():+.1f} dB  max={h3_rel_arr.max():+.1f} dB")
    if h2_rel_arr.max() > 0:
        print(f"    ⚠️  H2 exceeds fundamental — harmonics dominate!")

    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fc", type=float, default=60)
    p.add_argument("--h2", type=float, default=5.0)
    p.add_argument("--h3", type=float, default=5.0)
    args = p.parse_args()
    run_evenness_test(fc=args.fc, h2=args.h2, h3=args.h3)
