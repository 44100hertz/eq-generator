"""
Bass evenness test: measures perceived output flatness using the C enhancer.

Runs sine sweeps through the actual C DSP (via enhancer_ffi) and
reports fundamental + harmonic levels, verifying that the enhancer
produces even output across bass frequencies.

Replaces the old model + Python-DSP approach with direct C pipeline
measurement.
"""

import sys
import struct
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen import enhancer_ffi as effi
from eqgen.pipeline import build_default_eq_coeffs
from eqgen.analysis import generate_test_sine, goertzel_magnitude


def run():
    """Entry point for the test runner."""
    run_evenness_test()


def run_evenness_test(fc=60.0, h2=0.5, h3=1.0, fs=44100.0):
    """Test output flatness across bass frequencies through the C enhancer."""
    print("=" * 70)
    print(f"  BASS EVENNESS: C enhancer sweep, fc={fc}, h2={h2}, h3={h3}")
    print("=" * 70)

    eq_coeffs = build_default_eq_coeffs(fs)
    enh = effi.create_enhancer(
        cutoff_hz=fc, h2_amp=h2, h3_amp=h3,
        release_secs=0.2, fs=fs, coeffs=eq_coeffs,
    )

    bass_freqs = list(range(25, 185, 5))
    amp_in = 0.001  # very low to avoid clipping
    duration = 1.0
    pad = int(0.2 * fs)

    print(f"\n  Input: {amp_in:.4f} ({20*np.log10(amp_in):+.0f} dBFS)")
    print(f"\n  {'Freq':>6s}  {'Fund dB':>8s}  {'H2 rel':>8s}  "
          f"{'H3 rel':>8s}  {'RMS dB':>8s}")
    print(f"  {'-'*42}")

    results = []
    for freq in bass_freqs:
        signal = generate_test_sine([freq], [amp_in], duration, fs, stereo=False)

        total_n = len(signal) + 2 * pad
        pcm = bytearray(total_n * 4)
        for i in range(total_n):
            v = 0.0
            if pad <= i < pad + len(signal):
                v = signal[i - pad]
            sv = int(np.clip(v * 32767, -32768, 32767))
            struct.pack_into('<hh', pcm, i * 4, sv, sv)

        effi.reset_enhancer(enh)
        for i in range(0, len(pcm), 4):
            l = struct.unpack_from('<h', pcm, i)[0]
            r = struct.unpack_from('<h', pcm, i + 2)[0]
            l_out, r_out = effi.process_stereo_frame(enh, l, r)
            struct.pack_into('<hh', pcm, i, int(l_out), int(r_out))

        start = pad + len(signal) // 4
        end = pad + len(signal)
        mono = np.array([
            struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
            for i in range(start, end)
        ])

        rms_out = float(np.sqrt(np.mean(mono ** 2)))
        fund = float(goertzel_magnitude(mono, freq, fs))
        h2_amp = float(goertzel_magnitude(mono, 2 * freq, fs)) if 2 * freq < fs / 2 else 0.0
        h3_amp = float(goertzel_magnitude(mono, 3 * freq, fs)) if 3 * freq < fs / 2 else 0.0

        fund_db = 20.0 * np.log10(max(fund, 1e-12))
        h2_rel = 20.0 * np.log10(max(h2_amp, 1e-12) / max(fund, 1e-12))
        h3_rel = 20.0 * np.log10(max(h3_amp, 1e-12) / max(fund, 1e-12))
        rms_db = 20.0 * np.log10(max(rms_out, 1e-12))

        flag = "✓" if abs(rms_db + 60) < 3 else "⚠️"
        print(f"  {freq:6.0f}  {fund_db:+8.2f}  {h2_rel:+8.2f}  "
              f"{h3_rel:+8.2f}  {rms_db:+8.2f}  {flag}")

        results.append({
            "freq": freq, "fund_db": fund_db,
            "h2_rel": h2_rel, "h3_rel": h3_rel,
            "rms_db": rms_db,
        })

    effi.destroy_enhancer(enh)

    # Summary
    rms_vals = np.array([r["rms_db"] for r in results])
    valid = rms_vals[rms_vals > -80]
    if len(valid) > 0:
        print(f"\n  ── Summary ──")
        print(f"  RMS output: mean={valid.mean():+.1f} dBFS  "
              f"std={valid.std():.1f} dB  range={valid.max()-valid.min():.1f} dB")
        if valid.std() > 3:
            print(f"  ⚠️  UNEVEN — >3 dB standard deviation across bass")
        else:
            print(f"  ✅ Even — <3 dB variation across bass")

        # Harmonic dominance
        h2r = np.array([r["h2_rel"] for r in results])
        h3r = np.array([r["h3_rel"] for r in results])
        print(f"  H2/H1: mean={h2r.mean():+.1f} dB  max={h2r.max():+.1f} dB")
        print(f"  H3/H1: mean={h3r.mean():+.1f} dB  max={h3r.max():+.1f} dB")

    return results


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fc", type=float, default=60)
    p.add_argument("--h2", type=float, default=0.5)
    p.add_argument("--h3", type=float, default=1.0)
    args = p.parse_args()
    run_evenness_test(fc=args.fc, h2=args.h2, h3=args.h3)
