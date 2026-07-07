"""
Harmonic purity and linearity tests using the C enhancer pipeline.

Runs sine tones through the actual C DSP (via enhancer_ffi) and
measures fundamental, 2nd, and 3rd harmonic levels via Goertzel analysis.

Covers:
  - Full pipeline sweep (C enhancer output across bass frequencies)
  - Harmonic linearity vs input amplitude
  - 2nd/3rd harmonic measurement across the bass range
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.dsp import BassEnhancerConfig, butterworth_lp_mag, butterworth_hp_mag
from eqgen.analysis import generate_test_sine, goertzel_magnitude, measure_tones, format_tone_report
from eqgen import sweep
from eqgen.pipeline import build_default_eq_coeffs


# ═══════════════════════════════════════════════════════════════════════════════
# Full pipeline sweep: C enhancer output across bass frequencies
# ═══════════════════════════════════════════════════════════════════════════════

def pipeline_sweep(cutoff=60.0, h2=0.5, h3=1.0, fs=44100.0):
    """Run sine sweeps through the C enhancer and report harmonic output."""
    print("=" * 70)
    print(f"  C ENHANCER SWEEP: cutoff={cutoff} Hz, h2={h2}, h3={h3}")
    print("=" * 70)

    eq_coeffs = build_default_eq_coeffs(fs)

    test_freqs = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 100, 120,
                   150, 200, 300, 500, 1000, 2000]

    results = sweep.run_sine_sweep(
        freqs_hz=test_freqs,
        eq_coeffs_q28=eq_coeffs,
        fc=cutoff, h2=h2, h3=h3, fs=fs,
        amplitude=0.001,
        duration_sec=1.0,
    )

    report = sweep.sweep_report(results, amplitude_in=0.001)
    print(report)

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Harmonic linearity: vary input amplitude at a single frequency
# ═══════════════════════════════════════════════════════════════════════════════

def harmonic_linearity(freq=50.0, cutoff=60.0, h2=0.5, h3=1.0, fs=44100.0):
    """Measure harmonic levels vs input amplitude to verify linearity."""
    print("\n" + "=" * 70)
    print(f"  HARMONIC LINEARITY: {freq} Hz, cutoff={cutoff}, h2={h2}, h3={h3}")
    print("=" * 70)

    eq_coeffs = build_default_eq_coeffs(fs)
    amplitudes = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09]

    results = sweep.measure_harmonics_vs_amplitude(
        freq=freq, fc=cutoff, h2=h2, h3=h3, fs=fs,
        amplitudes=amplitudes,
        eq_coeffs_q28=eq_coeffs,
    )

    print(f"\n  {'In dBFS':>8s}  {'Fund dBFS':>10s}  "
          f"{'H2 rel':>8s}  {'H3 rel':>8s}  {'RMS dBFS':>10s}")
    print(f"  {'-'*50}")

    h2_rel_vals, h3_rel_vals = [], []
    for r in results:
        amp = r["amplitude"]
        fund = max(r["fundamental"], 1e-12)
        h2_amp = max(r["h2"], 1e-12)
        h3_amp = max(r["h3"], 1e-12)
        rms = max(r["rms"], 1e-12)

        in_db = 20.0 * np.log10(amp)
        fund_db = 20.0 * np.log10(fund)
        h2_rel = 20.0 * np.log10(h2_amp / fund)
        h3_rel = 20.0 * np.log10(h3_amp / fund)
        rms_db = 20.0 * np.log10(rms)

        h2_rel_vals.append(h2_rel)
        h3_rel_vals.append(h3_rel)

        print(f"  {in_db:+8.2f}  {fund_db:+10.2f}  "
              f"{h2_rel:+8.2f}  {h3_rel:+8.2f}  {rms_db:+10.2f}")

    if h2_rel_vals:
        h2_std = float(np.std(h2_rel_vals))
        h3_std = float(np.std(h3_rel_vals))
        linear = "LINEAR" if h2_std < 1.5 and h3_std < 1.5 else "NONLINEAR"
        print(f"\n  H2 σ={h2_std:.2f} dB  H3 σ={h3_std:.2f} dB  → {linear}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Detailed tone purity at a single frequency
# ═══════════════════════════════════════════════════════════════════════════════

def tone_purity(freq=50.0, amp=0.5, cutoff=60.0, h2=0.5, h3=1.0, fs=44100.0):
    """Detailed Goertzel measurement of a single tone through the C enhancer."""
    import struct
    from eqgen import enhancer_ffi as effi

    print("\n" + "=" * 70)
    print(f"  TONE PURITY: {freq} Hz, A={amp} ({20*np.log10(amp):+.1f} dBFS)")
    print(f"  cutoff={cutoff}, h2={h2}, h3={h3}")
    print("=" * 70)

    eq_coeffs = build_default_eq_coeffs(fs)
    duration = 0.5
    pad = int(0.3 * fs)
    n_samples = int(duration * fs)

    signal = generate_test_sine([freq], [amp], duration, fs, stereo=False)

    # Pack as int16 stereo
    total = n_samples + 2 * pad
    pcm = bytearray(total * 4)
    for i in range(total):
        v = 0.0
        if pad <= i < pad + n_samples:
            v = signal[i - pad]
        sv = int(np.clip(v * 32767, -32768, 32767))
        struct.pack_into('<hh', pcm, i * 4, sv, sv)

    enh = effi.create_enhancer(
        cutoff_hz=cutoff, h2_amp=h2, h3_amp=h3,
        release_secs=0.2, fs=fs, coeffs_q28=eq_coeffs,
    )

    for i in range(0, len(pcm), 4):
        l = struct.unpack_from('<h', pcm, i)[0]
        r = struct.unpack_from('<h', pcm, i + 2)[0]
        l_out, r_out = effi.process_stereo_frame(enh, l, r)
        struct.pack_into('<hh', pcm, i, l_out, r_out)

    effi.destroy_enhancer(enh)

    # Decode steady-state mono
    start = pad + n_samples // 2
    end = pad + n_samples
    mono = np.array([
        struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
        for i in range(start, end)
    ])

    m = measure_tones(mono, freq, fs, harmonics=(1, 2, 3))
    format_tone_report(m, "C Enhancer Output", fundamental=freq)

    return m


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all harmonic analysis tests through the C enhancer."""
    pipeline_sweep(cutoff=60.0, h2=0.5, h3=1.0)
    harmonic_linearity(freq=50.0, cutoff=60.0, h2=0.5, h3=1.0)
    tone_purity(freq=50.0, amp=0.5, cutoff=60.0, h2=0.5, h3=1.0)
    # Also test at 40 Hz
    tone_purity(freq=40.0, amp=0.5, cutoff=60.0, h2=0.5, h3=1.0)


if __name__ == "__main__":
    run()
