"""
Sweep-based analysis through the C enhancer pipeline.

Instead of a psychoacoustic model predicting enhancer behavior, this
module runs actual sine sweeps through the C DSP (via enhancer_ffi)
and measures the output — fundamental, 2nd harmonic, 3rd harmonic —
at the output of the real pipeline.

This replaces the model.py preprocess_eq_curve approach: the pipeline
IS the response, and sweeps measure what actually comes out.
"""

import struct
import numpy as np
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from eqgen import enhancer_ffi as effi
from eqgen.analysis import goertzel_magnitude


def run_sine_sweep(
    freqs_hz: List[float],
    eq_coeffs: List[float],
    fc: float = 60.0,
    h2: float = 0.5,
    h3: float = 1.0,
    fs: float = 44100.0,
    amplitude: float = 0.001,
    duration_sec: float = 1.0,
    release_secs: float = 0.2,
) -> Dict[float, Dict[str, float]]:
    """Run sine tones through the C enhancer and measure output harmonics.

    For each frequency in freqs_hz:
      1. Generate a pure sine at `amplitude`
      2. Pass through the C enhancer (with EQ + bass enhancement)
      3. Measure fundamental, H2, H3 amplitudes via Goertzel

    Returns:
        dict mapping freq → {
            "fundamental": linear amplitude at f,
            "h2": linear amplitude at 2f,
            "h3": linear amplitude at 3f,
            "rms": total RMS of output,
        }
    """
    n_samples = int(duration_sec * fs)
    steady_start = n_samples // 4  # skip startup transient

    enh = effi.create_enhancer(
        cutoff_hz=fc, h2_amp=h2, h3_amp=h3,
        release_secs=release_secs, fs=fs,
        coeffs=eq_coeffs,
    )

    results = {}
    for f_test in freqs_hz:
        t = np.arange(n_samples) / fs
        sine = amplitude * np.sin(2.0 * np.pi * f_test * t)

        # Convert to int16 stereo PCM
        pcm = bytearray(n_samples * 4)
        for i in range(n_samples):
            v = int(np.clip(sine[i] * 32767, -32768, 32767))
            struct.pack_into('<hh', pcm, i * 4, v, v)

        effi.reset_enhancer(enh)
        for i in range(0, len(pcm), 4):
            l = struct.unpack_from('<h', pcm, i)[0]
            r = struct.unpack_from('<h', pcm, i + 2)[0]
            l_out, r_out = effi.process_stereo_frame(enh, l, r)
            struct.pack_into('<hh', pcm, i, l_out, r_out)

        # Decode steady-state
        out_float = np.array([
            struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
            for i in range(steady_start, n_samples)
        ])

        rms_out = float(np.sqrt(np.mean(out_float ** 2)))
        fund_amp = float(goertzel_magnitude(out_float, f_test, fs))
        h2_amp = float(goertzel_magnitude(out_float, 2 * f_test, fs)) if 2 * f_test < fs / 2 else 0.0
        h3_amp = float(goertzel_magnitude(out_float, 3 * f_test, fs)) if 3 * f_test < fs / 2 else 0.0

        results[f_test] = {
            "fundamental": fund_amp,
            "h2": h2_amp,
            "h3": h3_amp,
            "rms": rms_out,
        }

    effi.destroy_enhancer(enh)
    return results


def sweep_report(
    results: Dict[float, Dict[str, float]],
    amplitude_in: float = 0.001,
) -> str:
    """Generate a human-readable report from sweep results.

    Returns a multi-line string showing per-frequency harmonic levels
    and flatness metrics.
    """
    lines = []
    amp_db = 20.0 * np.log10(amplitude_in)

    lines.append(f"  Input amplitude: {amplitude_in:.4f} ({amp_db:+.0f} dBFS)")
    lines.append(f"")
    header = (f"  {'Freq':>6s}  {'RMS':>8s}  {'dBFS':>7s}  "
              f"{'Fund':>9s}  {'H2 rel':>9s}  {'H3 rel':>9s}")
    lines.append(header)
    lines.append(f"  {'':->6s}  {'':->8s}  {'':->7s}  "
                 f"{'':->9s}  {'':->9s}  {'':->9s}")

    freqs = sorted(results.keys())
    for f in freqs:
        r = results[f]
        fund = r["fundamental"]
        h2 = r["h2"]
        h3 = r["h3"]
        rms = r["rms"]

        fund_db = 20.0 * np.log10(max(fund, 1e-12))
        h2_rel = 20.0 * np.log10(max(h2, 1e-12) / max(fund, 1e-12))
        h3_rel = 20.0 * np.log10(max(h3, 1e-12) / max(fund, 1e-12))
        rms_db = 20.0 * np.log10(max(rms, 1e-12))

        lines.append(f"  {f:6.0f}  {rms:8.4f}  {rms_db:+6.1f}  "
                      f"{fund_db:+8.1f}  {h2_rel:+8.1f}  {h3_rel:+8.1f}")

    # Flatness summary
    rmss = np.array([results[f]["rms"] for f in freqs])
    if len(rmss) > 0:
        ref = float(np.mean(rmss))
        spread = 20.0 * np.log10(max(np.max(rmss), 1e-6) / max(np.min(rmss), 1e-6))
        lines.append(f"")
        lines.append(f"  ── Flatness ──")
        lines.append(f"  Mean output RMS: {ref:.4f} ({20*np.log10(ref):+.1f} dBFS)")
        lines.append(f"  Min/Max spread:  {spread:.1f} dB")

        # Per-region
        bass = [results[f]["rms"] for f in freqs if f <= 120]
        treble = [results[f]["rms"] for f in freqs if f > 120]
        if bass:
            lines.append(f"  Bass (≤120 Hz):   {20*np.log10(max(bass)/max(min(bass),1e-6)):.1f} dB spread")
        if treble:
            lines.append(f"  Treble (>120 Hz): {20*np.log10(max(treble)/max(min(treble),1e-6)):.1f} dB spread")

    return "\n".join(lines)


def measure_harmonics_vs_amplitude(
    freq: float,
    fc: float = 60.0,
    h2: float = 0.5,
    h3: float = 1.0,
    fs: float = 44100.0,
    amplitudes: Optional[List[float]] = None,
    eq_coeffs_q28: Optional[List[int]] = None,
    duration_sec: float = 0.5,
) -> List[Dict[str, float]]:
    """Measure 2nd and 3rd harmonic levels vs input amplitude at a single frequency.

    Sweeps input amplitude and records the fundamental, H2, and H3 levels
    at the enhancer output.  Useful for verifying harmonic linearity.

    Returns list of dicts with keys: amplitude, fundamental, h2, h3, rms.
    """
    if amplitudes is None:
        amplitudes = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09]

    if eq_coeffs_q28 is None:
        eq_coeffs_q28 = []

    n_samples = int(duration_sec * fs)
    steady_start = n_samples // 4

    enh = effi.create_enhancer(
        cutoff_hz=fc, h2_amp=h2, h3_amp=h3,
        release_secs=0.2, fs=fs, coeffs_q28=eq_coeffs_q28,
    )

    results = []
    for amp in amplitudes:
        t = np.arange(n_samples) / fs
        sine = amp * np.sin(2.0 * np.pi * freq * t)

        pcm = bytearray(n_samples * 4)
        for i in range(n_samples):
            v = int(np.clip(sine[i] * 32767, -32768, 32767))
            struct.pack_into('<hh', pcm, i * 4, v, v)

        effi.reset_enhancer(enh)
        for i in range(0, len(pcm), 4):
            l = struct.unpack_from('<h', pcm, i)[0]
            r = struct.unpack_from('<h', pcm, i + 2)[0]
            l_out, r_out = effi.process_stereo_frame(enh, l, r)
            struct.pack_into('<hh', pcm, i, l_out, r_out)

        out_float = np.array([
            struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
            for i in range(steady_start, n_samples)
        ])

        rms_out = float(np.sqrt(np.mean(out_float ** 2)))
        fund_amp = float(goertzel_magnitude(out_float, freq, fs))
        h2_amp = float(goertzel_magnitude(out_float, 2 * freq, fs)) if 2 * freq < fs / 2 else 0.0
        h3_amp = float(goertzel_magnitude(out_float, 3 * freq, fs)) if 3 * freq < fs / 2 else 0.0

        results.append({
            "amplitude": amp,
            "fundamental": fund_amp,
            "h2": h2_amp,
            "h3": h3_amp,
            "rms": rms_out,
        })

    effi.destroy_enhancer(enh)
    return results
