"""
Sweep-based analysis through the C enhancer pipeline — full speaker path.

Runs actual audio through the real C DSP stack including volume LUT (overboost),
smart volume pre-gain, and loudness shelf — the exact pipeline that runs on
the ESP32 firmware and the desktop PipeWire filter.

This is the single source of truth for end-to-end pipeline analysis.
Replace ad-hoc sine→pack→enhancer→unpack code in tests with calls here.
"""

import struct
from typing import Dict, List, Optional, Tuple

import numpy as np

from eqgen import enhancer_ffi as effi
from eqgen.analysis import goertzel_magnitude


# ── Pipeline runner: configure DSP → run tone(s) → measure ────────────

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
    vol_gain: float = 1.0,
    pre_gain: Optional[float] = None,
    loudness_boost: float = 0.0,
) -> Dict[float, Dict[str, float]]:
    """Run sine tones through the C enhancer and measure output harmonics.

    The FULL pipeline:
      1. Multiply by vol_gain (BT volume LUT — includes overboost on speaker)
      2. enhancer_process_stereo: DC blocker → pre_gain → EQ → LP/HP split →
         crossfade → loudness shelf → tanh clamp

    For each frequency in freqs_hz, measures fundamental, H2, H3 amplitudes.

    Returns:
        dict mapping freq → {"fundamental", "h2", "h3", "rms"}
    """
    n_samples = int(duration_sec * fs)
    steady_start = n_samples // 4

    enh = effi.create_enhancer(
        cutoff_hz=fc, h2_amp=h2, h3_amp=h3,
        release_secs=release_secs, fs=fs,
        pre_gain=pre_gain if pre_gain is not None else 1.0,
        coeffs=eq_coeffs,
    )

    results = {}
    for f_test in freqs_hz:
        t = np.arange(n_samples) / fs
        sine = amplitude * np.sin(2.0 * np.pi * f_test * t)

        # int16 stereo PCM (matches real audio path)
        pcm = bytearray(n_samples * 4)
        for i in range(n_samples):
            v = int(np.clip(sine[i] * 32767, -32768, 32767))
            struct.pack_into('<hh', pcm, i * 4, v, v)

        effi.reset_enhancer(enh)
        for i in range(0, len(pcm), 4):
            l = struct.unpack_from('<h', pcm, i)[0]
            r = struct.unpack_from('<h', pcm, i + 2)[0]

            # Full signal chain: vol_gain → enhancer
            lf = (l / 32768.0) * vol_gain
            rf = (r / 32768.0) * vol_gain
            lf, rf = effi.process_stereo_frame(enh, lf, rf)

            struct.pack_into('<hh', pcm, i, int(np.clip(lf * 32767, -32768, 32767)), int(np.clip(rf * 32767, -32768, 32767)))

        out_float = np.array([
            struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
            for i in range(steady_start, n_samples)
        ])

        rms_out = float(np.sqrt(np.mean(out_float ** 2)))
        fund_amp = float(goertzel_magnitude(out_float, f_test, fs))
        h2_amp = float(goertzel_magnitude(out_float, 2 * f_test, fs)) if 2 * f_test < fs / 2 else 0
        h3_amp = float(goertzel_magnitude(out_float, 3 * f_test, fs)) if 3 * f_test < fs / 2 else 0

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
    """Generate a human-readable report from sweep results."""
    lines = []
    amp_db = 20.0 * np.log10(amplitude_in)

    lines.append(f"  Input amplitude: {amplitude_in:.4f} ({amp_db:+.0f} dBFS)")
    lines.append("")
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

    # Flatness
    rmss = np.array([results[f]["rms"] for f in freqs])
    if len(rmss) > 0:
        ref = float(np.mean(rmss))
        spread = 20.0 * np.log10(max(np.max(rmss), 1e-6) / max(np.min(rmss), 1e-6))
        lines.append("")
        lines.append(f"  ── Flatness ──")
        lines.append(f"  Mean RMS: {ref:.4f} ({20*np.log10(ref):+.1f} dBFS)  spread: {spread:.1f} dB")

    return "\n".join(lines)


def measure_harmonics_vs_amplitude(
    freq: float,
    fc: float = 60.0,
    h2: float = 0.5,
    h3: float = 1.0,
    fs: float = 44100.0,
    amplitudes: Optional[List[float]] = None,
    eq_coeffs: Optional[List[float]] = None,
    duration_sec: float = 0.5,
    vol_gain: float = 1.0,
    pre_gain: Optional[float] = None,
) -> List[Dict[str, float]]:
    """Measure harmonic levels vs input amplitude through full pipeline.

    Returns list of dicts: {"amplitude", "fundamental", "h2", "h3", "rms"}.
    """
    if amplitudes is None:
        amplitudes = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09]
    if eq_coeffs is None:
        eq_coeffs = []

    n_samples = int(duration_sec * fs)
    steady_start = n_samples // 4

    enh = effi.create_enhancer(
        cutoff_hz=fc, h2_amp=h2, h3_amp=h3,
        release_secs=0.2, fs=fs,
        pre_gain=pre_gain if pre_gain is not None else 1.0,
        coeffs=eq_coeffs,
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
            lf = (l / 32768.0) * vol_gain
            rf = (r / 32768.0) * vol_gain
            lf, rf = effi.process_stereo_frame(enh, lf, rf)
            struct.pack_into('<hh', pcm, i,
                int(np.clip(lf * 32767, -32768, 32767)),
                int(np.clip(rf * 32767, -32768, 32767)))

        out_float = np.array([
            struct.unpack_from('<h', pcm, i * 4)[0] / 32768.0
            for i in range(steady_start, n_samples)
        ])

        rms_out = float(np.sqrt(np.mean(out_float ** 2)))
        fund_amp = float(goertzel_magnitude(out_float, freq, fs))
        h2_amp = float(goertzel_magnitude(out_float, 2 * freq, fs)) if 2 * freq < fs / 2 else 0
        h3_amp = float(goertzel_magnitude(out_float, 3 * freq, fs)) if 3 * freq < fs / 2 else 0

        results.append({
            "amplitude": amp,
            "fundamental": fund_amp,
            "h2": h2_amp,
            "h3": h3_amp,
            "rms": rms_out,
        })

    effi.destroy_enhancer(enh)
    return results


# ── Volume LUT: exact mirror of smart_volume.h / firmware ─────────────

def build_vol_lut(
    vol: int,
    speaker_level_db: int = 60,
    overboost_db: float = 0.0,
    compensation_db: float = 0.0,
) -> float:
    """Compute vol_lut[vol] exactly as the C firmware does."""
    if vol <= 0:
        return 0.0

    sv_db_floor = float(-speaker_level_db)
    if sv_db_floor > -24.0:
        sv_db_floor = -24.0
    if sv_db_floor < -80.0:
        sv_db_floor = -80.0

    db = (float(vol) / 127.0) * (-sv_db_floor) + sv_db_floor + compensation_db + overboost_db
    if db > overboost_db:
        db = overboost_db
    return 10.0 ** (db / 20.0)


def build_full_vol_lut(
    speaker_level_db: int = 60,
    overboost_db: float = 0.0,
    compensation_db: float = 0.0,
) -> np.ndarray:
    """Build the full 128-entry vol_lut matching firmware."""
    lut = np.zeros(128, dtype=np.float64)
    for v in range(1, 128):
        lut[v] = build_vol_lut(v, speaker_level_db, overboost_db, compensation_db)
    return lut


def compute_smart_volume(
    vol: int,
    pg_loud: float = 1.0,
    quiet_shelf_db: float = 8.0,
) -> dict:
    """Compute pre_gain/loudness_boost for a volume level (matches smart_volume.h)."""
    t = float(vol) / 127.0

    atten_norm = 1.0 - t
    shelf_db = quiet_shelf_db * (atten_norm ** 0.33)
    shelf_linear = 10.0 ** (shelf_db / 20.0)
    boost = shelf_linear - 1.0

    max_shelf_linear = 10.0 ** (quiet_shelf_db / 20.0)
    pg_quiet = pg_loud / max_shelf_linear
    pre_gain = pg_quiet + t * (pg_loud - pg_quiet)

    return {
        "pre_gain": pre_gain,
        "shelf_db": shelf_db,
        "boost": boost,
    }
