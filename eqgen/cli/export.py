#!/usr/bin/env python3
"""
Export fitted IIR EQ biquads as a C header for ESP32 firmware.

Generates src/eq_coeffs.h with Q4.28 fixed-point coefficients.

Usage:
    python -m eqgen.cli.export --speaker small --h2 0.33 --h3 0.33 --fc 60 \
        --max-bands 6 -o src/eq_coeffs.h
    python -m eqgen.cli.export --speaker flat --fc 60 -o src/eq_coeffs.h
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.eq_fit import fit_eq_curve
from eqgen.quantize import BiquadQ28, quantize_biquads_q28
from eqgen.pipeline import run_pipeline, pre_gain_from_max_gain
from eqgen.presets import MAX_IIR_BANDS


SPEAKERS = {
    "flat": lambda f: 1.0,
    "small": lambda f: 0.25 if f <= 50 else (1.0 if f >= 100 else 0.25 + 0.75 * (f - 50) / 50),
    "steep": lambda f: 0.063 if f <= 40 else (1.0 if f >= 80
        else 0.063 + 0.937 * (f - 40) / 40),
}


def _fit_header_bands(freqs, shifted_db, fs, max_bands, f_min, f_max):
    """Fit biquads at a single sample rate, return (bq_q28, bands, fit)."""
    fit = fit_eq_curve(freqs, shifted_db, fs, max_bands=max_bands,
                       min_freq=f_min, max_freq=f_max)
    bq_q28 = quantize_biquads_q28(fit.biquads)
    return bq_q28, fit.bands, fit.n_bands


def generate_header(speaker_name: str, speaker_fn, fs: float, fc: float,
                    h2: float, h3: float, max_bands: int) -> str:
    """Generate C header content for a speaker model's EQ coefficients
    at both 44100 Hz and 48000 Hz."""
    f_min = fc / 2.0
    f_max = min(16000.0, fs * 0.49)
    freqs = np.logspace(np.log10(f_min), np.log10(f_max), 200)

    # Target: flat corrected output → correction = 1/speaker
    target_linear = np.array([1.0 / max(speaker_fn(f), 1e-12) for f in freqs])
    target_db = np.clip(20.0 * np.log10(np.maximum(target_linear, 1e-12)), -24.0, 24.0)

    max_gain_db = max(0.0, float(np.max(target_db)))
    pre_gain = pre_gain_from_max_gain(max_gain_db)
    pre_gain_q16 = int(round(pre_gain * 65536.0))

    # Fit at both 44100 and 48000
    bq_q28_44, bands_44, n_bands_44 = _fit_header_bands(
        freqs, target_db, 44100.0, max_bands, f_min, f_max)
    bq_q28_48, bands_48, n_bands_48 = _fit_header_bands(
        freqs, target_db, 48000.0, max_bands, f_min, f_max)

    lines = []
    lines.append(f"// Auto-generated EQ coefficients for: {speaker_name}")
    lines.append(f"// Generated: python -m eqgen.cli.export --speaker {speaker_name} --max-bands {max_bands}")
    lines.append(f"// Two coefficient arrays — one per sample rate — because biquad coefficients")
    lines.append(f"// are sample-rate-dependent.  At runtime, eqgen_get_coeffs(rate) selects")
    lines.append(f"// the correct array for 44100 Hz or 48000 Hz.")
    lines.append(f"//")
    lines.append(f"// 44.1k: {n_bands_44} biquads, 48k: {n_bands_48} biquads — Q4.28 fixed-point")
    lines.append(f"// fc={fc:.0f} Hz, h2_amp={h2:.2f}, h3_amp={h3:.2f}")
    lines.append("")
    lines.append("#pragma once")
    lines.append("#include <stdint.h>")
    lines.append("")
    lines.append(f"#define EQGEN_CUTOFF_HZ            {fc:.1f}f")
    lines.append(f"#define EQGEN_H2_AMP               {h2:.3f}f")
    lines.append(f"#define EQGEN_H3_AMP               {h3:.3f}f")
    lines.append(f"#define EQGEN_PRE_GAIN_Q16          {pre_gain_q16}  // {pre_gain:.2f}x = {20*np.log10(pre_gain):+.1f} dB")
    lines.append("")
    lines.append("/* ── Smart volume (AVRCP-based loudness compensation) ────────────── */")
    lines.append("#define EQGEN_LOUDNESS_FC_HZ         200.0f  /* one-pole shelf corner freq  */")
    lines.append("#define EQGEN_QUIET_SHELF_DB           8.0f  /* boost at DC when vol → 0   */")
    lines.append("#define EQGEN_QUIET_FUNDAMENTAL_BLEED  0.25f /* 25% LP bleed when h2 halved */")
    lines.append("")
    lines.append("/* Derived: quiet harmonic amplitudes (half of configured) */")
    lines.append(f"#define EQGEN_QUIET_H2_AMP       ({h2:.3f}f * 0.5f)")
    lines.append(f"#define EQGEN_QUIET_H3_AMP       ({h3:.3f}f * 0.5f)")
    lines.append(f"#define EQGEN_FS_44100              44100")
    lines.append(f"#define EQGEN_FS_48000              48000")
    lines.append(f"#define EQGEN_N_BIQUADS            {n_bands_44}")
    lines.append("")

    # 44100 Hz coefficients
    lines.append(f"static const int32_t eqgen_coeffs_q28_44100[{n_bands_44 * 5}] = {{")
    for i, bq in enumerate(bq_q28_44):
        band = bands_44[i]
        lines.append(
            f"    {bq.b0:>11d}, {bq.b1:>11d}, {bq.b2:>11d}, "
            f"{bq.a1:>11d}, {bq.a2:>11d},"
            f"  // [{i}] f0={band['f0']:6.1f} Hz gain={band['gain_db']:+5.1f} dB"
        )
    lines.append("};")
    lines.append("")

    # 48000 Hz coefficients
    lines.append(f"static const int32_t eqgen_coeffs_q28_48000[{n_bands_48 * 5}] = {{")
    for i, bq in enumerate(bq_q28_48):
        band = bands_48[i]
        lines.append(
            f"    {bq.b0:>11d}, {bq.b1:>11d}, {bq.b2:>11d}, "
            f"{bq.a1:>11d}, {bq.a2:>11d},"
            f"  // [{i}] f0={band['f0']:6.1f} Hz gain={band['gain_db']:+5.1f} dB"
        )
    lines.append("};")
    lines.append("")

    lines.append("/** Select coefficient array for the given sample rate. */")
    lines.append("static inline const int32_t *eqgen_get_coeffs(int sample_rate) {")
    lines.append("    return (sample_rate == 48000) ? eqgen_coeffs_q28_48000 : eqgen_coeffs_q28_44100;")
    lines.append("}")
    lines.append("")
    lines.append("/** Select nominal Fs for the given sample rate. */")
    lines.append("static inline int eqgen_get_fs(int sample_rate) {")
    lines.append("    return (sample_rate == 48000) ? EQGEN_FS_48000 : EQGEN_FS_44100;")
    lines.append("}")
    lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Export EQ biquad coefficients as C header for ESP32")
    parser.add_argument("--speaker", choices=list(SPEAKERS.keys()),
                        default="small", help="Speaker model [default: small]")
    parser.add_argument("--fs", type=float, default=44100.0,
                        help="Sample rate in Hz")
    parser.add_argument("--fc", type=float, default=60.0,
                        help="Bass enhancement cutoff Hz")
    parser.add_argument("--h2", type=float, default=0.33,
                        help="2nd harmonic amplitude")
    parser.add_argument("--h3", type=float, default=0.33,
                        help="3rd harmonic amplitude")
    parser.add_argument("--max-bands", type=int, default=MAX_IIR_BANDS,
                        help="Max biquad bands")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output file (default: stdout)")

    args = parser.parse_args()

    speaker_fn = SPEAKERS[args.speaker]
    header = generate_header(args.speaker, speaker_fn, args.fs, args.fc,
                             args.h2, args.h3, args.max_bands)

    if args.output:
        p = Path(args.output)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(header)
        print(f"Wrote {args.output}")
    else:
        print(header)


if __name__ == "__main__":
    main()
