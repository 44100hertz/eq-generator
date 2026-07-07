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
from eqgen.pipeline import run_pipeline


SPEAKERS = {
    "flat": lambda f: 1.0,
    "small": lambda f: 0.25 if f <= 50 else (1.0 if f >= 100 else 0.25 + 0.75 * (f - 50) / 50),
    "steep": lambda f: 0.063 if f <= 40 else (1.0 if f >= 80
        else 0.063 + 0.937 * (f - 40) / 40),
}


def generate_header(speaker_name: str, speaker_fn, fs: float, fc: float,
                    h2: float, h3: float, max_bands: int) -> str:
    """Generate C header content for a speaker model's EQ coefficients.

    Uses the pipeline (measurement → correction) approach. For synthetic
    speaker models, we compute correction = 1/speaker_response."""
    f_min = fc / 2.0
    f_max = min(16000.0, fs * 0.49)
    freqs = np.logspace(np.log10(f_min), np.log10(f_max), 200)

    # Target: flat perceived output → correction = 1/speaker
    target_linear = np.array([1.0 / max(speaker_fn(f), 1e-12) for f in freqs])
    target_db = np.clip(20.0 * np.log10(np.maximum(target_linear, 1e-12)), -24.0, 24.0)

    # Pre-gain: shift the target down so the EQ only needs cuts.
    max_gain_db = max(0.0, float(np.max(target_db)))
    pre_gain = 10.0 ** (max_gain_db / 20.0) if max_gain_db > 0.0 else 1.0
    pre_gain_q16 = int(round(pre_gain * 65536.0))
    shifted_db = target_db - max_gain_db if max_gain_db > 0.0 else target_db

    fit = fit_eq_curve(freqs, shifted_db, fs, max_bands=max_bands,
                       min_freq=f_min, max_freq=f_max)

    bq_q28 = quantize_biquads_q28(fit.biquads)

    lines = []
    lines.append(f"// Auto-generated EQ coefficients for: {speaker_name}")
    lines.append(f"// fs={fs:.0f} Hz, fc={fc:.0f} Hz, h2_amp={h2:.2f}, h3_amp={h3:.2f}")
    lines.append(f"// {fit.n_bands} biquads, Q4.28 format")
    lines.append(f"// Regenerate: python -m eqgen.cli.export --speaker {speaker_name} -o <path>")
    lines.append("")
    lines.append("#pragma once")
    lines.append("#include <stdint.h>")
    lines.append("")
    lines.append(f"#define EQGEN_CUTOFF_HZ            {fc:.1f}f")
    lines.append(f"#define EQGEN_H2_AMP               {h2:.3f}f")
    lines.append(f"#define EQGEN_H3_AMP               {h3:.3f}f")
    lines.append(f"#define EQGEN_PRE_GAIN_Q16          {pre_gain_q16}  // {pre_gain:.2f}x = {max_gain_db:+.1f} dB")
    lines.append(f"#define EQGEN_FS                   {fs:.0f}")
    lines.append(f"#define EQGEN_N_BIQUADS            {fit.n_bands}")
    lines.append("")
    lines.append(f"static const int32_t eqgen_coeffs_q28[{fit.n_bands * 5}] = {{")
    for i, bq in enumerate(bq_q28):
        band = fit.bands[i]
        lines.append(
            f"    {bq.b0:>11d}, {bq.b1:>11d}, {bq.b2:>11d}, "
            f"{bq.a1:>11d}, {bq.a2:>11d},"
            f"  // [{i}] f0={band['f0']:6.1f} Hz gain={band['gain_db']:+5.1f} dB"
        )
    lines.append("};")
    lines.append("")
    lines.append(f"// Target sample rate for these coefficients")
    lines.append(f"#define EQGEN_FS {fs:.0f}")
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
    parser.add_argument("--max-bands", type=int, default=6,
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
