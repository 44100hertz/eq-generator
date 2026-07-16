#!/usr/bin/env python3
"""
Export fitted IIR EQ biquads as a C header for ESP32 firmware.

Generates src/eq_coeffs.h with float coefficients.

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
from eqgen.model import perceptual_weight
from eqgen.dsp import pre_gain_from_max_gain
from eqgen.presets import MAX_IIR_BANDS
from eqgen.smart_volume import SV_LOUDNESS_FC, SV_LOUDNESS_Q, FM_SLOPE
from eqgen.header_gen import generate_c_header


SPEAKERS = {
    "flat": lambda f: 1.0,
    "small": lambda f: 0.25 if f <= 50 else (1.0 if f >= 100 else 0.25 + 0.75 * (f - 50) / 50),
    "steep": lambda f: 0.063 if f <= 40 else (1.0 if f >= 80
        else 0.063 + 0.937 * (f - 40) / 40),
}


def _fit_header_bands(freqs, shifted_db, fs, max_bands, f_min, f_max):
    """Fit biquads at a single sample rate, return (biquads, bands, n_bands)."""
    fit = fit_eq_curve(freqs, shifted_db, fs, max_bands=max_bands,
                       min_freq=f_min, max_freq=f_max,
                       error_weight=perceptual_weight)
    return fit.biquads, fit.bands, fit.n_bands


def _fit_export_biquads(speaker_fn, fc, max_bands):
    """Run the fitting pipeline for export, return (biquads, bands, pre_gain).

    Performs fitting at both 44100 Hz and 48000 Hz and returns the
    values needed by generate_header.
    """
    fs = 44100.0
    f_min = fc / 2.0
    f_max = min(16000.0, fs * 0.49)
    freqs = np.logspace(np.log10(f_min), np.log10(f_max), 200)

    # Target: flat corrected output → correction = 1/speaker
    target_linear = np.array([1.0 / max(speaker_fn(f), 1e-12) for f in freqs])
    target_db = np.clip(20.0 * np.log10(np.maximum(target_linear, 1e-12)), -24.0, 24.0)

    max_gain_db = max(0.0, float(np.max(target_db)))
    pre_gain = pre_gain_from_max_gain(max_gain_db)

    biquads_44, bands_44, _ = _fit_header_bands(
        freqs, target_db, 44100.0, max_bands, f_min, f_max)
    biquads_48, bands_48, _ = _fit_header_bands(
        freqs, target_db, 48000.0, max_bands, f_min, f_max)

    return biquads_44, bands_44, biquads_48, bands_48, pre_gain


def generate_header(speaker_name: str, biquads_44, bands_44,
                    biquads_48, bands_48, pre_gain: float,
                    fc: float, h2: float, h3: float, max_bands: int) -> str:
    """Generate C header content from pre-computed biquads."""
    n_bands_44 = len(bands_44)
    n_bands_48 = len(bands_48)

    header_lines = [
        f"Auto-generated EQ coefficients for: {speaker_name}",
        f"Generated: python -m eqgen.cli.export --speaker {speaker_name} --max-bands {max_bands}",
        "Two coefficient arrays — one per sample rate — because biquad coefficients",
        "are sample-rate-dependent.  At runtime, eqgen_get_coeffs(rate) selects",
        "the correct array for 44100 Hz or 48000 Hz.",
        "",
        f"44.1k: {n_bands_44} biquads, 48k: {n_bands_48} biquads — float coefficients",
        f"fc={fc:.0f} Hz, h2_amp={h2:.2f}, h3_amp={h3:.2f}",
    ]

    define_lines = [
        f"#define EQGEN_CUTOFF_HZ            {fc:.1f}f",
        f"#define EQGEN_H2_AMP               {h2:.3f}f",
        f"#define EQGEN_H3_AMP               {h3:.3f}f",
        f"#define EQGEN_PRE_GAIN             {pre_gain:.6f}f  // {pre_gain:.2f}x = {20*np.log10(pre_gain):+.1f} dB",
        "",
        "/* ── Smart volume (AVRCP-based loudness compensation) ────────────── */",
        "/* 2nd-order low shelf fitted to ISO 226:2023 equal-loudness contours. */",
        f"#define FM_SLOPE                    {FM_SLOPE:.4f}f  /* dB shelf per dB SPL drop */",
        f"#define EQGEN_LOUDNESS_FC_HZ        {SV_LOUDNESS_FC:.1f}f  /* shelf corner freq        */",
        f"#define EQGEN_LOUDNESS_Q            {SV_LOUDNESS_Q:.3f}f  /* shelf Q                  */",
        "#define EQGEN_QUIET_SHELF_DB           8.0f  /* boost at DC when vol → 0   */",
        "#define EQGEN_QUIET_FUNDAMENTAL_BLEED  0.40f /* max LP bleed when harmonics fully cut */",
        "",
        "/* Harmonic→bleed crossfade thresholds (t = vol/127) */",
        "#define EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T  0.20f  /* vol≤25: h2/h3=0, bleed=max */",
        "#define EQGEN_HARMONIC_BLEED_CROSSFADE_HI_T  0.50f  /* vol≥63: h2/h3=full, bleed=0 */",
        f"#define EQGEN_FS_44100              44100",
        f"#define EQGEN_FS_48000              48000",
        f"#define EQGEN_N_BIQUADS            {n_bands_44}",
        "",
        "/* ── Firmware config ────────────────────────────────────────────── */",
        "#define EQGEN_RELEASE_SECS           0.200f",
        "#define EQGEN_LIMITER_RELEASE_SECS   0.049f",
        '#define EQGEN_BT_DEVICE_NAME        "eqgen"',
        "#define EQGEN_SPEAKER_LEVEL_DB        60",
        "",
    ]

    return generate_c_header(
        biquads_44, bands_44,
        biquads_48, bands_48,
        header_lines=header_lines,
        define_lines=define_lines,
    )


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
    bq44, bn44, bq48, bn48, pg = _fit_export_biquads(
        speaker_fn, args.fc, args.max_bands)
    header = generate_header(args.speaker, bq44, bn44, bq48, bn48, pg,
                             args.fc, args.h2, args.h3, args.max_bands)

    if args.output:
        p = Path(args.output)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(header)
        print(f"Wrote {args.output}")
    else:
        print(header)


if __name__ == "__main__":
    main()
