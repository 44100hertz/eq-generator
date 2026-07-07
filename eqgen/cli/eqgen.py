#!/usr/bin/env python3
"""
EQ curve designer: measurements → EQ curve data.

Replaces the old JamesDSP CSV output with plain JSON so the EQ curve
can be fed to the audition / wire / export pipelines.

Usage:
    python -m eqgen.cli.eqgen -m meas.wav -t target.wav -o curve.json
    python -m eqgen.cli.eqgen -m meas.wav -t target.wav --fc 50 --h2 0.5
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import run_pipeline, curve_to_json, db_to_ratio, ratio_to_db


def parse_rolloff(spec: str):
    """Parse a 'FREQ,DB_OCTAVE' rolloff spec."""
    parts = spec.split(",")
    return float(parts[0]), float(parts[1])


def main():
    ap = argparse.ArgumentParser(
        description="Speaker EQ correction suite: measurements → EQ curve (JSON)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m eqgen.cli.eqgen -m meas1.wav meas2.wav -t target.wav -o eq.json
    python -m eqgen.cli.eqgen -m meas.wav -t target.wav --noise noise.wav \\
        --bass-enhancer-cutoff 50.0 --low-rolloff 220,2 --low-rolloff 50,-2
""",
    )

    ap.add_argument("-m", "--measurement", nargs="+", required=True,
                    help="One or more measurement WAV files")
    ap.add_argument("-t", "--target", required=True,
                    help="Target WAV file (reference system recording)")
    ap.add_argument("-n", "--noise", default=None,
                    help="Background noise WAV for spectral subtraction")
    ap.add_argument("-o", "--output", default=None,
                    help="Output JSON file (default: stdout)")
    ap.add_argument("--high-rolloff", action="append", default=[],
                    metavar="FREQ,DB_OCTAVE",
                    help="High-frequency rolloff (e.g. 4000,-6). Repeatable.")
    ap.add_argument("--low-rolloff", action="append", default=[],
                    metavar="FREQ,DB_OCTAVE",
                    help="Low-frequency rolloff (e.g. 100,-12). Repeatable.")
    ap.add_argument("--bass-enhancer-cutoff", "--fc", type=float, default=None,
                    dest="fc", help="Bass enhancer cutoff Hz")
    ap.add_argument("--smooth-exponent", type=float, default=1.0,
                    help="CV smoothing aggressiveness. 0=fixed, 1=linear, 2=amplified [1.0]")
    ap.add_argument("--h2", type=float, default=1.0,
                    help="2nd harmonic amplitude")
    ap.add_argument("--h3", type=float, default=1.0,
                    help="3rd harmonic amplitude")

    args = ap.parse_args()

    high_rolloffs = [parse_rolloff(s) for s in args.high_rolloff]
    low_rolloffs = [parse_rolloff(s) for s in args.low_rolloff]

    freqs, gains_db, sample_rate = run_pipeline(
        args.measurement, args.target, args.noise,
        bass_enhancer_cutoff=args.fc,
        h2=args.h2, h3=args.h3, smooth_exponent=args.smooth_exponent,
        high_rolloffs=high_rolloffs,
        low_rolloffs=low_rolloffs,
    )

    output = {
        "sample_rate": sample_rate,
        "points": curve_to_json(freqs, gains_db),
    }

    json_str = json.dumps(output, indent=2)

    if args.output:
        with open(args.output, "w") as f:
            f.write(json_str)
            f.write("\n")
        print(f"Wrote {args.output}", file=sys.stderr)
    else:
        print(json_str)


if __name__ == "__main__":
    main()
