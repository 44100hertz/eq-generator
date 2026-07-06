#!/usr/bin/env python3
"""
Standalone enhancer: process audio through the C DSP (flat EQ).

No speaker measurement needed — uses a default 3-HP-cascade EQ.

Usage:
    python -m eqgen.cli.enhance ~/Music /tmp/out
"""

import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import build_default_eq_coeffs, process_track


def main():
    music_dir = sys.argv[1] if len(sys.argv) > 1 else None
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "/tmp/eqgen_out"

    if not music_dir:
        print("Usage: python -m eqgen.cli.enhance <music_dir> [output_dir]")
        print("  Processes 3 random tracks through the C bass enhancer")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)

    music_dir_abs = os.path.abspath(music_dir)
    result = subprocess.run(
        ["find", music_dir_abs, "-type", "f",
         "(", "-name", "*.flac", "-o", "-name", "*.mp3",
         "-o", "-name", "*.wav", "-o", "-name", "*.m4a", ")",
         "-not", "-path", "*/.Trash*", "-not", "-path", "*/._*"],
        capture_output=True, text=True)
    all_files = [f.strip() for f in result.stdout.split("\n") if f.strip()]

    if not all_files:
        print(f"No audio files found in {music_dir}")
        sys.exit(1)

    random.shuffle(all_files)
    print(f"Found {len(all_files)} audio files. Processing 3 random tracks...\n")

    eq_coeffs = build_default_eq_coeffs()
    processed = 0
    for path in all_files[:10]:
        if processed >= 3:
            break
        name = os.path.basename(path)
        safe_name = os.path.splitext(name)[0].replace(" ", "_")[:60]
        out_path = os.path.join(out_dir, safe_name + "_enhanced.wav")

        print(f"[{processed+1}/3] {os.path.basename(os.path.dirname(path))}/{name}")
        if process_track(path, out_path, eq_coeffs, len(eq_coeffs) // 5):
            processed += 1
        print()

    print(f"Done! Output in {out_dir}/")


if __name__ == "__main__":
    main()
