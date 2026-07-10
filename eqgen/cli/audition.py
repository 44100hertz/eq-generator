#!/usr/bin/env python3
"""
Offline audition: speaker measurement → IIR fit → C DSP → processed WAV.

Usage:
    python -m eqgen.cli.audition technics/standing /tmp/out
    python -m eqgen.cli.audition technics/standing /tmp/out --tracks song.flac
    python -m eqgen.cli.audition cardboard /tmp/out --fc 50 --h2 0.5 --h3 1.0
"""

import argparse
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import (
    run_pipeline, design_eq, process_track, pre_gain_from_max_gain,
)
from eqgen.presets import MAX_IIR_BANDS

MEAS_DIR = ROOT / "measurements"


def run_speaker(speaker_name: str, out_dir: str, music_dir: str = None,
                fc: float = 60.0, h2: float = 0.5, h3: float = 1.0,
                max_bands: int = MAX_IIR_BANDS, smooth_exponent: float = 1.0,
                 tracks: list = None):
    """Design EQ for a speaker and process music tracks through it."""
    meas_dir = MEAS_DIR / speaker_name
    if not meas_dir.exists():
        print(f"ERROR: measurement dir not found: {meas_dir}")
        return

    meas_wavs = sorted(meas_dir.glob("measurement*.wav"))
    if not meas_wavs:
        meas_wavs = sorted(meas_dir.glob("response*.wav"))
    target_path = meas_dir / "target.wav"
    if not meas_wavs or not target_path.exists():
        print(f"ERROR: need measurement*.wav and target.wav in {meas_dir}")
        return

    meas_path = str(meas_wavs[-1])
    target = str(target_path)

    print("=" * 70)
    print(f"  SPEAKER: {speaker_name}  |  fc={fc} h2={h2} h3={h3}")
    print(f"  Measurement: {os.path.basename(meas_path)}")
    print("=" * 70)

    # ── EQ pipeline ────────────────────────────────────────────────
    print(f"\n── EQ pipeline (Welch + adaptive points + model)...")
    eq_freqs, target_db, fs, max_gain_db = run_pipeline(
        [meas_path], target,
        bass_enhancer_cutoff=fc, h2=h2, h3=h3)

    print(f"  {len(eq_freqs)} adaptive EQ points, {eq_freqs[0]:.0f}–{eq_freqs[-1]:.0f} Hz")
    pre_gain = pre_gain_from_max_gain(max_gain_db)
    if max_gain_db > 0:
        print(f"  Max correction gain: {max_gain_db:+.1f} dB → pre-gain: {20*np.log10(pre_gain):+.1f} dB")

    # ── IIR fit ────────────────────────────────────────────────────
    print(f"\n── Fitting IIR biquads...")
    coeffs, bands, eq_freqs, target_db, fitted_db = design_eq(
        eq_freqs, target_db, fs, max_bands=max_bands,
        min_peaking_freq=fc)

    print(f"  {len(bands)} bands")
    err = fitted_db - target_db
    for label, lo, hi in [("Bass ≤250Hz", 0, 250),
                           ("Mid 250-2k", 250, 2000),
                           ("Treble >2k", 2000, 99999)]:
        m = (eq_freqs >= lo) & (eq_freqs <= hi)
        if m.any():
            print(f"    {label:<16s}  max err = {np.max(np.abs(err[m])):+.1f} dB")

    # ── Find music ─────────────────────────────────────────────────
    if tracks:
        all_files = [t for t in tracks if os.path.exists(t)]
        if not all_files:
            print("ERROR: none of the specified tracks exist")
            return
        max_tracks = len(all_files)
    else:
        if not music_dir:
            music_dir = os.path.expanduser("~/Music")
        print(f"\n── Finding music in {music_dir}...")
        result = subprocess.run(
            ["find", music_dir, "-type", "f",
             "(", "-name", "*.flac", "-o", "-name", "*.mp3",
             "-o", "-name", "*.wav", "-o", "-name", "*.m4a", ")",
             "-not", "-path", "*/.Trash*"],
            capture_output=True, text=True)
        all_files = [f.strip() for f in result.stdout.split("\n") if f.strip()]
        max_tracks = 3
        random.shuffle(all_files)

    os.makedirs(out_dir, exist_ok=True)

    # ── Process tracks ─────────────────────────────────────────────
    print(f"\n── Processing {min(max_tracks, len(all_files))} track(s)...")
    processed = 0
    for path in (all_files if tracks else all_files[:10]):
        if processed >= max_tracks:
            break
        name = os.path.basename(path)
        safe = os.path.splitext(name)[0].replace(" ", "_")[:50]
        out_path = os.path.join(
            out_dir, f"{speaker_name.replace('/', '_')}_{safe}.wav")
        print(f"  [{processed+1}/{max_tracks}] {name}")
        if process_track(path, out_path, coeffs, len(bands),
                         cutoff_hz=fc, h2=h2, h3=h3, pre_gain=pre_gain):
            processed += 1
        print()

    print(f"  ✅ Output in {out_dir}/")


def main():
    ap = argparse.ArgumentParser(
        description="End-to-end: measure → EQ → C DSP → WAV")
    ap.add_argument("speaker", help="Speaker name (technics/standing, cardboard)")
    ap.add_argument("out_dir", help="Output directory for WAV files")
    ap.add_argument("--music-dir", default=None,
                    help="Directory to scan for music")
    ap.add_argument("--tracks", nargs="*", default=None,
                    help="Specific audio files to process")
    ap.add_argument("--fc", type=float, default=60.0,
                    help="Bass enhancer cutoff Hz [60]")
    ap.add_argument("--h2", type=float, default=0.5,
                    help="2nd harmonic amplitude [0.5]")
    ap.add_argument("--h3", type=float, default=1.0,
                    help="3rd harmonic amplitude [1.0]")
    ap.add_argument("--smooth-exponent", type=float, default=1.0,
                    help="CV smoothing aggressiveness [1.0]")
    ap.add_argument("--max-bands", type=int, default=MAX_IIR_BANDS,
                    help="Max IIR biquad bands [40]")
    args = ap.parse_args()

    run_speaker(args.speaker, args.out_dir,
                music_dir=args.music_dir,
                fc=args.fc, h2=args.h2, h3=args.h3,
                max_bands=args.max_bands, tracks=args.tracks)


if __name__ == "__main__":
    main()
