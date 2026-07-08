#!/usr/bin/env python3
"""
Preset CLI — manage EQGen presets from the command line.

Usage:
    python -m eqgen.cli.preset list
    python -m eqgen.cli.preset show <name>
    python -m eqgen.cli.preset create <name> [--from-dir DIR]
    python -m eqgen.cli.preset delete <name>
    python -m eqgen.cli.preset set <name> --fc 50 --h2 0.5 --h3 1.0
    python -m eqgen.cli.preset run <name>
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.presets import PresetManager, Preset


def cmd_list(pm: PresetManager):
    """List all presets."""
    names = pm.list_presets()
    if not names:
        print("No presets found.")
        return
    print(f"{'NAME':<30} {'FC':>6} {'H2':>6} {'H3':>6}  DESCRIPTION")
    print("-" * 70)
    for name in names:
        try:
            p = pm.load(name)
            fc = f"{p.fc:.0f}" if p.fc else ""
            print(f"{p.name:<30} {fc:>6} {p.h2:>6.2f} {p.h3:>6.2f}  {p.description or ''}")
        except Exception as e:
            print(f"{name:<30} {'ERROR':>6}  {e}")


def cmd_show(pm: PresetManager, name: str):
    """Show a preset's details."""
    p = pm.load(name)
    print(json.dumps(p.to_dict(), indent=2))


def cmd_create(pm: PresetManager, name: str, from_dir: str = None):
    """Create a new preset."""
    if from_dir:
        preset = pm.guess_preset_from_measurement_dir(from_dir)
        if preset is None:
            print(f"ERROR: no WAV files found in measurements/{from_dir}")
            sys.exit(1)
        preset.name = name
    else:
        preset = Preset(name=name, description=f"Preset: {name}")

    path = pm.save(preset)
    print(f"Created: {path}")


def cmd_delete(pm: PresetManager, name: str):
    """Delete a preset."""
    if pm.delete(name):
        print(f"Deleted: {name}")
    else:
        print(f"Preset '{name}' not found")
        sys.exit(1)


def cmd_set(pm: PresetManager, name: str, **kwargs):
    """Update preset parameters."""
    try:
        p = pm.load(name)
    except FileNotFoundError:
        print(f"Preset '{name}' not found")
        sys.exit(1)

    # Update fields from kwargs
    for key, value in kwargs.items():
        if value is not None and hasattr(p, key):
            setattr(p, key, value)

    pm.save(p)
    print(f"Updated: {name}")
    print(json.dumps(p.to_dict(), indent=2))


def cmd_run(pm: PresetManager, name: str):
    """Run the pipeline for a preset and print results."""
    try:
        p = pm.load(name)
    except FileNotFoundError:
        print(f"Preset '{name}' not found")
        sys.exit(1)

    meas_paths = p.resolve_measurements()
    target_path = p.resolve_target()
    noise_path = p.resolve_noise()

    if not meas_paths or not target_path:
        print("ERROR: preset has no measurement or target paths")
        sys.exit(1)

    from eqgen.pipeline import run_pipeline

    print(f"Running pipeline for '{name}'...")
    freqs, gains_db, sample_rate, max_gain_db = run_pipeline(
        meas_paths, target_path, noise_path,
        bass_enhancer_cutoff=p.fc, h2=p.h2, h3=p.h3,
        smooth_exponent=p.smooth_exponent,
    )

    print(f"\nResults:")
    print(f"  Sample rate: {sample_rate:.0f} Hz")
    print(f"  EQ points:   {len(freqs)} ({freqs[0]:.0f} – {freqs[-1]:.0f} Hz)")
    print(f"  Max gain:    {max_gain_db:+.1f} dB")
    print(f"\nTo export: python -m eqgen.cli.eqgen -m {' '.join(meas_paths)} -t {target_path} -o eq.json")


def main():
    ap = argparse.ArgumentParser(
        description="EQGen preset management",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    eqgen-preset list
    eqgen-preset show technics-standing
    eqgen-preset create my-speaker --from-dir technics/standing
    eqgen-preset set my-speaker --fc 50 --h2 0.5
    eqgen-preset run my-speaker
        """,
    )

    sub = ap.add_subparsers(dest="command")

    # list
    sub.add_parser("list", help="List all presets")

    # show
    p_show = sub.add_parser("show", help="Show preset details")
    p_show.add_argument("name", help="Preset name")

    # create
    p_create = sub.add_parser("create", help="Create a new preset")
    p_create.add_argument("name", help="Preset name")
    p_create.add_argument("--from-dir", default=None,
                          help="Guess preset from measurements/<dir>")

    # delete
    p_delete = sub.add_parser("delete", help="Delete a preset")
    p_delete.add_argument("name", help="Preset name")

    # set
    p_set = sub.add_parser("set", help="Update preset parameters")
    p_set.add_argument("name", help="Preset name")
    p_set.add_argument("--fc", type=float, default=None, help="Bass enhancer cutoff")
    p_set.add_argument("--h2", type=float, default=None, help="2nd harmonic amplitude")
    p_set.add_argument("--h3", type=float, default=None, help="3rd harmonic amplitude")
    p_set.add_argument("--max-bands", type=int, default=None, help="Max IIR biquad bands")
    p_set.add_argument("--smooth-exponent", type=float, default=None, help="CV smoothing")
    p_set.add_argument("--release", type=float, default=None, help="Envelope release (s)")
    p_set.add_argument("--description", type=str, default=None, help="Preset description")

    # run
    p_run = sub.add_parser("run", help="Run pipeline for a preset")
    p_run.add_argument("name", help="Preset name")

    args = ap.parse_args()
    pm = PresetManager()

    if args.command == "list":
        cmd_list(pm)
    elif args.command == "show":
        cmd_show(pm, args.name)
    elif args.command == "create":
        cmd_create(pm, args.name, args.from_dir)
    elif args.command == "delete":
        cmd_delete(pm, args.name)
    elif args.command == "set":
        cmd_set(pm, args.name,
                fc=args.fc, h2=args.h2, h3=args.h3,
                max_bands=getattr(args, 'max_bands', None),
                smooth_exponent=args.smooth_exponent,
                release=args.release,
                description=args.description)
    elif args.command == "run":
        cmd_run(pm, args.name)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
