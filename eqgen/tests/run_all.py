#!/usr/bin/env python3
"""
Run all Python analysis scripts in sequence.

Usage:
    python run_all.py              # run everything
    python run_all.py --skip MODEL  # skip model-based analyses
    python run_all.py --skip DSP    # skip DSP tests
    python run_all.py --only harmonics  # run only a specific module
"""

import sys
import argparse
import time
from pathlib import Path

# Add the project root to sys.path so modules can import from eqgen.*
ROOT = str(Path(__file__).resolve().parent.parent.parent)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# ── Analysis modules (in run order) ────────────────────────────────────
# Each module exposes a run() function.

MODULES = {
    # DSP tests — verify the enhancer's signal-processing behavior
    "harmonics":   ("eqgen.tests.test_harmonics",   "Harmonic purity & linearity"),
    "compressor":  ("eqgen.tests.test_compressor",  "Compressor & dynamics"),
    "gain_safety": ("eqgen.tests.test_gain_safety", "Gain safety on real measurements"),
    "chebyshev":   ("eqgen.tests.test_chebyshev",   "Chebyshev math analysis"),

    # Model tests — verify the psychoacoustic model
    "model":       ("eqgen.tests.test_model",       "Model verification (DSP vs formula)"),
    "loudness":    ("eqgen.model",                  "Equal-loudness weighting"),
    "model_gain":  ("eqgen.model",                  "Model gain analysis"),
}

# ── Group definitions ──────────────────────────────────────────────────
DSP_MODULES    = {"harmonics", "compressor", "gain_safety"}
MODEL_MODULES  = {"chebyshev", "model", "loudness", "model_gain"}


def main():
    parser = argparse.ArgumentParser(
        description="Run all Python analysis scripts for eqgen."
    )
    parser.add_argument(
        "--skip", nargs="*", default=[],
        choices=["DSP", "MODEL"] + list(MODULES.keys()),
        help="Skip specific groups or modules (DSP, MODEL, or module names)"
    )
    parser.add_argument(
        "--only", nargs="*", default=None,
        choices=list(MODULES.keys()),
        help="Run only the specified modules"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available modules and exit"
    )
    args = parser.parse_args()

    if args.list:
        print("Available modules:")
        print(f"\n  DSP tests (verify signal processing):")
        for name in DSP_MODULES:
            print(f"    {name:20s}  {MODULES[name][1]}")
        print(f"\n  Model tests (verify psychoacoustic model):")
        for name in MODEL_MODULES:
            print(f"    {name:20s}  {MODULES[name][1]}")
        return

    # Determine which modules to run
    if args.only:
        selected = set(args.only)
    else:
        selected = set(MODULES.keys())

    # Apply skips
    for skip in args.skip:
        if skip == "DSP":
            selected -= DSP_MODULES
        elif skip == "MODEL":
            selected -= MODEL_MODULES
        else:
            selected.discard(skip)

    # Sort in the order defined above
    ordered = [name for name in MODULES if name in selected]

    if not ordered:
        print("No modules selected to run.")
        return

    # ── Run each module ────────────────────────────────────────────────
    total_start = time.time()
    passed = 0
    failed = 0

    for i, name in enumerate(ordered):
        module_path, description = MODULES[name]
        header = f" [{i+1}/{len(ordered)}] {description} "

        # Special case: the model.py library has multiple run functions
        if module_path == "eqgen.model":
            import eqgen.model as _mod
            if name == "loudness":
                fn = _mod.run_loudness_analysis
            elif name == "model_gain":
                fn = _mod.run_model_gain_analysis
            else:
                fn = _mod.run_loudness_analysis  # shouldn't happen
        else:
            _mod = __import__(module_path, fromlist=["run"])
            fn = _mod.run

        print()
        print("=" * 70)
        print(f"{header:=^70}")
        print("=" * 70)

        try:
            mod_start = time.time()
            fn()
            mod_elapsed = time.time() - mod_start
            passed += 1
            print(f"\n  ✅ {description} completed in {mod_elapsed:.1f}s")
        except Exception as e:
            failed += 1
            print(f"\n  ❌ {description} FAILED: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()

    total_elapsed = time.time() - total_start

    # ── Summary ────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print(f"{' RESULTS ':=^70}")
    print("=" * 70)
    print(f"  Total:   {passed + failed} modules")
    print(f"  Passed:  {passed}")
    print(f"  Failed:  {failed}")
    print(f"  Time:    {total_elapsed:.1f}s")
    print("=" * 70)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
