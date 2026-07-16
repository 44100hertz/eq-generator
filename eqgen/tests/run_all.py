#!/usr/bin/env python3
"""
Run all eqgen tests in sequence.

Usage:
    python run_all.py              # run everything
    python run_all.py --skip HEURISTIC  # skip heuristic tests
    python run_all.py --skip INTEG      # skip integration tests
    python run_all.py --only roundtrip  # run only a specific module
"""

import sys
import argparse
import time
from pathlib import Path

ROOT = str(Path(__file__).resolve().parent.parent.parent)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# ── Test modules (in run order) ────────────────────────────────────────

MODULES = {
    # Heuristic tests — verify algorithm correctness independent of the DSP pipeline
    "eq_pipeline":        ("eqgen.tests.test_eq_pipeline",       "EQ curve → IIR fit"),
    "fullband_fit":       ("eqgen.tests.test_fullband_fit",      "Full-band IIR fit (synthetic)"),
    "real_fullband_fit":  ("eqgen.tests.test_real_fullband_fit", "Full-band IIR fit (real data)"),

    # Integration tests — full pipeline end-to-end through the C DSP
    "roundtrip":          ("eqgen.tests.test_roundtrip",         "Round-trip pipeline"),
    "input_signals":      ("eqgen.tests.test_input_signals",     "Input signal tones through C DSP"),
}

HEURISTIC_MODULES = {"eq_pipeline", "fullband_fit", "real_fullband_fit"}
INTEG_MODULES     = {"roundtrip", "input_signals"}


def main():
    parser = argparse.ArgumentParser(
        description="Run all eqgen tests."
    )
    parser.add_argument(
        "--skip", nargs="*", default=[],
        choices=["HEURISTIC", "INTEG"] + list(MODULES.keys()),
        help="Skip specific groups or modules (HEURISTIC, INTEG, or module names)"
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
        print("Heuristic tests (algorithm correctness):")
        for name in HEURISTIC_MODULES:
            print(f"    {name:20s}  {MODULES[name][1]}")
        print(f"\nIntegration tests (full pipeline):")
        for name in INTEG_MODULES:
            print(f"    {name:20s}  {MODULES[name][1]}")
        return

    if args.only:
        selected = set(args.only)
    else:
        selected = set(MODULES.keys())

    for skip in args.skip:
        if skip == "HEURISTIC":
            selected -= HEURISTIC_MODULES
        elif skip == "INTEG":
            selected -= INTEG_MODULES
        else:
            selected.discard(skip)

    ordered = [name for name in MODULES if name in selected]

    if not ordered:
        print("No modules selected to run.")
        return

    total_start = time.time()
    passed = 0
    failed = 0

    for i, name in enumerate(ordered):
        module_path, description = MODULES[name]
        header = f" [{i+1}/{len(ordered)}] {description} "

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
