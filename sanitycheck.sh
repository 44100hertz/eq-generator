#!/usr/bin/env bash
# Run the speaker EQ correction suite.
#
# Requires numpy, scipy (available via `nix develop` or system Python).
# If nix is available, automatically enters the dev shell.
set -euo pipefail
cd "$(dirname "$0")"

if command -v nix &>/dev/null && [ -f flake.nix ]; then
    PYTHON="nix develop --command python"
else
    PYTHON=python
fi

# With bass enhancer preprocessing (standard)
$PYTHON eqgen.py \
  -m measurements/technics/standing/measurement2.wav \
  -t measurements/technics/standing/target.wav \
  --noise measurements/technics/standing/noise2.wav \
  -b 0 \
  --bass-enhancer-cutoff 50.0 \
  --h2 0.5 --h3 0.5 \
  -o sanitycheck

# Without enhancer (for comparison)
$PYTHON eqgen.py \
  -m measurements/technics/standing/measurement2.wav \
  -t measurements/technics/standing/target.wav \
  --noise measurements/technics/standing/noise2.wav \
  -b 0 \
  -o sanitycheck_noenhancer

echo "Done. Outputs: sanitycheck_desktop.csv, sanitycheck_mobile.csv"
echo "               sanitycheck_noenhancer_desktop.csv, sanitycheck_noenhancer_mobile.csv"
