#!/usr/bin/env bash
# Run all python analysis scripts.
#
# Requires numpy (available via `nix develop` or your system python).
# If nix is available, automatically enters the dev shell.
#
# Usage:
#   ./run_analysis.sh              # run everything
#   ./run_analysis.sh --skip MODEL  # skip model analyses
#   ./run_analysis.sh --only linearity  # run a single module
#   ./run_analysis.sh --list        # list available modules
set -euo pipefail
cd "$(dirname "$0")"

if command -v nix &>/dev/null && [ -f flake.nix ]; then
    exec nix develop --command python python/run_all.py "$@"
else
    exec python python/run_all.py "$@"
fi
