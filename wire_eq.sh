#!/usr/bin/env bash
# wire_eq.sh — End-to-end: measure speaker → build LADSPA plugin → wire to output.
#
# Usage:
#   ./wire_eq.sh setup <speaker_name>
#   ./wire_eq.sh setup -m response.wav -n noise.wav -t target.wav
#   ./wire_eq.sh teardown
#   ./wire_eq.sh build <speaker_name>
#
# Auto-enters nix develop if nix is available and numpy is missing.
set -euo pipefail
cd "$(dirname "$0")"

if command -v nix &>/dev/null && [ -f flake.nix ]; then
    exec nix develop --command python3 wire_eq.py "$@"
else
    exec python3 wire_eq.py "$@"
fi
