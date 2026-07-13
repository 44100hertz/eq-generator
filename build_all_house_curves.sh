#!/usr/bin/env bash
# Build ESP32 firmware .bin for every house curve based on lunchbox-morebass preset.
# Output goes to /run/media/samp/BBD7-44A9/lunchbox-presets/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="/run/media/samp/BBD7-44A9/lunchbox-presets"
PRESETS_DIR="$SCRIPT_DIR/presets"
TEMPLATE="$PRESETS_DIR/lunchbox-morebass.json"
HOUSE_CURVES="$PRESETS_DIR/house_curves.json"

mkdir -p "$OUTDIR"

# Extract all house curve names from house_curves.json
CURVES=$(python3 -c "
import json
with open('$HOUSE_CURVES') as f:
    curves = json.load(f)
for name in curves:
    print(name)
")

COUNT=0
for curve in $CURVES; do
    COUNT=$((COUNT + 1))
    PRESET_NAME="lunchbox-morebass-$curve"
    PRESET_FILE="$PRESETS_DIR/$PRESET_NAME.json"

    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  [$COUNT] Building for house curve: $curve"
    echo "  Preset: $PRESET_NAME"
    echo "════════════════════════════════════════════════════════════════"

    # Create temp preset with this house curve
    python3 -c "
import json
with open('$TEMPLATE') as f:
    preset = json.load(f)
preset['name'] = '$PRESET_NAME'
preset['house_curve'] = '$curve'
with open('$PRESET_FILE', 'w') as f:
    json.dump(preset, f, indent=2)
    f.write('\n')
"

    # Step 1: Run wire build (generates eq_coeffs.h + builds DSP)
    echo "  → Running wire build..."
    cd "$SCRIPT_DIR"
    python -m eqgen.cli.wire build "$PRESET_NAME"

    # Step 2: Build ESP32 firmware
    echo "  → Building ESP32 firmware..."
    cd "$SCRIPT_DIR/firmware"
    idf.py build

    # Step 3: Copy .bin to output
    BIN_FILE="$SCRIPT_DIR/firmware/build/eqgen.bin"
    OUT_FILE="$OUTDIR/$PRESET_NAME.bin"
    cp "$BIN_FILE" "$OUT_FILE"
    echo "  ✅ Copied to $OUT_FILE"

    # Clean up temp preset
    rm -f "$PRESET_FILE"
done

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Done! Built $COUNT firmware images:"
echo "════════════════════════════════════════════════════════════════"
ls -lh "$OUTDIR"/*.bin
