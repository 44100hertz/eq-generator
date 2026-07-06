# eqgen firmware

ESP32 Bluetooth A2DP sink with DSP (harmonic bass enhancer + EQ biquads).

## Build

Requires ESP-IDF v4.4 or later.

```bash
# 1. Generate EQ coefficients from the desktop tuning tools
cd ..
python -m eqgen.cli.export --speaker small --fc 60 --h2 0.33 --h3 0.33 \
    -o src/eq_coeffs.h

# 2. Build and flash
cd firmware
idf.py set-target esp32
idf.py build
idf.py flash
```

After flashing, the device appears as **"eqgen"** in your phone's Bluetooth menu.
Pair and play audio — it flows through the DSP and out the I2S DAC.

## Hardware

| ESP32 pin | PCM5102A pin | Signal |
|---|---|---|
| GPIO 26 | BCK | Bit clock |
| GPIO 25 | LRCK / WS | Word select |
| GPIO 22 | DIN | Data |
| GND | GND | Ground |

The PCM5102A module gets 5V from the ESP32 dev board's VIN pin (or USB).
No MCLK, no mute pins needed — the DAC auto-starts on the first BCK edge.

## Changing the EQ

1. Tune on desktop with the PipeWire plugin or offline audition tools
2. Re-run `python -m eqgen.cli.export -o src/eq_coeffs.h`
3. `idf.py build flash` — done

The DSP code never changes. Only the coefficients header gets regenerated.
