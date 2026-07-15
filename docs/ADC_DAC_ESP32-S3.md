# ADC → DSP → DAC on ESP32-S3

Wired I2S audio pipeline coexisting with the existing Bluetooth A2DP path.

## Architecture

```
PCM1808 ──I2S──→ ESP32-S3 ──I2S──→ PCM5102A
  (ADC)      RX │ DSP │ TX       (DAC)
                └─ enhancer + EQ
```

Full-duplex on I2S_NUM_0: one `i2s_new_channel` call creates both TX and RX handles. Shared BCLK/WS, separate data lines. MCLK output from the ESP32 feeds the PCM1808's SCK — no external oscillator.

## Build-time source selection

Kconfig `choice EQGEN_AUDIO_SOURCE`:

| Value | Target | Behavior |
|-------|--------|----------|
| `BT_A2DP` | esp32 (default) | Existing BT sink → StreamBuffer → DSP → I2S TX |
| `I2S_ADC` | esp32s3 (default) | I2S RX DMA → DSP → I2S TX DMA |

```kconfig
# Kconfig.projbuild (new)
choice EQGEN_AUDIO_SOURCE
    prompt "Audio source"
    default EQGEN_AUDIO_SOURCE_BT_A2DP if IDF_TARGET_ESP32
    default EQGEN_AUDIO_SOURCE_I2S_ADC if IDF_TARGET_ESP32S3

config EQGEN_AUDIO_SOURCE_BT_A2DP
    bool "Bluetooth A2DP"
    select BT_ENABLED
config EQGEN_AUDIO_SOURCE_I2S_ADC
    bool "I2S ADC (PCM1808)"
endchoice
```

Per-target overrides (layered on shared `sdkconfig.defaults`):

```
sdkconfig.defaults              # common: partition, perf, toolchain
sdkconfig.defaults.esp32        # BT enabled, BTDM classic-only
sdkconfig.defaults.esp32s3      # S3 CPU freq, no BT
```

## New files

### `main/i2s_in.c` + `main/i2s_in.h`
PCM1808 RX driver. Gist:
- `i2s_in_init(sample_rate)` — RX channel on I2S_NUM_0, slave, standard Philips, 24-bit → 32-bit slot
- `i2s_in_read(buf, frame_count, &bytes_read)` — blocking DMA read

### `main/i2s_duplex.c` + `main/i2s_duplex.h`
Combined init replacing the separate TX-only init when `I2S_ADC`:
```c
i2s_chan_config_t chan_cfg = { .id = I2S_NUM_0, .role = I2S_ROLE_MASTER };
i2s_new_channel(&chan_cfg, &tx_handle, &rx_handle);
// tx_handle → PCM5102A, rx_handle → PCM1808
```
MCLK on TX GPIO config (both channels share the clock tree).

### `main/audio_source.h`
Thin abstraction:
```c
typedef struct {
    int  (*init)(int sample_rate);
    int  (*read)(int16_t *buf, uint32_t frames, size_t *bytes_read);
    int  (*get_volume)(void);
    void (*deinit)(void);
} AudioSource;
extern const AudioSource *audio_source;
```

### `main/audio_source_bt.c`
Wraps existing `xStreamBufferReceive` + `bt_a2dp_get_volume()`.

### `main/audio_source_i2s.c`
Wraps `i2s_in_read`. Fixed volume (compile-time constant, or ADC pot later).

## Modified files

| File | Change |
|------|--------|
| `CMakeLists.txt` | Document both targets |
| `main/CMakeLists.txt` | `SRCS` conditional on `CONFIG_EQGEN_AUDIO_SOURCE_*`; link `i2s_in.c` or `bt_a2dp.c`; `REQUIRES bt` only for A2DP |
| `main/main.c` | `audio_data_handler` → `audio_source->read()`; same DSP, same batching, same I2S TX write. Drop BT-specific volume LUT when `I2S_ADC`. |
| `README.md` | Wiring tables for both modes |

## Pins (I2S_ADC mode)

| Signal | GPIO | PCM1808 pin | PCM5102A pin |
|--------|------|-------------|--------------|
| MCLK   | 1    | SCK         | —            |
| BCLK   | 4    | BCK         | BCK          |
| LRCK   | 5    | LRC         | LRCK         |
| DOUT   | 6    | —           | DIN          |
| DIN    | 7    | OUT         | —            |

PCM1808 mode pins: FMT=GND, MD1=GND, MD0=GND (I2S, slave, 24-bit).

## DSP loop (unchanged)

```
audio_source->read(buf)
  → int16→float (÷32768)
  → BassEnhancer_process_stereo
  → float→int16 (×32768, error-feedback dither)
  → i2s_out_write
```

Same batching (BATCH_SIZE=256), same rate-change reinit, same enhancer reset on stream start.

## What stays untouched

- `src/enhancer.c`, `src/biquad.h`, `src/eq_coeffs.h`
- `sfx_player.c` (used in BT mode)
- All existing Kconfig defaults for esp32 target
