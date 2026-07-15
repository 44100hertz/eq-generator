# Multi-Driver Active Crossover

## Current state

eqgen processes one full-range speaker: measurement → one correction curve → one biquad cascade, shared
across L/R channels. The enhancer's internal 60 Hz LR4 split is for bass harmonic synthesis only — it has
nothing to do with speaker crossovers.

## Target

Two-way active speaker: EQ'd tweeter on L channel, EQ'd woofer on R channel, mono-summed input.
Crossover and per-driver EQ run on the ESP32.

```
L_in ───→ DC → pre-gain → [LR4 HP] → [EQ_tweeter] → L_out  (tweeter)
R_in ─┘
        → DC → pre-gain → [LR4 LP] → [EQ_woofer]  → [bass enhancer] → R_out  (woofer)
```

No phase modeling. LR4 crossover guarantees identical phase on both outputs; the acoustic sum at
crossover is flat. Per-driver phase is baked into the measurement — measure *through* the crossover
filter and EQ what comes out.

## Measurement protocol

1. Flash firmware with **crossover biquads only** (no EQ), L/R outputs active
2. Run a sweep: physically disconnect or mute the **woofer** → record tweeter-only response
3. Run a sweep: disconnect/mute the **tweeter** → record woofer-only response
4. Run eqgen pipeline **twice** — two independent correction curves
5. Flash merged firmware (crossover + both EQ cascades)

Alternatively: split a single measurement WAV into two pipelines by applying the crossover filters in
software post-hoc. Avoids multiple firmware flashes. Slightly less accurate because the physical amp +
driver chain isn't captured identically, but close enough for initial tuning.

## Phase 1 — Python (pipeline, no firmware changes)

### 1a. Crossover biquad design (`eq_fit.py`)

Add `design_butter_lp(fc, fs)` → `BiquadCoeffs`. Already have `design_butter_hp`. Then:

```python
def design_lr4_crossover(fc, fs):
    """Two 2nd-order Butterworth sections each → 4th-order LR."""
    lp1 = design_butter_lp(fc, fs)
    lp2 = design_butter_lp(fc, fs)
    hp1 = design_butter_hp(fc, fs)
    hp2 = design_butter_hp(fc, fs)
    return [lp1, lp2], [hp1, hp2]
```

### 1b. Multi-channel pipeline (`pipeline.py`)

Replace single `run_pipeline(measurement_paths, target_path, ...)` with:

```python
@dataclass
class ChannelConfig:
    name: str              # "tweeter" | "woofer"
    measurement_paths: list[str]
    target_path: str | None   # None → flat
    house_curve: list | None
    bass_enhancer_cutoff: float | None
    max_bands: int

def run_crossover_pipeline(
    crossover_fc: float,
    channels: list[ChannelConfig],
    ...
) -> CrossoverResult:
```

Internally: runs the existing Welch → CV smoothing → correction → IIR fit pipeline for each channel
independently.

### 1c. Merged output

```python
@dataclass
class CrossoverResult:
    crossover_fc: float
    lp_biquads: list[BiquadCoeffs]   # 2x 2nd-order LP
    hp_biquads: list[BiquadCoeffs]   # 2x 2nd-order HP
    tweeter_biquads: list[BiquadCoeffs]
    woofer_biquads: list[BiquadCoeffs]
    tweeter_pre_gain: float
    woofer_pre_gain: float
    sample_rate: float
```

Export as C header with separate coefficient arrays for HP path, LP path, tweeter EQ, woofer EQ.
Single `tweeter.measurements/` and `woofer.measurements/` directories in the preset structure.

## Phase 2 — Firmware

### 2a. Restructure `BassEnhancerCfg` / `BassEnhancer`

Current: one `eq_coeffs` array, shared across L/R. Enhancer runs on both channels.

New: separate EQ for tweeter (L) and woofer (R). Enhancer runs on woofer channel only (bass synthesis
makes no sense above the speaker crossover).

```c
typedef struct {
    /* ... existing fields ... */
    int eq_n_biquads;                /* max of tweeter/woofer */
    const float *eq_tweeter_coeffs;  /* tweeter EQ (L channel, HP path) */
    const float *eq_woofer_coeffs;   /* woofer EQ  (R channel, LP path) */
} BassEnhancerCfg;
```

### 2b. Signal path rewrite

```
L channel (tweeter):
  1. DC blocker
  2. pre_gain
  3. [speaker LR4 HP crossover] (2 biquads)
  4. [tweeter EQ biquads]
  5. Output (no enhancer — skip bass section)

R channel (woofer):
  1. DC blocker
  2. pre_gain
  3. [speaker LR4 LP crossover] (2 biquads)
  4. [woofer EQ biquads]
  5. [bass enhancer: LP/HP split, harmonics, mix, loudness]
  6. Output
```

The enhancer's internal LR4 (60 Hz) stays — it sits *after* the woofer EQ in the R channel.

### 2c. Crossover biquad initialization

Add to `BassEnhancerCfg`:

```c
float speaker_lp_coeffs[LR4_SECTIONS * 5];  /* speaker crossover LP */
float speaker_hp_coeffs[LR4_SECTIONS * 5];  /* speaker crossover HP */
```

Per-channel state gets `Biquad speaker_cross[LR4_SECTIONS]`.

`BassEnhancerCfg_init` grows a `speaker_crossover_hz` parameter. Zero or negative → disabled (full-range
mode, backward compatible).

## Phase 3 — Preset format (`presets.py`)

Extend `Preset` to optionally include:

```python
speaker_config: dict | None  # {"crossover_hz": 2500, "channels": {"tweeter": ..., "woofer": ...}}
```

When `speaker_config` is present, `run_crossover_pipeline` is used. When absent, existing single-channel
path — full backward compatibility.

## Verification

1. Unit test: LR4 HP + LP biquad cascade response sums to 0 dB at all frequencies (±0.1 dB)
2. End-to-end: measure a known 2-way speaker (e.g., cardboard box test) → run crossover pipeline →
   verify combined acoustic response is flat within ±2 dB across the crossover region
3. Firmware: run standalone with sine sweeps on each channel, verify HP/LP rolloff at crossover freq
