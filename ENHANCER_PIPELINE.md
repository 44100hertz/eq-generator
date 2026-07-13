# Enhancer DSP Pipeline

The bass enhancer (`src/enhancer.c` / `src/enhancer.h`) is the core of EQGen's real-time processing — it runs on the ESP32 firmware and the desktop PipeWire filter. This document describes the signal flow and the rationale behind each stage.

## Signal Chain (per channel)

```
input float → DC blocker → pre-gain → EQ cascade → LR4 crossover → crossfade → loudness shelf → tanh clamp
                                                         │
                                                    lookahead HP
                                                    (ring buffer)
```

### 1. DC Blocker
First-order IIR HP at 5 Hz. Removes DC offset before any gain or filtering.

### 2. Pre-gain
Applied before the EQ to prevent the biquad cascade from overflowing. Set to `1 / max(EQ gain)` so that the signal *after* the EQ at the maximum-gain frequency is roughly at the level of the original signal before vol-scaling.

### 3. EQ Cascade
32 cascaded biquads (Direct Form I) implementing the speaker correction curve. Designed by `eqgen.eq_fit` using iterative greedy matching. Linear — generates no distortion of its own.

### 4. LR4 Crossover
**Replaced the original 1st-order complementary LP/HP split.** The old split (`HP = 1 − LP`) produced a 40% HP component at 31 Hz due to phase shift, which the lookahead interpreted as actual high-frequency content. This falsely reduced the headroom budget and triggered harmonic generation on pure bass tones.

The LR4 uses two cascaded 2nd-order Butterworth sections each for LP and HP at the crossover frequency (`cutoff_hz`, typically 70 Hz). Properties:

| Frequency | LP bleed | HP bleed | Sum magnitude |
|-----------|----------|----------|---------------|
| Below fc  | ~100%    | ~4%      | 1.0 (allpass) |
| At fc     | 50%      | 50%      | 1.0 (allpass) |
| Above fc  | ~4%      | ~100%    | 1.0 (allpass) |

The sum `LP + HP = H_allpass` with `|H_allpass| = 1` at all frequencies — flat response, no notch.

* `lp_cross[2]` — 2 cascaded LP biquads for the bass path (delayed)
* `hp_cross[2]` — 2 cascaded HP biquads for the highs path (delayed)
* `hp_lookahead[2]` — 2 cascaded HP biquads for the lookahead (undelayed)

### 5. Lookahead & Headroom Budget

The lookahead runs the current sample through an independent HP cascade (`hp_lookahead`) and stores the absolute value in a 240-sample ring buffer (~5 ms). The peak over the window becomes `upcoming`:

```
room = push_gain × (1 − upcoming)
```

The headroom budget tracks upcoming HF peaks so the crossfade can pre-emptively shift bass energy to harmonics before the peak arrives through the 240-sample delay line.

### 6. Envelope & Normalization

The delayed LP signal feeds a peak-hold envelope follower (instant attack, exponential release τ = 200 ms). A single-pole smoother (α = 0.002) kills the 2f ripple that would otherwise AM-modulate the Chebyshev normalization, creating buzz.

**Adaptive blend:** During steady state the smoothed envelope kills 2f ripple (ripple amplitude ≈ 0.077·A < 0.1 threshold). During transients (sweeps, attacks) the envelope jumps by >0.3 and the blend switches to near-instant tracking, preventing `norm = lp_fund / env_norm` from spiking above 1.0. A norm spike would kick the Chebyshev polynomials into their unbounded region (|x| > 1), producing an impulse that rings through the HP harmonic filter as a crackle.

### 7. Crossfade

The crossfade weight `w` determines how much bass perceived loudness comes from harmonics vs fundamental:

```
target = min(env, 1.0)           // perceived bass level (capped at tanh ceiling)
w     = (target − room)          // fraction from harmonics
        / (target × (1 − h))
```

When `target ≤ room`: pure fundamental (`w = 0`).
When `target > room`: blend toward harmonics.

Both `w` and `harm_amp = target − room` are **slew-limited** at 0.002/sample (~500 samples for full transition) to prevent the Chebyshev from turning on/off in one sample during sweeps.

### 8. Chebyshev Harmonic Generator

* T2(x) = 2x² − 1 — generates 2nd harmonic from normalized bass signal
* T3(x) = 4x³ − 3x — generates 3rd harmonic, with analytic fundamental cancellation

The T3 cancellation is scaled by `harm_amp / env` so it tracks the slewed amplitude — otherwise the full cancellation fires even when `harm_amp ≈ 0`, injecting a phantom fundamental into the HP filter that rings as a crackle.

### 9. Output Clamp

`out = tanh(dry_hp + fund_out + harm_out)` — soft saturator. Produces benign odd-harmonic distortion instead of hard-clipping when the overboost pushes signal past unity.

## Known Limitations

- Very aggressive non-musical signals (e.g., 100% volume sines rapidly sweeping up and down through the crossover region) can produce minor transient artifacts. This is an inherent limitation of the lookahead ring buffer's 5 ms window — no realistic musical content triggers it.
- The LR4 allpass phase shift means the crossover has frequency-dependent group delay. This is inaudible for music but measurable in impulse responses.
