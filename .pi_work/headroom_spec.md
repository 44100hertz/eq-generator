# Headroom-based enhancer redesign

## Signal flow (per channel)

```
eq_out ──┬──► HP_lookahead ──► hp_ring[240] ──► upcoming_peak
         │
         └──► delay_line[240]
                   │
                   ▼   (delayed by 240 samples = 5ms @ 48kHz)
              ┌────┴────────────────────┐
              │  LP(fc) ──► env        │
              │  norm = lp/env         │
              │                        │
              │  bleed = clip(         │
              │    lp * bleed_cap,     │
              │    room)               │
              │                        │
              │  harm = clip(          │
              │    Chebyshev(norm,     │
              │      h2_amp, h3_amp),  │
              │    room - |bleed|)     │
              └────────────────────────┘
                   │
            ┌──────┴──────┐
            │ room = push_gain * max(0, 1.0 - upcoming_peak)
            └─────────────┘
                   ▼
            dry_hp + bleed + HP(harm) → loudness shelf → out
```

## Key design rules

1. **Volume BEFORE enhancer** (already done). Volume attenuation creates headroom naturally.
2. **LP runs on eq_out AFTER EQ correction** — already correct in current code.
3. **No reactive limiter** — replaced by lookahead headroom budget.
4. **No volume-driven harmonic↔bleed crossfade** — bleed + harmonics fill whatever room exists.

## What goes away

- AGC harmonic limiter (env_lim, limiter_release_secs, limiter_release_coeff)
- smart_volume harmonic↔bleed crossfade (t, crossfade_lo, crossfade_hi)
- pipeline.py step 8 model_gain_needed EQ modification
- EQGEN_H2_AMP/H3_AMP as eq_coeffs.h constants (they come from calculation now)

## What's new

- LOOKAHEAD_LEN = 240 samples (5ms @ 48kHz)
- delay line ring buffer per channel
- hp_ring ring buffer per channel (absolute values of dry_hp for lookahead window)
- upcoming_peak = max of hp_ring
- push_gain (0.0 = no fill, 1.0 = fill to 0dBFS, can go above 1.0)
- h2_amp, h3_amp, bleed_cap computed from measurement data
- h2_amp = Σ[ear(2f) × correction_ratio(2f)] / Σ[ear(f) × correction_ratio(f)] × human_factor
- h3_amp = same at 3f
- bleed_cap = computed from speaker headroom at bass frequencies

## Files to change

| File | Changes |
|------|---------|
| src/enhancer.h | Remove limiter fields, add lookahead state, push_gain |
| src/enhancer.c | Rewrite enhancer_process_channel with delay+ring+headroom |
| src/smart_volume.h | Strip crossfade; h2/h3/bleed become runtime-constant (passed from cfg) |
| eqgen/model.py | Add compute_harmonic_efficacy() |
| eqgen/pipeline.py | Replace step 8 with efficacy computation |
| eqgen/cli/wire.py | Update header generation constants |
