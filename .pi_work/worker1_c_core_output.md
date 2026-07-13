# Worker 1: C Core Headroom Enhancer Refactor

## Summary

Replaced the reactive AGC harmonic limiter with a 5ms lookahead headroom budget system. Volume scaling stays before the enhancer (not reverted).

## Changed files

| File | Change |
|------|--------|
| `src/enhancer.h` | New LOOKAHEAD_LEN constant (240). Removed limiter fields, lp_t3_alpha, lp_t3/env_t3/env_lim from channel. Added push_gain, lp_alpha, hp_lookahead, hp_ring, eq_delay, delay_pos. Updated function signatures. |
| `src/enhancer.c` | Rewrote enhancer_process_channel with delay line + HP ring + peak detector + headroom budget (bleed priority 1, harmonics priority 2). Simplified to single LP+env path. Removed limiter init/reset. Updated profiling label. |
| `src/enhancer_api.h` | Replaced limiter_release_secs with push_gain in EnhancerParams. |
| `src/enhancer_api.c` | Updated BassEnhancerCfg_init call site. |
| `src/smart_volume.h` | SmartVolumeParams already stripped of h2/h3/bleed — only shelf + pre_gain remain. |
| `src/eq_coeffs.h` | EQGEN_LIMITER_RELEASE_SECS → EQGEN_PUSH_GAIN. Removed crossfade LO_T. Updated comments. |
| `src/filter.c` | Uses EQGEN_PUSH_GAIN, update_params already matches new 2-arg signature. Volume-before-enhancer preserved. |
| `src/ladspa/ladspa_wrapper.c` | Both BassEnhancerCfg_init calls updated to EQGEN_PUSH_GAIN. |
| `src/test_enhancer.c` | Both calls updated to push_gain=1.0. |
| `firmware/main/main.c` | Eqgen_push_gain, update_params already matches. Volume-before-enhancer preserved. |
| `eqgen/enhancer_ffi.py` | EnhancerParams struct: limiter_release_secs → push_gain. create_enhancer signature updated. |

## Signal flow

```
eq_out ──┬──► hp_lookahead ──► hp_ring[240] ──► upcoming_peak
         │
         └──► eq_delay[240]
                   │
                   ▼   (delayed by 240 samples)
              ┌────┴────────────────────┐
              │  LP(fc) ──► env        │
              │  norm = lp/env         │
              │  bleed = clip(lp*bleed_cap, room)   │
              │  harm  = clip(Chebyshev(norm), room - |bleed|)
              │  out = dry_hp + bleed + harm │
              └────────────────────────┘
                   │
            room = push_gain * max(0, 1.0 - upcoming_peak)
```

## Build & test

```
$ make all     # clean compile of enhancer.so + ladspa + filter
$ ./test_enhancer
  All 5 tests PASS
  Peak output: 0.0784 (well within 0dBFS)
```

## Remaining

- Python pipeline (`eqgen/pipeline.py` step 8, `model.py`, `wire.py`) still needs updating to compute h2_amp/h3_amp/bleed_cap from measurement data and output new header format.
- Firmware ESP-IDF build not tested (needs ESP32 toolchain).
