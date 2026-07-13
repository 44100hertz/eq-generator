# Worker 3: Python-side headroom enhancer changes

## Changes summary

### 1. `eqgen/model.py` — Added `compute_harmonic_efficacy()`

Computes h2_amp, h3_amp, bleed_cap from the EQ correction curve and speaker measurement data.

Formula:
```
h2_amp = Σ ear(2f) × correction_ratio(2f) / Σ ear(f) × correction_ratio(f) × 0.35
```

Where `correction_ratio(f) = 10^(correction_db(f)/20)` and `0.35` = human hearing rolloff (~ -9 dB).

- Flat correction: h2 ≈ 1.81, h3 = 2.0 (clipped) — harmonics are very audible relative to the fundamental
- +30 dB bass correction: h2 ≈ 0.23 — the EQ handles the bass, harmonics barely needed
- fc=0: all zeros (no enhancer active)

### 2. `eqgen/pipeline.py` — Replaced step 8

**Was:** model_gain_needed() loop modifying corr[] to account for expected harmonic output, then holding bass flat at fc/2.

**Now:** compute_harmonic_efficacy() computes h2/h3/bleed without modifying the EQ curve. The EQ just flattens the raw speaker response.

**Return value changed:** `(freqs, gains_db, sample_rate, max_gain_db)` → `(freqs, gains_db, sample_rate, max_gain_db, efficacy)` where efficacy = `{"h2_amp": ..., "h3_amp": ..., "bleed_cap": ...}`.

**compute_smart_volume_curves()**: Removed `h2`, `h3` parameters — now only computes shelf + pre_gain curves (harmonic crossfade is gone).

### 3. `eqgen/cli/wire.py` — Updated header generation

- `run_full_pipeline()`: unpacks new 5th return value (efficacy)
- cfg dict now includes `h2_amp`, `h3_amp`, `bleed_cap`, `push_gain` from efficacy
- cfg dict no longer has `limiter_release_secs`
- `generate_eq_header()` outputs: `EQGEN_H2_AMP`, `EQGEN_H3_AMP`, `EQGEN_BLEED_CAP`, `EQGEN_PUSH_GAIN`
- Removed: `EQGEN_LIMITER_RELEASE_SECS`, `EQGEN_QUIET_FUNDAMENTAL_BLEED`, `EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T`, `EQGEN_HARMONIC_BLEED_CROSSFADE_HI_T`

### 4. All other callers — Mechanical 4→5 tuple unpack

| File | Change |
|------|--------|
| `cli/audition.py` | `_efficacy` added |
| `cli/preset.py` | `_efficacy` added |
| `cli/eqgen.py` | `_efficacy` added |
| `cli/graph_check.py` | No change (uses detailed=True dict) |
| `server.py` | Updated compute_smart_volume_curves call (removed h2/h3) |
| `tests/test_sine_sweep.py` | `_efficacy` added |
| `tests/test_real_fullband_fit.py` | `_efficacy` added |
| `cli/export.py` | No change (doesn't call run_pipeline) |

## Validation

- All 9 Python files parse with no syntax errors
- `compute_harmonic_efficacy()` verified with manual test cases:
  - Flat 0dB correction → h2≈1.81, h3=2.0, bleed=0.5
  - +30dB bass correction → h2≈0.23, h3≈0.22, bleed≈0.22
  - fc=0 → all zeros

## Open risks

1. **h2/h3 in presets are still used** as fallback when efficacy dict returns zeros (fc=0 case) — but in `run_full_pipeline`, h2/h3 from the preset are only passed to `run_pipeline()` which uses them via the `h2`/`h3` params that are now unused in the new step 8. The computed efficacy values are used for the header. This is a minor API dead spot — the h2/h3 pipeline params could be removed in a follow-up.
2. **server.py** still loads preset.h2/preset.h3 but no longer passes them to `compute_smart_volume_curves`. The preset class still has those fields — fine, they're just unused.
3. **process_track()** in pipeline.py still uses hardcoded h2/h3 defaults — the enhancer_ffi call signature will need updating when the C side changes.

## Next steps

- Worker 1 (C core) needs to implement the delay line + ring buffer + headroom budget
- Worker 2 (smart_volume.h) needs to strip the crossfade
- After both, the `enhancer_ffi.py` wrapper needs API alignment
