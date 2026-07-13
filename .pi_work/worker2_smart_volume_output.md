# Worker 2: smart_volume.h refactor

## Changes made

### src/smart_volume.h
- Removed `h2_amp`, `h3_amp`, `bleed` fields from `SmartVolumeParams` struct
- Removed the harmonic↔bleed crossfade block (entire t-based if/else chain using EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T/HI_T)
- Simplified `#error` guard to only check `EQGEN_QUIET_SHELF_DB`
- Updated doc comment to explain h2/h3/bleed are now compile-time constants
- Preserved: shelf boost computation, pre-gain computation, volume LUT rebuild, fft_pre_gain

### src/enhancer.c
- Updated `BassEnhancer_update_params` signature: removed `h2_amp`, `h3_amp`, `fundamental_bleed` parameters
- Now only takes `(enh, pre_gain, loudness_boost)`

### src/enhancer.h
- Already had the new signature from worker 1 (verified)

### src/filter.c
- Updated `BassEnhancer_update_params` call to pass only `(svp.pre_gain, svp.boost)`
- Updated fprintf log line to remove h2/h3/bleed references
- Volume-before-enhancer preserved (not reverted)

### firmware/main/main.c
- Updated `BassEnhancer_update_params` call to pass only `(svp.pre_gain, svp.boost)`
- Updated ESP_LOGI log line to remove h2/h3/bleed references
- Volume-before-enhancer preserved (not reverted)

## Validation
- `grep -rn "svp\.h2_amp\|svp\.h3_amp\|svp\.bleed" src/ firmware/` returns no results
- All callers updated to match new `BassEnhancer_update_params` signature
- No volume-before-enhancer revert in filter.c or main.c

## Residual risks
- `BassEnhancerCfg_init` in filter.c and main.c still pass `EQGEN_LIMITER_RELEASE_SECS` — this is worker 1's domain to update the init API. Code currently compiles since the constant still exists in eq_coeffs.h and the parameter is still accepted.

```acceptance-report
{
  "criteriaSatisfied": [
    {
      "id": "criterion-1",
      "status": "satisfied",
      "evidence": "smart_volume.h: removed crossfade + h2/h3/bleed from struct; enhancer.c: removed h2/h3/bleed from update_params; filter.c + main.c: updated callers, volume-before-enhancer preserved"
    }
  ],
  "changedFiles": [
    "src/smart_volume.h",
    "src/enhancer.c",
    "src/filter.c",
    "firmware/main/main.c"
  ],
  "testsAddedOrUpdated": [],
  "commandsRun": [
    {
      "command": "grep -rn 'svp\\.h2_amp\\|svp\\.h3_amp\\|svp\\.bleed' src/ firmware/",
      "result": "passed",
      "summary": "no residual references to removed struct fields"
    }
  ],
  "validationOutput": [
    "smart_volume.h: SmartVolumeParams has 4 fields (pre_gain, boost, shelf_db, fft_pre_gain) — h2/h3/bleed removed",
    "smart_volume_compute: only shelf + pre-gain computation remains",
    "enhancer.c BassEnhancer_update_params: 3-param signature (h2/h3/bleed removed)",
    "filter.c: update_params call passes 2 args; log line sanitized",
    "main.c: update_params call passes 2 args; log line sanitized"
  ],
  "residualRisks": [
    "BassEnhancerCfg_init still accepts limiter_release_secs — worker 1 will handle"
  ],
  "noStagedFiles": true,
  "diffSummary": "Stripped harmonic↔bleed volume crossfade from smart_volume; h2/h3/bleed are now compile-time constants set at enhancer init; BassEnhancer_update_params simplified to (enh, pre_gain, boost); both filter.c and main.c callers updated; volume-before-enhancer preserved.",
  "reviewFindings": [],
  "manualNotes": "enhancer.h already had the new update_params signature — worker 1 had already updated it. enhancer.c was out of sync and has been fixed."
}
```
