# Smart Volume

Loudness-compensated volume control driven by AVRCP volume changes. Fully implemented in firmware and mirrored in the desktop Python tools.

## What it does

At lower volumes the ear is less sensitive to bass (Fletcher-Munson / ISO 226 equal-loudness contours). Smart volume applies two corrections on every volume change:

1. **Loudness shelf** — a 2nd-order low shelf (fc=80 Hz, Q=0.332) boosts bass at quiet volumes. The shelf gain scales linearly with the dB drop from peak: ~2.5 dB per 10 dB of attenuation (FM_SLOPE). At vol=127 the shelf is 0 dB (flat).
2. **Pre-gain reduction** — offsets the shelf boost so the signal doesn't clip the enhancer input. Linear interpolation between `pg_quiet` and `pg_loud` as volume changes.

Both corrections are recomputed on every AVRCP volume event (0–127), keeping the IIR biquad frequencies fixed — only gain/Q change.

## Where the code lives

| Layer | File | Role |
|---|---|---|
| C shared header | `src/smart_volume.h` | Shelf boost + pre-gain math (`smart_volume_compute`, `smart_volume_rebuild_lut`) |
| C shared header | `src/volume_control.h` | Orchestration: `volume_set()` chains compute → update DSP → rebuild LUT |
| C DSP | `src/enhancer.c` / `enhancer.h` | `dsp_pipe_update_params()` applies pre_gain and loudness_boost at runtime |
| Firmware | `firmware/main/main.c` | Calls `volume_set()` on every AVRCP volume change; uses `vol_lut[128]` in the audio loop |
| Python analysis | `eqgen/smart_volume.py` | Desktop equivalent: models the effective correction curve at any volume |
| Python sweep | `eqgen/sweep.py` | `build_vol_lut()` / `build_full_vol_lut()` mirror the C LUT for end-to-end testing |

## Algorithm

```
vol ∈ [0, 127], t = vol/127

drop_db  = speaker_level × (1 − t)           # dB below peak
shelf_db = FM_SLOPE × drop_db                 # loudness shelf boost

max_shelf = 10^(FM_SLOPE × speaker_level / 20)
pg_quiet  = pg_loud / max_shelf
pre_gain  = pg_quiet + t × (pg_loud − pg_quiet)
```

The volume LUT maps vol → pre-enhancer gain (dB-linear) so the DSP receives an already-attenuated signal. The enhancer's headroom budget and crossfade operate relative to this attenuated level, naturally filling quieter levels with more harmonics.

The `speaker_level_db` parameter (typically 40–80) sets the total system gain range. Clamped to 24–80 dB. `overboost_db` adds extra gain at vol=127 to deliberately drive the enhancer into harmonic generation.
