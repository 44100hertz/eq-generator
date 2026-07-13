/**
 * volume_control.h — shared live volume-control orchestration
 *
 * Wraps the 4-step update (smart_volume_compute → BassEnhancer_update_params
 * → smart_volume_rebuild_lut → vol_lut lookup) so that both the desktop
 * PipeWire filter (src/filter.c) and the ESP32 firmware
 * (firmware/main/main.c) use identical code.
 *
 * Any edit to the box's volume behaviour must be made here — it
 * automatically translates to the desktop.
 *
 * Include AFTER eq_coeffs.h, enhancer.h, smart_volume.h.
 */
#pragma once

#include <stdint.h>
#include "smart_volume.h"
#include "enhancer.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Initialise the volume LUT with zero shelf compensation ──────── */

static inline void volume_init_lut(float vol_lut[128])
{
    smart_volume_rebuild_lut(vol_lut, 0.0f,
                             EQGEN_SPEAKER_LEVEL_DB,
                             EQGEN_OVERBOOST_DB);
}

/* ── Handle a volume change (0–127) ───────────────────────────────
 *
 * 1. Compute shelf boost & adjusted pre-gain for this volume
 * 2. Push boost/pre-gain into the enhancer
 * 3. Rebuild the LUT with shelf_db as compensation (so pre-enhancer
 *    gain already includes the loudness contour)
 *
 * Call this every time the volume index changes.  vol_lut[vol] is
 * then the effective pre-enhancer gain factor for the audio loop.
 * ──────────────────────────────────────────────────────────────── */

static inline SmartVolumeParams volume_set(uint8_t vol,
                                            float vol_lut[128],
                                            BassEnhancer *enh)
{
    float pg_loud = 1.0f;
#ifdef EQGEN_PRE_GAIN
    pg_loud = EQGEN_PRE_GAIN;
#endif

    SmartVolumeParams svp;
    smart_volume_compute(vol, pg_loud, &svp);

    BassEnhancer_update_params(enh, svp.pre_gain, svp.boost);

    smart_volume_rebuild_lut(vol_lut, svp.shelf_db,
                             EQGEN_SPEAKER_LEVEL_DB,
                             EQGEN_OVERBOOST_DB);
    return svp;
}

/* ── Get the pre-enhancer gain factor for the audio loop ────────── */

static inline float volume_gain(float vol_lut[128], uint8_t vol)
{
    return vol_lut[vol];
}

#ifdef __cplusplus
}
#endif
