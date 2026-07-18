/**
 * volume_control.h — shared live volume-control orchestration
 *
 * Wraps the 4-step update (smart_volume_compute → dsp_pipe_update_params
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
 * 1. Compute shelf boost & pre-gain for this volume
 * 2. Push boost/pre-gain into the enhancer
 * 3. Rebuild the LUT — no compensation needed since the shelf
 *    runs pre-enhancer (part of the target), not post-enhancer.
 *
 * Call this every time the volume index changes.  vol_lut[vol] is
 * then the effective pre-enhancer gain factor for the audio loop.
 * ──────────────────────────────────────────────────────────────── */

static inline SmartVolumeParams volume_set(uint8_t vol,
                                            float vol_lut[128],
                                            DspPipe *enh)
{
    float pg_loud = 1.0f;
#ifdef EQGEN_PRE_GAIN
    pg_loud = EQGEN_PRE_GAIN;
#endif

    SmartVolumeParams svp;
    smart_volume_compute(vol, pg_loud, &svp);

    /* push_gain ramps from floor (compile-time) to 1.0 at max volume.
     * At vol=0: push_gain = floor (max power saving).
     * At vol=127: push_gain = 1.0 (full fundamental, show off the speaker). */
    float t = (float)vol / 127.0f;
    float push_gain_v = EQGEN_PUSH_GAIN + (1.0f - EQGEN_PUSH_GAIN) * t;

    dsp_pipe_update_params(enh, svp.pre_gain, svp.boost, push_gain_v);

    smart_volume_rebuild_lut(vol_lut, 0.0f,
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
