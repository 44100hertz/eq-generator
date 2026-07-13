/**
 * smart_volume.h — shared smart volume math (float)
 *
 * Both the desktop PipeWire filter (src/filter.c) and the ESP32
 * firmware (firmware/main/main.c) include this to guarantee identical
 * loudness-contour behaviour at every volume step.
 *
 * eq_coeffs.h must be included BEFORE this header.
 */
#pragma once

#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Compile-time constants from eq_coeffs.h ──────────────────────── */

#if !defined(EQGEN_QUIET_SHELF_DB)
#  error "smart_volume.h requires eq_coeffs.h to be included first"
#endif

#ifndef EQGEN_OVERBOOST_DB
#  define EQGEN_OVERBOOST_DB  0.0f
#endif

/* ── Results of a smart-volume interpolation ─────────────────────── */

typedef struct {
    float   pre_gain;
    float   boost;      /* (G-1) for loudness shelf */
    float   shelf_db;
    float   fft_pre_gain;
} SmartVolumeParams;

/* ── Compute shelf boost and pre-gain for a given volume level ────
 *
 * Maps 0–127 to pre-enhancer gain via vol_lut[vol] (linear-in-dB).
 * Additionally applies a loudness shelf (cube-root power law) to
 * compensate for Fletcher-Munson equal-loudness contours at low
 * volumes — exactly as the target speaker box does.
 * ──────────────────────────────────────────────────────────────── */

static inline void smart_volume_compute(uint8_t vol,
                                        float pg_loud,
                                        SmartVolumeParams *out)
{
    float t = (float)vol / 127.0f;   /* 0=quietest, 1=loudest */

    /* 1. Shelf boost: power-law (cube root) tracks equal-loudness contours.
     * The ear's bass sensitivity drops steeply in the first ~15 dB of
     * attenuation, then saturates.  A linear mapping would put only 50%
     * shelf at -30 dB attenuation — far too slow. */
    float atten_norm = 1.0f - t;   /* 0=loud, 1=quiet */
    out->shelf_db    = EQGEN_QUIET_SHELF_DB * powf(atten_norm, 0.33f);
    float shelf_linear = powf(10.0f, out->shelf_db / 20.0f);
    out->boost       = shelf_linear - 1.0f;

    /* 3. Pre-gain: reduce at quiet to offset shelf boost */
    float max_shelf_linear = powf(10.0f, EQGEN_QUIET_SHELF_DB / 20.0f);
    float pg_quiet         = pg_loud / max_shelf_linear;
    out->pre_gain          = pg_quiet + t * (pg_loud - pg_quiet);
    out->fft_pre_gain      = out->pre_gain;
}

/* ── Rebuild the 128-entry volume LUT ───────────────────────────── */

static inline void smart_volume_rebuild_lut(float vol_lut[128],
                                            float compensation_db,
                                            float speaker_level_db,
                                            float overboost_db)
{
    /* System-gain floor: vol=1 sits at -speaker_level dB.
     * speaker_level_db is roughly total system gain (speaker+amp)
     * in dB.  Typical: 40 = sensitive/high-gain, 60 = quiet.
     *
     * overboost_db: additional gain at vol=127.  Drives the enhancer
     * input above 0 dBFS → crossfade converts excess into harmonics. */
    float sv_db_floor = -speaker_level_db;
    if (sv_db_floor > -24.0f) sv_db_floor = -24.0f;   /* clamp minimum range */
    if (sv_db_floor < -80.0f) sv_db_floor = -80.0f;   /* sanity ceiling   */

    for (int i = 1; i < 128; i++) {
        float db = ((float)i / 127.0f) * (-sv_db_floor) + sv_db_floor
                   + compensation_db + overboost_db;
        if (db > overboost_db) db = overboost_db;
        vol_lut[i] = powf(10.0f, db / 20.0f);
    }
    vol_lut[0] = 0.0f;
}

#ifdef __cplusplus
}
#endif
