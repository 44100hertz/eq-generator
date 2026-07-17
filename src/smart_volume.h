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

#if !defined(EQGEN_SPEAKER_LEVEL_DB)
#  error "smart_volume.h requires eq_coeffs.h (EQGEN_SPEAKER_LEVEL_DB) to be included first"
#endif

#ifndef EQGEN_OVERBOOST_DB
#  define EQGEN_OVERBOOST_DB  0.0f
#endif

/* Shelf slope: dB of bass boost per dB of SPL drop at 1 kHz.
 * Fitted to ISO 226:2023 equal-loudness contours (differential
 * compensation from 70 phon).  Defined in eq_coeffs.h (generated).
 * At vol=127 (drop=0) the shelf is always 0 dB. */
#ifndef FM_SLOPE
#  error "eq_coeffs.h must be included before smart_volume.h (defines FM_SLOPE)"
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
 * Additionally applies a loudness shelf (linear FM contour: 0.8 dB
 * per 10 dB of attenuation) to compensate for Fletcher-Munson
 * equal-loudness contours at low volumes.
 *
 * The shelf tracks the absolute listening level: the drop from peak
 * is speaker_level × (1 − vol/127) dB.  A speaker with 80 dB range
 * drops 80 dB to vol=0 → 6.4 dB shelf; a 40 dB range speaker drops
 * 40 dB → 3.2 dB shelf.  Peak volume (vol=127) always gets 0 shelf.
 * ──────────────────────────────────────────────────────────────── */

static inline void smart_volume_compute(uint8_t vol,
                                        float pg_loud,
                                        SmartVolumeParams *out)
{
    float t = (float)vol / 127.0f;   /* 0=quietest, 1=loudest */

    /* 1. Shelf boost: linear with dB drop from peak volume.
     *    0.8 dB / 10 dB is an adequate linear approximation of the
     *    FM equal-loudness contours over our gain range (24–80 dB). */
    float drop_db  = (1.0f - t) * (float)EQGEN_SPEAKER_LEVEL_DB;
    out->shelf_db  = FM_SLOPE * drop_db;
    float shelf_linear = powf(10.0f, out->shelf_db / 20.0f);
    out->boost     = shelf_linear - 1.0f;

    /* 2. Pre-gain: constant.  The loudness shelf is applied
     *    pre-enhancer (part of the psychoacoustic target), so
     *    the enhancer always operates at full signal level. */
    out->pre_gain     = pg_loud;
    out->fft_pre_gain = pg_loud;
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
