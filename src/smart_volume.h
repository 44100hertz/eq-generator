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

#if !defined(EQGEN_QUIET_H2_AMP) || !defined(EQGEN_QUIET_SHELF_DB)
#  error "smart_volume.h requires eq_coeffs.h to be included first"
#endif

#ifndef SV_DB_FLOOR
#  define SV_DB_FLOOR  -60.0f   /* dB gain at vol=1 (vol=0 is mute) */
#endif

/* ── Results of a smart-volume interpolation ─────────────────────── */

typedef struct {
    float   h2_amp;
    float   h3_amp;
    float   pre_gain;
    float   boost;      /* (G-1) for loudness shelf */
    float   bleed;      /* fundamental bleed amount */
    float   shelf_db;
    float   fft_pre_gain;
} SmartVolumeParams;

/* ── Interpolate between quiet and loud presets ──────────────────── */

static inline void smart_volume_compute(uint8_t vol,
                                        float pg_loud,
                                        SmartVolumeParams *out)
{
    float t = (float)vol / 127.0f;   /* 0=quietest, 1=loudest */

    /* 1. Harmonic amplitudes: half → full */
    out->h2_amp = EQGEN_QUIET_H2_AMP + t * (EQGEN_H2_AMP - EQGEN_QUIET_H2_AMP);
    out->h3_amp = EQGEN_QUIET_H3_AMP + t * (EQGEN_H3_AMP - EQGEN_QUIET_H3_AMP);

    /* 2. Shelf boost: linearly from full to 0 dB */
    out->shelf_db    = EQGEN_QUIET_SHELF_DB * (1.0f - t);
    float shelf_linear = powf(10.0f, out->shelf_db / 20.0f);
    out->boost       = shelf_linear - 1.0f;

    /* 3. Pre-gain: reduce at quiet to offset shelf boost */
    float max_shelf_linear = powf(10.0f, EQGEN_QUIET_SHELF_DB / 20.0f);
    float pg_quiet         = pg_loud / max_shelf_linear;
    out->pre_gain          = pg_quiet + t * (pg_loud - pg_quiet);
    out->fft_pre_gain      = out->pre_gain;

    /* 4. Fundamental bleed: proportional to harmonic reduction */
    out->bleed = EQGEN_QUIET_FUNDAMENTAL_BLEED * (1.0f - t);
}

/* ── Rebuild the 128-entry volume LUT ───────────────────────────── */

static inline void smart_volume_rebuild_lut(float vol_lut[128],
                                            float compensation_db)
{
    for (int i = 1; i < 128; i++) {
        float db = ((float)i / 127.0f) * (-SV_DB_FLOOR) + SV_DB_FLOOR
                   + compensation_db;
        if (db > 0.0f) db = 0.0f;
        vol_lut[i] = powf(10.0f, db / 20.0f);
    }
    vol_lut[0] = 0.0f;
}

#ifdef __cplusplus
}
#endif
