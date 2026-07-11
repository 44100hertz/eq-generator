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

#if !defined(EQGEN_QUIET_SHELF_DB) || !defined(EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T)
#  error "smart_volume.h requires eq_coeffs.h to be included first"
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

    /* 1. Harmonic → bleed crossfade.
     *
     * At low volumes (t ≤ LO_T, vol ≤ 25): harmonics fully cut,
     * replaced by fundamental bleed.  The bleed signal comes from
     * the same EQ-preprocessed LP(fc) as T2 harmonics, so the EQ
     * correction is still applied — we just skip Chebyshev shaping.
     * At -45 dBFS, EQ errors are inaudible anyway.
     *
     * At high volumes (t ≥ HI_T, vol ≥ 63): full preset, no bleed.
     * Bit-identical to the old linear-interpolation behaviour.
     *
     * Between: linear crossfade. */
    if (t >= EQGEN_HARMONIC_BLEED_CROSSFADE_HI_T) {
        out->h2_amp = EQGEN_H2_AMP;
        out->h3_amp = EQGEN_H3_AMP;
        out->bleed  = 0.0f;
    } else if (t <= EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T) {
        out->h2_amp = 0.0f;
        out->h3_amp = 0.0f;
        out->bleed  = EQGEN_QUIET_FUNDAMENTAL_BLEED;
    } else {
        float xfade = (t - EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T)
                    / (EQGEN_HARMONIC_BLEED_CROSSFADE_HI_T
                       - EQGEN_HARMONIC_BLEED_CROSSFADE_LO_T);
        out->h2_amp = xfade * EQGEN_H2_AMP;
        out->h3_amp = xfade * EQGEN_H3_AMP;
        out->bleed  = (1.0f - xfade) * EQGEN_QUIET_FUNDAMENTAL_BLEED;
    }

    /* 2. Shelf boost: power-law (cube root) tracks equal-loudness contours.
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
                                            float speaker_level_db)
{
    /* System-gain floor: vol=1 sits at -speaker_level dB.
     * speaker_level_db is roughly total system gain (speaker+amp)
     * in dB.  Typical: 40 = sensitive/high-gain, 60 = quiet. */
    float sv_db_floor = -speaker_level_db;
    if (sv_db_floor > -24.0f) sv_db_floor = -24.0f;   /* clamp minimum range */
    if (sv_db_floor < -80.0f) sv_db_floor = -80.0f;   /* sanity ceiling   */

    for (int i = 1; i < 128; i++) {
        float db = ((float)i / 127.0f) * (-sv_db_floor) + sv_db_floor
                   + compensation_db;
        if (db > 0.0f) db = 0.0f;
        vol_lut[i] = powf(10.0f, db / 20.0f);
    }
    vol_lut[0] = 0.0f;
}

#ifdef __cplusplus
}
#endif
