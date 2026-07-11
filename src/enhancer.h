/**
 * enhancer.h — Harmonic bass enhancer (float)
 *
 * Pipeline per channel:
 *   1. DC blocker (5 Hz HP)
 *   2. Pre-gain
 *   3. EQ preprocessing (cascaded biquads)
 *   4. LP(fc) → envelope → normalize → T2 → scale     (T2 harmonics)
 *   5. LP(fc/2) → envelope → normalize → T3 → scale   (T3 harmonics)
 *   6. HP(fc) dry + HP(fc) harmonics → mix
 *   7. Fundamental bleed (optional)
 *   8. Harmonic AGC limiter
 *   9. Loudness shelf (optional)
 *
 * Usage:
 *   BassEnhancer enh;
 *   BassEnhancer_init(&enh, &cfg, eq_bqs_left, eq_bqs_right);
 *   // Per stereo frame:
 *   float left = ...;
 *   float right = ...;
 *   BassEnhancer_process_stereo(&enh, &left, &right);
 */

#pragma once
#include <stdint.h>
#include "biquad.h"
#include "lp.h"
#include "envelope.h"
#include "dc_blocker.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Configuration ─────────────────────────────────────────────────── */
typedef struct {
    float   cutoff_hz;      /* bass crossover frequency (typ. 60 Hz)   */
    float   h2_amp;         /* 2nd harmonic amplitude (0..1, typ. 0.33)*/
    float   h3_amp;         /* 3rd harmonic amplitude (0..1, typ. 0.33)*/
    float   fs;             /* sample rate (typ. 44100)                */
    float   release_secs;   /* envelope release time (typ. 0.2)        */
    float   limiter_release_secs; /* limiter release time (typ. 0.049)   */

    /* Pre-computed coefficients */
    float   release_coeff;      /* exp(-1/(fs*release))                */
    float   limiter_release_coeff; /* exp(-1/(fs*lim_rel))             */
    float   pre_gain;           /* gain before EQ                      */

    /* Loudness compensation */
    float   loudness_alpha;      /* one-pole LP alpha (0 = disabled)   */
    float   loudness_boost;      /* (G-1), 0 = no shelf                */
    float   fundamental_bleed;   /* LP→mix bleed (0=none, ~0.25=max)   */

    /* Butterworth HP coefficients */
    float   hp_coeffs[5];       /* HP at cutoff_hz  (dry path)         */
    float   hp_harm_coeffs[5];  /* HP at cutoff_hz  (harmonics path)   */

    /* First-order LP coefficients */
    float   lp_t2_alpha;        /* LP at cutoff_hz  (for T2 path)      */
    float   lp_t3_alpha;        /* LP at cutoff_hz/2 (for T3 path)     */

    /* EQ settings */
    int     eq_n_biquads;       /* number of EQ biquads                */
    const float *eq_coeffs;     /* pointer to eq_coeffs array           */
} BassEnhancerCfg;

/* ── Per-channel state ─────────────────────────────────────────────── */
typedef struct {
    DCBlocker   dc_block;
    float       loudness_state;

    Biquad     *eq_bqs;          /* array of eq_n_biquads Biquad       */

    /* Enhancer filters */
    LP          lp_t2;            /* LP filter for T2 path              */
    LP          lp_t3;            /* LP filter for T3 path              */
    Biquad      hp;               /* HP filter for dry path             */
    Biquad      hp_harm;          /* HP filter for harmonics path       */

    /* Envelopes */
    Env         env_t2;           /* envelope for T2 path               */
    Env         env_t3;           /* envelope for T3 path               */
    Env         env_lim;          /* envelope for harmonic limiter      */
} BassEnhancerChan;

/* ── Stereo enhancer ───────────────────────────────────────────────── */
typedef struct {
    BassEnhancerCfg cfg;
    BassEnhancerChan left;
    BassEnhancerChan right;
} BassEnhancer;

/* ── Initialization ────────────────────────────────────────────────── */

/** Compute Butterworth HP coefficients.
 *  coeffs_out[5] = {b0, b1, b2, a1, a2}.
 */
void bass_design_butter_hp(float fc, float fs, float coeffs_out[5]);

/** Configure the loudness compensation shelf.
 *  fc: corner frequency (typ. 200 Hz)
 *  boost_db: max gain at DC (typ. 8 dB).  0 disables the shelf.
 */
void BassEnhancerCfg_set_loudness(BassEnhancerCfg *cfg,
                                  float fc, float fs, float boost_db);

/** Update runtime parameters without resetting filter state.
 *  Pass NaN for any parameter to leave it unchanged. */
void BassEnhancer_update_params(BassEnhancer *enh,
                                float h2_amp, float h3_amp,
                                float pre_gain,
                                float loudness_boost,
                                float fundamental_bleed);

/** Initialize BassEnhancerCfg from user-friendly parameters.
 *  Designs all LP/HP filters and pre-computes coefficients. */
void BassEnhancerCfg_init(BassEnhancerCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float limiter_release_secs,
                          float pre_gain,
                          int eq_n_biquads, const float *eq_coeffs);

/** Initialize a BassEnhancer instance. */
void BassEnhancer_init(BassEnhancer *enh,
                       const BassEnhancerCfg *cfg,
                       Biquad *eq_bqs_left,
                       Biquad *eq_bqs_right);

/** Reset all filter states. */
void BassEnhancer_reset(BassEnhancer *enh);

/* ── Processing ────────────────────────────────────────────────────── */

/** Process one stereo sample pair in-place.
 *  Input: float samples in range [-1.0, 1.0]
 *  Output: float samples in range [-1.0, 1.0] (written in-place)
 */
void BassEnhancer_process_stereo(BassEnhancer *enh,
                                 float *left, float *right);

/* ── Profiling counters ────────────────────────────────────────────── */

typedef struct {
    uint32_t frames;
    uint64_t cycles_total;
    uint64_t cycles_eq;
    uint64_t cycles_env;
    uint64_t cycles_harm;
    uint64_t cycles_mix;
    uint64_t cycles_i2s;
} EnhancerProfile;

extern EnhancerProfile enh_profile;

void enhancer_profile_report(void);

/* ── Low-level Chebyshev ──────────────────────────────────────────── */

/** Chebyshev T2: 2x² - 1. */
static inline float cheb_t2(float x) {
    return 2.0f * x * x - 1.0f;
}

/** Chebyshev T3: 4x³ - 3x. */
static inline float cheb_t3(float x) {
    return 4.0f * x * x * x - 3.0f * x;
}

#ifdef __cplusplus
}
#endif
