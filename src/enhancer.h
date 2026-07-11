/**
 * enhancer.h — Harmonic bass enhancer (process_fixed_v2) in Q16
 *
 * Pipeline per channel:
 *   1. EQ preprocessing (cascaded BiquadQ28 from eq_coeffs.h)
 *   2. LP(fc) → envelope → T2 → scale     (LPQ16, no Q28 truncation)
 *   3. LP(fc/2) → envelope → T3 → scale   (LPQ16, no Q28 truncation)
 *   4. HP(fc) dry + HP(fc) harmonics → mix
 *
 * LP filters use first-order Q16 (lp_q16.h) to avoid Q28 biquad
 * truncation at the tiny fc/fs ratios used for bass extraction.
 * HP coefficients are Butterworth, designed offline at Q4.28 precision.
 *
 * Usage:
 *   BassEnhancer enh;
 *   BassEnhancer_init(&enh, &lut, cfg);
 *   // Per stereo frame:
 *   int32_t left = ...;   // Q16 input sample
 *   int32_t right = ...;
 *   BassEnhancer_process_stereo(&enh, &left, &right);
 */

#pragma once
#include <stdint.h>
#include "biquad_q28.h"
#include "lp_q16.h"
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

    /* ── below: pre-computed Q16/Q28 values ──────────────────────── */
    int32_t release_coeff_q16;  /* exp(-1/(fs*release)) in Q16        */
    int32_t h2_amp_q16;         /* h2_amp in Q16                      */
    int32_t h3_amp_q16;         /* h3_amp in Q16                      */
    int32_t limiter_release_coeff_q16; /* exp(-1/(fs*lim_rel)) in Q16  */
    int32_t pre_gain_q16;       /* pre-gain before EQ (Q16, 1.0=65536)*/

    /* ── Loudness compensation ──────────────────────────────────── */
    int32_t loudness_alpha_q16;     /* one-pole LP alpha (0 = disabled)    */
    int32_t loudness_boost_q16;     /* (G-1) in Q16, 0 = no shelf         */
    int32_t fundamental_bleed_q16;  /* LP→mix bleed (0=none, ~16400=max)  */

    /* Butterworth HP coefficients in Q4.28 */
    int32_t hp_coeffs[5];       /* HP at cutoff_hz  (dry path)        */
    int32_t hp_harm_coeffs[5];  /* HP at cutoff_hz  (harmonics path)  */

    /* First-order LP coefficients in Q16 (no Q28 truncation) */
    int32_t lp_t2_alpha_q16;   /* LP at cutoff_hz  (for T2 path)     */
    int32_t lp_t3_alpha_q16;   /* LP at cutoff_hz/2 (for T3 path)    */

    /* ── EQ settings ─────────────────────────────────────────────── */
    int     eq_n_biquads;       /* number of EQ biquads               */
    const int32_t *eq_coeffs;   /* pointer to eq_coeffs_q28 array     */
} BassEnhancerCfg;

/* ── Per-channel state ─────────────────────────────────────────────── */
typedef struct {
    /* DC blocker (first stage, blocks subsonic DC) */
    DCBlocker dc_block;

    /* Loudness shelf state (one-pole LP, applied after limiter) */
    int32_t loudness_state;

    /* EQ biquads */
    BiquadQ28 *eq_bqs;          /* array of eq_n_biquads BiquadQ28    */

    /* Enhancer filters */
    LPQ16      lp_t2;            /* LP filter for T2 path              */
    LPQ16      lp_t3;            /* LP filter for T3 path              */
    BiquadQ28  hp;               /* HP filter for dry path             */
    BiquadQ28  hp_harm;          /* HP filter for harmonics path       */

    /* Envelopes */
    Env env_t2;                 /* envelope for T2 path               */
    Env env_t3;                 /* envelope for T3 path               */
    Env env_lim;                /* envelope for harmonic limiter      */
} BassEnhancerChan;

/* ── Stereo enhancer ───────────────────────────────────────────────── */
typedef struct {
    BassEnhancerCfg cfg;
    BassEnhancerChan left;
    BassEnhancerChan right;
} BassEnhancer;

/* ── Initialization ────────────────────────────────────────────────── */

/** Compute Butterworth HP coefficients in Q4.28.
 *  Designed offline; this is a runtime convenience using tan().
 *  coeffs_out[5] = {b0, b1, b2, a1, a2} in Q4.28.
 */
void bass_design_butter_hp_q28(float fc, float fs, int32_t coeffs_out[5]);

/** Configure the loudness compensation shelf.
 *  fc: corner frequency (typ. 200 Hz)
 *  boost_db: max gain at DC (typ. 8 dB).  0 disables the shelf.
 *
 *  Call after BassEnhancerCfg_init(). Safe to call at runtime;
 *  only updates cfg fields, does not touch filter state. */
void BassEnhancerCfg_set_loudness(BassEnhancerCfg *cfg,
                                  float fc, float fs, float boost_db);

/** Update runtime parameters without resetting filter state.
 *  All in Q16 format. Pass -1 for any parameter to leave it unchanged.
 *  fundamental_bleed controls how much of the LP (un-HP'd) signal
 *  is mixed directly into the output to compensate for reduced harmonics. */
void BassEnhancer_update_params(BassEnhancer *enh,
                                int32_t h2_q16, int32_t h3_q16,
                                int32_t pre_gain_q16,
                                int32_t loudness_boost_q16,
                                int32_t fundamental_bleed_q16);

/** Initialize BassEnhancerCfg from user-friendly float parameters.
 *  Designs all LP/HP filters and pre-computes Q16/Q28 values.
 */
void BassEnhancerCfg_init(BassEnhancerCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float limiter_release_secs,
                          float pre_gain,
                          int eq_n_biquads, const int32_t *eq_coeffs);

/** Initialize a BassEnhancer instance.
 *  Call BassEnhancerCfg_init() first, then allocate and pass BiquadQ28
 *  arrays for eq_bqs (needs eq_n_biquads entries per channel).
 */
void BassEnhancer_init(BassEnhancer *enh,
                       const BassEnhancerCfg *cfg,
                       const ReciprocalLUT *lut,
                       BiquadQ28 *eq_bqs_left,
                       BiquadQ28 *eq_bqs_right);

/** Reset all filter states (call on parameter change or stop). */
void BassEnhancer_reset(BassEnhancer *enh);

/* ── Processing ────────────────────────────────────────────────────── */

/** Process one stereo sample pair in-place.
 *  Each sample is int32_t in Q16 format (range approx [-32767, 32767]).
 *  left and right are modified in-place.
 */
void BassEnhancer_process_stereo(BassEnhancer *enh,
                                 int32_t *left, int32_t *right);

/* ── Profiling counters (compile-time zero-cost when PROFILE_CYCLES unset) ── */

typedef struct {
    uint32_t frames;        /* total frames processed */
    uint64_t cycles_total;  /* total cycles in process_stereo */
    uint64_t cycles_eq;     /* cycles in biquad cascade */
    uint64_t cycles_env;    /* cycles in envelope + normalize */
    uint64_t cycles_harm;   /* cycles in Chebyshev + scale */
    uint64_t cycles_mix;    /* cycles in HP + mix + limiter */
    uint64_t cycles_i2s;    /* cycles in I2S write (measured in caller) */
} EnhancerProfile;

extern EnhancerProfile enh_profile;

void enhancer_profile_report(void);

/* ── Low-level Chebyshev (exposed for testing) ────────────────────── */
int32_t cheb_t2(int32_t x);
int32_t cheb_t3(int32_t x);

#ifdef __cplusplus
}
#endif
