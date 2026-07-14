/**
 * enhancer.h — Harmonic bass enhancer (float)
 *
 * Pipeline per channel:
 *   1. DC blocker (5 Hz HP)
 *   2. Pre-gain
 *   3. EQ preprocessing (cascaded biquads)
 *   4. Lookahead: LR4 LP/HP split → ring buffer → upcoming peak
 *   5. Delay line (5 ms) → LR4 LP → envelope → normalize
 *   6. Crossfade: target = env, room = push_gain·(1−upcoming)
 *      w = (target−room)/(target·(1−h)), blend fundamental ↔ harmonics
 *   7. Mix: LR4 HP + (1−w)·lp_fund + harm
 *   8. Loudness shelf (optional)
 *   9. Soft output clamp (tanh)
 *
 * Volume scaling MUST run before this enhancer so the headroom
 * budget is relative to the already-attenuated signal.
 *
 * The LP/HP split uses an LR4 crossover (two cascaded 2nd-order
 * Butterworth sections each).  Unlike a complementary 1st-order split,
 * the LR4 has negligible phase bleed below cutoff — the lookahead
 * detects actual HF content, not crossover artifacts.
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
#include "envelope.h"
#include "dc_blocker.h"

#define LOOKAHEAD_LEN         240   /* 5 ms @ 48 kHz */
#define LR4_SECTIONS          2     /* cascaded 2nd-order sections */

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
    float   push_gain;       /* headroom fill strength (1.0 = 0 dBFS) */

    /* Pre-computed coefficients */
    float   release_coeff;      /* exp(-1/(fs*release))                */
    float   env_smooth_alpha;   /* smoothed envelope tracking rate     */
    float   lim_release_coeff;  /* limiter release (49ms) one-pole coeff */
    float   pre_gain;           /* gain before EQ                      */
    float   h2_scale;           /* h2_amp / (h2_amp + h3_amp)        */
    float   h3_scale;           /* h3_amp / (h2_amp + h3_amp)        */

    /* Loudness compensation */
    float   loudness_alpha;      /* one-pole LP alpha (0 = disabled)   */
    float   loudness_boost;      /* (G-1), 0 = no shelf                */

    /* 2nd-order Butterworth HP (harmonics cleanup path only) */
    float   hp_harm_coeffs[5];

    /* LR4 crossover: 2 cascaded 2nd-order Butterworth sections each */
    float   lp_cross_coeffs[LR4_SECTIONS * 5];  /* LP biquad coeffs */
    float   hp_cross_coeffs[LR4_SECTIONS * 5];  /* HP biquad coeffs */

    /* EQ settings */
    int     eq_n_biquads;       /* number of EQ biquads                */
    const float *eq_coeffs;     /* pointer to eq_coeffs array           */
} BassEnhancerCfg;

/* ── Per-channel state ─────────────────────────────────────────────── */
typedef struct {
    DCBlocker   dc_block;
    float       loudness_state;

    Biquad     *eq_bqs;          /* array of eq_n_biquads Biquad       */

    /* LR4 crossover filter states */
    Biquad      lp_cross[LR4_SECTIONS];  /* LP for fundamental extract */
    Biquad      hp_cross[LR4_SECTIONS];  /* HP for dry high-pass      */
    Biquad      hp_lookahead[LR4_SECTIONS]; /* HP for lookahead       */
    Biquad      hp_harm;              /* HP for harmonics cleanup       */

    /* Envelopes */
    Env         env;              /* instantaneous envelope             */
    float       env_smooth;       /* env low-passed to kill ripple     */

    /* Crossfade slew state */
    float       w_slew;           /* low-passed crossfade weight       */
    float       harm_amp_slew;    /* low-passed harmonic amplitude     */

    /* Bass-sum limiter */
    float       bass_gain_red;    /* gain reduction factor, 0..1       */

    /* Lookahead */
    float       hp_ring[LOOKAHEAD_LEN];  /* |dry_hp| window         */
    int         hp_ring_pos;
    float       eq_delay[LOOKAHEAD_LEN]; /* delay line for eq_out   */
    int         delay_pos;
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

/** Design LR4 crossover coefficients.
 *  lp_coeffs[10]: two cascaded 2nd-order LP biquads (5 coeffs each)
 *  hp_coeffs[10]: two cascaded 2nd-order HP biquads (5 coeffs each)
 */
void bass_design_lr4(float fc, float fs,
                      float lp_coeffs[LR4_SECTIONS * 5],
                      float hp_coeffs[LR4_SECTIONS * 5]);

/** Configure the loudness compensation shelf.
 *  fc: corner frequency (typ. 200 Hz)
 *  boost_db: max gain at DC (typ. 8 dB).  0 disables the shelf.
 */
void BassEnhancerCfg_set_loudness(BassEnhancerCfg *cfg,
                                  float fc, float fs, float boost_db);

/** Update runtime parameters without resetting filter state.
 *  h2_amp, h3_amp are set once via cfg and
 *  never changed at runtime.  Pass NaN to leave unchanged. */
void BassEnhancer_update_params(BassEnhancer *enh,
                                float pre_gain,
                                float loudness_boost);

/** Initialize BassEnhancerCfg from user-friendly parameters.
 *  Designs all LP/HP filters and pre-computes coefficients. */
void BassEnhancerCfg_init(BassEnhancerCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float push_gain,
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
