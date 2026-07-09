/**
 * enhancer.c — Harmonic bass enhancer implementation (process_fixed_v2)
 *
 * Port of dsp.py process_fixed_v2 to ESP32 Q16 fixed-point.
 *
 * Uses:
 *   - BiquadQ28 for all filtering (Q4.28 coefficients, Q16 state)
 *   - Env for peak-hold envelope following
 *   - FR_DIV for 1/envelope normalization
 *   - Chebyshev T2(x) = 2x² - 1, T3(x) = 4x³ - 3x
 *
 * Intermediate precision: 64-bit accumulator for poly eval, 32-bit clip.
 */

#include "enhancer.h"
#include <math.h>
#include <string.h>

#ifdef ESP_PLATFORM
#include "esp_cpu.h"
#include "esp_log.h"
#else
#include <stdio.h>
#define ESP_LOGI(tag, fmt, ...)  fprintf(stderr, "[" tag "] " fmt "\n", ##__VA_ARGS__)
static inline uint32_t esp_cpu_get_cycle_count(void) { return 0; }
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Q4.28 scale: 2^28 = 268435456 */
#define Q28 (268435456.0)
#define Q16 (65536.0)

static int32_t float2q28(float v) {
    if (v >= 7.999999f) return 0x7FFFFFFF;
    if (v <= -8.0f)     return (int32_t)0x80000000;
    return (int32_t)(v * Q28 + (v >= 0.0f ? 0.5f : -0.5f));
}

void bass_design_butter_hp_q28(float fc, float fs, int32_t coeffs_out[5]) {
    float omega = (float)tan((double)M_PI * fc / fs);
    float c = 1.0f + 1.41421356237f * omega + omega * omega;

    float b0 = 1.0f / c;
    float b1 = -2.0f / c;
    float b2 = 1.0f / c;
    float a1 = (2.0f * (omega * omega - 1.0f)) / c;
    float a2 = (1.0f - 1.41421356237f * omega + omega * omega) / c;

    coeffs_out[0] = float2q28(b0);
    coeffs_out[1] = float2q28(b1);
    coeffs_out[2] = float2q28(b2);
    coeffs_out[3] = float2q28(a1);
    coeffs_out[4] = float2q28(a2);
}

/* ── Envelope reciprocal LUT ───────────────────────────────────────── */

void ReciprocalLUT_init(ReciprocalLUT *lut) {
    for (int i = 0; i < ENV_LUT_SIZE; i++) {
        /* lut[i] = 1.0 / ((i + 1) / 256.0) = 256.0 / (i + 1) */
        /* Store as Q16.16: value * 65536 * 65536 / something? */
        /* Actually: for an input x with top 8 bits = i, we want 1/x in Q16.16 */
        /* x is in Q16, so x = (i * 2^N) for some normalization shift */
        /* We store reciprocal of normalized index for later use */
        float recip = 1.0f / ((float)(i + 1) / 256.0f);
        lut->entries[i] = (uint32_t)(recip * 65536.0f + 0.5f);
    }
}

/* ── Configuration initialization ──────────────────────────────────── */

void BassEnhancerCfg_init(BassEnhancerCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float limiter_release_secs,
                          float pre_gain,
                          int eq_n_biquads, const int32_t *eq_coeffs)
{
    memset(cfg, 0, sizeof(*cfg));
    cfg->cutoff_hz = cutoff_hz;
    cfg->h2_amp     = h2_amp;
    cfg->h3_amp     = h3_amp;
    cfg->fs         = fs;
    cfg->release_secs = release_secs;
    cfg->limiter_release_secs = limiter_release_secs;

    /* Pre-compute Q16 values */
    cfg->release_coeff_q16 = (int32_t)(exp(-1.0 / (fs * release_secs)) * Q16 + 0.5);
    cfg->h2_amp_q16 = (int32_t)(h2_amp * Q16 + 0.5);
    cfg->h3_amp_q16 = (int32_t)(h3_amp * Q16 + 0.5);

    /* Pre-gain clamp: [1/16, 16.0] = [-24, +24] dB range in Q16 */
    if (pre_gain < 0.0625f) pre_gain = 0.0625f;
    if (pre_gain > 16.0f)   pre_gain = 16.0f;
    cfg->pre_gain_q16 = (int32_t)(pre_gain * Q16 + (pre_gain >= 0.0f ? 0.5f : -0.5f));

    /* Limiter release coefficient */
    cfg->limiter_release_coeff_q16 = (int32_t)(
        exp(-1.0 / (fs * limiter_release_secs)) * Q16 + 0.5);

    /* Design LP/HP filters */
    /* LP filters use first-order Q16 (lp_q16.h) to avoid Q28 biquad
       truncation at the tiny fc/fs ratios used for bass extraction. */
    cfg->lp_t2_alpha_q16 = (int32_t)((1.0f - expf(-2.0f * (float)M_PI * cutoff_hz / fs)) * Q16 + 0.5f);
    cfg->lp_t3_alpha_q16 = (int32_t)((1.0f - expf(-2.0f * (float)M_PI * cutoff_hz * 0.5f / fs)) * Q16 + 0.5f);
    bass_design_butter_hp_q28(cutoff_hz, fs, cfg->hp_coeffs);
    bass_design_butter_hp_q28(cutoff_hz, fs, cfg->hp_harm_coeffs);

    /* EQ settings */
    cfg->eq_n_biquads = eq_n_biquads;
    cfg->eq_coeffs    = eq_coeffs;
}

/* ── Enhancer init / reset ─────────────────────────────────────────── */

void BassEnhancer_init(BassEnhancer *enh,
                       const BassEnhancerCfg *cfg,
                       const ReciprocalLUT *lut,
                       BiquadQ28 *eq_bqs_left,
                       BiquadQ28 *eq_bqs_right)
{
    memset(enh, 0, sizeof(*enh));
    enh->cfg = *cfg;

    /* Initialize EQ biquad arrays */
    enh->left.eq_bqs  = eq_bqs_left;
    enh->right.eq_bqs = eq_bqs_right;
    for (int i = 0; i < cfg->eq_n_biquads; i++) {
        BiquadQ28_init(&eq_bqs_left[i],  &cfg->eq_coeffs[i * 5]);
        BiquadQ28_init(&eq_bqs_right[i], &cfg->eq_coeffs[i * 5]);
    }

    /* Initialize LP filters (first-order Q16) */
    LPQ16_init(&enh->left.lp_t2,   enh->cfg.lp_t2_alpha_q16);
    LPQ16_init(&enh->left.lp_t3,   enh->cfg.lp_t3_alpha_q16);
    LPQ16_init(&enh->right.lp_t2,  enh->cfg.lp_t2_alpha_q16);
    LPQ16_init(&enh->right.lp_t3,  enh->cfg.lp_t3_alpha_q16);

    /* Initialize HP biquads */
    BiquadQ28_init(&enh->left.hp,      enh->cfg.hp_coeffs);
    BiquadQ28_init(&enh->left.hp_harm, enh->cfg.hp_harm_coeffs);
    BiquadQ28_init(&enh->right.hp,      enh->cfg.hp_coeffs);
    BiquadQ28_init(&enh->right.hp_harm, enh->cfg.hp_harm_coeffs);

    /* Initialize DC blockers (5 Hz cutoff) */
    float dc_R = (float)exp(-2.0 * M_PI * 5.0 / cfg->fs);
    int32_t dc_R_q16 = (int32_t)(dc_R * Q16 + 0.5f);
    DCBlocker_init(&enh->left.dc_block,  dc_R_q16);
    DCBlocker_init(&enh->right.dc_block, dc_R_q16);

    /* Initialize envelopes */
    Env_init(&enh->left.env_t2,  lut, cfg->release_coeff_q16);
    Env_init(&enh->left.env_t3,  lut, cfg->release_coeff_q16);
    Env_init(&enh->right.env_t2, lut, cfg->release_coeff_q16);
    Env_init(&enh->right.env_t3, lut, cfg->release_coeff_q16);

    Env_init(&enh->left.env_lim,  lut, cfg->limiter_release_coeff_q16);
    Env_init(&enh->right.env_lim, lut, cfg->limiter_release_coeff_q16);
}

void BassEnhancer_reset(BassEnhancer *enh) {
    BassEnhancerCfg *cfg = &enh->cfg;

    for (int i = 0; i < cfg->eq_n_biquads; i++) {
        BiquadQ28_reset(&enh->left.eq_bqs[i]);
        BiquadQ28_reset(&enh->right.eq_bqs[i]);
    }

    LPQ16_reset(&enh->left.lp_t2);
    LPQ16_reset(&enh->left.lp_t3);
    LPQ16_reset(&enh->right.lp_t2);
    LPQ16_reset(&enh->right.lp_t3);
    BiquadQ28_reset(&enh->left.hp);
    BiquadQ28_reset(&enh->left.hp_harm);
    BiquadQ28_reset(&enh->right.hp);
    BiquadQ28_reset(&enh->right.hp_harm);

    DCBlocker_reset(&enh->left.dc_block);
    DCBlocker_reset(&enh->right.dc_block);

    /* Only reset envelope peak — keep lut and release coeff set by init.
       memset(..., 0, sizeof(Env)) would wipe the lut pointer → NULL deref. */
    enh->left.env_t2.peak   = 0;
    enh->left.env_t3.peak   = 0;
    enh->right.env_t2.peak  = 0;
    enh->right.env_t3.peak  = 0;
    enh->left.env_lim.peak  = 0;
    enh->right.env_lim.peak = 0;
}

/* ── Single-channel processing ─────────────────────────────────────── */

/** Chebyshev T2: 2x² - 1. x in Q16 (±1 = 65536). Result in Q16. */
int32_t cheb_t2(int32_t x) {
    /* x_sq = x² in Q32 */
    int64_t x_sq = (int64_t)x * (int64_t)x;
    /* 2x²: shift Q32 -> Q16 (>> 16), then *2 */
    int64_t two_x2 = (x_sq >> 15);  /* = 2 * x² / 2^16 */

    /* T2 = 2x² - 1 = two_x2 - 65536 */
    return (int32_t)(two_x2 - 65536);
}

/** Chebyshev T3: 4x³ - 3x. x in Q16. Result in Q16. */
int32_t cheb_t3(int32_t x) {
    /* x² in Q32 */
    int64_t x_sq = (int64_t)x * (int64_t)x;
    /* x³ in Q48 */
    int64_t x_cu = x_sq * (int64_t)x;
    /* 4x³: Q48 -> Q16: >> 32 then *4, or >> 30 */
    int64_t four_x3 = (x_cu >> 30);  /* 4 * x³ / 2^48 */
    /* 3x in Q16 */
    int64_t three_x = (int64_t)x * 3;
    return (int32_t)(four_x3 - three_x);
}

/* ── Profiling ─────────────────────────────────────────────────────── */

EnhancerProfile enh_profile;

void enhancer_profile_report(void) {
    if (enh_profile.frames == 0) return;
    uint32_t n = enh_profile.frames;
    ESP_LOGI("enhancer", "=== %lu frames ===", (unsigned long)n);
    ESP_LOGI("enhancer", "  total: %llu cy/frame  (%llu cy total)",
             (unsigned long long)(enh_profile.cycles_total / n),
             (unsigned long long)enh_profile.cycles_total);
    ESP_LOGI("enhancer", "  EQ cascade: %llu cy/frame",
             (unsigned long long)(enh_profile.cycles_eq / n));
    ESP_LOGI("enhancer", "  env+normalize: %llu cy/frame",
             (unsigned long long)(enh_profile.cycles_env / n));
    ESP_LOGI("enhancer", "  Chebyshev+scale: %llu cy/frame",
             (unsigned long long)(enh_profile.cycles_harm / n));
    ESP_LOGI("enhancer", "  HP+mix+limiter: %llu cy/frame",
             (unsigned long long)(enh_profile.cycles_mix / n));
    ESP_LOGI("enhancer", "  I2S write: %llu cy/frame",
             (unsigned long long)(enh_profile.cycles_i2s / n));
    memset(&enh_profile, 0, sizeof(enh_profile));
}

/* ── Per-channel tick ──────────────────────────────────────────────── */

static int32_t enhancer_process_channel(BassEnhancerChan *ch,
                                        const BassEnhancerCfg *cfg,
                                        int32_t x)
{
    const ReciprocalLUT *lut = ch->env_t2.lut;
    uint32_t c0, c1, c2, c3;

    c0 = esp_cpu_get_cycle_count();

    /* ── Stage 0: DC blocker (first-order HP at ~5 Hz) ───────────── */
    x = DCBlocker_tick(&ch->dc_block, x);

    /* ── Stage 0a: Pre-gain (uniform gain before EQ, Q16) ────────── */
    if (cfg->pre_gain_q16 != 65536) {
        x = (int32_t)(((int64_t)x * (int64_t)cfg->pre_gain_q16) >> 16);
    }

    /* ── Stage 1: EQ preprocessing ────────────────────────────────── */
    int32_t eq_out = BiquadQ28_cascade(ch->eq_bqs, cfg->eq_n_biquads, x);

    c1 = esp_cpu_get_cycle_count();

    /* ── Stage 1: First-order LP at cutoff_hz for T2 path ──────── */
    int32_t lp_t2 = LPQ16_tick(&ch->lp_t2, eq_out);
    int32_t env_t2 = Env_tick(&ch->env_t2, lp_t2);

    int32_t norm_t2;
    if (env_t2 > 6) {
        uint32_t inv = ReciprocalLUT_lookup(lut, env_t2);
        norm_t2 = (int32_t)(((int64_t)lp_t2 * (int64_t)inv) >> 16);
    } else {
        norm_t2 = 0;
    }

    /* Scale by h2_amp and apply Chebyshev T2 */
    int32_t harm_scaled_t2;
    if (cfg->h2_amp_q16 == 0) {
        harm_scaled_t2 = 0;
    } else {
        int32_t cheb_in_t2 = (int32_t)(((int64_t)norm_t2 * (int64_t)cfg->h2_amp_q16) >> 16);
        int32_t harm_t2 = cheb_t2(cheb_in_t2);
        harm_scaled_t2 = (int32_t)(((int64_t)harm_t2 * (int64_t)env_t2) >> 16);
    }

    /* ── Stage 2: First-order LP at cutoff_hz/2 for T3 path ─────── */
    int32_t lp_t3 = LPQ16_tick(&ch->lp_t3, eq_out);
    int32_t env_t3 = Env_tick(&ch->env_t3, lp_t3);

    int32_t norm_t3;
    if (env_t3 > 6) {
        uint32_t inv = ReciprocalLUT_lookup(lut, env_t3);
        norm_t3 = (int32_t)(((int64_t)lp_t3 * (int64_t)inv) >> 16);
    } else {
        norm_t3 = 0;
    }

    int32_t harm_scaled_t3;
    if (cfg->h3_amp_q16 == 0) {
        harm_scaled_t3 = 0;
    } else {
        int32_t cheb_in_t3 = (int32_t)(((int64_t)norm_t3 * (int64_t)cfg->h3_amp_q16) >> 16);
        int32_t harm_t3 = cheb_t3(cheb_in_t3);
        harm_scaled_t3 = (int32_t)(((int64_t)harm_t3 * (int64_t)env_t3) >> 16);
    }

    c2 = esp_cpu_get_cycle_count();

    /* ── Stage 3: Mix ──────────────────────────────────────────────── */
    int32_t harm_sum = harm_scaled_t2 + harm_scaled_t3;
    int32_t harm_hp  = BiquadQ28_tick_q44(&ch->hp_harm, harm_sum);
    int32_t dry_hp   = BiquadQ28_tick_q44(&ch->hp, eq_out);
    int32_t out = dry_hp + harm_hp;

    /* ── Harmonic AGC limiter ──────────────────────────────────── */
    int32_t env_lim  = Env_tick(&ch->env_lim, out);
    int32_t env_peak = env_lim > 65536 ? env_lim : 65536;
    uint32_t inv = ReciprocalLUT_lookup(lut, env_peak);
    int32_t lim_gain = (int32_t)(inv >> 16);
    harm_hp = (int32_t)(((int64_t)harm_hp * (int64_t)lim_gain) >> 16);
    out = dry_hp + harm_hp;

    /* Brick-wall clamp (±1.0 in Q16) */
    if (out >  65536) out =  65536;
    if (out < -65536) out = -65536;

    c3 = esp_cpu_get_cycle_count();

    /* ── Accumulate profile ─────────────────────────────────────── */
    enh_profile.frames++;
    enh_profile.cycles_total += (c3 - c0);
    enh_profile.cycles_eq     += (c1 - c0);
    enh_profile.cycles_env    += (c2 - c1);
    enh_profile.cycles_mix    += (c3 - c2);
    /* cycles_harm is a subset of env, we'll approximate it */
    enh_profile.cycles_harm   += (c2 - c1) / 2;  /* rough: half of env block */

    return out;
}

/* ── Stereo processing ─────────────────────────────────────────────── */

void BassEnhancer_process_stereo(BassEnhancer *enh,
                                 int32_t *left, int32_t *right)
{
    *left  = enhancer_process_channel(&enh->left,  &enh->cfg, *left);
    *right = enhancer_process_channel(&enh->right, &enh->cfg, *right);
}
