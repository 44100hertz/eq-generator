/**
 * enhancer.c — Harmonic bass enhancer implementation (float)
 *
 * Port of the Q16 fixed-point enhancer to full float precision.
 *
 * Uses:
 *   - Biquad for all filtering (float coefficients and state)
 *   - Env for peak-hold envelope following
 *   - Chebyshev T2(x) = 2x² - 1, T3(x) = 4x³ - 3x
 */

#include "enhancer.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef ESP_PLATFORM
#include "esp_cpu.h"
#include "esp_log.h"
#else
#include <stdio.h>
#define ESP_LOGI(tag, fmt, ...)  fprintf(stderr, "[" tag "] " fmt "\n", ##__VA_ARGS__)
static inline uint32_t esp_cpu_get_cycle_count(void) { return 0; }
#endif

#define SQRT2 1.4142135623730951f

void bass_design_butter_hp(float fc, float fs, float coeffs_out[5]) {
    float omega = tanf((float)M_PI * fc / fs);
    float c = 1.0f + SQRT2 * omega + omega * omega;

    coeffs_out[0] = 1.0f / c;
    coeffs_out[1] = -2.0f / c;
    coeffs_out[2] = 1.0f / c;
    coeffs_out[3] = (2.0f * (omega * omega - 1.0f)) / c;
    coeffs_out[4] = (1.0f - SQRT2 * omega + omega * omega) / c;
}

/* ── Configuration initialization ──────────────────────────────────── */

void BassEnhancerCfg_init(BassEnhancerCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float limiter_release_secs,
                          float pre_gain,
                          int eq_n_biquads, const float *eq_coeffs)
{
    memset(cfg, 0, sizeof(*cfg));
    cfg->cutoff_hz        = cutoff_hz;
    cfg->h2_amp           = h2_amp;
    cfg->h3_amp           = h3_amp;
    cfg->fs               = fs;
    cfg->release_secs     = release_secs;
    cfg->limiter_release_secs = limiter_release_secs;
    cfg->pre_gain         = pre_gain;

    /* Pre-compute release coefficients */
    cfg->release_coeff = expf(-1.0f / (fs * release_secs));
    cfg->limiter_release_coeff = expf(-1.0f / (fs * limiter_release_secs));

    /* Design LP filters (first-order) */
    cfg->lp_t2_alpha = 1.0f - expf(-2.0f * (float)M_PI * cutoff_hz / fs);
    cfg->lp_t3_alpha = 1.0f - expf(-2.0f * (float)M_PI * cutoff_hz * 0.5f / fs);
    bass_design_butter_hp(cutoff_hz, fs, cfg->hp_coeffs);
    bass_design_butter_hp(cutoff_hz, fs, cfg->hp_harm_coeffs);

    /* EQ settings */
    cfg->eq_n_biquads = eq_n_biquads;
    cfg->eq_coeffs    = eq_coeffs;

    /* Loudness defaults: disabled */
    cfg->loudness_alpha    = 0.0f;
    cfg->loudness_boost    = 0.0f;
    cfg->fundamental_bleed = 0.0f;
}

/* ── Loudness shelf setup ──────────────────────────────────────── */

void BassEnhancerCfg_set_loudness(BassEnhancerCfg *cfg,
                                  float fc, float fs, float boost_db)
{
    if (boost_db <= 0.0f || fc <= 0.0f) {
        cfg->loudness_alpha = 0.0f;
        cfg->loudness_boost = 0.0f;
        return;
    }
    cfg->loudness_alpha = 1.0f - expf(-2.0f * (float)M_PI * fc / fs);
    cfg->loudness_boost = powf(10.0f, boost_db / 20.0f) - 1.0f;
}

/* ── Runtime parameter update ──────────────────────────────────── */

#ifndef NAN
#define NAN (0.0f/0.0f)
#endif

void BassEnhancer_update_params(BassEnhancer *enh,
                                float h2_amp, float h3_amp,
                                float pre_gain,
                                float loudness_boost,
                                float fundamental_bleed)
{
    if (!isnan(h2_amp))            enh->cfg.h2_amp            = h2_amp;
    if (!isnan(h3_amp))            enh->cfg.h3_amp            = h3_amp;
    if (!isnan(pre_gain))          enh->cfg.pre_gain          = pre_gain;
    if (!isnan(loudness_boost))    enh->cfg.loudness_boost    = loudness_boost;
    if (!isnan(fundamental_bleed)) enh->cfg.fundamental_bleed = fundamental_bleed;
}

/* ── Enhancer init / reset ─────────────────────────────────────────── */

void BassEnhancer_init(BassEnhancer *enh,
                       const BassEnhancerCfg *cfg,
                       Biquad *eq_bqs_left,
                       Biquad *eq_bqs_right)
{
    memset(enh, 0, sizeof(*enh));
    enh->cfg = *cfg;

    /* Initialize EQ biquad arrays */
    enh->left.eq_bqs  = eq_bqs_left;
    enh->right.eq_bqs = eq_bqs_right;
    for (int i = 0; i < cfg->eq_n_biquads; i++) {
        biquad_init(&eq_bqs_left[i],  &cfg->eq_coeffs[i * 5]);
        biquad_init(&eq_bqs_right[i], &cfg->eq_coeffs[i * 5]);
    }

    /* Initialize LP filters */
    lp_init(&enh->left.lp_t2,   cfg->lp_t2_alpha);
    lp_init(&enh->left.lp_t3,   cfg->lp_t3_alpha);
    lp_init(&enh->right.lp_t2,  cfg->lp_t2_alpha);
    lp_init(&enh->right.lp_t3,  cfg->lp_t3_alpha);

    /* Initialize HP biquads */
    biquad_init(&enh->left.hp,      enh->cfg.hp_coeffs);
    biquad_init(&enh->left.hp_harm, enh->cfg.hp_harm_coeffs);
    biquad_init(&enh->right.hp,      enh->cfg.hp_coeffs);
    biquad_init(&enh->right.hp_harm, enh->cfg.hp_harm_coeffs);

    /* Initialize DC blockers (5 Hz cutoff) */
    float dc_R = expf(-2.0f * (float)M_PI * 5.0f / cfg->fs);
    dc_blocker_init(&enh->left.dc_block,  dc_R);
    dc_blocker_init(&enh->right.dc_block, dc_R);

    /* Initialize envelopes */
    env_init(&enh->left.env_t2,   cfg->release_coeff);
    env_init(&enh->left.env_t3,   cfg->release_coeff);
    env_init(&enh->right.env_t2,  cfg->release_coeff);
    env_init(&enh->right.env_t3,  cfg->release_coeff);

    env_init(&enh->left.env_lim,  cfg->limiter_release_coeff);
    env_init(&enh->right.env_lim, cfg->limiter_release_coeff);
}

void BassEnhancer_reset(BassEnhancer *enh) {
    BassEnhancerCfg *cfg = &enh->cfg;

    for (int i = 0; i < cfg->eq_n_biquads; i++) {
        biquad_reset(&enh->left.eq_bqs[i]);
        biquad_reset(&enh->right.eq_bqs[i]);
    }

    lp_reset(&enh->left.lp_t2);
    lp_reset(&enh->left.lp_t3);
    lp_reset(&enh->right.lp_t2);
    lp_reset(&enh->right.lp_t3);
    biquad_reset(&enh->left.hp);
    biquad_reset(&enh->left.hp_harm);
    biquad_reset(&enh->right.hp);
    biquad_reset(&enh->right.hp_harm);

    dc_blocker_reset(&enh->left.dc_block);
    dc_blocker_reset(&enh->right.dc_block);

    enh->left.env_t2.peak   = 0.0f;
    enh->left.env_t3.peak   = 0.0f;
    enh->right.env_t2.peak  = 0.0f;
    enh->right.env_t3.peak  = 0.0f;
    enh->left.env_lim.peak  = 0.0f;
    enh->right.env_lim.peak = 0.0f;

    enh->left.loudness_state  = 0.0f;
    enh->right.loudness_state = 0.0f;
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

static float enhancer_process_channel(BassEnhancerChan *ch,
                                      const BassEnhancerCfg *cfg,
                                      float x)
{
    uint32_t c0, c1, c2, c3;
    (void)c0; (void)c1; (void)c2; (void)c3;

    c0 = esp_cpu_get_cycle_count();

    /* ── Stage 0: DC blocker ───────────────────────────────────── */
    x = dc_blocker_tick(&ch->dc_block, x);

    /* ── Stage 0a: Pre-gain ─────────────────────────────────────── */
    if (cfg->pre_gain != 1.0f) {
        x *= cfg->pre_gain;
    }

    /* ── Stage 1: EQ preprocessing ──────────────────────────────── */
    float eq_out = biquad_cascade(ch->eq_bqs, cfg->eq_n_biquads, x);

    /* Bypass enhancer when cutoff ≤ 0 */
    if (cfg->cutoff_hz <= 0.0f) {
        c1 = esp_cpu_get_cycle_count();
        c3 = c1;
        enh_profile.frames++;
        enh_profile.cycles_total += (c3 - c0);
        enh_profile.cycles_eq     += (c1 - c0);
        return eq_out;
    }

    c1 = esp_cpu_get_cycle_count();

    /* ── Stage 2: First-order LP at cutoff_hz for T2 path ──────── */
    float lp_t2  = lp_tick(&ch->lp_t2, eq_out);
    float env_t2 = env_tick(&ch->env_t2, lp_t2);

    /* Floor envelope to 0.25 → 1/env capped at 4.0 (+12 dB) */
    float env_t2_norm = env_t2 > 0.25f ? env_t2 : 0.25f;
    float norm_t2;
    if (env_t2 > 1e-6f) {
        norm_t2 = lp_t2 / env_t2_norm;
    } else {
        norm_t2 = 0.0f;
    }

    /* Scale by h2_amp and apply Chebyshev T2 */
    float harm_scaled_t2;
    if (cfg->h2_amp <= 0.0f) {
        harm_scaled_t2 = 0.0f;
    } else {
        float cheb_in_t2 = norm_t2 * cfg->h2_amp;
        float harm_t2 = cheb_t2(cheb_in_t2);
        harm_scaled_t2 = harm_t2 * env_t2;
    }

    /* ── Stage 3: First-order LP at cutoff_hz/2 for T3 path ────── */
    float lp_t3  = lp_tick(&ch->lp_t3, eq_out);
    float env_t3 = env_tick(&ch->env_t3, lp_t3);

    float env_t3_norm = env_t3 > 0.25f ? env_t3 : 0.25f;
    float norm_t3;
    if (env_t3 > 1e-6f) {
        norm_t3 = lp_t3 / env_t3_norm;
    } else {
        norm_t3 = 0.0f;
    }

    float harm_scaled_t3;
    if (cfg->h3_amp <= 0.0f) {
        harm_scaled_t3 = 0.0f;
    } else {
        float cheb_in_t3 = norm_t3 * cfg->h3_amp;
        float harm_t3 = cheb_t3(cheb_in_t3);
        harm_scaled_t3 = harm_t3 * env_t3;
    }

    /* ── Stage 3a: Fundamental bleed ────────────────────────────── */
    float fundamental = 0.0f;
    if (cfg->fundamental_bleed != 0.0f) {
        fundamental = lp_t2 * cfg->fundamental_bleed;
    }

    c2 = esp_cpu_get_cycle_count();

    /* ── Stage 4: Mix ───────────────────────────────────────────── */
    float harm_sum = harm_scaled_t2 + harm_scaled_t3 + fundamental;
    float harm_hp  = biquad_tick(&ch->hp_harm, harm_sum);
    float dry_hp   = biquad_tick(&ch->hp, eq_out);
    float out = dry_hp + harm_hp;

    /* ── Harmonic AGC limiter ──────────────────────────────────── */
    float env_lim  = env_tick(&ch->env_lim, out);
    float env_peak = env_lim > 1.0f ? env_lim : 1.0f;
    float lim_gain = 1.0f / env_peak;
    harm_hp *= lim_gain;
    out = dry_hp + harm_hp;

    /* ── Stage 5: Loudness shelf (one-pole low shelf, after limiter) ── */
    if (cfg->loudness_boost != 0.0f) {
        float diff = out - ch->loudness_state;
        ch->loudness_state += cfg->loudness_alpha * diff;
        out += ch->loudness_state * cfg->loudness_boost;
    }

    c3 = esp_cpu_get_cycle_count();

    /* ── Accumulate profile ─────────────────────────────────────── */
    enh_profile.frames++;
    enh_profile.cycles_total += (c3 - c0);
    enh_profile.cycles_eq     += (c1 - c0);
    enh_profile.cycles_env    += (c2 - c1);
    enh_profile.cycles_mix    += (c3 - c2);
    enh_profile.cycles_harm   += (c2 - c1) / 2;

    return out;
}

/* ── Stereo processing ─────────────────────────────────────────────── */

void BassEnhancer_process_stereo(BassEnhancer *enh,
                                 float *left, float *right)
{
    *left  = enhancer_process_channel(&enh->left,  &enh->cfg, *left);
    *right = enhancer_process_channel(&enh->right, &enh->cfg, *right);

    /* Safety net: clamp NaN/Inf to zero.  Once NaN enters biquad state
     * it propagates permanently — catching it here resets the output
     * but NOT the biquad state.  If NaN appears, dsp_init should also
     * reset biquads on next stream start. */
    if (!isfinite(*left))  *left  = 0.0f;
    if (!isfinite(*right)) *right = 0.0f;
}
