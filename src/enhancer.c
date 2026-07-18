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

/* ── LR4 crossover design ─────────────────────────────────────────── */

void bass_design_lr4(float fc, float fs,
                      float lp_coeffs[LR4_SECTIONS * 5],
                      float hp_coeffs[LR4_SECTIONS * 5])
{
    float k = tanf((float)M_PI * fc / fs);
    float k2 = k * k;
    float norm = 1.0f / (1.0f + SQRT2 * k + k2);

    /* Single 2nd-order Butterworth LP section */
    float lp_b0 = k2 * norm;
    float lp_b1 = 2.0f * k2 * norm;
    float lp_b2 = k2 * norm;
    float lp_a1 = 2.0f * (k2 - 1.0f) * norm;
    float lp_a2 = (1.0f - SQRT2 * k + k2) * norm;

    /* Single 2nd-order Butterworth HP section */
    float hp_b0 = norm;
    float hp_b1 = -2.0f * norm;
    float hp_b2 = norm;
    float hp_a1 = lp_a1;  /* same denominator as LP */
    float hp_a2 = lp_a2;

    /* Two cascaded sections for LR4 */
    for (int i = 0; i < LR4_SECTIONS; i++) {
        lp_coeffs[i * 5 + 0] = lp_b0;
        lp_coeffs[i * 5 + 1] = lp_b1;
        lp_coeffs[i * 5 + 2] = lp_b2;
        lp_coeffs[i * 5 + 3] = lp_a1;
        lp_coeffs[i * 5 + 4] = lp_a2;

        hp_coeffs[i * 5 + 0] = hp_b0;
        hp_coeffs[i * 5 + 1] = hp_b1;
        hp_coeffs[i * 5 + 2] = hp_b2;
        hp_coeffs[i * 5 + 3] = hp_a1;
        hp_coeffs[i * 5 + 4] = hp_a2;
    }
}


/* ── Configuration initialization ──────────────────────────────────── */

void dsp_pipe_cfg_init(DspPipeCfg *cfg,
                          float cutoff_hz, float h2_amp, float h3_amp,
                          float release_secs, float fs,
                          float push_gain,
                          float pre_gain,
                          int eq_n_biquads, const float *eq_coeffs)
{
    memset(cfg, 0, sizeof(*cfg));
    cfg->cutoff_hz        = cutoff_hz;
    cfg->h2_amp           = h2_amp;
    cfg->h3_amp           = h3_amp;
    cfg->fs               = fs;
    cfg->release_secs     = release_secs;
    cfg->push_gain        = push_gain;
    cfg->pre_gain         = pre_gain;

    /* Pre-compute release coefficient */
    cfg->release_coeff = expf(-1.0f / (fs * release_secs));
    cfg->env_smooth_alpha = 2.0f * (float)M_PI * cutoff_hz / (3.0f * fs);

    /* Pre-compute limiter release (49 ms) */
    cfg->lim_release_coeff = 1.0f - expf(-1.0f / (fs * 0.049f));

    /* Pre-compute full-band limiter release (3.0 s) */
    cfg->fb_release_coeff = 1.0f - expf(-1.0f / (fs * 3.0f));
    /* Pre-compute Chebyshev input scales */
    float h_sum = h2_amp + h3_amp;
    cfg->h2_scale = (h_sum > 1e-6f) ? h2_amp / h_sum : 0.0f;
    cfg->h3_scale = (h_sum > 1e-6f) ? h3_amp / h_sum : 0.0f;

    /* Design LR4 crossover LP/HP coefficients */
    bass_design_lr4(cutoff_hz, fs, cfg->lp_cross_coeffs, cfg->hp_cross_coeffs);

    /* Design HP for harmonics cleanup (2nd-order Butterworth) */
    bass_design_butter_hp(cutoff_hz, fs, cfg->hp_harm_coeffs);

    /* EQ settings */
    cfg->eq_n_biquads = eq_n_biquads;
    cfg->eq_coeffs    = eq_coeffs;

    /* Loudness defaults: flat (unity) biquad */
    cfg->loudness_fc = 80.0f;
    cfg->loudness_Q  = 0.332f;
    cfg->loudness_coeffs[0] = 1.0f;  /* b0 */
    cfg->loudness_coeffs[1] = 0.0f;  /* b1 */
    cfg->loudness_coeffs[2] = 0.0f;  /* b2 */
    cfg->loudness_coeffs[3] = 0.0f;  /* a1 */
    cfg->loudness_coeffs[4] = 0.0f;  /* a2 */
}

/* ── 2nd-order low-shelf biquad design (RBJ Audio EQ Cookbook) ── */

static void design_loudness_shelf(float fc, float Q, float boost_db,
                                   float fs, float coeffs[5])
{
    if (boost_db <= 0.001f || fc <= 0.0f) {
        /* Flat: unity gain at all frequencies */
        coeffs[0] = 1.0f;  /* b0 */
        coeffs[1] = 0.0f;  /* b1 */
        coeffs[2] = 0.0f;  /* b2 */
        coeffs[3] = 0.0f;  /* a1 */
        coeffs[4] = 0.0f;  /* a2 */
        return;
    }
    float A    = powf(10.0f, boost_db / 40.0f);
    float w0   = 2.0f * (float)M_PI * fc / fs;
    float cosw = cosf(w0);
    float sinw = sinf(w0);
    float alpha = sinw / (2.0f * Q);
    float sqrtA = sqrtf(A);

    float b0 = A * ((A + 1.0f) - (A - 1.0f) * cosw + 2.0f * sqrtA * alpha);
    float b1 = 2.0f * A * ((A - 1.0f) - (A + 1.0f) * cosw);
    float b2 = A * ((A + 1.0f) - (A - 1.0f) * cosw - 2.0f * sqrtA * alpha);
    float a0 = (A + 1.0f) + (A - 1.0f) * cosw + 2.0f * sqrtA * alpha;
    float a1 = -2.0f * ((A - 1.0f) + (A + 1.0f) * cosw);
    float a2 = (A + 1.0f) + (A - 1.0f) * cosw - 2.0f * sqrtA * alpha;

    /* Normalise to a0=1 */
    float inv_a0 = 1.0f / a0;
    coeffs[0] = b0 * inv_a0;
    coeffs[1] = b1 * inv_a0;
    coeffs[2] = b2 * inv_a0;
    coeffs[3] = a1 * inv_a0;
    coeffs[4] = a2 * inv_a0;
}

/* ── Loudness shelf setup ──────────────────────────────────────── */

void dsp_pipe_cfg_set_loudness(DspPipeCfg *cfg,
                                  float fc, float Q, float fs,
                                  float boost_db)
{
    cfg->loudness_fc = fc;
    cfg->loudness_Q  = Q;
    design_loudness_shelf(fc, Q, boost_db, fs, cfg->loudness_coeffs);
}

/* ── Runtime parameter update ──────────────────────────────────── */

#ifndef NAN
#define NAN (0.0f/0.0f)
#endif

void dsp_pipe_update_params(DspPipe *enh,
                                float pre_gain,
                                float loudness_boost,
                                float push_gain)
{
    if (!isnan(pre_gain)) {
        enh->cfg.pre_gain = pre_gain;
    }
    if (!isnan(push_gain)) {
        enh->cfg.push_gain = push_gain;
    }
    if (!isnan(loudness_boost)) {
        /* Convert (G-1) back to dB and recompute biquad coefficients */
        float G = loudness_boost + 1.0f;
        float boost_db = (G > 0.001f) ? 20.0f * log10f(G) : 0.0f;
        design_loudness_shelf(enh->cfg.loudness_fc,
                               enh->cfg.loudness_Q,
                               boost_db,
                               enh->cfg.fs,
                               enh->cfg.loudness_coeffs);
    }
}

/* ── Enhancer init / reset ─────────────────────────────────────────── */

void dsp_pipe_init(DspPipe *enh,
                       const DspPipeCfg *cfg,
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

    /* Initialize LR4 crossover LP and HP biquads.
     * Use enh->cfg (the persistent copy) for coefficient pointers,
     * NOT &cfg which is stack-local and becomes dangling when this
     * function returns (the biquad stores a pointer, not a copy). */
    for (int i = 0; i < LR4_SECTIONS; i++) {
        biquad_init(&enh->left.lp_cross[i],  &enh->cfg.lp_cross_coeffs[i * 5]);
        biquad_init(&enh->right.lp_cross[i], &enh->cfg.lp_cross_coeffs[i * 5]);
        biquad_init(&enh->left.hp_cross[i],  &enh->cfg.hp_cross_coeffs[i * 5]);
        biquad_init(&enh->right.hp_cross[i], &enh->cfg.hp_cross_coeffs[i * 5]);
        biquad_init(&enh->left.hp_lookahead[i],  &enh->cfg.hp_cross_coeffs[i * 5]);
        biquad_init(&enh->right.hp_lookahead[i], &enh->cfg.hp_cross_coeffs[i * 5]);
    }

    /* Initialize HP for harmonics cleanup */
    biquad_init(&enh->left.hp_harm,  enh->cfg.hp_harm_coeffs);
    biquad_init(&enh->right.hp_harm, enh->cfg.hp_harm_coeffs);

    enh->left.env_smooth = 0.0f;
    enh->left.w_slew = 0.0f;
    enh->right.w_slew = 0.0f;
    enh->left.harm_amp_slew = 0.0f;
    enh->right.harm_amp_slew = 0.0f;
    enh->left.bass_gain_red = 1.0f;
    enh->right.bass_gain_red = 1.0f;
    enh->left.fb_env = 0.0f;
    enh->right.fb_env = 0.0f;


    /* Initialize envelope */
    env_init(&enh->left.env,  cfg->release_coeff);
    env_init(&enh->right.env, cfg->release_coeff);
}

void dsp_pipe_reset(DspPipe *enh) {
    DspPipeCfg *cfg = &enh->cfg;

    for (int i = 0; i < cfg->eq_n_biquads; i++) {
        biquad_reset(&enh->left.eq_bqs[i]);
        biquad_reset(&enh->right.eq_bqs[i]);
    }

    for (int i = 0; i < LR4_SECTIONS; i++) {
        biquad_reset(&enh->left.lp_cross[i]);
        biquad_reset(&enh->right.lp_cross[i]);
        biquad_reset(&enh->left.hp_cross[i]);
        biquad_reset(&enh->right.hp_cross[i]);
        biquad_reset(&enh->left.hp_lookahead[i]);
        biquad_reset(&enh->right.hp_lookahead[i]);
    }

    biquad_reset(&enh->left.hp_harm);
    biquad_reset(&enh->right.hp_harm);


    enh->left.env.peak  = 0.0f;
    enh->left.w_slew  = 0.0f;
    enh->right.w_slew = 0.0f;
    enh->left.harm_amp_slew  = 0.0f;
    enh->right.harm_amp_slew = 0.0f;
    enh->left.bass_gain_red  = 1.0f;
    enh->right.bass_gain_red = 1.0f;
    enh->left.fb_env  = 0.0f;
    enh->right.fb_env = 0.0f;

    for (int i = 0; i < LOOKAHEAD_LEN; i++) {
        enh->left.hp_ring[i]  = 0.0f;
        enh->right.hp_ring[i] = 0.0f;
        enh->left.eq_delay[i]  = 0.0f;
        enh->right.eq_delay[i] = 0.0f;
    }
    enh->left.hp_ring_pos  = 0;
    enh->right.hp_ring_pos = 0;
    enh->left.delay_pos    = 0;
    enh->right.delay_pos   = 0;
    enh->left.hp_max      = 0.0f;
    enh->left.hp_max_pos  = -1;
    enh->right.hp_max     = 0.0f;
    enh->right.hp_max_pos = -1;


    enh->left.loudness_x1  = 0.0f;
    enh->left.loudness_x2  = 0.0f;
    enh->left.loudness_y1  = 0.0f;
    enh->left.loudness_y2  = 0.0f;
    enh->right.loudness_x1 = 0.0f;
    enh->right.loudness_x2 = 0.0f;
    enh->right.loudness_y1 = 0.0f;
    enh->right.loudness_y2 = 0.0f;
}

/* ── Profiling ─────────────────────────────────────────────────────── */

DspPipeProfile dsp_profile;

void dsp_pipe_profile_report(void) {
    if (dsp_profile.frames == 0) return;
    uint32_t n = dsp_profile.frames;
    ESP_LOGI("enhancer", "=== %lu frames ===", (unsigned long)n);
    ESP_LOGI("enhancer", "  total: %llu cy/frame  (%llu cy total)",
             (unsigned long long)(dsp_profile.cycles_total / n),
             (unsigned long long)dsp_profile.cycles_total);
    ESP_LOGI("enhancer", "  EQ cascade: %llu cy/frame",
             (unsigned long long)(dsp_profile.cycles_eq / n));
    ESP_LOGI("enhancer", "  env+normalize: %llu cy/frame",
             (unsigned long long)(dsp_profile.cycles_env / n));
    ESP_LOGI("enhancer", "  Chebyshev+scale: %llu cy/frame",
             (unsigned long long)(dsp_profile.cycles_harm / n));
    ESP_LOGI("enhancer", "  HP+mix+headroom: %llu cy/frame",
             (unsigned long long)(dsp_profile.cycles_mix / n));
    ESP_LOGI("enhancer", "  I2S write: %llu cy/frame",
             (unsigned long long)(dsp_profile.cycles_i2s / n));
    memset(&dsp_profile, 0, sizeof(dsp_profile));
}

/* ── Per-channel tick ──────────────────────────────────────────────── */

static float dsp_pipe_process_channel(DspPipeChan *ch,
                                      const DspPipeCfg *cfg,
                                      float x)
{
#ifdef EQGEN_PROFILE
    uint32_t c0, c1, c2, c3;
    (void)c0; (void)c1; (void)c2; (void)c3;
    c0 = esp_cpu_get_cycle_count();
#endif


    /* ── Stage 0a: Pre-gain ─────────────────────────────────────── */
    if (cfg->pre_gain != 1.0f) {
        x *= cfg->pre_gain;
    }

    /* ── Stage 1: EQ preprocessing ──────────────────────────────── */
    float eq_out = biquad_cascade(ch->eq_bqs, cfg->eq_n_biquads, x);

    /* ── Loudness shelf (2nd-order biquad low shelf) ─────────────
     * Pre-enhancer — part of the psychoacoustic target.  At low
     * listening volumes, a Fletcher-Munson shelf compensates for
     * the ear's reduced bass sensitivity so the enhancer operates
     * on the full perceptual target. */
    {
        const float *c = cfg->loudness_coeffs;
        float y = c[0] * eq_out
                + c[1] * ch->loudness_x1
                + c[2] * ch->loudness_x2
                - c[3] * ch->loudness_y1
                - c[4] * ch->loudness_y2;
        ch->loudness_x2 = ch->loudness_x1;
        ch->loudness_x1 = eq_out;
        ch->loudness_y2 = ch->loudness_y1;
        ch->loudness_y1 = y;
        eq_out = y;
    }

    /* Bypass enhancer when cutoff ≤ 0 */
    if (cfg->cutoff_hz <= 0.0f) {
#ifdef EQGEN_PROFILE
        c1 = esp_cpu_get_cycle_count();
        dsp_profile.frames++;
        dsp_profile.cycles_total += (c1 - c0);
        dsp_profile.cycles_eq     += (c1 - c0);
#endif
        return eq_out;
    }

#ifdef EQGEN_PROFILE
    c1 = esp_cpu_get_cycle_count();
#endif

    /* ── Lookahead: LR4 HP filter on current sample ──────────── */
    float hp_now = biquad_cascade(ch->hp_lookahead, LR4_SECTIONS, eq_out);
    float hp_abs = fabsf(hp_now);

    /* Incremental running max — avoids O(LOOKAHEAD_LEN) scan every sample.
     * The ring position we're about to write is the oldest sample in the
     * window; when the current max lives there and gets overwritten, we
     * rescan.  Amortized O(1). */
    if (hp_abs >= ch->hp_max || ch->hp_max_pos < 0) {
        ch->hp_max = hp_abs;
        ch->hp_max_pos = ch->hp_ring_pos;
    } else if (ch->hp_max_pos == ch->hp_ring_pos) {
        /* Overwriting the max — rescan */
        ch->hp_max = 0.0f;
        for (int i = 0; i < LOOKAHEAD_LEN; i++) {
            if (ch->hp_ring[i] > ch->hp_max) {
                ch->hp_max = ch->hp_ring[i];
                ch->hp_max_pos = i;
            }
        }
        /* Also scan the incoming sample (not yet in ring) */
        if (hp_abs > ch->hp_max) {
            ch->hp_max = hp_abs;
            ch->hp_max_pos = ch->hp_ring_pos;
        }
    }

    ch->hp_ring[ch->hp_ring_pos] = hp_abs;
    ch->hp_ring_pos = (ch->hp_ring_pos + 1) % LOOKAHEAD_LEN;

    float upcoming = ch->hp_max;

    /* ── Delay line: get delayed eq_out, write current ──────────── */
    float eq_delayed = ch->eq_delay[(ch->delay_pos + 1) % LOOKAHEAD_LEN];
    ch->eq_delay[ch->delay_pos] = eq_out;
    ch->delay_pos = (ch->delay_pos + 1) % LOOKAHEAD_LEN;

    /* ── Process delayed signal: LR4 LP/HP split ─────────────── */
    float lp_fund = biquad_cascade(ch->lp_cross, LR4_SECTIONS, eq_delayed);
    float dry_hp   = biquad_cascade(ch->hp_cross, LR4_SECTIONS, eq_delayed);

    /* Envelope for Chebyshev normalization */
    float env = env_tick(&ch->env, lp_fund);
    /* Smooth env to kill 2f ripple: ripple AM in norm contaminates
     * the Chebyshev output with modulation sidebands that
     * interfere with the dry_hp fundamental, creating buzz. */
    ch->env_smooth += cfg->env_smooth_alpha * (env - ch->env_smooth);
    /* Adaptive blend: during steady state the smoothed envelope kills
     * 2f ripple (ripple amplitude ≈ 0.077·A < 0.1 threshold).  During
     * transients (sweeps, attacks) env jumps by >0.3 and the blend
     * switches to near-instant tracking, preventing norm from spiking
     * above 1.0 and kicking the Chebyshev or HP biquad states. */
    float env_diff = fabsf(env - ch->env_smooth);
    float blend = fminf(fmaxf((env_diff - 0.1f) * 5.0f, 0.0f), 1.0f);
    float env_norm = ch->env_smooth * (1.0f - blend) + env * blend;
    if (env_norm < 0.25f) env_norm = 0.25f;
    float norm = (env > 1e-6f) ? lp_fund / env_norm : 0.0f;

    /* ── Headroom budget ────────────────────────────────────────── */
    float room = cfg->push_gain * (1.0f - upcoming);
    if (room < 0.0f) room = 0.0f;

    /* ── Crossfade: fundamental ↔ harmonics ──────────────────────
     *
     * target = env — smoothed bass level (envelope of lp_fund).
     * h2_amp, h3_amp — individual perceptual efficiencies.
     *   h2 at unit physical amplitude delivers ~1/h2_amp times the
     *   perceived loudness of a fundamental at unit amplitude;
     *   h3 is similarly ~1/h3_amp times.
     *   (Small hN_amp → harmonic N is efficient → less physical needed.)
     * h_sum = h2_amp + h3_amp.
     * w      = fraction of fundamental loudness replaced by harmonics.
     *
     * Derivation (preserve perceived loudness = target):
     *   perceived = (1-w)·target + h2_phys/h2_amp + h3_phys/h3_amp = target
     *   physical  = (1-w)·target + h2_phys + h3_phys  = room
     *
     * Split proportional to efficiency:
     *   h2_phys = harm_amp · h2_amp/h_sum
     *   h3_phys = harm_amp · h3_amp/h_sum
     *   → h2_phys/h2_amp = h3_phys/h3_amp = harm_amp/h_sum
     *   → perceived = (1-w)·target + 2·harm_amp/h_sum = target
     *
     * Solving:
     *   harm_amp = (target - room) · h_sum / (2 - h_sum)
     *   w        = (target - room) / (target · (1 - h_sum/2))
     *
     * When target ≤ room: pure fundamental (w=0).
     * When target > room: blend to preserve perceived flatness.
     *
     * Output:  dry_hp + (1-w)·lp_fund + harm_out
     *   dry_hp = LR4 HP (24 dB/oct)
     *   lp_fund: LR4 LP (24 dB/oct)
     *   Chebyshev amplitudes come from h2_amp/h3_amp (computed
     *   by harmonic efficacy) — no separate budget clamp.
     *   |fund_out| ≤ (1-w)·target
     *   bass sum ≤ room by construction.
     */
    /* Perceived loudness is bounded by the tanh ceiling (1.0).
     * When overboost pushes env >> 1, the tanh clamps the output
     * and the extra level can't be perceived — cap target so the
     * crossfade doesn't try to generate harmonics for inaudible headroom. */
    float target = fminf(env, 1.0f);
    float h = cfg->h2_amp + cfg->h3_amp;
    float w_target;
    if (target <= room || h >= 1.0f || target < 1e-6f) {
        w_target = 0.0f;
    } else {
        w_target = (target - room) / (target * (1.0f - h * 0.5f));
        if (w_target > 1.0f) w_target = 1.0f;
    }
    /* Slew the crossfade weight to prevent the Chebyshev from turning
     * on/off in one sample during sweeps — a sudden w jump kicks the
     * HP harmonic biquad and sounds like a crackle.  Slew rate of
     * 0.002 (~500 samples / 11 ms for full transition) is fast enough
     * to track amplitude changes but slow enough to prevent clicks. */
    float w = ch->w_slew;
    float w_diff = w_target - w;
    if (w_diff > 0.002f) {
        w += 0.002f;
    } else if (w_diff < -0.002f) {
        w -= 0.002f;
    } else {
        w = w_target;
    }
    ch->w_slew = w;

    float fund_out = (1.0f - w) * lp_fund;

    float harm_out = 0.0f;
    if (w > 0.0f && h > 0.0f) {
        /* Harmonics fill the headroom gap, scaled by h_sum/(2-h_sum)
         * so perceived loudness = target while physical stays
         * within room.  Without this scaling, efficient harmonics
         * (small h2_amp, h3_amp) would produce harm_amp ≈ target - room,
         * over-filling headroom.
         * Slew harm_amp at the same rate as w to prevent the
         * Chebyshev from turning on in one sample during sweeps. */
        float h_safe = fmaxf(h, 1e-6f);
        float harm_amp_target = (target - room) * h_safe / (2.0f - h_safe);
        if (harm_amp_target < 0.0f) harm_amp_target = 0.0f;
        float harm_amp = ch->harm_amp_slew;
        float ha_diff = harm_amp_target - harm_amp;
        if (ha_diff > 0.002f) {
            harm_amp += 0.002f;
        } else if (ha_diff < -0.002f) {
            harm_amp -= 0.002f;
        } else {
            harm_amp = harm_amp_target;
        }
        ch->harm_amp_slew = harm_amp;

        float harm_t2 = 0.0f;
        if (cfg->h2_amp > 0.0f) {
            harm_t2 = cheb_t2(norm * cfg->h2_scale) * harm_amp;
        }
        float harm_t3 = 0.0f;
        if (cfg->h3_amp > 0.0f) {
            harm_t3 = cheb_t3(norm * cfg->h3_scale) * harm_amp;
            /* T3(x) = 4x³-3x.  With x = α·sin(ωt), fundamental =
             * (3α³-3α)·sin(ωt).  Cancel this analytically.
             * Scale the cancellation by harm_amp/env so it tracks
             * the slewed harmonic amplitude — otherwise the full
             * cancellation fires even when harm_amp ~ 0, injecting
             * a phantom fundamental into the HP filter. */
            float s = cfg->h3_scale;
            float harm_ratio = harm_amp / fmaxf(env, 1e-10f);
            if (harm_ratio > 1.0f) harm_ratio = 1.0f;
            float t3_fund = 3.0f * s * (s * s - 1.0f) * lp_fund * harm_ratio;
            harm_t3 -= t3_fund;
        }

        harm_out = biquad_tick(&ch->hp_harm, harm_t2 + harm_t3);
    } else {
        ch->harm_amp_slew = 0.0f;
    }

#ifdef EQGEN_PROFILE
    c2 = esp_cpu_get_cycle_count();
#endif

    /* ── Bass-sum limiter (replaces tanh) ────────────────────────
     * When dry_hp + bass_sum would exceed 1.0, attenuate only the
     * bass (fundamental + harmonics).  Dry HF passes through clean.
     * Instant attack / 49ms release — catches kick drum transients
     * without pumping.                                        */
    float bass_sum = fund_out + harm_out;
    float total_abs = fabsf(dry_hp + bass_sum);

    float target_gain = 1.0f;
    if (total_abs > 1.0f) {
        float hp_abs = fabsf(dry_hp);
        if (hp_abs >= 1.0f) {
            target_gain = 0.0f;
        } else {
            target_gain = (1.0f - hp_abs) / fmaxf(fabsf(bass_sum), 1e-10f);
            if (target_gain < 0.0f) target_gain = 0.0f;
        }
    }

    /* Envelope: instant attack, smoothed release toward 1.0 */
    if (target_gain < ch->bass_gain_red) {
        ch->bass_gain_red = target_gain;
    } else {
        ch->bass_gain_red += (1.0f - ch->bass_gain_red) * cfg->lim_release_coeff;
    }

    bass_sum *= ch->bass_gain_red;
    float out = dry_hp + bass_sum;


    /* ── Full-band peak limiter ─────────────────────────────────
     * Instant attack / 3-second release.  Only engages when the
     * bass-sum limiter can't protect against treble transients.
     * At normal gains this never fires — it's a safety net for
     * aggressive overboost configurations.                    */
    {
        float peak = fabsf(out);
        if (peak > ch->fb_env) {
            ch->fb_env = peak;
        } else {
            ch->fb_env += (peak - ch->fb_env) * cfg->fb_release_coeff;
        }
        if (ch->fb_env > 1.0f) {
            out /= ch->fb_env;
        }
    }

#ifdef EQGEN_PROFILE
    c3 = esp_cpu_get_cycle_count();
    dsp_profile.frames++;
    dsp_profile.cycles_total += (c3 - c0);
    dsp_profile.cycles_eq     += (c1 - c0);
    dsp_profile.cycles_env    += (c2 - c1);
    dsp_profile.cycles_mix    += (c3 - c2);
    dsp_profile.cycles_harm   += (c2 - c1) / 2;
#endif

    return out;
}

/* ── Stereo processing ─────────────────────────────────────────────── */

void dsp_pipe_process_stereo(DspPipe *enh,
                                 float *left, float *right)
{
    *left  = dsp_pipe_process_channel(&enh->left,  &enh->cfg, *left);
    *right = dsp_pipe_process_channel(&enh->right, &enh->cfg, *right);

    /* Safety net: clamp NaN/Inf to zero.  Once NaN enters biquad state
     * it propagates permanently — catching it here resets the output
     * but NOT the biquad state.  If NaN appears, dsp_init should also
     * reset biquads on next stream start. */
    if (!isfinite(*left))  *left  = 0.0f;
    if (!isfinite(*right)) *right = 0.0f;
}
