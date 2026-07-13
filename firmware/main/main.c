/**
 * main.c — eqgen firmware entry point
 *
 * Wires Bluetooth A2DP sink → DSP (enhancer + EQ) → I2S DAC.
 *
 * The DSP chain runs at compile-time from eq_coeffs.h.
 * To change the EQ curve, re-run the desktop export tool and re-flash.
 *
 * All DSP is float now — no more Q16/Q28 integer math.
 */

/* ── Profiling: comment out for production builds ─────────────── */
// #define EQGEN_PROFILE

#include <math.h>
#include <stdbool.h>
#include <stdint.h>

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "esp_system.h"
#include "esp_timer.h"

#include "enhancer.h"
#include "eq_coeffs.h"
#include "envelope.h"
#include "smart_volume.h"
#include "volume_control.h"

#include "bt_a2dp.h"
#include "i2s_out.h"
#include "sfx_data.h"
#include "rom/ets_sys.h"

static const char *TAG = "eqgen";

/* ── DSP instance (globally accessible from the BT callback) ──── */

static BassEnhancer      enhancer;
static Biquad             eq_bqs_left[EQGEN_N_BIQUADS];
static Biquad             eq_bqs_right[EQGEN_N_BIQUADS];

/* ── Track the last I2S/DSP sample rate ───────────────────────── */

static int i2s_rate = 0;
static int dsp_rate = 0;

/* ── Stream-start detection for biquad state reset ────────────── */
static bool audio_streaming = false;

/* ── SFX event flags (set by Bluedroid task, cleared by DSP task) ── */

static volatile int sfx_pending = 0;  /* 1=connected, 2=disconnected */
static volatile uint32_t sfx_event_ms = 0;  /* last event timestamp (ms) */

/* ── Volume LUT: maps 0..127 → float gain (dB-linear) ────────── */
static float vol_lut[128];
static uint8_t last_vol = 255;  /* 255 = uninitialised, forces first update */

/* ── Forward decls ─────────────────────────────────────────────── */
static void update_smart_volume(int rate, uint8_t vol);

/* ────────────────────────────────────────────────────────────────
 *  DSP init
 * ──────────────────────────────────────────────────────────────── */

static void dsp_init(int rate)
{
    const float *coeffs = eqgen_get_coeffs(rate);
    float fs = (float)eqgen_get_fs(rate);

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_PUSH_GAIN,
#ifdef EQGEN_PRE_GAIN
                         EQGEN_PRE_GAIN,
#else
                         1.0f,
#endif
                         EQGEN_N_BIQUADS,
                         coeffs);

    BassEnhancer_init(&enhancer, &cfg,
                      eq_bqs_left, eq_bqs_right);

    /* Pre-compute loudness shelf alpha now so the one-pole filter
     * is ready.  Boost starts at 0 — update_smart_volume() will
     * push the correct boost on first volume change.
     *
     * Do NOT use BassEnhancerCfg_set_loudness(..., 0.0f) here:
     * that early-returns with alpha=0, permanently disabling the
     * shelf regardless of later boost updates. */
    cfg.loudness_alpha = 1.0f - expf(-2.0f * (float)M_PI * EQGEN_LOUDNESS_FC_HZ / fs);
    cfg.loudness_boost = 0.0f;
    enhancer.cfg = cfg;

    dsp_rate = rate;

    ESP_LOGI(TAG, "DSP ready @ %d Hz: %d biquads, fc=%.0f Hz, h2=%.2f h3=%.2f",
             rate, EQGEN_N_BIQUADS, (double)EQGEN_CUTOFF_HZ,
             (double)EQGEN_H2_AMP, (double)EQGEN_H3_AMP);
}

/* ── SFX playback helper ────────────────────────────────────────── */

#define SFX_BATCH 256  /* stereo frames per DMA write */
#define SFX_FADE   200  /* samples to fade in/out */

/* ── I2S batch write: accumulate into buffer, flush in chunks.
 *   Reduces i2s_channel_write calls from ~44k/sec to ~344/sec. ── */
#define BATCH_SIZE  256  /* stereo frames per I2S flush */

static void play_sfx(const int16_t *mono, uint32_t n_samples, int rate)
{
    int16_t stereo[SFX_BATCH * 2];

    uint8_t vol = bt_a2dp_get_volume();
    float vol_f = volume_gain(vol_lut, vol);   /* dB-linear LUT — same as audio path */

    /* Cap SFX to prevent blasting on high-gain speakers.
     * SFX never plays louder than -(speaker_level - 20) dB FS. */
    float sfx_cap_linear = powf(10.0f, -(EQGEN_SPEAKER_LEVEL_DB - 20.0f) / 20.0f);
    vol_f = fminf(vol_f, sfx_cap_linear);

    for (uint32_t off = 0; off < n_samples; off += SFX_BATCH) {
        uint32_t chunk = SFX_BATCH;
        if (off + chunk > n_samples) chunk = n_samples - off;

        for (uint32_t i = 0; i < chunk; i++) {
            float s = (float)mono[off + i] / 32768.0f;  /* int16 → float */

            /* Apply enhancer pre-gain so SFX level tracks audio level.
             * Audio goes through pre_gain → EQ → mix → LUT.
             * SFX bypasses the enhancer, so we apply both here. */
#ifdef EQGEN_PRE_GAIN
            s *= EQGEN_PRE_GAIN;
#else
            s *= 0.5f;   /* −6 dB headroom fallback */
#endif

            /* Fade in/out */
            uint32_t pos = off + i;
            if (pos < SFX_FADE) {
                s *= (float)pos / (float)SFX_FADE;
            } else if (pos >= n_samples - SFX_FADE) {
                s *= (float)(n_samples - 1 - pos) / (float)SFX_FADE;
            }

            /* Apply BT system volume */
            s *= vol_f;

            /* Clamp and convert to int16 */
            int32_t si = (int32_t)(s * 32768.0f);
            if (si >  32767) si =  32767;
            if (si < -32768) si = -32768;

            stereo[i * 2]     = (int16_t)si;
            stereo[i * 2 + 1] = (int16_t)si;
        }
        i2s_out_write(stereo, chunk, NULL);
    }
    (void)rate;
}

/* ────────────────────────────────────────────────────────────────
 *  Audio data callback — called from the dsp_i2s task
 * ──────────────────────────────────────────────────────────────── */

static void bt_event_handler(bt_event_t event)
{
    if (event == BT_EVENT_CONNECTED) {
        sfx_pending = 1;
    } else if (event == BT_EVENT_DISCONNECTED) {
        sfx_pending = 2;
    }
    sfx_event_ms = (uint32_t)(esp_timer_get_time() / 1000);
}

static void audio_data_handler(const uint8_t *data, uint32_t len, int rate)
{
    /* ── Play pending SFX before anything else ── */
    if (sfx_pending) {
        int64_t age_ms = ((int64_t)esp_timer_get_time() / 1000) - (int64_t)sfx_event_ms;
        if (age_ms >= 200) {
            int pending = sfx_pending;
            sfx_pending = 0;

            int play_rate = (i2s_rate > 0) ? i2s_rate : rate;
            if (play_rate != i2s_rate) {
                i2s_out_init(play_rate);
                i2s_rate = play_rate;
            }

#ifdef SFX_CONNECTED_SAMPLES
            if (pending & 1) {
                play_sfx(sfx_connected, SFX_CONNECTED_SAMPLES, play_rate);
            }
#endif
#ifdef SFX_DISCONNECTED_SAMPLES
            if (pending & 2) {
                play_sfx(sfx_disconnected, SFX_DISCONNECTED_SAMPLES, play_rate);
            }
#endif
        }
    }

    /* ── Error feedback state for float→int16 noise shaping ── */
    static float  l_err = 0.0f, r_err = 0.0f;

    /* No audio data — 50ms poll tick for events */
    if (len == 0) {
        audio_streaming = false;
        return;
    }

    /* Reset biquad/envelope state on stream start to flush any
     * stale state from a previous connection or silent gap. */
    if (!audio_streaming) {
        BassEnhancer_reset(&enhancer);
        l_err = 0.0f;
        r_err = 0.0f;
        audio_streaming = true;
        ESP_LOGI(TAG, "Stream start — enhancer state reset");
    }

    /* Re-init I2S and DSP if sample rate changed */
    if (rate != i2s_rate) {
        i2s_out_init(rate);
        i2s_rate = rate;
    }
    if (rate != dsp_rate) {
        dsp_init(rate);
    }

    /* ── DSP processing ── */
#ifdef EQGEN_PROFILE
    int64_t t0 = esp_timer_get_time();
#endif

    const int16_t *in = (const int16_t *)data;
    uint32_t frames = len / 4;   /* 2 channels × 2 bytes */

    uint8_t vol = bt_a2dp_get_volume();
    if (vol != last_vol) {
        update_smart_volume(rate, vol);
    }
    float vol_f = volume_gain(vol_lut, vol);

    /* ── Batch buffer for I2S writes: processes into local
     *   buffer, flushes in chunks to amortize DMA overhead. ── */
    static int16_t out_buf[BATCH_SIZE * 2];
    static uint32_t out_idx = 0;

    for (uint32_t i = 0; i < frames; i++) {
        /* int16 → float [-1, 1] */
        float l = (float)in[i * 2]     / 32768.0f;
        float r = (float)in[i * 2 + 1] / 32768.0f;

        /* Apply BT system volume before enhancer so the limiter
         * only engages at high volumes. */
        l *= vol_f;
        r *= vol_f;

        BassEnhancer_process_stereo(&enhancer, &l, &r);

        /* Safety net: if enhancer produced NaN/Inf, zero it.
         * This prevents undefined behavior in the float→int cast. */
        if (!isfinite(l)) l = 0.0f;
        if (!isfinite(r)) r = 0.0f;

        /* float → int16 with first-order error feedback.
         * Leaky integrator diffuses truncation error into HF
         * where it's inaudible — ~43 dB quieter in-band at -80 dB. */
        float lf = l * 32768.0f + 0.999f * l_err;
        float rf = r * 32768.0f + 0.999f * r_err;
        int32_t li = (int32_t)lf;
        int32_t ri = (int32_t)rf;
        if (li >= -32768 && li <= 32767) {
            l_err = lf - (float)li;
        } else {
            li = (li > 32767) ? 32767 : -32768;
            l_err = 0.0f;
        }
        if (ri >= -32768 && ri <= 32767) {
            r_err = rf - (float)ri;
        } else {
            ri = (ri > 32767) ? 32767 : -32768;
            r_err = 0.0f;
        }

        out_buf[out_idx++] = (int16_t)li;
        out_buf[out_idx++] = (int16_t)ri;

        /* Flush batch to I2S */
        if (out_idx >= BATCH_SIZE * 2) {
            i2s_out_write(out_buf, BATCH_SIZE, NULL);
            out_idx = 0;
        }
    }

    /* Flush remainder */
    if (out_idx > 0) {
        uint32_t rem_frames = out_idx / 2;
        i2s_out_write(out_buf, rem_frames, NULL);
        out_idx = 0;
    }

#ifdef EQGEN_PROFILE
    int64_t elapsed_us = esp_timer_get_time() - t0;
    static int64_t total_us = 0, total_frames = 0;
    total_us += elapsed_us;
    total_frames += frames;
    if (total_frames >= 44100) {
        float load_pct = (float)total_us / 10000.0f;
        float us_per_frame = (float)total_us / (float)total_frames;
        ESP_LOGI(TAG, "DSP: %.1f%% CPU, %.1f us/frame (%lu frames)",
                 load_pct, us_per_frame, (unsigned long)total_frames);
        enhancer_profile_report();
        total_us = 0;
        total_frames = 0;
    }
#endif
}

/* ── Smart volume: interpolate between quiet & loud preset ───────── */

static void update_smart_volume(int rate, uint8_t vol)
{
    if (dsp_rate <= 0) return;

    SmartVolumeParams svp = volume_set(vol, vol_lut, &enhancer);
    last_vol = vol;

    ESP_LOGI(TAG, "Smart vol: vol=%u t=%.2f shelf=%.1f dB pg=%.3f",
             (unsigned)vol,
             (double)(float)vol / 127.0,
             (double)svp.shelf_db,
             (double)svp.pre_gain);
}

/* ────────────────────────────────────────────────────────────────
 *  App entry
 * ──────────────────────────────────────────────────────────────── */

void app_main(void)
{
    ESP_LOGI(TAG, "eqgen firmware starting (float DSP)");

    bt_a2dp_set_event_callback(bt_event_handler);
    bt_a2dp_sink_init(EQGEN_BT_DEVICE_NAME, audio_data_handler);

    volume_init_lut(vol_lut);

    ESP_LOGI(TAG, "Firmware running — pair your phone now.");
}
