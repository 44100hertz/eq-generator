/**
 * main.c — eqgen firmware entry point
 *
 * Wires Bluetooth A2DP sink → DSP (enhancer + EQ) → I2S DAC.
 *
 * The DSP chain runs at compile-time from eq_coeffs.h.
 * To change the EQ curve, re-run the desktop export tool and re-flash.
 */

#include <stdint.h>

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "esp_system.h"
#include "esp_timer.h"

#include "enhancer.h"
#include "eq_coeffs.h"
#include "envelope.h"

#include "bt_a2dp.h"
#include "i2s_out.h"
#include "sfx_data.h"
#include "rom/ets_sys.h"

static const char *TAG = "eqgen";

/* ── DSP instance (globally accessible from the BT callback) ──── */

static BassEnhancer      enhancer;
static BiquadQ28         eq_bqs_left[EQGEN_N_BIQUADS];
static BiquadQ28         eq_bqs_right[EQGEN_N_BIQUADS];
static ReciprocalLUT     rec_lut;

/* ── Track the last I2S/DSP sample rate ───────────────────────── */

static int i2s_rate = 0;
static int dsp_rate = 0;

/* ── SFX event flags (set by Bluedroid task, cleared by DSP task) ── */

static volatile int sfx_pending = 0;  /* 1=connected, 2=disconnected */
static volatile uint32_t sfx_event_ms = 0;  /* last event timestamp (ms) */

/* ────────────────────────────────────────────────────────────────
 *  DSP init
 * ──────────────────────────────────────────────────────────────── */

static void dsp_init(int rate)
{
    const int32_t *coeffs = eqgen_get_coeffs(rate);
    float fs = (float)eqgen_get_fs(rate);

    ReciprocalLUT_init(&rec_lut);

    /* Set up DSP configuration from eq_coeffs.h defines */
    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_LIMITER_RELEASE_SECS,
#ifdef EQGEN_PRE_GAIN_Q16
                         (float)EQGEN_PRE_GAIN_Q16 / 65536.0f,
#else
                         1.0f,
#endif
                         EQGEN_N_BIQUADS,
                         coeffs);

    BassEnhancer_init(&enhancer, &cfg, &rec_lut,
                      eq_bqs_left, eq_bqs_right);

    dsp_rate = rate;

    ESP_LOGI(TAG, "DSP ready @ %d Hz: %d biquads, fc=%.0f Hz, h2=%.2f h3=%.2f",
             rate, EQGEN_N_BIQUADS, (double)EQGEN_CUTOFF_HZ,
             (double)EQGEN_H2_AMP, (double)EQGEN_H3_AMP);
}

/* ── SFX playback helper — batch-write to avoid per-sample DMA overhead ── */

#define SFX_BATCH 256  /* stereo frames per DMA write */
#define SFX_FADE   200  /* samples to fade in/out */

static void play_sfx(const int16_t *mono, uint32_t n_samples, int rate)
{
    int16_t stereo[SFX_BATCH * 2];

    /* Read AVRCP absolute volume once for the SFX duration.
     * Q16 gain: 0=mute, 65536=unity (vol 127). */
    uint8_t vol = bt_a2dp_get_volume();
    int32_t vol_q16 = vol ? ((int32_t)vol << 16) / 127 : 0;

    for (uint32_t off = 0; off < n_samples; off += SFX_BATCH) {
        uint32_t chunk = SFX_BATCH;
        if (off + chunk > n_samples) chunk = n_samples - off;

        for (uint32_t i = 0; i < chunk; i++) {
            int32_t s = mono[off + i];
            s = s >> 2;  /* −12 dB */

            /* Fade in/out to avoid amplitude discontinuity */
            uint32_t pos = off + i;
            if (pos < SFX_FADE) {
                s = (s * (int32_t)pos) / SFX_FADE;
            } else if (pos >= n_samples - SFX_FADE) {
                s = (s * (int32_t)(n_samples - 1 - pos)) / SFX_FADE;
            }

            /* Apply BT system volume */
            s = (s * vol_q16) >> 16;

            stereo[i * 2]     = (int16_t)s;
            stereo[i * 2 + 1] = (int16_t)s;
        }
        /* DMA backpressure (portMAX_DELAY) paces writes at the I2S rate.
         * No busy-wait needed — the DMA buffer acts as a FIFO. */
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
        /* Debounce: wait 200ms after last event for rapid-fire
         * connect/disconnect cycles to settle (common in A2DP init).
         * Audio data is NOT discarded during the wait — it flows through
         * normally.  The SFX is simply deferred until the debounce
         * timer expires. */
        int64_t age_ms = ((int64_t)esp_timer_get_time() / 1000) - (int64_t)sfx_event_ms;
        if (age_ms >= 200) {
            int pending = sfx_pending;
            sfx_pending = 0;

            /* On connect (audio may start soon): init I2S at BT rate.
             * On disconnect (audio already stopped): use last known rate. */
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

    /* No audio data — 50ms poll tick for events */
    if (len == 0) return;

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

    /* Process each stereo frame through the DSP */
    const int16_t *in = (const int16_t *)data;
    uint32_t frames = len / 4;   /* 2 channels × 2 bytes */

    /* Read AVRCP absolute volume once per chunk.
     * Q16 gain: 0=mute, 65536=unity (vol 127). */
    uint8_t vol = bt_a2dp_get_volume();
    int32_t vol_q16 = vol ? ((int32_t)vol << 16) / 127 : 0;

    for (uint32_t i = 0; i < frames; i++) {
        /* int16 → Q16 (left-shift by 1) */
        int32_t l = ((int32_t)in[i * 2])     << 1;
        int32_t r = ((int32_t)in[i * 2 + 1]) << 1;

        BassEnhancer_process_stereo(&enhancer, &l, &r);

        /* Apply BT system volume while still in Q16 for precision */
        l = (l * vol_q16) >> 16;
        r = (r * vol_q16) >> 16;

        /* Q16 → int16 (right-shift by 1, clamp) */
        l >>= 1;
        r >>= 1;
        if (l >  32767) l =  32767;
        if (l < -32768) l = -32768;
        if (r >  32767) r =  32767;
        if (r < -32768) r = -32768;

        i2s_out_write_stereo((int16_t)l, (int16_t)r);
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

/* ────────────────────────────────────────────────────────────────
 *  App entry
 * ──────────────────────────────────────────────────────────────── */

void app_main(void)
{
    ESP_LOGI(TAG, "eqgen firmware starting");

    /* DSP is initialized lazily when the first audio frame
     * arrives (so we know the negotiated sample rate).
     * This ensures the correct biquad coefficients (44.1k or 48k)
     * and HP/LP filter tuning are selected at runtime. */

    /* Start Bluetooth A2DP sink.
     * This also starts the DSP+I2S task internally.
     * I2S is initialized lazily when the first audio frame
     * arrives (so we know the negotiated sample rate). */
    bt_a2dp_set_event_callback(bt_event_handler);
    bt_a2dp_sink_init(EQGEN_BT_DEVICE_NAME, audio_data_handler);

    /* Nothing else — FreeRTOS scheduler runs the BT stack
     * and the dsp_i2s task.  app_main returns and the idle
     * task reclaims this stack. */
    ESP_LOGI(TAG, "Firmware running — pair your phone now.");
}
