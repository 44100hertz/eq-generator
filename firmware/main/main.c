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

#include "enhancer.h"
#include "eq_coeffs.h"
#include "envelope.h"

#include "bt_a2dp.h"
#include "i2s_out.h"

static const char *TAG = "eqgen";

/* ── DSP instance (globally accessible from the BT callback) ──── */

static BassEnhancer      enhancer;
static BiquadQ28         eq_bqs_left[EQGEN_N_BIQUADS];
static BiquadQ28         eq_bqs_right[EQGEN_N_BIQUADS];
static ReciprocalLUT     rec_lut;

/* ── Track the last I2S sample rate ───────────────────────────── */

static int i2s_rate = 0;

/* ────────────────────────────────────────────────────────────────
 *  DSP init
 * ──────────────────────────────────────────────────────────────── */

static void dsp_init(void)
{
    ReciprocalLUT_init(&rec_lut);

    /* Set up DSP configuration from eq_coeffs.h defines */
    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         0.2f,     /* release_secs */
                         (float)EQGEN_FS,
                         0.049f,   /* limiter_release_secs */
                         EQGEN_N_BIQUADS,
                         eqgen_coeffs_q28);

    BassEnhancer_init(&enhancer, &cfg, &rec_lut,
                      eq_bqs_left, eq_bqs_right);

    ESP_LOGI(TAG, "DSP ready: %d biquads, fc=%.0f Hz, h2=%.2f h3=%.2f",
             EQGEN_N_BIQUADS, (double)EQGEN_CUTOFF_HZ,
             (double)EQGEN_H2_AMP, (double)EQGEN_H3_AMP);
}

/* ────────────────────────────────────────────────────────────────
 *  Audio data callback — called from the dsp_i2s task
 * ──────────────────────────────────────────────────────────────── */

static void audio_data_handler(const uint8_t *data, uint32_t len, int rate)
{
    /* Re-init I2S if sample rate changed (e.g. first connect, or
     * phone negotiated a different rate than expected). */
    if (rate != i2s_rate) {
        i2s_out_init(rate);
        i2s_rate = rate;
    }

    /* Process each stereo frame through the DSP */
    const int16_t *in = (const int16_t *)data;
    uint32_t frames = len / 4;   /* 2 channels × 2 bytes */

    for (uint32_t i = 0; i < frames; i++) {
        /* int16 → Q16 (left-shift by 1) */
        int32_t l = ((int32_t)in[i * 2])     << 1;
        int32_t r = ((int32_t)in[i * 2 + 1]) << 1;

        BassEnhancer_process_stereo(&enhancer, &l, &r);

        /* Q16 → int16 (right-shift by 1, clamp) */
        l >>= 1;
        r >>= 1;
        if (l >  32767) l =  32767;
        if (l < -32768) l = -32768;
        if (r >  32767) r =  32767;
        if (r < -32768) r = -32768;

        i2s_out_write_stereo((int16_t)l, (int16_t)r);
    }
}

/* ────────────────────────────────────────────────────────────────
 *  App entry
 * ──────────────────────────────────────────────────────────────── */

void app_main(void)
{
    ESP_LOGI(TAG, "eqgen firmware starting");

    /* Initialize DSP (uses baked-in eq_coeffs.h) */
    dsp_init();

    /* Start Bluetooth A2DP sink.
     * This also starts the DSP+I2S task internally.
     * I2S is initialized lazily when the first audio frame
     * arrives (so we know the negotiated sample rate). */
    bt_a2dp_sink_init("eqgen", audio_data_handler);

    /* Nothing else — FreeRTOS scheduler runs the BT stack
     * and the dsp_i2s task.  app_main returns and the idle
     * task reclaims this stack. */
    ESP_LOGI(TAG, "Firmware running — pair your phone now.");
}
