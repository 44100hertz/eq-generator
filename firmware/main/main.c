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

#include <math.h>
#include <stdint.h>

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#include "esp_log.h"
#include "esp_system.h"
#include "esp_timer.h"

#include "enhancer.h"
#include "eq_coeffs.h"
#include "envelope.h"
#include "fft_eq.h"
#include "smart_volume.h"

#include "bt_a2dp.h"
#include "i2s_out.h"
#include "sfx_data.h"
#include "rom/ets_sys.h"

static const char *TAG = "eqgen";

/* ── DSP instance (globally accessible from the BT callback) ──── */

static BassEnhancer      enhancer;
static Biquad             eq_bqs_left[EQGEN_N_BIQUADS];
static Biquad             eq_bqs_right[EQGEN_N_BIQUADS];

#if EQGEN_FFT_BINS > 0
static FftEq            *fft_eq;
static float             fft_pre_gain;   /* pre-gain applied in float before biquads */
#endif

/* ── Track the last I2S/DSP sample rate ───────────────────────── */

static int i2s_rate = 0;
static int dsp_rate = 0;

/* ── SFX event flags (set by Bluedroid task, cleared by DSP task) ── */

static volatile int sfx_pending = 0;  /* 1=connected, 2=disconnected */
static volatile uint32_t sfx_event_ms = 0;  /* last event timestamp (ms) */

/* ── Volume LUT: maps 0..127 → float gain (dB-linear) ────────── */
static float vol_lut[128];
static uint8_t last_vol = 255;  /* 255 = uninitialised, forces first update */

/* ────────────────────────────────────────────────────────────────
 *  DSP init
 * ──────────────────────────────────────────────────────────────── */

static void dsp_init(int rate)
{
    const float *coeffs = eqgen_get_coeffs(rate);
    float fs = (float)eqgen_get_fs(rate);

#if EQGEN_FFT_BINS > 0
#ifdef EQGEN_PRE_GAIN
    fft_pre_gain = EQGEN_PRE_GAIN;
#else
    fft_pre_gain = 1.0f;
#endif

    const float *fft_gains = eqgen_get_fft_gains(rate);
    if (fft_eq) fft_eq_destroy(fft_eq);
    fft_eq = fft_eq_create(fft_gains);
    if (!fft_eq) {
        ESP_LOGE(TAG, "FFT EQ allocation failed — bypassing FFT path");
    }
#endif

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_LIMITER_RELEASE_SECS,
#if defined(EQGEN_PRE_GAIN) && EQGEN_FFT_BINS == 0
                         EQGEN_PRE_GAIN,
#else
                         1.0f,
#endif
                         EQGEN_N_BIQUADS,
                         coeffs);

    BassEnhancer_init(&enhancer, &cfg,
                      eq_bqs_left, eq_bqs_right);

    /* Configure loudness shelf (corner: 200 Hz, max boost: 8 dB). */
    BassEnhancerCfg_set_loudness(&cfg, EQGEN_LOUDNESS_FC_HZ, fs, 0.0f);
    enhancer.cfg = cfg;

    dsp_rate = rate;

#if EQGEN_FFT_BINS > 0
    ESP_LOGI(TAG, "DSP ready @ %d Hz: FFT=%d-pt + %d biquads, fc=%.0f Hz, h2=%.2f h3=%.2f",
             rate, EQGEN_FFT_N, EQGEN_N_BIQUADS, (double)EQGEN_CUTOFF_HZ,
             (double)EQGEN_H2_AMP, (double)EQGEN_H3_AMP);
#else
    ESP_LOGI(TAG, "DSP ready @ %d Hz: %d biquads, fc=%.0f Hz, h2=%.2f h3=%.2f",
             rate, EQGEN_N_BIQUADS, (double)EQGEN_CUTOFF_HZ,
             (double)EQGEN_H2_AMP, (double)EQGEN_H3_AMP);
#endif
}

/* ── SFX playback helper ────────────────────────────────────────── */

#define SFX_BATCH 256  /* stereo frames per DMA write */
#define SFX_FADE   200  /* samples to fade in/out */

static void play_sfx(const int16_t *mono, uint32_t n_samples, int rate)
{
    int16_t stereo[SFX_BATCH * 2];

    uint8_t vol = bt_a2dp_get_volume();
    float vol_f = vol ? (float)vol / 127.0f : 0.0f;

    for (uint32_t off = 0; off < n_samples; off += SFX_BATCH) {
        uint32_t chunk = SFX_BATCH;
        if (off + chunk > n_samples) chunk = n_samples - off;

        for (uint32_t i = 0; i < chunk; i++) {
            float s = (float)mono[off + i] / 32768.0f;  /* int16 → float */
            s *= 0.25f;  /* −12 dB */

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

#if EQGEN_FFT_BINS > 0
    static float  fft_buf_l[EQGEN_FFT_HOP];
    static float  fft_buf_r[EQGEN_FFT_HOP];
    static float  fft_out_l[EQGEN_FFT_HOP];
    static float  fft_out_r[EQGEN_FFT_HOP];
    static int    fft_buf_count = 0;
    static int    fft_buf_rate = 0;

    if (rate != fft_buf_rate) {
        fft_buf_count = 0;
        fft_buf_rate = rate;
    }
#endif

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
    float vol_f = vol_lut[vol];

#if EQGEN_FFT_BINS > 0
    for (uint32_t i = 0; i < frames; i++) {
        /* Accumulate: int16 → float [-1, 1] */
        fft_buf_l[fft_buf_count] = (float)in[i * 2]     / 32768.0f;
        fft_buf_r[fft_buf_count] = (float)in[i * 2 + 1] / 32768.0f;
        fft_buf_count++;

        if (fft_buf_count == EQGEN_FFT_HOP) {
            if (fft_eq) {
                fft_eq_process_frame(fft_eq, fft_buf_l, fft_buf_r,
                                     fft_out_l, fft_out_r);
            } else {
                for (int j = 0; j < EQGEN_FFT_HOP; j++) {
                    fft_out_l[j] = fft_buf_l[j];
                    fft_out_r[j] = fft_buf_r[j];
                }
            }

            /* Apply pre-gain in float before biquads */
            for (int j = 0; j < EQGEN_FFT_HOP; j++) {
                fft_out_l[j] *= fft_pre_gain;
                fft_out_r[j] *= fft_pre_gain;
            }

            /* Feed FFT output through IIR biquads + bass enhancer */
            for (int j = 0; j < EQGEN_FFT_HOP; j++) {
                float l = fft_out_l[j];
                float r = fft_out_r[j];

                BassEnhancer_process_stereo(&enhancer, &l, &r);

                /* Apply BT system volume */
                l *= vol_f;
                r *= vol_f;

                /* float → int16, clamp */
                int32_t li = (int32_t)(l * 32768.0f);
                int32_t ri = (int32_t)(r * 32768.0f);
                if (li >  32767) li =  32767;
                if (li < -32768) li = -32768;
                if (ri >  32767) ri =  32767;
                if (ri < -32768) ri = -32768;

                i2s_out_write_stereo((int16_t)li, (int16_t)ri);
            }

            fft_buf_count = 0;
        }
    }
#else
    /* ── No-FFT path: sample-by-sample IIR + enhancer ── */
    for (uint32_t i = 0; i < frames; i++) {
        /* int16 → float [-1, 1] */
        float l = (float)in[i * 2]     / 32768.0f;
        float r = (float)in[i * 2 + 1] / 32768.0f;

        BassEnhancer_process_stereo(&enhancer, &l, &r);

        /* Apply BT system volume */
        l *= vol_f;
        r *= vol_f;

        /* float → int16, clamp */
        int32_t li = (int32_t)(l * 32768.0f);
        int32_t ri = (int32_t)(r * 32768.0f);
        if (li >  32767) li =  32767;
        if (li < -32768) li = -32768;
        if (ri >  32767) ri =  32767;
        if (ri < -32768) ri = -32768;

        i2s_out_write_stereo((int16_t)li, (int16_t)ri);
    }
#endif

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

    float pg_loud = 1.0f;
#ifdef EQGEN_PRE_GAIN
    pg_loud = EQGEN_PRE_GAIN;
#endif

    SmartVolumeParams svp;
    smart_volume_compute(vol, pg_loud, &svp);

#if EQGEN_FFT_BINS > 0
    fft_pre_gain = svp.fft_pre_gain;
    BassEnhancer_update_params(&enhancer, svp.h2_amp, svp.h3_amp,
                               NAN, svp.boost, svp.bleed);
#else
    BassEnhancer_update_params(&enhancer, svp.h2_amp, svp.h3_amp,
                               svp.pre_gain, svp.boost, svp.bleed);
#endif

    smart_volume_rebuild_lut(vol_lut, svp.shelf_db);

    last_vol = vol;

    ESP_LOGI(TAG, "Smart vol: vol=%u t=%.2f h2=%.3f h3=%.3f shelf=%.1f dB pg=%.3f bleed=%.3f",
             (unsigned)vol,
             (double)(float)vol / 127.0,
             (double)svp.h2_amp,
             (double)svp.h3_amp,
             (double)svp.shelf_db,
             (double)svp.pre_gain,
             (double)svp.bleed);
}

/* ────────────────────────────────────────────────────────────────
 *  App entry
 * ──────────────────────────────────────────────────────────────── */

void app_main(void)
{
    ESP_LOGI(TAG, "eqgen firmware starting (float DSP)");

    bt_a2dp_set_event_callback(bt_event_handler);
    bt_a2dp_sink_init(EQGEN_BT_DEVICE_NAME, audio_data_handler);

    smart_volume_rebuild_lut(vol_lut, 0.0f);

    ESP_LOGI(TAG, "Firmware running — pair your phone now.");
}
