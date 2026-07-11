/**
 * quick_bench.c — Float biquad cycle-count sweep for ESP32
 *
 * Uses the real biquad_tick() from ../../src/biquad.h so numbers are
 * directly applicable to the DSP chain.
 *
 * Output: a table showing cycles/sample at each biquad count, so we can
 * pick the largest N that fits inside the 48 kHz stereo budget.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "esp_timer.h"
#include "esp_cpu.h"
#include "esp_log.h"

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

/* Pull in the real biquad implementation */
#include "biquad.h"

/* ── Config ───────────────────────────────────────────────────────── */

#define FS           48000        /* sample rate (Hz) */
#define N_SAMP        512         /* samples per benchmark run */
#define WARMUP_SAMP   128         /* warmup samples (separate pass) */
#define ITERS          64         /* repeats; best / worst reported */
#define MAX_BQ         80         /* highest sweep point */

#define CPU_MHZ       240

/* ── Helpers ──────────────────────────────────────────────────────── */

static inline float randf(void) {
    return (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
}

/* ── Benchmark a cascade of N biquads ──────────────────────────────── */

/**
 * Run ITERS timed cascades; return best & worst cycle counts.
 * Coefficients: peaking filter at 1 kHz, Q=1.414, +6 dB — a
 * representative "busy" EQ band.  (Peaking has the most non-zero
 * coefficients, so it's a worst-case multiply load.)
 */
static void bench_n(const float *coeffs, int n,
                    uint32_t *best, uint32_t *worst)
{
    /* Allocate biquads + sample buffers on heap */
    Biquad *bqs = calloc((size_t)n, sizeof(Biquad));
    float  *in  = malloc((size_t)N_SAMP * sizeof(float));
    float  *out = malloc((size_t)N_SAMP * sizeof(float));
    if (!bqs || !in || !out) {
        ESP_LOGE("bench", "alloc failed for n=%d", n);
        if (best)  *best  = 0;
        if (worst) *worst = 0;
        free(bqs); free(in); free(out);
        return;
    }

    for (int i = 0; i < n; i++) biquad_init(&bqs[i], coeffs);
    for (int i = 0; i < N_SAMP; i++) in[i] = randf();

    /* Warmup */
    for (int i = 0; i < WARMUP_SAMP; i++) {
        float y = randf();
        for (int j = 0; j < n; j++) y = biquad_tick(&bqs[j], y);
    }

    uint32_t best_cyc  = UINT32_MAX;
    uint32_t worst_cyc = 0;

    volatile float sink = 0.0f;
    for (int iter = 0; iter < ITERS; iter++) {
        uint32_t t0 = esp_cpu_get_cycle_count();
        for (int i = 0; i < N_SAMP; i++) {
            float y = in[i];
            for (int j = 0; j < n; j++) y = biquad_tick(&bqs[j], y);
            out[i] = y;
        }
        uint32_t t1 = esp_cpu_get_cycle_count();
        sink += out[N_SAMP - 1];  /* prevent dead-code elimination */
        uint32_t dt = t1 - t0;
        if (dt < best_cyc)  best_cyc  = dt;
        if (dt > worst_cyc) worst_cyc = dt;
    }
    (void)sink;

    *best  = best_cyc;
    *worst = worst_cyc;

    free(bqs);
    free(in);
    free(out);
}

/* ── Generate a representative peaking EQ coefficient set ─────────── */

static void make_coeffs(float fc, float Q, float gain_db, float fs,
                        float out[5])
{
    float A  = powf(10.0f, gain_db / 40.0f);
    float w0 = 2.0f * (float)M_PI * fc / fs;
    float alpha = sinf(w0) / (2.0f * Q);

    float b0 =  1.0f + alpha * A;
    float b1 = -2.0f * cosf(w0);
    float b2 =  1.0f - alpha * A;
    float a0 =  1.0f + alpha / A;
    float a1 = -2.0f * cosf(w0);
    float a2 =  1.0f - alpha / A;

    float norm = 1.0f / a0;
    out[0] = b0 * norm;
    out[1] = b1 * norm;
    out[2] = b2 * norm;
    out[3] = a1 * norm;
    out[4] = a2 * norm;
}

/* ── Main ─────────────────────────────────────────────────────────── */

void app_main(void)
{
    ESP_LOGI("bench", "=== Float biquad sweep, %d MHz ===", CPU_MHZ);

    /* one peaking EQ band at 1 kHz */
    float coeffs[5];
    make_coeffs(1000.0f, 1.414f, 6.0f, (float)FS, coeffs);

    /* Budget: 48k stereo = 96k samples/sec.
     * 240e6 / 96e3 = 2500 cycles/sample available. */
    float budget_us = 1e6f / (float)(FS * 2);  /* ≈ 10.42 µs per stereo pair */
    float budget_cy = (float)(CPU_MHZ * 1000000) / (float)(FS * 2);

    ESP_LOGI("bench", "Budget: %.1f cy/stereo-sample (%.2f us)",
             budget_cy, budget_us);
    ESP_LOGI("bench", "");
    ESP_LOGI("bench", " N_BQ | cy/bq | cy/samp(1ch) | us/samp(1ch) | cy/samp(2ch) | CPU%%(2ch)");
    ESP_LOGI("bench", "------+-------+--------------+--------------+--------------+----------");

    /* Sweep from 1 up to MAX_BQ */
    int sweep[] = {1,2,4,6,8,10,12,14,16,20,24,28,32,36,40,48,56,64,72,80};
    int n_pts = sizeof(sweep) / sizeof(sweep[0]);

    for (int si = 0; si < n_pts; si++) {
        int n = sweep[si];
        uint32_t best, worst;
        bench_n(coeffs, n, &best, &worst);
        if (best == 0 && worst == 0) continue;

        float cy_per_bq   = (float)best / (float)(N_SAMP * n);
        float cy_per_samp = (float)best / (float)N_SAMP;     /* 1 channel cascade */
        float us_per_samp = cy_per_samp / (float)CPU_MHZ;
        float cy_stereo   = cy_per_samp * 2.0f;               /* 2 channels */
        float cpu_stereo  = cy_stereo / budget_cy * 100.0f;
        (void)worst;

        ESP_LOGI("bench", " %4d | %5.1f | %12.0f | %12.2f | %12.0f | %7.1f%%",
                 n, cy_per_bq, cy_per_samp, us_per_samp, cy_stereo, cpu_stereo);
    }

    ESP_LOGI("bench", "");
    ESP_LOGI("bench", "Budget: %.0f cy / stereo sample (48 kHz)", budget_cy);

    ESP_LOGI("bench", "Done.");
    while (1) { vTaskDelay(pdMS_TO_TICKS(5000)); }
}
