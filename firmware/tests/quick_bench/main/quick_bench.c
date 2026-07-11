/**
 * quick_bench.c — Float biquad cycle counts on ESP32
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

/* ── Cycle counter ──────────────────────────────────────────────── */
static inline __attribute__((always_inline))
uint32_t cc(void) { uint32_t c; __asm__ __volatile__("rsr.ccount %0":"=a"(c)); return c; }

/* ── Float DF I biquad ─────────────────────────────────────────── */
typedef struct {
    float b0,b1,b2,a1,a2;
    float x1,x2,y1,y2;
} BqFlt;

static inline float bqflt_tick(BqFlt *b, float x) {
    float y = b->b0*x + b->b1*b->x1 + b->b2*b->x2 - b->a1*b->y1 - b->a2*b->y2;
    b->x2 = b->x1; b->x1 = x;
    b->y2 = b->y1; b->y1 = y;
    return y;
}

/* ── Init with realistic EQ params (peaking at 1 kHz, Q=1.4, +6 dB) ─ */
static void init_bqs(BqFlt *f, int n) {
    float fs = 48000, fc = 1000, Q = 1.414f, gain_db = 6.0f;
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
    float fb0 = b0 * norm, fb1 = b1 * norm, fb2 = b2 * norm;
    float fa1 = a1 * norm, fa2 = a2 * norm;
    for (int i = 0; i < n; i++) {
        f[i].b0 = fb0; f[i].b1 = fb1; f[i].b2 = fb2;
        f[i].a1 = fa1; f[i].a2 = fa2;
        f[i].x1 = f[i].x2 = f[i].y1 = f[i].y2 = 0;
    }
}

#define N_BQ    6
#define N_SAMP  256
#define ITERS   200
#define WARMUP  20

void app_main(void) {
    ESP_LOGI("bench", "=== Float biquad benchmark, %d biquads, 240 MHz ===", N_BQ);

    BqFlt  bqflt[N_BQ];
    init_bqs(bqflt, N_BQ);

    float f_in[N_SAMP], f_out[N_SAMP];
    srand(42);
    for (int i = 0; i < N_SAMP; i++) {
        f_in[i] = (float)rand()/(float)RAND_MAX * 2.0f - 1.0f;
    }

    /* Warmup */
    for (int w = 0; w < WARMUP; w++) {
        for (int i = 0; i < N_SAMP; i++) {
            float y = f_in[i];
            for (int j = 0; j < N_BQ; j++) y = bqflt_tick(&bqflt[j], y);
        }
    }

    /* Benchmark float cascade */
    {
        uint32_t best = UINT32_MAX;
        for (int iter = 0; iter < ITERS; iter++) {
            uint32_t t0 = cc();
            for (int i = 0; i < N_SAMP; i++) {
                float y = f_in[i];
                for (int j = 0; j < N_BQ; j++) y = bqflt_tick(&bqflt[j], y);
                f_out[i] = y;
            }
            uint32_t t1 = cc();
            if (t1 - t0 < best) best = t1 - t0;
        }
        float us = (float)best / 240.0f;
        float per_tick = (float)best / (N_SAMP * N_BQ);
        ESP_LOGI("bench", "Float %d-bq cascade: %6u cyc = %7.2f us  (%.1f cy/tick)",
                 N_BQ, best, us, per_tick);
    }

    ESP_LOGI("bench", "Done.");
    while (1) { vTaskDelay(pdMS_TO_TICKS(5000)); }
}
