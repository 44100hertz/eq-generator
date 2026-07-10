/**
 * bench_sweep.c — Sweep FFT sizes and biquad counts for CPU benchmarking.
 *
 * Build: gcc -O2 -Wall -o bench_sweep bench_sweep.c fft_eq.c enhancer.c \
 *        -I. -lm -DFFT_EQ_N=<n>
 *
 * Usage: ./bench_sweep <fft_n> <n_biquads> <n_frames>
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "biquad_q28.h"
#include "fft_eq.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Q16 65536
#define FRAME_SAMPLES (FFT_EQ_HOP)

static int f2q16(float x) {
    if (x >= 0.9999f) return 0x7FFFFFFF;
    if (x <= -1.0f)   return (int32_t)0x80000000;
    return (int32_t)(x * Q16 + (x >= 0.0f ? 0.5f : -0.5f));
}

static float randf(void) { return (float)rand() / (float)RAND_MAX; }

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void design_bq(int32_t coeffs[5], float fc, float fs) {
    float omega = tanf((float)M_PI * fc / fs);
    float c = 1.0f + 1.41421356237f * omega + omega * omega;
    float b0 = 1.0f / c;
    float b1 = -2.0f / c;
    float b2 = 1.0f / c;
    float a1 = (2.0f * (omega * omega - 1.0f)) / c;
    float a2 = (1.0f - 1.41421356237f * omega + omega * omega) / c;
    float q28 = 268435456.0f;
    coeffs[0] = (int32_t)(b0 * q28 + (b0 >= 0.0f ? 0.5f : -0.5f));
    coeffs[1] = (int32_t)(b1 * q28 + (b1 >= 0.0f ? 0.5f : -0.5f));
    coeffs[2] = (int32_t)(b2 * q28 + (b2 >= 0.0f ? 0.5f : -0.5f));
    coeffs[3] = (int32_t)(a1 * q28 + (a1 >= 0.0f ? 0.5f : -0.5f));
    coeffs[4] = (int32_t)(a2 * q28 + (a2 >= 0.0f ? 0.5f : -0.5f));
}

int main(int argc, char *argv[]) {
    int n_bq = (argc > 1) ? atoi(argv[1]) : 6;
    int n_frames = (argc > 2) ? atoi(argv[2]) : 2000;

    float fs = 48000.0f;

    /* Setup biquad cascade */
    BiquadQ28 *bqs_l = NULL, *bqs_r = NULL;
    int32_t *bq_coeffs = NULL;
    float bq_freqs[40];
    if (n_bq > 0) {
        bqs_l = calloc((size_t)n_bq, sizeof(BiquadQ28));
        bqs_r = calloc((size_t)n_bq, sizeof(BiquadQ28));
        bq_coeffs = calloc((size_t)(n_bq * 5), sizeof(int32_t));
        for (int i = 0; i < n_bq; i++) {
            bq_freqs[i] = 20.0f * powf(20000.0f / 20.0f, (float)i / (float)n_bq);
            design_bq(&bq_coeffs[i * 5], bq_freqs[i], fs);
            BiquadQ28_init(&bqs_l[i], &bq_coeffs[i * 5]);
            BiquadQ28_init(&bqs_r[i], &bq_coeffs[i * 5]);
        }
    }

    /* Setup FFT EQ (random gains) */
    float fft_gains[FFT_EQ_BINS];
    for (int i = 0; i < FFT_EQ_BINS; i++) {
        fft_gains[i] = 0.5f + randf() * 1.5f;
    }

    FftEq *fft_eq = fft_eq_create(fft_gains);
    if (!fft_eq) { fprintf(stderr, "FFT alloc failed\n"); return 1; }

    /* Generate pseudo-pink test signal */
    float *input_l  = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *input_r  = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *mid_l    = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *mid_r    = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_l = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_r = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float pink_l = 0.0f, pink_r = 0.0f;
    srand(42);
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        pink_l += (randf() - 0.5f) * 0.15f;
        pink_r += (randf() - 0.5f) * 0.15f;
        pink_l *= 0.999f;
        pink_r *= 0.999f;
        input_l[i] = pink_l;
        input_r[i] = pink_r;
    }

    /* ── Benchmark: combined biquads → FFT EQ ── */
    fft_eq_reset(fft_eq);
    double t0 = now_sec();
    for (int f = 0; f < n_frames; f++) {
        if (n_bq > 0) {
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                int32_t ql = f2q16(input_l[i]);
                int32_t qr = f2q16(input_r[i]);
                ql = BiquadQ28_cascade(bqs_l, n_bq, ql);
                qr = BiquadQ28_cascade(bqs_r, n_bq, qr);
                mid_l[i] = (float)ql / Q16;
                mid_r[i] = (float)qr / Q16;
            }
        } else {
            memcpy(mid_l, input_l, FRAME_SAMPLES * sizeof(float));
            memcpy(mid_r, input_r, FRAME_SAMPLES * sizeof(float));
        }
        fft_eq_process_frame(fft_eq, mid_l, mid_r, output_l, output_r);
    }
    double elapsed = now_sec() - t0;
    double us_per_frame = elapsed * 1e6 / (double)n_frames;
    double frame_budget_us = 1e6 / fs * (double)FRAME_SAMPLES;

    /* Also measure biquads alone and FFT alone for breakdown */
    fft_eq_reset(fft_eq);

    double t_bq = 0.0, t_fft = 0.0;
    if (n_bq > 0) {
        t0 = now_sec();
        for (int f = 0; f < n_frames; f++) {
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                int32_t ql = f2q16(input_l[i]);
                int32_t qr = f2q16(input_r[i]);
                ql = BiquadQ28_cascade(bqs_l, n_bq, ql);
                qr = BiquadQ28_cascade(bqs_r, n_bq, qr);
            }
        }
        t_bq = (now_sec() - t0) / (double)n_frames;
    }

    t0 = now_sec();
    for (int f = 0; f < n_frames; f++) {
        fft_eq_process_frame(fft_eq, input_l, input_r, output_l, output_r);
    }
    t_fft = (now_sec() - t0) / (double)n_frames;

    double us_bq_frame = t_bq * 1e6;
    double us_fft_frame = t_fft * 1e6;

    /* Print JSON result for easy parsing */
    printf("{\n");
    printf("  \"fft_n\": %d,\n", FFT_EQ_N);
    printf("  \"hop\": %d,\n", FFT_EQ_HOP);
    printf("  \"bins\": %d,\n", FFT_EQ_BINS);
    printf("  \"n_biquads\": %d,\n", n_bq);
    printf("  \"n_frames\": %d,\n", n_frames);
    printf("  \"total_ms\": %.6f,\n", elapsed * 1e3);
    printf("  \"us_per_frame\": %.2f,\n", us_per_frame);
    printf("  \"us_bq_frame\": %.2f,\n", us_bq_frame);
    printf("  \"us_fft_frame\": %.2f,\n", us_fft_frame);
    printf("  \"frame_budget_us\": %.2f,\n", frame_budget_us);
    printf("  \"cpu_pct_desktop\": %.2f,\n", us_per_frame / frame_budget_us * 100.0);
    printf("  \"desktop_speedup\": %.1f\n", frame_budget_us / us_per_frame);
    printf("}\n");

    fft_eq_destroy(fft_eq);
    free(input_l); free(input_r); free(mid_l); free(mid_r);
    free(output_l); free(output_r);
    free(bqs_l); free(bqs_r); free(bq_coeffs);
    return 0;
}
