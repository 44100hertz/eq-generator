/**
 * bench_hybrid.c — Performance benchmark: 6 biquads + 256-pt FFT EQ (float)
 *
 * Build:
 *   make bench  (from src/)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "biquad.h"
#include "fft_eq.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS            48000.0f
#define N_BIQUADS     6
#define N_FRAMES      2000
#define FRAME_SAMPLES 128

/* ── Timing helpers ────────────────────────────────────────────────── */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ── Generate random float in [0, 1] ──────────────────────────────── */

static float randf(void) {
    return (float)rand() / (float)RAND_MAX;
}

/* ── Make 6 biquad coeffs (HP at various freqs) ───────────────────── */

static void design_bq(float coeffs[5], float fc, float fs) {
    float omega = tanf((float)M_PI * fc / fs);
    float c = 1.0f + 1.41421356237f * omega + omega * omega;
    coeffs[0] = 1.0f / c;
    coeffs[1] = -2.0f / c;
    coeffs[2] = 1.0f / c;
    coeffs[3] = (2.0f * (omega * omega - 1.0f)) / c;
    coeffs[4] = (1.0f - 1.41421356237f * omega + omega * omega) / c;
}

/* ── Main ──────────────────────────────────────────────────────────── */

int main(void) {
    printf("=== Hybrid EQ Performance Benchmark (float) ===\n");
    printf("Sample rate: %.0f Hz\n", FS);
    printf("Frames: %d × %d samples = %d total\n\n",
           N_FRAMES, FRAME_SAMPLES, N_FRAMES * FRAME_SAMPLES);

    /* ── Setup biquad cascade ────────────────────────────────────── */
    float bq_freqs[N_BIQUADS] = {40, 100, 250, 500, 1200, 3000};
    float bq_coeffs[N_BIQUADS * 5];
    Biquad bqs_l[N_BIQUADS], bqs_r[N_BIQUADS];

    for (int i = 0; i < N_BIQUADS; i++) {
        design_bq(&bq_coeffs[i * 5], bq_freqs[i], FS);
        biquad_init(&bqs_l[i], &bq_coeffs[i * 5]);
        biquad_init(&bqs_r[i], &bq_coeffs[i * 5]);
    }

    /* ── Setup FFT EQ (random gains) ─────────────────────────────── */
    float fft_gains[FFT_EQ_BINS];
    for (int i = 0; i < FFT_EQ_BINS; i++) {
        fft_gains[i] = 0.5f + randf() * 1.5f;
    }

    FftEq *fft_eq = fft_eq_create(fft_gains);
    if (!fft_eq) {
        fprintf(stderr, "Failed to create FftEq\n");
        return 1;
    }

    /* ── Generate test signal (pink-ish noise) ──────────────────── */
    srand(42);
    float *input_l  = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *input_r  = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_l = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_r = calloc((size_t)FRAME_SAMPLES, sizeof(float));
    if (!input_l || !input_r || !output_l || !output_r) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }

    float pink_l = 0.0f, pink_r = 0.0f;
    for (int i = 0; i < FRAME_SAMPLES; i++) {
        pink_l += (randf() - 0.5f) * 0.15f;
        pink_r += (randf() - 0.5f) * 0.15f;
        pink_l *= 0.999f;
        pink_r *= 0.999f;
        input_l[i] = pink_l;
        input_r[i] = pink_r;
    }

    printf("Signal: pseudo-pink noise, amplitude ~0.3 RMS\n\n");

    /* ── Benchmark 1: Biquads only ───────────────────────────────── */
    {
        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                output_l[i] = biquad_cascade(bqs_l, N_BIQUADS, input_l[i]);
                output_r[i] = biquad_cascade(bqs_r, N_BIQUADS, input_r[i]);
            }
        }
        double elapsed = now_sec() - t0;
        double per_frame_us = elapsed * 1e6 / (double)N_FRAMES;
        double per_sample_ns = elapsed * 1e9 / (double)(N_FRAMES * FRAME_SAMPLES * 2);

        printf("── %d-biquad cascade ──\n", N_BIQUADS);
        printf("  Total:       %.3f ms\n", elapsed * 1e3);
        printf("  Per frame:   %.1f µs  (%d stereo samples)\n", per_frame_us, FRAME_SAMPLES);
        printf("  Per sample:  %.0f ns\n", per_sample_ns);
        printf("  Throughput:  %.1f× real-time\n\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);
    }

    /* ── Benchmark 2: FFT EQ only ────────────────────────────────── */
    {
        fft_eq_reset(fft_eq);

        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            fft_eq_process_frame(fft_eq, input_l, input_r, output_l, output_r);
        }
        double elapsed = now_sec() - t0;
        double per_frame_us = elapsed * 1e6 / (double)N_FRAMES;

        printf("── FFT EQ (%d-pt, %d-hop, Hann) ──\n", FFT_EQ_N, FFT_EQ_HOP);
        printf("  Total:       %.3f ms\n", elapsed * 1e3);
        printf("  Per frame:   %.1f µs  (%d stereo samples)\n", per_frame_us, FRAME_SAMPLES);
        printf("  Throughput:  %.1f× real-time\n\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);
    }

    /* ── Benchmark 3: Combined (FFT → biquads) ───────────────────── */
    {
        fft_eq_reset(fft_eq);
        for (int i = 0; i < N_BIQUADS; i++) {
            biquad_reset(&bqs_l[i]);
            biquad_reset(&bqs_r[i]);
        }

        float *mid_l = calloc((size_t)FRAME_SAMPLES, sizeof(float));
        float *mid_r = calloc((size_t)FRAME_SAMPLES, sizeof(float));

        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            /* FFT pass */
            fft_eq_process_frame(fft_eq, input_l, input_r, mid_l, mid_r);
            /* IIR pass */
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                output_l[i] = biquad_cascade(bqs_l, N_BIQUADS, mid_l[i]);
                output_r[i] = biquad_cascade(bqs_r, N_BIQUADS, mid_r[i]);
            }
        }
        double elapsed = now_sec() - t0;
        double per_frame_us = elapsed * 1e6 / (double)N_FRAMES;

        printf("── Combined: FFT EQ → %d biquads ──\n", N_BIQUADS);
        printf("  Total:       %.3f ms\n", elapsed * 1e3);
        printf("  Per frame:   %.1f µs  (%d stereo samples)\n", per_frame_us, FRAME_SAMPLES);
        printf("  Throughput:  %.1f× real-time\n\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);

        free(mid_l);
        free(mid_r);
    }

    /* ── Summary ──────────────────────────────────────────────────── */
    {
        fft_eq_reset(fft_eq);
        for (int i = 0; i < N_BIQUADS; i++) {
            biquad_reset(&bqs_l[i]);
            biquad_reset(&bqs_r[i]);
        }

        /* Biquad time */
        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                output_l[i] = biquad_cascade(bqs_l, N_BIQUADS, input_l[i]);
                output_r[i] = biquad_cascade(bqs_r, N_BIQUADS, input_r[i]);
            }
        }
        double t_bq = (now_sec() - t0) / (double)N_FRAMES;

        /* FFT time */
        t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            fft_eq_process_frame(fft_eq, input_l, input_r, output_l, output_r);
        }
        double t_fft = (now_sec() - t0) / (double)N_FRAMES;

        double frame_budget_us = 1e6 / FS * (double)FRAME_SAMPLES;
        double us_per_frame = (t_bq + t_fft) * 1e6;

        printf("── Summary ──\n");
        printf("  Biquads (%d):  %.1f µs/frame\n", N_BIQUADS, t_bq * 1e6);
        printf("  FFT EQ:        %.1f µs/frame\n", t_fft * 1e6);
        printf("  Combined:      %.1f µs/frame\n", us_per_frame);
        printf("  Budget:        %.0f µs/frame (128 samples @ %.0f Hz)\n",
               frame_budget_us, FS);
        printf("  CPU (desktop): %.1f%%\n\n",
               us_per_frame / frame_budget_us * 100.0);
    }

    fft_eq_destroy(fft_eq);
    free(input_l);
    free(input_r);
    free(output_l);
    free(output_r);

    return 0;
}
