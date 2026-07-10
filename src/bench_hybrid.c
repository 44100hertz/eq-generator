/**
 * bench_hybrid.c — Performance benchmark: 6 biquads + 256-pt FFT EQ
 *
 * Measures wall-clock time for each processing stage independently,
 * then combined. Uses random FFT gains — no real EQ curve needed.
 *
 * Build:
 *   gcc -O2 -Wall -o bench_hybrid bench_hybrid.c fft_eq.c enhancer.c -I. -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "biquad_q28.h"
#include "fft_eq.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define FS            48000.0f
#define N_BIQUADS     6
#define N_FRAMES      2000         /* frames to process per test */
#define FRAME_SAMPLES 128          /* FFT hop size */

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

/* ── Q16 conversion helpers ────────────────────────────────────────── */

#define Q16 65536

static int32_t f2q16(float x) {
    if (x >= 0.9999f) return 0x7FFFFFFF;
    if (x <= -1.0f)   return (int32_t)0x80000000;
    return (int32_t)(x * Q16 + (x >= 0.0f ? 0.5f : -0.5f));
}

static float q162f(int32_t x) {
    return (float)x / Q16;
}

/* ── Make 6 biquad coeffs (HP at various freqs, stable) ───────────── */

/* Borrowed from enhancer.c — designs Butterworth HP in Q4.28 */
static void design_bq(int32_t coeffs[5], float fc, float fs) {
    float omega = tanf((float)M_PI * fc / fs);
    float c = 1.0f + 1.41421356237f * omega + omega * omega;
    float b0 = 1.0f / c;
    float b1 = -2.0f / c;
    float b2 = 1.0f / c;
    float a1 = (2.0f * (omega * omega - 1.0f)) / c;
    float a2 = (1.0f - 1.41421356237f * omega + omega * omega) / c;

    float q28 = 268435456.0;
    coeffs[0] = (int32_t)(b0 * q28 + (b0 >= 0.0f ? 0.5f : -0.5f));
    coeffs[1] = (int32_t)(b1 * q28 + (b1 >= 0.0f ? 0.5f : -0.5f));
    coeffs[2] = (int32_t)(b2 * q28 + (b2 >= 0.0f ? 0.5f : -0.5f));
    coeffs[3] = (int32_t)(a1 * q28 + (a1 >= 0.0f ? 0.5f : -0.5f));
    coeffs[4] = (int32_t)(a2 * q28 + (a2 >= 0.0f ? 0.5f : -0.5f));
}

/* ── Main ──────────────────────────────────────────────────────────── */

int main(void) {
    printf("=== Hybrid EQ Performance Benchmark ===\n");
    printf("Sample rate: %.0f Hz\n", FS);
    printf("Frames: %d × %d samples = %d total\n\n",
           N_FRAMES, FRAME_SAMPLES, N_FRAMES * FRAME_SAMPLES);

    /* ── Setup biquad cascade (Q28) ──────────────────────────────── */
    float bq_freqs[N_BIQUADS] = {40, 100, 250, 500, 1200, 3000};
    int32_t bq_coeffs[N_BIQUADS * 5];
    BiquadQ28 bqs_l[N_BIQUADS], bqs_r[N_BIQUADS];

    for (int i = 0; i < N_BIQUADS; i++) {
        design_bq(&bq_coeffs[i * 5], bq_freqs[i], FS);
        BiquadQ28_init(&bqs_l[i], &bq_coeffs[i * 5]);
        BiquadQ28_init(&bqs_r[i], &bq_coeffs[i * 5]);
    }

    /* ── Setup FFT EQ (random gains) ─────────────────────────────── */
    float fft_gains[FFT_EQ_BINS];  /* 129 entries: DC through Nyquist */
    for (int i = 0; i < FFT_EQ_BINS; i++) {
        fft_gains[i] = 0.5f + randf() * 1.5f;  /* 0.5 .. 2.0 */
    }

    FftEq *fft_eq = fft_eq_create(fft_gains);
    if (!fft_eq) {
        fprintf(stderr, "Failed to create FftEq\n");
        return 1;
    }

    /* ── Generate test signal (pink-ish noise for realism) ───────── */
    srand(42);  /* deterministic */
    float *input_l  = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *input_r  = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_l = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));
    float *output_r = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));
    if (!input_l || !input_r || !output_l || !output_r) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }

    /* Pink-ish: integrate white noise (simple 1/f approx) */
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
                int32_t ql = f2q16(input_l[i]);
                int32_t qr = f2q16(input_r[i]);
                ql = BiquadQ28_cascade(bqs_l, N_BIQUADS, ql);
                qr = BiquadQ28_cascade(bqs_r, N_BIQUADS, qr);
                output_l[i] = q162f(ql);
                output_r[i] = q162f(qr);
            }
        }
        double elapsed = now_sec() - t0;
        double per_frame_us = elapsed * 1e6 / (double)N_FRAMES;
        double per_sample_ns = elapsed * 1e9 / (double)(N_FRAMES * FRAME_SAMPLES * 2);

        printf("── %d-biquad cascade ──\n", N_BIQUADS);
        printf("  Total:       %.3f ms\n", elapsed * 1e3);
        printf("  Per frame:   %.1f µs  (%d stereo samples)\n", per_frame_us, FRAME_SAMPLES);
        printf("  Per sample:  %.0f ns\n", per_sample_ns);
        printf("  Throughput:  %.1f× real-time\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);
        printf("\n");
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
        printf("  Throughput:  %.1f× real-time\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);
        printf("\n");
    }

    /* ── Benchmark 3: Combined (biquads → FFT) ───────────────────── */
    {
        fft_eq_reset(fft_eq);
        /* Reset biquad states */
        for (int i = 0; i < N_BIQUADS; i++) {
            BiquadQ28_reset(&bqs_l[i]);
            BiquadQ28_reset(&bqs_r[i]);
        }

        float *mid_l = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));
        float *mid_r = (float *)calloc((size_t)FRAME_SAMPLES, sizeof(float));

        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            /* IIR pass */
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                int32_t ql = f2q16(input_l[i]);
                int32_t qr = f2q16(input_r[i]);
                ql = BiquadQ28_cascade(bqs_l, N_BIQUADS, ql);
                qr = BiquadQ28_cascade(bqs_r, N_BIQUADS, qr);
                mid_l[i] = q162f(ql);
                mid_r[i] = q162f(qr);
            }
            /* FFT pass */
            fft_eq_process_frame(fft_eq, mid_l, mid_r, output_l, output_r);
        }
        double elapsed = now_sec() - t0;
        double per_frame_us = elapsed * 1e6 / (double)N_FRAMES;

        printf("── Combined: %d biquads → FFT EQ ──\n", N_BIQUADS);
        printf("  Total:       %.3f ms\n", elapsed * 1e3);
        printf("  Per frame:   %.1f µs  (%d stereo samples)\n", per_frame_us, FRAME_SAMPLES);
        printf("  Throughput:  %.1f× real-time\n",
               (double)(N_FRAMES * FRAME_SAMPLES) / FS / elapsed);
        printf("\n");

        free(mid_l);
        free(mid_r);
    }

    /* ── Summary ──────────────────────────────────────────────────── */
    /* Re-measure each independently for fair breakdown */
    {
        fft_eq_reset(fft_eq);
        for (int i = 0; i < N_BIQUADS; i++) {
            BiquadQ28_reset(&bqs_l[i]);
            BiquadQ28_reset(&bqs_r[i]);
        }

        /* Biquad time */
        double t0 = now_sec();
        for (int f = 0; f < N_FRAMES; f++) {
            for (int i = 0; i < FRAME_SAMPLES; i++) {
                int32_t ql = f2q16(input_l[i]);
                int32_t qr = f2q16(input_r[i]);
                ql = BiquadQ28_cascade(bqs_l, N_BIQUADS, ql);
                qr = BiquadQ28_cascade(bqs_r, N_BIQUADS, qr);
                output_l[i] = q162f(ql);
                output_r[i] = q162f(qr);
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
        /* Desktop is ~gigahertz, ESP32 is 240 MHz. Ratio depends on IPC.
           Rough estimate: desktop IPC ~2-3, Xtensa LX6 ~1. 240/4000*2 ≈ 12%.
           Conservative: 15× desktop slowdown. */
        double cpu_pct_desktop = us_per_frame / frame_budget_us * 100.0;
        double esp32_factor = 15.0;

        printf("── Summary ──\n");
        printf("  Biquads (%d):  %.1f µs/frame\n", N_BIQUADS, t_bq * 1e6);
        printf("  FFT EQ:        %.1f µs/frame\n", t_fft * 1e6);
        printf("  Combined:      %.1f µs/frame\n", us_per_frame);
        printf("  Budget:        %.0f µs/frame (128 samples @ %.0f Hz)\n",
               1e6 / FS * (double)FRAME_SAMPLES, FS);
        printf("  CPU @ 240 MHz: ~%.1f%% (est. %.0f× desktop slowdown)\n",
               cpu_pct_desktop * esp32_factor, esp32_factor);
        printf("\n");
    }

    /* Cleanup */
    fft_eq_destroy(fft_eq);
    free(input_l);
    free(input_r);
    free(output_l);
    free(output_r);

    return 0;
}
