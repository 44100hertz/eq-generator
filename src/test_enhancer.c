/**
 * test_enhancer.c — Verify C enhancer implementation (float)
 *
 * Compile:
 *   make test   (from src/)
 *
 * Tests:
 *   1. DC gain of a single biquad
 *   2. Impulse response of Butterworth HP cascade
 *   3. Chebyshev T2/T3 correctness
 *   4. Full enhancer pipeline with synthetic bass tone
 *   5. CPU budget estimate
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <assert.h>

#include "biquad.h"
#include "envelope.h"
#include "enhancer.h"
#include "dc_blocker.h"

#define FS 44100.0f

/* ── Test 1: Biquad DC gain ───────────────────────────────────────── */

static int test_biquad_dc_gain(void) {
    printf("Test 1: Biquad DC gain...\n");

    float hp_coeffs[5];
    bass_design_butter_hp(25.0f, FS, hp_coeffs);

    Biquad bq;
    biquad_init(&bq, hp_coeffs);

    /* Feed constant DC = 0.5 for 200k samples */
    float out = 0.0f;
    for (int i = 0; i < 200000; i++) {
        out = biquad_tick(&bq, 0.5f);
    }

    printf("  HP alone after 200k DC: %.6f (raw — expect non-zero)\n", out);

    /* Now test DC blocker + HP cascade */
    DCBlocker dcb;
    float dc_R = expf(-2.0f * (float)M_PI * 5.0f / FS);
    dc_blocker_init(&dcb, dc_R);
    biquad_reset(&bq);

    float chain_out = 0.0f;
    for (int i = 0; i < 200000; i++) {
        chain_out = biquad_tick(&bq, dc_blocker_tick(&dcb, 0.5f));
    }
    printf("  DC blocker + HP after 200k: %.6f (expect < 0.01)\n", chain_out);

    int pass = (fabsf(chain_out) < 0.01f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 2: Chebyshev T2/T3 ──────────────────────────────────────── */

static int test_chebyshev(void) {
    printf("Test 2: Chebyshev polynomials...\n");
    int failures = 0;

    /* Chebyshev T2(0.5) = 2*(0.5)² - 1 = -0.5 */
    {
        float t2 = cheb_t2(0.5f);
        printf("  T2(0.5): expected -0.5000, got %.4f\n", t2);
        if (fabsf(t2 - (-0.5f)) > 0.001f) failures++;
    }

    /* Chebyshev T2(1.0) = 2 - 1 = 1.0 */
    {
        float t2 = cheb_t2(1.0f);
        printf("  T2(1.0): expected 1.0000, got %.4f\n", t2);
        if (fabsf(t2 - 1.0f) > 0.01f) failures++;
    }

    /* Chebyshev T3(0.5) = 4*(0.5)³ - 3*0.5 = -1.0 */
    {
        float t3 = cheb_t3(0.5f);
        printf("  T3(0.5): expected -1.0000, got %.4f\n", t3);
        if (fabsf(t3 - (-1.0f)) > 0.002f) failures++;
    }

    printf("  %s\n\n", failures == 0 ? "PASS" : "FAIL");
    return failures;
}

/* ── Test 3: Envelope follower ────────────────────────────────────── */

static int test_envelope(void) {
    printf("Test 3: Envelope follower...\n");

    float release_coeff = expf(-1.0f / (FS * 0.2f));

    Env env;
    env_init(&env, release_coeff);

    /* Feed a sine wave at 50 Hz, amplitude 0.5 */
    float max_peak = 0.0f;
    for (int i = 0; i < 44100; i++) {
        float t = (float)i / FS;
        float sin_val = 0.5f * sinf(2.0f * (float)M_PI * 50.0f * t);
        float peak = env_tick(&env, sin_val);
        if (peak > max_peak) max_peak = peak;
    }

    printf("  Peak of 0.5 sine: %.4f (expect ~0.5)\n", max_peak);

    int pass = (fabsf(max_peak - 0.5f) < 0.05f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 4: Full enhancer pipeline ───────────────────────────────── */

static int test_enhancer_pipeline(void) {
    printf("Test 4: Full enhancer pipeline...\n");

    float eq_coeffs[5];
    bass_design_butter_hp(25.0f, FS, eq_coeffs);

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         60.0f,   /* cutoff_hz */
                         0.33f,   /* h2_amp */
                         0.33f,   /* h3_amp */
                         0.2f,    /* release_secs */
                         FS,
                         1.0f,    /* push_gain */
                         1.0f,    /* pre_gain */
                         1,       /* eq_n_biquads */
                         eq_coeffs);

    Biquad eq_left[1], eq_right[1];
    BassEnhancer enh;
    BassEnhancer_init(&enh, &cfg, eq_left, eq_right);

    /* Feed a 40 Hz sine wave (below cutoff) */
    float rms_in = 0.0f, rms_out = 0.0f;
    int n_samples = 44100 * 2;  /* 2 seconds */
    float peak_out = 0.0f;

    for (int i = 0; i < n_samples; i++) {
        float t = (float)i / FS;
        float x = 0.25f * sinf(2.0f * (float)M_PI * 40.0f * t);

        float left  = x;
        float right = x;

        BassEnhancer_process_stereo(&enh, &left, &right);

        rms_in  += x * x;
        rms_out += left * left;

        if (fabsf(left) > peak_out) peak_out = fabsf(left);
    }

    rms_in  = sqrtf(rms_in / n_samples);
    rms_out = sqrtf(rms_out / n_samples);

    printf("  Input RMS:  %.4f\n", rms_in);
    printf("  Output RMS: %.4f\n", rms_out);
    printf("  Gain:       %.2f dB\n", 20.0f * log10f(rms_out / rms_in + 1e-12f));
    printf("  Peak:       %.4f (no clipping if < 1.0)\n", peak_out);

    int pass = (rms_out > 1e-4f && rms_out < 1.0f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 5: CPU budget estimate ──────────────────────────────────── */

static int test_cpu_budget(void) {
    printf("Test 5: CPU budget estimate...\n");

    float eq_coeffs[5];
    bass_design_butter_hp(40.0f, FS, eq_coeffs);

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg, 60.0f, 0.33f, 0.33f, 0.2f, FS, 1.0f, 1.0f, 1, eq_coeffs);

    Biquad eq_left[1], eq_right[1];
    BassEnhancer enh;
    BassEnhancer_init(&enh, &cfg, eq_left, eq_right);

    float left = 0.0f, right = 0.0f;

    clock_t start = clock();
    int n_iter = 100000;
    for (int i = 0; i < n_iter; i++) {
        left = 0.1f;
        right = 0.1f;
        BassEnhancer_process_stereo(&enh, &left, &right);
    }
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double us_per_frame = elapsed * 1e6 / n_iter;
    double samples_per_sec = n_iter / elapsed;
    double cpu_pct = elapsed * FS / n_iter * 100.0;

    printf("  Processed %d stereo frames in %.3f s\n", n_iter, elapsed);
    printf("  %.1f µs per stereo frame\n", us_per_frame);
    printf("  %.0f samples/sec (%.1f× real-time at %d Hz)\n",
           samples_per_sec, samples_per_sec / FS, (int)FS);
    printf("  CPU usage: %.1f%%\n", cpu_pct);

    printf("  PASS (benchmark)\n\n");
    return 0;
}

/* ── Main ──────────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;

    printf("========================================\n");
    printf("  Bass Enhancer C Implementation Tests (float)\n");
    printf("========================================\n\n");

    failures += test_biquad_dc_gain();
    failures += test_chebyshev();
    failures += test_envelope();
    failures += test_enhancer_pipeline();
    failures += test_cpu_budget();

    printf("========================================\n");
    printf("  %d test(s) failed\n", failures);
    printf("========================================\n");

    return failures;
}
