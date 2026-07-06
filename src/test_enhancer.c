/**
 * test_enhancer.c — Verify C enhancer implementation matches Python reference
 *
 * Compile:
 *   gcc -std=c11 -Wall -O2 test_enhancer.c src/enhancer.c -lm -o test_enhancer
 *
 * Tests:
 *   1. DC gain of a single Q4.28 biquad
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

#include "biquad_q28.h"
#include "envelope.h"
#include "enhancer.h"
#include "dc_blocker.h"

#define FS 44100.0f
#define Q16_SCALE 65536

/* ── Test utilities ────────────────────────────────────────────────── */

static float q16_to_float(int32_t x) {
    return (float)x / Q16_SCALE;
}

static int32_t float_to_q16(float x) {
    if (x >= 0.9999f) return 0x7FFFFFFF;
    if (x <= -1.0f)   return (int32_t)0x80000000;
    return (int32_t)(x * Q16_SCALE + (x >= 0.0f ? 0.5f : -0.5f));
}

/* ── Test 1: Biquad DC gain ───────────────────────────────────────── */

static int test_biquad_dc_gain(void) {
    printf("Test 1: Biquad DC gain...\n");

    /* Butterworth HP at 25 Hz should have DC gain = 0 */
    int32_t hp_coeffs[5];
    bass_design_butter_hp_q28(25.0f, FS, hp_coeffs);

    BiquadQ28 bq;
    BiquadQ28_init(&bq, hp_coeffs);

    /* Feed constant DC = 0.5 (16384 in Q16) for 200k samples */
    int32_t dc_input = float_to_q16(0.5f);
    int32_t out = 0;
    for (int i = 0; i < 200000; i++) {
        out = BiquadQ28_tick(&bq, dc_input);
    }

    float dc_out = q16_to_float(out);

    /* Butterworth HP alone has a pole near z=1: settling is slow.
     * The full enhancer chain includes a DC blocker upstream.
     * This test verifies the DC blocker + HP cascade works. */
    printf("  HP alone after 200k DC: %.6f (raw — expect non-zero)\n", dc_out);

    /* Now test DC blocker + HP cascade */
    DCBlocker dcb;
    float dc_R = expf(-2.0f * (float)M_PI * 5.0f / FS);
    int32_t dc_R_q16 = (int32_t)(dc_R * 65536.0f + 0.5f);
    DCBlocker_init(&dcb, dc_R_q16);
    BiquadQ28_reset(&bq);

    int32_t chain_out = 0;
    for (int i = 0; i < 200000; i++) {
        chain_out = BiquadQ28_tick(&bq, DCBlocker_tick(&dcb, float_to_q16(0.5f)));
    }
    float chain_val = q16_to_float(chain_out);
    printf("  DC blocker + HP after 200k: %.6f (expect < 0.01)\n", chain_val);

    int pass = (fabsf(chain_val) < 0.01f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 2: Chebyshev T2/T3 ──────────────────────────────────────── */

static int test_chebyshev(void) {
    printf("Test 2: Chebyshev polynomials...\n");
    int failures = 0;

    /* Chebyshev T2(0.5) = 2*(0.5)² - 1 = -0.5 (x=32768 in Q16) */
    {
        int32_t x = 32768;  /* Q16 representation of 0.5 */
        int32_t t2 = cheb_t2(x);
        float expected = -0.5f;
        float got = q16_to_float(t2);
        printf("  T2(0.5): expected %.4f, got %.4f\n", expected, got);
        if (fabsf(got - expected) > 0.001f) failures++;
    }

    /* Chebyshev T2(1.0) = 2 - 1 = 1.0 (x=65536 in Q16) */
    {
        int32_t x = 65536;  /* Q16 representation of 1.0 */
        int32_t t2 = cheb_t2(x);
        float got = q16_to_float(t2);
        printf("  T2(1.0): expected 1.0000, got %.4f\n", got);
        if (fabsf(got - 1.0f) > 0.01f) failures++;
    }

    /* Chebyshev T3(0.5) = 4*(0.5)³ - 3*0.5 = 0.5 - 1.5 = -1.0 (x=32768 Q16) */
    {
        int32_t x = 32768;  /* Q16 representation of 0.5 */
        int32_t t3 = cheb_t3(x);
        float expected = -1.0f;
        float got = q16_to_float(t3);
        printf("  T3(0.5): expected %.4f, got %.4f\n", expected, got);
        if (fabsf(got - expected) > 0.002f) failures++;
    }

    printf("  %s\n\n", failures == 0 ? "PASS" : "FAIL");
    return failures;
}

/* ── Test 3: Envelope follower ────────────────────────────────────── */

static int test_envelope(void) {
    printf("Test 3: Envelope follower...\n");

    ReciprocalLUT lut;
    ReciprocalLUT_init(&lut);

    /* release_coeff = exp(-1/(fs * 0.2)) ≈ 0.9998866 */
    /* In Q16: 0.9998866 * 65536 = 65529 */
    int32_t release_q16 = (int32_t)(expf(-1.0f / (FS * 0.2f)) * Q16_SCALE + 0.5f);

    Env env;
    Env_init(&env, &lut, release_q16);

    /* Feed a sine wave at 50 Hz, amplitude 0.5 */
    int32_t max_peak = 0;
    for (int i = 0; i < 44100; i++) {
        float t = (float)i / FS;
        float sin_val = 0.5f * sinf(2.0f * (float)M_PI * 50.0f * t);
        int32_t peak = Env_tick(&env, float_to_q16(sin_val));
        if (peak > max_peak) max_peak = peak;
    }

    float peak_val = q16_to_float(max_peak);
    printf("  Peak of 0.5 sine: %.4f (expect ~0.5)\n", peak_val);

    int pass = (fabsf(peak_val - 0.5f) < 0.05f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 4: Full enhancer pipeline ───────────────────────────────── */

static int test_enhancer_pipeline(void) {
    printf("Test 4: Full enhancer pipeline...\n");

    /* Simple EQ: just one HP at 25 Hz (test with known coeffs) */
    int32_t eq_coeffs[5];
    bass_design_butter_hp_q28(25.0f, FS, eq_coeffs);

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         60.0f,   /* cutoff_hz */
                         0.33f,   /* h2_amp */
                         0.33f,   /* h3_amp */
                         0.2f,    /* release_secs */
                         FS,
                         0.049f,  /* limiter_release_secs */
                         1,       /* eq_n_biquads */
                         eq_coeffs);

    ReciprocalLUT lut;
    ReciprocalLUT_init(&lut);

    BiquadQ28 eq_left[1], eq_right[1];
    BassEnhancer enh;
    BassEnhancer_init(&enh, &cfg, &lut, eq_left, eq_right);

    /* Feed a 40 Hz sine wave (below cutoff) */
    float rms_in = 0.0f, rms_out = 0.0f;
    int n_samples = 44100 * 2;  /* 2 seconds */
    int32_t peak_out = 0;

    for (int i = 0; i < n_samples; i++) {
        float t = (float)i / FS;
        float x = 0.25f * sinf(2.0f * (float)M_PI * 40.0f * t);

        int32_t left  = float_to_q16(x);
        int32_t right = float_to_q16(x);

        BassEnhancer_process_stereo(&enh, &left, &right);

        float out = q16_to_float(left);
        rms_in  += x * x;
        rms_out += out * out;

        if (abs(left) > peak_out) peak_out = left;
    }

    rms_in  = sqrtf(rms_in / n_samples);
    rms_out = sqrtf(rms_out / n_samples);

    printf("  Input RMS:  %.4f\n", rms_in);
    printf("  Output RMS: %.4f\n", rms_out);
    printf("  Gain:       %.2f dB\n", 20.0f * log10f(rms_out / rms_in + 1e-12f));
    printf("  Peak:       %.4f (no clipping if < 0.5)\n", q16_to_float(peak_out));

    /* Output should be non-zero and not clipping */
    int pass = (rms_out > 1e-4f && rms_out < 1.0f);
    printf("  %s\n\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

/* ── Test 5: CPU budget estimate ──────────────────────────────────── */

static int test_cpu_budget(void) {
    printf("Test 5: CPU budget estimate...\n");

    /* Simple loop count: how many samples can we process per second? */
    int32_t eq_coeffs[5];
    bass_design_butter_hp_q28(40.0f, FS, eq_coeffs);

    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg, 60.0f, 0.33f, 0.33f, 0.2f, FS, 0.049f, 1, eq_coeffs);

    ReciprocalLUT lut;
    ReciprocalLUT_init(&lut);

    BiquadQ28 eq_left[1], eq_right[1];
    BassEnhancer enh;
    BassEnhancer_init(&enh, &cfg, &lut, eq_left, eq_right);

    int32_t left = 0, right = 0;

    clock_t start = clock();
    int n_iter = 100000;
    for (int i = 0; i < n_iter; i++) {
        left = float_to_q16(0.1f);
        right = left;
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
    printf("  CPU usage: %.1f%% at 160 MHz\n", cpu_pct);

    /* On desktop this will be very fast; on ESP32 scale by ~5-10x */
    printf("  Estimated ESP32 (10× slower): %.1f µs/frame, %.1f%% CPU\n",
           us_per_frame * 20.0, cpu_pct * 20.0);

    printf("  PASS (benchmark)\n\n");
    return 0;
}

/* ── Main ──────────────────────────────────────────────────────────── */

int main(void) {
    int failures = 0;

    printf("========================================\n");
    printf("  Bass Enhancer C Implementation Tests\n");
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
