/**
 * fft_eq.c — Overlap-add FFT EQ implementation
 *
 * Uses a straightforward radix-2 DIT complex FFT (configurable via EQGEN_FFT_N).
 * For real input, imag parts are zeroed — this wastes ~2x vs a proper
 * real FFT, but gives a conservative performance bound.
 */

#include "fft_eq.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* ── Radix-2 DIT complex FFT (in-place, interleaved re/im) ─────────── */

static void fft_cpx(float *data, int n, int inverse) {
    /* Bit-reversal permutation */
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (i < j) {
            float tr = data[2*i], ti = data[2*i+1];
            data[2*i]   = data[2*j];
            data[2*i+1] = data[2*j+1];
            data[2*j]   = tr;
            data[2*j+1] = ti;
        }
        int m = n >> 1;
        while (m >= 1 && j >= m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    /* Butterfly stages */
    float dir = inverse ? 1.0f : -1.0f;
    for (int len = 2; len <= n; len <<= 1) {
        float ang = 2.0f * (float)M_PI / (float)len * dir;
        float w_re = cosf(ang);
        float w_im = sinf(ang);
        for (int i = 0; i < n; i += len) {
            float cur_re = 1.0f, cur_im = 0.0f;
            int half = len >> 1;
            for (int k = 0; k < half; k++) {
                int a = i + k;
                int b = a + half;
                /* t = cur * data[b] */
                float t_re = cur_re * data[2*b] - cur_im * data[2*b+1];
                float t_im = cur_re * data[2*b+1] + cur_im * data[2*b];
                data[2*b]   = data[2*a] - t_re;
                data[2*b+1] = data[2*a+1] - t_im;
                data[2*a]   += t_re;
                data[2*a+1] += t_im;
                /* Twiddle factor update: cur *= w */
                float n_re = cur_re * w_re - cur_im * w_im;
                float n_im = cur_re * w_im + cur_im * w_re;
                cur_re = n_re;
                cur_im = n_im;
            }
        }
    }

    /* Scale for inverse transform */
    if (inverse) {
        float scale = 1.0f / (float)n;
        for (int i = 0; i < 2 * n; i++) {
            data[i] *= scale;
        }
    }
}

/* ── Hann window ───────────────────────────────────────────────────── */

static void make_hann_window(float *w, int n) {
    for (int i = 0; i < n; i++) {
        w[i] = 0.5f * (1.0f - cosf(2.0f * (float)M_PI * (float)i / (float)n));
    }
}

/* ── Create / destroy / reset ──────────────────────────────────────── */

FftEq *fft_eq_create(const float *gains) {
    FftEq *eq = (FftEq *)calloc(1, sizeof(FftEq));
    if (!eq) return NULL;

    eq->n   = FFT_EQ_N;
    eq->hop = FFT_EQ_HOP;

    eq->window     = (float *)calloc((size_t)eq->n, sizeof(float));
    eq->gains      = (float *)calloc((size_t)(eq->n / 2 + 1), sizeof(float));
    eq->overlap_l  = (float *)calloc((size_t)eq->n, sizeof(float));
    eq->overlap_r  = (float *)calloc((size_t)eq->n, sizeof(float));
    eq->olap_add_l = (float *)calloc((size_t)eq->hop, sizeof(float));
    eq->olap_add_r = (float *)calloc((size_t)eq->hop, sizeof(float));
    eq->fft_work   = (float *)calloc((size_t)(2 * eq->n), sizeof(float));

    if (!eq->window || !eq->gains || !eq->overlap_l || !eq->overlap_r ||
        !eq->olap_add_l || !eq->olap_add_r || !eq->fft_work) {
        fft_eq_destroy(eq);
        return NULL;
    }

    make_hann_window(eq->window, eq->n);
    memcpy(eq->gains, gains, (size_t)(eq->n / 2 + 1) * sizeof(float));

    return eq;
}

void fft_eq_destroy(FftEq *eq) {
    if (!eq) return;
    free(eq->window);
    free(eq->gains);
    free(eq->overlap_l);
    free(eq->overlap_r);
    free(eq->olap_add_l);
    free(eq->olap_add_r);
    free(eq->fft_work);
    free(eq);
}

void fft_eq_reset(FftEq *eq) {
    memset(eq->overlap_l,  0, (size_t)eq->n   * sizeof(float));
    memset(eq->overlap_r,  0, (size_t)eq->n   * sizeof(float));
    memset(eq->olap_add_l, 0, (size_t)eq->hop * sizeof(float));
    memset(eq->olap_add_r, 0, (size_t)eq->hop * sizeof(float));
}

/* ── Stereo processing (L↦real, R↦imag packed into one complex FFT) ── */

void fft_eq_process_frame(FftEq *eq,
                          const float *in_l, const float *in_r,
                          float *out_l, float *out_r) {
    int n   = eq->n;
    int hop = eq->hop;
    const float *window = eq->window;
    const float *gains  = eq->gains;
    float *work = eq->fft_work;

    /* 1. Shift overlap buffers left by hop */
    memmove(eq->overlap_l, eq->overlap_l + hop,
            (size_t)(n - hop) * sizeof(float));
    memmove(eq->overlap_r, eq->overlap_r + hop,
            (size_t)(n - hop) * sizeof(float));

    /* 2. Copy new samples into the right half */
    memcpy(eq->overlap_l + n - hop, in_l, (size_t)hop * sizeof(float));
    memcpy(eq->overlap_r + n - hop, in_r, (size_t)hop * sizeof(float));

    /* 3. Pack L (real) + R (imag) into one complex buffer, apply window.
     *    FFT is linear: FFT(L + jR) = FFT(L) + j·FFT(R).
     *    Applying the same real gain g[k] to both preserves this. */
    for (int i = 0; i < n; i++) {
        float w = window[i];
        work[2*i]     = eq->overlap_l[i] * w;
        work[2*i + 1] = eq->overlap_r[i] * w;
    }

    /* 4. Forward FFT (one instead of two — ~2× speedup) */
    fft_cpx(work, n, 0);

    /* 5. Multiply by per-bin gains (real gains preserve linearity) */
    work[0] *= gains[0];
    work[1] *= gains[0];
    for (int k = 1; k < n/2; k++) {
        float g = gains[k];
        work[2*k]         *= g;
        work[2*k + 1]     *= g;
        work[2*(n - k)]     *= g;
        work[2*(n - k) + 1] *= g;
    }
    work[n]     *= gains[n/2];   /* Nyquist real */
    work[n + 1] *= gains[n/2];   /* Nyquist imag */

    /* 6. Inverse FFT (one instead of two) */
    fft_cpx(work, n, 1);

    /* 7. Extract: real→L, imag→R. Overlap-add each channel independently. */
    for (int i = 0; i < hop; i++) {
        out_l[i] = work[2*i]     + eq->olap_add_l[i];
        out_r[i] = work[2*i + 1] + eq->olap_add_r[i];
    }

    /* 8. Save second half as overlap for next frame */
    for (int i = 0; i < hop; i++) {
        eq->olap_add_l[i] = work[2*(i + hop)];
        eq->olap_add_r[i] = work[2*(i + hop) + 1];
    }
}
