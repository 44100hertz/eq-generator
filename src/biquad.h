/**
 * biquad.h — Direct Form I biquad with float coefficients
 *
 * Coefficients are 5 floats: [b0, b1, b2, a1, a2].
 * State is 4 floats: x1, x2, y1, y2.
 *
 * Usage:
 *   Biquad bq[n];
 *   biquad_init(&bq[i], &coeffs[i * 5]);
 *   y = biquad_tick(&bq[i], x);
 */

#pragma once

#include "fpu.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Biquad state ──────────────────────────────────────────────────── */
typedef struct {
    const float *coeffs;  /* pointer to 5 floats: [b0, b1, b2, a1, a2] */
    float x1, x2;         /* input history */
    float y1, y2;         /* output history */
} Biquad;

/* ── Initialization ────────────────────────────────────────────────── */

/** Initialize a biquad from a pointer to 5 float coefficients. */
static inline void biquad_init(Biquad *bq, const float *co) {
    bq->coeffs = co;
    bq->x1 = 0.0f;
    bq->x2 = 0.0f;
    bq->y1 = 0.0f;
    bq->y2 = 0.0f;
}

/** Reset biquad state to zero without changing coefficients. */
static inline void biquad_reset(Biquad *bq) {
    bq->x1 = 0.0f;
    bq->x2 = 0.0f;
    bq->y1 = 0.0f;
    bq->y2 = 0.0f;
}

/* ── Tick ──────────────────────────────────────────────────────────── */

/** Process one sample through the biquad.
 *
 *  y = b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2
 */
__attribute__((always_inline))
static inline float biquad_tick(Biquad *bq, float x) {
    const float *c = bq->coeffs;
    float y = ftz(c[0] * x + c[1] * bq->x1 + c[2] * bq->x2
                - c[3] * bq->y1 - c[4] * bq->y2);

    bq->x2 = bq->x1;
    bq->x1 = x;
    bq->y2 = bq->y1;
    bq->y1 = y;

    return y;
}

/* ── Cascade convenience ───────────────────────────────────────────── */

/** Process one sample through a cascade of n biquads. */
static inline float biquad_cascade(Biquad *bqs, int n, float x) {
    float y = x;
    for (int i = 0; i < n; i++) {
        y = biquad_tick(&bqs[i], y);
    }
    return y;
}

#ifdef __cplusplus
}
#endif
