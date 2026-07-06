/**
 * biquad_q28.h — Direct Form I biquad with Q4.28 coefficients
 *
 * State variables (x1, x2, y1, y2) are Q16 (int32_t).
 * Coefficients are Q4.28 (int32_t, stored in eq_coeffs_q28[]).
 * Multiply: FR_MULK28(q16_sample, q28_coeff) → Q16 result with rounding.
 * Accumulator: int64_t, truncated to int32_t on output.
 *
 * Usage:
 *   BiquadQ28 bq[n];
 *   BiquadQ28_init(&bq[i], &eq_coeffs_q28[i * 5]);
 *   y = BiquadQ28_tick(&bq[i], x);
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Q4.28 multiply: Q16 * Q4.28 → Q16 ────────────────────────────── */
/* Uses 64-bit intermediate with round-to-nearest.                      */
/* If fr_math is available, this maps to FR_MULK28.                     */

#ifndef FR_MULK28
#include <stdint.h>
static inline int32_t FR_MULK28(int32_t x, int32_t k) {
    return (int32_t)((((int64_t)x * (int64_t)k) + ((int64_t)1 << 27)) >> 28);
}
#endif

/* ── Biquad state ──────────────────────────────────────────────────── */
typedef struct {
    const int32_t *coeffs;  /* pointer to 5 Q4.28 values: [b0, b1, b2, a1, a2] */
    int32_t x1, x2;         /* delay-line: input history (Q16) */
    int32_t y1, y2;         /* delay-line: output history (Q16) */
} BiquadQ28;

/* ── Initialization ────────────────────────────────────────────────── */

/** Initialize a biquad from a pointer to 5 Q4.28 coefficients.
 *  co must point to a static array of 5 int32_t: [b0, b1, b2, a1, a2].
 */
static inline void BiquadQ28_init(BiquadQ28 *bq, const int32_t *co) {
    bq->coeffs = co;
    bq->x1 = 0;
    bq->x2 = 0;
    bq->y1 = 0;
    bq->y2 = 0;
}

/** Reset biquad state to zero without changing coefficients. */
static inline void BiquadQ28_reset(BiquadQ28 *bq) {
    bq->x1 = 0;
    bq->x2 = 0;
    bq->y1 = 0;
    bq->y2 = 0;
}

/* ── Tick ──────────────────────────────────────────────────────────── */

/** Process one sample through the biquad.
 *
 *  Direct Form I, 64-bit accumulator in Q44, single round at end:
 *    y = round((b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2) / 2^28)
 *
 *  Accumulating in Q44 preserves exact coefficient cancellation
 *  (e.g. b0+b1+b2=0 for Butterworth HP), avoiding DC-bias artifacts.
 *
 *  Returns int32_t Q16 output (range approximately [-32767, 32767]).
 */
static inline int32_t BiquadQ28_tick(BiquadQ28 *bq, int32_t x) {
    const int32_t *c = bq->coeffs;
    int64_t acc;

    /* All 5 multiply-accumulates in Q44 (32*32=64 bit product) */
    acc  = (int64_t)(x)       * (int64_t)(c[0]);  /* b0 * x  */
    acc += (int64_t)(bq->x1)  * (int64_t)(c[1]);  /* b1 * x1 */
    acc += (int64_t)(bq->x2)  * (int64_t)(c[2]);  /* b2 * x2 */
    acc -= (int64_t)(bq->y1)  * (int64_t)(c[3]);  /* a1 * y1 */
    acc -= (int64_t)(bq->y2)  * (int64_t)(c[4]);  /* a2 * y2 */

    /* Single truncation from Q44 to Q16 (no rounding bias — avoids DC offset
     * amplification when 1+a1+a2 is tiny). Truncation error is mean-zero. */
    int32_t y = (int32_t)(acc >> 28);

    /* Shift delay line */
    bq->x2 = bq->x1;
    bq->x1 = x;
    bq->y2 = bq->y1;
    bq->y1 = y;

    return y;
}

/* ── Cascade convenience ───────────────────────────────────────────── */

/** Process one sample through a cascade of n biquads.
 *
 *  The output of each stage feeds the input of the next.
 *  Returns the final output.
 */
static inline int32_t BiquadQ28_cascade(BiquadQ28 *bqs, int n, int32_t x) {
    int32_t y = x;
    for (int i = 0; i < n; i++) {
        y = BiquadQ28_tick(&bqs[i], y);
    }
    return y;
}

#ifdef __cplusplus
}
#endif
