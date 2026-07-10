/**
 * biquad_q28.h — Direct Form I biquad with Q4.28 coefficients
 *
 * Input state (x1, x2) is Q16 (int32_t).
 * Output state (y1, y2) is Q32 (int64_t) — 16 extra fractional bits
 *   to preserve the slow decay of near-unity poles at low fc/fs.
 * Coefficients are Q4.28 (int32_t).
 * Accumulation: 64-bit Q(4.60) sum, single truncation to Q16.
 *
 * Usage:
 *   BiquadQ28 bq[n];
 *   BiquadQ28_init(&bq[i], &coeffs[i * 5]);
 *   y = BiquadQ28_tick(&bq[i], x);
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Biquad state ──────────────────────────────────────────────────── */
typedef struct {
    const int32_t *coeffs;  /* pointer to 5 Q4.28 values: [b0, b1, b2, a1, a2] */
    int32_t x1, x2;         /* delay-line: input history (Q16) */
    int64_t y1, y2;         /* delay-line: output history (Q32) — 16 extra bits */
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

/* ── Tick (precise, Q44 accumulation) ──────────────────────────────── */

/** Process one sample through the biquad with full Q44 accumulation.
 *
 *  All 5 products are summed in a 64-bit accumulator before a single
 *  truncation to Q16.  Preserves exact coefficient cancellation
 *  (e.g. b0+b1+b2=0 for Butterworth HP).  Used as the reference
 *  implementation on non-Xtensa platforms and for cancellation-
 *  sensitive HP filters on all platforms.
 */
/** Q16→Q32 upscale for matching precision with y1/y2 state. */
#define BQ28_Q16_TO_Q60(x)  (((int64_t)(x)) << 16)

static inline int32_t BiquadQ28_tick_q44(BiquadQ28 *bq, int32_t x) {
    const int32_t *c = bq->coeffs;
    int64_t acc;

    /* x-terms in Q60:  b(Q4.28) * x(Q16) * 2^16 = Q(4.60) */
    acc  = (int64_t)(x)       * (int64_t)(c[0]);          /* Q(4.44) → scale below */
    acc += (int64_t)(bq->x1)  * (int64_t)(c[1]);
    acc += (int64_t)(bq->x2)  * (int64_t)(c[2]);
    acc *= 65536LL;                                       /* Q(4.60) — avoid UB on signed << */

    /* y-terms already Q32:  coeff(Q4.28) * state(Q32) = Q(4.60)  */
    acc -= (int64_t)(bq->y1)  * (int64_t)(c[3]);
    acc -= (int64_t)(bq->y2)  * (int64_t)(c[4]);

    /* Store full-precision state directly from the accumulator.
     * acc is Q(4.60); >>28 gives Q32 (32 frac bits, preserving decay).
     * The Q16 output is the top 16 bits of the Q32 state.           */
    int64_t y_q32 = acc >> 28;       /* Q32: 16 frac extra over Q16 */
    int32_t y_q16 = (int32_t)(y_q32 >> 16);

    bq->x2 = bq->x1;
    bq->x1 = x;
    bq->y2 = bq->y1;
    bq->y1 = y_q32;

    return y_q16;
}

/* ── Tick (default) ────────────────────────────────────────────────── */

/** Process one sample through the biquad.
 *
 *  All 5 products are summed in a 64-bit accumulator before a single
 *  truncation to Q16.  The compiler generates carry-branch pairs on
 *  32-bit architectures — empirically ~48 cy/tick on Xtensa LX6.
 *
 *  For cancellation-sensitive filters (Butterworth HP) use
 *  BiquadQ28_tick_q44 which is this same implementation; the alias
 *  exists for clarity at call sites.
 */
static inline int32_t BiquadQ28_tick(BiquadQ28 *bq, int32_t x) {
    return BiquadQ28_tick_q44(bq, x);
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
