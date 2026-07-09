/**
 * lp_q16.h — First-order LP filter in Q16 (no Q28 truncation issue)
 *
 * y[n] = y[n-1] + α·(x[n] - y[n-1])  where α = 1 - exp(-2π·fc/fs)
 *
 * For fc=60 Hz, fs=44100: α ≈ 0.00848 → Q16 ≈ 556.
 * This avoids the Q28 biquad problem where b0 is so small (≈4875)
 * that normal signal levels truncate to zero.
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int32_t y1;    /* previous output (Q16) */
    int32_t alpha; /* smoothing coefficient (Q16): α = 1 - exp(-2π·fc/fs) */
} LPQ16;

/** Initialize with α coefficient in Q16. */
static inline void LPQ16_init(LPQ16 *lp, int32_t alpha_q16) {
    lp->y1    = 0;
    lp->alpha = alpha_q16;
}

/** Process one sample. x and return are Q16. */
static inline int32_t LPQ16_tick(LPQ16 *lp, int32_t x) {
    /* y = y1 + α·(x - y1) / 2^16 */
    int32_t diff = x - lp->y1;
    lp->y1 = lp->y1 + (int32_t)(((int64_t)lp->alpha * (int64_t)diff) >> 16);
    return lp->y1;
}

static inline void LPQ16_reset(LPQ16 *lp) {
    lp->y1 = 0;
}

#ifdef __cplusplus
}
#endif
