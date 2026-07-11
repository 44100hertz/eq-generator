/**
 * lp.h — First-order low-pass filter (float)
 *
 * y[n] = y[n-1] + α·(x[n] - y[n-1])  where α = 1 - exp(-2π·fc/fs)
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float y1;    /* previous output */
    float alpha; /* smoothing coefficient: α = 1 - exp(-2π·fc/fs) */
} LP;

/** Initialize with α coefficient. */
static inline void lp_init(LP *lp, float alpha) {
    lp->y1    = 0.0f;
    lp->alpha = alpha;
}

/** Process one sample. */
static inline float lp_tick(LP *lp, float x) {
    lp->y1 = lp->y1 + lp->alpha * (x - lp->y1);
    return lp->y1;
}

static inline void lp_reset(LP *lp) {
    lp->y1 = 0.0f;
}

#ifdef __cplusplus
}
#endif
