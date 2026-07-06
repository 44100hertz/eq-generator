/**
 * dc_blocker.h — First-order DC-blocking IIR filter (Q16)
 *
 * Transfer function: H(z) = (1 - z⁻¹) / (1 - R*z⁻¹)
 * R = exp(-2π*fc/fs) ≈ 1 - 2π*fc/fs for fc << fs
 *
 * No near-cancellation: coefficients are well-behaved in Q16.
 * Blocks DC completely (numerator is exactly 0 at z=1 in Q16 too
 * since 1 + (-1) = 0).
 *
 * Typical fc = 5 Hz at fs=44100 → R ≈ 0.99929 → Q16(R) ≈ 65489
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int32_t R;       /* pole coefficient in Q16: exp(-2π*fc/fs) * 65536 */
    int32_t x1;      /* delayed input */
    int32_t y1;      /* delayed output */
} DCBlocker;

/** Initialize DC blocker with R coefficient in Q16.
 *  R ≈ 65536 * exp(-2π*fc/fs).
 *  For fc=5 Hz, fs=44100: R ≈ 65489.
 */
static inline void DCBlocker_init(DCBlocker *dc, int32_t R_q16) {
    dc->R  = R_q16;
    dc->x1 = 0;
    dc->y1 = 0;
}

/** Process one sample.
 *  y = x - x1 + R*y1 (all in Q16, rounded).
 */
static inline int32_t DCBlocker_tick(DCBlocker *dc, int32_t x) {
    /* y = x - x1 + (R * y1) / 65536 */
    int32_t y = x - dc->x1 + (int32_t)(((int64_t)dc->R * (int64_t)dc->y1) >> 16);
    dc->x1 = x;
    dc->y1 = y;
    return y;
}

static inline void DCBlocker_reset(DCBlocker *dc) {
    dc->x1 = 0;
    dc->y1 = 0;
}

#ifdef __cplusplus
}
#endif
