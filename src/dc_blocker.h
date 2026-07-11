/**
 * dc_blocker.h — First-order DC-blocking IIR filter (float)
 *
 * Transfer function: H(z) = (1 - z⁻¹) / (1 - R*z⁻¹)
 * R = exp(-2π*fc/fs)
 *
 * Blocks DC completely (numerator is exactly 0 at z=1).
 * Typical fc = 5 Hz at fs=44100 → R ≈ 0.99929
 */

#pragma once

#include "fpu.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    float R;       /* pole coefficient: exp(-2π*fc/fs) */
    float x1;      /* delayed input */
    float y1;      /* delayed output */
} DCBlocker;

/** Initialize DC blocker.
 *  R = exp(-2π*fc/fs).
 *  For fc=5 Hz, fs=44100: R ≈ 0.99929.
 */
static inline void dc_blocker_init(DCBlocker *dc, float R) {
    dc->R  = R;
    dc->x1 = 0.0f;
    dc->y1 = 0.0f;
}

/** Process one sample. y = x - x1 + R*y1 */
static inline float dc_blocker_tick(DCBlocker *dc, float x) {
    float y = ftz(x - dc->x1 + dc->R * dc->y1);
    dc->x1 = x;
    dc->y1 = y;
    return y;
}

static inline void dc_blocker_reset(DCBlocker *dc) {
    dc->x1 = 0.0f;
    dc->y1 = 0.0f;
}

#ifdef __cplusplus
}
#endif
