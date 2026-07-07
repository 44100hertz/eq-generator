/**
 * enhancer_api.h — Clean opaque API for Python FFI
 *
 * Hides BassEnhancer internals behind a simple create/destroy/process API.
 * Suitable for ctypes FFI.
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Types ──────────────────────────────────────────────────────────── */

/** Opaque enhancer handle. */
typedef struct EnhancerHandle EnhancerHandle;

/** All coefficients needed to create an enhancer. */
typedef struct {
    float   cutoff_hz;
    float   h2_amp;
    float   h3_amp;
    float   release_secs;
    float   limiter_release_secs;
    float   pre_gain;      /* linear gain applied before EQ (typ. 1.0) */
    float   fs;
    int     eq_n_biquads;
    const int32_t *eq_coeffs_q28;  /* array of 5*n Q4.28 coefficients */
} EnhancerParams;

/* ── API ───────────────────────────────────────────────────────────── */

/** Create an enhancer. Returns NULL on failure. */
EnhancerHandle *enhancer_create(const EnhancerParams *params);

/** Destroy an enhancer. */
void enhancer_destroy(EnhancerHandle *enh);

/** Reset all filter/envelope state. */
void enhancer_reset(EnhancerHandle *enh);

/** Process one stereo frame.
 *  Input:  int16 samples (range [-32768, 32767])
 *  Output: int16 samples (range [-32768, 32767]), written in-place.
 */
void enhancer_process_stereo(EnhancerHandle *enh, int16_t *left, int16_t *right);

#ifdef __cplusplus
}
#endif
