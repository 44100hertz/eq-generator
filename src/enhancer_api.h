/**
 * enhancer_api.h — Clean opaque API for Python FFI (float)
 *
 * Hides BassEnhancer internals behind a simple create/destroy/process API.
 * Suitable for ctypes FFI.
 */

#pragma once

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
    const float *eq_coeffs;  /* array of 5*n float coefficients */
} EnhancerParams;

/* ── API ───────────────────────────────────────────────────────────── */

/** Create an enhancer. Returns NULL on failure. */
EnhancerHandle *enhancer_create(const EnhancerParams *params);

/** Destroy an enhancer. */
void enhancer_destroy(EnhancerHandle *enh);

/** Reset all filter/envelope state. */
void enhancer_reset(EnhancerHandle *enh);

/** Process one stereo frame.
 *  Input:  float samples (range [-1.0, 1.0])
 *  Output: float samples (range [-1.0, 1.0]), written in-place.
 */
void enhancer_process_stereo(EnhancerHandle *enh, float *left, float *right);

#ifdef __cplusplus
}
#endif
