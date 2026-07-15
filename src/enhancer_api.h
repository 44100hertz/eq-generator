/**
 * enhancer_api.h — Clean opaque API for Python FFI (float)
 *
 * Hides DspPipe internals behind a simple create/destroy/process API.
 * Suitable for ctypes FFI.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* ── Types ──────────────────────────────────────────────────────────── */

/** Opaque enhancer handle. */
typedef struct DspPipeHandle DspPipeHandle;

/** All coefficients needed to create an enhancer. */
typedef struct {
    float   cutoff_hz;
    float   h2_amp;
    float   h3_amp;
    float   release_secs;
    float   push_gain;     /* headroom fill strength (1.0 = 0 dBFS) */
    float   pre_gain;      /* linear gain applied before EQ (typ. 1.0) */
    float   fs;
    int     eq_n_biquads;
    const float *eq_coeffs;  /* array of 5*n float coefficients */
} DspPipeParams;

/* ── API ───────────────────────────────────────────────────────────── */

/** Create an enhancer. Returns NULL on failure. */
DspPipeHandle *dsp_pipe_create(const DspPipeParams *params);

/** Destroy an enhancer. */
void dsp_pipe_destroy(DspPipeHandle *enh);

/** Reset all filter/envelope state. */
void dsp_pipe_handle_reset(DspPipeHandle *enh);

/** Process one stereo frame.
 *  Input:  float samples (range [-1.0, 1.0])
 *  Output: float samples (range [-1.0, 1.0]), written in-place.
 */
void dsp_pipe_handle_process_stereo(DspPipeHandle *enh, float *left, float *right);

#ifdef __cplusplus
}
#endif
