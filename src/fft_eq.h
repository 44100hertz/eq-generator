/**
 * fft_eq.h — Overlap-add FFT EQ (float, 256-point, 50% overlap)
 *
 * Applies per-bin scalar gains via frequency-domain multiplication.
 * Input: 128 new float samples per tick (after IIR biquad pass).
 * Hann window + 50% overlap-add for perfect reconstruction when
 * gains are unity.
 *
 * State (per instance, stereo):
 *   overlap_l[256], overlap_r[256]  — input history + window target
 *   olap_add_l[128], olap_add_r[128] — overlap-add accumulator
 *   window[256]                      — pre-computed Hann window
 *   gains[129]                       — per-bin gains (DC..Nyquist)
 *   fft_work[512]                    — FFT workspace (256 complex, reusable)
 *
 * RAM: ~4.5 KB (2×1024 + 2×512 + 1024 + 516 + 2048 bytes)
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#define FFT_EQ_N    256
#define FFT_EQ_HOP  128
#define FFT_EQ_BINS 129   /* gains table: DC through Nyquist (n/2 + 1) */

typedef struct {
    int   n;          /* window size (256) */
    int   hop;        /* hop size (128) */

    float *window;    /* pre-computed Hann window [n] */
    float *gains;     /* per-bin scalar gains [n/2 + 1] */

    /* Per-channel overlap-add state */
    float *overlap_l;     /* [n] input history + window target */
    float *overlap_r;     /* [n] */
    float *olap_add_l;    /* [hop] accumulated overlap from prev frame */
    float *olap_add_r;    /* [hop] */

    /* Workspace (reused between channels) */
    float *fft_work;  /* [2*n] = [512] complex interleaved */
} FftEq;

/** Allocate and initialize an FftEq. gains must point to FFT_EQ_BINS+1 floats
 *  (indices 0..128, DC through Nyquist). Window and overlap buffers are
 *  allocated internally. Returns NULL on allocation failure. */
FftEq *fft_eq_create(const float *gains);

/** Free all memory associated with the FftEq. */
void fft_eq_destroy(FftEq *eq);

/** Reset overlap-add state to zero (call after parameter changes). */
void fft_eq_reset(FftEq *eq);

/** Process one frame (128 stereo samples) through the FFT EQ.
 *  in_l, in_r: 128 new float samples (from IIR output)
 *  out_l, out_r: 128 output samples.
 *  All buffers must be non-overlapping and 16-byte aligned for safety. */
void fft_eq_process_frame(FftEq *eq,
                          const float *in_l, const float *in_r,
                          float *out_l, float *out_r);

#ifdef __cplusplus
}
#endif
