/**
 * fft_eq.h — Overlap-add FFT EQ (configurable N, 50% overlap)
 *
 * Applies per-bin scalar gains via frequency-domain multiplication.
 * Input: N/2 new float samples per tick.
 * Hann window + 50% overlap-add for perfect reconstruction when
 * gains are unity.
 *
 * FFT_EQ_N can be overridden at compile time (-DFFT_EQ_N=512).
 * Default: 256-point FFT, 128-sample hop, 129 bins.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** EQGEN_FFT_N is the single source of truth (defined in eq_coeffs.h).
 *  Fall back to 256 if compiling without the generated header (e.g. bench). */
#if defined(EQGEN_FFT_N)
  #define FFT_EQ_N    EQGEN_FFT_N
#elif !defined(FFT_EQ_N)
  #define FFT_EQ_N    256
#endif
#define FFT_EQ_HOP  (FFT_EQ_N / 2)
#define FFT_EQ_BINS (FFT_EQ_N / 2 + 1)   /* gains table: DC through Nyquist (n/2 + 1) */

typedef struct {
    int   n;          /* window size (FFT_EQ_N) */
    int   hop;        /* hop size (FFT_EQ_HOP) */

    float *window;    /* pre-computed Hann window [n] */
    float *gains;     /* per-bin scalar gains [n/2 + 1] */

    /* Per-channel overlap-add state */
    float *overlap_l;     /* [n] input history + window target */
    float *overlap_r;     /* [n] */
    float *olap_add_l;    /* [hop] accumulated overlap from prev frame */
    float *olap_add_r;    /* [hop] */

    /* Workspace (reused between channels) */
    float *fft_work;  /* [2*n] complex interleaved */
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
