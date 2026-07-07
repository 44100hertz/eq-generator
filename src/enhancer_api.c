/**
 * enhancer_api.c — Opaque API for Python FFI
 *
 * Wraps BassEnhancer behind a simple create/destroy/process API.
 * Handles int16 ↔ Q16 conversion internally.
 */

#include "enhancer_api.h"
#include "enhancer.h"
#include "biquad_q28.h"
#include "dc_blocker.h"
#include "envelope.h"
#include <stdlib.h>
#include <string.h>

/* ── EnhancerHandle: the opaque struct ──────────────────────────────── */

struct EnhancerHandle {
    BassEnhancer      enh;
    BiquadQ28        *eq_bqs_left;
    BiquadQ28        *eq_bqs_right;
    ReciprocalLUT     lut;  /* stored in handle so env pointers stay valid */
};

/* ── Create / destroy ───────────────────────────────────────────────── */

EnhancerHandle *enhancer_create(const EnhancerParams *p) {
    EnhancerHandle *h = calloc(1, sizeof(EnhancerHandle));
    if (!h) return NULL;

    /* Allocate EQ biquad state arrays */
    int n = p->eq_n_biquads;
    h->eq_bqs_left  = calloc((size_t)n, sizeof(BiquadQ28));
    h->eq_bqs_right = calloc((size_t)n, sizeof(BiquadQ28));
    if ((n > 0) && (!h->eq_bqs_left || !h->eq_bqs_right)) {
        free(h->eq_bqs_left);
        free(h->eq_bqs_right);
        free(h);
        return NULL;
    }

    /* Configure and initialize */
    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
                         p->cutoff_hz, p->h2_amp, p->h3_amp,
                         p->release_secs, p->fs,
                         p->limiter_release_secs,
                         p->pre_gain,
                         p->eq_n_biquads, p->eq_coeffs_q28);

    ReciprocalLUT_init(&h->lut);

    BassEnhancer_init(&h->enh, &cfg, &h->lut, h->eq_bqs_left, h->eq_bqs_right);

    return h;
}

void enhancer_destroy(EnhancerHandle *h) {
    if (!h) return;
    free(h->eq_bqs_left);
    free(h->eq_bqs_right);
    free(h);
}

void enhancer_reset(EnhancerHandle *h) {
    if (h) BassEnhancer_reset(&h->enh);
}

/* ── Process ────────────────────────────────────────────────────────── */

void enhancer_process_stereo(EnhancerHandle *h, int16_t *left, int16_t *right) {
    /* int16 → Q16: shift left by 1 (multiply by 2).
     * int16 max = 32767. Q16 where 0.5 = 32768. So int16 * 2 = Q16. */
    int32_t l_q16 = ((int32_t)*left) << 1;
    int32_t r_q16 = ((int32_t)*right) << 1;

    BassEnhancer_process_stereo(&h->enh, &l_q16, &r_q16);

    /* Q16 → int16: shift right by 1, clamp */
    l_q16 >>= 1;
    r_q16 >>= 1;
    if (l_q16 >  32767) l_q16 =  32767;
    if (l_q16 < -32768) l_q16 = -32768;
    if (r_q16 >  32767) r_q16 =  32767;
    if (r_q16 < -32768) r_q16 = -32768;

    *left  = (int16_t)l_q16;
    *right = (int16_t)r_q16;
}
