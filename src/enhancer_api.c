/**
 * enhancer_api.c — Opaque API for Python FFI (float)
 *
 * Wraps BassEnhancer behind a simple create/destroy/process API.
 * Callers are responsible for int16 ↔ float conversion.
 */

#include "enhancer_api.h"
#include "enhancer.h"
#include "biquad.h"
#include "dc_blocker.h"
#include "envelope.h"
#include <stdlib.h>
#include <string.h>

/* ── EnhancerHandle: the opaque struct ──────────────────────────────── */

struct EnhancerHandle {
    BassEnhancer      enh;
    Biquad           *eq_bqs_left;
    Biquad           *eq_bqs_right;
};

/* ── Create / destroy ───────────────────────────────────────────────── */

EnhancerHandle *enhancer_create(const EnhancerParams *p) {
    EnhancerHandle *h = calloc(1, sizeof(EnhancerHandle));
    if (!h) return NULL;

    /* Allocate EQ biquad state arrays */
    int n = p->eq_n_biquads;
    h->eq_bqs_left  = calloc((size_t)n, sizeof(Biquad));
    h->eq_bqs_right = calloc((size_t)n, sizeof(Biquad));
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
                         p->push_gain,
                         p->pre_gain,
                         p->eq_n_biquads, p->eq_coeffs);

    BassEnhancer_init(&h->enh, &cfg, h->eq_bqs_left, h->eq_bqs_right);

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

void enhancer_process_stereo(EnhancerHandle *h, float *left, float *right) {
    BassEnhancer_process_stereo(&h->enh, left, right);
}
