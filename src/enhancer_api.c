/**
 * enhancer_api.c — Opaque API for Python FFI (float)
 *
 * Wraps DspPipe behind a simple create/destroy/process API.
 * Callers are responsible for int16 ↔ float conversion.
 */

#include "enhancer_api.h"
#include "enhancer.h"
#include "biquad.h"
#include "dc_blocker.h"
#include "envelope.h"
#include <stdlib.h>
#include <string.h>

/* ── DspPipeHandle: the opaque struct ──────────────────────────────── */

struct DspPipeHandle {
    DspPipe      enh;
    Biquad           *eq_bqs_left;
    Biquad           *eq_bqs_right;
};

/* ── Create / destroy ───────────────────────────────────────────────── */

DspPipeHandle *dsp_pipe_create(const DspPipeParams *p) {
    DspPipeHandle *h = calloc(1, sizeof(DspPipeHandle));
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
    DspPipeCfg cfg;
    dsp_pipe_cfg_init(&cfg,
                         p->cutoff_hz, p->h2_amp, p->h3_amp,
                         p->release_secs, p->fs,
                         p->push_gain,
                         p->pre_gain,
                         p->eq_n_biquads, p->eq_coeffs);

    dsp_pipe_init(&h->enh, &cfg, h->eq_bqs_left, h->eq_bqs_right);

    return h;
}

void dsp_pipe_destroy(DspPipeHandle *h) {
    if (!h) return;
    free(h->eq_bqs_left);
    free(h->eq_bqs_right);
    free(h);
}

void dsp_pipe_handle_reset(DspPipeHandle *h) {
    if (h) dsp_pipe_reset(&h->enh);
}

/* ── Process ────────────────────────────────────────────────────────── */

void dsp_pipe_handle_process_stereo(DspPipeHandle *h, float *left, float *right) {
    dsp_pipe_process_stereo(&h->enh, left, right);
}
