/**
 * eqgen_filter.c — PipeWire DSP filter wrapping enhancer.c.
 *
 * Follows the official pipewire tutorial (page_tutorial7.html).
 * Uses pw_filter_new_simple with PW_KEY_FORMAT_DSP port properties.
 */
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <spa/pod/builder.h>
#include <spa/param/latency-utils.h>

#include <pipewire/pipewire.h>
#include <pipewire/filter.h>

#include "enhancer.h"
#include "eq_coeffs.h"

struct eqgen_data {
    struct pw_main_loop *loop;
    struct pw_filter    *filter;
    void                *in_l_port;
    void                *in_r_port;
    void                *out_l_port;
    void                *out_r_port;

    BassEnhancer         enh;
    BiquadQ28            eq_l[EQGEN_N_BIQUADS];
    BiquadQ28            eq_r[EQGEN_N_BIQUADS];
    ReciprocalLUT        lut;
};

/* ── process (RT thread) ──────────────────────────────────────────── */

static void on_process(void *userdata, struct spa_io_position *position)
{
    struct eqgen_data *d = userdata;
    uint32_t n = (uint32_t)position->clock.duration;
    if (n == 0) return;

    float *il = pw_filter_get_dsp_buffer(d->in_l_port, n);
    float *ir = pw_filter_get_dsp_buffer(d->in_r_port, n);
    float *ol = pw_filter_get_dsp_buffer(d->out_l_port, n);
    float *or_ = pw_filter_get_dsp_buffer(d->out_r_port, n);
    if (!il || !ir || !ol || !or_) return;

    for (uint32_t i = 0; i < n; i++) {
        int32_t l = (int32_t)(il[i] * 65536.0f);
        int32_t r = (int32_t)(ir[i] * 65536.0f);
        BassEnhancer_process_stereo(&d->enh, &l, &r);
        ol[i] = (float)l / 65536.0f;
        or_[i] = (float)r / 65536.0f;
    }
}

static const struct pw_filter_events filter_events = {
    PW_VERSION_FILTER_EVENTS,
    .process = on_process,
};

/* ── signal ───────────────────────────────────────────────────────── */

static void do_quit(void *userdata, int signal_number)
{
    (void)signal_number;
    struct eqgen_data *d = userdata;
    pw_main_loop_quit(d->loop);
}

/* ── main ─────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
    struct eqgen_data data = { 0 };
    const struct spa_pod *params[1];
    uint32_t n_params = 0;
    uint8_t buffer[1024];
    struct spa_pod_builder b = SPA_POD_BUILDER_INIT(buffer, sizeof(buffer));

    pw_init(&argc, &argv);

    data.loop = pw_main_loop_new(NULL);
    pw_loop_add_signal(pw_main_loop_get_loop(data.loop),
                       SIGINT, do_quit, &data);
    pw_loop_add_signal(pw_main_loop_get_loop(data.loop),
                       SIGTERM, do_quit, &data);

    /* detect actual sample rate (CLI arg or default to EQGEN_FS) */
    float actual_fs = (float)EQGEN_FS;
    if (argc > 1) {
        int arg_fs = atoi(argv[1]);
        if (arg_fs > 0) actual_fs = (float)arg_fs;
    }

    /* init DSP at the actual PipeWire sample rate */
    BassEnhancerCfg cfg;
    BassEnhancerCfg_init(&cfg,
        EQGEN_CUTOFF_HZ, EQGEN_H2_AMP, EQGEN_H3_AMP,
        EQGEN_RELEASE_SECS, actual_fs,
        EQGEN_LIMITER_RELEASE_SECS,
        EQGEN_N_BIQUADS, eqgen_coeffs_q28);
    ReciprocalLUT_init(&data.lut);
    BassEnhancer_init(&data.enh, &cfg, &data.lut, data.eq_l, data.eq_r);

    /* create filter (pw_filter_new_simple — manages core/context) */
    data.filter = pw_filter_new_simple(
        pw_main_loop_get_loop(data.loop),
        "eqgen_filter",
        pw_properties_new(
            PW_KEY_MEDIA_TYPE,     "Audio",
            PW_KEY_MEDIA_CATEGORY, "Filter",
            PW_KEY_MEDIA_ROLE,     "DSP",
            PW_KEY_NODE_PASSIVE,   "follow",
            NULL),
        &filter_events,
        &data);

    /* ports with PW_KEY_FORMAT_DSP (the canonical way per tutorial7) */
    data.in_l_port = pw_filter_add_port(data.filter,
        PW_DIRECTION_INPUT,
        PW_FILTER_PORT_FLAG_MAP_BUFFERS,
        0,
        pw_properties_new(
            PW_KEY_FORMAT_DSP,  "32 bit float mono audio",
            PW_KEY_PORT_NAME,   "input_L",
            NULL),
        NULL, 0);

    data.in_r_port = pw_filter_add_port(data.filter,
        PW_DIRECTION_INPUT,
        PW_FILTER_PORT_FLAG_MAP_BUFFERS,
        0,
        pw_properties_new(
            PW_KEY_FORMAT_DSP,  "32 bit float mono audio",
            PW_KEY_PORT_NAME,   "input_R",
            NULL),
        NULL, 0);

    data.out_l_port = pw_filter_add_port(data.filter,
        PW_DIRECTION_OUTPUT,
        PW_FILTER_PORT_FLAG_MAP_BUFFERS,
        0,
        pw_properties_new(
            PW_KEY_FORMAT_DSP,  "32 bit float mono audio",
            PW_KEY_PORT_NAME,   "output_L",
            NULL),
        NULL, 0);

    data.out_r_port = pw_filter_add_port(data.filter,
        PW_DIRECTION_OUTPUT,
        PW_FILTER_PORT_FLAG_MAP_BUFFERS,
        0,
        pw_properties_new(
            PW_KEY_FORMAT_DSP,  "32 bit float mono audio",
            PW_KEY_PORT_NAME,   "output_R",
            NULL),
        NULL, 0);

    /* latency: 10ms */
    params[n_params++] = spa_process_latency_build(&b,
        SPA_PARAM_ProcessLatency,
        &SPA_PROCESS_LATENCY_INFO_INIT(.ns = 10 * SPA_NSEC_PER_MSEC));

    if (pw_filter_connect(data.filter,
            PW_FILTER_FLAG_RT_PROCESS,
            params, n_params) < 0) {
        fprintf(stderr, "eqgen_filter: can't connect\n");
        return -1;
    }

    fprintf(stderr, "eqgen_filter: %d biquads, fc=%.0f fs=%.0f pid=%d\n",
            EQGEN_N_BIQUADS, EQGEN_CUTOFF_HZ, actual_fs, getpid());

    pw_main_loop_run(data.loop);

    fprintf(stderr, "eqgen_filter: done\n");
    pw_filter_destroy(data.filter);
    pw_main_loop_destroy(data.loop);
    pw_deinit();
    return 0;
}
