/**
 * filter.c — stdin→stdout EQGen DSP filter (float)
 *
 * Reads raw float32le stereo interleaved samples from stdin,
 * processes through BassEnhancer_process_stereo, writes to stdout.
 * Uses baked-in eq_coeffs.h coefficients and supports sample-rate
 * independent processing.
 *
 * Usage:
 *   parec -d eqgen_sink.monitor --format=float32le --rate=48000 \
 *     | ./filter [sample_rate] \
 *     | pacat -d <sink> --format=float32le --rate=48000
 *
 * Build:
 *   make filter  (from src/)
 */

#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "eq_coeffs.h"
#include "enhancer.h"
#include "fft_eq.h"
#include "smart_volume.h"

static volatile int g_running = 1;

static void sig_handler(int sig) {
    (void)sig;
    g_running = 0;
}

int main(int argc, char *argv[]) {
    float fs = 48000.0f;
    if (argc > 1) {
        fs = (float)atof(argv[1]);
        if (fs < 1000.0f) fs = 48000.0f;
    }

    signal(SIGINT,  sig_handler);
    signal(SIGTERM, sig_handler);
    signal(SIGPIPE, sig_handler);

    /* Allocate enhancer with baked-in EQ coefficients. */
    int n_bq             = EQGEN_N_BIQUADS;
    Biquad *eq_bqs_l     = calloc((size_t)n_bq, sizeof(Biquad));
    Biquad *eq_bqs_r     = calloc((size_t)n_bq, sizeof(Biquad));
    if (!eq_bqs_l || !eq_bqs_r) {
        fprintf(stderr, "filter: out of memory\n");
        return 1;
    }

    const float *coeffs = eqgen_get_coeffs((int)fs);

    BassEnhancer enh;
    BassEnhancerCfg cfg;

    float pre_gain_f = 1.0f;
#ifdef EQGEN_PRE_GAIN
    pre_gain_f = EQGEN_PRE_GAIN;
#endif
    const float pg_ref = pre_gain_f;  /* never mutate — used as loud reference */

    /* Volume LUT — mirrors ESP32 rebuild_vol_lut() exactly. */
    float vol_lut[128];
    smart_volume_rebuild_lut(vol_lut, 0.0f);

    /* ── FFT hybrid EQ (if enabled) ─────────────────────────────── */
    FftEq *fft_eq = NULL;
#if EQGEN_FFT_BINS > 0
    {
        const float *fft_gains = eqgen_get_fft_gains((int)fs);
        fft_eq = fft_eq_create(fft_gains);
        if (!fft_eq) {
            fprintf(stderr, "filter: FFT EQ allocation failed\n");
        } else {
            fprintf(stderr, "filter: FFT EQ enabled (%d bins)\n", EQGEN_FFT_BINS);
        }
    }
#endif

    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_LIMITER_RELEASE_SECS,
                         fft_eq ? 1.0f : pre_gain_f,
                         n_bq,
                         coeffs);

    /* Initialise loudness shelf alpha (boost stays 0 until first control msg) */
    BassEnhancerCfg_set_loudness(&cfg, EQGEN_LOUDNESS_FC_HZ, fs, 0.0f);

    BassEnhancer_init(&enh, &cfg, eq_bqs_l, eq_bqs_r);

    /* Open control FIFO for live smart-volume adjustments (optional). */
    const char *ctrl_path = "/tmp/eqgen_sv_fifo";
    int ctrl_fd = open(ctrl_path, O_RDWR | O_NONBLOCK);
    if (ctrl_fd >= 0) {
        fprintf(stderr, "filter: smart-volume control on %s\n", ctrl_path);
    }

    fprintf(stderr, "filter: running at %.0f Hz\n", fs);

    float buf_interleaved[256]; /* 128 stereo frames: L,R,L,R,... */
    float buf_l[128];           /* left channel for FFT pass */
    float buf_r[128];           /* right channel for FFT pass */
    size_t frame_sz = 2 * sizeof(float);
    unsigned long frame_count = 0;
    unsigned long partial_count = 0;

    while (g_running) {
        /* ── Check smart-volume control FIFO (non-blocking) ────── */
        if (ctrl_fd >= 0) {
            unsigned char vol_byte;
            ssize_t nr;
            while ((nr = read(ctrl_fd, &vol_byte, 1)) > 0) {
                uint8_t vol = vol_byte;
                if (vol <= 127) {
                    SmartVolumeParams svp;
                    smart_volume_compute(vol, pg_ref, &svp);

                    smart_volume_rebuild_lut(vol_lut, svp.shelf_db);

                    if (fft_eq) {
                        pre_gain_f = svp.fft_pre_gain;
                        BassEnhancer_update_params(&enh, svp.h2_amp, svp.h3_amp,
                                                   NAN, svp.boost, svp.bleed);
                    } else {
                        BassEnhancer_update_params(&enh, svp.h2_amp, svp.h3_amp,
                                                   svp.pre_gain, svp.boost, svp.bleed);
                    }

                    fprintf(stderr, "filter: vol=%u h2=%.3f h3=%.3f pg=%.3f shelf=%.1f dB bleed=%.3f\n",
                            (unsigned)vol,
                            (double)svp.h2_amp,
                            (double)svp.h3_amp,
                            (double)svp.pre_gain,
                            (double)svp.shelf_db,
                            (double)svp.bleed);
                }
            }
            if (nr == 0 || (nr < 0 && errno != EAGAIN && errno != EWOULDBLOCK && errno != EINTR)) {
                close(ctrl_fd);
                ctrl_fd = open(ctrl_path, O_RDWR | O_NONBLOCK);
            }
        }

        size_t nread = fread(buf_interleaved, frame_sz, 128, stdin);
        if (nread == 0) {
            if (feof(stdin)) break;
            if (ferror(stdin)) { perror("filter: read error"); break; }
            continue;
        }
        frame_count++;
        if (nread != 128) {
            partial_count++;
            for (size_t i = nread; i < 128; i++) {
                buf_interleaved[i * 2]     = 0.0f;
                buf_interleaved[i * 2 + 1] = 0.0f;
            }
            if (partial_count <= 5 || partial_count % 1000 == 0) {
                fprintf(stderr, "filter: partial read %zu/128 (frame %lu, count %lu)\n",
                        nread, frame_count, partial_count);
            }
        }

        /* Deinterleave input */
        for (size_t i = 0; i < nread; i++) {
            buf_l[i] = buf_interleaved[i * 2];
            buf_r[i] = buf_interleaved[i * 2 + 1];
        }

        /* ── Pass 1: FFT EQ (broad correction) ────────────────── */
        if (fft_eq) {
            fft_eq_process_frame(fft_eq, buf_l, buf_r, buf_l, buf_r);
            for (size_t i = 0; i < nread; i++) {
                buf_l[i] *= pre_gain_f;
                buf_r[i] *= pre_gain_f;
            }
        }

        /* ── Pass 2: IIR biquad cascade + enhancer (surgical) ──── */
        for (size_t i = 0; i < nread; i++) {
            float l = buf_l[i];
            float r = buf_r[i];

            BassEnhancer_process_stereo(&enh, &l, &r);

            buf_l[i] = l;
            buf_r[i] = r;
        }

        /* ── Interleave output ─────────────────────────────────── */
        float vol_gain = vol_lut[127];
        for (size_t i = 0; i < nread; i++) {
            buf_interleaved[i * 2]     = buf_l[i] * vol_gain;
            buf_interleaved[i * 2 + 1] = buf_r[i] * vol_gain;
        }

        size_t nwritten = fwrite(buf_interleaved, frame_sz, nread, stdout);
        if (nwritten < nread) {
            perror("filter: write error");
            break;
        }
        fflush(stdout);
    }

    if (fft_eq) fft_eq_destroy(fft_eq);
    if (ctrl_fd >= 0) close(ctrl_fd);
    free(eq_bqs_l);
    free(eq_bqs_r);
    fprintf(stderr, "filter: done (%lu frames, %lu partial)\n",
            frame_count, partial_count);
    return 0;
}
