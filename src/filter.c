/**
 * filter.c — stdin→stdout EQGen DSP filter
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
 *   cc -O2 -Wall -o filter filter.c enhancer.c \
 *      -I. -lm -DBRIDGE_STANDALONE
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

    /* Allocate enhancer with baked-in EQ coefficients.
     * Select the coefficient array matching the runtime sample rate. */
    int n_bq             = EQGEN_N_BIQUADS;
    BiquadQ28 *eq_bqs_l  = calloc((size_t)n_bq, sizeof(BiquadQ28));
    BiquadQ28 *eq_bqs_r  = calloc((size_t)n_bq, sizeof(BiquadQ28));
    if (!eq_bqs_l || !eq_bqs_r) {
        fprintf(stderr, "filter: out of memory\n");
        return 1;
    }

    const int32_t *coeffs = eqgen_get_coeffs((int)fs);

    BassEnhancer enh;
    BassEnhancerCfg cfg;
    ReciprocalLUT lut;

    float pre_gain_f = 1.0f;
#ifdef EQGEN_PRE_GAIN_Q16
    pre_gain_f = (float)EQGEN_PRE_GAIN_Q16 / 65536.0f;
#endif
    const float pg_ref = pre_gain_f;  /* never mutate — used as loud reference */

    /* Volume LUT — mirrors ESP32 rebuild_vol_lut() exactly.
     * Desktop always simulates vol=127; LUT entry carries shelf compensation. */
    int32_t vol_lut[128];
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
        /* When FFT is active: apply pre_gain in float BEFORE Q16 conversion
         * to prevent clipping from FFT boost.  Bypass it in the enhancer. */
    }
#endif

    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_LIMITER_RELEASE_SECS,
                         fft_eq ? 1.0f : pre_gain_f,  /* if FFT: float path handles it */
                         n_bq,
                         coeffs);

    /* Initialise loudness shelf alpha (boost stays 0 until first control msg) */
    BassEnhancerCfg_set_loudness(&cfg, EQGEN_LOUDNESS_FC_HZ, fs, 0.0f);

    ReciprocalLUT_init(&lut);
    BassEnhancer_init(&enh, &cfg, &lut, eq_bqs_l, eq_bqs_r);

    /* Open control FIFO for live smart-volume adjustments (optional).
     * O_RDWR always succeeds — the filter only reads, but being a "writer"
     * means any O_WRONLY client (web UI) can connect without blocking. */
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
                    /* ── Shared interpolation (identical to ESP32) ── */
                    SmartVolumeParams svp;
                    smart_volume_compute(vol, pg_ref, &svp);

                    /* ── Rebuild LUT so vol_lut[127] carries shelf boost ── */
                    smart_volume_rebuild_lut(vol_lut, svp.shelf_db);

                    /* ── Apply to enhancer ── */
                    if (fft_eq) {
                        pre_gain_f = svp.fft_pre_gain;
                        BassEnhancer_update_params(&enh, svp.h2_q16, svp.h3_q16,
                                                   -1, svp.boost_q16, svp.bleed_q16);
                    } else {
                        BassEnhancer_update_params(&enh, svp.h2_q16, svp.h3_q16,
                                                   svp.pg_q16, svp.boost_q16, svp.bleed_q16);
                    }

                    fprintf(stderr, "filter: vol=%u h2=%.3f h3=%.3f pg=%.3f shelf=%.1f dB bleed=%.3f\n",
                            (unsigned)vol,
                            (double)(float)svp.h2_q16 / 65536.0,
                            (double)(float)svp.h3_q16 / 65536.0,
                            (double)(float)svp.pg_q16 / 65536.0,
                            (double)svp.shelf_db,
                            (double)(float)svp.bleed_q16 / 65536.0);
                }
            }
            /* Reopen on EOF (writer closed) or unexpected error */
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
            /* Zero-fill the rest so FFT EQ doesn't see stale data */
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
            /* Apply pre_gain before Q16 to prevent clipping from FFT boost.
             * When FFT EQ is active, enhancer pre_gain is 1.0 (no-op). */
            for (size_t i = 0; i < nread; i++) {
                buf_l[i] *= pre_gain_f;
                buf_r[i] *= pre_gain_f;
            }
        }

        /* ── Pass 2: IIR biquad cascade + enhancer (surgical) ──── */
        for (size_t i = 0; i < nread; i++) {
            float fl = buf_l[i];
            float fr = buf_r[i];

            /* float [-1,1] → Q16 */
            int32_t ql = (int32_t)(fl * 65536.0f);
            int32_t qr = (int32_t)(fr * 65536.0f);

            /* Clamp to Q16 range */
            if (ql >  65535) ql =  65535;
            if (ql < -65536) ql = -65536;
            if (qr >  65535) qr =  65535;
            if (qr < -65536) qr = -65536;

            BassEnhancer_process_stereo(&enh, &ql, &qr);

            /* Q16 → float */
            buf_l[i] = (float)ql / 65536.0f;
            buf_r[i] = (float)qr / 65536.0f;
        }

        /* ── Interleave output ───────────────────────────────────
         * Apply vol_lut[127] — always at "max volume", but the LUT
         * carries the shelf-dB compensation so midrange stays at the
         * loud reference level.  Identical to ESP32 vol=127 path.   */
        float vol_gain = (float)vol_lut[127] / 65536.0f;
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
