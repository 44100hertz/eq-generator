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

#include "eq_coeffs.h"
#include "enhancer.h"

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

    BassEnhancerCfg_init(&cfg,
                         EQGEN_CUTOFF_HZ,
                         EQGEN_H2_AMP,
                         EQGEN_H3_AMP,
                         EQGEN_RELEASE_SECS,
                         fs,
                         EQGEN_LIMITER_RELEASE_SECS,
#ifdef EQGEN_PRE_GAIN_Q16
                         (float)EQGEN_PRE_GAIN_Q16 / 65536.0f,
#else
                         1.0f,
#endif
                         n_bq,
                         coeffs);
    ReciprocalLUT_init(&lut);
    BassEnhancer_init(&enh, &cfg, &lut, eq_bqs_l, eq_bqs_r);

    fprintf(stderr, "filter: running at %.0f Hz\n", fs);

    float buf[512]; /* stereo interleaved: L,R,L,R,...  256 frames = 5.8ms @ 44100 */
    size_t frame_sz = 2 * sizeof(float);

    while (g_running) {
        size_t nread = fread(buf, frame_sz, 256, stdin);
        if (nread == 0) {
            if (feof(stdin)) break;
            if (ferror(stdin)) { perror("filter: read error"); break; }
            continue;
        }

        for (size_t i = 0; i < nread; i++) {
            float fl = buf[i * 2];
            float fr = buf[i * 2 + 1];

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
            buf[i * 2]     = (float)ql / 65536.0f;
            buf[i * 2 + 1] = (float)qr / 65536.0f;
        }

        size_t nwritten = fwrite(buf, frame_sz, nread, stdout);
        if (nwritten < nread) {
            perror("filter: write error");
            break;
        }
        fflush(stdout);
    }

    free(eq_bqs_l);
    free(eq_bqs_r);
    fprintf(stderr, "filter: done\n");
    return 0;
}
