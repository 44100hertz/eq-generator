#!/usr/bin/env python3
"""Test FFT EQ overlap-add for AM modulation artifacts.

Feeds a 1 kHz sine through the FFT EQ with non-unity gains
and checks for periodic amplitude variation in the output.
"""
import ctypes
import math
import os
import struct
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

# Compile test harness
test_c = SRC.parent / "prototypes" / "fft_eq_test.c"

# Write inline test
test_code = r"""
#include "fft_eq.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    /* Non-trivial gains: gentle tilt, will reveal any OLA issues */
    float gains[FFT_EQ_BINS];
    for (int i = 0; i < FFT_EQ_BINS; i++) {
        /* +6 dB/octave tilt from 0 dB at DC */
        float db = (float)i / (float)(FFT_EQ_BINS - 1) * 24.0f;
        gains[i] = powf(10.0f, db / 20.0f);
    }

    FftEq *eq = fft_eq_create(gains);
    if (!eq) return 1;

    /* Generate 1 kHz sine at 44100 Hz, 2 seconds, float [-1,1] */
    float fs = 44100.0f;
    float freq = 1000.0f;
    int total_samples = (int)(fs * 2.0f);
    int hop = FFT_EQ_HOP;

    float *samples = (float *)malloc((size_t)total_samples * sizeof(float));
    for (int i = 0; i < total_samples; i++) {
        samples[i] = sinf(2.0f * (float)M_PI * freq * (float)i / fs);
    }

    /* Process through FFT EQ */
    float buf_l[FFT_EQ_HOP];
    float buf_r[FFT_EQ_HOP];
    float *output = (float *)malloc((size_t)total_samples * sizeof(float));

    for (int offset = 0; offset < total_samples; offset += hop) {
        int n = hop;
        if (offset + n > total_samples) n = total_samples - offset;

        for (int i = 0; i < n; i++) {
            buf_l[i] = samples[offset + i];
            buf_r[i] = 0.0f;  /* silent right channel */
        }
        /* Pad partial final frame */
        for (int i = n; i < hop; i++) {
            buf_l[i] = 0.0f;
            buf_r[i] = 0.0f;
        }
        fft_eq_process_frame(eq, buf_l, buf_r, buf_l, buf_r);
        for (int i = 0; i < n; i++) {
            output[offset + i] = buf_l[i];
        }
    }

    /* Skip startup transient: first 4 frames (~11.6 ms) */
    int skip = hop * 4;

    /* Compute envelope (absolute value) and check for AM */
    /* AM would show as periodic variation in the running max */
    float max_env = 0.0f, min_env = 1e9f;
    for (int i = skip; i < total_samples - hop; i++) {
        float env = fabsf(output[i]);
        if (env > max_env) max_env = env;
        if (env < min_env) min_env = env;
    }

    /* Running RMS over hop-sized windows to detect periodicity */
    float prev_rms = 0.0f;
    int am_count = 0;
    for (int i = skip; i < total_samples - hop; i += hop) {
        float sum_sq = 0.0f;
        for (int j = 0; j < hop; j++) {
            sum_sq += output[i + j] * output[i + j];
        }
        float rms = sqrtf(sum_sq / (float)hop);
        if (prev_rms > 0.0f) {
            float ratio = rms / prev_rms;
            float ratio_db = 20.0f * log10f(ratio);
            if (fabsf(ratio_db) > 0.5f) am_count++;
        }
        prev_rms = rms;
    }

    int total_windows = (total_samples - skip) / hop;
    float am_pct = (float)am_count / (float)total_windows * 100.0f;

    printf("Env range: %.4f .. %.4f (%.2f dB)\n",
           min_env, max_env, 20.0f * log10f(max_env / min_env));
    printf("AM windows (>0.5 dB jump): %d / %d (%.1f%%)\n",
           am_count, total_windows, am_pct);

    free(samples);
    free(output);
    fft_eq_destroy(eq);

    if (am_pct > 5.0f) {
        printf("FAIL: excessive AM modulation\n");
        return 1;
    }
    printf("PASS: no significant AM modulation\n");
    return 0;
}
"""

test_c_path = SRC.parent / "prototypes" / "fft_eq_test.c"
test_c_path.write_text(test_code)

# Build
result = subprocess.run(
    ["gcc", "-O2", "-Wall", "-Wextra", "-I", str(SRC),
     "-o", str(SRC.parent / "prototypes" / "fft_eq_test"),
     str(test_c_path), str(SRC / "fft_eq.c"), "-lm"],
    capture_output=True, text=True)
if result.returncode != 0:
    print("BUILD FAILED:", result.stderr)
    sys.exit(1)

# Run
result = subprocess.run(
    [str(SRC.parent / "prototypes" / "fft_eq_test")],
    capture_output=True, text=True, timeout=10)
print(result.stdout)
if result.returncode != 0:
    print("TEST FAILED")
    sys.exit(1)
