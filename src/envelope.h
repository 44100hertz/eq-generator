/**
 * envelope.h — Q16 envelope follower with reciprocal LUT
 *
 * Peak-hold detector with instant attack and exponential release.
 * The reciprocal LUT allows fast 1/env for normalization.
 *
 * Usage:
 *   Env env;
 *   Env_init(&env, release_coeff_q16);
 *   int32_t env_val = Env_tick(&env, x);
 *   int32_t inv_env = Env_reciprocal(&env, &env->lut, env_val);
 */

#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Reciprocal LUT ────────────────────────────────────────────────── */
/* Maps 8-bit index to 1/index in Q16. Used for fast normalization.     */
/* Generated offline by quantize.py.                                    */

#define ENV_LUT_SIZE 256

typedef struct {
    uint32_t entries[ENV_LUT_SIZE];  /* 1/(index+1) scaled to Q16.16 */
} ReciprocalLUT;

/** Initialize reciprocal LUT. Call once at startup. */
void ReciprocalLUT_init(ReciprocalLUT *lut);

/** Look up 1/x in Q16.16 using the high 8 fractional bits of x.
 *  x must be > 0 (positive int32_t in Q16). */
static inline uint32_t ReciprocalLUT_lookup(const ReciprocalLUT *lut, int32_t x) {
    if (x <= 0) return 0xFFFFFFFF;  /* clamp */
    /* x is Q16: bits 15..8 are the top 8 fractional bits.
       Values >= 65536 (≥1.0) clamp to max index 255. */
    uint32_t ux = (uint32_t)x;
    if (ux >= 65536) return lut->entries[255];
    uint32_t idx = (ux >> 8) & 0xFF;
    return lut->entries[idx];
}

/* ── Envelope follower ─────────────────────────────────────────────── */
typedef struct {
    int32_t peak;       /* current envelope value (Q16) */
    int32_t release;    /* release coefficient (Q16): per-sample multiplier */
                        /* release = float_to_q16(exp(-1/(fs * release_secs))) */
    const ReciprocalLUT *lut;
} Env;

/** Initialize envelope follower.
 *  release_coeff = exp(-1 / (fs * release_time_secs)) in Q16.
 *  Example: for fs=44100, release=200ms:
 *    release_coeff = 0.9998866 * 65536 = 65529
 */
static inline void Env_init(Env *env, const ReciprocalLUT *lut, int32_t release_coeff) {
    env->peak = 0;
    env->release = release_coeff;
    env->lut = lut;
}

/** Process one sample.
 *  Peak-hold with instant attack, exponential release.
 *  x is Q16, output is Q16. */
static inline int32_t Env_tick(Env *env, int32_t x) {
    /* Instant attack: take abs, hold peak */
    int32_t ax = (x >= 0) ? x : -x;
    if (ax > env->peak) {
        env->peak = ax;
    } else {
        /* Exponential release: peak *= release_coeff (scaled) */
        env->peak = (int32_t)(((int64_t)env->peak * (int64_t)env->release) >> 16);
    }
    return env->peak;
}

#ifdef __cplusplus
}
#endif
