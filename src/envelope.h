/**
 * envelope.h — Float envelope follower (peak-hold with exponential release)
 *
 * Peak-hold detector with instant attack and exponential release.
 *
 * Usage:
 *   Env env;
 *   env_init(&env, release_coeff);
 *   float env_val = env_tick(&env, x);
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/* ── Envelope follower ─────────────────────────────────────────────── */
typedef struct {
    float peak;       /* current envelope value */
    float release;    /* release coefficient: exp(-1/(fs * release_secs)) */
} Env;

/** Initialize envelope follower.
 *  release_coeff = exp(-1 / (fs * release_time_secs)).
 *  Example: for fs=44100, release=200ms:
 *    release_coeff = exp(-1/(44100*0.2)) ≈ 0.9998866
 */
static inline void env_init(Env *env, float release_coeff) {
    env->peak    = 0.0f;
    env->release = release_coeff;
}

/** Process one sample.
 *  Peak-hold with instant attack, exponential release. */
static inline float env_tick(Env *env, float x) {
    float ax = (x >= 0.0f) ? x : -x;
    if (ax > env->peak) {
        env->peak = ax;
    } else {
        env->peak *= env->release;
    }
    return env->peak;
}

#ifdef __cplusplus
}
#endif
