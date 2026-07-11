/**
 * fpu.h — Xtensa LX6 FPU denormal flush-to-zero helper
 *
 * The ESP32 Xtensa LX6 FPU has no hardware denormal support.
 * Subnormal floats (< ~1.18e-38) trap to software emulation
 * costing ~300+ cycles instead of single-digit cycles.
 *
 * Audio DSP feedback loops (biquads, envelope releases, LP
 * filters) naturally decay toward zero and will produce
 * denormals if left unchecked.
 *
 * Usage: wrap any stored float that participates in a feedback
 * loop with ftz().
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/** Flush subnormals to zero using integer bit inspection.
 *  Exponent=0 means the float is either zero or subnormal.
 *  Return value is always a normal float or zero — never
 *  subnormal. */
static inline float ftz(float x) {
    union { float f; unsigned u; } u = { .f = x };
    if ((u.u & 0x7F800000u) == 0) return 0.0f;
    return x;
}

#ifdef __cplusplus
}
#endif
