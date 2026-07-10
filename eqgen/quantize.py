"""
Fixed-point quantization for Q4.28 biquad coefficients.

The actual DSP processing is done in C (src/enhancer.so via ctypes).
This module provides the offline design/validation tools:
  - BiquadQ28: Q4.28 coefficient quantization (mirrors C implementation)
  - ReciprocalLUT: 1/x lookup table for envelope normalization
  - Utility helpers: float↔Q16 conversion
"""

import numpy as np
from dataclasses import dataclass
from typing import List

from eqgen.eq_fit import BiquadCoeffs


# ─────────────────────────────────────────────────────────────────────────────
# Q16 helpers
# ─────────────────────────────────────────────────────────────────────────────

Q16 = 16
Q16_SCALE = 1 << Q16  # 65536
Q16_MAX = (1 << 31) - 1
Q16_MIN = -(1 << 31)


def float_to_q16(x: float) -> int:
    """Convert float to 16.16 fixed-point with round-to-nearest."""
    scaled = x * Q16_SCALE
    if scaled >= 0:
        result = int(scaled + 0.5)
    else:
        result = int(scaled - 0.5)
    return max(Q16_MIN, min(Q16_MAX, result))


def q16_to_float(q: int) -> float:
    """Convert 16.16 fixed-point back to float."""
    return q / Q16_SCALE


# ─────────────────────────────────────────────────────────────────────────────
# Q4.28 biquad (mirrors src/biquad_q28.h)
# ─────────────────────────────────────────────────────────────────────────────

Q28_SCALE = 1 << 28   # 268435456
Q28_SHIFT = 28


@dataclass
class BiquadQ28:
    """Biquad with Q4.28 coefficients. Mirrors C BiquadQ28 struct.

    Q4.28 range: [-8, 8-2^-28], resolution: 3.7e-9.
    Tick: Q44 accumulator with truncation (no rounding bias).
    """
    b0: int = 0; b1: int = 0; b2: int = 0
    a1: int = 0; a2: int = 0
    x1: int = 0; x2: int = 0
    y1: int = 0; y2: int = 0

    @classmethod
    def from_float(cls, bc: BiquadCoeffs) -> "BiquadQ28":
        """Quantize float coefficients to Q4.28."""
        Q28_INT32_MAX = 0x7FFFFFFF   # 2^31 - 1 = 7.999999996
        Q28_INT32_MIN = -0x80000000  # -2^31    = -8.0

        def q28(v: float) -> int:
            v_clamped = max(-8.0 + 1e-9, min(8.0 - 1e-9, v))
            raw = int(round(v_clamped * Q28_SCALE)) if v_clamped >= 0 \
                else -int(round(-v_clamped * Q28_SCALE))
            # bc.float2q28 returns 0x7FFFFFFF for v >= 7.999999; raw may
            # round up to 0x80000000, which overflows c_int32 → -8.0.
            if raw > Q28_INT32_MAX:
                return Q28_INT32_MAX
            if raw < Q28_INT32_MIN:
                return Q28_INT32_MIN
            return raw
        return cls(b0=q28(bc.b0), b1=q28(bc.b1), b2=q28(bc.b2),
                   a1=q28(bc.a1), a2=q28(bc.a2))

    def reset(self):
        self.x1 = self.x2 = self.y1 = self.y2 = 0

    def tick(self, x: int) -> int:
        """Q16 input → Q16 output. Q44 accumulator, truncation (no rounding)."""
        acc = 0
        acc += int(self.b0) * x
        acc += int(self.b1) * self.x1
        acc += int(self.b2) * self.x2
        acc -= int(self.a1) * self.y1
        acc -= int(self.a2) * self.y2

        y = acc >> Q28_SHIFT
        y = y & 0xFFFFFFFF
        if y >= 0x80000000:
            y -= 0x100000000

        self.x2 = self.x1
        self.x1 = x
        self.y2 = self.y1
        self.y1 = y
        return y


def q28_to_float(v: int) -> float:
    return v / Q28_SCALE


def quantize_biquads_q28(biquads: List[BiquadCoeffs]) -> List[BiquadQ28]:
    """Quantize float biquad list to Q4.28."""
    return [BiquadQ28.from_float(bc) for bc in biquads]


# ─────────────────────────────────────────────────────────────────────────────
# Reciprocal lookup table for envelope normalization
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ReciprocalLUT:
    """Lookup table for fast 1/x in Q16.16.

    Mirrors C ReciprocalLUT in envelope.h.
    Maps 256 entries: lut[i] = round(65536.0 * 256.0 / (i + 1))
    The C lookup uses the top 8 bits of x to index.
    """
    entries: List[int] = None

    @classmethod
    def build(cls) -> "ReciprocalLUT":
        """Build 256-entry LUT matching C ReciprocalLUT_init."""
        entries = []
        for i in range(256):
            # lut[i] = 65536.0 * 256.0 / (i + 1) → Q16.16 value
            entries.append(int(round(65536.0 * 256.0 / (i + 1))))
        return cls(entries=entries)

    def lookup(self, x_q16: int) -> int:
        """Look up 1/x in Q16.16 using top 8 bits of x_q16.
        Matches C ReciprocalLUT_lookup."""
        if x_q16 <= 0:
            return 0x7FFFFFFF
        ux = x_q16 & 0xFFFFFFFF
        # Count leading zeros, then shift to get top 8 bits
        clz = 0
        v = ux
        while v and clz < 32:
            v >>= 1
            clz += 1
        clz = 32 - clz  # count of significant bits
        shift = clz - 8
        if shift >= 0:
            idx = ux >> shift
        else:
            idx = ux << (-shift)
        idx = min(idx, 255)
        return self.entries[idx]

    def max_error(self, n_test: int = 10000) -> float:
        """Return max relative error over Q16 range [256, 32768]."""
        max_err = 0.0
        for i in range(n_test):
            x_q16 = 256 + (32768 - 256) * i // (n_test - 1)
            exact_recip = 65536.0 * 65536.0 / x_q16
            lut_val = self.lookup(x_q16)
            err = abs(lut_val - exact_recip) / exact_recip
            if err > max_err:
                max_err = err
        return max_err
