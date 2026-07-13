"""
Deprecated compatibility shim — all DSP is now float.

All functions here are identity passthroughs. Import from eqgen.eq_fit
(BiquadCoeffs) or eqgen.pipeline (pre_gain_from_max_gain) directly
instead. This module will be removed in a future cleanup."""

import numpy as np
from dataclasses import dataclass
from typing import List

from eqgen.eq_fit import BiquadCoeffs


# ─────────────────────────────────────────────────────────────────────────────
# BiquadQ28 — now a float passthrough (backward compatible)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class BiquadQ28:
    """Biquad coefficients — float passthrough (was Q4.28 fixed-point).

    Retained for backward compatibility.  Coefficients are now float values,
    identical to BiquadCoeffs from eq_fit.
    """
    b0: float = 0.0; b1: float = 0.0; b2: float = 0.0
    a1: float = 0.0; a2: float = 0.0

    @classmethod
    def from_float(cls, bc: BiquadCoeffs) -> "BiquadQ28":
        """Passthrough — coefficients are already float."""
        return cls(b0=bc.b0, b1=bc.b1, b2=bc.b2, a1=bc.a1, a2=bc.a2)


def q28_to_float(v: float) -> float:
    """Passthrough — values are already float."""
    return float(v)


def quantize_biquads_q28(biquads: List[BiquadCoeffs]) -> List[BiquadQ28]:
    """Passthrough — biquads stay float."""
    return [BiquadQ28.from_float(bc) for bc in biquads]


# ─────────────────────────────────────────────────────────────────────────────
# Q16 helpers — retained for backward compat (now identity)
# ─────────────────────────────────────────────────────────────────────────────

def float_to_q16(x: float) -> float:
    """Passthrough — float DSP has no Q16."""
    return float(x)


def q16_to_float(q: float) -> float:
    """Passthrough."""
    return float(q)
