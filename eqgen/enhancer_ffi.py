"""
C FFI wrapper for the harmonic bass enhancer shared library (src/enhancer.so).

Provides opaque-handle creation, destruction, reset, and stereo processing
wrapped behind a clean Python API.  Import this instead of copy-pasting
ctypes boilerplate into every script.
"""

import ctypes
from ctypes import (
    c_int32,
    c_float,
    c_int16,
    c_void_p,
    POINTER,
    Structure,
    byref,
    cast,
)
from pathlib import Path
from typing import List

# ── Library loading (lazy, cached) ────────────────────────────────────

_lib = None


def _load_lib() -> ctypes.CDLL:
    """Locate and load enhancer.so relative to this file's location."""
    so = Path(__file__).resolve().parent.parent / "src" / "enhancer.so"
    if not so.exists():
        raise FileNotFoundError(
            f"enhancer.so not found at {so}.  Build with:  make -C src"
        )
    return ctypes.CDLL(str(so))


def get_lib() -> ctypes.CDLL:
    """Return the cached enhancer shared library handle."""
    global _lib
    if _lib is None:
        _lib = _load_lib()
    return _lib


# ── Parameter struct ──────────────────────────────────────────────────


class EnhancerParams(Structure):
    _fields_ = [
        ("cutoff_hz", c_float),
        ("h2_amp", c_float),
        ("h3_amp", c_float),
        ("release_secs", c_float),
        ("limiter_release_secs", c_float),
        ("fs", c_float),
        ("eq_n_biquads", c_int32),
        ("eq_coeffs_q28", c_void_p),
    ]


# ── Function binding (idempotent) ─────────────────────────────────────

_bound = False
_eq_arr_keepalive = None  # prevent GC of ctypes coeff arrays


def _bind() -> None:
    """Set argument/return types on enhancer functions.  Safe to call repeatedly."""
    global _bound
    if _bound:
        return
    lib = get_lib()
    lib.enhancer_create.argtypes = [POINTER(EnhancerParams)]
    lib.enhancer_create.restype = c_void_p
    lib.enhancer_destroy.argtypes = [c_void_p]
    lib.enhancer_destroy.restype = None
    lib.enhancer_reset.argtypes = [c_void_p]
    lib.enhancer_reset.restype = None
    lib.enhancer_process_stereo.argtypes = [
        c_void_p,
        POINTER(c_int16),
        POINTER(c_int16),
    ]
    lib.enhancer_process_stereo.restype = None
    _bound = True


# ── Convenience API ───────────────────────────────────────────────────


def create_enhancer(
    cutoff_hz: float = 60.0,
    h2_amp: float = 0.33,
    h3_amp: float = 0.33,
    release_secs: float = 0.2,
    limiter_release_secs: float = 0.049,
    fs: float = 44100.0,
    coeffs_q28: List[int] | None = None,
) -> ctypes.c_void_p:
    """Create an enhancer instance.

    Args:
        cutoff_hz: bass crossover frequency (Hz).
        h2_amp, h3_amp: 2nd/3rd harmonic amplitudes (0…1).
        release_secs: envelope release for harmonic linearisation (s).
        limiter_release_secs: AGC limiter release time (s, typ. 0.049).
            The limiter applies 1.0/max(env,1.0) to harmonics only —
            always active, matching the EEL reference behaviour.
        fs: sample rate (Hz).
        coeffs_q28: flat list of Q4.28 biquad coefficients.
    """
    _bind()

    if coeffs_q28 is None:
        coeffs_q28 = []
    n_biquads = len(coeffs_q28) // 5

    eq_arr = (c_int32 * len(coeffs_q28))(*coeffs_q28)
    global _eq_arr_keepalive
    _eq_arr_keepalive = eq_arr  # prevent garbage collection
    params = EnhancerParams(
        cutoff_hz=cutoff_hz,
        h2_amp=h2_amp,
        h3_amp=h3_amp,
        release_secs=release_secs,
        limiter_release_secs=limiter_release_secs,
        fs=fs,
        eq_n_biquads=n_biquads,
        eq_coeffs_q28=cast(eq_arr, c_void_p),
    )
    lib = get_lib()
    enh = lib.enhancer_create(byref(params))
    if not enh:
        raise RuntimeError("enhancer_create returned NULL")
    return enh


def destroy_enhancer(enh: ctypes.c_void_p) -> None:
    """Free an enhancer instance."""
    get_lib().enhancer_destroy(enh)


def reset_enhancer(enh: ctypes.c_void_p) -> None:
    """Reset all filter states in an enhancer instance."""
    get_lib().enhancer_reset(enh)


def process_stereo_frame(enh: ctypes.c_void_p, left: int, right: int) -> tuple[int, int]:
    """Process a single stereo sample pair (int16 values) in-place."""
    l = c_int16(left)
    r = c_int16(right)
    get_lib().enhancer_process_stereo(enh, byref(l), byref(r))
    return l.value, r.value
