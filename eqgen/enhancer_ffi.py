"""
C FFI wrapper for the harmonic bass enhancer shared library (src/enhancer.so).

Provides opaque-handle creation, destruction, reset, and stereo processing
wrapped behind a clean Python API.  Import this instead of copy-pasting
ctypes boilerplate into every script.

All coefficients are now float (no more Q4.28 integer math).
"""

import ctypes
from ctypes import (
    c_int32,
    c_float,
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


class DspPipeParams(Structure):
    _fields_ = [
        ("cutoff_hz", c_float),
        ("h2_amp", c_float),
        ("h3_amp", c_float),
        ("release_secs", c_float),
        ("push_gain", c_float),
        ("pre_gain", c_float),
        ("fs", c_float),
        ("eq_n_biquads", c_int32),
        ("eq_coeffs", c_void_p),  # pointer to float array
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
    lib.dsp_pipe_create.argtypes = [POINTER(DspPipeParams)]
    lib.dsp_pipe_create.restype = c_void_p
    lib.dsp_pipe_destroy.argtypes = [c_void_p]
    lib.dsp_pipe_destroy.restype = None
    lib.dsp_pipe_handle_reset.argtypes = [c_void_p]
    lib.dsp_pipe_handle_reset.restype = None
    lib.dsp_pipe_handle_process_stereo.argtypes = [
        c_void_p,
        POINTER(c_float),
        POINTER(c_float),
    ]
    lib.dsp_pipe_handle_process_stereo.restype = None
    _bound = True


# ── Convenience API ───────────────────────────────────────────────────


def create_enhancer(
    cutoff_hz: float = 60.0,
    h2_amp: float = 0.33,
    h3_amp: float = 0.33,
    release_secs: float = 0.2,
    push_gain: float = 1.0,
    pre_gain: float = 1.0,
    fs: float = 44100.0,
    coeffs: List[float] | None = None,
) -> ctypes.c_void_p:
    """Create an enhancer instance.

    Args:
        cutoff_hz: bass crossover frequency (Hz).
        h2_amp, h3_amp: 2nd/3rd harmonic amplitudes (0…1).
        release_secs: envelope release for harmonic linearisation (s).
        push_gain: headroom fill strength (1.0 = fill to 0 dBFS).
        pre_gain: linear gain applied before EQ (default 1.0).
        fs: sample rate (Hz).
        coeffs: flat list of float biquad coefficients [b0,b1,b2,a1,a2,...].
    """
    _bind()

    if coeffs is None:
        coeffs = []
    n_biquads = len(coeffs) // 5

    eq_arr = (c_float * len(coeffs))(*coeffs)
    global _eq_arr_keepalive
    _eq_arr_keepalive = eq_arr  # prevent garbage collection
    params = DspPipeParams(
        cutoff_hz=cutoff_hz,
        h2_amp=h2_amp,
        h3_amp=h3_amp,
        release_secs=release_secs,
        push_gain=push_gain,
        pre_gain=pre_gain,
        fs=fs,
        eq_n_biquads=n_biquads,
        eq_coeffs=cast(eq_arr, c_void_p),
    )
    lib = get_lib()
    enh = lib.dsp_pipe_create(byref(params))
    if not enh:
        raise RuntimeError("dsp_pipe_create returned NULL")
    return enh


def destroy_enhancer(enh: ctypes.c_void_p) -> None:
    """Free an enhancer instance."""
    get_lib().dsp_pipe_destroy(enh)


def reset_enhancer(enh: ctypes.c_void_p) -> None:
    """Reset all filter states in an enhancer instance."""
    get_lib().dsp_pipe_handle_reset(enh)


def process_stereo_frame(enh: ctypes.c_void_p, left: float, right: float) -> tuple[float, float]:
    """Process a single stereo sample pair (float values in [-1.0, 1.0]) in-place."""
    l = c_float(left)
    r = c_float(right)
    get_lib().dsp_pipe_handle_process_stereo(enh, byref(l), byref(r))
    return l.value, r.value
