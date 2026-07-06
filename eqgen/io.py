"""
Audio I/O: WAV loading, measurement loading, and wave utilities.

Canonical versions — imported by end_to_end.py, test scripts, and eqgen.py.
"""

import struct
import wave
from pathlib import Path
from typing import Callable, Tuple

import numpy as np


def read_wav(path: str) -> Tuple[np.ndarray, float]:
    """Read a mono or stereo WAV, returning mono float64 samples in [-1, 1].

    Handles 8/16/24/32-bit integer and 32-bit float formats via scipy.
    """
    from scipy.io import wavfile

    rate, data = wavfile.read(path)

    # Convert to float64 in [-1, 1]
    if data.dtype == np.int16:
        data = data.astype(np.float64) / 32768.0
    elif data.dtype == np.int32:
        data = data.astype(np.float64) / 2147483648.0
    elif data.dtype == np.uint8:
        data = (data.astype(np.float64) - 128.0) / 128.0
    elif data.dtype == np.float32:
        data = data.astype(np.float64)
    else:
        raise ValueError(f"Unsupported WAV dtype: {data.dtype}")

    # Convert to mono
    if data.ndim == 1:
        pass
    elif data.ndim == 2:
        data = data.mean(axis=1)
    else:
        raise ValueError(f"Unexpected WAV shape: {data.shape}")

    return np.ascontiguousarray(data, dtype=np.float64), float(rate)


def load_measurement(
    meas_path: str, target_path: str
) -> Tuple[Callable[[float], float], float]:
    """Load brown-noise measurement and target WAVs.

    Uses stdlib wave + numpy (no scipy dependency).  Computes the speaker's
    frequency response relative to the target recording and returns a
    callable speaker_fn(f_hz) → linear magnitude.

    Returns (speaker_fn, sample_rate).
    """
    with wave.open(str(meas_path), "rb") as mf, wave.open(str(target_path), "rb") as tf:
        fs = mf.getframerate()
        meas_raw = mf.readframes(mf.getnframes())
        target_raw = tf.readframes(tf.getnframes())

    n_meas = len(meas_raw) // 2
    meas = np.array(
        [struct.unpack_from("<h", meas_raw, i * 2)[0] / 32768.0 for i in range(n_meas)],
        dtype=np.float64,
    )

    n_target = len(target_raw) // 2
    target = np.array(
        [struct.unpack_from("<h", target_raw, i * 2)[0] / 32768.0 for i in range(n_target)],
        dtype=np.float64,
    )

    n = min(len(meas), len(target))
    window = np.hanning(n)
    meas_fft = np.abs(np.fft.rfft(meas[:n] * window)) / n
    target_fft = np.abs(np.fft.rfft(target[:n] * window)) / n
    fft_freqs = np.fft.rfftfreq(n, 1.0 / fs)

    eps = np.max(target_fft) * 0.01
    speaker_mag = meas_fft / np.maximum(target_fft, eps)
    speaker_mag /= np.max(speaker_mag[(fft_freqs >= 20) & (fft_freqs <= 200)])

    def speaker_fn(f: float) -> float:
        if f >= fs / 2:
            return 0.0
        return float(np.interp(f, fft_freqs, speaker_mag))

    return speaker_fn, float(fs)
