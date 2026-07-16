"""
Audio file processing: decode audio via ffmpeg, process through C enhancer, write WAV.
"""

import os
import struct
import subprocess
from typing import List, Optional

import numpy as np
from scipy.io import wavfile

from eqgen import enhancer_ffi


def process_track(
    input_path: str,
    output_path: str,
    coeffs: List[float],
    n_biquads: int,
    cutoff_hz: float = 60.0,
    h2_amp: float = 0.5,
    h3_amp: float = 1.0,
    pre_gain: float = 1.0,
    start_sec: float = 30.0,
    duration_sec: float = 40.0,
    release_secs: float = 0.2,
) -> bool:
    """Decode audio via ffmpeg, process through C enhancer, write WAV.

    Returns True on success.
    """
    cmd = ["ffmpeg", "-y", "-v", "error",
           "-ss", str(start_sec), "-t", str(duration_sec),
           "-i", input_path,
           "-f", "s16le", "-acodec", "pcm_s16le",
           "-ar", "44100", "-ac", "2", "pipe:1"]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0 or len(result.stdout) < 4:
        print(f"  ffmpeg failed")
        return False

    pcm = result.stdout
    n_frames = len(pcm) // 4
    print(f"  44100Hz stereo, {n_frames} frames ({n_frames/44100:.1f}s)")

    enh = enhancer_ffi.create_enhancer(
        cutoff_hz=cutoff_hz, h2_amp=h2_amp, h3_amp=h3_amp,
        release_secs=release_secs,
        pre_gain=pre_gain,
        fs=44100.0, coeffs=coeffs)

    out_data = bytearray(len(pcm))
    peak_in = 0
    peak_out = 0
    for i in range(0, len(pcm), 4):
        l = struct.unpack_from('<h', pcm, i)[0]
        r = struct.unpack_from('<h', pcm, i+2)[0]
        l_out, r_out = enhancer_ffi.process_stereo_frame(enh, l, r)
        struct.pack_into('<hh', out_data, i, l_out, r_out)
        if abs(l_out) > peak_out: peak_out = abs(l_out)
        if abs(r_out) > peak_out: peak_out = abs(r_out)
        if abs(l) > peak_in: peak_in = abs(l)
        if abs(r) > peak_in: peak_in = abs(r)

    enhancer_ffi.destroy_enhancer(enh)

    int16_data = np.frombuffer(out_data, dtype=np.int16).reshape(-1, 2)
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    wavfile.write(output_path, 44100, int16_data)

    gain_db = 20.0 * np.log10(max(peak_out, 1) / max(peak_in, 1))
    print(f"  Peak: in={peak_in} out={peak_out} ({gain_db:+.1f} dB) \u2192 {output_path}")
    if peak_out >= 32767:
        print("  \u26a0\ufe0f  CLIPPING")
    return True
