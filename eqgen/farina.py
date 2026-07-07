"""
Thin wrapper around pyrirtool (https://github.com/maj4e/pyrirtool) for
Farina exponential sine-sweep harmonic distortion analysis.

pyrirtool provides sweep generation + deconvolution.  We add the
harmonic IR extraction and frequency-response helpers on top.
"""

import numpy as np

# pyrirtool (Python 2) uses scipy.signal.tukey — relocated to .windows
# in modern scipy.  Inject the alias before importing stimulus.
import scipy.signal
scipy.signal.tukey = scipy.signal.windows.tukey

from stimulus import stimulus


def generate_sweep(fs, f1, f2, duration, amplitude, silence_start=0.5,
                   silence_end=0.5):
    """Generate an exponential sine sweep and its inverse filter.

    Returns (signal, invfilter, sweep_len) — see pyrirtool.stimulus.
    """
    fs_int = int(fs)
    s = stimulus('sinesweep', fs_int)
    # pyrirtool (Python 2) multiplies float args by fs → float dims that
    # numpy rejects in Python 3.  Pass silence as 0 and pad ourselves.
    s.generate(fs_int, int(duration), amplitude, repetitions=1,
               silenceAtStart=0, silenceAtEnd=0,
               sweeprange=[f1, f2])
    # Manually pad with silence
    pre = np.zeros((int(silence_start * fs_int), 1))
    post = np.zeros((int(silence_end * fs_int), 1))
    signal = np.concatenate([pre, s.signal, post], axis=0)
    total_len = len(signal)
    return signal, s.invfilter, total_len


def deconvolve(output, invfilter, sweep_len, repetitions=1):
    """Deconvolve system output with the inverse filter.

    Returns a 1-D impulse response array.
    """
    s = stimulus('sinesweep', 1)  # fs unused for deconvolve only
    s.invfilter = invfilter
    s.Lp = sweep_len
    s.repetitions = repetitions
    output = np.squeeze(output)
    if output.ndim == 1:
        output = output.reshape(-1, 1)
    rir = s.deconvolve(output)
    return np.squeeze(rir)


def harmonic_ir(rir, fs, f1, f2, duration, order):
    """Extract the windowed impulse response for a given harmonic order.

    In the Farina method the k-th harmonic distortion IR appears
    *before* the linear IR by Δt = T · ln(k) / ln(f₂/f₁).  We locate
    the linear peak, offset back, and apply a 10 ms Hann window.

    Returns (ir_win, peak_idx).
    """
    from math import log as _log

    if f2 <= 0:
        f2 = fs / 2.0
    delta_t = duration * _log(order) / _log(f2 / f1) if order > 1 else 0.0
    delta_samples = int(round(delta_t * fs))

    lin_peak = np.argmax(np.abs(rir))
    h_peak = max(0, lin_peak - delta_samples)

    win_half = int(round(0.005 * fs))
    start = max(0, h_peak - win_half)
    end = min(len(rir), h_peak + win_half + 1)
    win_len = end - start
    win = np.hanning(win_len)

    ir_win = np.zeros_like(rir)
    ir_win[start:end] = rir[start:end] * win
    return ir_win, h_peak


def harmonic_response(ir_win, fs, eval_freqs):
    """Compute the dB magnitude response of a windowed harmonic IR."""
    n_fft = len(ir_win)
    H = np.fft.rfft(ir_win)
    freqs = np.fft.rfftfreq(n_fft, 1.0 / fs)
    mag = np.abs(H)
    mag_interp = np.interp(eval_freqs, freqs, mag)
    with np.errstate(divide='ignore'):
        return 20.0 * np.log10(np.maximum(mag_interp, 1e-20))
