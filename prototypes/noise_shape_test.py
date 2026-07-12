#!/usr/bin/env python3
"""
Test error-feedback noise shaping for float→int16 conversion.
Simulates: generate low-frequency sine at very low gain,
truncate to int16, compare with/without leaky integrator error diffusion.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram

FS = 48000
DURATION = 2.0  # seconds
FREQ = 70.0     # Hz — low enough that harmonics ≠ inaudible
N = int(FS * DURATION)

def quantize_truncate(x):
    """Simple truncation (current behavior)."""
    scaled = x * 32768.0
    i = np.floor(scaled).astype(np.int32)
    i = np.clip(i, -32768, 32767)
    return i.astype(np.float64) / 32768.0

def quantize_error_feedback(x, k=0.999):
    """First-order error feedback with leaky integrator.
    Error `err` is in *LSB units* (i.e., int16 domain).
    """
    n = len(x)
    out = np.zeros(n, dtype=np.float64)
    err = 0.0
    for i in range(n):
        # Convert to LSB domain, add leaked previous error
        scaled = x[i] * 32768.0 + k * err
        q = int(np.floor(scaled))
        if -32768 <= q <= 32767:
            clipped = q
        elif q > 32767:
            clipped = 32767
        else:
            clipped = -32768
        # Error in LSB units: what we wanted minus what we got
        err = scaled - float(clipped)
        out[i] = float(clipped) / 32768.0
    return out

def analyze(label, x_original, x_quantized, color):
    """Compute noise spectrum and print stats."""
    noise = x_quantized - x_original
    f, Pxx = periodogram(noise, FS, scaling='density')
    # RMS noise
    rms = np.sqrt(np.mean(noise**2))
    # Noise in the signal band (DC-200Hz)
    mask = f <= 200
    inband = np.sqrt(np.sum(Pxx[mask]) * (f[1] - f[0]))

    print(f"  {label}:")
    print(f"    RMS noise  = {rms:.2e} ( {20*np.log10(rms):.1f} dB FS )")
    print(f"    In-band    = {inband:.2e} ( {20*np.log10(inband):.1f} dB FS )")

    # Average noise floor from 1kHz to Nyquist (where ears are less sensitive)
    mask_hf = f >= 1000
    hf_avg = np.mean(Pxx[mask_hf])
    print(f"    HF floor   = {hf_avg:.2e} ( {10*np.log10(hf_avg):.1f} dB/Hz )")

    return f, Pxx

def main():
    t = np.arange(N) / FS
    # Generate sine at multiple levels
    levels_db = [-60, -80]

    fig, axes = plt.subplots(len(levels_db), 2, figsize=(14, 4 * len(levels_db)))
    if len(levels_db) == 1:
        axes = axes.reshape(1, -1)

    print(f"=== Noise shaping test @ {FREQ} Hz, {FS} Hz sample rate ===\n")

    for row, db in enumerate(levels_db):
        amplitude = 10.0 ** (db / 20.0)
        x = amplitude * np.sin(2.0 * np.pi * FREQ * t)

        x_trunc = quantize_truncate(x)
        x_fb    = quantize_error_feedback(x, k=0.999)

        print(f"--- {db} dB FS ---")
        f, Pxx_trunc = analyze("Truncate",       x, x_trunc, 'tab:red')
        f, Pxx_fb    = analyze("Error feedback", x, x_fb,    'tab:blue')

        # ── Time domain (first 500 samples) ──
        zoom = 500
        ax_t = axes[row, 0]
        ax_t.plot(t[:zoom] * 1000, x_trunc[:zoom], '.-', ms=2, alpha=0.5,
                  color='tab:red', label='truncate')
        ax_t.plot(t[:zoom] * 1000, x_fb[:zoom], '.-', ms=2, alpha=0.5,
                  color='tab:blue', label='error feedback (k=0.999)')
        ax_t.axhline(0, color='gray', lw=0.5)
        ax_t.set_xlabel('Time (ms)')
        ax_t.set_ylabel('Amplitude')
        ax_t.set_title(f'Time domain — {db} dB FS sine')
        ax_t.legend(fontsize=8)

        # ── Frequency domain (noise spectrum) ──
        ax_f = axes[row, 1]
        ax_f.semilogx(f, 10*np.log10(Pxx_trunc + 1e-40), color='tab:red',
                      alpha=0.7, lw=0.8, label='truncate noise')
        ax_f.semilogx(f, 10*np.log10(Pxx_fb + 1e-40), color='tab:blue',
                      alpha=0.7, lw=1.2, label='error feedback noise')
        # Mark the signal band
        ax_f.axvspan(0, 200, alpha=0.1, color='green')
        ax_f.axvline(FREQ, color='orange', lw=0.5, ls='--', label=f'{FREQ} Hz')
        ax_f.set_xlabel('Frequency (Hz)')
        ax_f.set_ylabel('Noise PSD (dB/Hz)')
        ax_f.set_title(f'Noise spectrum — {db} dB FS sine')
        ax_f.legend(fontsize=8)
        ax_f.set_ylim(bottom=-180)

    plt.tight_layout()
    outpath = '/home/cyan/code/eqgen/prototypes/noise_shape_test.png'
    plt.savefig(outpath, dpi=150)
    print(f"\nSaved plot: {outpath}")

    # ── Also test different leak coefficients ──
    print("\n=== Leak coefficient sweep @ -80 dB ===\n")
    amplitude = 10.0 ** (-80.0 / 20.0)
    x = amplitude * np.sin(2.0 * np.pi * FREQ * t)

    for k in [0.9, 0.99, 0.999, 0.9999]:
        x_fb = quantize_error_feedback(x, k=k)
        noise = x_fb - x
        f, Pxx = periodogram(noise, FS, scaling='density')
        rms = np.sqrt(np.mean(noise**2))
        inband = np.sqrt(np.sum(Pxx[f <= 200]) * (f[1] - f[0]))
        hf = np.mean(Pxx[f >= 1000])
        print(f"  k={k:.4f}:  RMS={20*np.log10(rms):.1f} dB  inband={20*np.log10(inband):.1f} dB  HF={10*np.log10(hf):.1f} dB/Hz")


if __name__ == '__main__':
    main()
