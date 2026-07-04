#!/usr/bin/env python3
"""
Analyze the harmonic bass enhancer output for pure sine wave inputs.

Usage:
    python analyze.py                          # run all tests
    python analyze.py --freq 50 --dur 0.5      # single sine test
    python analyze.py --compare                # original vs linear comparison
"""

import numpy as np
import argparse
import sys
from pathlib import Path

# Add parent dir for imports
sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    process_original,
    process_linear,
    generate_test_sine,
    magnitude_spectrum,
    find_peaks,
    describe_harmonics,
)


def analyze_sine(freq: float, amplitude: float, duration: float, cfg: BassEnhancerConfig):
    """Analyze a pure sine through the linear enhancer, printing all spectral peaks."""
    fs = cfg.fs

    # Generate test signal, pad with silence to let filters settle
    pad_samples = int(0.1 * fs)  # 100ms padding
    signal = generate_test_sine([freq], [amplitude], duration, fs, stereo=False)
    signal_2ch = np.zeros((2, signal.shape[0] + 2 * pad_samples))
    signal_2ch[0, pad_samples:pad_samples + signal.shape[0]] = signal
    signal_2ch[1, :] = signal_2ch[0, :] * 0.95

    # Process through both variants
    out_orig = process_original(signal_2ch, cfg)
    out_lin = process_linear(signal_2ch, cfg)

    # Analyze steady-state portion (last half, mono sum)
    start = pad_samples + len(signal) // 2
    end = pad_samples + len(signal)
    mono_orig = (out_orig[0, start:end] + out_orig[1, start:end]) / 2.0
    mono_lin = (out_lin[0, start:end] + out_lin[1, start:end]) / 2.0

    # Also analyze the LP output and harmonic signal (for linear)
    # We'll recompute to extract intermediates
    lpsig, harm_sig, env_sig = extract_linear_intermediates(signal_2ch, cfg)
    mono_lp = lpsig[0, start:end]
    mono_harm = harm_sig[0, start:end]

    # FFT
    fft_size = min(16384, end - start)
    freqs_o, mags_o = magnitude_spectrum(mono_orig, fs, fft_size)
    freqs_l, mags_l = magnitude_spectrum(mono_lin, fs, fft_size)
    freqs_lp, mags_lp = magnitude_spectrum(mono_lp, fs, fft_size)
    freqs_h, mags_h = magnitude_spectrum(mono_harm, fs, fft_size)

    peaks_o = find_peaks(freqs_o, mags_o)
    peaks_l = find_peaks(freqs_l, mags_l)
    peaks_lp = find_peaks(freqs_lp, mags_lp)
    peaks_h = find_peaks(freqs_h, mags_h)

    print(f"\n{'='*70}")
    print(f"  Pure sine: {freq} Hz, amplitude {amplitude} ({20*np.log10(amplitude):+.1f} dBFS)")
    print(f"  Duration: {duration}s, cutoff: {cfg.cutoff} Hz, fs: {fs} Hz")
    print(f"  h2_amp={cfg.h2_amp}, h3_amp={cfg.h3_amp}")
    print(f"{'='*70}")

    # ── LP filter output ───────────────────────────────────────────────
    print(f"\n  ── LP filter output (bass extracted from input) ──")
    print(f"  {'Freq (Hz)':>10s}  {'Magnitude':>10s}  {'dBFS':>8s}")
    print(f"  {'-'*32}")
    for f, m in peaks_lp:
        db = 20 * np.log10(m) if m > 0 else -np.inf
        print(f"  {f:10.2f}  {m:10.6f}  {db:+8.2f}")

    # ── Harmonic signal (before HP) ────────────────────────────────────
    print(f"\n  ── Harmonic signal (raw Chebyshev output, before output HP) ──")
    print(f"  {'Freq (Hz)':>10s}  {'Magnitude':>10s}  {'dBFS':>8s}  {'Note':>30s}")
    print(f"  {'-'*62}")
    for f, m in peaks_h:
        db = 20 * np.log10(m) if m > 0 else -np.inf
        note = ""
        if abs(f - 2 * freq) < 2:
            note = f"← 2nd harmonic ({2*freq:.1f} Hz)"
        elif abs(f - 3 * freq) < 2:
            note = f"← 3rd harmonic ({3*freq:.1f} Hz)"
        elif abs(f - freq) < 2:
            note = f"← FUNDAMENTAL LEAKAGE"
        else:
            note = f"← UNEXPECTED (not harmonic of {freq} Hz)"
        print(f"  {f:10.2f}  {m:10.6f}  {db:+8.2f}  {note}")

    # ── Original plugin output ─────────────────────────────────────────
    print(f"\n  ── ORIGINAL plugin output ──")
    print(f"  {'Freq (Hz)':>10s}  {'Magnitude':>10s}  {'dBFS':>8s}  {'Note':>30s}")
    print(f"  {'-'*62}")
    for f, m in peaks_o:
        db = 20 * np.log10(m) if m > 0 else -np.inf
        note = ""
        if abs(f - freq) < 2:
            note = f"← Fundamental"
        elif abs(f - 2 * freq) < 2:
            note = f"← 2nd harmonic"
        elif abs(f - 3 * freq) < 2:
            note = f"← 3rd harmonic"
        else:
            note = f"← UNEXPECTED"
        print(f"  {f:10.2f}  {m:10.6f}  {db:+8.2f}  {note}")

    # ── LINEAR plugin output ───────────────────────────────────────────
    print(f"\n  ── LINEAR plugin output ──")
    print(f"  {'Freq (Hz)':>10s}  {'Magnitude':>10s}  {'dBFS':>8s}  {'Note':>30s}")
    print(f"  {'-'*62}")
    for f, m in peaks_l:
        db = 20 * np.log10(m) if m > 0 else -np.inf
        note = ""
        if abs(f - freq) < 2:
            note = f"← Fundamental"
        elif abs(f - 2 * freq) < 2:
            note = f"← 2nd harmonic"
        elif abs(f - 3 * freq) < 2:
            note = f"← 3rd harmonic"
        else:
            note = f"← UNEXPECTED"
        print(f"  {f:10.2f}  {m:10.6f}  {db:+8.2f}  {note}")

    # ── Comparison summary ─────────────────────────────────────────────
    print(f"\n  ── Comparison (relative to fundamental) ──")
    print(f"  {'Harmonic':>12s}  {'Original':>10s}  {'Linear':>10s}  {'Delta':>10s}")
    print(f"  {'-'*46}")

    # Find fundamental magnitude in each
    fund_o = next((m for f, m in peaks_o if abs(f - freq) < 2), None)
    fund_l = next((m for f, m in peaks_l if abs(f - freq) < 2), None)

    for h, hf in [(2, 2 * freq), (3, 3 * freq)]:
        m_o = next((m for f, m in peaks_o if abs(f - hf) < 2), None)
        m_l = next((m for f, m in peaks_l if abs(f - hf) < 2), None)
        if fund_o and m_o:
            db_o = 20 * np.log10(m_o / fund_o)
        else:
            db_o = float('-inf')
        if fund_l and m_l:
            db_l = 20 * np.log10(m_l / fund_l)
        else:
            db_l = float('-inf')
        delta = db_l - db_o if db_o > -200 and db_l > -200 else float('nan')
        print(f"  {f'H{h} ({hf:.0f} Hz)':>12s}  {db_o:+9.2f}  {db_l:+9.2f}  {delta:+9.2f}")

    # RMS levels
    rms_o = np.sqrt(np.mean(mono_orig**2))
    rms_l = np.sqrt(np.mean(mono_lin**2))
    print(f"\n  ── RMS levels ──")
    print(f"  Original output RMS:  {rms_o:.6f}  ({20*np.log10(rms_o):+.2f} dBFS)")
    print(f"  Linear output RMS:    {rms_l:.6f}  ({20*np.log10(rms_l):+.2f} dBFS)")
    print(f"  Delta:                {20*np.log10(rms_l/rms_o):+.2f} dB")

    return {
        "freq": freq,
        "amplitude": amplitude,
        "rms_original": rms_o,
        "rms_linear": rms_l,
        "peaks_original": peaks_o,
        "peaks_linear": peaks_l,
        "peaks_harm": peaks_h,
    }


def extract_linear_intermediates(samples, cfg):
    """Run the linear enhancer and return intermediate signals for debugging."""
    import numpy as np
    from harmonic_bass import (
        design_butter_lp,
        design_butter_hp,
        biquad_tick,
        EnvFollower,
    )

    n = samples.shape[1]
    lp_out = np.zeros_like(samples)
    harm_out = np.zeros_like(samples)
    env_out = np.zeros((2, n))

    lp_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    lp_L = np.zeros(4)
    lp_R = np.zeros(4)

    env_L2 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R2 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_L3 = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R3 = EnvFollower.from_params(cfg.env_release, cfg.fs)

    max_corr = 10.0 ** (cfg.corr_limit / 20.0)
    floor = 0.0001

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        r_lp, lp_R = biquad_tick(sR, lp_coeffs, lp_R)
        lp_out[0, i] = l_lp
        lp_out[1, i] = r_lp

        env_L2.tick(l_lp)
        env_R2.tick(r_lp)
        env_L3.tick(l_lp)
        env_R3.tick(r_lp)

        envL2_s = max(env_L2.read(), floor)
        envR2_s = max(env_R2.read(), floor)
        envL3_s = max(env_L3.read(), floor)
        envR3_s = max(env_R3.read(), floor)

        h2_corrL = min(1.0 / envL2_s, max_corr)
        h3_corrL = min(1.0 / (envL3_s * envL3_s), max_corr)
        h2_corrR = min(1.0 / envR2_s, max_corr)
        h3_corrR = min(1.0 / (envR3_s * envR3_s), max_corr)

        harmL = cfg.h2_amp * h2_corrL * (2.0 * l_lp * l_lp - 1.0) \
              + cfg.h3_amp * h3_corrL * (4.0 * l_lp * l_lp * l_lp - 3.0 * l_lp)
        harmR = cfg.h2_amp * h2_corrR * (2.0 * r_lp * r_lp - 1.0) \
              + cfg.h3_amp * h3_corrR * (4.0 * r_lp * r_lp * r_lp - 3.0 * r_lp)

        harm_out[0, i] = harmL
        harm_out[1, i] = harmR
        env_out[0, i] = env_L2.read()
        env_out[1, i] = env_R2.read()

    return lp_out, harm_out, env_out


def sweep_analysis(cfg: BassEnhancerConfig):
    """Analyze frequency response with a slow sine sweep."""
    import time

    fs = cfg.fs
    duration = 2.0
    n = int(duration * fs)
    t = np.arange(n) / fs

    # Log sweep from 25 to 120 Hz
    f_start, f_end = 25.0, 120.0
    phase = 2.0 * np.pi * f_start * (f_end / f_start) ** (t / duration) * t
    sweep_signal = 0.5 * np.sin(phase)

    signal_2ch = np.zeros((2, n))
    signal_2ch[0, :] = sweep_signal
    signal_2ch[1, :] = sweep_signal * 0.95

    print(f"\n{'='*70}")
    print(f"  Sweep analysis: {f_start}-{f_end} Hz, {duration}s")
    print(f"{'='*70}")

    out_orig = process_original(signal_2ch, cfg)
    out_lin = process_linear(signal_2ch, cfg)

    # Short-time RMS
    window_ms = 50
    window_samp = int(window_ms * 0.001 * fs)
    hop = window_samp // 4

    times = []
    rms_orig = []
    rms_lin = []

    for start in range(0, n - window_samp, hop):
        end = start + window_samp
        t_sec = (start + window_samp / 2) / fs
        freq_est = f_start * (f_end / f_start) ** (t_sec / duration)

        mo = (out_orig[0, start:end] + out_orig[1, start:end]) / 2.0
        ml = (out_lin[0, start:end] + out_lin[1, start:end]) / 2.0

        times.append(freq_est)  # use frequency as x-axis
        rms_orig.append(20 * np.log10(np.sqrt(np.mean(mo**2)) + 1e-12))
        rms_lin.append(20 * np.log10(np.sqrt(np.mean(ml**2)) + 1e-12))

    # Summarize
    orig_arr = np.array(rms_orig)
    lin_arr = np.array(rms_lin)

    print(f"\n  ── RMS over sweep ──")
    print(f"  Original:  min={orig_arr.min():+.1f} dBFS  max={orig_arr.max():+.1f} dBFS  mean={orig_arr.mean():+.1f} dBFS")
    print(f"  Linear:    min={lin_arr.min():+.1f} dBFS  max={lin_arr.max():+.1f} dBFS  mean={lin_arr.mean():+.1f} dBFS")
    print(f"  Delta:     mean={np.mean(lin_arr - orig_arr):+.1f} dB  max={np.max(lin_arr - orig_arr):+.1f} dB")

    return {"times": times, "rms_orig": rms_orig, "rms_lin": rms_lin}


def level_analysis(cfg: BassEnhancerConfig):
    """Analyze how output level varies with input amplitude."""
    print(f"\n{'='*70}")
    print(f"  Level linearity analysis: varying input amplitude at 50 Hz")
    print(f"{'='*70}")

    freq = 50.0
    duration = 0.5
    fs = cfg.fs
    amplitudes = np.logspace(-2, 0, 20)  # 0.01 to 1.0

    print(f"\n  {'Input (dBFS)':>12s}  {'Original (dBFS)':>16s}  {'Linear (dBFS)':>16s}  {'Delta (dB)':>12s}")
    print(f"  {'-'*60}")

    results = []
    for amp in amplitudes:
        pad = int(0.2 * fs)
        signal = generate_test_sine([freq], [amp], duration, fs, stereo=False)
        signal_2ch = np.zeros((2, len(signal) + 2 * pad))
        signal_2ch[0, pad:pad + len(signal)] = signal
        signal_2ch[1, :] = signal_2ch[0, :] * 0.95

        out_o = process_original(signal_2ch, cfg)
        out_l = process_linear(signal_2ch, cfg)

        start = pad + len(signal) // 4
        end = pad + len(signal)
        mo = (out_o[0, start:end] + out_o[1, start:end]) / 2.0
        ml = (out_l[0, start:end] + out_l[1, start:end]) / 2.0

        rms_o = np.sqrt(np.mean(mo**2))
        rms_l = np.sqrt(np.mean(ml**2))
        db_in = 20 * np.log10(amp)
        db_o = 20 * np.log10(rms_o) if rms_o > 0 else -200
        db_l = 20 * np.log10(rms_l) if rms_l > 0 else -200
        delta = db_l - db_o

        print(f"  {db_in:+12.2f}  {db_o:+16.2f}  {db_l:+16.2f}  {delta:+12.2f}")
        results.append((db_in, db_o, db_l, delta))

    return results


def main():
    parser = argparse.ArgumentParser(description="Analyze harmonic bass enhancer")
    parser.add_argument("--freq", type=float, default=50.0, help="Test sine frequency (Hz)")
    parser.add_argument("--amp", type=float, default=0.5, help="Test sine amplitude")
    parser.add_argument("--dur", type=float, default=0.5, help="Test sine duration (s)")
    parser.add_argument("--cutoff", type=float, default=60.0, help="Cutoff frequency (Hz)")
    parser.add_argument("--h2", type=float, default=0.13, help="2nd harmonic amplitude")
    parser.add_argument("--h3", type=float, default=0.10, help="3rd harmonic amplitude")
    parser.add_argument("--fs", type=float, default=44100.0, help="Sample rate")
    parser.add_argument("--corr-limit", type=float, default=20.0, help="Max correction (dB)")
    parser.add_argument("--sweep", action="store_true", help="Run sweep analysis")
    parser.add_argument("--levels", action="store_true", help="Run level linearity analysis")
    parser.add_argument("--all", action="store_true", help="Run all analyses")
    args = parser.parse_args()

    cfg = BassEnhancerConfig(
        cutoff=args.cutoff,
        h2_amp=args.h2,
        h3_amp=args.h3,
        fs=args.fs,
        corr_limit=args.corr_limit,
    )

    if args.all or (not args.sweep and not args.levels):
        analyze_sine(args.freq, args.amp, args.dur, cfg)

    if args.all or args.sweep:
        sweep_analysis(cfg)

    if args.all or args.levels:
        level_analysis(cfg)


def run():
    """Run the analysis with default settings."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10, fs=44100.0, corr_limit=20.0,
    )
    analyze_sine(50.0, 0.5, 0.5, cfg)
    sweep_analysis(cfg)
    level_analysis(cfg)


if __name__ == "__main__":
    main()
