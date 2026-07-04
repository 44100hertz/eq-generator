"""
Harmonic purity and linearity tests.

Verifies that the DSP variants produce clean harmonics (no leakage,
no intermodulation) and that harmonic levels are linear with input
amplitude.

Covers:
  - FFT-based single-frequency analysis (peaks, comparison)
  - Sweep response (frequency-domain)
  - Level-linearity sweep (amplitude-domain)
  - Goertzel-based linearity across all variants
  - Variant comparison: original vs linear vs fixed vs fixed_v2
  - T₃ rejection: single-LP vs dual-LP
"""

import numpy as np
from dsp import (
    BassEnhancerConfig,
    process_original,
    process_linear,
    process_fixed,
    process_fixed_v2,
    design_butter_lp,
    design_butter_hp,
    biquad_tick,
    EnvFollower,
)
from analysis import (
    magnitude_spectrum,
    find_peaks,
    generate_test_sine,
    goertzel_magnitude,
    measure_tones,
    format_tone_report,
)
from model import (
    ear_sensitivity,
    butterworth_lp_mag,
    butterworth_hp_mag,
    small_speaker,
)


# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
# Full pipeline: EQ → enhancer → speaker (harmonic transparency)
# ═══════════════════════════════════════════════════════════════════════════════

def pipeline_harmonic_control():
    """Test the full pipeline: EQ preprocessor → enhancer → speaker.

    A small speaker with heavy bass roll-off requires large EQ boosts
    below cutoff.  Those boosted bass frequencies feed the Chebyshev
    polynomials and produce 3rd harmonics at 3f — frequencies where the
    speaker is already flat.  Without a narrower LP on T₃, the 3rd
    harmonic can become louder at the ear than the fundamental itself.

    This test models the complete chain and measures whether excess,
    disconnected 3rd harmonics are controlled.
    """
    fc = 60.0
    h2, h3 = 0.33, 0.33
    cfg = BassEnhancerConfig(
        cutoff=fc, h2_amp=h2, h3_amp=h3,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  FULL PIPELINE: EQ → enhancer → speaker")
    print(f"  Speaker: -12 dB @ 50 Hz, flat > 100 Hz")
    print(f"  Enhancer: cutoff={fc} Hz, h2={h2}, h3={h3}")
    print("=" * 70)

    print(f"\n  EQ gain G solves sqrt(H1²+H2²+H3²) at ear = 0.5")
    print(f"  (matching the Rust bass_enhancer_preprocess logic)")

    freqs = [25, 30, 35, 40, 45, 50, 55, 60, 70, 80]
    target = 0.5

    print(f"\n  {'Freq':>6s}  {'Spkr':>6s}  {'G':>8s}  "
          f"{'Variant':>7s}  {'H1 ear':>8s}  {'H2 ear':>8s}  {'H3 ear':>8s}  "
          f"{'H3/H1':>8s}  {'H3 freq':>9s}")
    print(f"  {'-'*76}")

    results = []
    for freq in freqs:
        Sf  = small_speaker(freq)
        S2f = small_speaker(2 * freq)
        S3f = small_speaker(3 * freq)

        hp_f  = butterworth_hp_mag(freq, fc)
        hp_2f = butterworth_hp_mag(2 * freq, fc)
        hp_3f = butterworth_hp_mag(3 * freq, fc)
        lp_f  = butterworth_lp_mag(freq, fc)
        lp_f2 = butterworth_lp_mag(freq, fc / 2.0)

        a = (hp_f  * Sf)**2
        b = h2**2 * lp_f**2  * hp_2f**2 * S2f**2
        c_fixed = h3**2 * lp_f**2  * hp_3f**2 * S3f**2
        c_v2    = h3**2 * lp_f2**2 * hp_3f**2 * S3f**2

        G_fixed = target / np.sqrt(a + b + c_fixed) if (a + b + c_fixed) > 1e-12 else 1.0
        G_v2    = target / np.sqrt(a + b + c_v2)    if (a + b + c_v2)    > 1e-12 else 1.0

        pad = int(0.3 * cfg.fs)
        sig = generate_test_sine([freq], [1.0], 0.5, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(sig) + 2 * pad))
        stereo[0, pad:pad + len(sig)] = sig
        stereo[1, :] = stereo[0, :] * 0.95
        ss = slice(pad + len(sig) // 2, pad + len(sig))

        for label, proc, G in [("fixed", process_fixed, G_fixed),
                                ("v2",    process_fixed_v2, G_v2)]:
            stereo_gained = stereo * G
            out = proc(stereo_gained, cfg)
            mono = (out[0, ss] + out[1, ss]) / 2.0

            h1_amp = goertzel_magnitude(mono, freq, cfg.fs)
            h2_amp = goertzel_magnitude(mono, 2 * freq, cfg.fs)
            h3_amp = goertzel_magnitude(mono, 3 * freq, cfg.fs)

            h1_ear = h1_amp * Sf
            h2_ear = h2_amp * S2f
            h3_ear = h3_amp * S3f

            h1_db = 20 * np.log10(h1_ear) if h1_ear > 0 else -np.inf
            h2_db = 20 * np.log10(h2_ear) if h2_ear > 0 else -np.inf
            h3_db = 20 * np.log10(h3_ear) if h3_ear > 0 else -np.inf
            h3_rel = 20 * np.log10(h3_ear / max(h1_ear, 1e-12))

            h3_freq = 3 * freq
            region = "ABOVE 2fc" if h3_freq > 2 * fc else "in bass"

            print(f"  {freq:6.0f}  {20*np.log10(Sf):+5.0f}dB  {20*np.log10(G):+7.2f}dB  "
                  f"{label:>7s}  {h1_db:+7.2f}  {h2_db:+7.2f}  {h3_db:+7.2f}  "
                  f"{h3_rel:+7.2f}  {h3_freq:5.0f} Hz {region}")

            results.append({
                "freq": freq, "variant": label,
                "h1_ear": h1_db, "h2_ear": h2_db, "h3_ear": h3_db,
                "h3_rel_h1": h3_rel, "h3_freq": h3_freq,
                "above_2fc": h3_freq > 2 * fc,
            })

        print()

    above = [r for r in results if r["above_2fc"]]
    fixed_above = [r["h3_rel_h1"] for r in above if r["variant"] == "fixed"]
    v2_above = [r["h3_rel_h1"] for r in above if r["variant"] == "v2"]

    print(f"  ── Summary: H3/H1 ratio for H3 above 2fc ──")
    print(f"  fixed:  mean={np.mean(fixed_above):+.1f} dB  "
          f"worst={max(fixed_above):+.1f} dB")
    print(f"  v2:     mean={np.mean(v2_above):+.1f} dB  "
          f"worst={max(v2_above):+.1f} dB")
    improvement = np.mean([f - v for f, v in zip(fixed_above, v2_above)])
    print(f"  v2 reduces H3/H1 by {improvement:.1f} dB on average above 2fc")

    worst_fixed = max([r for r in above if r["variant"] == "fixed"],
                      key=lambda r: r["h3_rel_h1"])
    worst_v2 = [r for r in above if r["variant"] == "v2"
                 and abs(r["freq"] - worst_fixed["freq"]) < 0.1][0]
    print(f"\n  Worst case ({worst_fixed['freq']:.0f} Hz):")
    print(f"    fixed: H3 at {worst_fixed['h3_freq']:.0f} Hz is {worst_fixed['h3_rel_h1']:+.1f} dB rel H1")
    print(f"    v2:    H3 at {worst_v2['h3_freq']:.0f} Hz is {worst_v2['h3_rel_h1']:+.1f} dB rel H1")

    return results


# FFT-based sine analysis (from analyze.py)
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_sine(freq: float, amplitude: float, duration: float, cfg: BassEnhancerConfig):
    """Analyze a pure sine through the linear enhancer, printing all spectral peaks."""
    fs = cfg.fs

    pad_samples = int(0.1 * fs)
    signal = generate_test_sine([freq], [amplitude], duration, fs, stereo=False)
    signal_2ch = np.zeros((2, signal.shape[0] + 2 * pad_samples))
    signal_2ch[0, pad_samples:pad_samples + signal.shape[0]] = signal
    signal_2ch[1, :] = signal_2ch[0, :] * 0.95

    out_orig = process_original(signal_2ch, cfg)
    out_lin = process_linear(signal_2ch, cfg)

    start = pad_samples + len(signal) // 2
    end = pad_samples + len(signal)
    mono_orig = (out_orig[0, start:end] + out_orig[1, start:end]) / 2.0
    mono_lin = (out_lin[0, start:end] + out_lin[1, start:end]) / 2.0

    fft_size = min(16384, end - start)
    freqs_o, mags_o = magnitude_spectrum(mono_orig, fs, fft_size)
    freqs_l, mags_l = magnitude_spectrum(mono_lin, fs, fft_size)

    peaks_o = find_peaks(freqs_o, mags_o)
    peaks_l = find_peaks(freqs_l, mags_l)

    print(f"\n{'='*70}")
    print(f"  Pure sine: {freq} Hz, amplitude {amplitude} ({20*np.log10(amplitude):+.1f} dBFS)")
    print(f"  Duration: {duration}s, cutoff: {cfg.cutoff} Hz, fs: {fs} Hz")
    print(f"  h2_amp={cfg.h2_amp}, h3_amp={cfg.h3_amp}")
    print(f"{'='*70}")

    for label, peaks in [("ORIGINAL plugin output", peaks_o), ("LINEAR plugin output", peaks_l)]:
        print(f"\n  ── {label} ──")
        print(f"  {'Freq (Hz)':>10s}  {'Magnitude':>10s}  {'dBFS':>8s}  {'Note':>30s}")
        print(f"  {'-'*62}")
        for f, m in peaks:
            db = 20 * np.log10(m) if m > 0 else -np.inf
            note = ""
            if abs(f - freq) < 2:
                note = "← Fundamental"
            elif abs(f - 2 * freq) < 2:
                note = "← 2nd harmonic"
            elif abs(f - 3 * freq) < 2:
                note = "← 3rd harmonic"
            else:
                note = "← UNEXPECTED"
            print(f"  {f:10.2f}  {m:10.6f}  {db:+8.2f}  {note}")

    # Comparison summary
    print(f"\n  ── Comparison (relative to fundamental) ──")
    print(f"  {'Harmonic':>12s}  {'Original':>10s}  {'Linear':>10s}")
    print(f"  {'-'*36}")

    fund_o = next((m for f, m in peaks_o if abs(f - freq) < 2), None)
    fund_l = next((m for f, m in peaks_l if abs(f - freq) < 2), None)

    for h, hf in [(2, 2 * freq), (3, 3 * freq)]:
        m_o = next((m for f, m in peaks_o if abs(f - hf) < 2), None)
        m_l = next((m for f, m in peaks_l if abs(f - hf) < 2), None)
        db_o = 20 * np.log10(m_o / fund_o) if fund_o and m_o else float('-inf')
        db_l = 20 * np.log10(m_l / fund_l) if fund_l and m_l else float('-inf')
        print(f"  {f'H{h} ({hf:.0f} Hz)':>12s}  {db_o:+9.2f}  {db_l:+9.2f}")


def sweep_analysis(cfg: BassEnhancerConfig):
    """Analyze frequency response with a slow sine sweep."""
    fs = cfg.fs
    duration = 2.0
    n = int(duration * fs)
    t = np.arange(n) / fs

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

    window_ms = 50
    window_samp = int(window_ms * 0.001 * fs)
    hop = window_samp // 4

    rms_orig = []
    rms_lin = []

    for start in range(0, n - window_samp, hop):
        end = start + window_samp
        mo = (out_orig[0, start:end] + out_orig[1, start:end]) / 2.0
        ml = (out_lin[0, start:end] + out_lin[1, start:end]) / 2.0
        rms_orig.append(20 * np.log10(np.sqrt(np.mean(mo**2)) + 1e-12))
        rms_lin.append(20 * np.log10(np.sqrt(np.mean(ml**2)) + 1e-12))

    oa = np.array(rms_orig)
    la = np.array(rms_lin)
    print(f"\n  ── RMS over sweep ──")
    print(f"  Original:  min={oa.min():+.1f} dBFS  max={oa.max():+.1f} dBFS  mean={oa.mean():+.1f} dBFS")
    print(f"  Linear:    min={la.min():+.1f} dBFS  max={la.max():+.1f} dBFS  mean={la.mean():+.1f} dBFS")
    print(f"  Delta:     mean={np.mean(la - oa):+.1f} dB  max={np.max(la - oa):+.1f} dB")


def level_analysis(cfg: BassEnhancerConfig):
    """Analyze how output level varies with input amplitude."""
    print(f"\n{'='*70}")
    print(f"  Level linearity analysis: varying input amplitude at 50 Hz")
    print(f"{'='*70}")

    freq = 50.0
    duration = 0.5
    fs = cfg.fs
    amplitudes = np.logspace(-2, 0, 20)

    print(f"\n  {'Input (dBFS)':>12s}  {'Original (dBFS)':>16s}  {'Linear (dBFS)':>16s}  {'Delta (dB)':>12s}")
    print(f"  {'-'*60}")

    for amp in amplitudes:
        pad = int(0.2 * fs)
        signal = generate_test_sine([freq], [amp], duration, fs, stereo=False)
        signal_2ch = np.zeros((2, len(signal) + 2 * pad))
        signal_2ch[0, pad:pad + len(signal)] = signal
        signal_2ch[1, :] = signal_2ch[0, :] * 0.95

        out_o = process_original(signal_2ch, cfg)
        out_l = process_linear(signal_2ch, cfg)

        ss_slice = slice(pad + len(signal) // 4, pad + len(signal))
        mo = (out_o[0, ss_slice] + out_o[1, ss_slice]) / 2.0
        ml = (out_l[0, ss_slice] + out_l[1, ss_slice]) / 2.0

        rms_o = np.sqrt(np.mean(mo**2))
        rms_l = np.sqrt(np.mean(ml**2))
        db_in = 20 * np.log10(amp)
        db_o = 20 * np.log10(rms_o) if rms_o > 0 else -200
        db_l = 20 * np.log10(rms_l) if rms_l > 0 else -200
        delta = db_l - db_o

        print(f"  {db_in:+12.2f}  {db_o:+16.2f}  {db_l:+16.2f}  {delta:+12.2f}")


# ═══════════════════════════════════════════════════════════════════════════════
# Goertzel-based linearity (from linearity_verify.py)
# ═══════════════════════════════════════════════════════════════════════════════

def analyze_variant(label, process_fn, freq, cfg, amplitudes=None):
    """Full linearity sweep for one variant using Goertzel analysis."""
    if amplitudes is None:
        amplitudes = [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.125, 0.09]

    duration = 0.5
    pad = int(0.3 * cfg.fs)

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  f={freq} Hz, cutoff={cfg.cutoff} Hz, h2={cfg.h2_amp}, h3={cfg.h3_amp}")
    print(f"{'='*70}")
    print(f"\n  {'In dBFS':>8s}  {'Fund dBFS':>10s}  "
          f"{'H2 rel':>8s}  {'H3 rel':>8s}  "
          f"{'Noise dB':>9s}  {'SNR dB':>8s}")
    print(f"  {'-'*56}")

    h2_rel_list, h3_rel_list, snr_list = [], [], []

    for amp in amplitudes:
        signal = generate_test_sine([freq], [amp], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out = process_fn(stereo, cfg)
        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))
        mono = (out[0, ss_slice] + out[1, ss_slice]) / 2.0

        m = measure_tones(mono, freq, cfg.fs, harmonics=(1, 2, 3))
        tones = m["tones"]

        fund_db = 20 * np.log10(max(tones.get(1, 0), 1e-12))
        h2_rel = 20 * np.log10(max(tones.get(2, 0), 1e-12) / max(tones.get(1, 0), 1e-12))
        h3_rel = 20 * np.log10(max(tones.get(3, 0), 1e-12) / max(tones.get(1, 0), 1e-12))
        noise_db = 20 * np.log10(max(m["noise_rms"], 1e-12))
        snr = m.get("snr_db", -np.inf)

        h2_rel_list.append(h2_rel)
        h3_rel_list.append(h3_rel)
        snr_list.append(snr)

        print(f"  {20*np.log10(amp):+8.2f}  {fund_db:+10.2f}  "
              f"{h2_rel:+8.2f}  {h3_rel:+8.2f}  "
              f"{noise_db:+9.2f}  {snr:+8.1f}")

    valid_h2 = [v for v in h2_rel_list if v > -200]
    valid_h3 = [v for v in h3_rel_list if v > -200]

    if valid_h2:
        h2_mean, h2_std, h2_range = np.mean(valid_h2), np.std(valid_h2), max(valid_h2) - min(valid_h2)
        linear = "✅ LINEAR" if h2_range < 1.5 else f"❌ NONLINEAR (±{h2_range:.1f} dB spread)"
        print(f"\n  H2:  mean={h2_mean:+.2f} dB  σ={h2_std:.2f} dB  range={h2_range:.2f} dB  {linear}")

    if valid_h3:
        h3_mean, h3_std, h3_range = np.mean(valid_h3), np.std(valid_h3), max(valid_h3) - min(valid_h3)
        linear = "✅ LINEAR" if h3_range < 1.5 else f"❌ NONLINEAR (±{h3_range:.1f} dB spread)"
        print(f"  H3:  mean={h3_mean:+.2f} dB  σ={h3_std:.2f} dB  range={h3_range:.2f} dB  {linear}")

    print(f"  SNR: mean={np.mean(snr_list):+.1f} dB  "
          f"min={min(snr_list):+.1f} dB  (higher = cleaner)")


# ═══════════════════════════════════════════════════════════════════════════════
# Variant comparison (from fix_verify.py)
# ═══════════════════════════════════════════════════════════════════════════════

def compare_all_variants(freq=50.0, amp=0.5, duration=0.3):
    """Compare original, linear, and fixed variants on a pure sine."""
    cfg = BassEnhancerConfig(cutoff=60.0, h2_amp=0.13, h3_amp=0.10, fs=44100.0)

    pad = int(0.2 * cfg.fs)
    signal = generate_test_sine([freq], [amp], duration, cfg.fs, stereo=False)
    stereo = np.zeros((2, len(signal) + 2 * pad))
    stereo[0, pad:pad + len(signal)] = signal
    stereo[1, :] = stereo[0, :] * 0.95

    out_orig = process_original(stereo, cfg)
    out_lin = process_linear(stereo, cfg)
    out_fix = process_fixed(stereo, cfg)

    start = pad + len(signal) // 2
    end = pad + len(signal)
    fft_size = min(16384, end - start)

    def get_mono(out):
        return (out[0, start:end] + out[1, start:end]) / 2.0

    results = {}
    for name, out in [("original", out_orig), ("linear", out_lin), ("fixed", out_fix)]:
        mono = get_mono(out)
        rms = np.sqrt(np.mean(mono**2))
        freqs, mags = magnitude_spectrum(mono, cfg.fs, fft_size)
        peaks = find_peaks(freqs, mags)
        results[name] = {"rms": rms, "peaks": peaks}

    print(f"\n{'='*70}")
    print(f"  VARIANT COMPARISON: {freq} Hz sine, A={amp} ({20*np.log10(amp):+.1f} dBFS)")
    print(f"  h2_amp={cfg.h2_amp}, h3_amp={cfg.h3_amp}, cutoff={cfg.cutoff} Hz")
    print(f"{'='*70}")

    print(f"\n  {'Freq (Hz)':>10s}  {'Original':>10s}  {'Linear':>10s}  {'Fixed':>10s}  {'Note':>30s}")
    print(f"  {'-'*72}")

    all_freqs = sorted(set(int(round(p[0])) for r in results.values() for p in r["peaks"]))
    for f_round in all_freqs:
        vals = {}
        for name in ["original", "linear", "fixed"]:
            match = next((m for f, m in results[name]["peaks"] if abs(f - f_round) < 2), None)
            vals[name] = 20 * np.log10(match) if match else float('-inf')

        note = ""
        if abs(f_round - freq) < 3:
            note = "← Fundamental"
        elif abs(f_round - 2 * freq) < 3:
            note = f"← 2nd harmonic ({2*freq:.0f} Hz)"
        elif abs(f_round - 3 * freq) < 3:
            note = f"← 3rd harmonic ({3*freq:.0f} Hz)"
        else:
            note = "← UNEXPECTED"

        def fmt(v):
            return f"{v:+8.2f}" if v > -200 else "       --"

        print(f"  {f_round:10.0f}  {fmt(vals['original']):>10s}  {fmt(vals['linear']):>10s}  {fmt(vals['fixed']):>10s}  {note}")

    print(f"\n  ── RMS levels ──")
    for name in ["original", "linear", "fixed"]:
        rms = results[name]["rms"]
        print(f"  {name:>10s}: {rms:.6f}  ({20*np.log10(rms):+.2f} dBFS)")

    f_peaks = results["fixed"]["peaks"]
    fund_peak = next((m for f, m in f_peaks if abs(f - freq) < 2), None)
    if fund_peak:
        unexpected = [(f, m) for f, m in f_peaks
                      if abs(f - freq) > 3 and abs(f - 2 * freq) > 3 and abs(f - 3 * freq) > 3]
        if unexpected:
            print(f"  ⚠️  Unexpected peaks: {[(round(f,1), f'{20*np.log10(m/fund_peak):+.1f} dB') for f,m in unexpected]}")
        else:
            print(f"  ✅ No unexpected peaks — clean harmonic generation!")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# T₃ rejection comparison (from fix_v2.py)
# ═══════════════════════════════════════════════════════════════════════════════

def compare_t3_rejection():
    """Compare T₃ output between single-LP (fixed) and dual-LP (fixed_v2)."""
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )

    print("=" * 70)
    print("  T₃ REJECTION ABOVE CUTOFF: fixed (single LP) vs fixed_v2 (fc/2 LP)")
    print("=" * 70)

    for freq in [50, 60, 80, 100, 120, 180, 300]:
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [0.5], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out_v1 = process_fixed(stereo, cfg)
        out_v2 = process_fixed_v2(stereo, cfg)

        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))
        mono_v1 = (out_v1[0, ss_slice] + out_v1[1, ss_slice]) / 2.0
        mono_v2 = (out_v2[0, ss_slice] + out_v2[1, ss_slice]) / 2.0

        h3_v1 = goertzel_magnitude(mono_v1, 3 * freq, cfg.fs)
        h3_v2 = goertzel_magnitude(mono_v2, 3 * freq, cfg.fs)
        fund_v1 = goertzel_magnitude(mono_v1, freq, cfg.fs)
        fund_v2 = goertzel_magnitude(mono_v2, freq, cfg.fs)

        h3_rel_v1 = 20 * np.log10(max(h3_v1, 1e-12) / max(fund_v1, 1e-12))
        h3_rel_v2 = 20 * np.log10(max(h3_v2, 1e-12) / max(fund_v2, 1e-12))
        delta = h3_rel_v2 - h3_rel_v1

        rms_v1 = np.sqrt(np.mean(mono_v1**2))
        rms_v2 = np.sqrt(np.mean(mono_v2**2))
        rms_delta = 20 * np.log10(rms_v2 / rms_v1) if rms_v1 > 0 else 0

        status = "✅ quieter" if delta < -1 else ("⚠️ same" if abs(delta) < 2 else "❌ LOUDER")
        print(f"  {freq:3d} Hz:  H3 rel: fixed={h3_rel_v1:+.1f} dB  v2={h3_rel_v2:+.1f} dB  "
              f"Δ={delta:+.1f} dB  {status}  (RMS Δ={rms_delta:+.1f} dB)")

    # Detailed comparison at 180 Hz
    print(f"\n  ── Detailed T₃ output at 180 Hz ──")
    freq = 180.0
    pad = int(0.3 * cfg.fs)
    signal = generate_test_sine([freq], [0.5], 0.5, cfg.fs, stereo=False)
    stereo = np.zeros((2, len(signal) + 2 * pad))
    stereo[0, pad:pad + len(signal)] = signal
    stereo[1, :] = stereo[0, :] * 0.95

    for label, fn in [("fixed (single LP at fc)", process_fixed),
                       ("fixed_v2 (separate fc/2 LP)", process_fixed_v2)]:
        out = fn(stereo, cfg)
        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))
        mono = (out[0, ss_slice] + out[1, ss_slice]) / 2.0
        m = measure_tones(mono, freq, cfg.fs, harmonics=(1, 2, 3))
        format_tone_report(m, label, fundamental=freq)


# ═══════════════════════════════════════════════════════════════════════════════
# Below-cutoff sweep: fixed vs fixed_v2
# ═══════════════════════════════════════════════════════════════════════════════

def below_cutoff_sweep():
    """Sweep through bass frequencies, comparing H2/H3 output from both variants.

    Quantifies the tradeoff: v2's narrower fc/2 LP for T₃ means less 3rd harmonic
    below cutoff. How much do we lose?
    """
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )
    amplitude = 0.5
    freqs = [25, 30, 35, 40, 45, 50, 55, 60]

    print("=" * 70)
    print("  BELOW-CUTOFF SWEEP: fixed vs fixed_v2")
    print(f"  cutoff={cfg.cutoff} Hz, h2={cfg.h2_amp}, h3={cfg.h3_amp}, A={amplitude}")
    print("=" * 70)

    # Header
    print(f"\n  {'Freq':>6s}  "
          f"{'H1 fixed':>10s}  {'H2 fixed':>10s}  {'H3 fixed':>10s}  "
          f"{'H1 v2':>10s}  {'H2 v2':>10s}  {'H3 v2':>10s}  "
          f"{'ΔH2':>8s}  {'ΔH3':>8s}")
    print(f"  {'-'*90}")

    results = []
    for freq in freqs:
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [amplitude], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out_fixed = process_fixed(stereo, cfg)
        out_v2 = process_fixed_v2(stereo, cfg)

        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))

        m_fixed = measure_tones(
            (out_fixed[0, ss_slice] + out_fixed[1, ss_slice]) / 2.0,
            freq, cfg.fs, harmonics=(1, 2, 3),
        )
        m_v2 = measure_tones(
            (out_v2[0, ss_slice] + out_v2[1, ss_slice]) / 2.0,
            freq, cfg.fs, harmonics=(1, 2, 3),
        )

        t_f = m_fixed["tones"]
        t_v = m_v2["tones"]

        fund_f = max(t_f.get(1, 0), 1e-12)
        fund_v = max(t_v.get(1, 0), 1e-12)

        def rel(tones, h):
            amp = max(tones.get(h, 0), 1e-12)
            return amp / max(tones.get(1, 0), 1e-12)

        h2_f_rel = 20 * np.log10(rel(t_f, 2))
        h3_f_rel = 20 * np.log10(rel(t_f, 3))
        h2_v_rel = 20 * np.log10(rel(t_v, 2))
        h3_v_rel = 20 * np.log10(rel(t_v, 3))

        dh2 = h2_v_rel - h2_f_rel
        dh3 = h3_v_rel - h3_f_rel

        print(f"  {freq:6.0f}  "
              f"{20*np.log10(t_f[1]):+10.2f}  {h2_f_rel:+10.2f}  {h3_f_rel:+10.2f}  "
              f"{20*np.log10(t_v[1]):+10.2f}  {h2_v_rel:+10.2f}  {h3_v_rel:+10.2f}  "
              f"{dh2:+8.2f}  {dh3:+8.2f}")

        results.append({
            "freq": freq,
            "h2_fixed": h2_f_rel, "h3_fixed": h3_f_rel,
            "h2_v2": h2_v_rel, "h3_v2": h3_v_rel,
            "dh2": dh2, "dh3": dh3,
            "fund_fixed": t_f[1], "fund_v2": t_v[1],
        })

    # ── Summary ───────────────────────────────────────────────────────
    print(f"\n  ── Summary ──")
    dh2_vals = [r["dh2"] for r in results]
    dh3_vals = [r["dh3"] for r in results]
    print(f"  ΔH2 (v2 − fixed): mean={np.mean(dh2_vals):+.2f} dB  "
          f"min={min(dh2_vals):+.2f} dB  max={max(dh2_vals):+.2f} dB")
    print(f"  ΔH3 (v2 − fixed): mean={np.mean(dh3_vals):+.2f} dB  "
          f"min={min(dh3_vals):+.2f} dB  max={max(dh3_vals):+.2f} dB")

    # Show where the T₃ LP (fc/2) starts cutting in
    print(f"\n  T₃ LP at fc/2={cfg.cutoff/2:.0f} Hz starts attenuating above ~{cfg.cutoff/2:.0f} Hz.")
    print(f"  This reduces H3 level in v2 vs fixed, increasingly above that point.")
    worst = min(results, key=lambda r: r["dh3"])
    print(f"  Worst-case H3 loss: {worst['dh3']:+.2f} dB at {worst['freq']:.0f} Hz")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Loudness-weighted comparison
# ═══════════════════════════════════════════════════════════════════════════════

def loudness_weighted_comparison():
    """Compare fixed vs fixed_v2 using equal-loudness-weighted perceived level.

    The ear is much more sensitive at 2f/3f than at f for bass frequencies.
    At 40 Hz the ear is 4.6× more sensitive at H2 (80 Hz) and 8.3× more
    sensitive at H3 (120 Hz) than at the fundamental.

    A flat-RMS model ignores this. The loudness-weighted model tells us
    whether losing H3 content actually matters perceptually.
    """
    cfg = BassEnhancerConfig(
        cutoff=60.0, h2_amp=0.13, h3_amp=0.10,
        fs=44100.0, env_release=200.0,
    )
    amplitude = 0.5
    freqs = [25, 30, 35, 40, 45, 50, 55, 60]

    print("=" * 70)
    print("  LOUDNESS-WEIGHTED COMPARISON: fixed vs fixed_v2")
    print(f"  cutoff={cfg.cutoff} Hz, h2={cfg.h2_amp}, h3={cfg.h3_amp}, A={amplitude}")
    print("=" * 70)

    # Show ear sensitivity at relevant frequencies
    print(f"\n  ── Ear sensitivity (ISO 226:2023 ≈60 phon) ──")
    print(f"  {'Freq':>6s}  {'s(f)':>8s}  {'s(2f)':>8s}  {'s(3f)':>8s}  "
          f"{'s(2f)/s(f)':>12s}  {'s(3f)/s(f)':>12s}")
    print(f"  {'-'*60}")
    for f in freqs:
        sf = ear_sensitivity(f)
        s2f = ear_sensitivity(2 * f)
        s3f = ear_sensitivity(3 * f)
        print(f"  {f:5.0f} Hz  {sf:8.3f}  {s2f:8.3f}  {s3f:8.3f}  "
              f"{s2f/sf:12.1f}×  {s3f/sf:12.1f}×")

    # ── Per-frequency comparison ──────────────────────────────────────
    print(f"\n  ── Per-frequency loudness comparison ──")
    print(f"  {'Freq':>6s}  "
          f"{'RMS fixed':>10s}  {'RMS v2':>10s}  {'ΔRMS':>8s}  "
          f"{'Loud fixed':>10s}  {'Loud v2':>10s}  {'ΔLoud':>8s}")
    print(f"  {'-'*68}")

    results = []
    for freq in freqs:
        duration = 0.5
        pad = int(0.3 * cfg.fs)
        signal = generate_test_sine([freq], [amplitude], duration, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95

        out_fixed = process_fixed(stereo, cfg)
        out_v2 = process_fixed_v2(stereo, cfg)

        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))
        mono_f = (out_fixed[0, ss_slice] + out_fixed[1, ss_slice]) / 2.0
        mono_v = (out_v2[0, ss_slice] + out_v2[1, ss_slice]) / 2.0

        # Physical RMS
        rms_f = 20 * np.log10(np.sqrt(np.mean(mono_f**2)))
        rms_v = 20 * np.log10(np.sqrt(np.mean(mono_v**2)))
        drms = rms_v - rms_f

        # Goertzel amplitudes for each harmonic
        h1_f = goertzel_magnitude(mono_f, freq, cfg.fs)
        h2_f = goertzel_magnitude(mono_f, 2 * freq, cfg.fs)
        h3_f = goertzel_magnitude(mono_f, 3 * freq, cfg.fs)
        h1_v = goertzel_magnitude(mono_v, freq, cfg.fs)
        h2_v = goertzel_magnitude(mono_v, 2 * freq, cfg.fs)
        h3_v = goertzel_magnitude(mono_v, 3 * freq, cfg.fs)

        # Ear sensitivity at f, 2f, 3f
        sf = ear_sensitivity(freq)
        s2f = ear_sensitivity(2 * freq)
        s3f = ear_sensitivity(3 * freq)

        # Perceived loudness: sqrt(Σ s² · A²)
        loud_f = np.sqrt((sf * h1_f)**2 + (s2f * h2_f)**2 + (s3f * h3_f)**2)
        loud_v = np.sqrt((sf * h1_v)**2 + (s2f * h2_v)**2 + (s3f * h3_v)**2)
        loud_f_db = 20 * np.log10(loud_f) if loud_f > 0 else -np.inf
        loud_v_db = 20 * np.log10(loud_v) if loud_v > 0 else -np.inf
        dloud = loud_v_db - loud_f_db

        print(f"  {freq:6.0f}  "
              f"{rms_f:+10.2f}  {rms_v:+10.2f}  {drms:+8.2f}  "
              f"{loud_f_db:+10.2f}  {loud_v_db:+10.2f}  {dloud:+8.2f}")

        results.append({
            "freq": freq,
            "rms_fixed": rms_f, "rms_v2": rms_v, "drms": drms,
            "loud_fixed": loud_f_db, "loud_v2": loud_v_db, "dloud": dloud,
            "h1_f": h1_f, "h2_f": h2_f, "h3_f": h3_f,
            "h1_v": h1_v, "h2_v": h2_v, "h3_v": h3_v,
        })

    # ── Harmonic contribution breakdown ───────────────────────────────
    print(f"\n  ── Harmonic contribution to perceived loudness ──")
    print(f"  {'Freq':>6s}  "
          f"{'H1% fixed':>10s}  {'H2% fixed':>10s}  {'H3% fixed':>10s}  "
          f"{'H1% v2':>10s}  {'H2% v2':>10s}  {'H3% v2':>10s}")
    print(f"  {'-'*68}")

    for r in results:
        freq = r["freq"]
        sf = ear_sensitivity(freq)
        s2f = ear_sensitivity(2 * freq)
        s3f = ear_sensitivity(3 * freq)

        def contributions(h1, h2, h3):
            w1 = (sf * h1)**2
            w2 = (s2f * h2)**2
            w3 = (s3f * h3)**2
            total = w1 + w2 + w3
            if total == 0:
                return 0, 0, 0
            return 100 * w1 / total, 100 * w2 / total, 100 * w3 / total

        c1f, c2f, c3f = contributions(r["h1_f"], r["h2_f"], r["h3_f"])
        c1v, c2v, c3v = contributions(r["h1_v"], r["h2_v"], r["h3_v"])

        print(f"  {freq:6.0f}  "
              f"{c1f:9.1f}%  {c2f:9.1f}%  {c3f:9.1f}%  "
              f"{c1v:9.1f}%  {c2v:9.1f}%  {c3v:9.1f}%")

    # ── Summary ───────────────────────────────────────────────────────
    dloud_vals = [r["dloud"] for r in results]
    drms_vals = [r["drms"] for r in results]
    print(f"\n  ── Summary ──")
    print(f"  ΔRMS  (v2 − fixed): mean={np.mean(drms_vals):+.2f} dB  "
          f"min={min(drms_vals):+.2f} dB  max={max(drms_vals):+.2f} dB")
    print(f"  ΔLoud (v2 − fixed): mean={np.mean(dloud_vals):+.2f} dB  "
          f"min={min(dloud_vals):+.2f} dB  max={max(dloud_vals):+.2f} dB")

    # Explain the difference
    avg_gap = np.mean(dloud_vals) - np.mean(drms_vals)
    print(f"\n  Flat-RMS says v2 is {'quieter' if np.mean(drms_vals) < 0 else 'louder'} by {abs(np.mean(drms_vals)):.2f} dB.")
    print(f"  Loudness-weighted says v2 is {'quieter' if np.mean(dloud_vals) < 0 else 'louder'} by {abs(np.mean(dloud_vals)):.2f} dB.")
    if abs(avg_gap) > 0.1:
        direction = "less" if avg_gap > 0 else "more"
        print(f"  → The loudness penalty is {abs(avg_gap):.2f} dB {direction} than RMS suggests,")
        print(f"    because harmonics contribute disproportionately to perceived loudness at bass frequencies.")
    else:
        print(f"  → Flat-RMS and loudness-weighted agree — the ear weighting doesn't change the story.")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all harmonic and linearity tests."""
    cfg_old = BassEnhancerConfig(cutoff=60.0, h2_amp=0.13, h3_amp=0.10, fs=44100.0,
                                  corr_limit=20.0)
    cfg_v2 = BassEnhancerConfig(cutoff=60.0, h2_amp=0.33, h3_amp=0.33, fs=44100.0,
                                 env_release=200.0)

    # ── FFT-based analysis ────────────────────────────────────────────
    print("=" * 70)
    print("  FFT-BASED SINE ANALYSIS")
    print("=" * 70)
    analyze_sine(50.0, 0.5, 0.5, cfg_old)
    sweep_analysis(cfg_old)
    level_analysis(cfg_old)

    # ── Goertzel-based linearity ──────────────────────────────────────
    print("=" * 70)
    print("  GOERTZEL LEVEL LINEARITY")
    print("=" * 70)
    analyze_variant("ORIGINAL (harmonic_bass_enhancer.eel)", process_original, 50.0, cfg_old)
    analyze_variant("LINEAR (broken — divide after)", process_linear, 50.0, cfg_old)
    analyze_variant("FIXED (normalize before)", process_fixed, 50.0, cfg_old)

    # ── Detailed tone purity at single level ──────────────────────────
    print(f"\n{'='*70}")
    print(f"  DETAILED TONE PURITY: 50 Hz, A=0.5 (−6 dBFS)")
    print(f"{'='*70}")
    for label, fn in [("Original", process_original),
                       ("Linear (broken)", process_linear),
                       ("Fixed", process_fixed)]:
        pad = int(0.3 * cfg_old.fs)
        signal = generate_test_sine([50.0], [0.5], 0.5, cfg_old.fs, stereo=False)
        stereo = np.zeros((2, len(signal) + 2 * pad))
        stereo[0, pad:pad + len(signal)] = signal
        stereo[1, :] = stereo[0, :] * 0.95
        out = fn(stereo, cfg_old)
        ss_slice = slice(pad + len(signal) // 2, pad + len(signal))
        mono = (out[0, ss_slice] + out[1, ss_slice]) / 2.0
        m = measure_tones(mono, 50.0, cfg_old.fs, harmonics=(1, 2, 3))
        format_tone_report(m, label, fundamental=50.0)

    # ── Variant comparison ────────────────────────────────────────────
    for freq in [40, 50, 55]:
        for amp in [0.3, 0.5, 0.7]:
            compare_all_variants(freq=freq, amp=amp)

    # ── T₃ rejection ──────────────────────────────────────────────────
    compare_t3_rejection()

    # ── Full pipeline (EQ → enhancer → speaker) ──────────────────────
    pipeline_harmonic_control()

    # ── Below-cutoff sweep ────────────────────────────────────────────
    below_cutoff_sweep()

    # ── Loudness-weighted comparison ──────────────────────────────────
    loudness_weighted_comparison()

    # ── V2 linearity quick check ──────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  V2 LINEARITY: 50 Hz, varying amplitude")
    print(f"{'='*70}")
    print(f"  {'In dBFS':>8s}  {'H2 rel f':>10s}  {'H3 rel f':>10s}")
    print(f"  {'-'*32}")
    h2_vals, h3_vals = [], []
    for amp in [1.0, 0.5, 0.25, 0.125, 0.09]:
        pad = int(0.3 * cfg_v2.fs)
        sig = generate_test_sine([50], [amp], 0.5, cfg_v2.fs, stereo=False)
        stereo = np.zeros((2, len(sig) + 2*pad))
        stereo[0, pad:pad+len(sig)] = sig; stereo[1,:] = stereo[0,:]*0.95
        out = process_fixed_v2(stereo, cfg_v2)
        ss_slice = slice(pad+len(sig)//2, pad+len(sig))
        mono = (out[0, ss_slice] + out[1, ss_slice])/2
        f = goertzel_magnitude(mono, 50, cfg_v2.fs)
        h2 = goertzel_magnitude(mono, 100, cfg_v2.fs)
        h3 = goertzel_magnitude(mono, 150, cfg_v2.fs)
        h2r = 20*np.log10(max(h2,1e-12)/max(f,1e-12))
        h3r = 20*np.log10(max(h3,1e-12)/max(f,1e-12))
        h2_vals.append(h2r); h3_vals.append(h3r)
        print(f"  {20*np.log10(amp):+8.2f}  {h2r:+10.2f}  {h3r:+10.2f}")
    print(f"  H2 σ={np.std(h2_vals):.2f} dB  H3 σ={np.std(h3_vals):.2f} dB  "
          f"{'✅' if np.std(h2_vals)<0.01 and np.std(h3_vals)<0.01 else '❌'}")


if __name__ == "__main__":
    run()
