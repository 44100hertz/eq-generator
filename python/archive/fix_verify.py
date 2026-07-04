"""
Fixed variant: normalize BEFORE Chebyshev, then scale back.

Key insight:
    Tₙ(cos θ) = cos(nθ)  is EXACT when the input has amplitude 1.

    For input x = A·cos(ωt), if we first normalize x̂ = x/A ≈ cos(ωt),
    then Tₙ(x̂) = cos(nωt) — NO leakage terms at all.

    Then multiply by A to restore level: A·cos(nωt).

Compare to the broken approach:
    Tₙ(x)/Aⁿ has leakage because Tₙ(A·cos θ) ≠ Aⁿ·cos(nθ) for A < 1.
    The 1/Aⁿ factor amplifies the leftover (A²−1)·DC and 3A(A²−1)·cosθ terms.

This is mathematically the same final result that level linearisation
was TRYING to achieve, but done correctly.
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from harmonic_bass import (
    BassEnhancerConfig,
    process_original,
    process_linear,
    design_butter_lp,
    design_butter_hp,
    biquad_tick,
    EnvFollower,
    generate_test_sine,
    magnitude_spectrum,
    find_peaks,
)


def process_fixed(
    samples: np.ndarray,
    cfg: BassEnhancerConfig,
) -> np.ndarray:
    """FIXED variant: normalize input BEFORE Chebyshev, scale back after.

    Core fix:
        _norm = _lp / max(env, floor)   # normalize to amplitude ≈ 1
        T₂(_norm) * env                  # clean 2nd harmonic, no DC offset
        T₃(_norm) * env                  # clean 3rd harmonic, no fund. leakage

    Only ONE envelope per channel needed (T₂ and T₃ share it).
    No correction limit needed (no divisions by tiny numbers that blow up).
    """
    assert samples.ndim == 2 and samples.shape[0] == 2, "Expected (2, N) stereo input"
    n = samples.shape[1]
    out = samples.copy()

    lp_coeffs = design_butter_lp(cfg.cutoff, cfg.fs)
    hp_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)
    hp_harm_coeffs = design_butter_hp(cfg.cutoff, cfg.fs)

    lp_L = np.zeros(4); lp_R = np.zeros(4)
    hp_L = np.zeros(4); hp_R = np.zeros(4)
    hp_harm_L = np.zeros(4); hp_harm_R = np.zeros(4)

    env_L = EnvFollower.from_params(cfg.env_release, cfg.fs)
    env_R = EnvFollower.from_params(cfg.env_release, cfg.fs)
    floor = 0.0001

    for i in range(n):
        sL = samples[0, i]
        sR = samples[1, i]

        l_lp, lp_L = biquad_tick(sL, lp_coeffs, lp_L)
        r_lp, lp_R = biquad_tick(sR, lp_coeffs, lp_R)

        env_L.tick(l_lp)
        env_R.tick(r_lp)

        envL_s = max(env_L.read(), floor)
        envR_s = max(env_R.read(), floor)

        # Normalize to amplitude ~1 before Chebyshev
        normL = l_lp / envL_s
        normR = r_lp / envR_s

        # T₂(norm) = 2·norm² − 1 = cos(2ωt)  [exact for |norm| ≤ 1]
        # T₃(norm) = 4·norm³ − 3·norm = cos(3ωt)  [exact for |norm| ≤ 1]
        # Scale back by envelope to restore amplitude
        harmL = envL_s * (
            cfg.h2_amp * (2.0 * normL * normL - 1.0)
            + cfg.h3_amp * (4.0 * normL * normL * normL - 3.0 * normL)
        )
        harmR = envR_s * (
            cfg.h2_amp * (2.0 * normR * normR - 1.0)
            + cfg.h3_amp * (4.0 * normR * normR * normR - 3.0 * normR)
        )

        wetL, hp_L = biquad_tick(sL, hp_coeffs, hp_L)
        harmL_hp, hp_harm_L = biquad_tick(harmL, hp_harm_coeffs, hp_harm_L)
        wetR, hp_R = biquad_tick(sR, hp_coeffs, hp_R)
        harmR_hp, hp_harm_R = biquad_tick(harmR, hp_harm_coeffs, hp_harm_R)

        out[0, i] = wetL + harmL_hp
        out[1, i] = wetR + harmR_hp

    return out


def compare_all(freq: float = 50.0, amp: float = 0.5, duration: float = 0.3):
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

    # Steady-state analysis
    start = pad + len(signal) // 2
    end = pad + len(signal)
    fft_size = min(16384, end - start)

    def monorms(out):
        mono = (out[0, start:end] + out[1, start:end]) / 2.0
        return mono, np.sqrt(np.mean(mono**2))

    results = {}
    for name, out in [("original", out_orig), ("linear", out_lin), ("fixed", out_fix)]:
        mono, rms = monorms(out)
        freqs, mags = magnitude_spectrum(mono, cfg.fs, fft_size)
        peaks = find_peaks(freqs, mags)
        results[name] = {"rms": rms, "peaks": peaks, "freqs": freqs, "mags": mags}

    print(f"\n{'='*70}")
    print(f"  COMPARISON: {freq} Hz sine, A={amp} ({20*np.log10(amp):+.1f} dBFS)")
    print(f"  h2_amp={cfg.h2_amp}, h3_amp={cfg.h3_amp}, cutoff={cfg.cutoff} Hz")
    print(f"{'='*70}")

    # Spectral peaks table
    print(f"\n  {'Freq (Hz)':>10s}  {'Original':>10s}  {'Linear':>10s}  {'Fixed':>10s}  {'Note':>30s}")
    print(f"  {'-'*72}")
    all_freqs = sorted(set(
        int(round(p[0])) for r in results.values() for p in r["peaks"]
    ))
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

    # RMS comparison
    print(f"\n  ── RMS levels ──")
    for name in ["original", "linear", "fixed"]:
        rms = results[name]["rms"]
        print(f"  {name:>10s}: {rms:.6f}  ({20*np.log10(rms):+.2f} dBFS)")
    print(f"  {'Delta lin':>10s}: {20*np.log10(results['linear']['rms']/results['original']['rms']):+.2f} dB")
    print(f"  {'Delta fix':>10s}: {20*np.log10(results['fixed']['rms']/results['original']['rms']):+.2f} dB")

    # Harmonic purity check for fixed variant
    print(f"\n  ── Fixed variant harmonic purity ──")
    f_peaks = results["fixed"]["peaks"]
    fund_peak = next((m for f, m in f_peaks if abs(f - freq) < 2), None)
    if fund_peak:
        for h, hf in [(2, 2 * freq), (3, 3 * freq)]:
            h_peak = next((m for f, m in f_peaks if abs(f - hf) < 2), None)
            if h_peak:
                db = 20 * np.log10(h_peak / fund_peak)
                print(f"  H{h} ({hf:.0f} Hz): {db:+.2f} dB relative to fundamental")
        # Check for unexpected peaks
        unexpected = [(f, m) for f, m in f_peaks
                      if abs(f - freq) > 3 and abs(f - 2 * freq) > 3 and abs(f - 3 * freq) > 3]
        if unexpected:
            print(f"  ⚠️  Unexpected peaks: {[(round(f,1), f'{20*np.log10(m/fund_peak):+.1f} dB') for f,m in unexpected]}")
        else:
            print(f"  ✅ No unexpected peaks — clean harmonic generation!")

    return results


def run():
    """Run the fix verification analysis."""
    print("=" * 70)
    print("  FIXED VARIANT: normalize before Chebyshev")
    print("=" * 70)

    for freq in [40, 50, 55]:
        for amp in [0.3, 0.5, 0.7]:
            compare_all(freq=freq, amp=amp)


if __name__ == "__main__":
    run()
