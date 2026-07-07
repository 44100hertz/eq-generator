"""
Pure math verification of the Chebyshev harmonic generation and level linearisation.

For a pure sine x = A·cos(ωt), the Chebyshev polynomials produce:

    T₂(x) = 2x² − 1 = A²·cos(2ωt) + (A² − 1)           [2nd harmonic + DC offset]
    T₃(x) = 4x³ − 3x = A³·cos(3ωt) + 3A(A²−1)·cos(ωt)  [3rd harmonic + fundamental leakage]

Level linearisation (original intent: cancel A²/A³ nonlinearity so output ∝ A):

    T₂ / A  = A·cos(2ωt) + (A − 1/A)                     [DC offset for A≠1]
    T₃ / A² = A·cos(3ωt) + 3(A − 1/A)·cos(ωt)            [FUNDAMENTAL LEAKAGE for A≠1]

The fundamental leakage term in T₃/A² is:
    3(A − 1/A)·cos(ωt)
    For A = 0.5:  3(0.5 − 2.0) = −4.5  →  inverted, 4.5× larger than desired
    For A = 0.1:  3(0.1 − 10.0) = −29.7 →  absolutely massive (before correction limit)

This leakage passes through the output HP filter. Even though the HP attenuates
frequencies near the fundamental, the leakage amplitude can be so large that the
attenuated residue is audible as an "overtone" or tonal coloration.

The correction limit (corr_limit dB) caps 1/A and 1/A², limiting the worst case.
But it also creates a discontinuity where the correction clips.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.dsp import design_butter_lp, design_butter_hp, biquad_tick
from eqgen.dsp import butterworth_hp_mag


def analyze_chebyshev_mix(
    freq: float = 50.0,
    amplitude: float = 0.5,
    h2_amp: float = 0.13,
    h3_amp: float = 0.10,
    cutoff: float = 60.0,
    fs: float = 44100.0,
    corr_limit_db: float = 20.0,
    duration: float = 0.2,
):
    """Analyze the harmonic generator math on a pure sine with LP pre-filtering.

    Simulates the exact signal chain:
      1. LP filter at `cutoff` (attenuates if freq ≈ cutoff)
      2. T₂ and T₃ Chebyshev polynomials
      3. Level linearisation: T₂/A + T₃/A² (with correction limit)
      4. HP filter at `cutoff` on harmonic output

    Reports all spectral components and their amplitudes.
    """
    # ── 1. Generate pure sine ──────────────────────────────────────────
    n = int(duration * fs)
    t = np.arange(n) / fs
    x_in = amplitude * np.sin(2.0 * np.pi * freq * t)

    # ── 2. LP filter (2nd order Butterworth, bilinear) ─────────────────
    lp_coeffs = design_butter_lp(cutoff, fs)
    hp_coeffs = design_butter_hp(cutoff, fs)

    lp_state = np.zeros(4)
    x_lp = np.zeros(n)
    for i in range(n):
        x_lp[i], lp_state = biquad_tick(x_in[i], lp_coeffs, lp_state)

    # Steady-state LP output: after filter settling
    settle = int(0.05 * fs)  # skip first 50ms
    x_lp_ss = x_lp[settle:]

    # Effective amplitude after LP
    lp_gain = np.sqrt(2.0 * np.mean(x_lp_ss**2)) / (amplitude / np.sqrt(2))
    A_eff = amplitude * lp_gain

    print(f"  LP filter gain at {freq} Hz (cutoff={cutoff}): {20*np.log10(lp_gain):+.2f} dB")
    print(f"  Effective amplitude into Chebyshev: A_eff = {A_eff:.4f}")

    # ── 3. Chebyshev polynomials ───────────────────────────────────────
    t2_raw = 2.0 * x_lp_ss**2 - 1.0        # T₂(x)
    t3_raw = 4.0 * x_lp_ss**3 - 3.0 * x_lp_ss  # T₃(x)

    # ── 4. Theoretical spectral components ─────────────────────────────
    print(f"\n  ── Theoretical Chebyshev output (A_eff = {A_eff:.4f}) ──")
    print(f"  {'Component':>20s}  {'Amplitude':>10s}  {'dB rel A':>10s}  {'Note':>30s}")

    # T₂: A²·cos(2ωt) + (A²−1)
    t2_h2_amp = A_eff**2
    t2_dc = A_eff**2 - 1.0
    print(f"  {'T₂ 2nd harmonic':>20s}  {t2_h2_amp:10.6f}  {20*np.log10(t2_h2_amp/(A_eff+1e-12)):+10.2f}  {'(at 2f) desired'}")
    print(f"  {'T₂ DC offset':>20s}  {abs(t2_dc):10.6f}  {'':>10s}  {'(removed by HP)'}")

    # T₃: A³·cos(3ωt) + 3A(A²−1)·cos(ωt)
    t3_h3_amp = A_eff**3
    t3_leak = 3.0 * A_eff * (A_eff**2 - 1.0)
    print(f"  {'T₃ 3rd harmonic':>20s}  {t3_h3_amp:10.6f}  {20*np.log10(t3_h3_amp/(A_eff+1e-12)):+10.2f}  {'(at 3f) desired'}")
    print(f"  {'T₃ fundamental leak':>20s}  {abs(t3_leak):10.6f}  {20*np.log10(abs(t3_leak)/(A_eff+1e-12)):+10.2f}  {'(at f) LEAKAGE!' if abs(t3_leak) > 1e-6 else '(at f) negligible'}")

    # ── 5. Level linearisation ─────────────────────────────────────────
    env = A_eff  # perfect envelope tracking in steady state
    max_corr = 10.0 ** (corr_limit_db / 20.0)

    corr_2 = min(1.0 / max(env, 1e-4), max_corr)
    corr_3 = min(1.0 / max(env * env, 1e-4), max_corr)

    print(f"\n  ── Level linearisation (env = A_eff = {env:.4f}) ──")
    print(f"  T₂ correction: 1/A = {1.0/env:.4f} (capped at {max_corr:.1f}) → effective: {corr_2:.4f}")
    print(f"  T₃ correction: 1/A² = {1.0/(env*env):.4f} (capped at {max_corr:.1f}) → effective: {corr_3:.4f}")

    # Corrected amplitudes
    t2_h2_corrected = h2_amp * corr_2 * t2_h2_amp
    t2_dc_corrected = h2_amp * corr_2 * t2_dc
    t3_h3_corrected = h3_amp * corr_3 * t3_h3_amp
    t3_leak_corrected = h3_amp * corr_3 * t3_leak

    print(f"\n  ── Corrected harmonic mix (h2={h2_amp}, h3={h3_amp}) ──")
    print(f"  {'Component':>20s}  {'Amplitude':>10s}  {'dB rel A':>10s}")
    print(f"  {'T₂ 2nd harmonic':>20s}  {t2_h2_corrected:10.6f}  {20*np.log10(abs(t2_h2_corrected)/(A_eff+1e-12)):+10.2f}")
    print(f"  {'T₃ 3rd harmonic':>20s}  {t3_h3_corrected:10.6f}  {20*np.log10(abs(t3_h3_corrected)/(A_eff+1e-12)):+10.2f}")
    print(f"  {'T₃ fund. leakage':>20s}  {abs(t3_leak_corrected):10.6f}  {20*np.log10(abs(t3_leak_corrected)/(A_eff+1e-12)):+10.2f}")
    print(f"  {'T₂ DC offset':>20s}  {abs(t2_dc_corrected):10.6f}  {'(removed by HP)':>10s}")

    # ── 6. HP filter effect ────────────────────────────────────────────
    hp_gain_f_db = 20 * np.log10(butterworth_hp_mag(freq, cutoff))
    hp_gain_2f_db = 20 * np.log10(butterworth_hp_mag(2 * freq, cutoff))
    hp_gain_3f_db = 20 * np.log10(butterworth_hp_mag(3 * freq, cutoff))

    hp_gain_f = 10 ** (hp_gain_f_db / 20)
    hp_gain_2f = 10 ** (hp_gain_2f_db / 20)
    hp_gain_3f = 10 ** (hp_gain_3f_db / 20)

    print(f"\n  ── Output HP filter (cutoff={cutoff} Hz) ──")
    print(f"  At {freq:.0f} Hz (fundamental): {hp_gain_f_db:+.1f} dB  →  gain = {hp_gain_f:.4f}")
    print(f"  At {2*freq:.0f} Hz (2nd harmonic):  {hp_gain_2f_db:+.1f} dB  →  gain = {hp_gain_2f:.4f}")
    print(f"  At {3*freq:.0f} Hz (3rd harmonic):  {hp_gain_3f_db:+.1f} dB  →  gain = {hp_gain_3f:.4f}")

    # After HP
    leak_after_hp = t3_leak_corrected * hp_gain_f
    h2_after_hp = t2_h2_corrected * hp_gain_2f
    h3_after_hp = t3_h3_corrected * hp_gain_3f

    print(f"\n  ── After HP filter ──")
    print(f"  T₃ fundamental leak after HP: {abs(leak_after_hp):.6f}  ({20*np.log10(abs(leak_after_hp)/(A_eff+1e-12)):+.2f} dB rel A)")
    print(f"  T₂ 2nd harmonic after HP:      {h2_after_hp:.6f}  ({20*np.log10(abs(h2_after_hp)/(A_eff+1e-12)):+.2f} dB rel A)")
    print(f"  T₃ 3rd harmonic after HP:      {h3_after_hp:.6f}  ({20*np.log10(abs(h3_after_hp)/(A_eff+1e-12)):+.2f} dB rel A)")

    # ── 7. Summary ─────────────────────────────────────────────────────
    print(f"\n  ════════════════════════════════════════════════════════════")
    print(f"  SUMMARY: Two signals at {freq:.0f} Hz in the final output")
    print(f"  ════════════════════════════════════════════════════════════")
    print(f"  1. Original (HP'd):                     {A_eff * hp_gain_f:.6f}")
    print(f"  2. T₃ fundamental leakage (HP'd):      {leak_after_hp:.6f}")
    print(f"     → Combined at {freq:.0f} Hz:         {(A_eff + t3_leak_corrected) * hp_gain_f:.6f}")
    print(f"     → Net change in fundamental level:   {20*np.log10(abs((A_eff + t3_leak_corrected) / A_eff)):+.2f} dB")

    if abs(t3_leak_corrected) > A_eff * 0.1:
        print(f"\n  ⚠️  FUNDAMENTAL LEAKAGE IS SIGNIFICANT!")
        if t3_leak_corrected < 0:
            print(f"  The leakage is INVERTED (negative), partially cancelling the original fundamental.")
        print(f"  The 2nd/3rd harmonics become RELATIVELY louder = perceived 'overtone'.")
        print(f"\n  Root cause: T₃(x) = 4x³−3x produces fundamental leakage for |x| < 1.")
        print(f"  Level linearisation (÷A²) amplifies this leakage by 1/A².")
        print(f"  Fix options:")
        print(f"    1. Apply a separate LP at fc/2 for T₃ (as in original plugin)")
        print(f"    2. Subtract the fundamental from the T₃ output before linearisation")
        print(f"    3. Use a different polynomial that is orthogonal for all amplitudes")
        print(f"    4. Apply the correction DIVISION to the FILTERED harmonic only (post-HP)")

    return {
        "A_eff": A_eff,
        "t2_h2_corrected": t2_h2_corrected,
        "t3_h3_corrected": t3_h3_corrected,
        "t3_leak_corrected": t3_leak_corrected,
        "leak_after_hp": leak_after_hp,
        "h2_after_hp": h2_after_hp,
        "h3_after_hp": h3_after_hp,
    }


def param_sweep():
    """Sweep amplitude and frequency to find worst-case leakage."""
    print("=" * 70)
    print("  PARAMETER SWEEP: fundamental leakage vs amplitude & frequency")
    print("=" * 70)

    freqs = [30, 40, 50, 55, 60, 70]
    amps = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    cutoff = 60.0
    corr_limit_db = 20.0
    max_corr = 10.0 ** (corr_limit_db / 20.0)

    print(f"\n  {'Freq':>6s}  {'Amp':>6s}  {'A_eff':>8s}  {'Corr₂':>8s}  {'Corr₃':>8s}  "
          f"{'Leak(dB)':>10s}  {'H2(dB)':>10s}  {'H3(dB)':>10s}")

    for f in freqs:
        lp_gain = 1.0 / np.sqrt(1.0 + (f / cutoff)**4)

        for amp in amps:
            A_eff = amp * lp_gain
            corr_2 = min(1.0 / max(A_eff, 1e-4), max_corr)
            corr_3 = min(1.0 / max(A_eff * A_eff, 1e-4), max_corr)

            h2 = 0.13 * corr_2 * A_eff**2
            h3 = 0.10 * corr_3 * A_eff**3
            leak = 0.10 * corr_3 * 3.0 * A_eff * (A_eff**2 - 1.0)

            hp_f = butterworth_hp_mag(f, cutoff)
            hp_2f = butterworth_hp_mag(2 * f, cutoff)
            hp_3f = butterworth_hp_mag(3 * f, cutoff)

            leak_db = 20 * np.log10(abs(leak * hp_f) / (A_eff * hp_f + 1e-12)) if abs(leak) > 1e-12 else -200
            h2_db = 20 * np.log10(abs(h2 * hp_2f) / (A_eff * hp_f + 1e-12))
            h3_db = 20 * np.log10(abs(h3 * hp_3f) / (A_eff * hp_f + 1e-12))

            print(f"  {f:6.0f}  {amp:6.2f}  {A_eff:8.4f}  {corr_2:8.2f}  {corr_3:8.2f}  "
                  f"{leak_db:+10.2f}  {h2_db:+10.2f}  {h3_db:+10.2f}")

    # ── Fix analysis ───────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print("  FIX ANALYSIS: what if we skip T₃ level linearisation?")
    print("=" * 70)

    print(f"\n  Without level linearisation on T₃:")
    print(f"  T₃ output = h3_amp * (4x³ - 3x) = h3_amp * [A³·cos(3ωt) + 3A(A²−1)·cos(ωt)]")
    print(f"  Fundamental leak = h3_amp * 3A(A²−1)  (NOT amplified by 1/A²)")
    print(f"")

    for amp in [0.1, 0.3, 0.5, 0.7, 1.0]:
        leak = 0.10 * 3.0 * amp * (amp**2 - 1.0)
        h3 = 0.10 * amp**3
        leak_db = 20 * np.log10(abs(leak) / (amp + 1e-12)) if abs(leak) > 1e-12 else -200
        h3_db = 20 * np.log10(h3 / (amp + 1e-12))
        print(f"  A={amp:.1f}:  leak={leak_db:+6.1f} dB  h3={h3_db:+6.1f} dB")


def run():
    """Run the chebyshev math analysis."""
    print("=" * 70)
    print("  CHEBYSHEV MATH ANALYSIS")
    print("=" * 70)
    analyze_chebyshev_mix()
    param_sweep()


if __name__ == "__main__":
    run()
