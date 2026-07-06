"""
Gain safety tests.

Verify that the EQ preprocessing pipeline doesn't produce dangerously
large gains or unstable curves from real measurement data.

Catches:
  - Excessive correction gain (speaker can't reproduce that frequency)
  - Sharp peak-to-trough variations (measurement noise amplified into resonance)
  - fc/3 truncation creating artificial cliffs
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.model import (
    preprocess_eq_curve,
    butterworth_hp_mag,
    butterworth_lp_mag,
)


def analyze_gain_safety(
    meas_path: str,
    target_path: str,
    fc: float = 50.0,
    max_gain_db: float = 24.0,
    max_variation_db_per_octave: float = 12.0,
) -> dict:
    """Analyze an EQ curve for unsafe gain levels."""
    import scipy.io.wavfile as wav

    rate, meas_raw = wav.read(meas_path)
    _, target_raw = wav.read(target_path)
    if meas_raw.ndim == 2:
        meas_raw = meas_raw.mean(axis=1)
    if target_raw.ndim == 2:
        target_raw = target_raw.mean(axis=1)
    meas = meas_raw.astype(float) / 32768.0
    target = target_raw.astype(float) / 32768.0

    h2, h3 = 0.33, 0.33

    print("=" * 70)
    print("  GAIN SAFETY ANALYSIS")
    print(f"  fc={fc} Hz, h2={h2}, h3={h3}")
    print(f"  Limits: max gain = {max_gain_db} dB, max variation = {max_variation_db_per_octave} dB/octave")
    print("=" * 70)

    # ── FFT both signals once ────────────────────────────────────────
    n = min(len(meas), len(target))
    window = np.hanning(n)
    meas_fft = np.abs(np.fft.rfft(meas[:n] * window)) / n
    target_fft = np.abs(np.fft.rfft(target[:n] * window)) / n
    freqs_fft = np.fft.rfftfreq(n, 1.0 / rate)

    def mag_at(freq):
        """Get magnitude at a specific frequency via nearest FFT bin."""
        idx = np.argmin(np.abs(freqs_fft - freq))
        return meas_fft[idx], target_fft[idx], freqs_fft[idx]

    # ── 1. Check measurement SNR at bass frequencies ──────────────────
    print(f"\n  ── Measurement SNR below cutoff ──")
    print(f"  {'Freq':>8s}  {'meas RMS':>12s}  {'target RMS':>12s}  {'correction':>12s}  {'Note':>20s}")
    print(f"  {'-'*68}")

    issues = []
    for f in range(16, 55):
        m, t, actual_f = mag_at(f)
        comp = t / m if m > 1e-12 else 1e6
        comp_db = 20 * np.log10(max(comp, 1e-12))

        note = ""
        if comp_db > max_gain_db:
            note = f"⚠️  EXCESSIVE ({comp_db:+.0f} dB)"
            issues.append(("excessive_gain", actual_f, comp_db))
        elif m < 1e-7:
            note = "speaker near silence"

        print(f"  {actual_f:8.1f}  {m:12.8f}  {t:12.8f}  {comp_db:+11.1f} dB  {note:>20s}")

    # ── 2. Compute enhancer EQ curve at FFT bin resolution ───────────
    print(f"\n  ── Enhancer EQ curve shape ──")

    points = []
    for i, freq in enumerate(freqs_fft):
        if freq < 15 or freq > fc * 2.5:
            continue
        m_f  = meas_fft[i]
        t_f  = target_fft[i]
        m_2f = np.interp(2*freq, freqs_fft, meas_fft) if 2*freq < rate/2 else 0
        m_3f = np.interp(3*freq, freqs_fft, meas_fft) if 3*freq < rate/2 else 0

        comp = t_f / m_f if m_f > 1e-12 else 1e6

        hp_f  = butterworth_hp_mag(freq, fc)
        hp_2f = butterworth_hp_mag(2*freq, fc) if 2*freq < rate/2 else 1.0
        hp_3f = butterworth_hp_mag(3*freq, fc) if 3*freq < rate/2 else 1.0
        lp_f  = butterworth_lp_mag(freq, fc)
        lp_f2 = butterworth_lp_mag(freq, fc / 2)

        a = (hp_f  * m_f)**2
        b = h2**2 * lp_f**2  * hp_2f**2 * m_2f**2
        c = h3**2 * lp_f2**2 * hp_3f**2 * m_3f**2

        denom = np.sqrt(a + b + c)
        target_amp = comp * m_f
        G = target_amp / denom if denom > 1e-12 else 0.01
        G = min(G, comp)
        points.append((freq, 20 * np.log10(max(G, 1e-12))))

    freqs_arr = np.array([p[0] for p in points])
    gains_arr = np.array([p[1] for p in points])

    max_g = np.max(gains_arr)
    gains_norm = gains_arr - max_g

    # ── 3. Check peak-to-trough variation ─────────────────────────────
    print(f"\n  ── Peak-to-trough analysis (sub-{fc*2:.0f} Hz) ──")
    mask = freqs_arr < fc * 2
    sub_freqs = freqs_arr[mask]
    sub_gains = gains_norm[mask]

    # Find all local peaks and adjacent troughs
    peak_trough_ratios = []
    for i in range(2, len(sub_gains) - 2):
        if sub_gains[i] > sub_gains[i-1] and sub_gains[i] > sub_gains[i+1]:
            # This is a local peak. Find nearest troughs on each side
            left_min = np.min(sub_gains[max(0,i-10):i])
            right_min = np.min(sub_gains[i+1:min(len(sub_gains),i+11)])
            trough = max(left_min, right_min)
            if sub_gains[i] - trough > 6:  # significant peak
                octaves = np.log2(sub_freqs[i] / sub_freqs[max(0,i-10)])
                if octaves > 0:
                    variation = (sub_gains[i] - trough) / octaves
                    peak_trough_ratios.append((sub_freqs[i], sub_gains[i] - trough, variation))

    if peak_trough_ratios:
        print(f"  {'Freq':>8s}  {'Peak-Trough':>14s}  {'Var/octave':>14s}  {'Status':>12s}")
        print(f"  {'-'*52}")
        max_var = 0
        for f, pk, var in sorted(peak_trough_ratios, key=lambda x: -x[1])[:10]:
            status = "OK" if var <= max_variation_db_per_octave else "⚠️  UNSTABLE"
            max_var = max(max_var, var)
            print(f"  {f:8.1f}  {pk:+13.1f} dB  {var:+13.1f} dB  {status:>12s}")
            if var > max_variation_db_per_octave:
                issues.append(("unstable_curve", f, var))

    # ── 4. Check fc/3 boundary ────────────────────────────────────────
    third = fc / 3.0
    print(f"\n  ── fc/3 boundary (f={third:.1f} Hz) ──")
    near_third = [(f, g) for f, g in zip(freqs_arr, gains_norm) if abs(f - third) < 3]
    if near_third:
        third_gain = min(near_third, key=lambda x: abs(x[0] - third))
        nearby_above = [(f, g) for f, g in near_third if f > third + 1]
        if nearby_above:
            closest_above = min(nearby_above, key=lambda x: x[0])
            cliff = third_gain[1] - closest_above[1]
            print(f"  Gain at fc/3 ({third:.1f} Hz):  {third_gain[1]:+.1f} dB")
            print(f"  Gain above fc/3 ({closest_above[0]:.1f} Hz): {closest_above[1]:+.1f} dB")
            print(f"  Cliff: {cliff:+.1f} dB")
            if abs(cliff) > 10:
                print(f"  ⚠️  Sharp fc/3 boundary — sub-20 Hz will be loud")
                issues.append(("fc3_cliff", third, cliff))

    # ── 5. Summary ────────────────────────────────────────────────────
    print(f"\n  ── Summary ──")
    if issues:
        print(f"  ❌ {len(issues)} safety issue(s) found:")
        for kind, freq, val in issues:
            if kind == "excessive_gain":
                print(f"     Excessive gain at {freq:.1f} Hz: {val:+.0f} dB")
            elif kind == "unstable_curve":
                print(f"     Unstable curve at {freq:.1f} Hz: {val:.1f} dB/oct variation")
            elif kind == "fc3_cliff":
                print(f"     fc/3 cliff: {val:+.0f} dB at {freq:.1f} Hz boundary")
    else:
        print(f"  ✅ All safety checks passed")

    return {
        "issues": issues,
        "max_gain": float(np.max(gains_arr)),
        "max_variation": float(max([v for _, _, v in peak_trough_ratios]) if peak_trough_ratios else 0),
    }


def run():
    """Run gain safety analysis on the sanity check data."""
    # Path relative to project root
    import os
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    meas_path = os.path.join(root, "measurements", "technics", "standing", "measurement2.wav")
    target_path = os.path.join(root, "measurements", "technics", "standing", "target.wav")
    analyze_gain_safety(meas_path, target_path, fc=50.0)


if __name__ == "__main__":
    run()
