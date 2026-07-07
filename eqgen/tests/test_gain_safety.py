"""
Gain safety tests.

Analyze raw measurement data for excessive correction gain and unstable
EQ curves.  Uses only FFT analysis of measurement/target WAVs — no
psychoacoustic model.  The pipeline IS the response; we check
measurement quality directly.

Catches:
  - Excessive correction gain (speaker can't reproduce that frequency)
  - Sharp peak-to-trough variations (measurement noise amplified into resonance)
"""

import sys
import os
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))


def analyze_gain_safety(
    meas_path: str,
    target_path: str,
    max_gain_db: float = 24.0,
    max_variation_db_per_octave: float = 12.0,
) -> dict:
    """Analyze measurement vs target for unsafe correction gain levels."""
    import scipy.io.wavfile as wav

    rate, meas_raw = wav.read(meas_path)
    _, target_raw = wav.read(target_path)
    if meas_raw.ndim == 2:
        meas_raw = meas_raw.mean(axis=1)
    if target_raw.ndim == 2:
        target_raw = target_raw.mean(axis=1)
    meas = meas_raw.astype(float) / 32768.0
    target = target_raw.astype(float) / 32768.0

    print("=" * 70)
    print("  GAIN SAFETY ANALYSIS (raw measurement → correction)")
    print(f"  Limits: max gain = {max_gain_db} dB, "
          f"max variation = {max_variation_db_per_octave} dB/octave")
    print("=" * 70)

    # FFT both signals
    n = min(len(meas), len(target))
    window = np.hanning(n)
    meas_fft = np.abs(np.fft.rfft(meas[:n] * window)) / n
    target_fft = np.abs(np.fft.rfft(target[:n] * window)) / n
    freqs_fft = np.fft.rfftfreq(n, 1.0 / rate)

    def mag_at(freq):
        idx = np.argmin(np.abs(freqs_fft - freq))
        return meas_fft[idx], target_fft[idx], freqs_fft[idx]

    # Check measurement SNR at bass frequencies
    print(f"\n  ── Correction gain below 55 Hz ──")
    print(f"  {'Freq':>8s}  {'meas RMS':>12s}  {'target RMS':>12s}  "
          f"{'correction':>12s}  {'Note':>20s}")
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

        print(f"  {actual_f:8.1f}  {m:12.8f}  {t:12.8f}  "
              f"{comp_db:+11.1f} dB  {note:>20s}")

    # Compute correction curve = target / measurement (in dB)
    print(f"\n  ── Correction curve shape ──")

    mask = (freqs_fft >= 15) & (freqs_fft <= rate / 2 * 0.95)
    f_sub = freqs_fft[mask]
    m_sub = meas_fft[mask]
    t_sub = target_fft[mask]

    correction_db = 20.0 * np.log10(np.maximum(t_sub / np.maximum(m_sub, 1e-12), 1e-12))
    correction_db -= np.max(correction_db)  # normalize

    # Peak-to-trough analysis
    print(f"\n  ── Peak-to-trough analysis ──")
    peak_trough_ratios = []
    gains_sub = correction_db
    freqs_sub = f_sub

    for i in range(2, len(gains_sub) - 2):
        if gains_sub[i] > gains_sub[i-1] and gains_sub[i] > gains_sub[i+1]:
            left_min = np.min(gains_sub[max(0, i-10):i])
            right_min = np.min(gains_sub[i+1:min(len(gains_sub), i+11)])
            trough = max(left_min, right_min)
            if gains_sub[i] - trough > 6:
                octaves = np.log2(freqs_sub[i] / freqs_sub[max(0, i-10)])
                if octaves > 0:
                    variation = (gains_sub[i] - trough) / octaves
                    peak_trough_ratios.append((freqs_sub[i], gains_sub[i] - trough, variation))

    if peak_trough_ratios:
        print(f"  {'Freq':>8s}  {'Peak-Trough':>14s}  {'Var/octave':>14s}  {'Status':>12s}")
        print(f"  {'-'*52}")
        for f, pk, var in sorted(peak_trough_ratios, key=lambda x: -x[1])[:10]:
            status = "OK" if var <= max_variation_db_per_octave else "⚠️  UNSTABLE"
            print(f"  {f:8.1f}  {pk:+13.1f} dB  {var:+13.1f} dB  {status:>12s}")
            if var > max_variation_db_per_octave:
                issues.append(("unstable_curve", f, var))

    # Summary
    print(f"\n  ── Summary ──")
    if issues:
        print(f"  ❌ {len(issues)} safety issue(s) found:")
        for kind, freq, val in issues:
            if kind == "excessive_gain":
                print(f"     Excessive gain at {freq:.1f} Hz: {val:+.0f} dB")
            elif kind == "unstable_curve":
                print(f"     Unstable curve at {freq:.1f} Hz: {val:.1f} dB/oct variation")
    else:
        print(f"  ✅ All safety checks passed")

    return {
        "issues": issues,
        "max_gain": float(np.max(correction_db)),
        "max_variation": float(
            max([v for _, _, v in peak_trough_ratios]) if peak_trough_ratios else 0
        ),
    }


def run():
    """Run gain safety analysis on real measurement data."""
    root = Path(__file__).resolve().parent.parent.parent
    meas_path = root / "measurements" / "technics" / "standing" / "measurement2.wav"
    target_path = root / "measurements" / "technics" / "standing" / "target.wav"
    if meas_path.exists() and target_path.exists():
        analyze_gain_safety(str(meas_path), str(target_path))
    else:
        print("  ⚠️  Measurement files not found — skipping gain safety test.")


if __name__ == "__main__":
    run()
