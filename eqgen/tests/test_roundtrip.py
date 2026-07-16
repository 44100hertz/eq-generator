"""
Round-trip integration test: synthetic speaker → pipeline → enhancer → verify.

Covers the full signal chain:
  1. Generate brown-noise measurement through a known speaker rolloff
  2. Run the EQ pipeline (Welch FFT → CV smoothing → correction curve)
  3. Fit IIR biquads via greedy peaking fitter
  4. Compute h2/h3 harmonic efficacy
  5. Feed sine tones through the C enhancer with fitted EQ
  6. Verify output fundamental amplitudes are consistent (±3 dB across bass)

This is the only test that exercises measurement → EQ design → enhancer
in one shot, catching regressions that component tests miss.
"""

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.io import wavfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen.pipeline import run_pipeline, design_eq
from eqgen.model import compute_harmonic_efficacy, small_speaker
from eqgen.sweep import run_sine_sweep


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: write a mono float64 signal as a 16-bit WAV
# ═══════════════════════════════════════════════════════════════════════════════

def _write_wav(path: str, samples: np.ndarray, rate: int = 44100):
    """Write mono float64 samples in [-1,1] as 16-bit PCM WAV."""
    samples = np.clip(np.asarray(samples, dtype=np.float64), -1.0, 1.0)
    int16 = (samples * 32767).astype(np.int16)
    wavfile.write(path, rate, int16)


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: generate brown (1/f²) noise
# ═══════════════════════════════════════════════════════════════════════════════

def _brown_noise(n_samples: int, seed: int = 42) -> np.ndarray:
    """Brown noise: cumulative sum of white noise → 1/f² spectrum."""
    rng = np.random.RandomState(seed)
    white = rng.randn(n_samples).astype(np.float64)
    brown = np.cumsum(white)
    brown /= np.max(np.abs(brown)) + 1e-10
    return brown * 0.9


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: apply a frequency-dependent gain (speaker model) via FFT filtering
# ═══════════════════════════════════════════════════════════════════════════════

def _apply_speaker(signal: np.ndarray, speaker_fn, rate: float) -> np.ndarray:
    """Apply speaker frequency response to a signal via FFT overlap-add."""
    n = len(signal)
    freqs = np.fft.rfftfreq(n, 1.0 / rate)
    H = np.array([speaker_fn(f) for f in freqs], dtype=np.float64)
    spectrum = np.fft.rfft(signal)
    spectrum *= H
    result = np.fft.irfft(spectrum, n=n)
    peak = np.max(np.abs(result))
    if peak > 0.95:
        result /= peak / 0.95
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Speaker models
# ═══════════════════════════════════════════════════════════════════════════════

def _steep_speaker(f: float) -> float:
    """-20 dB at 40 Hz, flat above 150 Hz."""
    if f <= 40:
        return 0.1
    elif f >= 150:
        return 1.0
    else:
        return 0.1 + 0.9 * (f - 40) / 110


# ═══════════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════════

def test_roundtrip_small_speaker():
    """Full pipeline → enhancer round-trip for a -12 dB @ 50 Hz speaker."""
    print("=" * 70)
    print("  ROUND-TRIP: Small speaker (-12 dB @ 50 Hz)")
    print("=" * 70)

    fs = 44100.0
    duration = 5.0       # seconds of measurement noise
    n_samples = int(fs * duration)
    fc = 60.0            # bass enhancer cutoff

    # 1. Generate brown-noise "target" (flat reference) and "measurement"
    #    (through the speaker rolloff)
    target_noise = _brown_noise(n_samples, seed=1)
    meas_noise = _apply_speaker(_brown_noise(n_samples, seed=2), small_speaker, fs)

    with tempfile.TemporaryDirectory() as tmpdir:
        target_path = os.path.join(tmpdir, "target.wav")
        meas_path = os.path.join(tmpdir, "measurement.wav")
        _write_wav(target_path, target_noise, int(fs))
        _write_wav(meas_path, meas_noise, int(fs))

        # 2. Run the pipeline
        freqs, gains_db, rate, max_gain_db, efficacy = run_pipeline(
            measurement_paths=[meas_path],
            target_path=target_path,
            bass_enhancer_cutoff=fc,
            n_eval=256,
        )

    h2 = efficacy["h2_amp"]
    h3 = efficacy["h3_amp"]
    print(f"\n  Pipeline output:")
    print(f"    Max gain:      {max_gain_db:+.1f} dB")
    print(f"    h2 efficacy:   {h2:.4f}")
    print(f"    h3 efficacy:   {h3:.4f}")
    print(f"    Correction at 50 Hz: {gains_db[np.argmin(np.abs(freqs - 50))]:+.1f} dB")

    # 3. Fit IIR biquads
    coeffs, bands, fit_freqs, fit_target, fit_response = design_eq(
        freqs, gains_db, fs, max_bands=12,
    )

    n_bands = len(bands)
    print(f"    IIR bands:      {n_bands}")
    for b in bands[:3]:
        print(f"      {b['type']:>10s} f0={b['f0']:6.1f} Hz  gain={b['gain_db']:+.1f} dB  Q={b['Q']:.1f}")

    # Verify the fit is reasonable
    fit_error = np.max(np.abs(fit_response - fit_target))
    print(f"    IIR fit max error: {fit_error:.2f} dB")
    assert fit_error < 3.0, f"IIR fit error {fit_error:.1f} dB exceeds 3 dB"

    # 4. Compute pre-gain
    from eqgen.dsp import pre_gain_from_max_gain
    pre_gain = float(pre_gain_from_max_gain(max_gain_db))
    print(f"    Pre-gain:       {pre_gain:.4f} ({20*np.log10(pre_gain):+.1f} dB)")

    # 5. Run enhancer on sine tones across the bass band
    test_freqs = [30, 40, 50, 60, 80, 100, 150, 200, 300, 500]
    print(f"\n  Enhancer output fundamentals:")
    print(f"    {'Freq':>6s}  {'Fund':>8s}  {'H2':>8s}  {'H3':>8s}  {'Peak':>8s}")
    print(f"    {'-'*46}")

    results = run_sine_sweep(
        freqs_hz=test_freqs, eq_coeffs=coeffs, fc=fc, h2=h2, h3=h3,
        fs=fs, pre_gain=pre_gain, amplitude=0.5, duration_sec=0.5,
    )

    fundamentals = {}
    for f in test_freqs:
        r = results[f]
        fundamentals[f] = r["fundamental"]
        print(f"    {f:6.0f}  {r['fundamental']:8.4f}  {r['h2']:8.4f}  {r['h3']:8.4f}  {r['rms']:8.4f}")

    # 6. Verify correction: enhancer output is pre-emphasized for the
    #    speaker rolloff.  Multiply by speaker gain to get the "acoustic"
    #    output — this should be approximately flat across frequencies.
    corrected = np.array([fundamentals[f] * small_speaker(f) for f in test_freqs])
    mid_corrected = np.median(corrected)
    corrected_db = 20.0 * np.log10(corrected / mid_corrected)
    max_dev = float(np.max(np.abs(corrected_db)))
    print(f"\n  Acoustic flatness (enhancer output × speaker): max dev = {max_dev:+.1f} dB")
    for f, db in zip(test_freqs, corrected_db):
        flag = "✅" if abs(db) < 3.0 else "⚠️"
        print(f"    {flag} {f:4.0f} Hz: {db:+.1f} dB")
    assert max_dev < 3.0, \
        f"Corrected deviation {max_dev:.1f} dB exceeds 3 dB — pipeline not correcting speaker"

    # 7. Verify the enhancer didn't clip
    for f in test_freqs:
        assert results[f]["rms"] <= 0.95, f"RMS {results[f]['rms']:.4f} > 0.95 at {f} Hz — likely clipping"

    print(f"\n  ✅ Round-trip passed — speaker corrected to ±{max_dev:.1f} dB flatness")

    # ── Harmonic stress: find maximum safe overboost ──
    # Fix push_gain low enough that the crossfade engages harmonics
    # for all bass signals, then sweep input amplitude (simulating
    # overboost_db) to find where the output peaks at 0 dBFS.
    # The amplitude at the clipping boundary, divided by 0.5 (the
    # nominal level used in the flatness test), gives the overboost
    # headroom in linear ratio → 20*log10 = overboost_db ceiling.
    print(f"\n  {'='*60}")
    print(f"  HARMONIC STRESS: finding max safe overboost")
    print(f"  {'='*60}")

    stress_freq = 50.0
    stress_freqs = [stress_freq]
    # Low push_gain forces crossfade engagement — harmonics always active.
    stress_push_gain = 0.3
    nominal_amp = 0.5  # reference amplitude from the flatness test

    # Binary search for amplitude where peak crosses 1.0
    lo_amp, hi_amp = 0.125, 8.0
    safe_amp = 0.125
    clip_amp = 8.0
    n_iters = 10
    for it in range(n_iters):
        mid_amp = np.sqrt(lo_amp * hi_amp)  # geometric: even dB steps
        r = run_sine_sweep(
            freqs_hz=stress_freqs, eq_coeffs=coeffs, fc=fc, h2=h2, h3=h3,
            fs=fs, pre_gain=pre_gain, push_gain=stress_push_gain,
            amplitude=float(mid_amp), duration_sec=0.5,
        )
        rms = r[stress_freq]["rms"]
        fund = r[stress_freq]["fundamental"]
        h2m = r[stress_freq]["h2"]
        h3m = r[stress_freq]["h3"]
        clips = rms >= 0.7
        marker = "CLIP" if clips else "safe"
        overboost_eq = 20.0 * np.log10(mid_amp / nominal_amp)
        print(f"    amp={mid_amp:.3f} (+{overboost_eq:+.1f} dB)  "
              f"rms={rms:.4f}  fund={fund:.4f}  h2={h2m:.4f}  h3={h3m:.4f}  [{marker}]")
        if clips:
            clip_amp = mid_amp
            hi_amp = mid_amp
        else:
            safe_amp = mid_amp
            lo_amp = mid_amp

    overboost_safe_db = 20.0 * np.log10(safe_amp / nominal_amp)
    overboost_clip_db = 20.0 * np.log10(clip_amp / nominal_amp)
    print(f"\n  Safe amplitude:     {safe_amp:.3f} (overboost ≤ {overboost_safe_db:+.1f} dB)")
    print(f"  Clipping amplitude: {clip_amp:.3f} (overboost ≥ {overboost_clip_db:+.1f} dB)")

    # Full profile at the safe boundary
    r_safe = run_sine_sweep(
        freqs_hz=stress_freqs, eq_coeffs=coeffs, fc=fc, h2=h2, h3=h3,
        fs=fs, pre_gain=pre_gain, push_gain=stress_push_gain,
        amplitude=safe_amp, duration_sec=0.5,
    )
    sr = r_safe[stress_freq]
    if sr["fundamental"] > 1e-10:
        h2_db = 20.0 * np.log10(max(sr["h2"], 1e-10) / sr["fundamental"])
        h3_db = 20.0 * np.log10(max(sr["h3"], 1e-10) / sr["fundamental"])
        print(f"  At boundary:  h2 rel fund = {h2_db:+.1f} dB  "
              f"h3 rel fund = {h3_db:+.1f} dB")

    assert safe_amp > nominal_amp, \
        f"Safe amplitude {safe_amp:.3f} ≤ nominal {nominal_amp} — no overboost headroom"
    print(f"\n  ✅ Harmonic stress passed — overboost ceiling ≈ {overboost_safe_db:+.1f} dB")
    print(f"     (use overboost_db ≤ {overboost_safe_db:+.1f} to avoid clipping)")


def test_roundtrip_flat_speaker():
    """Round-trip for a flat speaker: pipeline should produce near-zero correction."""
    print("\n" + "=" * 70)
    print("  ROUND-TRIP: Flat speaker (no correction needed)")
    print("=" * 70)

    fs = 44100.0
    duration = 4.0
    n_samples = int(fs * duration)
    fc = 60.0

    # Both target and measurement are flat (same brown noise with different seeds)
    target_noise = _brown_noise(n_samples, seed=10)
    meas_noise = _brown_noise(n_samples, seed=20)

    with tempfile.TemporaryDirectory() as tmpdir:
        target_path = os.path.join(tmpdir, "target.wav")
        meas_path = os.path.join(tmpdir, "measurement.wav")
        _write_wav(target_path, target_noise, int(fs))
        _write_wav(meas_path, meas_noise, int(fs))

        freqs, gains_db, rate, max_gain_db, efficacy = run_pipeline(
            measurement_paths=[meas_path],
            target_path=target_path,
            bass_enhancer_cutoff=fc,
            n_eval=256,
        )

    print(f"\n  Pipeline output:")
    print(f"    Max gain:      {max_gain_db:+.1f} dB")

    # For a flat speaker, correction should be small
    mid = (freqs >= 500) & (freqs <= 2000)
    mid_gains = gains_db[mid]
    max_dev = float(np.max(np.abs(mid_gains)))
    print(f"    Midrange deviation: ±{max_dev:.1f} dB")

    # With only 4s of brown noise, Welch variance is significant.
    # Allow up to 6 dB ripple — this test mainly verifies the pipeline
    # doesn't explode on a flat speaker.
    assert max_dev < 6.0, \
        f"Midrange deviation {max_dev:.1f} dB exceeded — pipeline may be over-correcting"

    # Fit biquads and verify fit error is low (biquads should be near-unity)
    coeffs, bands, fit_freqs, fit_target, fit_response = design_eq(
        freqs, gains_db, fs, max_bands=4,
    )
    fit_error = np.max(np.abs(fit_response - fit_target))
    print(f"    IIR bands:      {len(bands)}")
    print(f"    IIR fit error:  {fit_error:.2f} dB")
    assert fit_error < 2.0, f"Flat speaker IIR fit error {fit_error:.1f} dB is excessive"

    print(f"\n  ✅ Flat speaker round-trip passed")


# ═══════════════════════════════════════════════════════════════════════════════
# Runner
# ═══════════════════════════════════════════════════════════════════════════════

def run():
    """Run all round-trip tests."""
    test_roundtrip_small_speaker()
    test_roundtrip_flat_speaker()

    print("\n" + "=" * 70)
    print("  All round-trip tests passed.")
    print("=" * 70)


if __name__ == "__main__":
    run()
