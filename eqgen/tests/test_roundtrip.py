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
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from eqgen import enhancer_ffi as effi
from eqgen.analysis import goertzel_magnitude
from eqgen.dsp import build_default_eq_coeffs
from eqgen.pipeline import run_pipeline, design_eq
from eqgen.model import compute_harmonic_efficacy


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: write a mono float64 signal as a 16-bit WAV
# ═══════════════════════════════════════════════════════════════════════════════

def _write_wav(path: str, samples: np.ndarray, rate: int = 44100):
    """Write mono float64 samples in [-1,1] as 16-bit PCM WAV."""
    samples = np.asarray(samples, dtype=np.float64)
    samples = np.clip(samples, -1.0, 1.0)
    int16 = (samples * 32767).astype(np.int16)

    import struct as _struct
    datasize = len(int16) * 2
    with open(path, 'wb') as f:
        f.write(b'RIFF')
        f.write(_struct.pack('<I', 36 + datasize))
        f.write(b'WAVEfmt ')
        f.write(_struct.pack('<I', 16))        # fmt chunk size
        f.write(_struct.pack('<H', 1))          # PCM
        f.write(_struct.pack('<H', 1))          # mono
        f.write(_struct.pack('<I', rate))
        f.write(_struct.pack('<I', rate * 2))   # byte rate
        f.write(_struct.pack('<H', 2))          # block align
        f.write(_struct.pack('<H', 16))         # bits per sample
        f.write(b'data')
        f.write(_struct.pack('<I', datasize))
        f.write(int16.tobytes())


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
# Helper: process sine tones through C enhancer with given EQ
# ═══════════════════════════════════════════════════════════════════════════════

def _measure_enhancer_output(
    freqs_hz: list,
    coeffs: list,
    fc: float, h2: float, h3: float,
    fs: float = 44100.0,
    pre_gain: float = 1.0,
    amplitude: float = 0.5,
    duration: float = 0.5,
) -> dict:
    """Feed sine tones through C enhancer, return {freq: {fund, h2, h3, peak}}."""
    results = {}
    eq_coeffs = list(coeffs) if coeffs else []

    for f in freqs_hz:
        n_samples = int(duration * fs)
        t = np.arange(n_samples) / fs
        sine = amplitude * np.sin(2.0 * np.pi * f * t)

        # Create enhancer
        enh = effi.create_enhancer(
            cutoff_hz=fc, h2_amp=float(h2), h3_amp=float(h3),
            release_secs=0.2, fs=fs, pre_gain=pre_gain,
            coeffs=eq_coeffs,
        )

        # Process stereo (mono duplicated)
        out = np.zeros(len(sine), dtype=np.float64)
        for i, x in enumerate(sine):
            lf, rf = effi.process_stereo_frame(enh, float(x), float(x))
            out[i] = float(lf)

        effi.destroy_enhancer(enh)

        # Measure (skip first 50ms for filter settling)
        settle = int(0.05 * fs)
        fund = float(goertzel_magnitude(out[settle:], f, fs))
        h2m = float(goertzel_magnitude(out[settle:], 2 * f, fs)) if 2 * f < fs / 2 else 0.0
        h3m = float(goertzel_magnitude(out[settle:], 3 * f, fs)) if 3 * f < fs / 2 else 0.0
        peak_out = float(np.max(np.abs(out[settle:])))

        results[f] = {"fund": fund, "h2": h2m, "h3": h3m, "peak": peak_out}

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Speaker models
# ═══════════════════════════════════════════════════════════════════════════════

def _small_speaker(f: float) -> float:
    """-12 dB at 50 Hz, flat above 100 Hz."""
    if f <= 50:
        return 0.25
    elif f >= 100:
        return 1.0
    else:
        return 0.25 + 0.75 * (f - 50) / 50


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
    meas_noise = _apply_speaker(_brown_noise(n_samples, seed=2), _small_speaker, fs)

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

    results = _measure_enhancer_output(
        test_freqs, coeffs, fc, h2, h3, fs=fs, pre_gain=pre_gain,
        amplitude=0.5, duration=0.5,
    )

    fundamentals = {}
    for f in test_freqs:
        r = results[f]
        fundamentals[f] = r["fund"]
        print(f"    {f:6.0f}  {r['fund']:8.4f}  {r['h2']:8.4f}  {r['h3']:8.4f}  {r['peak']:8.4f}")

    # 6. Verify correction: enhancer output is pre-emphasized for the
    #    speaker rolloff.  Multiply by speaker gain to get the "acoustic"
    #    output — this should be approximately flat across frequencies.
    corrected = np.array([fundamentals[f] * _small_speaker(f) for f in test_freqs])
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
        assert results[f]["peak"] <= 1.0, f"Peak {results[f]['peak']:.3f} > 1.0 at {f} Hz — clipping"

    print(f"\n  ✅ Round-trip passed — speaker corrected to ±{max_dev:.1f} dB flatness")


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
