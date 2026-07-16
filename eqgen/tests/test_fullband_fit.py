"""Full-band IIR fit test: 40 peaking filters on synthetic data."""
import numpy as np

from eqgen.eq_fit import fit_eq_curve, cascade_response_db
from eqgen.model import small_speaker


def run():
    fs = 44100.0
    freqs = np.logspace(np.log10(20), np.log10(20000), 300)

    # Target: inverse of small speaker (flat perceived output)
    speaker_linear = np.array([small_speaker(f) for f in freqs])
    target_db = 20.0 * np.log10(1.0 / np.maximum(speaker_linear, 1e-12))

    print("=" * 70)
    print("  40 PEAKING FILTERS — FULL-BAND FIT (20 Hz – 20 kHz)")
    print("  Target: flat output → inverse of speaker (-12 dB @ 50 Hz)")
    print("=" * 70)

    fit = fit_eq_curve(
        freqs, target_db, fs, max_bands=40,
        min_freq=20.0, max_freq=20000.0, min_peaking_freq=40.0,
    )

    fitted_db = cascade_response_db(fit.biquads, freqs, fs)
    err = fitted_db - target_db

    bands = [
        ("Sub-bass  (20-40 Hz)",    20, 40),
        ("Bass      (40-100 Hz)",   40, 100),
        ("Upper bass(100-250 Hz)",  100, 250),
        ("Low mid   (250-500 Hz)",  250, 500),
        ("Mid       (500-2k Hz)",   500, 2000),
        ("Upper mid (2k-8k Hz)",    2000, 8000),
        ("Treble    (8k-20k Hz)",   8000, 20000),
    ]

    print(f"\n{'Octave band':>25s}  {'Max err':>8s}  {'RMS err':>8s}")
    print(f"{'':->25s}  {'':->8s}  {'':->8s}")
    for label, lo, hi in bands:
        m = (freqs >= lo) & (freqs <= hi)
        if m.any():
            max_e = np.max(np.abs(err[m]))
            rms_e = np.sqrt(np.mean(err[m] ** 2))
            print(f"  {label:>23s}:  {max_e:+7.2f} dB  {rms_e:+7.2f} dB")

    overall_max = np.max(np.abs(err))
    overall_rms = np.sqrt(np.mean(err ** 2))
    print(f"\n  {'OVERALL':>23s}:  {overall_max:+7.2f} dB  {overall_rms:+7.2f} dB")

    print(f"\n  Band distribution ({len(fit.bands)} total):")
    for i, b in enumerate(fit.bands):
        print(f"    [{i:2d}] f0={b['f0']:6.1f} Hz  gain={b['gain_db']:+5.1f} dB  Q={b['Q']:.2f}")

    assert overall_max < 20.0, f"Full-band fit error {overall_max:.1f} dB exceeds 20 dB"
    print(f"\n  ✅ test_fullband_fit passed")
