"""
Fit a 2nd-order low shelf for loudness compensation.

KEY INSIGHT: The loudness shelf compensates for the *change* in ear
sensitivity as volume drops from peak (80 phon → quieter levels).
It does NOT provide the full absolute ISO 226 compensation (which would
require 50+ dB of bass boost even at peak volume — physically impossible).

Differential compensation from 80 phon → phon P:
    delta_comp(f, P) = comp(f, P) - comp(f, 80)
    where comp(f, P) = SPL(f, P) - SPL(1kHz, P)

This gives manageable gains (0-15 dB) and a clean linear relationship
with SPL drop at 1kHz.
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from eqgen.iso226 import contour_spl_db
from eqgen.eq_fit import design_low_shelf, biquad_response


def shelf_mag_db(freqs, fc, gain_db, Q, fs=48000):
    """Magnitude response in dB of a 2nd-order low shelf."""
    from eqgen.dsp import second_order_low_shelf_mag
    mag = np.array([second_order_low_shelf_mag(f, fc, gain_db, Q, fs) for f in freqs])
    return 20.0 * np.log10(np.maximum(mag, 1e-12))


def one_pole_shelf_mag_db(freqs, fc, gain_db):
    """Magnitude response in dB of the current one-pole shelf.
    H(s) = (s + G·ωc) / (s + ωc). DC gain = G, HF gain = 1.
    |H(f)|² = (G² + w²) / (1 + w²)  where w = f/fc."""
    w2 = (freqs / fc) ** 2
    G_lin = 10.0 ** (gain_db / 20.0)
    mag2 = (G_lin**2 + w2) / (1.0 + w2)
    return 10.0 * np.log10(np.maximum(mag2, 1e-12))


def abs_compensation(freqs, phon):
    """Absolute dB boost needed vs 1kHz at given phon level."""
    spl_1k = contour_spl_db(1000.0, phon)
    return np.array([contour_spl_db(f, phon) - spl_1k for f in freqs])


def differential_compensation(freqs, phon, ref_phon=80):
    """Additional compensation needed vs reference phon level.
    This is the change in ear sensitivity as volume drops from ref_phon.
    """
    comp_ref = abs_compensation(freqs, ref_phon)
    comp_phon = abs_compensation(freqs, phon)
    return comp_phon - comp_ref


def main():
    fs = 48000.0
    fit_freqs = np.logspace(np.log10(20), np.log10(1000), 200)
    display_freqs = np.array([20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200,
                              250, 315, 400, 500, 630, 800, 1000])
    mask = fit_freqs <= 1000.0
    REF_PHON = 80

    print("=" * 72)
    print("  LOUDNESS SHELF: DIFFERENTIAL ISO 226 FIT")
    print(f"  Reference: {REF_PHON} phon (peak volume, 0 shelf)")
    print("=" * 72)

    # ── Show the differential compensation curves ──
    print(f"\n  ── Differential compensation vs {REF_PHON} phon (dB) ──")
    header = f"  {'Freq':>7s}"
    for phon in [0, 20, 40, 60, 70]:
        header += f"  {'Δ'+str(phon):>7s}"
    print(header)
    print(f"  {'-'*47}")
    for f in display_freqs:
        if f > 1000:
            break
        line = f"  {f:5.0f} Hz"
        for phon in [0, 20, 40, 60, 70]:
            dc = differential_compensation(np.array([f]), phon, REF_PHON)[0]
            line += f"  {dc:+6.1f}"
        print(line)

    # ── Step 1: Fit fc, Q to minimize shape error across phon levels ──
    # Compute the average differential shape (normalized to unit gain).
    # The shape of the differential compensation is nearly identical across
    # phon levels — only the magnitude changes.
    fit_phons = list(range(0, 75, 10))
    diff_shapes = []
    for phon in fit_phons:
        diff = differential_compensation(fit_freqs, phon, REF_PHON)
        mx = float(np.max(diff))
        if mx > 0.01:
            diff_shapes.append(diff / mx)

    avg_diff_shape = np.mean(diff_shapes, axis=0)

    # Fit fc, Q where a shelf with unit gain best matches avg_diff_shape
    best_fc, best_Q = 60.0, 0.5
    best_err = float('inf')

    for fc in np.logspace(np.log10(25), np.log10(400), 70):
        for Q in np.logspace(np.log10(0.25), np.log10(4.0), 40):
            # Use 6 dB shelf (in the practical range) and normalize
            shelf_db = shelf_mag_db(fit_freqs, fc, 6.0, Q, fs)
            shelf_norm = shelf_db / float(np.max(shelf_db))

            err = np.sqrt(np.mean((shelf_norm[mask] - avg_diff_shape[mask])**2))
            if err < best_err:
                best_err = err
                best_fc, best_Q = fc, Q

    print(f"\n  ── Best (fc, Q) for differential shape ──")
    print(f"  fc = {best_fc:.1f} Hz   Q = {best_Q:.3f}   shape RMS err = {best_err:.4f}")

    # ── Step 2: With fixed fc/Q, fit gain at each phon level ──
    print(f"\n  ── Gain fit: fc={best_fc:.1f} Hz, Q={best_Q:.3f} ──")
    print(f"  {'Phon':>5s}  {'1kHz SPL':>9s}  {'Drop dB':>8s}  "
          f"{'Gain dB':>8s}  {'RMS err':>8s}")

    gains = []
    drops = []
    spl_1k_ref = contour_spl_db(1000.0, REF_PHON)

    for phon in range(0, 85, 10):
        if phon > 80:
            continue
        diff_target = differential_compensation(fit_freqs, phon, REF_PHON)
        max_t = float(np.max(diff_target[mask]))
        spl_1k = contour_spl_db(1000.0, phon)
        drop = spl_1k_ref - spl_1k

        if max_t <= 0.1:
            gain_db = 0.0
            err = 0.0
        else:
            # Ternary search for gain
            lo, hi = 0.0, max_t * 1.5
            for _ in range(30):
                m1 = lo + (hi - lo) / 3.0
                m2 = hi - (hi - lo) / 3.0
                e1 = np.sqrt(np.mean(
                    (shelf_mag_db(fit_freqs[mask], best_fc, m1, best_Q, fs)
                     - diff_target[mask])**2
                ))
                e2 = np.sqrt(np.mean(
                    (shelf_mag_db(fit_freqs[mask], best_fc, m2, best_Q, fs)
                     - diff_target[mask])**2
                ))
                if e1 < e2:
                    hi = m2
                else:
                    lo = m1
            gain_db = (lo + hi) / 2.0
            err = np.sqrt(np.mean(
                (shelf_mag_db(fit_freqs[mask], best_fc, gain_db, best_Q, fs)
                 - diff_target[mask])**2
            ))

        gains.append(gain_db)
        drops.append(drop)
        print(f"  {phon:5.0f}  {spl_1k:9.1f}  {drop:8.1f}  "
              f"{gain_db:8.1f}  {err:8.2f}")

    # ── Step 3: Linear fit gain = slope * drop ──
    drops_arr = np.array(drops)
    gains_arr = np.array(gains)

    # Through-origin fit
    slope = np.sum(gains_arr * drops_arr) / np.sum(drops_arr**2)

    # Also fit with intercept
    A = np.vstack([drops_arr, np.ones_like(drops_arr)]).T
    slope_w, intercept = np.linalg.lstsq(A, gains_arr, rcond=None)[0]

    print(f"\n  ── Linear relationships ──")
    print(f"  Through-origin:  shelf_db = {slope:.4f} × drop_db")
    print(f"  With intercept:  shelf_db = {slope_w:.4f} × drop_db + {intercept:.1f}")

    print(f"\n  {'Drop':>8s}  {'Gain':>8s}  {'Slope*drop':>11s}  {'Slope*drop+I':>13s}")
    for d, g in zip(drops, gains):
        print(f"  {d:8.1f}  {g:8.1f}  {slope*d:11.1f}  {slope_w*d+intercept:13.1f}")

    # ── Step 4: Comparison at key phon levels ──
    for phon_ref in [70, 60, 50, 40]:
        drop = contour_spl_db(1000.0, REF_PHON) - contour_spl_db(1000.0, phon_ref)
        diff_target = differential_compensation(display_freqs, phon_ref, REF_PHON)
        gain_ref = slope * drop  # using through-origin slope

        new_db = shelf_mag_db(display_freqs, best_fc, gain_ref, best_Q, fs)
        old_db = one_pole_shelf_mag_db(display_freqs, 200.0, gain_ref)

        print(f"\n  ── {phon_ref} phon (drop={drop:.1f} dB, shelf_gain={gain_ref:.1f} dB) ──")
        print(f"  {'Freq':>7s}  {'Target':>7s}  {'Old 1-pole':>10s}  "
              f"{'New 2nd':>10s}  {'Old err':>8s}  {'New err':>8s}")
        print(f"  {'-'*61}")

        old_rms, new_rms = 0.0, 0.0
        old_mid, new_mid = 0.0, 0.0
        n, nm = 0, 0
        for f, t, o, n2 in zip(display_freqs, diff_target, old_db, new_db):
            if f > 1000:
                break
            oe = o - t
            ne = n2 - t
            old_rms += oe**2
            new_rms += ne**2
            n += 1
            if 200 <= f <= 500:
                old_mid += oe**2
                new_mid += ne**2
                nm += 1
            print(f"  {f:6.0f} Hz  {t:+6.1f} dB  {o:+10.1f} dB  "
                  f"{n2:+10.1f} dB  {oe:+8.1f}  {ne:+8.1f}")

        old_rms = np.sqrt(old_rms / n)
        new_rms = np.sqrt(new_rms / n)
        old_mid = np.sqrt(old_mid / nm) if nm > 0 else 0
        new_mid = np.sqrt(new_mid / nm) if nm > 0 else 0
        print(f"  {'-'*61}")
        print(f"  RMS 20-1000Hz:  Old={old_rms:.1f}  New={new_rms:.1f}  "
              f"({(old_rms-new_rms)/old_rms*100:.0f}% better)")
        print(f"  RMS 200-500Hz:  Old={old_mid:.1f}  New={new_mid:.1f}  "
              f"({(old_mid-new_mid)/old_mid*100:.0f}% better)")

    # ── Summary ──
    print(f"\n  ═══════════════════════════════════════════════════════════")
    print(f"  RECOMMENDED PARAMETERS")
    print(f"  ═══════════════════════════════════════════════════════════")
    print(f"  fc          = {best_fc:.0f} Hz")
    print(f"  Q           = {best_Q:.3f}")
    print(f"  Slope       = {slope:.4f} dB shelf per dB SPL drop at 1kHz")
    print(f"               = {slope*10:.2f} dB shelf per 10 dB SPL drop")
    print(f"")
    print(f"  Old: fc=200 Hz, Q=N/A (one-pole), slope=0.08")
    print(f"  New: fc={best_fc:.0f} Hz, Q={best_Q:.3f}, slope={slope:.4f}")
    print(f"")
    print(f"  Given speaker_level_db and vol (0-127):")
    print(f"    drop_db  = speaker_level_db * (1 - vol/127)")
    print(f"    shelf_db = {slope:.4f} * drop_db")
    print(f"  ═══════════════════════════════════════════════════════════")

    # ── Shelf shape comparison at fixed gain ──
    print(f"\n  ── Shape overlay (6 dB gain, normalized) ──")
    freqs_viz = np.logspace(np.log10(20), np.log10(1000), 100)
    new_norm = shelf_mag_db(freqs_viz, best_fc, 6.0, best_Q, fs)
    new_norm = new_norm / float(np.max(new_norm))
    old_norm = one_pole_shelf_mag_db(freqs_viz, 200.0, 6.0)
    old_norm = old_norm / float(np.max(old_norm))
    diff_norm = avg_diff_shape  # from fit_freqs

    print(f"  {'Freq':>7s}  {'Target':>7s}  {'Old':>7s}  {'New':>7s}  "
          f"{'Old err':>8s}  {'New err':>8s}")
    for f in display_freqs:
        if f > 1000:
            break
        t = float(np.interp(f, fit_freqs, diff_norm))
        o = float(np.interp(f, freqs_viz, old_norm))
        n = float(np.interp(f, freqs_viz, new_norm))
        print(f"  {f:6.0f} Hz  {t:7.4f}  {o:7.4f}  {n:7.4f}  "
              f"{o-t:+8.4f}  {n-t:+8.4f}")


if __name__ == "__main__":
    main()
