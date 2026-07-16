"""Test ISO 226:2023 contour loading, interpolation, and ear_sensitivity."""

import csv
import sys
import numpy as np
from pathlib import Path

ROOT = str(Path(__file__).resolve().parent.parent.parent)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from eqgen.iso226 import contour_spl_db, sensitivity


# ── Helpers ─────────────────────────────────────────────────────────────

def _read_csv_raw() -> dict[float, dict[int, float]]:
    """Re-read the CSV independently so tests don't trust the module loader."""
    csv_path = Path(__file__).resolve().parent.parent / "iso226_table.csv"
    result: dict[float, dict[int, float]] = {}
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        phon_cols: dict[int, int] = {}
        for i, h in enumerate(header):
            if h.startswith("phon_"):
                p = int(h.split("_")[1])
                if p <= 80:
                    phon_cols[i] = p
        for row in reader:
            freq = float(row[0])
            result[freq] = {}
            for col, phon in phon_cols.items():
                val = row[col].strip() if col < len(row) else ""
                if val:
                    result[freq][phon] = float(val)
    return result


_RAW = _read_csv_raw()
_FREQ_LIST = sorted(_RAW.keys())


# ── Anchor tests ────────────────────────────────────────────────────────

def test_1khz_anchor():
    """At 1 kHz, each phon contour should read close to the phon value."""
    for phon in [0, 10, 20, 30, 40, 50, 60, 70, 80]:
        spl = contour_spl_db(1000.0, phon)
        # The CSV 1 kHz values are within ~0.1 dB of phon (by design).
        expected = _RAW[1000.0][phon]
        assert abs(spl - expected) < 1e-6, f"1 kHz {phon} phon: {spl} vs CSV {expected}"


def test_sensitivity_1khz_is_one():
    """sensitivity(1000) must equal 1.0 for any phon."""
    for phon in [0, 20, 40, 60, 80]:
        s = sensitivity(1000.0, phon)
        assert abs(s - 1.0) < 1e-9, f"sensitivity(1000, {phon}) = {s} != 1.0"


# ── Spot-checks against CSV ─────────────────────────────────────────────

def test_csv_exact_match_at_table_frequencies():
    """At every table row, contour_spl_db must match the CSV value exactly."""
    for freq in _FREQ_LIST:
        for phon, csv_val in _RAW[freq].items():
            val = contour_spl_db(freq, phon)
            assert abs(val - csv_val) < 1e-9, f"{freq} Hz, {phon} phon: {val} vs {csv_val}"


def test_interpolation_between_phon_contours():
    """Mid-phon interpolation should split the difference linearly."""
    # At 1 kHz: phon_30 = 29.9, phon_40 = 40.0
    val_35 = contour_spl_db(1000.0, 35.0)
    assert abs(val_35 - 34.95) < 1e-6, f"35 phon at 1 kHz should be ~34.95, got {val_35}"


def test_interpolation_between_frequencies():
    """Mid-frequency interpolation in log space."""
    # 50 Hz, 60 phon = 90.1  (from CSV row 50)
    # 63 Hz, 60 phon = 86.0
    # 56.5 Hz should be roughly halfway in log between 50 and 63.
    val = contour_spl_db(56.5, 60.0)
    csv_50 = _RAW[50.0][60]   # 90.1
    csv_63 = _RAW[63.0][60]   # 86.0
    # log-space t
    t = (np.log10(56.5) - np.log10(50)) / (np.log10(63) - np.log10(50))
    expected = csv_50 + t * (csv_63 - csv_50)
    assert abs(val - expected) < 1e-6, f"56.5 Hz phon_60: {val} vs expected {expected}"


# ── Monotonicity ────────────────────────────────────────────────────────

def test_monotonic_phon():
    """Higher phon → higher SPL at every frequency."""
    for freq in _FREQ_LIST:
        vals = [(p, contour_spl_db(freq, p)) for p in sorted(_RAW[freq].keys())]
        for (p1, s1), (p2, s2) in zip(vals, vals[1:]):
            assert s2 > s1, f"{freq} Hz: phon {p2} SPL ({s2}) not > phon {p1} ({s1})"


def test_sensitivity_monotonic():
    """sensitivity decreases as frequency drops below 1 kHz."""
    s1000 = sensitivity(1000.0, 60.0)
    s400 = sensitivity(400.0, 60.0)
    s100 = sensitivity(100.0, 60.0)
    s50 = sensitivity(50.0, 60.0)
    assert s1000 > s400 > s100 > s50, (
        f"sensitivity not monotonic: 1k={s1000:.4f} 400={s400:.4f} "
        f"100={s100:.4f} 50={s50:.4f}"
    )


# ── Edge cases ──────────────────────────────────────────────────────────

def test_clamp_below_table():
    """Frequencies below 20 Hz clamp to the 20 Hz value."""
    spl20 = contour_spl_db(20.0, 60.0)
    spl10 = contour_spl_db(10.0, 60.0)
    spl5  = contour_spl_db(5.0, 60.0)
    assert abs(spl10 - spl20) < 1e-9
    assert abs(spl5 - spl20) < 1e-9


def test_clamp_above_table():
    """Frequencies above 12.5 kHz clamp to the 12.5 kHz value."""
    spl12500 = contour_spl_db(12500.0, 60.0)
    spl16000 = contour_spl_db(16000.0, 60.0)
    assert abs(spl16000 - spl12500) < 1e-9


def test_clamp_phon():
    """Phon outside [0, 80] is clamped."""
    s0 = contour_spl_db(1000.0, 0.0)
    s_neg = contour_spl_db(1000.0, -10.0)
    s80 = contour_spl_db(1000.0, 80.0)
    s_high = contour_spl_db(1000.0, 200.0)
    assert abs(s_neg - s0) < 1e-9
    assert abs(s_high - s80) < 1e-9


def test_no_phon_above_80_loaded():
    """Columns phon_90 and phon_100 must not be in the loaded contours."""
    from eqgen import iso226
    for phon in iso226._CONTOURS:
        assert phon <= 80, f"phon_{phon} should not be loaded"


# ── ear_sensitivity backward compat ─────────────────────────────────────

def test_ear_sensitivity_delegates():
    """eqgen.model.ear_sensitivity must match eqgen.iso226.sensitivity(f, 60)."""
    from eqgen.model import ear_sensitivity
    from eqgen.iso226 import sensitivity as iso_sens
    for f in [20, 50, 100, 400, 1000, 4000, 8000, 12500]:
        assert abs(ear_sensitivity(f) - iso_sens(f, 60.0)) < 1e-12


def test_sensitivity_above_1khz():
    """Ear is more sensitive at 3-4 kHz than at 1 kHz (sensitivity > 1)."""
    s = sensitivity(3150.0, 60.0)
    assert s > 1.0, f"sensitivity at 3150 Hz should be > 1, got {s:.4f}"


# ── Benchmark-related sanity ────────────────────────────────────────────

def test_sensitivity_range():
    """sensitivity stays in a reasonable range."""
    for f in np.logspace(np.log10(20), np.log10(12500), 100):
        s = sensitivity(f, 60.0)
        assert s > 0.0, f"sensitivity at {f:.1f} Hz is zero"
        assert s < 3.0, f"sensitivity at {f:.1f} Hz = {s:.4f} (suspiciously high)"


# ── Runner ──────────────────────────────────────────────────────────────

def run():
    tests = [
        test_1khz_anchor,
        test_sensitivity_1khz_is_one,
        test_csv_exact_match_at_table_frequencies,
        test_interpolation_between_phon_contours,
        test_interpolation_between_frequencies,
        test_monotonic_phon,
        test_sensitivity_monotonic,
        test_clamp_below_table,
        test_clamp_above_table,
        test_clamp_phon,
        test_no_phon_above_80_loaded,
        test_ear_sensitivity_delegates,
        test_sensitivity_above_1khz,
        test_sensitivity_range,
    ]
    failed = 0
    for t in tests:
        try:
            t()
        except Exception as e:
            print(f"FAIL {t.__name__}: {e}")
            failed += 1
        else:
            print(f"PASS {t.__name__}")
    print(f"\n{failed}/{len(tests)} failed")
    return failed


if __name__ == "__main__":
    sys.exit(run())
