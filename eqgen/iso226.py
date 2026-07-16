"""
ISO 226:2023 equal-loudness contours, loaded from iso226_table.csv.

Provides:
    contour_spl_db(freq_hz, phon)  — dB SPL at a given frequency and phon level.
    sensitivity(freq_hz, phon=60)  — relative ear sensitivity (1.0 at 1 kHz).

Phon columns above 80 are ignored (user instruction).
"""

import csv
import numpy as np
from pathlib import Path


_TABLE_PATH = Path(__file__).resolve().parent / "iso226_table.csv"
_MAX_PHON = 80


def _load_table() -> dict[int, tuple[list[float], list[float]]]:
    """Return {phon: (freqs[], spl_db[])} for each phon column ≤ 80.

    Missing cells (empty CSV fields) are skipped, so each contour may have
    a different set of frequency points.
    """
    contours: dict[int, list[float | None]] = {}
    with open(_TABLE_PATH, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        col_map: dict[int, int] = {}
        for i, h in enumerate(header):
            if h.startswith("phon_"):
                p = int(h.split("_")[1])
                if p <= _MAX_PHON:
                    col_map[i] = p
                    contours[p] = []

        freqs: list[float] = []
        for row in reader:
            freqs.append(float(row[0]))
            for col, phon in col_map.items():
                val = row[col].strip() if col < len(row) else ""
                contours[phon].append(float(val) if val else None)

    # Strip missing entries from each contour.
    result: dict[int, tuple[list[float], list[float]]] = {}
    for phon, vals in contours.items():
        f_ok, s_ok = [], []
        for f_, s_ in zip(freqs, vals):
            if s_ is not None:
                f_ok.append(f_)
                s_ok.append(s_)
        result[phon] = (f_ok, s_ok)
    return result


_CONTOURS = _load_table()
_AVAILABLE_PHONS = sorted(_CONTOURS.keys())


def _interp_log_freq(freq_hz: float, fs: list[float], ss: list[float]) -> float:
    """Log-frequency interpolation within a single phon contour. Flat outside range."""
    if freq_hz <= fs[0]:
        return ss[0]
    if freq_hz >= fs[-1]:
        return ss[-1]
    # Binary search the bracket.
    lo, hi = 0, len(fs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if fs[mid] <= freq_hz:
            lo = mid
        else:
            hi = mid
    t = (np.log10(freq_hz) - np.log10(fs[lo])) / (np.log10(fs[hi]) - np.log10(fs[lo]))
    return ss[lo] + t * (ss[hi] - ss[lo])


def contour_spl_db(freq_hz: float, phon: float) -> float:
    """SPL in dB SPL at *freq_hz* for the given *phon* contour.

    Interpolation: log₁₀(f) between adjacent rows of the table, linear in
    phon between the bracketing 10-phon contours.

    Phon is clamped to [0, 80].  Frequency is clamped to the table extent
    (20 Hz – 12.5 kHz).
    """
    freq = float(freq_hz)
    p = max(0.0, min(float(phon), float(_MAX_PHON)))

    # Collect each contour's SPL at this frequency.
    spl_by_phon: dict[int, float] = {}
    for ph, (fs, ss) in _CONTOURS.items():
        spl_by_phon[ph] = _interp_log_freq(freq, fs, ss)

    # Linear interpolation in phon.
    p_list = sorted(spl_by_phon.keys())
    s_list = [spl_by_phon[p_] for p_ in p_list]
    return float(np.interp(p, p_list, s_list))


def sensitivity(freq_hz: float, phon: float = 60.0) -> float:
    """Ear sensitivity at *freq_hz* relative to 1 kHz = 1.0.

    Values < 1 mean the ear is LESS sensitive (needs more SPL to achieve
    the same perceived loudness).  Based on ISO 226:2023 contours.
    """
    spl_f = contour_spl_db(freq_hz, phon)
    spl_1k = contour_spl_db(1000.0, phon)
    return float(10.0 ** ((spl_1k - spl_f) / 20.0))
