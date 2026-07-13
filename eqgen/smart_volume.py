"""
Smart volume simulation: model the correction curve at different volume levels.

Matches the C firmware (smart_volume.h / filter.c) exactly:
  - Shelf boost: linear with dB drop from peak (0.8 dB per 10 dB, FM contours)
  - Pre-gain: linear interpolation pg_quiet → pg_loud
  - Volume LUT: attenuation mapped through speaker_level

Volume attenuation happens BEFORE the enhancer — this creates
headroom that the enhancer fills with bleed/harmonics.  The
effective curve shows the net perceived response.

The shelf follows FM equal-loudness contours linearly:
  shelf_db = 0.08 × drop_from_peak_db
where drop_from_peak_db = speaker_level × (1 − vol/127).

This keeps compensation proportional to the absolute listening level:
a speaker with 80 dB range drops 80 dB to vol=0 (≈20 dB SPL → 6.4 dB
shelf), while a 40 dB range speaker drops only 40 dB (≈20 dB SPL →
3.2 dB shelf).  The peak (vol=127) always gets 0 shelf.
"""

import numpy as np

from eqgen.dsp import first_order_lp_mag, pre_gain_from_max_gain

# Constants matching smart_volume.h / firmware
SV_LOUDNESS_FC = 200.0   # Hz — one-pole shelf corner
FM_SLOPE = 0.08           # dB shelf per dB of attenuation (0.8 dB / 10 dB)


def compute_smart_volume_curves(
    freqs: np.ndarray,
    correction_db: np.ndarray,
    pre_gain_db_loud: float = 0.0,
    speaker_level_db: int = 60,
    overboost_db: float = 0.0,
) -> dict:
    """Compute effective curves at different smart-volume levels.

    Returns dict with keys: curves, lut_db, params.
    """
    pg_loud_linear = 10.0 ** (pre_gain_db_loud / 20.0)

    # Volume LUT: linear-in-dB mapping from sv_db_floor → overboost_db
    sv_db_floor = float(-speaker_level_db)
    if sv_db_floor > -24.0:
        sv_db_floor = -24.0
    if sv_db_floor < -80.0:
        sv_db_floor = -80.0

    lut_range = -sv_db_floor  # 24–80 dB — the speaker's total dynamic range

    # Peak shelf at vol=0 = FM_SLOPE × lut_range (linear with total range)
    max_shelf_db = FM_SLOPE * lut_range
    max_shelf_linear = 10.0 ** (max_shelf_db / 20.0)

    # Precompute per-step data (0, 32, 64, 96, 127)
    steps = [0, 32, 64, 96, 127]
    curves = []
    lut_db = []

    for vol in steps:
        # Volume attenuation (matches smart_volume_rebuild_lut with compensation_db=0)
        vol_db = sv_db_floor + (float(vol) / 127.0) * lut_range + overboost_db
        if vol_db > overboost_db:
            vol_db = overboost_db

        # Smart volume: linear FM shelf (matches smart_volume_compute)
        t = float(vol) / 127.0
        drop_db = lut_range * (1.0 - t)  # dB below peak
        shelf_db_val = FM_SLOPE * drop_db

        # Pre-gain (matches smart_volume_compute)
        pg_quiet = pg_loud_linear / max_shelf_linear
        pre_gain_linear = pg_quiet + t * (pg_loud_linear - pg_quiet)
        pre_gain_db = 20.0 * np.log10(max(pre_gain_linear, 1e-12))

        # Shelf response: one-pole low shelf at 200 Hz
        if shelf_db_val > 0.01:
            shelf_gain = np.array([
                1.0 + (10.0 ** (shelf_db_val / 20.0) - 1.0) * first_order_lp_mag(f, SV_LOUDNESS_FC)
                for f in freqs
            ])
            shelf_db_arr = 20.0 * np.log10(np.maximum(shelf_gain, 1e-12))
        else:
            shelf_db_arr = np.zeros_like(freqs)

        # Effective: correction + shelf + pre_gain + volume attenuation.
        # Volume is applied BEFORE the enhancer — it creates headroom.
        effective_db = correction_db + shelf_db_arr + pre_gain_db + vol_db

        label = f"vol={vol}" if vol > 0 else "vol\u21920"
        curves.append({
            "label": label,
            "vol": vol,
            "t": round(t, 2),
            "vol_db": round(vol_db, 1),
            "pre_gain_db": round(pre_gain_db, 1),
            "shelf_max_db": round(shelf_db_val, 1),
            "freqs": freqs.tolist(),
            "effective_db": effective_db.tolist(),
        })
        lut_db.append(round(vol_db, 1))

    return {
        "curves": curves,
        "vol_steps": steps,
        "lut_db": lut_db,
        "params": {
            "fc": float(SV_LOUDNESS_FC),
            "fm_slope": float(FM_SLOPE),
            "shelf_max_db": round(max_shelf_db, 2),
            "pre_gain_db_loud": round(pre_gain_db_loud, 1),
            "speaker_level_db": speaker_level_db,
            "sv_db_floor": round(sv_db_floor, 1),
            "overboost_db": round(overboost_db, 1),
        },
    }
