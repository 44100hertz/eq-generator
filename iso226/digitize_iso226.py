#!/usr/bin/env python3
"""Digitize "ISO226 curves.png" (equal-loudness contours) into a dB-SPL table.

Method:
  * Y axis (dB SPL) calibrated from the horizontal gridlines (130..0 dB) and
    the black x-axis line (-10 dB).
  * X axis (log frequency) calibrated from the centers of the x-axis tick
    labels (20, 40, ..., 10000 Hz).
  * The 11 contours (0,10,...,100 phon; 0 is the dashed T_f threshold curve)
    are seeded near 20 Hz, where they are stacked in order and well
    separated, then tracked across the plot by following dark near-neutral
    pixels with slope prediction. Dashed curves are coasted through gaps.
  * Red annotation text and the pink watermark are rejected by requiring
    dark pixels to be near-neutral (small max-min channel spread).

Outputs: printed table, ../eqgen/iso226_table.csv, iso226_verify.png (overlay).

Verified by: 1 kHz anchor (N phon must read N dB), strict curve ordering,
max tracking step, and visual inspection of the overlay.
"""

import csv

import numpy as np
from PIL import Image, ImageDraw

IMG_PATH = "ISO226 curves.png"
PHONS = list(range(0, 101, 10))
FREQS = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400,
         500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000,
         6300, 8000, 10000, 12500]

ROW_TOP, ROW_BOT = 100, 1004  # plot interior; the x-axis line sits at ~1009
WIN = 6          # search half-window (px) around predicted y while tracking
MAX_MISS = 130   # columns to coast through dash gaps / text overlaps


def clusters_1d(idxs, gap=2):
    """Group sorted indices into clusters, return their centers."""
    if len(idxs) == 0:
        return []
    out, start, prev = [], idxs[0], idxs[0]
    for i in idxs[1:]:
        if i > prev + gap:
            out.append((start + prev) / 2)
            start = i
        prev = i
    out.append((start + prev) / 2)
    return out


def main():
    im = np.array(Image.open(IMG_PATH).convert("RGB")).astype(int)
    H, W, _ = im.shape
    mx, mn = im.max(2), im.min(2)
    dark = (mx < 150) & (mx - mn < 60)          # black curves/axes/text
    gray = (mx - mn < 15) & (mx > 170) & (mx < 248)  # light-gray gridlines

    # --- axes -----------------------------------------------------------
    yaxis_col = int(np.argmax(dark.sum(0)[: W // 3]))
    xaxis_row = int(np.argmax(dark.sum(1)[H // 2:]) + H // 2)
    print(f"y-axis col {yaxis_col}, x-axis row {xaxis_row}")

    # --- Y calibration: gridlines are 130,120,...,0 dB ------------------
    row_hits = gray[:, yaxis_col + 5:].sum(1)
    rows = np.where(row_hits > 0.35 * (W - yaxis_col))[0]
    glines = clusters_1d(rows)
    assert len(glines) == 14, f"expected 14 horizontal gridlines, got {len(glines)}"
    db_vals = np.array([130 - 10 * i for i in range(14)])
    m, c = np.polyfit(db_vals, glines, 1)       # row = m*dB + c
    resid = np.abs(np.polyval([m, c], db_vals) - glines).max()
    axis_pred = np.polyval([m, c], -10)
    print(f"Y fit: {m:.4f} px/dB, max residual {resid:.2f} px, "
          f"predicted -10 dB row {axis_pred:.1f} vs axis {xaxis_row}")
    assert resid < 1.0 and abs(axis_pred - xaxis_row) < 2.0
    row_to_db = lambda r: (r - c) / m

    # --- X calibration: label centers are 20,40,...,10000 Hz ------------
    band = dark[xaxis_row + 16: xaxis_row + 56, yaxis_col:]
    cols = np.where(band.any(0))[0] + yaxis_col
    blobs = clusters_1d(cols, gap=15)[:10]
    label_freqs = np.array([20, 40, 80, 160, 315, 630, 1250, 2500, 5000, 10000])
    diffs = np.diff(blobs)
    assert len(blobs) == 10 and np.all((diffs > 120) & (diffs < 140)), \
        f"bad x labels: {blobs}"
    a, b = np.polyfit(np.log10(label_freqs), blobs, 1)  # col = a*log10(f)+b
    resid_x = np.abs(np.polyval([a, b], np.log10(label_freqs)) - blobs).max()
    print(f"X fit: {a:.2f} px/decade, max residual {resid_x:.2f} px")
    assert resid_x < 1.5
    f_to_col = lambda f: a * np.log10(f) + b
    col_to_f = lambda x: 10 ** ((x - b) / a)

    # --- per-column dark clusters (the curve pixels) --------------------
    col_clusters = [clusters_1d(np.where(dark[ROW_TOP:ROW_BOT, x])[0] + ROW_TOP)
                    for x in range(W)]

    # --- seed: first column near 20 Hz with all 11 contours in order ----
    seed = None
    for x in range(int(f_to_col(20)), int(f_to_col(26))):
        cl = col_clusters[x]
        if (len(cl) == 11 and min(np.diff(cl)) >= 12
                and cl[0] < 200 and cl[-1] > 400):
            seed = (x, cl)
            break
    assert seed, "no seed column with 11 ordered contours near 20 Hz"
    sx, scl = seed
    print(f"seed at col {sx} ({col_to_f(sx):.1f} Hz): "
          f"{[round(row_to_db(r), 1) for r in scl]} dB")

    # --- track each contour left and right from the seed ----------------
    def track(x0, y0, step):
        pts = []
        y, slope, miss = y0, 0.0, 0
        last_x, last_y = x0, y0
        x = x0 + step
        while yaxis_col + 2 <= x < W - 2 and miss < MAX_MISS:
            pred = y + slope
            win = min(WIN + 0.5 * miss, 14)  # widen while coasting dash gaps
            near = [cy for cy in col_clusters[x] if abs(cy - pred) <= win]
            if near:
                ynew = min(near, key=lambda cy: abs(cy - pred))
                raw = np.clip((ynew - last_y) / (x - last_x), -3, 3)
                slope = 0.6 * slope + 0.4 * raw
                y, last_x, last_y = ynew, x, ynew
                pts.append((x, ynew))
                miss = 0
            else:
                miss += 1
                slope *= 0.95  # coast asymptotically flat: a runaway slope
                y += slope     # can capture the contour below at a curve end
            x += step
        return pts

    tracks = {}  # phon -> (xs, ys) sorted by x
    for phon, y0 in zip(PHONS[::-1], scl):  # top cluster = 100 phon
        left = track(sx, y0, -1)
        right = track(sx, y0, +1)
        pts = left[::-1] + [(sx, y0)] + right
        xs = np.array([p[0] for p in pts])
        ys = np.array([p[1] for p in pts])
        tracks[phon] = (xs, ys)
        step_max = (np.abs(np.diff(ys) / np.diff(xs)).max()
                    if len(ys) > 1 else 0)
        print(f"{phon:3d} phon: {col_to_f(xs[0]):7.1f} Hz .. "
              f"{col_to_f(xs[-1]):7.1f} Hz, {len(xs)} pts, "
              f"max step {step_max:.2f} px/col")
        assert step_max < 4, f"{phon} phon: tracking jump {step_max} px/col"

    # --- truncate captures ----------------------------------------------
    # A contour that really ends leaves no pixels; if a long found-gap is
    # followed by a merge with another contour, everything past the gap is
    # a capture, not the curve. (Contours never touch: min separation in
    # this chart is ~25 px, so <3 px for many columns means capture.)
    for phon in PHONS:
        xs, ys = tracks[phon]
        for g in np.where(np.diff(xs) > 45)[0]:
            x2 = xs[g + 1]
            for other in PHONS:
                if other == phon:
                    continue
                xo, yo = tracks[other]
                lo, hi = max(x2, xo[0]), min(x2 + 80, xo[-1])
                if lo > hi:
                    continue
                xc = np.arange(lo, hi + 1)
                if np.abs(np.interp(xc, xs, ys) -
                          np.interp(xc, xo, yo)).min() < 3:
                    tracks[phon] = (xs[:g + 1], ys[:g + 1])
                    print(f"{phon} phon: curve ends at "
                          f"{col_to_f(xs[g]):.0f} Hz (ignored capture "
                          f"by {other} phon after a {xs[g+1]-xs[g]:.0f} px gap)")
                    break
            else:
                continue
            break

    # --- sample at standard frequencies ---------------------------------
    table = {}  # phon -> list of dB (nan if curve not drawn there)
    for phon in PHONS:
        xs, ys = tracks[phon]
        vals = []
        for f in FREQS:
            fx = f_to_col(f)
            if xs[0] <= fx <= xs[-1]:
                vals.append(row_to_db(np.interp(fx, xs, ys)))
            elif abs(fx - xs[0]) <= 4:  # curve starts/ends a hair outside
                vals.append(row_to_db(ys[0]))
            elif abs(fx - xs[-1]) <= 4:
                vals.append(row_to_db(ys[-1]))
            else:
                vals.append(float("nan"))
        table[phon] = vals

    # --- checks -----------------------------------------------------------
    i1k = FREQS.index(1000)
    # N phon must read N dB at 1 kHz. Exempt: 0 phon (it is the measured T_f
    # threshold, ~2-4 dB at 1 kHz, not the definition) and 100 phon (the
    # chart only draws it at low frequencies).
    for phon in PHONS[1:-1]:
        err = abs(table[phon][i1k] - phon)
        assert err < 1.0, f"1 kHz anchor: {phon} phon reads {table[phon][i1k]}"
    worst = 99.0
    for i, f in enumerate(FREQS):
        col_vals = [(p, table[p][i]) for p in PHONS if not np.isnan(table[p][i])]
        for (p1, v1), (p2, v2) in zip(col_vals, col_vals[1:]):
            worst = min(worst, v2 - v1)
    print(f"1 kHz anchor OK; min gap between adjacent contours {worst:.2f} dB")
    assert worst > 0.3, "contours cross or touch"

    # --- output -----------------------------------------------------------
    header = ["freq_hz"] + [f"phon_{p}" for p in PHONS]
    with open("../eqgen/iso226_table.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for i, f in enumerate(FREQS):
            w.writerow([f] + ["" if np.isnan(table[p][i])
                              else f"{table[p][i]:.1f}" for p in PHONS])

    print("\nISO 226 equal-loudness contours digitized from chart "
          "(dB SPL, ~±1 dB reading accuracy)\n")
    print(f"{'Hz':>7} " + " ".join(f"{p:>6}" for p in PHONS))
    for i, f in enumerate(FREQS):
        cells = ["    — " if np.isnan(table[p][i]) else f"{table[p][i]:6.1f}"
                 for p in PHONS]
        print(f"{f:>7} " + " ".join(cells))

    # --- verification overlay ---------------------------------------------
    out = Image.open(IMG_PATH).convert("RGB")
    dr = ImageDraw.Draw(out)
    for phon in PHONS:
        xs, ys = tracks[phon]
        dr.line(list(zip(xs, ys)), fill=(0, 200, 0), width=1)
        for i, f in enumerate(FREQS):
            if not np.isnan(table[phon][i]):
                fx, fy = f_to_col(f), np.interp(f_to_col(f), xs, ys)
                dr.ellipse([fx - 3, fy - 3, fx + 3, fy + 3], fill=(255, 0, 0))
    out.save("iso226_verify.png")
    print("\nwrote ../eqgen/iso226_table.csv and iso226_verify.png")


if __name__ == "__main__":
    main()
