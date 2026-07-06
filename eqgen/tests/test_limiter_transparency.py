#!/usr/bin/env python3
"""
Limiter transparency benchmark.

Processes track-2 from each album through the enhancer and measures:
  - gain_reduction_pct: % of samples where limiter is active (gain < 0.99)
  - mean_attenuation_db: average gain reduction in dB
  - modulation_depth_db: std dev of limiter gain (gain pumping)
  - heavy_limit_pct: % of samples with >6dB attenuation (gain < 0.5)

Low modulation depth = transparent. High heavy_limit_pct = audible.
"""

import sys, os, subprocess, time
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy.signal import lfilter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from eqgen.dsp import design_butter_lp, design_butter_hp, EnvFollower, biquad_tick


def ba(c):
    return (np.array([c[0], c[1], c[2]]), np.array([1.0, c[3], c[4]]))


def decode(path, sr=44100):
    p = subprocess.run(
        ["ffmpeg", "-y", "-v", "error", "-i", path,
         "-ar", str(sr), "-ac", "2", "-f", "f32le", "pipe:1"],
        capture_output=True,
    )
    return np.frombuffer(p.stdout, dtype=np.float32).reshape(-1, 2).T.astype(np.float64), sr


def measure_limiter(chunk, fs, cutoff, h):
    """Returns (gain_trace, reduction_pct, mean_att_db, mod_db, heavy_pct)."""
    b_lp, a_lp = ba(design_butter_lp(cutoff, fs))
    b_lp3, a_lp3 = ba(design_butter_lp(cutoff / 2, fs))
    b_hp, a_hp = ba(design_butter_hp(cutoff, fs))
    hp_c = design_butter_hp(cutoff, fs)

    L, R = chunk[0], chunk[1]
    N = len(L)

    lp2_L = lfilter(b_lp, a_lp, L)
    lp3_L = lfilter(b_lp3, a_lp3, L)
    wet_L = lfilter(b_hp, a_hp, L)

    env2 = EnvFollower.from_params(200, fs)
    env3 = EnvFollower.from_params(200, fs)
    lim_env = EnvFollower.from_params(200, fs)
    lim_gain = 1.0
    floor = 0.0001
    hp_st = np.zeros(4)
    gains = np.empty(N)

    for i in range(N):
        env2.tick(lp2_L[i])
        env3.tick(lp3_L[i])
        e2 = max(env2.read(), floor)
        e3 = max(env3.read(), floor)
        t2 = e2 * h * (2.0 * (lp2_L[i] / e2) ** 2 - 1.0)
        t3 = e3 * h * (4.0 * (lp3_L[i] / e3) ** 3 - 3.0 * (lp3_L[i] / e3))
        harm_hp, hp_st = biquad_tick(t2 + t3, hp_c, hp_st)
        total = (wet_L[i] + harm_hp) * lim_gain
        lim_env.tick(abs(total))
        lim_gain = 1.0 / max(lim_env.read(), 1.0)
        gains[i] = lim_gain

    active = gains < 0.99
    heavy = gains < 0.5
    att_db = -20 * np.log10(np.maximum(gains, 1e-10))

    return {
        "active_pct": np.mean(active) * 100,
        "mean_att_db": np.mean(att_db),
        "modulation_db": np.std(att_db),
        "heavy_pct": np.mean(heavy) * 100,
        "min_gain": gains.min(),
        "mean_gain": gains.mean(),
    }


def find_track2(music_root):
    exts = {".flac", ".mp3", ".wav", ".opus", ".ogg", ".m4a", ".aac", ".wv", ".ape"}
    albums = defaultdict(list)
    for root, dirs, files in os.walk(music_root):
        af = [f for f in files if Path(f).suffix.lower() in exts]
        if af:
            albums[root].extend(sorted(af))
    return [os.path.join(d, fs[1]) for d, fs in sorted(albums.items()) if len(fs) >= 2]


def main():
    music_root = sys.argv[1] if len(sys.argv) > 1 else (
        "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/"
    )
    cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
    h = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

    t2 = find_track2(music_root)
    print(f"Testing {len(t2)} albums at fc={cutoff:.0f} h={h:.2f}", file=sys.stderr)

    all_results = []
    t0 = time.time()

    for idx, path in enumerate(t2):
        try:
            samples, fs = decode(path)
            if samples.shape[1] < fs:
                continue
            mid = samples.shape[1] // 2
            start = max(0, mid - fs // 2)
            chunk = samples[:, start:start + fs]
            pk = max(abs(chunk[0]).max(), abs(chunk[1]).max())
            if pk < 0.01:
                continue
            chunk *= 1.0 / pk  # normalize to peak 1.0

            r = measure_limiter(chunk, fs, cutoff, h)
            all_results.append(r)

            elapsed = time.time() - t0
            if idx % 10 == 0:
                print(f"  [{idx+1}/{len(t2)} {elapsed:.0f}s] "
                      f"active={r['active_pct']:.0f}% "
                      f"att={r['mean_att_db']:.2f}dB "
                      f"mod={r['modulation_db']:.3f}dB "
                      f"heavy={r['heavy_pct']:.1f}%",
                      file=sys.stderr)

            if elapsed > 30:
                print(f"Stopping at {len(all_results)} tracks, {elapsed:.0f}s", file=sys.stderr)
                break

        except Exception as e:
            print(f"  SKIP: {e}", file=sys.stderr)

    # Report
    if not all_results:
        print("No data")
        return

    active = np.array([r["active_pct"] for r in all_results])
    att = np.array([r["mean_att_db"] for r in all_results])
    mod = np.array([r["modulation_db"] for r in all_results])
    heavy = np.array([r["heavy_pct"] for r in all_results])

    print(f"\n{'='*60}")
    print(f"  LIMITER TRANSPARENCY BENCHMARK  (fc={cutoff:.0f} h={h:.2f})")
    print(f"  {len(all_results)} tracks × 1s mid-chunk, peak-normalized")
    print(f"{'='*60}")
    print(f"  {'':20s}  {'mean':>8s}  {'median':>8s}  {'p90':>8s}  {'max':>8s}")
    print(f"  {'active %':20s}  {active.mean():8.1f}  {np.median(active):8.1f}  {np.percentile(active,90):8.1f}  {active.max():8.1f}")
    print(f"  {'attenuation dB':20s}  {att.mean():8.2f}  {np.median(att):8.2f}  {np.percentile(att,90):8.2f}  {att.max():8.2f}")
    print(f"  {'modulation dB':20s}  {mod.mean():8.3f}  {np.median(mod):8.3f}  {np.percentile(mod,90):8.3f}  {mod.max():8.3f}")
    print(f"  {'heavy (>6dB) %':20s}  {heavy.mean():8.1f}  {np.median(heavy):8.1f}  {np.percentile(heavy,90):8.1f}  {heavy.max():8.1f}")

    # Interpretation
    print(f"\n  Interpretation:")
    if np.median(heavy) < 1:
        print(f"    ✅ Heavy limiting is rare (<1% of samples)")
    else:
        print(f"    ⚠️  Heavy limiting on {np.median(heavy):.0f}% of samples — reduce h or raise cutoff")
    if np.median(mod) < 0.5:
        print(f"    ✅ Low modulation ({np.median(mod):.2f} dB) — gain changes are smooth")
    else:
        print(f"    ⚠️  High modulation ({np.median(mod):.2f} dB) — gain pumping may be audible")


if __name__ == "__main__":
    main()
