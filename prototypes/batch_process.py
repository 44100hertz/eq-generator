"""Quick batch processor — runs end_to_end.py design once, then processes specific tracks."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'python'))

from end_to_end import load_measurement, design_eq, process_track
import numpy as np

TRACKS = [
    ("Daft Punk - One More Time",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/(2022) Daft Punk - Discovery (Ada - 0190296617164)/A1 One More Time.flac"),
    ("Radiohead - 15 Step",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/2007 - In Rainbows (UK 1st press 45RPM - _Xurbia_Xendless Limited _X_X001) (PBTHAL)/01 - 15 Step.flac"),
    ("Radiohead - Nude",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/2007 - In Rainbows (UK 1st press 45RPM - _Xurbia_Xendless Limited _X_X001) (PBTHAL)/03 - Nude.flac"),
    ("Radiohead - Reckoner",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/2007 - In Rainbows (UK 1st press 45RPM - _Xurbia_Xendless Limited _X_X001) (PBTHAL)/07 - Reckoner.flac"),
    ("Daft Punk - Something About Us",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/(2022) Daft Punk - Discovery (Ada - 0190296617164)/C2 Something About Us.flac"),
    ("Radiohead - Weird Fishes",
     "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/2007 - In Rainbows (UK 1st press 45RPM - _Xurbia_Xendless Limited _X_X001) (PBTHAL)/04 - Weird FishesArpeggi.flac"),
]

SPEAKERS = ["technics/standing", "cardboard"]

fc, h2, h3 = 60.0, 0.33, 0.33
start_sec, duration = 60, 45

for spk in SPEAKERS:
    print(f"\n{'='*60}")
    print(f"  SPEAKER: {spk}")
    print(f"{'='*60}")

    meas_dir = os.path.join(os.path.dirname(__file__), "measurements", spk)
    import glob
    meas_wavs = sorted(glob.glob(os.path.join(meas_dir, "measurement*.wav")))
    if not meas_wavs:
        meas_wavs = sorted(glob.glob(os.path.join(meas_dir, "response*.wav")))
    target_path = os.path.join(meas_dir, "target.wav")

    if not meas_wavs or not os.path.exists(target_path):
        print(f"  SKIP: missing measurements")
        continue

    # Design EQ once
    speaker_fn, fs, _, _, _ = load_measurement(meas_wavs[-1], target_path)
    coeffs, bands, _, _, _ = design_eq(speaker_fn, fs, fc=fc, h2=h2, h3=h3, max_bands=8)
    print(f"  EQ: {len(bands)} bands")

    out_dir = f"/tmp/eqgen_{spk.replace('/', '_')}"
    os.makedirs(out_dir, exist_ok=True)

    for label, path in TRACKS:
        if not os.path.exists(path):
            print(f"  SKIP: {label} — missing")
            continue
        safe = label.replace(" ", "_").replace("/", "-")[:50]
        out = os.path.join(out_dir, f"{safe}.wav")
        print(f"  {label}...")
        process_track(path, out, coeffs, len(bands),
                      cutoff_hz=fc, h2=h2, h3=h3,
                      start_sec=start_sec, duration_sec=duration)

    print(f"  → {out_dir}/")
