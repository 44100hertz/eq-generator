"""Quick linearity check for v2."""
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from harmonic_bass import BassEnhancerConfig, generate_test_sine
from fix_v2 import process_fixed_v2
from goertzel_analyzer import goertzel_magnitude


def run():
    """Run the v2 linearity check."""
    cfg = BassEnhancerConfig(cutoff=60.0, h2_amp=0.33, h3_amp=0.33, fs=44100.0, env_release=200.0)

    print("  V2 LINEARITY: 50 Hz, varying amplitude")
    print(f"  {'In dBFS':>8s}  {'H2 rel f':>10s}  {'H3 rel f':>10s}")
    print(f"  {'-'*32}")

    h2_vals, h3_vals = [], []
    for amp in [1.0, 0.5, 0.25, 0.125, 0.09]:
        pad = int(0.3 * cfg.fs)
        sig = generate_test_sine([50], [amp], 0.5, cfg.fs, stereo=False)
        stereo = np.zeros((2, len(sig) + 2*pad))
        stereo[0, pad:pad+len(sig)] = sig; stereo[1,:] = stereo[0,:]*0.95
        out = process_fixed_v2(stereo, cfg)
        ss = pad+len(sig)//2
        mono = (out[0,ss:pad+len(sig)] + out[1,ss:pad+len(sig)])/2
        f = goertzel_magnitude(mono, 50, cfg.fs)
        h2 = goertzel_magnitude(mono, 100, cfg.fs)
        h3 = goertzel_magnitude(mono, 150, cfg.fs)
        h2r = 20*np.log10(max(h2,1e-12)/max(f,1e-12))
        h3r = 20*np.log10(max(h3,1e-12)/max(f,1e-12))
        h2_vals.append(h2r); h3_vals.append(h3r)
        print(f"  {20*np.log10(amp):+8.2f}  {h2r:+10.2f}  {h3r:+10.2f}")

    print(f"  H2 σ={np.std(h2_vals):.2f} dB  H3 σ={np.std(h3_vals):.2f} dB  "
          f"{'✅' if np.std(h2_vals)<0.01 and np.std(h3_vals)<0.01 else '❌'}")


if __name__ == "__main__":
    run()
