"""
Show how the H2_AMP/H3_AMP change inverts the model's gain computation.
The model solves G = target / sqrt(a + b + c). Smaller h2/h3 → smaller
denom → LARGER G → hotter bass. But the enhancer itself got quieter,
so the old model constants were compensating for a LOUDER enhancer.
"""

import numpy as np

def hp(f, fc): w=f/fc; return w*w/np.sqrt(1+w**4)
def lp(f, fc): w=f/fc; return 1/np.sqrt(1+w**4)

def model_gain(f, fc, h2, h3, S_at_f, S_at_2f, S_at_3f, target=1.0):
    a = hp(f,fc)**2 * S_at_f**2
    b = h2**2 * lp(f,fc)**2 * hp(2*f,fc)**2 * S_at_2f**2
    c = h3**2 * lp(f,fc/2)**2 * hp(3*f,fc)**2 * S_at_3f**2
    G = target / np.sqrt(a+b+c) if a+b+c>1e-12 else 1.0
    comp = target / S_at_f if S_at_f>1e-12 else 1.0
    return min(G, comp)

def run():
    """Run the model gain analysis."""
    fc = 60.0
    freqs = [30, 40, 50, 60, 80, 100, 120]

    print("Speaker: -12 dB at 50 Hz, flat > 100 Hz")
    print(f"Cutoff: {fc} Hz")
    print()
    print(f"{'Freq':>6s}  {'Old G':>8s} {'Old dB':>8s}  {'New G':>8s} {'New dB':>8s}  {'Delta':>8s}")
    print(f"{'':->50s}")

    for f in freqs:
        Sf = speaker(f)
        S2f = speaker(2*f) if 2*f < 20000 else 1.0
        S3f = speaker(3*f) if 3*f < 20000 else 1.0
        G_old = model_gain(f, fc, 1.0, 2.0, Sf, S2f, S3f)
        G_new = model_gain(f, fc, 0.33, 0.33, Sf, S2f, S3f)
        dB_old = 20*np.log10(G_old)
        dB_new = 20*np.log10(G_new)
        print(f"  {f:4.0f} Hz  {G_old:8.3f} {dB_old:+7.2f}  {G_new:8.3f} {dB_new:+7.2f}  {dB_new-dB_old:+7.2f}")

    # What h2/h3 give the SAME gains as the old model?
    print(f"\n{'='*50}")
    print("What h2/h3 produce the SAME gain as old model?")
    print("(Old model was tuned with h2=1.0, h3=2.0)")
    print(f"{'='*50}")
    print(f"{'Freq':>6s}  {'Old G':>8s}  ", end="")
    for ratio in [0.33, 0.5, 0.7, 1.0]:
        print(f"{'h='+str(ratio):>8s}", end="  ")
    print()
    for f in [30, 40, 50, 60, 80]:
        Sf = speaker(f)
        S2f = speaker(2*f) if 2*f < 20000 else 1.0
        S3f = speaker(3*f) if 3*f < 20000 else 1.0
        G_old = model_gain(f, fc, 1.0, 2.0, Sf, S2f, S3f)
        print(f"  {f:4.0f} Hz  {G_old:8.3f}  ", end="")
        for ratio in [0.33, 0.5, 0.7, 1.0]:
            G = model_gain(f, fc, ratio*2.0, ratio*2.0, Sf, S2f, S3f)
            print(f"{G:8.3f}  ", end="")
        print()

    # But wait - the actual enhancer now uses h2=0.33, h3=0.33.
    # The correct model constants ARE 0.33, 0.33 - the model math is right.
    # The issue is: the OLD model (1.0, 2.0) was WRONG for the old enhancer
    # but happened to produce gains that worked well with real music.
    # Now the model is CORRECT but the gains are larger.
    #
    # Solution: either retune h2/h3 in the EEL plugin to be louder
    # (so the model denominator gets larger again), or accept that the
    # model needs different target behavior below cutoff.


# Simulate a small speaker: -12 dB at 50 Hz, flat above 100 Hz
def speaker(f): return 0.25 if f <= 50 else (1.0 if f >= 100 else 0.25 + 0.75*(f-50)/50)


if __name__ == "__main__":
    run()
