#!/usr/bin/env python3
"""Quick: test fixed h values on 30s of tracks. Write results as they come."""
import sys, os, subprocess, time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.signal import lfilter

sys.path.insert(0, str(Path(__file__).parent / "python"))
from dsp import design_butter_lp, design_butter_hp, EnvFollower, biquad_tick

def ba(c): return (np.array([c[0],c[1],c[2]]), np.array([1.0,c[3],c[4]]))

def decode(path, sr=44100):
    p = subprocess.run(["ffmpeg","-y","-v","error","-i",path,"-ar",str(sr),"-ac","2","-f","f32le","pipe:1"], capture_output=True)
    return np.frombuffer(p.stdout,dtype=np.float32).reshape(-1,2).T.astype(np.float64), sr

def frac_limited(chunk, fs, fc, h):
    b_lp,a_lp=ba(design_butter_lp(fc,fs))
    b_lp3,a_lp3=ba(design_butter_lp(fc/2,fs))
    b_hp,a_hp=ba(design_butter_hp(fc,fs))
    hp_c=design_butter_hp(fc,fs)
    L,R=chunk[0],chunk[1]; N=len(L)
    lp2_L=lfilter(b_lp,a_lp,L); lp2_R=lfilter(b_lp,a_lp,R)
    lp3_L=lfilter(b_lp3,a_lp3,L); lp3_R=lfilter(b_lp3,a_lp3,R)
    wet_L=lfilter(b_hp,a_hp,L); wet_R=lfilter(b_hp,a_hp,R)
    env2L=EnvFollower.from_params(200,fs); env2R=EnvFollower.from_params(200,fs)
    env3L=EnvFollower.from_params(200,fs); env3R=EnvFollower.from_params(200,fs)
    lim_env=EnvFollower.from_params(200,fs); lim_gain=1.0; floor=0.0001
    hpL=np.zeros(4); hpR=np.zeros(4); cnt=0
    for i in range(N):
        env2L.tick(lp2_L[i]); env2R.tick(lp2_R[i])
        env3L.tick(lp3_L[i]); env3R.tick(lp3_R[i])
        e2L=max(env2L.read(),floor); e2R=max(env2R.read(),floor)
        e3L=max(env3L.read(),floor); e3R=max(env3R.read(),floor)
        t2L=e2L*h*(2.0*(lp2_L[i]/e2L)**2-1.0); t2R=e2R*h*(2.0*(lp2_R[i]/e2R)**2-1.0)
        t3L=e3L*h*(4.0*(lp3_L[i]/e3L)**3-3.0*(lp3_L[i]/e3L)); t3R=e3R*h*(4.0*(lp3_R[i]/e3R)**3-3.0*(lp3_R[i]/e3R))
        hhL,hpL=biquad_tick(t2L+t3L,hp_c,hpL); hhR,hpR=biquad_tick(t2R+t3R,hp_c,hpR)
        totalL=(wet_L[i]+hhL)*lim_gain; totalR=(wet_R[i]+hhR)*lim_gain
        lim_env.tick(max(abs(totalL),abs(totalR)))
        lim_gain=1.0/max(lim_env.read(),1.0)
        if lim_gain<0.5: cnt+=1
    return cnt/N

def find_track2(mr):
    exts={".flac",".mp3",".wav",".opus",".ogg",".m4a",".aac",".wv",".ape"}
    al=defaultdict(list)
    for r,d,f in os.walk(mr):
        af=[x for x in f if Path(x).suffix.lower() in exts]
        if af: al[r].extend(sorted(af))
    return [os.path.join(d,fs[1]) for d,fs in sorted(al.items()) if len(fs)>=2]

def main():
    mr=sys.argv[1] if len(sys.argv)>1 else "/run/media/samp/787be337-88e4-4b95-92f9-45d37615cd02/music/mu/"
    t2=find_track2(mr); print(f"{len(t2)} albums",file=sys.stderr)
    cutoffs=[40,60,80]; h_vals=[0.5,1.0,1.5,2.0,3.0,5.0]
    # results[fc][h] = list of fractions
    results={fc:{h:[] for h in h_vals} for fc in cutoffs}
    t0=time.time()
    out=open("/tmp/h_sweep.csv","w")
    out.write("track,fc,h,frac_limited\n")
    for idx,path in enumerate(t2):
        try:
            samples,fs=decode(path)
            if samples.shape[1]<fs: continue
            mid=samples.shape[1]//2; start=max(0,mid-2098)
            chunk=samples[:,start:start+4196]
            pk=max(abs(chunk[0]).max(),abs(chunk[1]).max())
            if pk<0.01: continue
            chunk*=1.0/pk
            name=Path(path).stem[:40]
            for fc in cutoffs:
                for h in h_vals:
                    f=frac_limited(chunk,fs,fc,h)
                    results[fc][h].append(f)
                    out.write(f"{name},{fc},{h},{f:.6f}\n")
            out.flush()
            elapsed=time.time()-t0
            print(f"[{idx+1} {elapsed:.0f}s]",file=sys.stderr)
            if elapsed > 30 or idx >= 50:
                print(f"Stopping at {idx+1} tracks, {elapsed:.0f}s",file=sys.stderr)
                break
        except Exception as e:
            print(f"SKIP {Path(path).name[:30]}: {e}",file=sys.stderr)
    out.close()

    # Report
    print("\n"+"="*60)
    print("  Highest h where ≥90% of tracks have <10% lim_gain<0.5")
    print("="*60)
    for fc in cutoffs:
        print(f"\n  fc={fc:.0f} Hz:")
        best=0
        for h in h_vals:
            arr=np.array(results[fc][h])
            if len(arr)==0: continue
            passes=np.mean(arr<=0.10)
            print(f"    h={h:.2f}: pass={passes*100:.0f}% (n={len(arr)})")
            if passes>=0.90: best=h
        print(f"    → max safe: h={best:.2f}")

if __name__=="__main__":
    main()
