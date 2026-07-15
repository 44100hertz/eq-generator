# eqgen — Speaker EQ Correction with Psychoacoustic Bass Enhancement

A rapid prototyping suite for designing fixed‑point IIR EQ curves with harmonic
bass enhancement.  Measure a speaker on the desktop, tune parameters, audition
instantly, then **export Q4.28 coefficients** for deployment on an ESP32
running DSP + Bluetooth A2DP → I2S DAC → line out.

The DSP core in `src/` is the same C code that ships on the embedded target —
the desktop tools exist so you never have to tune in C.

---

## How it works

Small speakers can't reproduce deep bass.  Instead of wasting power trying, the
harmonic bass enhancer **cuts** real sub‑bass from the signal and **synthesizes**
2nd and 3rd harmonics from what remains.

The ear is ~5× more sensitive at 100–150 Hz than at 50 Hz (ISO 226 equal‑loudness
contours), so the perceived loudness can match the original — with less cone
excursion, less distortion, and lower power draw.

EQ correction is applied *before* the enhancer so the combined output (dry +
harmonics) matches the target frequency response.

### Pipeline (desktop)

```
measurement.wav ─┐
                 ├─→ Welch FFT (multi‑window, per‑bin CV) ─→ adaptive EQ points
target.wav ──────┘         │                                      │
                           │  noise.wav (spectral subtraction)     │
                           ▼                                      ▼
                      pooled stats ────────────────────→ resample to EQ grid
                                                                    │
                                                     ┌──────────────┤
                                                     ▼              ▼
                                               target curve   measurement curve
                                                     │              │
                                                     ▼              ▼
                                              correction = target / measurement
                                                     │
                                                     ▼ (optional)
                                              bass enhancer preprocess
                                              (solve G s.t. perceived output = target)
                                                     │
                                                     ▼
                                              IIR biquad fit (greedy golden‑section)
                                                     │
                                                     ▼
                                              Q4.28 quantization
                                                    ╱          ╲
                                                   ▼            ▼
                                            JamesDSP CSV     enhancer.so
                                            (desktop/mobile)  (ctypes FFI)
                                                                 │
                                                  ┌──────────────┴──────────────┐
                                                  ▼                             ▼
                                          process_music.py               wire_eq.py
                                          (offline WAV)              (live PipeWire LADSPA)
```

### DSP core (runs identically on desktop and ESP32)

```
  input (Q16)
     │
     ▼
  DC blocker
     │
     ▼
  cascaded EQ biquads (Q4.28)
     │
     ├──→ LP(fc) ──→ envelope ──→ normalize ──→ T₂ ──→ scale ──┐
     │   (T₂ path: 2nd harmonic from all bass below fc)           │
     │                                                            │
     └──→ LP(fc/2) ──→ envelope ──→ normalize ──→ T₃ ──→ scale ─┤
         (T₃ path: 3rd harmonic from ultra‑bass below fc/2)      │
                                                                  ▼
                                                            harmonics mix
                                                                  │
     HP(fc) ←── original dry signal ──────────────────────────────┤
                                                                  ▼
                                                            output (Q16)
     Limited to 1.0/max(env,1.0) on harmonics only — dry path untouched.
```

---

## Directory overview

```
eqgen/
├── eqgen.py                 # EQ design → JamesDSP CSV (desktop + mobile)
├── end_to_end.py            # Full pipeline: speaker → IIR → C DSP → WAV
├── process_music.py         # Standalone C enhancer (no measurement, flat EQ)
├── wire_eq.py / wire_eq.sh  # Live system‑wide EQ via PipeWire LADSPA
├── export_coeffs.py         # Export fitted biquads as C header for ESP32
├── sanitycheck.sh           # Quick smoke test
├── batch_process.py         # Batch process hardcoded tracks
├── compare_fitters.py       # Compare IIR fitting strategies
├── test_fullband_fit.py     # Full‑band (20 Hz–20 kHz) IIR fit tests
├── test_real_fullband_fit.py
├── test_sine_sweep.py
│
├── python/                  # Pure Python modules
│   ├── model.py             # Psychoacoustic model (sine → perceived amplitude)
│   ├── eq_fit.py            # Greedy IIR peaking‑filter fitter
│   ├── quantize.py          # Q4.28 fixed‑point quantization + ReciprocalLUT
│   ├── dsp.py               # Python simulation of the C enhancer (4 variants)
│   ├── analysis.py          # Goertzel, FFT analysis, test‑signal generation
│   ├── enhancer_ffi.py      # ctypes wrapper for src/enhancer.so
│   ├── audio_io.py          # WAV I/O and measurement loading
│   └── run_all.py           # Test runner
│
├── src/                     # C DSP — same code that ships on ESP32
│   ├── enhancer.c / .h      # Harmonic bass enhancer (Q16 fixed‑point)
│   ├── enhancer_api.c / .h  # C API used by Python via ctypes
│   ├── biquad.h             # Float biquad filter
│   ├── envelope.h           # Envelope follower with reciprocal LUT
│   ├── eq_coeffs.h          # Generated EQ coefficients (autogenerated)
│   ├── eqgen_filter.c       # Native PipeWire filter process
│   ├── ladspa_wrapper.c     # LADSPA plugin wrapper
│   ├── test_enhancer.c      # Standalone C test
│   └── Makefile
│
├── measurements/            # Per‑speaker recording directories
│   └── <speaker_name>/
│       ├── measurement*.wav # Brown‑noise recording of the speaker
│       ├── target.wav       # Brown‑noise recording of reference system
│       ├── noise.wav        # Ambient noise floor (optional, for subtraction)
│       └── *.txt            # Legacy FFT‑dump measurements
│
├── flake.nix                # Nix dev shell (Python + numpy/scipy/matplotlib)
├── harmonic_bass_enhancer*.eel  # Original REAPER/JSFX prototypes
└── README.md
```

---

## Taking measurements

### Equipment
- A measurement microphone (flat response, or one you've compensated for)
- The speaker under test
- A **reference system** — flat‑response speakers or known‑good headphones
- Audacity (or `sox`) to generate brown noise

### Procedure

1. **Generate brown noise** — 30–60 seconds at a comfortable listening level.
   In Audacity: Generate → Noise → Brownian, export as 16‑bit WAV.

2. **Record the reference system** playing the noise → `target.wav`.
   Place the mic where your ears would be.  This captures the room + reference
   system transfer function.

3. **Record the speaker under test** playing the same noise → `measurement.wav`
   (or `measurement2.wav`, etc. for multiple takes).  Same mic position, same
   gain settings.

4. **Record ambient silence** → `noise.wav` (optional).  Used for spectral
   subtraction if your room has audible background noise.

5. **Place files** in `measurements/<speaker_name>/`.

### What the pipeline does with them

The target recording represents what you *want* the speaker to sound like.
The measurement recording is what the speaker actually sounds like.
The correction curve = `target / measurement` in linear magnitude.

If you skip the target and just want a flat EQ, record the reference system
once and reuse that target for all speakers.

---

## Usage — typical workflow

### 1. Design the EQ curve

```bash
python eqgen.py \
  -m measurements/technics/standing/measurement2.wav \
  -t measurements/technics/standing/target.wav \
  --noise measurements/technics/standing/noise2.wav \
  --bass-enhancer-cutoff 50 \
  --h2 0.5 --h3 1.0 \
  -o my_eq
```

Produces `my_eq_desktop.csv` (tab‑separated) and `my_eq_mobile.csv` (GraphicEQ
string).  Both are normalised so the maximum gain is 0 dB — ready to import
into JamesDSP.

Without bass enhancement:

```bash
python eqgen.py -m meas.wav -t target.wav -o flat_eq
```

### 2. Audition offline

Process real music through the full pipeline:

```bash
python end_to_end.py technics/standing /tmp/out --tracks song1.flac song2.mp3
```

Or use the C enhancer standalone (flat EQ, no measurement):

```bash
python process_music.py ~/Music /tmp/out
```

### 3. Audition live (Linux + PipeWire)

Wire the enhancer into your system audio path:

```bash
./wire_eq.sh setup technics/standing
```

All system audio now flows through:  app → null sink → native DSP filter → hardware output.

Tear down:

```bash
./wire_eq.sh teardown
```

### 4. Iterate

Tweak cutoff, harmonic amplitudes, noise threshold, number of bands — re‑run,
re‑listen until it sounds right.

### 5. Export for ESP32

When satisfied, bake the coefficients into a C header:

```bash
python export_coeffs.py --speaker small --fc 60 --h2 0.33 --h3 0.33 -o src/eq_coeffs.h
make -C src
```

The DSP code in `src/enhancer.c` + `src/biquad.h` + `src/envelope.h`
compiles directly under ESP‑IDF — no changes needed.

---

## Entry point reference

| Script | Purpose |
|---|---|
| `eqgen.py` | Measurements → JamesDSP CSV (desktop + mobile).  Core tuning tool. |
| `end_to_end.py` | Speaker → IIR fit → C DSP → WAV.  Full offline audition. |
| `process_music.py` | Run the C enhancer on audio files (flat EQ, no measurements). |
| `wire_eq.py` / `wire_eq.sh` | Live PipeWire LADSPA system‑wide EQ.  Real‑time audition. |
| `export_coeffs.py` | Export fitted biquads as `src/eq_coeffs.h` for ESP32 firmware. |
| `sanitycheck.sh` | Quick smoke test with bundled technics/standing measurements. |
| `batch_process.py` | Batch‑process hardcoded tracks across multiple speakers. |
| `compare_fitters.py` | Compare IIR fitting strategies (historical/tuning tool). |

---

## Parameter guide

| Parameter | CLI flag | Default | Range | Description |
|---|---|---|---|---|
| **Cutoff** | `--bass-enhancer-cutoff` / `--fc` | none | 20–120 Hz | Crossover frequency. Below this, harmonics are synthesized instead of the fundamental. |
| **h2 amp** | `--h2` | 0.5–1.0† | 0.0–1.0 | 2nd‑harmonic amplitude. Controls perceived bass warmth. |
| **h3 amp** | `--h3` | 0.33–1.0† | 0.0–1.0 | 3rd‑harmonic amplitude. Controls aggressive bass edge. |
| **Max noise** | `--max-noise` | 0.65 | 0.05–2.0 | Per‑span CV threshold. Lower = more EQ points (needs cleaner data). Higher = fewer, broader points. |
| **Boost** | `-b` / `--boost` | 0 | −6–+6 dB | Global EQ gain applied to all bands. |
| **Max bands** | `--max-bands` | 24 | 1–40 | Maximum IIR peaking biquads. Fewer = simpler, more = tighter fit. |
| **High rolloff** | `--high-rolloff FREQ,DB` | none | e.g. `4000,-6` | Shelf rolloff on the target curve above FREQ. Repeatable. |
| **Low rolloff** | `--low-rolloff FREQ,DB` | none | e.g. `100,-12` | Shelf rolloff on the target curve below FREQ. Repeatable. |
| **Release** | `--release` | 0.2 s | 0.01–2.0 s | Envelope follower release time for harmonic linearisation. |
| **Limiter release** | (in enhancer_ffi) | 0.049 s | 0.01–0.5 s | Harmonic limiter release. Shorter = tighter peak control. |

†  Defaults vary between scripts — `end_to_end.py` uses `(0.5, 1.0)`,
   `process_music.py` uses `(0.33, 0.33)`.  Tune for your speaker.

### What (h2, h3) do

- **h2 only**: clean octave‑up harmonic.  Smooth, warm bass that doesn't sound
  artificial.
- **h3 only**: octave‑and‑a‑fifth up.  Adds edge and punch but can sound buzzy
  on sustained notes.
- **Both**: the default.  `h2` dominates at higher bass frequencies (where LP(fc)
  passes more energy), `h3` dominates at sub‑bass (where LP(fc/2) is the only
  signal reaching the T₃ path).
- **Higher cutoff**: moves the effect up.  At 120 Hz you're enhancing low‑mids
  rather than sub‑bass — useful for tiny phone speakers.

---

## Build

### Desktop

```bash
# Python environment (one of)
nix develop                              # if you use Nix
pip install numpy scipy matplotlib       # otherwise

# C DSP shared library
make -C src                              # builds enhancer.so + eqgen_ladspa.so
```

### ESP32 (firmware)

1. Run `python -m eqgen.cli.export ... -o src/eq_coeffs.h` to generate coefficients.
2. Add these DSP source files to your ESP‑IDF project:
   - `src/enhancer.c`, `src/enhancer.h`
   - `src/biquad.h`
   - `src/lp.h`
   - `src/envelope.h`
   - `src/eq_coeffs.h` (generated)
3. In your audio callback (A2DP sink → I2S), convert int16 samples to float,
   call `dsp_pipe_process_stereo()`, and convert back.  The DSP runs at
   44.1 kHz or 48 kHz — whatever sample rate the coefficients were designed at.

   With 24 EQ biquads the DSP consumes ~20% of one core on a 240 MHz ESP32,
   leaving the second core free for the Bluetooth stack and I2S DMA.

---

## The psychoacoustic model (`python/model.py`)

The desktop pipeline uses a frequency‑domain model to predict what the enhancer
does *before* running audio through it.  For a sine wave at frequency `f` with
linear gain `G`:

```
A(f, G) = G · √[ HP(f,fc)² · S(f)²
                 + h2² · LP(f,fc)²  · HP(2f,fc)² · S(2f)²
                 + h3² · LP(f,fc/2)² · HP(3f,fc)² · S(3f)² ]
```

Where:
- `S(f)` is the speaker's frequency response
- `HP` / `LP` are 2nd‑order Butterworth high‑pass / low‑pass at the given cutoff
- Harmonic terms are weighted by 5× to match ISO 226 equal‑loudness contours
  (the ear is ~5× more sensitive at 100 Hz than 50 Hz)

The EQ preprocess step solves `A(f, G) = target(f)` for `G` at each frequency —
so the EQ + enhancer combined produce the desired perceived level.

---

## Python modules

| Module | Purpose |
|---|---|
| `model.py` | Psychoacoustic model: computes perceived output amplitude and the EQ gain needed to hit a target. |
| `eq_fit.py` | Greedy IIR biquad fitter with golden‑section gain/Q search. Cascades peaking filters at the frequencies with largest residual error. |
| `quantize.py` | Q4.28 fixed‑point quantization, ReciprocalLUT (1/x lookup for envelope normalisation). |
| `dsp.py` | Python simulation of all 4 enhancer variants (original, linear, fixed, fixed_v2). Used for debugging and offline analysis. |
| `analysis.py` | Goertzel single‑frequency DFT, FFT peak finding, harmonic description, test signal generator. |
| `enhancer_ffi.py` | ctypes wrapper for `src/enhancer.so` — create/destroy/process instances from Python. |
| `audio_io.py` | WAV reading (via scipy) and measurement loading. |

---

## Measurements directory

```
measurements/
├── technics/standing/    # Technics speakers, standing mic position
├── anker/                # Anker Bluetooth speaker
├── cardboard/            # Cardboard-box speaker prototype
├── jamo/                 # Jamo bookshelf speakers
├── kenrad/               # Vintage Kenrad radio
├── phone/                # Phone speaker
├── laptop/               # Laptop speakers
├── corolla/              # Toyota Corolla car stereo
├── edifier-wall/         # Edifier speakers near wall
├── opus1/                # Opus 1 DIY speaker
└── ALDI/                 # ALDI budget speaker
```

Each contains `.wav` recordings (for the current pipeline) and possibly legacy
`.txt` FFT dumps (recorded with an older measurement tool).

---

## Hardware

The standalone Bluetooth DSP box runs on:

- **ESP32-WROOM-32** dev board (38‑pin generic) — the original ESP32, not C3/S2/S3.
  The C3/S3 lines lack Bluetooth Classic hardware and cannot run A2DP.
- **PCM5102A I2S DAC module** — 112 dB SNR, line‑out via 3.5 mm jack, no MCLK needed.
- **3 wires**: BCK, LRCK (WS), DIN between ESP32 and DAC.  No codec, no driver IC.

Total BOM: ~$8–12.

| Component | Search keywords | Notes |
|---|---|---|
| ESP32 dev board | "ESP32 WROOM-32 38-pin" | Must be original ESP32, not C3/S2/S3 |
| I2S DAC module | "PCM5102 I2S DAC" | Line‑out on 3.5 mm jack, ~$2 |
| Power | USB‑C breakout or 5V barrel jack | Dev board already has 3.3 V regulator |
| Enclosure | "aluminum project box" + panel‑mount 3.5 mm jack | Optional |

Avoid: ESP32-A1S, ESP32-LyraT, or any "ESP32 audio board" with an onboard codec —
these are harder to hack than bare I2S.

## License

[to be determined]
