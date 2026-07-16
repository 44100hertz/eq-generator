# eqgen — Speaker EQ Correction with Psychoacoustic Bass Enhancement

A rapid prototyping suite for designing IIR EQ curves with harmonic
bass enhancement.  Measure a speaker on the desktop, tune parameters, audition
instantly, then **export float coefficients** for deployment on an ESP32
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
                                                    ╱          ╲
                                                   ▼            ▼
                                              JSON curve    enhancer.so
                                              (EQ points)   (ctypes FFI)
                                                                 │
                                                  ┌──────────────┴──────────────┐
                                                  ▼                             ▼
                                          eqgen.cli.enhance              eqgen.cli.wire
                                          (offline WAV)              (live PipeWire LADSPA)
```

### DSP core (runs identically on desktop and ESP32)

```
  input (float)
     │
     ▼
  DC blocker
     │
     ▼
  pre‑gain
     │
     ▼
  cascaded EQ biquads (float)
     │
     ▼
  LR4 crossover → LP(fc) ──→ envelope ──→ Chebyshev T₂/T₃ ──→ mix
     │                                                           │
     └──→ HP(fc) dry path ───────────────────────────────────────┤
                                                                  ▼
                                                            loudness shelf
                                                                  │
                                                                  ▼
                                                            full‑band limiter
                                                                  │
                                                                  ▼
                                                            output (float)

     Harmonics generated via Chebyshev polynomials from the LP path.
     Crossfade between dry HP and harmonics controlled by headroom budget.
     Full‑band peak limiter with 3 s release catches overboost transients.
```

---

## Directory overview

```
eqgen/
├── eqgen/                  # Python package
│   ├── cli/                # CLI entry points
│   │   ├── eqgen.py        # Design EQ curve → JSON
│   │   ├── audition.py     # Full pipeline: speaker → IIR → C DSP → WAV
│   │   ├── enhance.py      # Standalone C enhancer (flat EQ, no measurements)
│   │   ├── wire.py         # Live PipeWire system‑wide EQ
│   │   ├── export.py       # Export fitted biquads as C header
│   │   └── preset.py       # Preset management CLI
│   ├── pipeline.py         # Shared pipeline: WAV → EQ curve
│   ├── eq_fit.py           # Greedy IIR biquad fitter
│   ├── model.py            # Psychoacoustic model
│   ├── dsp.py              # Filter design utilities
│   ├── analysis.py         # Goertzel, FFT analysis, test‑signal generation
│   ├── sweep.py            # Sweep‑based analysis through C enhancer
│   ├── process.py          # Audio file processing through C enhancer
│   ├── smart_volume.py     # Smart volume curve modeling
│   ├── enhancer_ffi.py     # ctypes wrapper for src/enhancer.so
│   ├── io.py               # WAV I/O and measurement loading
│   ├── iso226.py           # ISO 226:2023 equal‑loudness contours
│   ├── presets.py          # Preset system (JSON‑based)
│   ├── server.py           # Web UI server
│   ├── web_ui/             # Web UI assets
│   ├── analysis_tools/     # Loudness analysis utilities
│   └── tests/              # Test suite
│
├── src/                    # C DSP — same code ships on ESP32
│   ├── enhancer.c / .h     # Harmonic bass enhancer (float DSP)
│   ├── enhancer_api.c / .h # C API for Python ctypes
│   ├── biquad.h            # Float biquad filter
│   ├── envelope.h          # Envelope follower
│   ├── lp.h                # Low‑pass filter
│   ├── smart_volume.h      # Smart volume math (shared C header)
│   ├── volume_control.h    # Volume orchestration wrapper
│   ├── filter.c            # Native PipeWire filter process
│   ├── sfx_data.h          # Sound effects data
│   ├── fpu.h               # ESP32 FPU workaround
│   ├── eq_coeffs.h         # Generated EQ coefficients (autogenerated)
│   ├── ladspa/             # LADSPA plugin wrapper
│   ├── test_enhancer.c     # Standalone C test
│   └── Makefile
│
├── firmware/               # ESP32 firmware (ESP‑IDF)
│   └── main/               # bt_a2dp, i2s_out, main, sfx_player
│
├── measurements/           # Per‑speaker recording directories
│   └── <speaker_name>/
│       ├── measurement*.wav # Brown‑noise recording of the speaker
│       ├── target.wav       # Brown‑noise recording of reference system
│       └── noise.wav        # Ambient noise floor (optional, for subtraction)
│
├── presets/                # Preset JSON files + house_curves.json
├── output/                 # Build artifacts + reports
├── prototypes/             # Historical / tuning scripts
├── docs/                   # Technical documentation
└── flake.nix               # Nix dev shell
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
python -m eqgen.cli.eqgen \
  -m measurements/technics/standing/measurement2.wav \
  -t measurements/technics/standing/target.wav \
  -n measurements/technics/standing/noise2.wav \
  --bass-enhancer-cutoff 50 \
  -o my_eq.json
```

Produces `my_eq.json` with the EQ curve (freq/gain_db points, sample rate, max gain).

Without bass enhancement:

```bash
python -m eqgen.cli.eqgen -m meas.wav -t target.wav -o flat_eq.json
```

### 2. Audition offline

Process real music through the full pipeline (measurement → IIR fit → C DSP → WAV):

```bash
python -m eqgen.cli.audition technics/standing /tmp/out --fc 60 --tracks song1.flac song2.mp3
```

Or use the C enhancer standalone (flat EQ, no measurement):

```bash
python -m eqgen.cli.enhance ~/Music /tmp/out
```

### 3. Audition live (Linux + PipeWire)

Wire the enhancer into your system audio path:

```bash
python -m eqgen.cli.wire setup technics/standing
```

All system audio now flows through:  app → null sink → native DSP filter → hardware output.

Tear down:

```bash
python -m eqgen.cli.wire teardown
```

### 4. Iterate

Tweak cutoff, harmonic amplitudes, max bands — re‑run, re‑listen until it sounds right.

### 5. Export for ESP32

When satisfied, bake the coefficients into a C header:

```bash
python -m eqgen.cli.export --speaker small --fc 60 --h2 0.33 --h3 0.33 \
    --max-bands 6 -o src/eq_coeffs.h
make -C src
```

The DSP code in `src/enhancer.c` + `src/biquad.h` + `src/envelope.h` +
`src/smart_volume.h` + `src/volume_control.h` compiles directly under
ESP‑IDF — no changes needed.

---

## Entry point reference

| Command | Purpose |
|---|---|
| `python -m eqgen.cli.eqgen` | Measurements → EQ curve JSON.  Core tuning tool. |
| `python -m eqgen.cli.audition` | Speaker → IIR fit → C DSP → WAV.  Full offline audition. |
| `python -m eqgen.cli.enhance` | Run the C enhancer on audio files (flat EQ, no measurements). |
| `python -m eqgen.cli.wire setup/teardown` | Live PipeWire LADSPA system‑wide EQ.  Real‑time audition. |
| `python -m eqgen.cli.export` | Export fitted biquads as `src/eq_coeffs.h` for ESP32 firmware. |
| `python -m eqgen.cli.preset` | Preset management CLI (list, show, guess, etc.). |
| `python -m eqgen.server` | Web UI server for preset management and visualization. |
| `make test` | Run all Python tests (`python -m eqgen.tests.run_all`). |

---

## Parameter guide

| Parameter | CLI flag | Default | Range | Description |
|---|---|---|---|---|
| **Cutoff** | `--bass-enhancer-cutoff` / `--fc` | 60 Hz | 20–120 Hz | Crossover frequency. Below this, harmonics are synthesized instead of the fundamental. |
| **h2 amp** | `--h2` | 0.5† | 0.0–1.0 | 2nd‑harmonic amplitude. Controls perceived bass warmth. |
| **h3 amp** | `--h3` | 0.33† | 0.0–1.0 | 3rd‑harmonic amplitude. Controls aggressive bass edge. |
| **Max bands** | `--max-bands` | 24 | 1–40 | Maximum IIR peaking biquads. Fewer = simpler, more = tighter fit. |
| **Smooth exponent** | `--smooth-exponent` | 1.0 | 0–2 | CV smoothing aggressiveness. 0=fixed bandwidth, 1=linear, 2=amplified. |
| **High rolloff** | `--high-rolloff FREQ,DB` | none | e.g. `4000,-6` | Shelf rolloff on the target curve above FREQ. Repeatable. |
| **Low rolloff** | `--low-rolloff FREQ,DB` | none | e.g. `100,-12` | Shelf rolloff on the target curve below FREQ. Repeatable. |
| **Release** | `--release` | 0.2 s | 0.01–2.0 s | Envelope follower release time for harmonic linearisation. |
| **Overboost** | `--overboost` | 0 dB | 0–6 dB | Extra gain at max volume to drive harmonics (used in export). |
| **Speaker level** | `--speaker-level` | 60 dB | 24–80 dB | Total system gain (speaker + amp) for smart volume LUT. |

†  Defaults vary between scripts — e.g. `eqgen.cli.audition` and `eqgen.cli.enhance`
   use `(0.5, 1.0)`.  Tune for your speaker.

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

1. Run `python -m eqgen.cli.export --speaker <name> -o src/eq_coeffs.h` to generate
   coefficients.
2. Add these DSP source files to your ESP‑IDF project:
   - `src/enhancer.c`, `src/enhancer.h`
   - `src/biquad.h`
   - `src/lp.h`
   - `src/envelope.h`
   - `src/smart_volume.h`
   - `src/volume_control.h`
   - `src/eq_coeffs.h` (generated)
3. In your audio callback (A2DP sink → I2S), convert int16 samples to float,
   call `dsp_pipe_process_stereo()`, and convert back.  The DSP runs at
   44.1 kHz or 48 kHz — whatever sample rate the coefficients were designed at.

   With 24 EQ biquads the DSP consumes ~20% of one core on a 240 MHz ESP32,
   leaving the second core free for the Bluetooth stack and I2S DMA.

### Makefile targets

| Target | Description |
|---|---|
| `make` / `make all` | Build `enhancer.so` (Python FFI) + `eqgen_ladspa.so` |
| `make test` | Run all Python tests |
| `make clean` | Remove build artifacts |
| `make eqgen ARGS="..."` | Design EQ curve → JSON |
| `make audition ARGS="..."` | Full pipeline audition |
| `make wire-setup ARGS="..."` | Live PipeWire system‑wide EQ |
| `make wire-teardown` | Remove PipeWire wiring |
| `make export ARGS="..."` | Export coefficients for ESP32 |
| `make server` | Start web UI server |
| `make flash ARGS=<speaker>` | Export + build + flash ESP32 firmware |

---

## Smart volume

The firmware implements a smart volume system that adjusts the loudness contour
and pre‑gain at every volume step (0–127, controlled via AVRCP from the paired
phone).  At low volumes, a loudness shelf compensates for the ear's reduced
bass sensitivity (Fletcher‑Munson / ISO 226), and pre‑gain is reduced to
maintain headroom.

The core math lives in `src/smart_volume.h` and `src/volume_control.h` —
shared C headers used identically on desktop and ESP32.  The Python simulation
is in `eqgen/smart_volume.py`.

See `SMART_VOLUME.md` for details.

---

## The psychoacoustic model (`eqgen/model.py`)

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
- Harmonic terms are weighted to match ISO 226 equal‑loudness contours
  (the ear is ~5× more sensitive at 100 Hz than 50 Hz)

The EQ preprocess step solves `A(f, G) = target(f)` for `G` at each frequency —
so the EQ + enhancer combined produce the desired perceived level.

---

## Python modules

| Module | Purpose |
|---|---|
| `pipeline.py` | Shared pipeline: Welch FFT → CV smoothing → correction → efficacy.  Single source of truth. |
| `eq_fit.py` | Greedy IIR biquad fitter with golden‑section gain/Q search. Cascades peaking filters at the frequencies with largest residual error. |
| `model.py` | Psychoacoustic model: computes perceived output amplitude and the EQ gain needed to hit a target. |
| `dsp.py` | Filter design utilities: biquad response, pre‑gain computation, default EQ coeffs. |
| `analysis.py` | Goertzel single‑frequency DFT, Welch FFT stats, test signal generation. |
| `enhancer_ffi.py` | ctypes wrapper for `src/enhancer.so` — create/destroy/process instances from Python. |
| `io.py` | WAV reading (via scipy) and measurement loading. |
| `iso226.py` | ISO 226:2023 equal‑loudness contours (loaded from digitized CSV table). |
| `sweep.py` | Sweep‑based analysis through the real C enhancer pipeline. |
| `process.py` | Audio file processing (ffmpeg → C enhancer → WAV). |
| `smart_volume.py` | Smart volume curve modeling — matches `smart_volume.h` exactly. |
| `presets.py` | Preset system: JSON‑based preset files, house curves, measurement directory scanning. |
| `server.py` | Web UI server for preset management and pipeline visualization. |

---

## Measurements directory

```
measurements/
├── technics/standing/    # Technics speakers, standing mic position
├── cardboard/            # Cardboard‑box speaker prototype
├── lunchbox/             # Lunchbox speaker build
├── lunchbox-noboost/     # Lunchbox without overboost
└── ...                   # (add your own speakers here)
```

Each contains `.wav` recordings (for the current pipeline).

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
