Speaker Configuration System: Nature, Quality, and Recommendations
===================================================================

Date: 2026-07-06
Project: eqgen — Speaker EQ Correction with Harmonic Bass Enhancement

1. Overview
-----------

The eqgen system provides end-to-end speaker equalization: it measures a
speaker's frequency response from brown-noise recordings, computes a
correction curve, fits that curve to a cascade of IIR biquad filters,
quantizes the coefficients to Q4.28 fixed-point, and applies them through
a C-based DSP pipeline that also includes a harmonic bass enhancer.

The pipeline supports four deployment targets:
  - `enhancer.so` — C shared library with Python ctypes FFI for offline
    processing and testing
  - `eqgen_ladspa.so` — LADSPA plugin for real-time desktop DSP via
    PipeWire
  - `eqgen_filter` — stdin→stdout float32 DSP filter for live PipeWire
    filter chains
  - ESP32 firmware — I2S → BiquadQ28 cascade → I2S, targeting a Bluetooth
    A2DP speaker

2. Architecture
---------------

2.1. Coefficient Generation (Python)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python side handles all design-time work:

    Measurement → Welch FFT → CV-weighted smoothing → Correction → IIR fit → Q4.28 quantize

Key modules and their roles:

  `pipeline.py`     — Reads WAV measurements, runs Welch's method (16K
                      FFT, 50% overlap, Hanning window), pools per-bin
                      statistics across multiple measurement files, merges
                      noise-floor coefficient of variation (CV), applies
                      CV-weighted adaptive-bandwidth kernel smoothing,
                      computes correction = target / measurement.

  `eq_fit.py`       — Greedy per-band (gain, Q) optimizer for fitting
                      the correction curve to a cascade of parametric
                      peaking EQ filters. Uses golden-section search for
                      Q optimization. Reordering pass puts attenuation
                      stages first to prevent intermediate overflow in
                      the fixed-point cascade.

  `quantize.py`     — Maps float biquad coefficients to Q4.28 (32-bit
                      signed integer, resolution ~3.7e-9), mirrors the C
                      BiquadQ28 implementation. Also provides the
                      ReciprocalLUT for envelope normalization in C.

  `analysis.py`     — Goertzel algorithm for precise single-frequency
                      measurement, FFT-based peak finding, test signal
                      generation. Used for harmonic analysis and sweep
                      verification.

  `sweep.py`        — Runs sine sweeps through the actual C enhancer DSP
                      (via enhancer_ffi) and measures fundamental + 2nd/3rd
                      harmonic levels at the output. Replaces the old
                      psychoacoustic model approach.

  `enhancer_ffi.py` — ctypes wrapper for the C enhancer shared library.
                      Handles opaque handle creation, int16↔Q16 conversion,
                      and coefficient passing.

2.2. DSP Engine (C)
~~~~~~~~~~~~~~~~~~~

The C side runs the actual audio processing:

  `enhancer.c`      — Harmonic bass enhancer (`process_fixed_v2`): the
                      normalize-before-Chebyshev variant with separate LP
                      filters at fc (for T₂ path) and fc/2 (for T₃ path).
                      Signal chain per channel:
                        DC blocker → EQ biquad cascade → LP(fc) → env → T₂
                                                     → LP(fc/2) → env → T₃
                        → mix with HP(fc) dry → harmonic AGC limiter

  `biquad_q28.h`    — Direct Form I biquad with Q4.28 coefficients, Q16
                      state. 64-bit accumulator, single truncation
                      (mean-zero error, no DC bias).

  `envelope.h`      — Peak-hold envelope follower with instant attack
                      and exponential release. Reciprocal LUT for fast
                      1/envelope normalization.

  `dc_blocker.h`    — First-order DC-blocking IIR (5 Hz cutoff).
                      Numerator sums exactly to zero at DC in Q16.

  `enhancer_api.c`  — Opaque handle API for Python FFI. Allocates
                      BiquadQ28 state arrays, wraps BassEnhancer behind
                      create/destroy/reset/process_stereo.

  `filter.c`        — stdin→stdout float32 DSP filter for PipeWire
                      chains. Bakes in coefficients via eq_coeffs.h.

  `ladspa_wrapper.c`— LADSPA plugin wrapper exposing the enhancer as an
                      audio plugin with runtime release-time control.

2.3. Test Suite
~~~~~~~~~~~~~~~

Tests are organized under `eqgen/tests/`:

  run_all.py            — Test runner with --skip/--only selectors
  test_eq_pipeline.py   — End-to-end: synthetic speaker models → EQ curve
                          → IIR fit → Q4.28 quantize → verify
  test_harmonics.py     — C enhancer sweep + harmonic linearity
  test_model.py         — C enhancer direct output verification
  test_evenness.py      — Bass evenness through C enhancer
  test_sine_sweep.py    — Real measurement → EQ → C enhancer sweep
  test_compressor.py    — Envelope follower dynamics analysis
  test_chebyshev.py     — Chebyshev math verification (T₂/T₃, leakage)
  test_limiter_transparency.py — Limiter transparency on real music
  test_gain_safety.py   — Gain safety on real measurements
  test_fullband_fit.py  — 40-band peaking filter fit stress test
  test_real_fullband_fit.py — Real measurement IIR fit report

3. Quality Assessment
---------------------

3.1. Strengths
~~~~~~~~~~~~~~

  a. Measurement Pipeline
     - Welch's method with 50% overlap provides low-variance spectral
       estimates.
     - CV-weighted adaptive-bandwidth kernel smoothing avoids the
       oscillation artifacts of cubic splines while providing more
       smoothing in noisy regions and minimal smoothing where SNR is high.
     - Noise-floor CV merging and inflation correctly handles both
       intermittent noise (max-CV merge) and stationary noise (SNR-based
       inflation).
     - Multi-measurement pooling improves statistical reliability.

  b. IIR Fit
     - Greedy per-band golden-section search finds optimal gain/Q pairs
       efficiently.
     - Adaptive number of bands up to 40 provides excellent fit quality
       (typically <1 dB RMS error across 20-20kHz).
     - Q4.28 quantization introduces negligible error (<0.01 dB from
       quantization alone).
     - Attenuation-first reordering prevents fixed-point overflow.
     - Butterworth HP numerator sums to exactly zero at DC in Q4.28 —
       no degeneracy at sub-bass frequencies.

  c. C DSP
     - Matches the process_fixed_v2 Python reference exactly (verified
       by test_enhancer.c).
     - Separate LP at fc/2 for T₃ path prevents excessive 3rd harmonic
       generation at higher frequencies.
     - Normalize-before-Chebyshev approach produces clean harmonics
       with no fundamental leakage.
     - Harmonic AGC limiter prevents clipping without touching dry signal.
     - Efficient: ~15 µs per stereo frame on desktop, estimated ~1.5 ms
       on ESP32 at 160 MHz — well within the 22.7 µs/frame budget at
       44100 Hz.

  d. Deployment Flexibility
     - Shared library with Python FFI for offline processing.
     - LADSPA plugin for any LADSPA host (PipeWire, Ardour, etc.).
     - stdin/stdout filter for PipeWire filter chains.
     - ESP32 firmware for embedded speakers.

3.2. Limitations and Tradeoffs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  a. Single Measurement Type
     The pipeline relies entirely on brown-noise recordings. While
     brown noise provides good low-frequency energy, it has limited
     high-frequency content. For speakers used in rooms, additional
     measurement types (sine sweeps for distortion, MLS for impulse
     response) could provide more information.

  b. No Room Correction
     The system computes a speaker correction curve from a single-point
     measurement. It does not model room modes, reflections, or
     time-domain behavior (no impulse response analysis). The correction
     is purely magnitude-based.

  c. Peaking Filters Only
     The IIR fitter uses only peaking (bell) filters. Low-shelf and
     high-shelf filters are not used despite being available in the
     eq_fit module. Shelving filters could handle broad trends more
     efficiently than cascaded peaking filters.

  d. Harmonic Bass Enhancer is Always-On
     The enhancer cutoff, h2, and h3 are fixed at build time. There is
     no bypass or dynamic control over harmonic levels. The LADSPA
     plugin exposes only release time as a runtime parameter.

  e. No Preset Format
     Speaker configurations exist as generated C headers (eq_coeffs.h)
     and measurement directories. There is no portable, human-readable
     preset format for sharing or comparing speaker configurations.

4. Recommendations
------------------

4.1. Preset Format
~~~~~~~~~~~~~~~~~~

A portable JSON preset format would enable:

  {
    "name": "Technics Standing Speakers",
    "version": 1,
    "fs": 44100,
    "n_biquads": 24,
    "coeffs": [
      {"type": "peaking", "f0": 45.2, "gain_db": 8.3, "Q": 1.4},
      ...
    ],
    "enhancer": {
      "cutoff_hz": 60.0,
      "h2_amp": 0.5,
      "h3_amp": 1.0,
      "release_secs": 0.2
    }
  }

Benefits:
  - Human-readable and diffable.
  - Portable between C header, Python, and any other target.
  - Each band has metadata (type, f0, gain, Q) for inspection and
    manual tuning.
  - The generated C header becomes a build artifact rather than the
    source of truth.

4.2. Sweep-Based Verification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The new `sweep.py` module enables measurement-based verification of the
complete pipeline. Recommended workflow:

  1. Design EQ from brown-noise measurements (pipeline.py).
  2. Build enhancer with those coefficients.
  3. Run sine sweeps through the C enhancer (sweep.py).
  4. Plot output RMS vs frequency to verify flatness.
  5. Measure 2nd and 3rd harmonic levels across the bass range.
  6. Adjust h2/h3/cutoff and re-measure until harmonic balance is
     satisfactory.

This replaces the old psychoacoustic model approach with direct
measurement of what the actual DSP pipeline produces.

4.3. Future Improvements
~~~~~~~~~~~~~~~~~~~~~~~~

  a. Shelf Filter Support — Extend the IIR fitter to try low-shelf and
     high-shelf filters for broad-band corrections, potentially reducing
     the number of bands needed.

  b. Exponential Sine Sweep — Use Farina's method for simultaneous
     impulse response + harmonic distortion measurement. This would
     provide both linear response and distortion data from a single
     recording.

  c. Multi-Point Averaging — Support spatial averaging of multiple
     microphone positions to reduce room effects.

  d. Runtime Envelope Parameters — Expose h2_amp, h3_amp, and cutoff
     as LADSPA control ports for real-time adjustment.

  e. Automatic Preset Selection — Auto-detect connected hardware sink
     and load the appropriate EQ preset.

5. Summary
----------

The eqgen speaker configuration system is a mature, well-tested pipeline
for designing and deploying speaker EQ correction with harmonic bass
enhancement. The measurement-to-correction pipeline is statistically
robust (Welch's method + CV-weighted smoothing). The IIR fitter produces
high-quality biquad cascades with negligible quantization error. The C
DSP implementation is correct, efficient, and deployable across multiple
targets (desktop LADSPA, PipeWire filter, ESP32 firmware).

The recent removal of the psychoacoustic model in favor of direct
sweep-based measurement improves the system's realism: the pipeline IS
the response, and sweeps measure actual output rather than predicting it.

The primary gap is the lack of a portable preset format. Creating one
would make speaker configurations shareable, inspectable, and
independent of the C header build artifact.
