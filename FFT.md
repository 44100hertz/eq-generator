# FFT + IIR Hybrid EQ

Hybrid architecture where an FFT overlap-add EQ provides broad correction across
the full band, and IIR biquads surgically correct the narrow peaks and dips the
FFT can't resolve at its coarse bin resolution. The bass enhancer runs last,
generating harmonics below the cutoff frequency.

FFT resolution at 256 points is 187.5 Hz @ 48 kHz — too coarse for narrow bass
resonances, but perfectly adequate for the broad shape of most speaker corrections.

## Architecture

```
  ┌──────────┐     ┌───────────────┐     ┌──────────────┐
  │   FFT    │────▶│  IIR biquads  │────▶│ Bass Enhancer│────▶ output
  │ (broad)  │     │ (surgical)    │     │ (harmonics)  │
  └──────────┘     └───────────────┘     └──────────────┘
```

The FFT applies per-bin scalar gains at bin center frequencies (sampled from the
target correction curve). The IIR biquads are fitted to the *residual* — what the
FFT missed. Both paths are linear and multiplicative (cascaded).

No crossfade or crossover filter is needed. Two sample rates are supported (44100
and 48000 Hz) with rate-specific coefficient arrays generated at design time.

## FFT EQ runtime

| Parameter | Value | Notes |
|-----------|-------|-------|
| Window size | 256 samples | 5.3 ms @ 48 kHz |
| Hop size | 128 samples | 50% overlap, 2.7 ms |
| Window | Hann | Standard overlap-add for perfect reconstruction at unity gain |
| Frequency resolution | 187.5 Hz | fs / 256 |
| Latency | ~5 ms | Frame + overlap-add overhead |

**Implementation**: `src/fft_eq.c` — custom radix-2 DIT complex FFT (in-place,
interleaved real/imag). Uses a straightforward complex FFT for simplicity;
a real FFT would save ~2×, but this gives a conservative performance bound.

## Processing order (per 128-sample frame)

```
  1. Shift overlap buffer: samples [128..255] → [0..127]
  2. Copy 128 new samples into overlap[128..255]
  3. Apply Hann window → FFT input buffer
  4. Forward FFT (256-point complex)
  5. Multiply each complex bin by its scalar gain from table (preserves
     conjugate symmetry)
  6. Inverse FFT → time domain
  7. Overlap-add: FFT output[0..127] + saved overlap from previous frame
  8. Save FFT output[128..255] as overlap for next frame
  9. Output 128 samples → IIR biquad cascade → bass enhancer
```

**State (per channel)**: overlap buffer [256 floats], overlap-add accumulator
[128 floats]. Workspace [512 floats] is reused between left and right channels.
Total RAM: ~4.5 KB for stereo.

Per-bin gains are clamped to ±12 dB (0.25× to 4.0×) for safety.

## Python pipeline

`compute_fft_residual()` in `eqgen/pipeline.py`:

1. Sample the target correction curve at FFT bin center frequencies
2. Clamp to safe gain range (±12 dB)
3. Cubic-spline interpolate the piecewise-constant bin gains back to the
   evaluation grid, approximating the Hann window's ~3-bin blending
4. Subtract from the target to get the residual
5. Threshold residual at ±1 dB — below this, FFT is "close enough" and IIR
   biquads aren't wasted on noise

```python
from eqgen.pipeline import compute_fft_residual, FFT_N

fft_bin_freqs, fft_gains_linear, fft_approx_db, residual_db = \
    compute_fft_residual(freqs, target_db, fs, FFT_N)
```

The residual is self-normalizing — it oscillates around 0 dB because the FFT
handles the broad shape. IIR biquads are fitted to the residual via
`fit_eq_curve()`, which handles both positive and negative gain at every band.

## Design-time flow

`eqgen.cli.wire.run_full_pipeline()` and `eqgen.cli.eqgen.design_eq()`:

1. Run the full measurement pipeline (Welch FFT → CV-weighted smoothing →
   correction → pre-gain → bass enhancer preprocessing)
2. Compute FFT residual at both 44100 and 48000 Hz
3. Fit IIR biquads to the residual at each rate via greedy peaking filter fitter
4. Export both FFT gain tables and IIR coefficients to `src/eq_coeffs.h`:
   - `EQGEN_FFT_BINS` (129 bins, DC through Nyquist)
   - `eqgen_fft_gains_44100[129]` and `eqgen_fft_gains_48000[129]` float tables
   - `EQGEN_N_BIQUADS` + Q4.28 biquad arrays at both sample rates
   - Runtime selector: `eqgen_get_fft_gains(rate)` and `eqgen_get_coeffs(rate)`

When `EQGEN_FFT_BINS == 0`, the FFT path is bypassed at runtime:
```c
#if EQGEN_FFT_BINS > 0
    FftEq *fft_eq = fft_eq_create(fft_gains);
#endif
```

## C runtime integration

`src/filter.c` (the PipeWire filter chain entry point):

```
  deinterleave 128 stereo samples
  → fft_eq_process_frame()          [FFT EQ, broad correction]
  → IIR biquad cascade (Q28)        [surgical residual fit]
  → BassEnhancer_process_stereo()   [harmonics below cutoff]
  → interleave → stdout
```

## Performance

Benchmark in `src/bench_hybrid.c`. On desktop hardware the combined FFT + 6-biquad
chain comfortably exceeds real-time by >100×. Conservative ESP32 estimate (15×
desktop slowdown) puts CPU usage in the low single-digit percent range for the
FFT portion. The custom radix-2 DIT FFT avoids external dependencies (no KissFFT,
no esp-dsp) at a modest compute cost — a 256-point real FFT would be ~2× faster
but the complex FFT is simpler and sufficient.

## Web visualization

The web server (`eqgen/server.py`) generates these traces per analysis run:

| Trace | Description |
|-------|-------------|
| `correction` | Full target correction curve (dB) |
| `fft_response` | FFT per-bin gain approximation (dB, at bin centers) |
| `iir_fit` | IIR biquad response fitted to the residual (dB, hi-res grid) |
| `combined` | FFT + IIR combined response (dB, hi-res grid) |
| `prediction_error` | Combined − target error (dB) |

Rendered client-side via Canvas with no additional JS dependencies.
