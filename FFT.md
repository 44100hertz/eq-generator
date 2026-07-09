# FFT + IIR Hybrid EQ

Hybrid architecture where an FFT provides unlimited treble bands, and IIR biquads
handle the narrow bass bands where FFT resolution is poor.

## Architecture

```
  ┌──────────┐     ┌──────────┐     ┌──────────┐
  │   IIR    │────▶│   FFT    │────▶│ Overlap  │────▶ I2S
  │ (bass)   │     │ (treble) │     │  Add     │
  └──────────┘     └──────────┘     └──────────┘
```

No crossover filter needed. IIR covers the full band; FFT applies the residual —
whatever the IIR missed above the crossover frequency. Both paths are linear and
additive.

## FFT parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Window size | 256 samples | 5.3 ms @ 48 kHz |
| Hop size | 128 samples | 50% overlap, 2.7 ms |
| Window | Hann | Standard overlap-add |
| Frequency resolution | 187.5 Hz | 48k / 256 |
| Latency | ~5 ms | Frame + overlap-add overhead |

## CPU budget

~12M MAC/s (256-point real FFT + per-bin gain multiply + inverse FFT + window +
overlap-add). Comparable to ~25 biquads on ESP32 at 240 MHz.

Two options for the FFT library:

1. **esp-dsp** (ESP-IDF component) — `dsps_fft2r_fc32`, radix-2 fixed-point,
   ~1500 cycles for 256-point. Path of least resistance if ESP-IDF is the target.
2. **KissFFT** (single `.c` file, no deps) — ~6000 cycles for 256-point float.
   Portable, works on desktop too.

## Generator changes (~30 lines Python)

After the IIR fit, split the frequency range:

```python
# Freqs above crossover (e.g. 300 Hz)
fcross = 300.0
mask = fft_bin_freqs >= fcross

# IIR response at FFT bin centers
iir_db = cascade_response_db(iir_bq, fft_bin_freqs, design_fs)

# Residual error (FFT must correct what IIR missed)
residual_db = ideal_correction_db - iir_db
residual_db[~mask] = 0.0  # FFT only handles treble

# Per-bin scalar gains (linear, not dB)
fft_gains = 10.0 ** (residual_db / 20.0)
fft_gains = np.clip(fft_gains, 0.0, 4.0)  # safety clamp

# Export both to eq_coeffs.h:
#   EQGEN_N_BIQUADS + biquad arrays (unchanged)
#   EQGEN_FFT_BINS (128 real bins)
#   eqgen_fft_gains[128] float32 table
```

## Runtime changes (~120 lines C)

New file `src/fft_eq.c`, called from `audio_data_handler` after the IIR pass:

```
Input: 128 new stereo samples (from IIR output)
State: float overlap_buf[256] per channel  (2 × 1 KB)
       float fft_in[256], fft_out[256]      (2 × 2 KB reusable workspace)

Per 128-sample tick:
  1. Shift overlap buffer: copy samples [128..255] → [0..127]
  2. Copy 128 new samples into overlap[128..255]
  3. Apply Hann window to overlap[0..255] → fft_in
  4. Real FFT: fft_in[256] → fft_out[128 complex bins]
  5. Multiply each complex bin by its scalar gain from table
  6. Inverse FFT: fft_out → fft_in (now time-domain)
  7. Overlap-add: fft_in[0..127] += overlap[0..127]
  8. Overlap-add: fft_in[128..255] → overlap[128..255] (save for next tick)
  9. Output fft_in[0..127] to I2S
```

RAM: ~6 KB (two overlap buffers + two FFT workspaces, reusable per tick).

## Web app visualization (~2 extra traces)

The pipeline result dict gets two new traces:

```python
result["fft_residual"] = [{"freq": f, "db": d} for f, d in ...]
result["combined"]     = [{"freq": f, "db": d} for f, d in ...]
```

Chart renderer handles them with zero new JS — just pass them to `drawChart()`:

```javascript
drawChart('cIir', [
  {name:'Ideal',      data: r.correction,     color: COLORS[1]},
  {name:'IIR fit',    data: r.iir_fit,        color: COLORS[4], dash: [6,3]},
  {name:'+ FFT',      data: r.fft_residual,   color: COLORS[2], dash: [2,2]},
  {name:'Combined',   data: r.combined,       color: COLORS[0]},
], 'dB');
```

## When to implement

12 biquads covers the vast majority of speaker corrections at 44.1/48 kHz.
The FFT pass is gated on `EQGEN_FFT_BINS > 0` in the header — zero bins means
bypass the FFT path entirely. Only flip the switch when biquad count becomes
the limiting factor.
