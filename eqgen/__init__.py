# eqgen — Speaker EQ Correction with Harmonic Bass Enhancement
#
# Modules:
#   pipeline     — WAV → EQ curve (Welch FFT, kernel smoothing)
#   eq_fit       — IIR biquad design (greedy peaking filter fit)
#   quantize     — Q4.28 fixed-point biquad quantization
#   dsp          — Filter design utilities and analysis helpers
#   analysis     — Goertzel, FFT, test signal generation
#   sweep        — Sweep-based analysis through the C enhancer
#   enhancer_ffi — ctypes wrapper for src/enhancer.so
#   io           — WAV I/O and measurement loading
