# eqgen — Speaker EQ Correction with Harmonic Bass Enhancement
#
# Modules:
#   pipeline     — WAV → EQ curve (Welch FFT, kernel smoothing)
#   eq_fit       — IIR biquad design (greedy peaking filter fit)
#   dsp          — Filter design utilities and analysis helpers
#   analysis     — Goertzel, FFT, test signal generation
#   sweep        — Sweep-based analysis through the C enhancer
#   enhancer_ffi — ctypes wrapper for src/enhancer.so
#   io           — WAV I/O and measurement loading
#   iso226       — ISO 226:2023 equal-loudness contour data
#   presets      — Preset system for bundling measurements + parameters
#   process      — Audio file processing through C enhancer
#   smart_volume — Smart volume curve modeling (matches firmware)
#   server       — Web UI server for preset management and visualization
