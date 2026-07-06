# TODO

## 1. Extract shared modules (deduplicate copy-pasted code) ✅

- [x] `python/io.py` — canonical `load_measurement()` pulled from `end_to_end.py` / `test_sine_sweep.py` / `test_real_fullband_fit.py`
- [x] `python/enhancer_ffi.py` — canonical C FFI wrapper pulled from `end_to_end.py` / `process_music.py` / `test_sine_sweep.py`
- [x] All callers updated to import from these shared modules

## 2. Pick and document the canonical IIR fitter ✅

- [x] `fit_eq_curve_optimized` (golden‑section gain/Q search) is now the **only** fitter, renamed to `fit_eq_curve`
- [x] Deleted `fit_eq_curve` (old basic) and `_refine_fit` (unused scipy refinement)
- [x] Deleted `fit_eq_curve_annealed`
- [x] All callers updated

## 3. Unify the two parallel pipelines ✅

- [x] `Graph.bass_enhancer_preprocess()` delegates to `model.model_gain_needed()`
- [x] Removed duplicate imports

## 4. Write a README ✅

- [x] Project purpose and architecture diagram
- [x] How to take measurements (brown‑noise recording procedure)
- [x] Usage examples for each entry point
- [x] Parameter guide (cutoff, h2, h3, max‑noise, boost)
- [x] Two‑layer architecture (desktop prototyping / ESP32‑C3 deployment)
- [x] Psychoacoustic model explanation

## 5. Canonicalize enhancer parameters

- [ ] Settle on tuned `(h2, h3)` values — `audition.py` defaults to `(0.5, 1.0)`, `enhance.py` to `(0.33, 0.33)`
- [ ] Document any speaker‑class‑specific defaults

## 6. Clean up dead weight ✅

- [x] Removed `python/archive/`
- [x] Removed `fr_math/`
- [x] Removed `.eel` files (JSFX prototypes — superseded by C DSP + LADSPA)
- [x] Removed JamesDSP CSV output (superseded by PipeWire LADSPA)
- [x] Added `src/Makefile`

## 7. Remove hardcoded paths

- [ ] `prototypes/batch_process.py` — absolute paths to external drive mount
- [ ] `eqgen/cli/audition.py` — hardcoded music directory fallback (now `~/Music`)
- [ ] Replace with CLI args or environment variables

## 8. Handle old .txt measurements

- [ ] Decide whether to write a converter or archive them
- [ ] Many speakers (technics, anker, kenrad, corolla, etc.) only have `.txt` FFT dumps — the pipeline only reads `.wav`

## 9. Project restructuring ✅

- [x] Moved Python modules into `eqgen/` package
- [x] Extracted shared pipeline logic into `eqgen/pipeline.py`
- [x] Created `eqgen/cli/` for entry points (`eqgen`, `audition`, `enhance`, `wire`, `export`)
- [x] Consolidated tests under `eqgen/tests/`
- [x] Moved non‑integrated scripts to `prototypes/`
- [x] Moved LADSPA sources to `src/ladspa/`
- [x] Moved `test_enhancer.c` into `src/`
- [x] Removed thin shell wrappers — replaced with Makefile targets
- [x] Added root `Makefile` with targets for all operations
- [x] Added `pyproject.toml` with entry points
- [x] Updated `.gitignore`, `flake.nix`

## 10. Python test suite integration

- [ ] Direct root test scripts (`test_fullband_fit.py` etc.) now merged into `eqgen/tests/` — imports updated
- [ ] Consider migrating from `run_all.py` to `pytest`
- [ ] Add test coverage for the new pipeline module

## 11. Standalone C test

- [ ] `src/test_enhancer.c` — decide whether to integrate, archive, or remove (currently buildable via `make -C src test`)
