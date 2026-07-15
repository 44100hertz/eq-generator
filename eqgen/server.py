#!/usr/bin/env python3
"""
EQGen Web Server — Preset management, pipeline execution, and visualization.

Serves a single-page web UI for assembling, tweaking, and visualizing presets.

Usage:
    python -m eqgen.server
    python -m eqgen.server --port 8080 --no-open
"""

import argparse
import atexit
import io
import json
import logging
import os
import signal
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import webbrowser
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
PID_FILE = ROOT / "output" / "server.pid"

# ── Bottle (lightweight, no extra deps beyond stdlib) ────────────────
# We use Python's built-in http.server + a tiny routing layer so there
# are zero new dependencies.  If FastAPI is installed, prefer it.
try:
    from fastapi import FastAPI, Request, HTTPException
    from fastapi.responses import HTMLResponse, JSONResponse, FileResponse, Response
    from fastapi.staticfiles import StaticFiles
    import uvicorn
    HAS_FASTAPI = True
except ImportError:
    HAS_FASTAPI = False

sys.path.insert(0, str(ROOT))

from eqgen.presets import Preset, PresetManager, PRESETS_DIR, MAX_IIR_BANDS, load_house_curves, list_house_curve_names, get_house_curve
from eqgen.pipeline import run_pipeline
from eqgen.eq_fit import cascade_response_db, BiquadCoeffs, FitResult, fit_eq_curve
from eqgen.web_ui import load_html


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline runner (runs in background thread, writes result to file)
# ═══════════════════════════════════════════════════════════════════════════════

_pipeline_results: dict = {}
_pipeline_lock = threading.Lock()

REPORTS_DIR = ROOT / "output" / "reports"
REPORTS_DIR.mkdir(parents=True, exist_ok=True)


def _run_pipeline_for_preset(preset_dict: dict, task_id: str):
    """Run the full pipeline for a preset and write results."""
    with _pipeline_lock:
        _pipeline_results[task_id] = {"status": "running", "started": time.time()}

    try:
        preset = Preset.from_dict(preset_dict)

        meas_paths = preset.resolve_measurements()
        target_path = preset.resolve_target()
        noise_path = preset.resolve_noise()

        if not meas_paths or (not target_path and not preset.house_curve):
            raise ValueError("Missing measurement or target paths")

        # Resolve house curve adjustment if specified
        curve_data = None
        if preset.house_curve:
            curve_data = get_house_curve(preset.house_curve)

        detailed = run_pipeline(
            meas_paths, target_path or "", noise_path,
            bass_enhancer_cutoff=preset.fc,
            smooth_exponent=preset.smooth_exponent,
            detailed=True,
            house_curve=curve_data,
        )

        freqs = np.array(detailed["freqs"])
        gains_db = np.array(detailed["gains_db"])
        fs = detailed["sample_rate"]
        max_gain_db = detailed.get("max_gain_db", 0.0)

        # ── IIR-only EQ: fit biquads directly to correction curve ─
        from eqgen.dsp import pre_gain_from_max_gain, compute_overboost_ceiling

        if preset.max_bands > 0:
            freqs_hires = np.logspace(np.log10(freqs[0]), np.log10(freqs[-1]), 500)

            # Fit IIR biquads to the full target curve
            fit_result = fit_eq_curve(freqs, gains_db, fs,
                                      max_bands=preset.max_bands,
                                      min_freq=freqs[0], max_freq=freqs[-1],
                                      min_peaking_freq=freqs[0])
            q28_floats = fit_result.biquads
            iir_response_db = cascade_response_db(q28_floats, freqs_hires, fs)

            # Prediction error: IIR fit vs target on same hi-res grid
            target_hires = np.interp(freqs_hires, freqs, gains_db)
            pred_err_db = iir_response_db - target_hires

            pre_gain_lin = float(pre_gain_from_max_gain(max_gain_db))
            coeffs_flat = [v for bc in fit_result.biquads
                           for v in [bc.b0, bc.b1, bc.b2, bc.a1, bc.a2]]
        else:
            freqs_hires = freqs
            fit_result = FitResult()
            iir_response_db = np.zeros(len(freqs_hires))
            target_hires = np.interp(freqs_hires, freqs, gains_db)
            pred_err_db = iir_response_db - target_hires

            pre_gain_lin = 1.0
            coeffs_flat = []
        overboost_ceiling_db, overboost_ceiling_freq = compute_overboost_ceiling(
            h2_amp=h2_amp, h3_amp=h3_amp,
            fc=preset.fc, coeffs=coeffs_flat,
            pre_gain=pre_gain_lin, fs=fs,
        )

        # ── Smart volume: compute effective curves at different volume levels ─
        from eqgen.smart_volume import compute_smart_volume_curves
        pg_db_loud = 20.0 * np.log10(pre_gain_lin)
        sv_result = compute_smart_volume_curves(freqs, gains_db,
                                                  pre_gain_db_loud=float(pg_db_loud),
                                                  speaker_level_db=preset.speaker_level,
                                                  overboost_db=preset.overboost_db)

        result = {
            "status": "complete",
            "task_id": task_id,
            "metrics": {
                "n_bands": len(fit_result.bands),
                "rms_err": float(np.sqrt(np.mean(err_arr ** 2))),
                "bass_err": float(np.max(np.abs(err_arr[err_freqs <= 250]))),
                "mid_err": float(np.max(np.abs(err_arr[(err_freqs > 250) & (err_freqs <= 2000)]))),
                "treble_err": float(np.max(np.abs(err_arr[err_freqs > 2000]))),
                "fc": preset.fc,
                "fs": fs, "max_gain_db": max_gain_db,
                "efficacy": detailed.get("efficacy", {}),
                "overboost_ceiling_db": round(overboost_ceiling_db, 1),
                "overboost_ceiling_freq": round(overboost_ceiling_freq, 1),
            },
            "raw_measurement": raw_meas,
            "raw_target": targ_db,
            "meas_smoothed": detailed["meas_subtracted"],
            "targ_smoothed": detailed["target_resampled"],
            "correction": [{"freq": float(f), "db": float(d)} for f, d in zip(freqs, gains_db)],
            "error": error,
            "iir_fit": [{"freq": float(f), "db": float(d)}
                        for f, d in zip(freqs_hires, iir_response_db)],
            "combined": [{"freq": float(f), "db": float(d)}
                         for f, d in zip(freqs_hires, iir_response_db)],
            "prediction_error": [{"freq": float(f), "db": float(d)}
                                 for f, d in zip(freqs_hires, pred_err_db)],
            "smart_volume": sv_result,
        }

    except Exception as e:
        result = {
            "status": "error",
            "task_id": task_id,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }

    with _pipeline_lock:
        _pipeline_results[task_id] = result


# ── Apply to System (PipeWire) ────────────────────────────────────────

def _apply_preset_to_system(preset_name: str, task_id: str):
    """Run pipeline, generate eq_coeffs.h, rebuild, restart PipeWire chain.

    Loads the preset by name from PresetManager, then delegates pipeline,
    header generation, and wiring to eqgen.cli.wire.
    """
    with _pipeline_lock:
        _pipeline_results[task_id] = {"status": "running", "started": time.time()}

    try:
        from eqgen.cli.wire import (
            run_full_pipeline, generate_eq_header, generate_sfx_header,
            setup_wiring, SRC_DIR,
        )

        preset = pm.load(preset_name)

        curve_data = get_house_curve(preset.house_curve) if preset.house_curve else None

        # 1-2. Run pipeline + IIR fit
        (eq_freqs, bq_q28_44, bq_q28_48, bands_44, bands_48,
         cfg, coeffs_flat) = run_full_pipeline(
            meas_paths=preset.resolve_measurements(),
            target_path=preset.resolve_target() or "",
            noise_path=preset.resolve_noise(),
            fc=preset.fc,
            max_bands=preset.max_bands,
            smooth_exponent=preset.smooth_exponent,
            bluetooth_id=preset.bluetooth_id or None,
            house_curve=curve_data,
            overboost_db=preset.overboost_db,
        )
        cfg["release_secs"] = preset.release

        # 3. Generate eq_coeffs.h + sfx_data.h
        generate_eq_header(bq_q28_44, bq_q28_48, bands_44, bands_48, cfg,
                          bluetooth_id=preset.bluetooth_id or None,
                          speaker_level=preset.speaker_level)
        generate_sfx_header()
        logger.info("Wrote eq_coeffs.h (%d biquads)", len(bands_44))

        # 4. Rebuild filter binary
        for f in ["filter.c", "ladspa/ladspa_wrapper.c"]:
            (SRC_DIR / f).touch()
        r = subprocess.run(
            ["make", "-C", str(SRC_DIR), "filter"],
            capture_output=True, text=True, timeout=60)
        if r.returncode != 0:
            raise RuntimeError(f"Build failed:\n{r.stderr}")
        logger.info("Rebuilt eqgen_filter")

        # 5. Wire — kill stale processes, then setup_wiring
        subprocess.run(["pkill", "-9", "-f", "pw-cat.*--target"], capture_output=True)
        subprocess.run(["pkill", "-9", "-f", "eqgen_filter"], capture_output=True)
        time.sleep(0.3)
        setup_wiring(cfg, coeffs_flat)
        logger.info("PipeWire chain restarted")

        result = {
            "status": "complete",
            "task_id": task_id,
            "message": f"Applied — {len(bands_44)} bands wired to PipeWire",
            "output": {
                "n_bands": len(bands_44),
                "pre_gain": round(cfg["pre_gain"], 3),
                "pre_gain_db": round(20.0 * np.log10(cfg["pre_gain"]), 1),
                "fs": float(cfg["fs"]),
            },
        }

    except Exception as e:
        result = {
            "status": "error",
            "task_id": task_id,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }

    with _pipeline_lock:
        _pipeline_results[task_id] = result


# ── Flash ESP32 ───────────────────────────────────────────────────────

FIRMWARE_DIR = ROOT / "firmware"


def _find_esp32_port() -> Optional[str]:
    """Auto-detect ESP32 USB-serial port.  Returns device path or None."""
    # Prefer stable by-id symlinks (survive replug port renumbering)
    by_id = Path("/dev/serial/by-id")
    if by_id.exists():
        for p in by_id.iterdir():
            name = p.name.lower()
            if any(chip in name for chip in ("cp210", "ch340", "ch341", "ft232")):
                return str(p.resolve())  # resolve symlink
    # Fallback: glob /dev/ttyUSB* or /dev/ttyACM*
    import glob
    for pattern in ("/dev/ttyUSB*", "/dev/ttyACM*"):
        ports = sorted(glob.glob(pattern))
        if ports:
            return ports[0]
    return None


def _flash_esp32(preset_name: str, task_id: str):
    """Run pipeline, generate eq_coeffs.h, idf.py build flash."""
    with _pipeline_lock:
        _pipeline_results[task_id] = {"status": "running", "started": time.time()}

    try:
        from eqgen.cli.wire import run_full_pipeline, generate_eq_header, generate_sfx_header, SRC_DIR

        preset = pm.load(preset_name)

        curve_data = get_house_curve(preset.house_curve) if preset.house_curve else None

        # 1. Run pipeline + IIR fit
        (eq_freqs, bq_q28_44, bq_q28_48, bands_44, bands_48,
         cfg, coeffs_flat) = run_full_pipeline(
            meas_paths=preset.resolve_measurements(),
            target_path=preset.resolve_target() or "",
            noise_path=preset.resolve_noise(),
            fc=preset.fc,
            max_bands=preset.max_bands,
            smooth_exponent=preset.smooth_exponent,
            bluetooth_id=preset.bluetooth_id or None,
            house_curve=curve_data,
            overboost_db=preset.overboost_db,
        )
        cfg["release_secs"] = preset.release

        # 2. Generate eq_coeffs.h + sfx_data.h
        generate_eq_header(bq_q28_44, bq_q28_48, bands_44, bands_48, cfg,
                          bluetooth_id=preset.bluetooth_id or None,
                          speaker_level=preset.speaker_level)
        generate_sfx_header()
        logger.info("Wrote eq_coeffs.h (%d biquads)", len(bands_44))

        # 3. idf.py build
        r = subprocess.run(
            ["idf.py", "build"],
            cwd=str(FIRMWARE_DIR),
            capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            raise RuntimeError(f"idf.py build failed:\n{r.stderr[-1000:]}")
        logger.info("ESP32 firmware built")

        # 4. idf.py flash — retry once on transient serial errors.
        max_attempts = 2
        last_err = ""
        for attempt in range(1, max_attempts + 1):
            flash_cmd = ["idf.py", "flash"]
            port = _find_esp32_port()
            if port:
                flash_cmd.extend(["-p", port])
                logger.info("Flashing on %s (attempt %d/%d)", port, attempt, max_attempts)
            else:
                logger.warning("No ESP32 port auto-detected — using idf.py default")
            r = subprocess.run(
                flash_cmd,
                cwd=str(FIRMWARE_DIR),
                capture_output=True, text=True, timeout=120)
            if r.returncode == 0:
                break
            last_err = r.stderr[-500:] if r.stderr else "unknown error"
            if "Permission denied" in last_err or "not readable" in last_err:
                raise RuntimeError(
                    f"Serial port permission denied. "
                    f"Use 'make flash ARGS={preset.name}' from a terminal, "
                    f"or add yourself to the dialout group.")
            if "serial noise" in last_err.lower() and attempt < max_attempts:
                logger.warning("Transient serial error, retrying…")
                continue
            raise RuntimeError(f"idf.py flash failed:\n{last_err}")
        else:
            raise RuntimeError(f"idf.py flash failed after {max_attempts} attempts:\n{last_err}")

        result = {
            "status": "complete",
            "task_id": task_id,
            "message": f"Flashed — {len(bands_44)} bands",
            "output": {
                "n_bands": len(bands_44),
                "pre_gain": round(cfg["pre_gain"], 3),
                "pre_gain_db": round(20.0 * np.log10(cfg["pre_gain"]), 1),
                "fs": float(cfg["fs"]),
            },
        }

    except Exception as e:
        result = {
            "status": "error",
            "task_id": task_id,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }

    with _pipeline_lock:
        _pipeline_results[task_id] = result


# ═══════════════════════════════════════════════════════════════════════════
# API helpers
# ═══════════════════════════════════════════════════════════════════════════

pm = PresetManager()


def api_list_presets():
    """List all presets with metadata."""
    names = pm.list_presets()
    presets = []
    for name in names:
        try:
            p = pm.load(name)
            presets.append(p.to_dict())
        except Exception:
            presets.append({"name": name, "error": "failed to load"})
    return presets


def api_get_preset(name: str):
    """Get a single preset."""
    p = pm.load(name)
    return p.to_dict()


def api_save_preset(data: dict):
    """Create or update a preset."""
    name = data.get("name")
    if not name:
        raise ValueError("Preset requires a 'name' field")
    preset = Preset.from_dict(data)
    preset.name = name
    pm.save(preset)
    return preset.to_dict()


def api_delete_preset(name: str):
    """Delete a preset."""
    if pm.delete(name):
        return {"deleted": name}
    raise FileNotFoundError(f"Preset '{name}' not found")


def api_scan_measurements():
    """Scan measurement directories."""
    return pm.scan_measurement_dirs()


def api_set_smart_volume(vol: int):
    """Write a volume value (0-127) to the smart-volume control FIFO."""
    ctrl_path = "/tmp/eqgen_sv_fifo"
    try:
        fd = os.open(ctrl_path, os.O_WRONLY)
        os.write(fd, bytes([max(0, min(127, vol))]))
        os.close(fd)
        return {"volume": vol, "status": "ok"}
    except (FileNotFoundError, OSError) as e:
        return {"volume": vol, "status": "no_filter", "error": str(e)}


def api_run_pipeline(data: dict):
    """Start a pipeline run in the background. Returns a task_id."""
    import uuid
    task_id = data.get("task_id") or str(uuid.uuid4())[:8]
    preset_dict = data.get("preset", data)

    thread = threading.Thread(
        target=_run_pipeline_for_preset,
        args=(preset_dict, task_id),
        daemon=True,
    )
    thread.start()
    return {"task_id": task_id, "status": "queued"}


def api_apply_pipeline(data: dict):
    """Start a full apply (build + PipeWire) in the background."""
    import uuid
    task_id = data.get("task_id") or str(uuid.uuid4())[:8]
    preset_name = data.get("preset_name")
    if not preset_name:
        raise ValueError("preset_name required")

    thread = threading.Thread(
        target=_apply_preset_to_system,
        args=(preset_name, task_id),
        daemon=True,
    )
    thread.start()
    return {"task_id": task_id, "status": "queued"}


def api_flash_esp32(data: dict):
    """Start an ESP32 flash (pipeline + header + idf.py build flash)."""
    import uuid
    task_id = data.get("task_id") or str(uuid.uuid4())[:8]
    preset_name = data.get("preset_name")
    if not preset_name:
        raise ValueError("preset_name required")

    thread = threading.Thread(
        target=_flash_esp32,
        args=(preset_name, task_id),
        daemon=True,
    )
    thread.start()
    return {"task_id": task_id, "status": "queued"}


def api_pipeline_status(task_id: str):
    """Get pipeline run status/results."""
    with _pipeline_lock:
        result = _pipeline_results.get(task_id)
    if result is None:
        return {"status": "unknown", "task_id": task_id}
    return result


# ═══════════════════════════════════════════════════════════════════════════
# FastAPI app
# ═══════════════════════════════════════════════════════════════════════════

if HAS_FASTAPI:
    app = FastAPI(title="EQGen", description="Speaker EQ Correction Web UI")

    @app.get("/")
    async def index():
        return HTMLResponse(content=WEB_UI_HTML, status_code=200)

    @app.get("/api/presets")
    async def list_presets():
        return api_list_presets()

    @app.get("/api/presets/{name}")
    async def get_preset(name: str):
        try:
            return api_get_preset(name)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail=f"Preset '{name}' not found")

    @app.post("/api/presets")
    async def save_preset(request: Request):
        data = await request.json()
        try:
            return api_save_preset(data)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

    @app.delete("/api/presets/{name}")
    async def delete_preset(name: str):
        try:
            return api_delete_preset(name)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail=f"Preset '{name}' not found")

    @app.get("/api/measurements")
    async def scan_measurements():
        return api_scan_measurements()

    @app.post("/api/run")
    async def run_pipeline_ep(request: Request):
        data = await request.json()
        return api_run_pipeline(data)

    @app.post("/api/analyze")
    async def analyze_pipeline_ep(request: Request):
        data = await request.json()
        return api_run_pipeline(data)

    @app.post("/api/apply")
    async def apply_pipeline_ep(request: Request):
        data = await request.json()
        return api_apply_pipeline(data)

    @app.post("/api/flash")
    async def flash_esp32_ep(request: Request):
        data = await request.json()
        return api_flash_esp32(data)

    @app.get("/api/house_curves")
    async def list_house_curves_ep():
        curves = load_house_curves()
        return {"curves": curves, "names": list_house_curve_names()}

    @app.get("/api/status/{task_id}")
    async def pipeline_status(task_id: str):
        return api_pipeline_status(task_id)

    @app.post("/api/smart-volume")
    async def set_smart_volume_ep(request: Request):
        data = await request.json()
        vol = int(data.get("volume", 127))
        return api_set_smart_volume(vol)


# ═══════════════════════════════════════════════════════════════════════════
# Bottle / stdlib fallback
# ═══════════════════════════════════════════════════════════════════════════

if not HAS_FASTAPI:
    from http.server import HTTPServer, BaseHTTPRequestHandler
    from urllib.parse import urlparse

    class _Handler(BaseHTTPRequestHandler):
        def _send_json(self, data, status=200):
            body = json.dumps(data).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _read_json(self):
            length = int(self.headers.get("Content-Length", 0))
            if length == 0:
                print(f"[DEBUG] {self.command} {self.path} — empty body", file=sys.stderr, flush=True)
                return {}
            raw = self.rfile.read(length)
            try:
                data = json.loads(raw.decode("utf-8"))
            except json.JSONDecodeError as e:
                print(f"[DEBUG] {self.command} {self.path} — JSON error: {e}", file=sys.stderr, flush=True)
                print(f"[DEBUG] raw body ({len(raw)} bytes): {raw[:300]}", file=sys.stderr, flush=True)
                raise
            if 'preset' in data:
                p = data['preset']
                print(f"[DEBUG] preset keys: {list(p.keys())}", file=sys.stderr, flush=True)
                print(f"[DEBUG] preset.name={p.get('name')!r}", file=sys.stderr, flush=True)
            else:
                print(f"[DEBUG] body keys: {list(data.keys())} (no 'preset' key!)", file=sys.stderr, flush=True)
            return data

        def _route(self, method, path):
            parsed = urlparse(path)
            p = parsed.path

            # Static
            if p == "/" or p == "/index.html":
                body = WEB_UI_HTML.encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return

            # API
            try:
                if p == "/api/presets" and method == "GET":
                    self._send_json(api_list_presets())
                elif p == "/api/presets" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_save_preset(data))
                elif p.startswith("/api/presets/") and method == "GET":
                    name = p.split("/api/presets/", 1)[1]
                    self._send_json(api_get_preset(name))
                elif p.startswith("/api/presets/") and method == "DELETE":
                    name = p.split("/api/presets/", 1)[1]
                    self._send_json(api_delete_preset(name))
                elif p == "/api/measurements" and method == "GET":
                    self._send_json(api_scan_measurements())
                elif p == "/api/house_curves" and method == "GET":
                    curves = load_house_curves()
                    self._send_json({"curves": curves, "names": list_house_curve_names()})
                elif p == "/api/run" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_run_pipeline(data))
                elif p == "/api/analyze" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_run_pipeline(data))
                elif p == "/api/apply" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_apply_pipeline(data))
                elif p == "/api/flash" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_flash_esp32(data))
                elif p == "/api/smart-volume" and method == "POST":
                    data = self._read_json()
                    self._send_json(api_set_smart_volume(int(data.get("volume", 127))))
                elif p.startswith("/api/status/") and method == "GET":
                    task_id = p.split("/api/status/", 1)[1]
                    self._send_json(api_pipeline_status(task_id))
                else:
                    self._send_json({"error": "not found"}, 404)
            except FileNotFoundError as e:
                self._send_json({"error": str(e)}, 404)
            except ValueError as e:
                self._send_json({"error": str(e)}, 400)
            except Exception as e:
                self._send_json({"error": str(e), "traceback": traceback.format_exc()}, 500)

        def do_GET(self):
            self._route("GET", self.path)

        def do_POST(self):
            self._route("POST", self.path)

        def do_DELETE(self):
            self._route("DELETE", self.path)

        def log_message(self, format, *args):
            pass  # suppress logs


# ═══════════════════════════════════════════════════════════════════════════
# Web UI (single-page app, embedded as a string)
# ═══════════════════════════════════════════════════════════════════════════

WEB_UI_HTML = load_html()

# Inject single-source-of-truth constants into the embedded JS
WEB_UI_HTML = WEB_UI_HTML.replace("__MAX_IIR_BANDS__", str(MAX_IIR_BANDS))


# ═══════════════════════════════════════════════════════════════════════════
# CLI entry point
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="EQGen Web Server — preset management and visualization")
    ap.add_argument("--port", type=int, default=8080, help="Server port [8080]")
    ap.add_argument("--host", default="127.0.0.1", help="Bind address [127.0.0.1]")
    ap.add_argument("--no-open", action="store_true", help="Don't open browser")
    ap.add_argument("--reload", action="store_true", help="Auto-reload on code changes (FastAPI only)")
    ap.add_argument("--pid-file", default=str(PID_FILE), help=f"PID file path [{PID_FILE}]")

    args = ap.parse_args()
    pid_path = Path(args.pid_file)
    pid_path.parent.mkdir(parents=True, exist_ok=True)

    # Write PID file + register cleanup
    pid_path.write_text(str(os.getpid()))
    atexit.register(_cleanup_pid, pid_path)

    if HAS_FASTAPI:
        if not args.no_open:
            threading.Timer(1.5, lambda: webbrowser.open(f"http://{args.host}:{args.port}")).start()
        uvicorn.run("eqgen.server:app", host=args.host, port=args.port, log_level="info", reload=args.reload)
    else:
        server = HTTPServer((args.host, args.port), _Handler)
        url = f"http://{args.host}:{args.port}"
        print(f"EQGen server running at {url}")
        print("Press Ctrl+C to stop.")
        if not args.no_open:
            threading.Timer(1.0, lambda: webbrowser.open(url)).start()
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nShutting down...")
            server.shutdown()


def _cleanup_pid(pid_path: Path):
    """Remove PID file if it belongs to this process."""
    try:
        if pid_path.exists():
            pid_path.unlink()
    except Exception:
        pass


if __name__ == "__main__":
    main()
