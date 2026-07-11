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
from eqgen.eq_fit import cascade_response_db, BiquadCoeffs, fit_eq_curve
from eqgen.quantize import q28_to_float, quantize_biquads_q28


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
            h2=preset.h2, h3=preset.h3,
            smooth_exponent=preset.smooth_exponent,
            detailed=True,
            house_curve=curve_data,
        )

        freqs = np.array(detailed["freqs"])
        gains_db = np.array(detailed["gains_db"])
        fs = detailed["sample_rate"]
        max_gain_db = detailed.get("max_gain_db", 0.0)

        # ── IIR-only EQ: fit biquads directly to correction curve ─
        freqs_hires = np.logspace(np.log10(freqs[0]), np.log10(freqs[-1]), 500)

        # Fit IIR biquads to the full target curve
        fit_result = fit_eq_curve(freqs, gains_db, fs,
                                  max_bands=preset.max_bands,
                                  min_freq=freqs[0], max_freq=freqs[-1],
                                  min_peaking_freq=freqs[0])
        bq_q28 = quantize_biquads_q28(fit_result.biquads)
        q28_floats = [BiquadCoeffs(b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                                   b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                                   a2=q28_to_float(b.a2)) for b in bq_q28]
        iir_response_db = cascade_response_db(q28_floats, freqs_hires, fs)

        # Prediction error: IIR fit vs target on same hi-res grid
        target_hires = np.interp(freqs_hires, freqs, gains_db)
        pred_err_db = iir_response_db - target_hires

        # Raw measurement (as-is, no level shift) and fitted target line.
        # The target is always pink or brown noise — a straight line in dB
        # vs log-freq.  Fit it, then shift to measurement’s midrange level.
        raw_meas = detailed["raw_measurement"]
        raw_targ = detailed["raw_target"]
        targ_arr = np.array([(p["freq"], p["db"]) for p in raw_targ])
        log_f_targ = np.log10(targ_arr[:, 0])
        slope, intercept = np.polyfit(log_f_targ, targ_arr[:, 1], 1)
        targ_fit = slope * log_f_targ + intercept

        meas_arr = np.array([(p["freq"], p["db"]) for p in raw_meas])
        mm = (meas_arr[:, 0] >= 500) & (meas_arr[:, 0] <= 2000)
        tm = (targ_arr[:, 0] >= 500) & (targ_arr[:, 0] <= 2000)
        meas_mid = float(np.mean(meas_arr[mm, 1])) if mm.any() else 0.0
        targ_mid = float(np.mean(targ_fit[tm])) if tm.any() else 0.0
        offset = meas_mid - targ_mid

        targ_db = [{"freq": float(targ_arr[i, 0]),
                     "db": float(targ_fit[i] + offset)}
                    for i in range(len(targ_arr))]

        # Error: fitted target - measurement (at raw measurement resolution)
        meas_log_f = np.log10(meas_arr[:, 0])
        targ_at_meas = slope * meas_log_f + intercept + offset
        error = [{"freq": float(meas_arr[i, 0]),
                  "db": float(targ_at_meas[i] - meas_arr[i, 1])}
                 for i in range(len(meas_arr))]

        err_arr = np.array([e["db"] for e in error])
        err_freqs = np.array([e["freq"] for e in error])

        # ── Smart volume: compute effective curves at different volume levels ─
        from eqgen.pipeline import compute_smart_volume_curves, pre_gain_from_max_gain
        pg_db_loud = 20.0 * np.log10(pre_gain_from_max_gain(max_gain_db))
        sv_result = compute_smart_volume_curves(freqs, gains_db, preset.h2, preset.h3,
                                                  pre_gain_db_loud=float(pg_db_loud))

        result = {
            "status": "complete",
            "task_id": task_id,
            "metrics": {
                "n_bands": len(fit_result.bands),
                "rms_err": float(np.sqrt(np.mean(err_arr ** 2))),
                "bass_err": float(np.max(np.abs(err_arr[err_freqs <= 250]))),
                "mid_err": float(np.max(np.abs(err_arr[(err_freqs > 250) & (err_freqs <= 2000)]))),
                "treble_err": float(np.max(np.abs(err_arr[err_freqs > 2000]))),
                "fc": preset.fc, "h2": preset.h2, "h3": preset.h3,
                "fs": fs, "max_gain_db": max_gain_db,
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

def _apply_preset_to_system(preset_dict: dict, task_id: str):
    """Run pipeline, generate eq_coeffs.h, rebuild, restart PipeWire chain.

    Delegates pipeline, header generation, and wiring to eqgen.cli.wire.
    setup_wiring() correctly tears down any existing chain before recreating
    it — no orphaned pw-cat or filter processes.
    """
    with _pipeline_lock:
        _pipeline_results[task_id] = {"status": "running", "started": time.time()}

    try:
        from eqgen.cli.wire import (
            run_full_pipeline, generate_eq_header, generate_sfx_header,
            setup_wiring, SRC_DIR,
        )

        preset = Preset.from_dict(preset_dict)
        meas_paths = preset.resolve_measurements()
        target_path = preset.resolve_target()
        noise_path = preset.resolve_noise()

        if not meas_paths or (not target_path and not preset.house_curve):
            raise ValueError("Missing measurement or target paths")

        curve_data = None
        if preset.house_curve:
            curve_data = get_house_curve(preset.house_curve)

        # 1-2. Run pipeline + IIR fit (delegates to wire.py)
        (eq_freqs, bq_q28_44, bq_q28_48, bands_44, bands_48,
         cfg, coeffs_flat) = run_full_pipeline(
            meas_paths=meas_paths,
            target_path=target_path or "",
            noise_path=noise_path,
            fc=preset.fc or 60.0,
            h2=preset.h2,
            h3=preset.h3,
            max_bands=preset.max_bands,
            smooth_exponent=preset.smooth_exponent,
            bluetooth_id=preset.bluetooth_id or None,
            house_curve=curve_data,
        )
        cfg["release_secs"] = preset.release
        cfg["limiter_release_secs"] = preset.limiter_release

        # 3. Generate eq_coeffs.h + sfx_data.h (delegates to wire.py)
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

        # 5. Wire — kill any stale processes first, then let setup_wiring
        #    handle the full teardown + recreation cycle.
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


def _flash_esp32(preset_dict: dict, task_id: str):
    """Run pipeline, generate eq_coeffs.h, idf.py build flash."""
    with _pipeline_lock:
        _pipeline_results[task_id] = {"status": "running", "started": time.time()}

    try:
        from eqgen.cli.wire import run_full_pipeline, generate_eq_header, generate_sfx_header, SRC_DIR

        preset = Preset.from_dict(preset_dict)
        meas_paths = preset.resolve_measurements()
        target_path = preset.resolve_target()
        noise_path = preset.resolve_noise()

        if not meas_paths or (not target_path and not preset.house_curve):
            raise ValueError("Missing measurement or target paths")

        curve_data = None
        if preset.house_curve:
            curve_data = get_house_curve(preset.house_curve)

        # 1. Run pipeline + IIR fit
        (eq_freqs, bq_q28_44, bq_q28_48, bands_44, bands_48,
         cfg, coeffs_flat) = run_full_pipeline(
            meas_paths=meas_paths,
            target_path=target_path or "",
            noise_path=noise_path,
            fc=preset.fc or 60.0,
            h2=preset.h2,
            h3=preset.h3,
            max_bands=preset.max_bands,
            smooth_exponent=preset.smooth_exponent,
            bluetooth_id=preset.bluetooth_id or None,
            house_curve=curve_data,
        )
        cfg["release_secs"] = preset.release
        cfg["limiter_release_secs"] = preset.limiter_release

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

        # 4. idf.py flash
        flash_cmd = ["idf.py", "flash"]
        port = _find_esp32_port()
        if port:
            flash_cmd.extend(["-p", port])
            logger.info("Flashing on %s", port)
        else:
            logger.warning("No ESP32 port auto-detected — using idf.py default")
        r = subprocess.run(
            flash_cmd,
            cwd=str(FIRMWARE_DIR),
            capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            err = r.stderr[-500:] if r.stderr else "unknown error"
            if "Permission denied" in err or "not readable" in err:
                raise RuntimeError(
                    f"Serial port permission denied. "
                    f"Use 'make flash ARGS={preset.name}' from a terminal, "
                    f"or add yourself to the dialout group.")
            raise RuntimeError(f"idf.py flash failed:\n{err}")

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
    preset_dict = data.get("preset", data)

    thread = threading.Thread(
        target=_apply_preset_to_system,
        args=(preset_dict, task_id),
        daemon=True,
    )
    thread.start()
    return {"task_id": task_id, "status": "queued"}


def api_flash_esp32(data: dict):
    """Start an ESP32 flash (pipeline + header + idf.py build flash)."""
    import uuid
    task_id = data.get("task_id") or str(uuid.uuid4())[:8]
    preset_dict = data.get("preset", data)

    thread = threading.Thread(
        target=_flash_esp32,
        args=(preset_dict, task_id),
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

WEB_UI_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>EQGen — Preset Manager</title>
<style>
  :root {
    --bg: #0d1117; --surface: #161b22; --border: #30363d;
    --text: #c9d1d9; --muted: #8b949e; --accent: #58a6ff;
    --red: #f85149; --green: #3fb950; --yellow: #d2991d;
    --purple: #bc8cff;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: var(--bg); color: var(--text); height: 100vh; display: flex; overflow: hidden; }
  #sidebar { width: 280px; min-width: 280px; background: var(--surface);
             border-right: 1px solid var(--border); display: flex; flex-direction: column; overflow: hidden; }
  #sidebar .header { padding: 16px; border-bottom: 1px solid var(--border); }
  #sidebar .header h1 { font-size: 18px; color: var(--accent); }
  #sidebar .header .sub { font-size: 11px; color: var(--muted); }
  #sidebar .actions { padding: 12px; display: flex; gap: 6px; border-bottom: 1px solid var(--border); }
  #sidebar .preset-list { flex: 1; overflow-y: auto; padding: 4px 0; }
  #sidebar .preset-item { padding: 10px 16px; cursor: pointer; border-left: 3px solid transparent;
                          display: flex; justify-content: space-between; align-items: center;
                          transition: background .15s; }
  #sidebar .preset-item:hover { background: rgba(255,255,255,0.04); }
  #sidebar .preset-item.active { background: rgba(88,166,255,0.1); border-left-color: var(--accent); }
  #sidebar .preset-item .name { font-weight: 600; font-size: 13px; }
  #sidebar .preset-item .desc { font-size: 11px; color: var(--muted); margin-top: 2px; }
  #sidebar .preset-item .del { opacity: 0; color: var(--red); font-size: 16px; padding: 2px 6px;
                                border: none; background: none; cursor: pointer; border-radius: 4px; }
  #sidebar .preset-item:hover .del { opacity: 0.6; }
  #sidebar .preset-item .del:hover { opacity: 1; background: rgba(248,81,73,0.15); }
  #main { flex: 1; display: flex; flex-direction: column; overflow: hidden; }
  #main .toolbar { padding: 12px 20px; border-bottom: 1px solid var(--border);
                   display: flex; align-items: center; gap: 12px; background: var(--surface); }
  #main .toolbar h2 { font-size: 16px; flex: 1; }
  #main .content { flex: 1; overflow-y: auto; padding: 20px; }
  .form-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; max-width: 800px; }
  .form-group { display: flex; flex-direction: column; gap: 4px; }
  .form-group.full { grid-column: 1 / -1; }
  .form-group label { font-size: 11px; color: var(--muted); text-transform: uppercase;
                      letter-spacing: .5px; font-weight: 600; }
  .form-group input, .form-group select, .form-group textarea {
    background: var(--bg); border: 1px solid var(--border); border-radius: 6px;
    padding: 8px 10px; color: var(--text); font-size: 13px; font-family: inherit;
    transition: border-color .15s; }
  .form-group input:focus, .form-group select:focus, .form-group textarea:focus {
    outline: none; border-color: var(--accent); }
  .form-group textarea { resize: vertical; min-height: 60px; }
  button { padding: 8px 16px; border-radius: 6px; border: 1px solid var(--border);
           background: var(--surface); color: var(--text); font-size: 13px; cursor: pointer;
           font-family: inherit; transition: all .15s; }
  button:hover { background: rgba(255,255,255,0.06); }
  button.primary { background: var(--accent); color: #000; border-color: var(--accent); font-weight: 600; }
  button.primary:hover { opacity: 0.9; }
  button.danger { color: var(--red); border-color: var(--red); }
  button.danger:hover { background: rgba(248,81,73,0.15); }
  button:disabled { opacity: 0.4; cursor: not-allowed; }
  .chart { background: var(--surface); border: 1px solid var(--border);
           border-radius: 8px; padding: 16px; margin-bottom: 16px; }
  .chart h3 { font-size: 13px; font-weight: 600; color: var(--muted); margin-bottom: 8px; }
  canvas { display: block; width: 100%; height: 280px; }
  .legend { display: flex; gap: 14px; flex-wrap: wrap; margin-top: 8px; font-size: 12px; }
  .legend span { display: flex; align-items: center; gap: 4px; }
  .legend .swatch { display: inline-block; width: 14px; height: 3px; border-radius: 2px; }
  #status-bar { padding: 4px 16px; font-size: 11px; color: var(--muted);
                border-top: 1px solid var(--border); background: var(--surface); }
  #status-bar.running { color: var(--yellow); }
  #status-bar.error { color: var(--red); }
  .tooltip { position: fixed; background: rgba(22,27,34,0.95); border: 1px solid var(--border);
             border-radius: 6px; padding: 6px 10px; font-size: 12px; pointer-events: none;
             display: none; z-index: 100; color: var(--text); }
  .empty { text-align: center; color: var(--muted); padding: 60px 20px; }
  .empty h3 { font-size: 16px; margin-bottom: 8px; }
  .empty p { font-size: 13px; max-width: 400px; margin: 0 auto 20px; }
</style>
</head>
<body>
<div id="sidebar">
  <div class="header">
    <h1>EQGen</h1>
    <div class="sub">Speaker EQ &amp; Bass Enhancement</div>
  </div>
  <div class="actions">
    <button onclick="newPreset()" style="flex:1">+ New Preset</button>
    <button onclick="scanMeasurements()" style="flex:1">Scan Mics</button>
  </div>
  <div class="preset-list" id="presetList"></div>
</div>

<div id="main">
  <div class="toolbar">
    <h2 id="mainTitle">Select a preset or create a new one</h2>
    <div style="display:flex;gap:6px">
      <button class="primary" id="btnApply" onclick="applyToSystem()" disabled>&#x26A1; Apply to System</button>
      <button id="btnFlash" onclick="flashESP32()" disabled>&#x1F4E1; Flash ESP32</button>
      <button id="btnSave" onclick="savePreset()" disabled>&#x1F4BE; Save</button>
    </div>
  </div>
  <div id="liveBar" style="padding:8px 20px;border-bottom:1px solid var(--border);background:var(--bg);display:flex;align-items:center;gap:12px;font-size:12px">
    <span style="color:var(--muted);font-weight:600">&#x1F50A; Live Smart Volume</span>
    <input type="range" id="liveSlider" min="0" max="127" value="127"
      style="flex:1;max-width:200px;accent-color:var(--accent)"
      oninput="onLiveSlider(this.value)">
    <span id="liveLabel" style="color:var(--accent);min-width:24px;text-align:right">127</span>
    <span id="liveStatus" style="color:var(--muted);font-size:11px"></span>
  </div>
  <div class="content" id="content">
    <div id="editor"></div>
    <div id="results"></div>
  </div>
  <div id="status-bar">Ready</div>
</div>

<script>
// ═══════════════════════════════════════════════════════════════════════
// State
// ═══════════════════════════════════════════════════════════════════════
let state = {
  presets: [],
  activePreset: null,
  pipelineResult: null,
  analyzeTaskId: null,
  applyTaskId: null,
  analyzePending: null,
  measurements: [],
  houseCurves: null,
};

const COLORS = ['#58a6ff','#3fb950','#f85149','#d2991d','#bc8cff','#ff7b72','#79c0ff','#a5d6ff'];
const DPR = 2;
const DEBOUNCE_MS = 600;

// ═══════════════════════════════════════════════════════════════════════
// API helpers
// ═══════════════════════════════════════════════════════════════════════
async function api(url, opts = {}) {
  const resp = await fetch(url, opts);
  if (!resp.ok) {
    const err = await resp.json().catch(() => ({}));
    throw new Error(err.error || err.detail || `HTTP ${resp.status}`);
  }
  return resp.json();
}

async function loadPresets() {
  state.presets = await api('/api/presets');
  renderPresetList();
}

async function loadPreset(name) {
  const p = await api(`/api/presets/${name}`);
  state.activePreset = p;
  state.pipelineResult = null;
  document.getElementById('results').innerHTML = '';
  renderEditor();
  setStatus('Ready');
  // Auto-trigger analysis
  scheduleAnalyze();
}

async function savePreset() {
  const p = readForm();
  if (!p || !p.name) { alert('Preset name is required'); return; }
  await api('/api/presets', { method: 'POST', headers: {'Content-Type':'application/json'}, body: JSON.stringify(p) });
  setStatus(`Saved "${p.name}"`);
  await loadPresets();
  // Don't reload — preserve current form state
  state.activePreset = p;
  renderPresetList();
}

async function deletePreset(name) {
  if (!confirm(`Delete preset "${name}"?`)) return;
  await api(`/api/presets/${name}`, { method: 'DELETE' });
  state.activePreset = null;
  state.pipelineResult = null;
  document.getElementById('results').innerHTML = '';
  await loadPresets();
  renderEditor();
  setStatus(`Deleted "${name}"`);
}

async function scanMeasurements() {
  state.measurements = await api('/api/measurements');
  if (state.activePreset) renderEditor();
}

// ═══════════════════════════════════════════════════════════════════════
// Auto-analyze (debounced, reactive)
// ═══════════════════════════════════════════════════════════════════════

function scheduleAnalyze() {
  // Debounce: reset timer on every call
  if (state.analyzePending) clearTimeout(state.analyzePending);
  state.analyzePending = setTimeout(doAnalyze, DEBOUNCE_MS);
}

async function doAnalyze() {
  state.analyzePending = null;
  const p = readForm();
  if (!p || !p.name) return;
  if (!p.measurements || p.measurements.length === 0 || (!p.target && !p.house_curve)) {
    document.getElementById('results').innerHTML = '';
    return;
  }

  setStatus('Analyzing...', 'running');
  try {
    const resp = await api('/api/analyze', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({preset: p})
    });
    state.analyzeTaskId = resp.task_id;
    pollAnalyzeStatus();
  } catch(e) {
    setStatus(`Analyze failed: ${e.message}`, 'error');
  }
}

async function pollAnalyzeStatus() {
  const tid = state.analyzeTaskId;
  if (!tid) return;
  try {
    const result = await api(`/api/status/${tid}`);
    if (result.status === 'complete') {
      state.pipelineResult = result;
      renderResults();
      const m = result.metrics || {};
      setStatus(`${m.n_bands||0} bands, ${(m.rms_err||0).toFixed(1)} dB RMS error`);
      state.analyzeTaskId = null;
    } else if (result.status === 'error') {
      state.pipelineResult = {error: result.error};
      renderResults();
      setStatus(`Analyze error: ${result.error}`, 'error');
      state.analyzeTaskId = null;
    } else {
      const elapsed = result.started ? ((Date.now() / 1000 - result.started)).toFixed(0) : '?';
      setStatus(`Analyzing... ${elapsed}s`, 'running');
      setTimeout(pollAnalyzeStatus, 1000);
    }
  } catch(e) {
    setStatus(`Poll error: ${e.message}`, 'error');
    state.analyzeTaskId = null;
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Apply to System (PipeWire)
// ═══════════════════════════════════════════════════════════════════════

async function applyToSystem() {
  const p = readForm();
  if (!p || !p.name) { setStatus('Preset name is required', 'error'); return; }
  if (!p.measurements || p.measurements.length === 0 || (!p.target && !p.house_curve)) {
    setStatus('Measurement and target paths required', 'error'); return;
  }

  const btn = document.getElementById('btnApply');
  btn.disabled = true;
  btn.textContent = '\u26A1 Applying...';
  setStatus('Building DSP and wiring PipeWire...', 'running');

  try {
    const resp = await api('/api/apply', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({preset: p})
    });
    state.applyTaskId = resp.task_id;
    pollApplyStatus();
  } catch(e) {
    setStatus(`Apply failed: ${e.message}`, 'error');
    btn.disabled = false;
    btn.textContent = '\u26A1 Apply to System';
  }
}

// ═══════════════════════════════════════════════════════════════════════
// Flash ESP32
// ═══════════════════════════════════════════════════════════════════════

async function flashESP32() {
  const p = readForm();
  if (!p || !p.name) { setStatus('Preset name is required', 'error'); return; }
  if (!p.measurements || p.measurements.length === 0 || (!p.target && !p.house_curve)) {
    setStatus('Measurement and target paths required', 'error'); return;
  }

  if (!confirm('Flash ESP32 with this EQ curve? This will rebuild the firmware and flash via USB.')) return;

  const btn = document.getElementById('btnFlash');
  btn.disabled = true;
  btn.textContent = '\u1F4E1 Flashing...';
  setStatus('Building firmware and flashing ESP32...', 'running');

  try {
    const resp = await api('/api/flash', {
      method: 'POST',
      headers: {'Content-Type':'application/json'},
      body: JSON.stringify({preset: p})
    });
    state.flashTaskId = resp.task_id;
    pollFlashStatus();
  } catch(e) {
    setStatus(`Flash failed: ${e.message}`, 'error');
    btn.disabled = false;
    btn.textContent = '\u1F4E1 Flash ESP32';
  }
}

async function pollFlashStatus() {
  const tid = state.flashTaskId;
  if (!tid) return;
  const btn = document.getElementById('btnFlash');
  try {
    const result = await api(`/api/status/${tid}`);
    if (result.status === 'complete') {
      btn.disabled = false;
      btn.textContent = '\u1F4E1 Flash ESP32';
      const m = result.output || {};
      setStatus(`\u2705 Flashed — ${m.n_bands||'?'} bands, ${result.message||''}`);
      state.flashTaskId = null;
    } else if (result.status === 'error') {
      btn.disabled = false;
      btn.textContent = '\u1F4E1 Flash ESP32';
      setStatus(`Flash error: ${result.error}`, 'error');
      state.flashTaskId = null;
    } else {
      const elapsed = result.started ? ((Date.now() / 1000 - result.started)).toFixed(0) : '?';
      setStatus(`Building / flashing... ${elapsed}s`, 'running');
      setTimeout(pollFlashStatus, 2000);
    }
  } catch(e) {
    btn.disabled = false;
    btn.textContent = '\u1F4E1 Flash ESP32';
    setStatus(`Flash poll error: ${e.message}`, 'error');
    state.flashTaskId = null;
  }
}

async function pollApplyStatus() {
  const tid = state.applyTaskId;
  if (!tid) return;
  const btn = document.getElementById('btnApply');
  try {
    const result = await api(`/api/status/${tid}`);
    if (result.status === 'complete') {
      btn.disabled = false;
      btn.textContent = '\u26A1 Apply to System';
      const m = result.output || {};
      setStatus(`\u2705 Applied — ${m.n_bands||'?'} bands, ${result.message||''}`);
      state.applyTaskId = null;
    } else if (result.status === 'error') {
      btn.disabled = false;
      btn.textContent = '\u26A1 Apply to System';
      setStatus(`Apply error: ${result.error}`, 'error');
      state.applyTaskId = null;
    } else {
      const elapsed = result.started ? ((Date.now() / 1000 - result.started)).toFixed(0) : '?';
      setStatus(`Building / wiring... ${elapsed}s`, 'running');
      setTimeout(pollApplyStatus, 1500);
    }
  } catch(e) {
    btn.disabled = false;
    btn.textContent = '\u26A1 Apply to System';
    setStatus(`Apply poll error: ${e.message}`, 'error');
    state.applyTaskId = null;
  }
}

function setStatus(msg, cls = '') {
  const bar = document.getElementById('status-bar');
  bar.textContent = typeof msg === 'string' ? msg : JSON.stringify(msg);
  bar.className = cls;
}

// ── Live smart-volume slider (sends to eqgen_filter via FIFO) ────

let _svDebounce = null;
function onLiveSlider(vol) {
  const label = document.getElementById('liveLabel');
  const status = document.getElementById('liveStatus');
  if (label) label.textContent = vol;

  // Debounce: send only after slider stops for 80ms
  if (_svDebounce) clearTimeout(_svDebounce);
  _svDebounce = setTimeout(async () => {
    try {
      const r = await api('/api/smart-volume', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({volume: parseInt(vol)}),
      });
      if (status) {
        status.textContent = r.status === 'ok' ? '✓ live' : '✗ ' + (r.error || 'no filter');
        status.style.color = r.status === 'ok' ? 'var(--green)' : 'var(--red)';
      }
    } catch(e) {
      if (status) { status.textContent = '✗ ' + e.message; status.style.color = 'var(--red)'; }
    }
  }, 80);
}

// ═══════════════════════════════════════════════════════════════════════
// Rendering
// ═══════════════════════════════════════════════════════════════════════

function renderPresetList() {
  const el = document.getElementById('presetList');
  if (state.presets.length === 0) {
    el.innerHTML = '<div class="empty" style="padding:40px 16px"><p style="font-size:12px">No presets yet.<br>Click "+ New Preset" to create one.</p></div>';
    return;
  }
  el.innerHTML = state.presets.map(p => {
    const active = state.activePreset && state.activePreset.name === p.name ? ' active' : '';
    const fcStr = p.fc != null ? `fc=${esc(String(p.fc))}` : '';
    const h2Str = p.h2 != null ? `h2=${esc(String(p.h2))}` : '';
    const h3Str = p.h3 != null ? `h3=${esc(String(p.h3))}` : '';
    const paramStr = [fcStr, h2Str, h3Str].filter(Boolean).join(' ');
    return `<div class="preset-item${active}" onclick="loadPreset('${esc(p.name)}')">
      <div><div class="name">${esc(p.name)}</div>
        ${paramStr ? `<div class="desc">${paramStr}</div>` : '<div class="desc">no parameters</div>'}
      </div>
      <button class="del" onclick="event.stopPropagation();deletePreset('${esc(p.name)}')" title="Delete">&times;</button>
    </div>`;
  }).join('');
}

function readForm() {
  const name = document.getElementById('fName')?.value.trim();
  if (!name) return null;
  return {
    name,
    description: document.getElementById('fDesc')?.value || '',
    measurements: (document.getElementById('fMeas')?.value || '').split('\n').map(s=>s.trim()).filter(Boolean),
    target: document.getElementById('fTarget')?.value || '',
    house_curve: document.getElementById('fHouseCurve')?.value || '',
    noise: document.getElementById('fNoise')?.value || null,
    fc: numVal('fFc'),
    h2: numVal('fH2', 0.5),
    h3: numVal('fH3', 1.0),
    max_bands: intVal('fMaxBands', __MAX_IIR_BANDS__),
    smooth_exponent: numVal('fSmooth', 1.0),
    release: numVal('fRelease', 0.2),
    limiter_release: numVal('fLimiter', 0.049),
    bluetooth_id: document.getElementById('fBtId')?.value || '',
    speaker_level: intVal('fSpkLevel', 60),
    high_rolloffs: [],
    low_rolloffs: [],
  };
}

function numVal(id, fallback) {
  const el = document.getElementById(id);
  if (!el || el.value === '') return fallback !== undefined ? fallback : null;
  const v = parseFloat(el.value);
  return isNaN(v) ? (fallback !== undefined ? fallback : null) : v;
}
function intVal(id, fallback) {
  const el = document.getElementById(id);
  if (!el || el.value === '') return fallback !== undefined ? fallback : null;
  const v = parseInt(el.value);
  return isNaN(v) ? (fallback !== undefined ? fallback : null) : v;
}

function renderEditor() {
  const p = state.activePreset;
  const editor = document.getElementById('editor');
  const title = document.getElementById('mainTitle');
  const btnSave = document.getElementById('btnSave');
  const btnApply = document.getElementById('btnApply');
  const btnFlash = document.getElementById('btnFlash');

  if (!p) {
    title.textContent = 'Select a preset or create a new one';
    btnSave.disabled = true;
    btnApply.disabled = true;
    btnFlash.disabled = true;
    editor.innerHTML = `<div class="empty">
      <h3>Welcome to EQGen</h3>
      <p>Load a preset to see EQ curves. Tweak parameters — charts update automatically. Click "Apply to System" to push the EQ to your live PipeWire chain.</p>
      <button class="primary" onclick="newPreset()">+ New Preset</button>
    </div>`;
    document.getElementById('presetList').querySelectorAll('.active').forEach(el => el.classList.remove('active'));
    return;
  }

  title.textContent = `Preset: ${esc(p.name)}`;
  btnSave.disabled = false;
  btnApply.disabled = false;
  btnFlash.disabled = false;
  document.getElementById('presetList').querySelectorAll('.active').forEach(el => el.classList.remove('active'));
  document.querySelector(`.preset-item[onclick*="${p.name}"]`)?.classList.add('active');

  const measStr = (p.measurements || []).join('\n');
  const measOpts = buildMeasOptions();

  editor.innerHTML = `
    <div class="form-grid">
      <div class="form-group"><label>Preset Name</label><input id="fName" value="${esc(p.name)}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Description</label><input id="fDesc" value="${esc(p.description || '')}"></div>
      <div class="form-group full"><label>Measurement WAV files (one per line)</label>
        <textarea id="fMeas" rows="3" oninput="scheduleAnalyze()">${esc(measStr)}</textarea>
        ${measOpts}
      </div>
      <div class="form-group"><label>House Curve</label><select id="fHouseCurve" onchange="onHouseCurveChange()">${buildHouseCurveOptions(p.house_curve || 'flat')}</select></div>
      <div class="form-group"><label>Target WAV</label><input id="fTarget" value="${esc(p.target || '')}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Noise WAV (optional)</label><input id="fNoise" value="${esc(p.noise || '')}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Cutoff fc (Hz)</label><input id="fFc" type="number" step="1" min="20" max="200" value="${p.fc ?? ''}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>H2 Amp</label><input id="fH2" type="number" step="0.01" min="0" max="2" value="${p.h2 ?? 0.5}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>H3 Amp</label><input id="fH3" type="number" step="0.01" min="0" max="2" value="${p.h3 ?? 1.0}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Max Bands</label><input id="fMaxBands" type="number" step="1" min="1" max="8" value="${p.max_bands ?? __MAX_IIR_BANDS__}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Smooth Exponent</label><input id="fSmooth" type="number" step="0.1" min="0" max="4" value="${p.smooth_exponent ?? 1.0}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Release (s)</label><input id="fRelease" type="number" step="0.01" min="0.01" max="2" value="${p.release ?? 0.2}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Limiter Release (s)</label><input id="fLimiter" type="number" step="0.001" min="0.01" max="0.5" value="${p.limiter_release ?? 0.049}" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Bluetooth ID</label><input id="fBtId" value="${esc(p.bluetooth_id || '')}" placeholder="e.g. Living Room Speaker" oninput="scheduleAnalyze()"></div>
      <div class="form-group"><label>Speaker Level (dB)</label><input id="fSpkLevel" type="number" step="1" min="20" max="100" value="${p.speaker_level ?? 60}" oninput="scheduleAnalyze()"></div>
    </div>
  `;
}

function buildMeasOptions() {
  if (!state.measurements.length) return '';
  let html = '<div style="margin-top:6px;font-size:11px;color:var(--muted)">Available: ';
  state.measurements.forEach(d => {
    html += `<span style="cursor:pointer;color:var(--accent);margin-right:8px" title="${d.path}" onclick="fillFromDir('${d.path}')">${esc(d.name)}</span>`;
  });
  html += '</div>';
  return html;
}

window.fillFromDir = function(dir) {
  const dirEntry = state.measurements.find(d => d.path === dir);
  if (!dirEntry) return;
  const meas = dirEntry.wavs.filter(w => w.kind === 'measurement').map(w => w.path);
  const targ = dirEntry.wavs.find(w => w.kind === 'target');
  const noise = dirEntry.wavs.filter(w => w.kind === 'noise');
  if (meas.length) document.getElementById('fMeas').value = meas.join('\n');
  if (targ) document.getElementById('fTarget').value = targ.path;
  if (noise.length) document.getElementById('fNoise').value = noise[noise.length-1].path;
  scheduleAnalyze();
};

async function loadHouseCurves() {
  try {
    const resp = await api('/api/house_curves');
    state.houseCurves = resp;
    if (state.activePreset) renderEditor();
  } catch(e) {
    console.error('Failed to load house curves:', e);
  }
}

function buildHouseCurveOptions(selected) {
  if (!state.houseCurves || !state.houseCurves.names) {
    return '<option value="">(loading...)</option>';
  }
  return state.houseCurves.names.map(n =>
    `<option value="${esc(n)}" ${n === selected ? 'selected' : ''}>${esc(n)}</option>`
  ).join('');
}

function onHouseCurveChange() {
  scheduleAnalyze();
}

function newPreset() {
  state.activePreset = { name: 'new-preset', measurements: [], target: '', house_curve: 'flat', fc: 60, h2: 0.5, h3: 1.0, max_bands: __MAX_IIR_BANDS__, smooth_exponent: 1.0, release: 0.2, limiter_release: 0.049, bluetooth_id: '', speaker_level: 60 };
  state.pipelineResult = null;
  document.getElementById('results').innerHTML = '';
  document.getElementById('presetList').querySelectorAll('.active').forEach(el => el.classList.remove('active'));
  renderEditor();
}

// ═══════════════════════════════════════════════════════════════════════
// Results / Charts
// ═══════════════════════════════════════════════════════════════════════

function renderResults() {
  const r = state.pipelineResult;
  const resultsDiv = document.getElementById('results');

  // Show loading state while waiting
  if (!r && state.analyzePending !== undefined) {
    return; // keep whatever was there
  }

  if (!r || r.status === 'error') {
    const errMsg = r?.error ? esc(String(r.error)) : 'Waiting for analysis...';
    resultsDiv.innerHTML = r?.error
      ? `<div class="chart"><h3>Analysis Error</h3><p style="color:var(--red);white-space:pre-wrap;font-family:monospace;font-size:11px">${errMsg}</p></div>`
      : '<div class="chart"><p style="color:var(--muted);text-align:center;padding:40px">Enter measurement paths and target. Charts will appear here.</p></div>';
    return;
  }
  const n_bands = (r.metrics && r.metrics.n_bands) || '—';
  const sv = r.smart_volume || null;

  let svHtml = '';
  if (sv && sv.curves) {
    svHtml = `
    <div class="chart">
      <h3>4. Smart Volume — Shelf Boost &amp; Pre-Gain by Volume
        <span style="font-weight:400;font-size:12px;color:var(--muted);margin-left:12px">
          Vol: <input type="range" id="svSlider" min="0" max="127" value="127"
            style="vertical-align:middle;width:160px;accent-color:var(--accent)"
            oninput="onSvSlider(this.value)">
          <span id="svLabel" style="display:inline-block;width:50px;text-align:right">127</span>
        </span>
        <span id="svPgLabel" style="font-weight:400;font-size:12px;color:var(--yellow);margin-left:16px"></span>
      </h3>
      <canvas id="cSmartVol"></canvas>
      <div class="legend">
        <span><i class="swatch" style="background:rgba(63,185,80,0.55)"></i> Correction (loud, vol→127)</span>
        <span><i class="swatch" style="background:rgba(88,166,255,0.7);height:0;border-bottom:2px dashed rgba(88,166,255,0.7)"></i> + Shelf boost (current vol)</span>
      </div>
    </div>`;
  }

  resultsDiv.innerHTML = `
    <div class="chart"><h3>1. Measurement &amp; Target</h3><canvas id="cMeas"></canvas>
      <div class="legend"><span><i class="swatch" style="background:#58a6ff"></i> Raw measurement</span><span><i class="swatch" style="background:#58a6ff;height:0;border-bottom:1px dashed #58a6ff"></i> Meas. smoothed</span><span><i class="swatch" style="background:#d2991d"></i> Target (fit)</span></div></div>
    <div class="chart"><h3>2. IIR EQ Fit (${n_bands} biquad bands)</h3><canvas id="cFft"></canvas>
      <div class="legend"><span><i class="swatch" style="background:rgba(63,185,80,0.55)"></i> Correction curve</span><span><i class="swatch" style="background:rgba(188,140,255,0.55);height:0;border:1px dashed rgba(188,140,255,0.55)"></i> IIR fit</span><span><i class="swatch" style="background:rgba(255,120,80,0.8);height:3px"></i> Combined</span></div></div>
    <div class="chart"><h3>3. Predicted Response Error (IIR fit − correction, ideally flat at 0 dB)</h3><canvas id="cError"></canvas>
      <div class="legend"><span><i class="swatch" style="background:#f85149"></i> Residual error</span></div></div>
    ${svHtml}
  `;

  requestAnimationFrame(() => {
    requestAnimationFrame(() => {
      drawChart('cMeas', [
        {name:'Raw meas',data:r.raw_measurement,color:'#58a6ff',width:1.8},
        {name:'Meas. smoothed',data:r.meas_smoothed,color:'#58a6ff',dash:[4,2],width:1.2},
        {name:'Target (fit)',data:r.targ_smoothed,color:'#d2991d',width:1.8},
      ], 'dB');
      drawChart('cFft', [
        {name:'Correction',data:r.correction,color:'rgba(63,185,80,0.55)',width:1.8},
        {name:'IIR fit',data:r.iir_fit,color:'rgba(188,140,255,0.55)',dash:[6,3],width:1.2},
        {name:'Combined',data:r.combined,color:'rgba(255,120,80,0.8)',width:2.4},
      ], 'dB');
      drawChart('cError', [
        {name:'Residual error',data:r.prediction_error,color:COLORS[2],width:1.4},
      ], 'dB', {yMin: -3, yMax: 3});
      if (sv && sv.curves) {
        state._svData = sv;
        onSvSlider(127);
      }
    });
  });
}

// ── Smart volume slider ──────────────────────────────────────────────

function onSvSlider(vol) {
  const sv = state._svData;
  if (!sv || !sv.curves) return;

  const label = document.getElementById('svLabel');
  if (label) label.textContent = vol;

  const pgLabel = document.getElementById('svPgLabel');

  const t = vol / 127;

  const curves = sv.curves;
  let idxLo = 0, idxHi = sv.curves.length - 1;
  for (let i = 0; i < curves.length - 1; i++) {
    if (t >= curves[i].t && t <= curves[i + 1].t) {
      idxLo = i; idxHi = i + 1; break;
    }
  }
  const a = (t - curves[idxLo].t) / (curves[idxHi].t - curves[idxLo].t + 0.0001);
  const pgDb = curves[idxLo].pre_gain_db + a * (curves[idxHi].pre_gain_db - curves[idxLo].pre_gain_db);
  const sdb = curves[idxLo].shelf_max_db + a * (curves[idxHi].shelf_max_db - curves[idxLo].shelf_max_db);

  if (pgLabel) {
    pgLabel.textContent = 'pre-gain: ' + pgDb.toFixed(1) + ' dB  |  shelf boost: +' + sdb.toFixed(1) + ' dB';
  }

  const freqs = curves[idxLo].freqs;
  const interp = [];
  for (let i = 0; i < freqs.length; i++) {
    const vLo = curves[idxLo].effective_db[i];
    const vHi = curves[idxHi].effective_db[i];
    interp.push({freq: freqs[i], db: vLo + a * (vHi - vLo)});
  }

  const loudCurve = sv.curves[sv.curves.length-1].effective_db.map(
    (d, i) => ({freq: sv.curves[0].freqs[i], db: d}));

  drawChart('cSmartVol', [
    {name:'Correction (loud)', data: loudCurve, color:'rgba(63,185,80,0.55)', width:1.6},
    {name:'+ Shelf (vol=' + vol + ')', data: interp, color:'rgba(88,166,255,0.7)', dash:[5,3], width:2.0},
  ], 'dB');
}

// ═══════════════════════════════════════════════════════════════════════
// Canvas chart drawing
// ═══════════════════════════════════════════════════════════════════════

function drawChart(canvasId, traces, yLabel, opts = {}) {
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  const dpr = DPR;
  const parent = canvas.parentElement;
  const W = parent.clientWidth - 32;
  const H = 280;
  canvas.width = W * dpr;
  canvas.height = H * dpr;
  canvas.style.width = W + 'px';
  canvas.style.height = H + 'px';
  ctx.scale(dpr, dpr);

  const pad = {l:48, r:16, t:10, b:28};
  const pw = W - pad.l - pad.r;
  const ph = H - pad.t - pad.b;

  let yMin = opts.yMin !== undefined ? opts.yMin : Infinity;
  let yMax = opts.yMax !== undefined ? opts.yMax : -Infinity;
  let xMin = Infinity, xMax = -Infinity;
  traces.forEach(t => {
    t.data.forEach(p => {
      if (p.freq < xMin) xMin = p.freq;
      if (p.freq > xMax) xMax = p.freq;
      const v = p.db !== undefined ? p.db : p.cv;
      if (opts.yMin === undefined && v < yMin) yMin = v;
      if (opts.yMax === undefined && v > yMax) yMax = v;
    });
  });
  if (!isFinite(yMin)) yMin = -30;
  if (!isFinite(yMax)) yMax = 10;
  const span = yMax - yMin || 1;
  yMin -= span * 0.05;
  yMax += span * 0.05;

  const toX = f => pad.l + (Math.log10(f) - Math.log10(xMin)) / (Math.log10(xMax) - Math.log10(xMin)) * pw;
  const toY = v => pad.t + (yMax - v) / (yMax - yMin) * ph;

  ctx.fillStyle = '#161b22';
  ctx.fillRect(0, 0, W, H);

  ctx.strokeStyle = 'rgba(255,255,255,0.04)';
  ctx.lineWidth = 1;
  const yStep = niceStep(yMax - yMin, 6);
  for (let y = Math.ceil(yMin / yStep) * yStep; y <= yMax; y += yStep) {
    const py = toY(y);
    ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(W - pad.r, py); ctx.stroke();
    ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'right';
    ctx.fillText(y.toFixed(1), pad.l - 6, py + 3);
  }
  const decades = [20,30,40,50,60,80,100,200,300,400,500,600,800,1000,2000,3000,4000,5000,6000,8000,10000,14000,20000];
  ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
  decades.forEach(f => {
    if (f >= xMin && f <= xMax) {
      const px = toX(f);
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, H - pad.b); ctx.stroke();
      ctx.fillText(f >= 1000 ? (f/1000).toFixed(f%1000===0?0:1)+'k' : f, px, H - pad.b + 14);
    }
  });

  if (yMin < 0 && yMax > 0) {
    ctx.strokeStyle = 'rgba(255,255,255,0.12)';
    ctx.setLineDash([4,4]);
    ctx.beginPath(); ctx.moveTo(pad.l, toY(0)); ctx.lineTo(W - pad.r, toY(0)); ctx.stroke();
    ctx.setLineDash([]);
  }

  traces.forEach(t => {
    ctx.strokeStyle = t.color || COLORS[0];
    ctx.lineWidth = t.width || 1.8;
    if (t.dash) ctx.setLineDash(t.dash);
    ctx.beginPath();
    let first = true;
    t.data.forEach(p => {
      const px = toX(p.freq);
      const py = toY(p.db !== undefined ? p.db : p.cv);
      if (first) { ctx.moveTo(px, py); first = false; }
      else ctx.lineTo(px, py);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  });

  ctx.save();
  ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
  ctx.translate(10, pad.t + ph/2);
  ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  // Hover tooltip
  const tip = document.getElementById('tooltip');
  if (!tip) return;
  canvas.addEventListener('mousemove', function(e) {
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    if (mx < pad.l || mx > W - pad.r || my < pad.t || my > H - pad.b) {
      tip.style.display = 'none'; return;
    }
    const logF = Math.log10(xMin) + (mx - pad.l) / pw * (Math.log10(xMax) - Math.log10(xMin));
    const freq = Math.pow(10, logF);
    let best = null, bestDist = Infinity;
    traces.forEach(t => {
      t.data.forEach(p => {
        const d = Math.abs(p.freq - freq);
        if (d < bestDist) { bestDist = d; best = {freq: p.freq, value: p.db !== undefined ? p.db : p.cv, name: t.name}; }
      });
    });
    if (best) {
      const vu = best.value !== undefined ? best.value : 0;
      tip.style.display = 'block';
      tip.style.left = (e.clientX + 12) + 'px';
      tip.style.top = (e.clientY - 10) + 'px';
      const fLabel = best.freq >= 1000 ? (best.freq/1000).toFixed(1)+'kHz' : best.freq.toFixed(0)+'Hz';
      const vLabel = vu.toFixed(1);
      tip.innerHTML = '<b>' + fLabel + '</b> &nbsp; ' + (best.name||'') + ' &nbsp; ' + vLabel + ' ' + yLabel;
    }
  });
  canvas.addEventListener('mouseleave', function() { tip.style.display = 'none'; });
}

function niceStep(span, target) {
  const rough = span / target;
  const exp = Math.pow(10, Math.floor(Math.log10(rough)));
  const mant = rough / exp;
  const nice = mant <= 1.5 ? 1 : mant <= 3 ? 2 : mant <= 7 ? 5 : 10;
  return nice * exp;
}

function esc(s) {
  if (s === null || s === undefined) return '';
  if (typeof s !== 'string') return esc(JSON.stringify(s));
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}

window.onerror = function(msg, url, line, col, err) {
  setStatus('JS Error: ' + msg + ' (line ' + line + ')', 'error');
  console.error('Unhandled error:', err || msg);
  return false;
};

// ═══════════════════════════════════════════════════════════════════════
// Init
// ═══════════════════════════════════════════════════════════════════════
renderEditor();
loadPresets().catch(e => setStatus('Failed to load presets: ' + e.message, 'error'));
scanMeasurements().catch(e => setStatus('Failed to scan measurements: ' + e.message, 'error'));
loadHouseCurves().catch(e => setStatus('Failed to load house curves: ' + e.message, 'error'));
</script>
<div class="tooltip" id="tooltip"></div>
</body>
</html>

"""

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
