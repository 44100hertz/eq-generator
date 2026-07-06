#!/usr/bin/env python3
"""
Graph check: run the full EQ pipeline and render a self-contained HTML report
with canvas-based charts of every analysis stage.

Usage:
    python -m eqgen.cli.graph_check technics/standing
    python -m eqgen.cli.graph_check -m meas.wav -t target.wav -n noise.wav
    python -m eqgen.cli.graph_check technics/standing -o /tmp/report.html
"""

import argparse
import json
import os
import subprocess
import sys
import webbrowser
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

from eqgen.pipeline import run_pipeline
from eqgen.eq_fit import cascade_response_db, BiquadCoeffs, fit_eq_curve
from eqgen.quantize import q28_to_float, quantize_biquads_q28

MEAS_DIR = ROOT / "measurements"


# ── Pipeline wrapper ──────────────────────────────────────────────────

def gather_all(speaker_name=None, meas_paths=None, noise_path=None,
               target_path=None, fc=60.0, h2=0.5, h3=1.0,
               max_noise=0.65, max_bands=40):
    """Run pipeline + IIR fit + harmonic model; return everything as a dict for JSON."""
    if speaker_name:
        meas_dir = MEAS_DIR / speaker_name
        if not meas_dir.exists():
            sys.exit(f"ERROR: measurement dir not found: {meas_dir}")
        meas_wavs = sorted(meas_dir.glob("measurement*.wav"))
        if not meas_wavs:
            meas_wavs = sorted(meas_dir.glob("response*.wav"))
        target_wav = meas_dir / "target.wav"
        noise_wav = None
        for nf in ["noise.wav", "noise2.wav", "noise1.wav"]:
            if (meas_dir / nf).exists():
                noise_wav = str(meas_dir / nf)
                break
        if not meas_wavs or not target_wav.exists():
            sys.exit(f"ERROR: need measurement*.wav + target.wav in {meas_dir}")
        meas_paths = [str(p) for p in meas_wavs]
        target_path = str(target_wav)
        noise_path = noise_wav

    if not meas_paths or not target_path:
        sys.exit("ERROR: need --measurement and --target (or a speaker name)")

    # ── Run pipeline with detailed intermediate data ──────────────────
    detailed = run_pipeline(
        meas_paths, target_path, noise_path,
        max_noise=max_noise, bass_enhancer_cutoff=fc, h2=h2, h3=h3,
        detailed=True,
    )

    freqs = np.array(detailed["freqs"])
    gains_db = np.array(detailed["gains_db"])
    fs = detailed["sample_rate"]

    # ── IIR fit ───────────────────────────────────────────────────────
    fit = fit_eq_curve(freqs, gains_db, fs, max_bands=max_bands,
                       min_freq=freqs[0], max_freq=freqs[-1],
                       min_peaking_freq=freqs[0])
    bq_q28 = quantize_biquads_q28(fit.biquads)
    q28_floats = [BiquadCoeffs(b0=q28_to_float(b.b0), b1=q28_to_float(b.b1),
                               b2=q28_to_float(b.b2), a1=q28_to_float(b.a1),
                               a2=q28_to_float(b.a2)) for b in bq_q28]
    fitted_db = cascade_response_db(q28_floats, freqs, fs)

    # ── High-resolution IIR response for smoother graph 4 ────────────
    freqs_hires = np.logspace(np.log10(freqs[0]), np.log10(freqs[-1]), 500)
    fitted_db_hires = cascade_response_db(q28_floats, freqs_hires, fs)
    iir_fit = [{"freq": float(f), "db": float(d)}
               for f, d in zip(freqs_hires, fitted_db_hires)]

    correction = [{"freq": float(f), "db": float(d)} for f, d in zip(freqs, gains_db)]

    # ── Normalize: align measurement to target in midrange (500-2000 Hz) ──
    raw_meas = detailed["raw_measurement"]
    raw_targ = detailed["raw_target"]

    def _midrange_offset(meas, targ):
        """Mean dB difference between target and measurement in 500-2000 Hz."""
        targ_by_freq = {p["freq"]: p["db"] for p in targ}
        diffs = []
        for p in meas:
            f = p["freq"]
            if 500 <= f <= 2000 and f in targ_by_freq:
                diffs.append(targ_by_freq[f] - p["db"])
        return float(np.mean(diffs)) if diffs else 0.0

    align_db = _midrange_offset(raw_meas, raw_targ)
    raw_meas_aligned = [{"freq": p["freq"], "db": p["db"] + align_db} for p in raw_meas]

    # ── Error (full-resolution, after midrange alignment) ────────────
    error = []
    targ_by_freq = {p["freq"]: p["db"] for p in raw_targ}
    for p in raw_meas_aligned:
        t_db = targ_by_freq.get(p["freq"])
        if t_db is not None:
            error.append({"freq": p["freq"], "db": p["db"] - t_db})

    # ── Predicted error: IIR fit + error → residual after correction ──
    #    error = meas − target  (graph 2),  IIR_fit ≈ target − meas  (graph 4)
    #    so adding them should cancel to ~0 dB.
    iir_freqs = freqs
    iir_dbs = fitted_db
    predicted_error = []
    for pt in error:
        iir_val = float(np.interp(pt["freq"], iir_freqs, iir_dbs))
        predicted_error.append({"freq": pt["freq"], "db": iir_val + pt["db"]})

    # ── Error metrics (computed from full-resolution error) ───────────
    err_arr = np.array([e["db"] for e in error])
    err_freqs = np.array([e["freq"] for e in error])
    metrics = {
        "n_bands": len(fit.bands),
        "rms_err": float(np.sqrt(np.mean(err_arr ** 2))),
        "bass_err": float(np.max(np.abs(err_arr[err_freqs <= 250]))),
        "mid_err": float(np.max(np.abs(err_arr[(err_freqs > 250) & (err_freqs <= 2000)]))),
        "treble_err": float(np.max(np.abs(err_arr[err_freqs > 2000]))),
        "fc": fc, "h2": h2, "h3": h3,
        "fs": fs, "speaker": speaker_name or "custom",
        "align_db": round(align_db, 1),
    }

    return {
        "metrics": metrics,
        "raw_measurement": raw_meas_aligned,
        "raw_target": detailed["raw_target"],
        "noise_cv": detailed["noise_cv"],
        "noise_floor": detailed["noise_floor"],
        "correction": correction,
        "error": error,
        "predicted_error": predicted_error,
        "iir_fit": iir_fit,
    }


# ── HTML report generator ─────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>EQGen — Pipeline Report</title>
<style>
  :root {
    --bg: #0d1117; --surface: #161b22; --border: #30363d;
    --text: #c9d1d9; --muted: #8b949e; --accent: #58a6ff;
    --red: #f85149; --green: #3fb950; --yellow: #d2991d;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font: 14px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: var(--bg); color: var(--text); padding: 24px; }
  h1 { font-size: 22px; margin-bottom: 4px; }
  .subtitle { color: var(--muted); margin-bottom: 24px; font-size: 13px; }
  .metrics { display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 32px; }
  .metric { background: var(--surface); border: 1px solid var(--border);
            border-radius: 8px; padding: 12px 16px; min-width: 130px; }
  .metric .label { font-size: 11px; color: var(--muted); text-transform: uppercase;
                   letter-spacing: .5px; }
  .metric .value { font-size: 18px; font-weight: 600; margin-top: 2px; }
  .metric .value.good { color: var(--green); }
  .metric .value.warn { color: var(--yellow); }
  .metric .value.bad { color: var(--red); }
  .chart { background: var(--surface); border: 1px solid var(--border);
           border-radius: 8px; padding: 16px; margin-bottom: 20px; }
  .chart h2 { font-size: 14px; font-weight: 600; color: var(--muted);
              margin-bottom: 8px; }
  canvas { display: block; width: 100%; height: 300px; }
  .legend { display: flex; gap: 16px; flex-wrap: wrap; margin-top: 8px; font-size: 12px; }
  .legend span { display: flex; align-items: center; gap: 4px; }
  .legend .swatch { display: inline-block; width: 12px; height: 3px; border-radius: 2px; }
  .tooltip { position: fixed; background: rgba(22,27,34,0.95); border: 1px solid var(--border);
             border-radius: 6px; padding: 6px 10px; font-size: 12px; pointer-events: none;
             display: none; z-index: 100; color: var(--text); }
</style>
</head>
<body>
<h1>EQGen Pipeline Report</h1>
<div class="subtitle">__SPEAKER__ &mdash; fc=__FC__Hz h2=__H2__ h3=__H3__ &mdash; __N_BANDS__ biquads @ __FS__Hz &mdash; alignment __ALIGN_DB__ dB</div>

<div class="metrics">
  <div class="metric"><div class="label">RMS Error</div><div class="value __RMS_ERR_CLS__">__RMS_ERR__ dB</div></div>
  <div class="metric"><div class="label">Bass Error ≤250Hz</div><div class="value __BASS_ERR_CLS__">__BASS_ERR__ dB</div></div>
  <div class="metric"><div class="label">Mid Error 250-2k</div><div class="value __MID_ERR_CLS__">__MID_ERR__ dB</div></div>
  <div class="metric"><div class="label">Treble Error >2k</div><div class="value __TREBLE_ERR_CLS__">__TREBLE_ERR__ dB</div></div>
  <div class="metric"><div class="label">Biquad Bands</div><div class="value">__N_BANDS__</div></div>
</div>

<div class="chart"><h2>1. Raw Measurement vs Target</h2>
  <canvas id="c_raw"></canvas>
  <div class="legend"><span><i class="swatch" style="background:#58a6ff"></i> Measurement</span>
    <span><i class="swatch" style="background:#d2991d"></i> Target</span></div>
</div>

<div class="chart"><h2>2. Error (Measurement &minus; Target)</h2>
  <canvas id="c_error"></canvas>
  <div class="legend"><span><i class="swatch" style="background:#f85149"></i> Error dB</span>
    <span><i class="swatch" style="background:rgba(255,255,255,0.1);border:1px dashed #30363d;height:0"></i> 0 dB line</span></div>
</div>

<div class="chart"><h2>3. Noise Floor (per-bin CV)</h2>
  <canvas id="c_noise"></canvas>
  <div class="legend"><span><i class="swatch" style="background:#8b949e"></i> CV (coefficient of variation)</span></div>
</div>

<div class="chart"><h2>4. IIR Fit vs Ideal Correction</h2>
  <canvas id="c_iir_fit"></canvas>
  <div class="legend"><span><i class="swatch" style="background:#3fb950"></i> Ideal correction</span>
    <span><i class="swatch" style="background:#bc8cff;height:1px;border:1px dashed #bc8cff"></i> IIR biquad fit</span></div>
</div>

<div class="chart"><h2>5. Predicted Error (corrected output &minus; target)</h2>
  <canvas id="c_predicted"></canvas>
  <div class="legend"><span><i class="swatch" style="background:#58a6ff"></i> Predicted output</span>
    <span><i class="swatch" style="background:rgba(255,255,255,0.12);border:1px dashed #30363d;height:0"></i> 0 dB (flat)</span></div>
</div>

<div class="tooltip" id="tooltip"></div>

<script>
const DATA = __JSON_DATA__;
const COLORS = ['#58a6ff','#3fb950','#f85149','#d2991d','#bc8cff','#ff7b72','#79c0ff','#a5d6ff'];
const CANVAS_DPI = 2;

function initChart(id, traces, yLabel, opts) {
  opts = opts || {};
  const canvas = document.getElementById(id);
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  const dpr = CANVAS_DPI;
  const parent = canvas.parentElement;
  const W = parent.clientWidth - 32;
  const H = 300;
  canvas.width = W * dpr;
  canvas.height = H * dpr;
  canvas.style.width = W + 'px';
  canvas.style.height = H + 'px';
  ctx.scale(dpr, dpr);

  const pad = {l: 52, r: 20, t: 12, b: 32};
  const pw = W - pad.l - pad.r;
  const ph = H - pad.t - pad.b;

  // value ranges
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
  const ySpan = yMax - yMin || 1;
  yMin -= ySpan * 0.05;
  yMax += ySpan * 0.05;

  const toX = f => pad.l + (Math.log10(f) - Math.log10(xMin)) / (Math.log10(xMax) - Math.log10(xMin)) * pw;
  const toY = v => pad.t + (yMax - v) / (yMax - yMin) * ph;

  // background
  ctx.fillStyle = '#161b22';
  ctx.fillRect(0, 0, W, H);

  // grid
  ctx.strokeStyle = 'rgba(255,255,255,0.04)';
  ctx.lineWidth = 1;
  const yStep = _niceStep(yMax - yMin, 6);
  const yStart = Math.ceil(yMin / yStep) * yStep;
  for (let y = yStart; y <= yMax; y += yStep) {
    const py = toY(y);
    ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(W - pad.r, py); ctx.stroke();
    ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'right';
    ctx.fillText(y.toFixed(1), pad.l - 6, py + 3);
  }
  // x axis labels (log decades)
  const decades = [20, 30, 40, 50, 60, 80, 100, 200, 300, 400, 500, 600, 800,
                   1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 14000, 20000];
  ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
  decades.forEach(f => {
    if (f >= xMin && f <= xMax) {
      const px = toX(f);
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, H - pad.b); ctx.stroke();
      const label = f >= 1000 ? (f/1000).toFixed(f%1000===0?0:1)+'k' : f;
      ctx.fillText(label, px, H - pad.b + 14);
    }
  });

  // zero line
  if (yMin < 0 && yMax > 0) {
    ctx.strokeStyle = 'rgba(255,255,255,0.12)';
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(pad.l, toY(0)); ctx.lineTo(W - pad.r, toY(0)); ctx.stroke();
    ctx.setLineDash([]);
  }

  // traces
  traces.forEach((t, i) => {
    const color = t.color || COLORS[i % COLORS.length];
    ctx.strokeStyle = color;
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

  // y axis label
  ctx.save();
  ctx.fillStyle = '#8b949e'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
  ctx.translate(10, pad.t + ph/2);
  ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  // hover tooltip
  const tip = document.getElementById('tooltip');
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
      tip.innerHTML = '<b>' + fLabel + '</b> &nbsp; ' + (best.name||'') + ' &nbsp; ' + vLabel + ' dB';
    }
  });
  canvas.addEventListener('mouseleave', function() { tip.style.display = 'none'; });
}

function _niceStep(span, target) {
  const rough = span / target;
  const exp = Math.pow(10, Math.floor(Math.log10(rough)));
  const mant = rough / exp;
  const nice = mant <= 1.5 ? 1 : mant <= 3 ? 2 : mant <= 7 ? 5 : 10;
  return nice * exp;
}

// ── Render charts ──
initChart('c_raw', [
  {name:'Measurement', data:DATA.raw_measurement, color:COLORS[0]},
  {name:'Target', data:DATA.raw_target, color:COLORS[3]},
], 'dB');

initChart('c_error', [
  {name:'Error', data:DATA.error, color:COLORS[2]},
], 'dB', {yMin: -30, yMax: 30});

initChart('c_noise', [
  {name:'CV', data:DATA.noise_cv, color:COLORS[5]||'#8b949e'},
], 'CV', {yMin: 0, yMax: 1.2});

initChart('c_iir_fit', [
  {name:'Ideal', data:DATA.correction, color:COLORS[1]},
  {name:'IIR fit', data:DATA.iir_fit, color:COLORS[4], dash:[6,3], width:1.4},
], 'dB');

initChart('c_predicted', [
  {name:'Predicted', data:DATA.predicted_error, color:COLORS[0]},
], 'dB');
</script>
</body>
</html>"""


def _err_cls(val):
    if val < 1.0:
        return "good"
    elif val < 3.0:
        return "warn"
    return "bad"


def generate_html(data: dict) -> str:
    """Fill the HTML template with pipeline data."""
    m = data["metrics"]
    html = HTML_TEMPLATE
    html = html.replace("__JSON_DATA__", json.dumps(data, indent=2))
    html = html.replace("__SPEAKER__", str(m.get("speaker", "custom")))
    html = html.replace("__FC__", str(m["fc"]))
    html = html.replace("__H2__", str(m["h2"]))
    html = html.replace("__H3__", str(m["h3"]))
    html = html.replace("__N_BANDS__", str(m["n_bands"]))
    html = html.replace("__FS__", f'{m["fs"]:.0f}')
    html = html.replace("__ALIGN_DB__", f'{m.get("align_db", 0):+.1f}')
    for key in ["rms_err", "bass_err", "mid_err", "treble_err"]:
        v = m[key]
        html = html.replace(f"__{key.upper()}__", f"{v:.1f}")
        html = html.replace(f"__{key.upper()}_CLS__", _err_cls(v))
    return html


# ── CLI ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Run full pipeline → self-contained HTML graph report.")
    ap.add_argument("speaker", nargs="?", default=None,
                    help="Speaker name (under measurements/)")
    ap.add_argument("-m", "--measurement", nargs="+", default=None,
                    help="Measurement WAV file(s)")
    ap.add_argument("-n", "--noise", default=None, help="Noise WAV file")
    ap.add_argument("-t", "--target", default=None, help="Target WAV file")
    ap.add_argument("-o", "--output", default=None, help="Output HTML path")
    ap.add_argument("--fc", type=float, default=60.0, help="Bass enhancer cutoff Hz")
    ap.add_argument("--h2", type=float, default=0.5, help="2nd harmonic amplitude")
    ap.add_argument("--h3", type=float, default=1.0, help="3rd harmonic amplitude")
    ap.add_argument("--max-noise", type=float, default=0.65)
    ap.add_argument("--max-bands", type=int, default=40)
    ap.add_argument("--no-open", action="store_true", help="Don't open in browser")
    args = ap.parse_args()

    if not args.speaker and not (args.measurement and args.target):
        ap.error("Need speaker name or -m/-t")

    print("Running pipeline...")
    data = gather_all(
        speaker_name=args.speaker,
        meas_paths=args.measurement,
        noise_path=args.noise,
        target_path=args.target,
        fc=args.fc, h2=args.h2, h3=args.h3,
        max_noise=args.max_noise, max_bands=args.max_bands,
    )

    out_path = args.output or f"output/{args.speaker or 'custom'}_report.html"
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    html = generate_html(data)
    with open(out_path, "w") as f:
        f.write(html)

    print(f"Wrote {out_path}  ({len(data['correction'])} EQ points, "
          f"{data['metrics']['n_bands']} biquads)")

    if not args.no_open:
        url = "file://" + os.path.abspath(out_path)
        webbrowser.open(url)


if __name__ == "__main__":
    main()
