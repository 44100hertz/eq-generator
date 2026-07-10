"""
Preset system for EQGen — bundles measurement files, target, noise, and all
tuning parameters into a single JSON file for reproducible EQ curves.

Presets are stored in the `presets/` directory as .json files.

Usage:
    from eqgen.presets import PresetManager
    pm = PresetManager()
    pm.list_presets()
    pm.load("my-speaker")
    pm.save(preset_dict)
"""

import json
import os
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parent.parent
PRESETS_DIR = ROOT / "presets"

# ── Global limit: maximum IIR biquad bands across the entire program ──
MAX_IIR_BANDS = 8


@dataclass
class Preset:
    """A complete EQGen preset that drives `run_pipeline()`.

    All path fields are relative to the project root (eqgen/..).
    """
    name: str
    description: str = ""
    measurements: List[str] = field(default_factory=list)
    target: str = ""
    noise: Optional[str] = None
    fc: Optional[float] = None          # bass enhancer cutoff Hz
    h2: float = 0.5                      # 2nd harmonic amplitude
    h3: float = 1.0                      # 3rd harmonic amplitude
    max_bands: int = MAX_IIR_BANDS       # max IIR biquad bands
    smooth_exponent: float = 1.0         # CV smoothing aggressiveness
    release: float = 0.2                 # envelope release (s)
    limiter_release: float = 0.049       # harmonic limiter release (s)
    bluetooth_id: str = ""               # Bluetooth device ID for speaker identity
    default_volume: int = 32              # AVRCP default volume on boot (0–127)
    high_rolloffs: List[List[float]] = field(default_factory=list)
    low_rolloffs: List[List[float]] = field(default_factory=list)
    pre_gain: Optional[float] = None     # computed pre-gain override

    def to_dict(self) -> dict:
        d = asdict(self)
        # Omit None values for cleaner JSON
        return {k: v for k, v in d.items() if v is not None}

    @classmethod
    def from_dict(cls, d: dict) -> "Preset":
        # Filter to known fields
        known = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in d.items() if k in known}
        # Convert list of lists for rolloffs
        if "high_rolloffs" in filtered:
            filtered["high_rolloffs"] = [list(r) if isinstance(r, (list, tuple)) else r
                                          for r in filtered["high_rolloffs"]]
        if "low_rolloffs" in filtered:
            filtered["low_rolloffs"] = [list(r) if isinstance(r, (list, tuple)) else r
                                         for r in filtered["low_rolloffs"]]
        return cls(**filtered)

    def resolve_measurements(self) -> List[str]:
        """Return absolute paths to measurement WAV files."""
        return [str(ROOT / p) for p in self.measurements]

    def resolve_target(self) -> Optional[str]:
        """Return absolute path to target WAV file."""
        if self.target:
            return str(ROOT / self.target)
        return None

    def resolve_noise(self) -> Optional[str]:
        """Return absolute path to noise WAV file."""
        if self.noise:
            return str(ROOT / self.noise)
        return None


class PresetManager:
    """Manage preset JSON files in the presets/ directory."""

    def __init__(self, presets_dir: Optional[Path] = None):
        self.presets_dir = Path(presets_dir) if presets_dir else PRESETS_DIR
        self.presets_dir.mkdir(parents=True, exist_ok=True)

    def list_presets(self) -> List[str]:
        """Return list of preset names (without .json extension)."""
        names = []
        for f in sorted(self.presets_dir.glob("*.json")):
            names.append(f.stem)
        return names

    def load(self, name: str) -> Preset:
        """Load a preset by name. Raises FileNotFoundError if missing."""
        path = self.presets_dir / f"{name}.json"
        with open(path, "r") as f:
            data = json.load(f)
        preset = Preset.from_dict(data)
        preset.name = name  # ensure name matches filename
        return preset

    def save(self, preset: Preset) -> Path:
        """Save a preset to disk. Returns the file path."""
        path = self.presets_dir / f"{preset.name}.json"
        data = preset.to_dict()
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
            f.write("\n")
        return path

    def delete(self, name: str) -> bool:
        """Delete a preset by name. Returns True if deleted, False if missing."""
        path = self.presets_dir / f"{name}.json"
        if path.exists():
            path.unlink()
            return True
        return False

    def exists(self, name: str) -> bool:
        """Check if a preset exists."""
        return (self.presets_dir / f"{name}.json").exists()

    def scan_measurement_dirs(self) -> List[dict]:
        """Scan measurements/ for directories containing WAV files.

        Returns list of {name, path, wavs: [{name, path, kind}]}.
        """
        meas_root = ROOT / "measurements"
        if not meas_root.exists():
            return []

        results = []
        for d in sorted(meas_root.rglob("*")):
            if not d.is_dir():
                continue
            wavs = sorted(d.glob("*.wav"))
            if not wavs:
                continue

            rel_dir = str(d.relative_to(ROOT))
            entries = []
            for w in wavs:
                rel = str(w.relative_to(ROOT))
                kind = "measurement"
                if w.stem == "target":
                    kind = "target"
                elif w.stem == "noise" or w.stem.startswith("noise"):
                    kind = "noise"
                elif w.stem.startswith("measurement") or w.stem.startswith("response"):
                    kind = "measurement"
                entries.append({
                    "name": w.name,
                    "path": rel,
                    "kind": kind,
                })
            results.append({
                "name": str(d.relative_to(meas_root)),
                "path": rel_dir,
                "wavs": entries,
            })
        return results

    def guess_preset_from_measurement_dir(self, dir_name: str) -> Optional[Preset]:
        """Create a preset by guessing from a measurement directory.

        dir_name is relative to measurements/, e.g. "technics/standing".
        """
        meas_root = ROOT / "measurements"
        d = meas_root / dir_name
        if not d.is_dir():
            return None

        meas_wavs = sorted(d.glob("measurement*.wav"))
        if not meas_wavs:
            meas_wavs = sorted(d.glob("response*.wav"))
        target_wav = d / "target.wav"
        noise_wav = None
        for nf in ["noise2.wav", "noise1.wav", "noise.wav"]:
            candidate = d / nf
            if candidate.exists():
                noise_wav = candidate
                break

        if not meas_wavs or not target_wav.exists():
            return None

        preset_name = dir_name.replace("/", "-")
        return Preset(
            name=preset_name,
            description=f"Preset for {dir_name}",
            measurements=[str(w.relative_to(ROOT)) for w in meas_wavs],
            target=str(target_wav.relative_to(ROOT)),
            noise=str(noise_wav.relative_to(ROOT)) if noise_wav else None,
            fc=60.0,
            h2=0.5,
            h3=1.0,
        )
