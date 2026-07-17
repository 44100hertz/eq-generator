#!/usr/bin/env python3
"""
Bridge evdev volume events → smart-volume FIFO.

Listens for KEY_VOLUMEUP / KEY_VOLUMEDOWN (or REL_DIAL / REL_HWHEEL
from a rotary encoder) and writes a 0-127 volume byte to
/tmp/eqgen_sv_fifo, driving the filter's smart-volume pipeline
identically to how AVRCP drives the ESP32 firmware.

Usage:
  # Auto-detect first suitable device (prefers rotary encoders):
  python tools/sv_knob.py

  # Specify a device explicitly:
  python tools/sv_knob.py --device /dev/input/by-id/usb-My_Knob

  # Adjust step size (default 2, range 1-16):
  python tools/sv_knob.py --step 4

  # Dry-run: print events without writing to FIFO:
  python tools/sv_knob.py --dry-run

Note: Disable your desktop's volume key shortcuts (e.g. KDE System Settings
→ Shortcuts → Volume Up/Down) so events reach this script.

Requires: python-evdev  (pip install evdev)
"""

import argparse
import errno
import os
import signal
import sys
import time

FIFO_PATH = "/tmp/eqgen_sv_fifo"

# Event codes we listen for
REL_DIAL = 0x07       # REL_DIAL — Griffin PowerMate, Surface Dial
REL_HWHEEL = 0x06     # REL_HWHEEL — horizontal scroll wheels
KEY_VOLUMEUP = 115
KEY_VOLUMEDOWN = 114


import subprocess as _sp

_NOTIFY_ID = None  # track notification to replace in-place
_LAST_NOTIFY = 0.0
_NOTIFY_COOLDOWN = 0.3  # seconds between notifications


def show_vol(vol: int, force: bool = False):
    """Show a transient volume OSD with progress bar via notify-send.
    Rate-limited to one per _NOTIFY_COOLDOWN seconds unless forced."""
    global _NOTIFY_ID, _LAST_NOTIFY
    now = time.monotonic()
    if not force and now - _LAST_NOTIFY < _NOTIFY_COOLDOWN:
        return
    _LAST_NOTIFY = now
    pct = round(vol / 127 * 100)
    bar = "█" * (pct // 5) + "░" * (20 - pct // 5)
    args = [
        "notify-send",
        "-a", "eqgen",
        "-h", f"int:value:{pct}",
        "-h", "string:synchronous:volume",
        "-t", "1500",
        "Volume",
        f"{bar}  {pct}%",
    ]
    _sp.Popen(args, stdout=_sp.DEVNULL, stderr=_sp.DEVNULL)
import json

def _read_knob_state(path):
    """Return (vol, knob_enabled) from state file, or defaults."""
    try:
        with open(path) as f:
            s = json.load(f)
            return s.get("vol", 63), s.get("knob_enabled", True)
    except (FileNotFoundError, json.JSONDecodeError):
        return 63, True


def _write_knob_state(path, vol, knob_enabled=True):
    """Write current vol and knob toggle to state file."""
    try:
        with open(path, "w") as f:
            json.dump({"vol": vol, "knob_enabled": knob_enabled}, f)
    except OSError:
        pass


def open_fifo():
    """Open the smart-volume FIFO for writing (non-blocking)."""
    try:
        fd = os.open(FIFO_PATH, os.O_WRONLY | os.O_NONBLOCK)
    except FileNotFoundError:
        sys.exit(
            f"ERROR: {FIFO_PATH} not found.\n"
            f"  Is the EQGen filter running?  (python -m eqgen.cli.wire setup ...)"
        )
    except PermissionError:
        sys.exit(f"ERROR: permission denied opening {FIFO_PATH}")
    return fd

def write_vol(fd, vol):
    """Write a single volume byte (0-127) to the FIFO."""
    vol = max(0, min(127, vol))
    try:
        os.write(fd, bytes([vol]))
    except OSError as e:
        if e.errno == errno.EAGAIN:
            pass  # pipe full, drop — next event will catch up
        else:
            raise


def find_devices(device_path=None):
    """Return a list of evdev device paths to listen on."""
    try:
        import evdev
    except ImportError:
        sys.exit("python-evdev not installed.  Run: pip install evdev")

    if device_path:
        return [evdev.InputDevice(device_path)]

    try:
        device_nodes = evdev.list_devices()
    except PermissionError:
        sys.exit(
            "Cannot list input devices — permission denied.\n"
            "  Add yourself to the 'input' group:  sudo usermod -a -G input $USER\n"
            "  Then log out and back in, or run:  newgrp input"
        )

    if not device_nodes:
        sys.exit(
            "No /dev/input/event* devices found.\n"
            "  (This can also mean you're not in the 'input' group —\n"
            "   run 'groups' to check, or try:  sudo usermod -a -G input $USER)"
        )

    candidates = []
    for dev in [evdev.InputDevice(p) for p in device_nodes]:
        caps = dev.capabilities()
        has_rel = caps.get(evdev.ecodes.EV_REL, [])
        has_key = caps.get(evdev.ecodes.EV_KEY, [])
        if REL_DIAL in has_rel or REL_HWHEEL in has_rel:
            candidates.append((0, dev))
        elif KEY_VOLUMEUP in has_key or KEY_VOLUMEDOWN in has_key:
            candidates.append((1, dev))

    candidates.sort(key=lambda x: x[0])
    if not candidates:
        sys.exit(
            "No volume-capable input device found.\n"
            "  (Checked for rotary encoders and media volume keys.)\n"
            "  Try:  python tools/sv_knob.py --list  to see all devices"
        )
    return [d for _, d in candidates]


def list_devices():
    """Print all evdev devices for diagnostics."""
    try:
        import evdev
    except ImportError:
        sys.exit("python-evdev not installed.  Run: pip install evdev")

    try:
        nodes = evdev.list_devices()
    except PermissionError:
        sys.exit("Permission denied.  Are you in the 'input' group?")

    if not nodes:
        print("No /dev/input/event* devices found.")
        return

    for p in nodes:
        d = evdev.InputDevice(p)
        caps = d.capabilities()
        rel = caps.get(evdev.ecodes.EV_REL, [])
        key = caps.get(evdev.ecodes.EV_KEY, [])
        has_dial = 0x07 in rel
        has_vol = 114 in key or 115 in key
        tag = ""
        if has_dial:
            tag = " [ROTARY]"
        elif has_vol:
            tag = " [VOLUME KEYS]"
        print(f"  {p}  {d.name!r}{tag}")


def main():
    ap = argparse.ArgumentParser(description="Bridge evdev → smart-volume FIFO")
    ap.add_argument("--list", "-l", action="store_true",
                    help="List all input devices and exit")
    ap.add_argument("--device", "-d", help="Input device path (auto-detect if omitted)")
    ap.add_argument("--step", "-s", type=int, default=2,
                    help="Volume steps per encoder tick (1-16, default 2)")
    ap.add_argument("--state-file",
                    default=os.path.expanduser("~/.config/eqgen/knob_state.json"),
                    help="Shared state file for server sync")
    ap.add_argument("--dry-run", "-n", action="store_true",
                    help="Print events, don't write to FIFO")
    args = ap.parse_args()

    if args.list:
        list_devices()
        return

    step = max(1, min(16, args.step))

    import evdev

    devices = find_devices(args.device)
    dev = devices[0]

    # Seed from state file, or default to 63
    vol, knob_enabled = _read_knob_state(args.state_file)
    if not knob_enabled:
        print(f"Knob disabled in state file — events will be ignored.  "
              f"Re-enable via web UI or by editing {args.state_file}")
    if not args.dry_run:
        fd = open_fifo()
        write_vol(fd, vol)   # sync initial state to filter
        _write_knob_state(args.state_file, vol, knob_enabled)
    else:
        fd = None

    print(f"Listening on {dev.path} ({dev.name}) — step={step}, vol={vol}")
    print(f"FIFO: {FIFO_PATH}" if not args.dry_run else "DRY RUN (no FIFO writes)")
    print("Press Ctrl+C to stop.")

    running = True

    def handler(sig, frame):
        nonlocal running
        running = False

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)

    try:
        for event in dev.read_loop():
            if not running:
                break

            if event.type == evdev.ecodes.EV_KEY:
                if event.value == 0:  # key release only
                    continue
                if event.code == KEY_VOLUMEUP:
                    vol = min(127, vol + step)
                elif event.code == KEY_VOLUMEDOWN:
                    vol = max(0, vol - step)
                else:
                    continue
            elif event.type == evdev.ecodes.EV_REL:
                delta = event.value
                if event.code == REL_DIAL:
                    vol = max(0, min(127, vol + delta * step))
                elif event.code == REL_HWHEEL:
                    vol = max(0, min(127, vol + delta * step))
                else:
                    continue
            else:
                continue

            # Re-read knob_enabled from state file (web UI can toggle it)
            _, knob_enabled = _read_knob_state(args.state_file)

            if not args.dry_run and fd is not None and knob_enabled:
                write_vol(fd, vol)
                _write_knob_state(args.state_file, vol, True)
                show_vol(vol)

            direction = "▲" if event.code == KEY_VOLUMEUP else \
                        "▼" if event.code == KEY_VOLUMEDOWN else \
                        "▲" if (event.type == evdev.ecodes.EV_REL and event.value > 0) else \
                        "▼" if (event.type == evdev.ecodes.EV_REL and event.value < 0) else "·"
            status = "" if knob_enabled else " [DISABLED]"
            print(f"vol={vol:3d}  {direction}{status}")

    except KeyboardInterrupt:
        pass
    finally:
        if fd is not None:
            os.close(fd)
        print("Stopped.")

if __name__ == "__main__":
    main()
