"""PyWebView launcher for the Novum Pipeline config editor.

Loads ``UI/Novum Pipeline.html`` in a native window and bridges it to the local
filesystem via an ``API`` instance exposed as ``window.pywebview.api`` on the JS
side (file pickers, config save, and the pipeline run/poll/cancel lifecycle).
See the README for env setup; run with ``python UI/launcher.py``.

``run_pipeline`` records the live ``Popen`` so the UI can poll ``pipeline_status``
and call ``cancel_pipeline``; Snakemake's ``.snakemake/locks/`` still guards against
overlapping runs. Subprocess output is written to ``Data/Results/.pipeline.log``,
which ``pipeline_status`` tails back into the UI on a non-clean exit.
"""
from __future__ import annotations

import os
import signal
import subprocess
import time
from pathlib import Path
from threading import Lock, Timer

# Qt high-DPI: must be set before pywebview imports PyQt5/QtWebEngine.
# QtWebEngine on WSLg ignores --force-device-scale-factor, so page scaling is
# handled in CSS (`html { zoom: ... }` in UI/src/tokens.css) instead.
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")
# WSLg has no GPU passthrough; skip Chromium's GPU init to silence the
# transient `GpuChannelMsg_CreateCommandBuffer` error before fallback.
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu --log-level=3")

import webview

from config_io import load_config as _load_config, save_config as _save_config


ROOT = Path(__file__).resolve().parent.parent
QUERIES_DIR = ROOT / "Data" / "Raw" / "RNAs"
TARGETS_DIR = ROOT / "Data" / "Raw" / "genomes"
CONFIG_PATH = ROOT / "Config" / "config.yaml"
HTML_PATH = ROOT / "UI" / "Novum Pipeline.html"
LOG_PATH = ROOT / "Data" / "Results" / ".pipeline.log"


#----- Lists files in `directory` as {name, path} dicts the React file pickers consume. The
#      path is spelled relative to ROOT, in POSIX form, because it is what lands in
#      config.yaml — deriving it here keeps it from drifting from `directory` -----#
def _list_dir(directory: Path) -> list[dict]:
    if not directory.is_dir():
        return []
    prefix = directory.relative_to(ROOT).as_posix()
    return [{"name": entry.name, "path": f"{prefix}/{entry.name}"}
            for entry in sorted(directory.iterdir(), key=lambda p: p.name)
            if entry.is_file()]


#----- The object exposed to JS as `window.pywebview.api`. Every public method here is
#      callable from the React side; `_lock` serialises the run-state mutations, which the
#      1 Hz status poll and the user's Cancel can otherwise reach concurrently -----#
class API:
    def __init__(self) -> None:
        self._lock = Lock()
        self._set_run(None)

    #----- The whole run state, set in one place: idle when `proc` is None, otherwise a fresh run
    #      with the clock started. `_finished_at` doubles as the "already harvested" latch -----#
    def _set_run(self, proc: "subprocess.Popen | None") -> None:
        self._proc = proc
        self._started_at: float | None = time.time() if proc else None
        self._finished_at: float | None = None
        self._cancelled = False

    #----- file pickers used by the React form -----#

    def list_queries(self) -> list[dict]:
        return _list_dir(QUERIES_DIR)

    def list_targets(self) -> list[dict]:
        return _list_dir(TARGETS_DIR)

    #----- config save -----#

    def load_config(self) -> dict:
        return _load_config(CONFIG_PATH)

    #----- Write config.yaml, turning an OS error into a result the UI can show -----#
    def save_config(self, config: dict) -> dict:
        try:
            return _save_config(CONFIG_PATH, config)
        except OSError as exc:
            return {"ok": False, "error": str(exc)}

    #----- pipeline lifecycle -----#

    #----- A run is live only while a process exists AND has not exited -----#
    def _is_running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    #----- Last `max_lines` of the run log, so a failure can be read in-app instead of
    #      sending the user back to the launcher terminal. Only the final `max_bytes` are
    #      read, so a long run's log is never slurped whole; the leading partial line is
    #      dropped by the line slice -----#
    @staticmethod
    def _read_log_tail(max_lines: int = 80, max_bytes: int = 65536) -> str:
        try:
            with LOG_PATH.open("rb") as fh:
                fh.seek(0, os.SEEK_END)
                fh.seek(max(0, fh.tell() - max_bytes))
                tail = fh.read().decode("utf-8", "replace")
        except OSError:
            return ""
        return "\n".join(tail.splitlines()[-max_lines:])

    #----- Saves the config, then launches snakemake against it -----#
    def run_pipeline(self, config: dict) -> dict:
        with self._lock:
            if self._is_running():
                return {
                    "ok": False,
                    "error": "Pipeline already running",
                    "pid": self._proc.pid,
                }
            save_result = self.save_config(config)
            if not save_result["ok"]:
                return save_result
            # Harvest a previous run that exited unpolled, so a failed launch reads as idle.
            self._set_run(None)
            try:
                LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
                # The child gets its own copy of the fd, so the parent's closes at launch.
                with LOG_PATH.open("w") as log:
                    proc = subprocess.Popen(
                        ["snakemake", "--use-conda", "--cores", "all"],
                        cwd=str(ROOT),
                        stdout=log,
                        stderr=subprocess.STDOUT,
                    )
            except FileNotFoundError as exc:
                return {"ok": False, "error": f"snakemake not found on PATH ({exc})"}
            self._set_run(proc)
            return {
                "ok": True,
                "pid": proc.pid,
                "log_path": str(LOG_PATH),
            }

    #----- Polled ~1/s by the UI while a run is alive; terminal states stay sticky until
    #      acknowledge_pipeline, so a fast finish cannot slip between two polls -----#
    def pipeline_status(self) -> dict:
        with self._lock:
            if self._proc is None:
                return {"state": "idle"}

            if self._proc.poll() is None:
                return {
                    "state": "running",
                    "pid": self._proc.pid,
                    "elapsed": time.time() - self._started_at,
                }

            if self._finished_at is None:          # first poll after the process exited
                self._finished_at = time.time()

            code = self._proc.returncode
            state = "cancelled" if self._cancelled else "succeeded" if code == 0 else "failed"
            result = {
                "state": state,
                "pid": self._proc.pid,
                "returncode": code,
                "elapsed": self._finished_at - self._started_at,
                "log_path": str(LOG_PATH),
            }
            # Surface the trace in-app for non-clean exits; success stays terse.
            if state in ("failed", "cancelled"):
                result["log_tail"] = self._read_log_tail()
            return result

    #----- SIGTERM the run; `_cancelled` is what makes the exit read as "cancelled" rather
    #      than "failed" when the status is next polled -----#
    def cancel_pipeline(self) -> dict:
        with self._lock:
            if not self._is_running():
                return {"ok": False, "error": "No pipeline running"}
            self._cancelled = True
            try:
                self._proc.send_signal(signal.SIGTERM)
            except OSError as exc:
                return {"ok": False, "error": str(exc)}
            return {"ok": True, "pid": self._proc.pid}

    #----- Drops the sticky terminal status so the UI returns to idle -----#
    def acknowledge_pipeline(self) -> dict:
        with self._lock:
            if self._is_running():
                return {"ok": False, "error": "Pipeline still running"}
            self._set_run(None)
            return {"ok": True}


#----- Builds the PyWebView window, wires the min-width snap-back, and starts the Qt loop -----#
def main() -> None:
    # Width floor only: the cards break at ~768 CSS px, so this is 768 x --ui-zoom
    # (tokens.css). Retune it whenever --ui-zoom changes.
    min_w = 960
    window = webview.create_window(
        "Novum Pipeline",
        str(HTML_PATH),
        js_api=API(),
        width=1280,
        height=820,
        min_size=(min_w, 1),
    )

    #----- Snaps sub-floor widths back to min_w, debounced: resizing mid-drag makes the
    #      window drift in the pull direction -----#
    timer: Timer | None = None

    def _snap_back(width, _height, *_args) -> None:
        nonlocal timer
        # Cancel any pending snap: the drag may have crossed back above the floor.
        if timer is not None:
            timer.cancel()
            timer = None
        if width >= min_w:
            return

        def _fire() -> None:
            # Read the live size at fire time so we restore the user's current
            # height, not the (possibly stale) height when the timer was armed.
            try:
                w = window.width
                h = window.height
            except Exception as exc:
                print(f"[launcher] min-width snap-back: live size read failed: {exc}")
                return
            if w >= min_w:
                return
            try:
                window.resize(min_w, h)
            except Exception as exc:
                print(f"[launcher] min-width snap-back failed: {exc}")

        timer = Timer(0.12, _fire)
        timer.daemon = True
        timer.start()

    window.events.resized += _snap_back

    # Qt, because GTK needs a system PyGObject conda envs can't see. private_mode=False keeps
    # localStorage (the theme toggle) and the code caches across launches. debug=True opens
    # DevTools.
    webview.start(gui="qt", debug=False, private_mode=False)


if __name__ == "__main__":
    main()
