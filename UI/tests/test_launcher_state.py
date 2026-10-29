"""Covers launcher.API's run-state machine — the sticky-until-acknowledged protocol the
React overlay depends on.

`app.jsx` polls `pipeline_status` about once a second and only leaves a terminal state when
the user dismisses the overlay, which calls `acknowledge_pipeline`. That stickiness is what
stops a fast run from finishing and resetting between two polls, so it is worth pinning:
the states are produced entirely here, and the UI merely renders them.

`launcher.py` imports pywebview at module scope, which the test env does not ship, so a stub
stands in — nothing under test touches it.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
sys.modules.setdefault("webview", types.SimpleNamespace(create_window=None, start=None))
sys.path.insert(0, str(ROOT / "UI"))

import launcher  # noqa: E402  (needs the stub above in place first)


class FakePopen:
    """Minimal Popen stand-in: `poll()` returns None until `finish()` sets a code."""

    def __init__(self, pid: int = 4242) -> None:
        self.pid = pid
        self.returncode = None
        self.signals: list[int] = []

    def poll(self):
        return self.returncode

    def finish(self, code: int = 0) -> None:
        self.returncode = code

    def send_signal(self, sig) -> None:
        self.signals.append(sig)


@pytest.fixture
def api():
    return launcher.API()


def _start(api, proc=None):
    proc = proc or FakePopen()
    api._set_run(proc)
    return proc


def test_starts_idle(api):
    assert api.pipeline_status() == {"state": "idle"}


def test_running_reports_pid_and_elapsed(api):
    proc = _start(api)
    status = api.pipeline_status()
    assert status["state"] == "running"
    assert status["pid"] == proc.pid
    assert status["elapsed"] >= 0


@pytest.mark.parametrize("code, expected", [(0, "succeeded"), (1, "failed"), (255, "failed")])
def test_exit_code_selects_terminal_state(api, code, expected):
    proc = _start(api)
    proc.finish(code)
    status = api.pipeline_status()
    assert status["state"] == expected
    assert status["returncode"] == code


def test_cancel_marks_the_exit_cancelled_not_failed(api):
    # SIGTERM makes snakemake exit non-zero, so without the `_cancelled` flag a user-driven
    # stop would be indistinguishable from a genuine failure.
    proc = _start(api)
    assert api.cancel_pipeline()["ok"] is True
    assert proc.signals, "cancel_pipeline should have signalled the process"
    proc.finish(-15)
    assert api.pipeline_status()["state"] == "cancelled"


def test_terminal_state_is_sticky_across_polls(api):
    # The property the 1 Hz poll relies on: a finished run keeps reporting its outcome
    # instead of falling back to idle.
    proc = _start(api)
    proc.finish(0)
    assert [api.pipeline_status()["state"] for _ in range(4)] == ["succeeded"] * 4


def test_elapsed_freezes_once_finished(api):
    proc = _start(api)
    proc.finish(0)
    first = api.pipeline_status()["elapsed"]
    assert api.pipeline_status()["elapsed"] == first


def test_acknowledge_returns_to_idle_only_when_stopped(api):
    proc = _start(api)
    assert api.acknowledge_pipeline() == {"ok": False, "error": "Pipeline still running"}
    proc.finish(0)
    assert api.acknowledge_pipeline() == {"ok": True}
    assert api.pipeline_status() == {"state": "idle"}


def test_cancel_without_a_run_is_refused(api):
    assert api.cancel_pipeline()["ok"] is False


def test_failed_run_carries_the_log_tail_and_success_does_not(api, monkeypatch):
    monkeypatch.setattr(launcher.API, "_read_log_tail", staticmethod(lambda **_: "boom"))
    proc = _start(api)
    proc.finish(2)
    assert api.pipeline_status()["log_tail"] == "boom"

    api2 = launcher.API()
    ok = _start(api2)
    ok.finish(0)
    assert "log_tail" not in api2.pipeline_status()


@pytest.fixture
def launchable(tmp_path, monkeypatch):
    # run_pipeline with the save, the log path and Popen stubbed; returns the procs launched.
    launched: list[FakePopen] = []

    def popen(*_a, stdout=None, **_k):
        launched.append(FakePopen(pid=len(launched) + 1))
        launched[-1].stdout = stdout
        return launched[-1]

    monkeypatch.setattr(launcher.API, "save_config", lambda self, cfg: {"ok": True})
    monkeypatch.setattr(launcher, "LOG_PATH", tmp_path / ".pipeline.log")
    monkeypatch.setattr(launcher.subprocess, "Popen", popen)
    return launched


def test_the_parent_holds_no_handle_on_the_run_log(api, launchable):
    # The child writes through its own copy of the fd; nothing is left for the parent to close.
    assert api.run_pipeline({})["ok"] is True
    assert launchable[0].stdout.closed


def test_relaunch_after_a_run_that_exited_unpolled(api, launchable):
    assert api.run_pipeline({})["ok"] is True
    launchable[0].finish(0)
    assert api.run_pipeline({})["ok"] is True
    assert api.pipeline_status()["pid"] == launchable[1].pid


def test_failed_relaunch_does_not_resurrect_the_previous_run(api, launchable, monkeypatch):
    api.run_pipeline({})
    launchable[0].finish(1)

    def missing(*_a, **_k):
        raise FileNotFoundError("snakemake")

    monkeypatch.setattr(launcher.subprocess, "Popen", missing)
    assert api.run_pipeline({})["ok"] is False
    assert api.pipeline_status() == {"state": "idle"}


def test_list_dir_paths_are_root_relative_posix(tmp_path, monkeypatch):
    # The `path` written here goes straight into config.yaml, so it must stay relative to
    # the repo root and POSIX-shaped whatever the platform.
    monkeypatch.setattr(launcher, "ROOT", tmp_path)
    directory = tmp_path / "Data" / "Raw" / "RNAs"
    directory.mkdir(parents=True)
    (directory / "b.fa").write_text(">b\n")
    (directory / "a.fa").write_text(">a\n")
    (directory / "sub").mkdir()

    assert launcher._list_dir(directory) == [
        {"name": "a.fa", "path": "Data/Raw/RNAs/a.fa"},
        {"name": "b.fa", "path": "Data/Raw/RNAs/b.fa"},
    ]


def test_list_dir_of_a_missing_directory_is_empty(tmp_path):
    assert launcher._list_dir(tmp_path / "nope") == []
