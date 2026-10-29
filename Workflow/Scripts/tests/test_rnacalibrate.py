from __future__ import annotations
import json
import math
import os
import subprocess
import sys
from functools import partial
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import pytest
import _gumbel_fit as gf
import rnacalibrate as rc

SCRIPTS = str(Path(__file__).resolve().parent.parent)


def _fasta(tmp_path, text):
    p = tmp_path / "t.fna"
    p.write_text(text)
    return p


def test_length_stats_population_stdev(tmp_path):
    f = _fasta(tmp_path, ">a\nAAAA\n>b\nAAAAAA\n")
    stats = rc.summarise_lengths(rc.record_lengths(f))
    assert stats == {"count": 2, "mean": 5.0, "std": 1.0}


def test_length_stats_single_record_has_zero_stdev(tmp_path):
    f = _fasta(tmp_path, ">a\nAAAA\n")
    assert rc.summarise_lengths(rc.record_lengths(f))["std"] == 0.0


def test_length_stats_rejects_empty_file(tmp_path):
    with pytest.raises(RuntimeError, match="No FASTA records"):
        rc.record_lengths(_fasta(tmp_path, ""))


def test_build_length_arg_rounds_to_ints():
    assert rc.build_length_arg({"count": 3, "mean": 986.4, "std": 792.2}) == {
        "mean": 986, "std": 792, "value": "986,792",
    }


LENGTH_ARG = {"mean": 986, "std": 792, "value": "986,792"}


def test_build_command_minimal_omits_optional_flags():
    assert rc.build_command("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG) == [
        "RNAcalibrate", "-k", "10000", "-q", "q.fa", "-t", "t.fna",
        "-m", "50000", "-l", "986,792",
    ]


def test_build_command_appends_every_optional_flag_in_order():
    assert rc.build_command(
        "RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG,
        randomize_targets=True, max_internal_loop=4, max_bulge_loop=5, seed="2,7",
    ) == [
        "RNAcalibrate", "-k", "10000", "-q", "q.fa", "-t", "t.fna",
        "-m", "50000", "-l", "986,792", "-u", "4", "-v", "5", "-f", "2,7", "-s",
    ]


def test_build_command_seed_is_forced_helix_not_rng():
    # Guards the naming trap: `seed` here means RNAhybrid's -f, never an RNG seed.
    cmd = rc.build_command("RNAcalibrate", "q.fa", "t.fna", 1, 1, LENGTH_ARG, seed="2,7")
    assert "-f" in cmd and "2,7" in cmd


def test_parse_output_returns_one_entry_per_row():
    parsed = rc.parse_rnacalibrate_output(
        "hsa-let-7a-5p 64 2.118843 0.166056\nhsa-miR-21-5p 63 2.025677 0.163574\n"
    )
    assert parsed == {"per_query": [
        {"query": "hsa-let-7a-5p", "sample_size": 64, "xi": 2.118843, "theta": 0.166056},
        {"query": "hsa-miR-21-5p", "sample_size": 63, "xi": 2.025677, "theta": 0.163574},
    ]}


def test_parse_output_rejects_wrong_column_count():
    with pytest.raises(RuntimeError, match="expected 4 columns"):
        rc.parse_rnacalibrate_output("hsa-let-7a-5p 64 2.118843\n")


def test_parse_output_rejects_empty():
    with pytest.raises(RuntimeError, match="did not produce any calibration rows"):
        rc.parse_rnacalibrate_output("\n\n")


def test_derive_seed_is_stable_against_known_values():
    # Hardcoded on purpose: these must not drift between runs, machines, or releases.
    assert rc.derive_seed(1, "UGAGGUAGUAGGUUGUAUAGUU") == 2093734617
    assert rc.derive_seed(1, "GGGGCCUGGCGGUGGGCGG") == 1920725615
    assert rc.derive_seed(1, "ACGU") == 600126125


def test_derive_seed_differs_per_sequence_and_per_base():
    assert rc.derive_seed(1, "ACGU") != rc.derive_seed(1, "ACGA")
    assert rc.derive_seed(1, "ACGU") != rc.derive_seed(2, "ACGU")


def test_derive_seed_stays_in_faketime_safe_range():
    for seq in ("ACGU", "UGAGGUAGUAGGUUGUAUAGUU", "G" * 200):
        value = rc.derive_seed(2**31 - 1, seq)
        assert 0 <= value < 2**31


def test_derive_seed_ignores_pythonhashseed():
    # The regression guard for using built-in hash(), which is salted per process.
    code = (
        f"import sys; sys.path.insert(0, {SCRIPTS!r}); "
        "import rnacalibrate; print(rnacalibrate.derive_seed(1, 'ACGU'))"
    )
    outputs = set()
    for hashseed in ("0", "1", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                              text=True, env=env, check=True)
        outputs.add(done.stdout.strip())
    assert outputs == {"600126125"}


def test_faketime_env_uses_ld_preload_on_linux():
    env = rc.build_faketime_env(42, "/opt/lib/libfaketime.so.1", "linux", {"PATH": "/bin"})
    assert env["LD_PRELOAD"] == "/opt/lib/libfaketime.so.1"
    assert env["FAKETIME_FMT"] == "%s"
    assert env["FAKETIME"] == "@42"
    assert env["PATH"] == "/bin"           # base env is preserved
    assert "DYLD_INSERT_LIBRARIES" not in env


def test_faketime_env_uses_dyld_on_darwin():
    env = rc.build_faketime_env(42, "/opt/lib/libfaketime.1.dylib", "darwin", {})
    assert env["DYLD_INSERT_LIBRARIES"] == "/opt/lib/libfaketime.1.dylib"
    assert env["DYLD_FORCE_FLAT_NAMESPACE"] == "1"
    assert env["FAKETIME"] == "@42"
    assert "LD_PRELOAD" not in env


def test_faketime_env_does_not_mutate_the_base_env():
    base = {"PATH": "/bin"}
    rc.build_faketime_env(1, "/x.so", "linux", base)
    assert base == {"PATH": "/bin"}


def test_faketime_env_rejects_unknown_platform():
    with pytest.raises(RuntimeError, match="Unsupported platform"):
        rc.build_faketime_env(1, "/x.so", "win32", {})


def test_faketime_env_seed_is_an_absolute_frozen_timestamp():
    # The '@' prefix freezes the clock; without it libfaketime treats the value as an offset.
    assert rc.build_faketime_env(0, "/x.so", "linux", {})["FAKETIME"] == "@0"


def test_resolve_faketime_library_reports_the_env_file_when_missing(monkeypatch):
    monkeypatch.setattr(rc, "which_required",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("not found")))
    with pytest.raises(RuntimeError, match="rnahybrid.yaml"):
        rc.resolve_faketime_library("linux")


def test_resolve_faketime_library_finds_the_object_beside_the_wrapper(tmp_path, monkeypatch):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "faketime").write_text("")
    (tmp_path / "lib" / "faketime").mkdir(parents=True)
    (tmp_path / "lib" / "faketime" / "libfaketime.so.1").write_text("")
    monkeypatch.setattr(rc, "which_required", lambda *a: str(tmp_path / "bin" / "faketime"))
    found = rc.resolve_faketime_library("linux")
    # .resolve() on both sides: macOS symlinks /tmp -> /private/tmp.
    assert found == (tmp_path / "lib" / "faketime" / "libfaketime.so.1").resolve()


def test_resolve_faketime_library_errors_when_wrapper_exists_but_object_does_not(tmp_path, monkeypatch):
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "faketime").write_text("")
    monkeypatch.setattr(rc, "which_required", lambda *a: str(tmp_path / "bin" / "faketime"))
    with pytest.raises(RuntimeError, match="libfaketime.so.1"):
        rc.resolve_faketime_library("linux")


@pytest.mark.parametrize("first, repeat, other, verdict", [
    ("A", "A", "B", "verified"),
    ("A", "B", "C", "blocked"),
    # Two identical runs alone would look like success; the third run is what proves
    # the seed actually controls the output.
    ("A", "A", "A", "blind"),
    ("A", "B", "A", "blocked"),     # blocked outranks blind
])
def test_classify_probe(first, repeat, other, verdict):
    assert rc.classify_probe(first, repeat, other) == verdict


@pytest.mark.parametrize("rng_seed, randomize_targets, expected", [
    (1, True, True),
    (None, True, False),
    # `-s` is the only consumer of the RNG. Without it RNAcalibrate is already
    # byte-reproducible, so pinning is inert — and probing would report a false "blind"
    # and fail a valid config.
    (1, False, False),
    (None, False, False),
])
def test_pinning_enabled(rng_seed, randomize_targets, expected):
    assert rc.pinning_enabled(rng_seed, randomize_targets) is expected


def test_query_key_takes_the_first_whitespace_token():
    assert rc.query_key("hsa-miR-21-5p") == "hsa-miR-21-5p"
    assert rc.query_key("hsa-miR-21-5p some description") == "hsa-miR-21-5p"


def test_iter_query_records_yields_header_and_sequence(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">hsa-let-7a-5p\nUGAGGUAG\nUAGGUU\n>hsa-miR-21-5p\nUAGCUUAU\n")
    assert list(rc.iter_query_records(q)) == [
        ("hsa-let-7a-5p", "UGAGGUAGUAGGUU"),
        ("hsa-miR-21-5p", "UAGCUUAU"),
    ]


def test_iter_query_records_rejects_duplicate_ids(tmp_path):
    # Today RNAcalibrate emits two rows and load_per_query_distributions silently
    # keeps the last. split_query_per_miRNA already guards this; match it.
    q = tmp_path / "q.fa"
    q.write_text(">hsa-let-7a-5p\nACGU\n>hsa-let-7a-5p\nUGCA\n")
    with pytest.raises(RuntimeError, match="Duplicate query ID"):
        list(rc.iter_query_records(q))


def test_iter_query_records_rejects_an_empty_file(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text("")
    with pytest.raises(RuntimeError, match="No FASTA records"):
        list(rc.iter_query_records(q))


def test_expect_single_row_returns_the_entry():
    parsed = {"per_query": [{"query": "m", "sample_size": 1, "xi": 0.1, "theta": 0.2}]}
    assert rc.expect_single_row(parsed, "m")["xi"] == 0.1


def test_expect_single_row_rejects_multiple_rows():
    parsed = {"per_query": [{"query": "a"}, {"query": "b"}]}
    with pytest.raises(RuntimeError, match="exactly 1 calibration row"):
        rc.expect_single_row(parsed, "a")


def test_sha256_file_matches_a_known_digest(tmp_path):
    f = tmp_path / "x.txt"
    f.write_bytes(b"ACGU\n")
    assert rc.sha256_file(f) == (
        "d57db2001ce2b4da9b363b703c19d7f7a94772c0c895388d06812de225366798"
    )


@pytest.mark.parametrize("variant", ["acgu", "ACGT", "acgt", " ACGU \n", "\tACGU"])
def test_derive_seed_and_sequence_sha256_share_one_normalisation(variant):
    # The contract is that both consumers hash the same string. Asserted without restating
    # derive_seed's mixing arithmetic, which would only ever agree with itself.
    assert rc.derive_seed(1, variant) == rc.derive_seed(1, "ACGU")
    assert rc.sequence_sha256(variant) == rc.sequence_sha256("ACGU")


def test_count_fasta_records(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">a\nACGU\n>b\nUGCA\n>c\nAAAA\n")
    assert rc.count_fasta_records(f) == 3


@pytest.mark.parametrize("rng_seed", [1, None])
def test_build_inputs_block_shape(tmp_path, rng_seed):
    q = tmp_path / "q.fa"
    q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"
    t.write_text(">g\nAAAA\n>h\nCCCC\n")
    block = rc.build_inputs_block(
        query=q, target=t, k=10000, max_target_length=50000,
        length_arg={"mean": 986, "std": 792, "value": "986,792"},
        forced_helix=None, max_internal_loop=None, max_bulge_loop=None, rng_seed=rng_seed,
    )
    assert block["query"]["n_records"] == 1
    assert block["target"]["n_records"] == 2
    assert block["query"]["sha256"] == rc.sha256_file(q)
    assert block["target"]["sha256"] == rc.sha256_file(t)
    assert block["params"] == {
        "k": 10000, "max_target_length": 50000, "length_arg": "986,792",
        "length_anchors": None,
        "forced_helix": None, "max_internal_loop": None, "max_bulge_loop": None,
    }
    assert block["rng_seed"] == rng_seed


# ---- verify_faketime: probe sequencing and the sleep between the two same-seed probes ----
# subprocess.run and time.sleep are faked, so these run instantly despite the real 1.1s sleep.


def _completed(stdout):
    return SimpleNamespace(stdout=stdout)


def test_verify_faketime_sleeps_between_first_and_second_probe_then_verifies(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        calls.append(("run", kwargs["env"]["FAKETIME"]))
        return _completed(kwargs["env"]["FAKETIME"])

    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: calls.append(("sleep", seconds)))

    rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG, 1,
                       "/opt/lib/libfaketime.so.1", "linux")

    # Same seed twice (A, B) around the sleep, then a different seed (C) with no sleep
    # before it — matching classify_probe's A==B, A!=C contract.
    assert [c[0] for c in calls] == ["run", "sleep", "run", "run"]
    assert calls[1][1] >= 1.0


@pytest.mark.parametrize("k, expected", [
    (10000, [5, 200, 10000]),   # the shipped case: both cheap rungs, then the real k
    (200, [5, 200]),            # k equal to a rung is not probed twice
    (100, [5, 100]),            # a rung at or above k is dropped, never probed above the run
    (5, [5]),
    (2, [2]),                   # k below every rung leaves only k itself
])
def test_probe_ladder_tops_out_at_the_configured_k(k, expected):
    # A rung below the fit's convergence threshold returns a constant `-nan`, which reads as
    # "the seed has no effect"; the run's own k is the only rung guaranteed to be in the
    # regime the calibration itself uses, so it is always last and always present.
    assert rc.probe_ladder(k) == expected


@pytest.mark.parametrize("blind_at, rungs", [
    ([rc.PROBE_K], 2),                      # rescued by k=200, as 27 of 43 samples are
    ([rc.PROBE_K, rc.PROBE_K_RETRY], 3),    # only the configured k can settle it
])
def test_verify_faketime_escalates_until_a_rung_discriminates(monkeypatch, blind_at, rungs):
    calls = []

    def fake_run(command, **kwargs):
        probe_k = int(command[2])
        calls.append(("run", probe_k))
        if probe_k in blind_at:
            return _completed("indistinguishable")
        return _completed(kwargs["env"]["FAKETIME"])

    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: calls.append(("sleep", seconds)))

    rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG, 1,
                       "/opt/lib/libfaketime.so.1", "linux")

    # 3 probes + 1 sleep per rung, and the last rung is the one that discriminated.
    assert [c[0] for c in calls] == ["run", "sleep", "run", "run"] * rungs
    assert calls[-1][1] == rc.probe_ladder(10000)[rungs - 1]


def test_verify_faketime_blind_at_every_rung_names_the_configured_k(monkeypatch):
    # The message must not name PROBE_K_RETRY: blind at k=200 is routine and recoverable,
    # while blind at the run's own k is the verdict that actually refuses to proceed.
    monkeypatch.setattr(rc.subprocess, "run",
                        lambda command, **kwargs: _completed("indistinguishable"))
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: None)

    with pytest.raises(RuntimeError) as excinfo:
        rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG, 1,
                           "/opt/lib/libfaketime.so.1", "linux")

    assert "k=10000" in str(excinfo.value)


# ---- error messages carry the full command, not just command[:2] ----

def test_verify_faketime_probe_error_includes_the_full_command(monkeypatch):
    def fake_run(command, **kwargs):
        raise subprocess.CalledProcessError(2, command, output="", stderr="blocked")

    monkeypatch.setattr(rc.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as excinfo:
        rc.verify_faketime("RNAcalibrate", "q.fa", "t.fna", 10000, 50000, LENGTH_ARG, 1,
                           "/opt/lib/libfaketime.so.1", "linux")

    message = str(excinfo.value)
    # The old bug rendered only command[:2] ("RNAcalibrate -k"); -t/-m only appear once
    # the full command is joined.
    assert "-t t.fna" in message
    assert "-m 50000" in message


def test_run_rnacalibrate_error_includes_the_full_command(tmp_path, monkeypatch):
    q = tmp_path / "q.fa"; q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"; t.write_text(">g\nAAAA\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")

    def fake_run(command, **kwargs):
        raise subprocess.CalledProcessError(1, command, output="", stderr="boom")

    monkeypatch.setattr(rc.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError) as excinfo:
        rc.run_rnacalibrate(q, t, tmp_path / "out.json", k=10, max_target_length=100)

    message = str(excinfo.value)
    assert "RNAcalibrate failed for miRNA m" in message
    assert "-t" in message and "-m 100" in message


@pytest.mark.parametrize("randomize_targets, rng_seed, warns", [
    (False, 42, True),          # set but inert: `-s` is never passed, so nothing to pin
    (False, None, False),
    (True, 42, False),          # the normal pinned case
])
def test_run_rnacalibrate_warns_only_when_the_seed_is_inert(
        tmp_path, monkeypatch, capsys, randomize_targets, rng_seed, warns):
    q = tmp_path / "q.fa"; q.write_text(">m\nACGU\n")
    t = tmp_path / "t.fna"; t.write_text(">g\nAAAA\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")
    monkeypatch.setattr(rc, "resolve_faketime_library", lambda *a: "/opt/lib/libfaketime.so.1")
    monkeypatch.setattr(rc, "verify_faketime", lambda *a, **k: None)  # probed above
    monkeypatch.setattr(rc.subprocess, "run", lambda *a, **k: _completed("m 300 2.4 0.18\n"))

    rc.run_rnacalibrate(q, t, tmp_path / "out.json", k=10, max_target_length=100,
                        randomize_targets=randomize_targets, rng_seed=rng_seed)

    err = capsys.readouterr().err
    if warns:
        assert f"rng_seed={rng_seed}" in err and "randomize_targets" in err
    else:
        assert err == ""


# --- the length-anchor ladder ---


def test_derive_seed_is_unchanged_when_no_anchor_is_given():
    # Pinned against the shipped code: the reference fit must keep its seed so per_query
    # stays comparable with the 96 runs on disk.
    assert rc.derive_seed(1, "GUGAGGACUCGGGAGGUGG") == 1932658044
    assert rc.derive_seed(0, "GUGAGGACUCGGGAGGUGG") == 1932658043
    assert rc.derive_seed(1, "ACGT") == 600126125

def test_derive_seed_still_normalises_t_to_u_and_case():
    assert rc.derive_seed(1, "ACGT") == rc.derive_seed(1, "acgu")
    assert rc.derive_seed(1, "ACGT", anchor=900) == rc.derive_seed(1, "acgu", anchor=900)

def test_derive_seed_differs_per_anchor_and_per_attempt():
    base = rc.derive_seed(1, "ACGU")
    seeds = {rc.derive_seed(1, "ACGU", anchor=a) for a in (76, 150, 300, 600, 900, 1500)}
    assert len(seeds) == 6
    assert base not in seeds
    attempts = {rc.derive_seed(1, "ACGU", anchor=900, attempt=i) for i in range(5)}
    assert len(attempts) == 5

def test_derive_seed_stays_inside_the_faketime_safe_range():
    for anchor in (76, 1500):
        for attempt in range(5):
            seed = rc.derive_seed(1, "ACGU", anchor=anchor, attempt=attempt)
            assert 0 <= seed < 2 ** 31

def test_derive_seed_is_immune_to_pythonhashseed(tmp_path):
    # sha256, never the salted built-in hash() — two interpreters must agree. Pinned
    # against the shipped code, like its pre-existing sibling above.
    code = ("import sys; sys.path.insert(0, %r); import rnacalibrate as rc; "
            "print(rc.derive_seed(1, 'ACGU', anchor=900, attempt=2))" % SCRIPTS)
    outs = set()
    for hashseed in ("0", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        outs.add(subprocess.run([sys.executable, "-c", code], capture_output=True,
                                text=True, env=env, check=True).stdout.strip())
    assert outs == {"1984738640"}


def test_anchor_length_arg_fixes_the_null_width_at_one_third():
    assert rc.anchor_length_arg(900) == {"mean": 900, "std": 300, "value": "900,300"}
    assert rc.anchor_length_arg(76) == {"mean": 76, "std": 25, "value": "76,25"}

def test_anchor_length_arg_widens_for_the_tier_two_repair():
    assert rc.anchor_length_arg(900, divisor=2)["value"] == "900,450"


def test_repair_ladder_is_two_redraws_then_widen_then_raise_k():
    assert rc.REPAIR_LADDER == (
        (3, 1, "fitted"),
        (3, 1, "fitted"),
        (3, 1, "fitted"),
        (2, 1, "refitted_widened"),
        (3, 3, "refitted_high_k"),
    )


def test_count_records_per_anchor_bins_by_the_geometric_midpoint():
    # 50 -> 76 (open bottom cell), 734 -> 600, 736 -> 900, 5000 -> 3000 (open top cell)
    counts = rc.count_records_per_anchor([50, 734, 736, 5000],
                                         (76, 150, 300, 600, 900, 1500, 3000))
    assert counts == {76: 1, 150: 0, 300: 0, 600: 1, 900: 1, 1500: 0, 3000: 1}


def test_record_lengths_feeds_both_the_l_stats_and_the_cell_counts(tmp_path):
    # The two derivations must describe the same population — one pass, one `if seq` filter.
    f = _fasta(tmp_path, "".join(
        f">r{i}\n{'A' * length}\n" for i, length in enumerate([50, 734, 736, 5000])))
    lengths = rc.record_lengths(f)
    assert lengths == [50, 734, 736, 5000]


def test_parse_output_lets_nan_through():
    # Parsing never judges a fit; classify_fit does, on both the reference and anchor paths.
    parsed = rc.parse_rnacalibrate_output("hsa-miR-1224-5p 222 -nan -nan\n")
    row = parsed["per_query"][0]
    assert row["sample_size"] == 222
    assert math.isnan(row["xi"]) and math.isnan(row["theta"])


# --- resolve_strata: Tier 3 and extrapolation ---

REFERENCE = {"query": "m", "sample_size": 240, "xi": 2.40, "theta": 0.19,
             "derived_seed": 7}

def _entry(anchor, theta, sample_size=250, status="fitted", xi=2.5):
    return {"anchor": anchor, "length_arg": f"{anchor},{anchor // 3}", "k": 10000,
            "sample_size": sample_size, "xi": xi, "theta": theta,
            "derived_seed": 1, "attempts": 1, "status": status,
            "reason": None, "residual": None}

def _rejected(anchor, reason="non_finite"):
    return {"anchor": anchor, "length_arg": None, "k": None, "sample_size": None,
            "xi": None, "theta": None, "derived_seed": None, "attempts": 5,
            "status": "rejected", "reason": reason, "residual": None}

def _on_curve(anchor, alpha=1.557, log_c=0.5, query_nt=22):
    return gf.predict_theta(alpha, log_c, anchor, query_nt)

def test_resolve_strata_keeps_sound_fits_and_records_their_residual():
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900, 1500)]
    strata, alpha, spread = rc.resolve_strata(entries, 22, REFERENCE)
    assert alpha == pytest.approx(1.557, abs=1e-9)
    assert spread == pytest.approx(0.0, abs=1e-9)
    assert [s["status"] for s in strata] == ["fitted"] * 6
    assert all(s["residual"] == pytest.approx(1.0) for s in strata)

def test_resolve_strata_extrapolates_a_rejected_cell_off_the_curve():
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900)]
    entries.append(_rejected(1500))
    strata, alpha, _ = rc.resolve_strata(entries, 22, REFERENCE)
    top = strata[-1]
    assert top["status"] == "extrapolated"
    assert top["reason"] == "non_finite"
    assert top["theta"] == pytest.approx(_on_curve(1500))
    assert top["xi"] == 2.5              # weighted median of the survivors' xi
    assert top["length_arg"] is None     # nothing was invoked for this cell

def test_extrapolated_xi_is_the_survivors_weighted_median():
    # The fixture above gives every anchor the same xi, so it cannot tell the weighted
    # median from any other summary. Here the five xi are distinct and anchor 600 carries
    # 40x the weight: 2.60 is its value, and no unweighted median (2.65), mean (2.69),
    # weighted mean (2.610), first survivor (2.90) or last (2.50) lands on it.
    xis = {76: 2.90, 150: 2.80, 300: 2.65, 600: 2.60, 900: 2.50}
    entries = [_entry(a, _on_curve(a), sample_size=400 if a == 600 else 10, xi=xis[a])
               for a in (76, 150, 300, 600, 900)]
    entries.append(_rejected(1500))
    strata, _, _ = rc.resolve_strata(entries, 22, REFERENCE)
    assert strata[-1]["status"] == "extrapolated"
    assert strata[-1]["xi"] == 2.60

REPAIRED_STATUSES = ("refitted_widened", "refitted_high_k")

def test_fitted_statuses_covers_every_repair_rung():
    # Spelled out rather than derived from FITTED_STATUSES: hand-listing it as ("fitted",)
    # would otherwise empty the parametrize below and turn that test into a silent skip.
    assert set(rc.FITTED_STATUSES) == {"fitted", *REPAIRED_STATUSES}

@pytest.mark.parametrize("status", REPAIRED_STATUSES)
def test_resolve_strata_fits_on_every_repaired_status_not_just_fitted(status):
    # FITTED_STATUSES is derived from REPAIR_LADDER so a renamed or added rung stays in the
    # least squares. Hand-listing ("fitted",) would silently re-label a REPAIRED cell as
    # extrapolated and drop its measured theta from the very curve it would be read off.
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900)]
    entries.append(_entry(1500, _on_curve(1500), status=status))
    strata, _, _ = rc.resolve_strata(entries, 22, REFERENCE)
    top = strata[-1]
    assert top["status"] == status
    assert top["residual"] == pytest.approx(1.0)
    assert top["length_arg"] is not None  # still the measured fit, not a modelled cell

def test_resolve_strata_extrapolates_an_above_ceiling_anchor():
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900, 1500)]
    entries.append(_rejected(3000, reason="above_ceiling"))
    strata, _, _ = rc.resolve_strata(entries, 22, REFERENCE)
    assert strata[-1]["anchor"] == 3000
    assert strata[-1]["status"] == "extrapolated"
    assert strata[-1]["reason"] == "above_ceiling"
    assert strata[-1]["theta"] == pytest.approx(_on_curve(3000))

def test_resolve_strata_extrapolates_a_three_sigma_outlier():
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900, 1500)]
    entries[3]["theta"] *= 1.5
    strata, alpha, _ = rc.resolve_strata(entries, 22, REFERENCE)
    outlier = next(s for s in strata if s["anchor"] == 600)
    assert outlier["status"] == "extrapolated"
    assert outlier["reason"] == "residual_outlier"
    assert alpha == pytest.approx(1.557, abs=1e-9)   # refitted without it

def test_resolve_strata_falls_back_to_the_reference_fit_when_alpha_is_undetermined():
    # Two sound anchors is below Tier 3's floor of three.
    entries = [_entry(76, _on_curve(76)), _entry(150, _on_curve(150)),
               _rejected(300), _rejected(600)]
    strata, alpha, spread = rc.resolve_strata(entries, 22, REFERENCE)
    assert alpha is None and spread is None
    assert [s["status"] for s in strata] == ["reference_fallback"] * 4
    assert all(s["xi"] == REFERENCE["xi"] and s["theta"] == REFERENCE["theta"]
               for s in strata)
    # "reason" is a closed, machine-read vocabulary: a sound entry with no reason of its
    # own reports the fallback itself, never the raw UnusableCurveError text; a rejected
    # entry keeps ITS OWN reason, since that is more specific than the fallback.
    assert [s["reason"] for s in strata] == [
        "unusable_curve", "unusable_curve", "non_finite", "non_finite",
    ]

def test_resolve_strata_returns_entries_in_ladder_order():
    entries = [_entry(a, _on_curve(a)) for a in (76, 150, 300, 600, 900, 1500)]
    strata, _, _ = rc.resolve_strata(entries, 22, REFERENCE)
    assert [s["anchor"] for s in strata] == [76, 150, 300, 600, 900, 1500]


# --- review fixes: _calibrate_anchor and the anchored run_rnacalibrate path, on real jobs ---

# Derived from the factory, not restated: the point of _stratum is that there is one spelling.
STRATUM_KEYS = set(rc._stratum(900))


def _bare_anchor_command():
    # No k / length_arg bound — _calibrate_anchor supplies both per attempt, exactly as
    # run_rnacalibrate's make_anchor_command does.
    return partial(rc.build_command, executable="RNAcalibrate", target="t.fna",
                  max_target_length=100, randomize_targets=True)


def _nan_then(n_failures, calls):
    # Degenerate for the first n_failures rungs, sound after; records (command, env) per call.
    def fake_run(command, **kwargs):
        calls.append((command, kwargs.get("env")))
        return _completed("m 2 -nan -nan\n" if len(calls) <= n_failures
                          else "m 300 2.4 0.18\n")
    return fake_run


def _run_anchor(tmp_path, monkeypatch, fake_run, **overrides):
    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    kwargs = dict(k=10000, rng_seed=1, library_path="/x.so",
                  platform_name="linux", pin_clock=True)
    kwargs.update(overrides)
    query_path = rc._write_query_file(tmp_path, 0, "m", "ACGU")
    return rc._calibrate_anchor(query_path, "m", "ACGU", 900, _bare_anchor_command(), **kwargs)


def test_calibrate_anchor_uses_a_fresh_derived_seed_per_repair_attempt(tmp_path, monkeypatch):
    calls = []
    result = _run_anchor(tmp_path, monkeypatch, _nan_then(2, calls))

    seeds_seen = [env["FAKETIME"] for _command, env in calls]
    expected = [f"@{rc.derive_seed(1, 'ACGU', anchor=900, attempt=i)}" for i in range(3)]
    assert seeds_seen == expected
    assert len(set(seeds_seen)) == 3          # a fresh draw every retry, never a repeat
    assert result["status"] == "fitted"
    assert result["attempts"] == 3
    assert result["reason"] is None


def test_calibrate_anchor_recovers_at_a_later_repair_rung(tmp_path, monkeypatch):
    # non_finite through all three "fitted" rungs, sound once widened.
    result = _run_anchor(tmp_path, monkeypatch, _nan_then(3, []))

    assert result["status"] == "refitted_widened"
    assert result["attempts"] == 4
    assert result["reason"] is None
    assert result["length_arg"] == "900,450"   # divisor=2, the widening rung
    assert result["k"] == 10000                # k multiplier is 1 at this rung


def test_calibrate_anchor_sleeps_past_a_clock_second_on_unpinned_retries(tmp_path, monkeypatch):
    sleeps = []
    monkeypatch.setattr(rc.time, "sleep", lambda seconds: sleeps.append(seconds))

    result = _run_anchor(tmp_path, monkeypatch, _nan_then(2, []),
                         rng_seed=None, library_path=None, pin_clock=False)

    # No sleep before attempt 0; one before each of the two retries that follow it. Unpinned
    # commands are otherwise identical, so without this an unpinned retry would collapse
    # into the same draw RNAcalibrate's time()-seeded RNG already made.
    assert sleeps == [rc.CLOCK_TICK_SLEEP_S, rc.CLOCK_TICK_SLEEP_S]
    assert result["status"] == "fitted"
    assert result["derived_seed"] is None


def test_run_rnacalibrate_anchored_path_writes_uniform_11_key_strata(tmp_path, monkeypatch):
    query = tmp_path / "q.fa"
    query.write_text(">mirA\nACGUACGUACGU\n")
    target = tmp_path / "t.fna"
    # Lengths land one per cell of anchors (76, 150, 300): 50 -> 76, 150 -> 150, 250 -> 300.
    target.write_text(">g1\n" + "A" * 50 + "\n>g2\n" + "A" * 150 + "\n>g3\n" + "A" * 250 + "\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")
    monkeypatch.setattr(rc, "resolve_faketime_library", lambda *a: "/opt/lib/libfaketime.so.1")
    monkeypatch.setattr(rc, "verify_faketime", lambda *a, **k: None)

    def fake_run(command, **kwargs):
        length_value = command[command.index("-l") + 1]
        header, sequence = next(rc.iter_fasta_records(command[command.index("-q") + 1]))
        name = rc.query_key(header)
        for anchor in (76, 150, 300):
            if length_value == rc.anchor_length_arg(anchor)["value"]:
                # Exactly on an alpha=1.5 curve, so fit_alpha_robust recovers it with zero
                # residual spread and every anchor succeeds on attempt 0.
                theta = gf.predict_theta(1.5, 0.3, anchor, len(sequence))
                return _completed(f"{name} 300 2.4 {theta}\n")
        return _completed(f"{name} 300 2.4 0.2\n")   # the reference job's own -l

    monkeypatch.setattr(rc.subprocess, "run", fake_run)

    out = tmp_path / "out.json"
    rc.run_rnacalibrate(query, target, out, k=10000, max_target_length=100000,
                        randomize_targets=True, rng_seed=1, length_anchors=(76, 150, 300))
    data = json.loads(out.read_text())

    assert data["calibration"]["anchors"] == [
        {"anchor": 76, "length_arg": "76,25", "n_records": 1, "status": "attempted"},
        {"anchor": 150, "length_arg": "150,50", "n_records": 1, "status": "attempted"},
        {"anchor": 300, "length_arg": "300,100", "n_records": 1, "status": "attempted"},
    ]

    [entry] = data["calibration"]["per_query"]
    assert entry["alpha"] == pytest.approx(1.5, abs=1e-6)
    assert entry["alpha_residual_spread"] == pytest.approx(0.0, abs=1e-6)
    strata = entry["strata"]
    assert [s["anchor"] for s in strata] == [76, 150, 300]
    assert [s["status"] for s in strata] == ["fitted", "fitted", "fitted"]
    assert [s["attempts"] for s in strata] == [1, 1, 1]
    assert all(set(s.keys()) == STRATUM_KEYS for s in strata)


# --- the reference (whole-file) fit: classified like an anchor fit, never emitted degenerate ---

# The invalid rows RNAcalibrate prints, each with classify_fit's verdict.
DEGENERATE_REFERENCE_ROWS = [
    ("105 -nan -nan", "non_finite"),
    ("0 0.000000 0.000000", "degenerate_sample"),
    ("-2 0.000000 0.000000", "degenerate_sample"),
    ("300 2.4 0.000000", "non_positive_theta"),
    ("300 13970.0 -2022.0", "non_positive_theta"),
    ("300 -3.33 0.803", "non_positive_xi"),
]


def _fake_calibration_run(tmp_path, monkeypatch, reference_row, anchor_row="300 2.4 0.18",
                          **run_kwargs):
    query = tmp_path / "q.fa"
    query.write_text(">mirA\nACGUACGUACGU\n")
    target = tmp_path / "t.fna"
    target.write_text(">g1\n" + "A" * 50 + "\n>g2\n" + "A" * 150 + "\n>g3\n" + "A" * 250 + "\n")

    monkeypatch.setattr(rc, "which_required", lambda *a: "RNAcalibrate")
    monkeypatch.setattr(rc, "resolve_faketime_library", lambda *a: "/opt/lib/libfaketime.so.1")
    monkeypatch.setattr(rc, "verify_faketime", lambda *a, **k: None)
    reference_l = rc.build_length_arg(rc.summarise_lengths(rc.record_lengths(target)))["value"]

    def fake_run(command, **kwargs):
        is_reference = command[command.index("-l") + 1] == reference_l
        return _completed(f"mirA {reference_row if is_reference else anchor_row}\n")

    monkeypatch.setattr(rc.subprocess, "run", fake_run)
    out = tmp_path / "out.json"
    rc.run_rnacalibrate(query, target, out, k=10000, max_target_length=100000,
                        randomize_targets=True, rng_seed=1, **run_kwargs)
    return out


@pytest.mark.parametrize("row, reason", DEGENERATE_REFERENCE_ROWS)
def test_a_ladderless_run_refuses_a_degenerate_reference_fit(tmp_path, monkeypatch, row, reason):
    # Without a ladder the reference fit IS RNAhybrid's -d; theta <= 0 there is no null at all.
    with pytest.raises(rc.DegenerateFitError) as excinfo:
        _fake_calibration_run(tmp_path, monkeypatch, row)
    message = str(excinfo.value)
    assert "mirA" in message and reason in message
    assert "derived_seed=" in message and "rnacalibrate.k" in message
    assert not (tmp_path / "out.json").exists()


@pytest.mark.parametrize("row, reason", DEGENERATE_REFERENCE_ROWS)
def test_reference_fallback_never_copies_a_degenerate_reference_fit(
        tmp_path, monkeypatch, row, reason):
    # Every anchor fails all five rungs, so the curve is unusable and each cell would fall
    # back to the reference fit — which must be refused rather than copied into all three.
    with pytest.raises(rc.DegenerateFitError, match=reason):
        _fake_calibration_run(tmp_path, monkeypatch, row, anchor_row="2 -nan -nan",
                              length_anchors=(76, 150, 300))
    assert not (tmp_path / "out.json").exists()


def test_a_sound_reference_fit_is_written_through(tmp_path, monkeypatch):
    out = _fake_calibration_run(tmp_path, monkeypatch, "300 2.4 0.18")
    [entry] = json.loads(out.read_text())["calibration"]["per_query"]
    assert (entry["sample_size"], entry["xi"], entry["theta"]) == (300, 2.4, 0.18)


# --- one binary-dependent test, following the test_intarna.py precedent ---

ECOLI = Path("Data/Raw/genomes/GCF_000005845.2_ASM584v2_cds_rna_from_genomic.fna")

# rnahybrid-2.1.2-h7b50bb2_4 at FAKETIME=@1000, k=200, -l 900,300; both reproduce 2/2.
# k=200 shifts N (B becomes 200) but not the mechanism.
PINNED_MODES = [
    ("hsa-miR-1224-5p", "GUGAGGACUCGGGAGGUGG", 85, "non_finite"),
    ("hsa-miR-2054", "CUGUAAUAUAAAUUUAAUUUAUU", -2, "degenerate_sample"),
]


@pytest.mark.skipif(not ECOLI.exists(), reason="E. coli target FASTA not present")
@pytest.mark.parametrize("mirna,sequence,expected_n,expected_reason", PINNED_MODES)
def test_the_binarys_two_failure_modes_are_still_what_the_guard_expects(
        tmp_path, mirna, sequence, expected_n, expected_reason):
    # The fit window is a hard-coded ABSOLUTE 2.0 cutoff on a scale that slides with
    # ln(m*n): a range entirely above it yields -nan, one entirely below caps `start` and
    # yields 0.000000 with a negative "sample size". Both are binary properties, so pin them.
    executable = rc.which_required("RNAcalibrate", "rnacalibrate")
    library = rc.resolve_faketime_library(sys.platform)
    query = tmp_path / "q.fa"
    query.write_text(f">{mirna}\n{sequence}\n")

    command = rc.build_command(
        executable=executable, query=str(query), target=str(ECOLI), k=200,
        max_target_length=50000, length_arg=rc.anchor_length_arg(900),
        randomize_targets=True,
    )
    env = rc.build_faketime_env(1000, library, sys.platform, os.environ)
    stdout = subprocess.run(command, check=True, capture_output=True, text=True,
                            env=env).stdout

    row = rc.expect_single_row(
        rc.parse_rnacalibrate_output(stdout), mirna)
    assert row["sample_size"] == expected_n
    assert gf.classify_fit(row["sample_size"], row["xi"], row["theta"]) == expected_reason
