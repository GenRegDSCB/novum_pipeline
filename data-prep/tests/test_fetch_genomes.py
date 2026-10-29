from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import contextlib
import errno
import gzip
import io
import tempfile
import zipfile
import zlib
from pathlib import Path
from unittest import mock

import requests

import fetch_genomes
from fetch_genomes import (assembly_stem, rna_ftp_url, target_dest, append_gz_fasta,
                           download_targets, _resolve_one, _API_BASE, resolve_assemblies)


def test_rna_ftp_url_splits_accession_digits():
    assert rna_ftp_url("GCF_000005845.2", "ASM584v2") == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/"
        "GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_rna_from_genomic.fna.gz")


def test_rna_ftp_url_handles_gca_prefix():
    assert "/GCA/014/131/755/" in rna_ftp_url("GCA_014131755.1", "ASM1413175v1")


def test_assembly_stem_sanitizes_unsafe_characters():
    assert assembly_stem("GCF_1.1", "ASM 1 (v2)") == "GCF_1.1_ASM_1__v2_"


def test_target_dest_names_both_branches():
    assert target_dest("g", "GCF_1.1", "ASM1", with_rna=True).name == \
        "GCF_1.1_ASM1_cds_rna_from_genomic.fna"
    assert target_dest("g", "GCF_1.1", "ASM1", with_rna=False).name == \
        "GCF_1.1_ASM1_cds_from_genomic.fna"


def _gz(path, text):
    Path(path).write_bytes(_rna_gz_bytes(text))
    return path


def test_append_gz_fasta_concatenates_cds_before_rna():
    d = Path(tempfile.mkdtemp())
    part = d / "x.part"
    part.write_bytes(b">lcl|X_cds_1 [gene=a]\nAAAA\n")
    added = append_gz_fasta(part, _gz(d / "r.gz", ">lcl|X_ncrna_1 [product=b]\nGGGG\n"))
    assert added == 1
    assert part.read_text() == (">lcl|X_cds_1 [gene=a]\nAAAA\n"
                                ">lcl|X_ncrna_1 [product=b]\nGGGG\n")


def test_append_gz_fasta_inserts_missing_newline():
    d = Path(tempfile.mkdtemp())
    part = d / "x.part"
    part.write_bytes(b">lcl|X_cds_1\nAAAA")          # no trailing newline
    append_gz_fasta(part, _gz(d / "r.gz", ">lcl|X_rrna_1\nGGGG\n"))
    assert part.read_text() == ">lcl|X_cds_1\nAAAA\n>lcl|X_rrna_1\nGGGG\n"


def test_append_gz_fasta_counts_every_record():
    d = Path(tempfile.mkdtemp())
    part = d / "x.part"
    part.write_bytes(b">lcl|X_cds_1\nAA\n")
    assert append_gz_fasta(part, _gz(d / "r.gz", ">a\nGG\n>b\nCC\n>c\nTT\n")) == 3


class _Resp:
    """Minimal stand-in for a requests.Response used in streaming mode."""
    def __init__(self, status_code=200, body=b""):
        self.status_code = status_code
        self._body = body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size=1):
        yield self._body


class _StubSession:
    """Serves one canned CDS zip for POST and RNA payload(s) for GET.

    `rna_gz` may be:
      - None: every GET 404s (no RNA available)
      - bytes: every GET returns that same body with status 200
      - dict[str, bytes | None]: looked up by whichever key appears as a
        substring of the GET url, so different accessions sharing one batch
        can each get their own RNA outcome (404 / ok / corrupt payload).
    """
    def __init__(self, zip_bytes=b"", rna_gz=None):
        self.zip_bytes = zip_bytes
        self.rna_gz = rna_gz
        self.posts = []
        self.gets = []

    def post(self, url, **kw):
        self.posts.append(url)
        return _Resp(200, self.zip_bytes)

    def get(self, url, **kw):
        self.gets.append(url)
        body = self.rna_gz
        if isinstance(body, dict):
            body = next((v for key, v in body.items() if key in url), None)
        return _Resp(404) if body is None else _Resp(200, body)


def _cds_zip_multi(pairs):
    """Bundles several (accession, text) entries into one ZIP — for tests where a
    single batched POST must serve more than one accession."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for acc, text in pairs:
            zf.writestr(f"ncbi_dataset/data/{acc}/cds_from_genomic.fna", text)
    return buf.getvalue()


def _cds_zip(acc, text):
    return _cds_zip_multi([(acc, text)])


def _rna_gz_bytes(text):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
        fh.write(text.encode())
    return buf.getvalue()


def _big_rna_text(n=500):
    """A multi-record RNA FASTA big enough that a truncated gzip of it forces
    the decompressor to yield many complete lines before the stream ends —
    unlike a single short record, where the whole payload can sit in one
    deflate block and a truncation raises before anything is emitted at all.
    Needed so a corrupt-gzip test actually exercises "partial data already
    written to disk", instead of accidentally looking clean by luck."""
    lines = []
    for i in range(n):
        lines.append(f">lcl|X_ncrna_{i} [product=b]")
        lines.append("GGGGCCCCAAAATTTT" * 4)
    return "\n".join(lines) + "\n"


def test_download_targets_merges_and_cleans_up_part_file():
    d = Path(tempfile.mkdtemp())
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1 [gene=a]\nAAAA\n"),
                     _rna_gz_bytes(">lcl|X_ncrna_1 [product=b]\nGGGG\n"))
    paths, status = download_targets({"ecoli": ("GCF_1.1", "ASM1")}, d, s, delay=0)
    assert status == {"ecoli": "merged"}
    out = Path(paths["ecoli"])
    assert out.name == "GCF_1.1_ASM1_cds_rna_from_genomic.fna"
    assert out.read_text().count(">") == 2
    assert list(d.glob("*.part")) == []


def test_download_targets_degrades_to_cds_only_on_404():
    d = Path(tempfile.mkdtemp())
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"), rna_gz=None)
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0)
    assert status == {"t": "missing"}
    assert Path(paths["t"]).name == "GCF_1.1_ASM1_cds_from_genomic.fna"
    assert list(d.glob("*.part")) == []


def test_download_targets_does_not_treat_cds_only_file_as_cached():
    d = Path(tempfile.mkdtemp())
    (d / "GCF_1.1_ASM1_cds_from_genomic.fna").write_text(">stale\nAAAA\n")
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nCCCC\n"), rna_gz=None)
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0)
    assert status == {"t": "missing"}                    # refetched, not reused
    assert Path(paths["t"]).read_text() == ">lcl|X_cds_1\nCCCC\n"


def test_download_targets_uses_merged_file_as_cache():
    d = Path(tempfile.mkdtemp())
    (d / "GCF_1.1_ASM1_cds_rna_from_genomic.fna").write_text(">cached\nAAAA\n")
    s = _StubSession()
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0)
    assert status == {"t": "cached"}
    assert s.posts == [] and s.gets == []                # no network at all


def test_download_targets_shares_one_download_across_taxa():
    d = Path(tempfile.mkdtemp())
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"),
                     _rna_gz_bytes(">lcl|X_ncrna_1\nGGGG\n"))
    hits = {"a": ("GCF_1.1", "ASM1"), "b": ("GCF_1.1", "ASM1"), "c": None}
    paths, status = download_targets(hits, d, s, delay=0)
    assert paths["a"] == paths["b"]
    assert "c" not in paths and "c" not in status
    assert len(s.gets) == 1                              # deduped by accession


def test_download_targets_corrupt_rna_errors_one_accession_without_killing_batch():
    # A truncated gzip still has a valid header, so the HTTP transfer "succeeds"
    # (no requests exception) and the failure only surfaces once gzip tries to
    # read past the missing end-of-stream marker (EOFError, not OSError).
    # Use a big multi-record payload (not a single short line): a short payload
    # can sit inside one deflate block and a truncation raises before anything
    # is emitted at all, which would make the corruption-on-disk scenario this
    # test exists to catch impossible to trigger.
    full = _rna_gz_bytes(_big_rna_text())
    truncated = full[: int(len(full) * 0.9)]

    d = Path(tempfile.mkdtemp())
    zip_bytes = _cds_zip_multi([
        ("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"),
        ("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n"),
    ])
    s = _StubSession(zip_bytes, rna_gz={
        "GCF_1.1": truncated,                             # corrupt: triggers EOFError
        "GCF_2.2": _rna_gz_bytes(">lcl|Y_ncrna_1\nGGGG\n"),  # healthy sibling in the same batch
    })
    hits = {"bad": ("GCF_1.1", "ASM1"), "good": ("GCF_2.2", "ASM2")}
    paths, status = download_targets(hits, d, s, delay=0)

    # The corrupt accession is isolated as its own "error" status; the CDS half
    # it already had still survives on disk (RNA merge failed, CDS did not).
    assert status["bad"] == "error"
    assert Path(paths["bad"]).name == "GCF_1.1_ASM1_cds_from_genomic.fna"
    # Critically: the second accession in the same batch still completed —
    # one accession's corrupt payload must not abort the others.
    assert status["good"] == "merged"
    assert Path(paths["good"]).name == "GCF_2.2_ASM2_cds_rna_from_genomic.fna"
    # No orphaned staging file for either accession.
    assert list(d.glob("*.part")) == []

    # Byte-content proof (Finding 1): the file named *_cds_from_genomic.fna must
    # actually contain ONLY the CDS record — never a partial RNA fragment that
    # leaked onto disk before the truncated gzip raised mid-stream.
    assert Path(paths["bad"]).read_text() == ">lcl|X_cds_1\nAAAA\n"


def _assert_rollback_on_truncated_gzip(original: bytes):
    # A big multi-record payload so the decompressor yields several complete lines
    # (which the un-fixed code wrote straight to disk) before the truncation raises —
    # a short single-record payload can fail before emitting anything at all, which
    # would make the assertion pass by luck.
    d = Path(tempfile.mkdtemp())
    part = d / "x.part"
    part.write_bytes(original)

    full = _rna_gz_bytes(_big_rna_text())
    truncated_gz = d / "r.gz"
    truncated_gz.write_bytes(full[: int(len(full) * 0.9)])

    raised = False
    try:
        append_gz_fasta(part, truncated_gz)
    except (EOFError, zlib.error):
        raised = True
    assert raised, "append_gz_fasta must propagate the gzip failure"
    assert part.read_bytes() == original      # rolled back to pristine pre-call bytes


def test_append_gz_fasta_rolls_back_part_on_corrupt_gzip():
    # Unit-level proof that append_gz_fasta itself is transactional.
    _assert_rollback_on_truncated_gzip(b">lcl|X_cds_1 [gene=a]\nAAAA\n")


def test_append_gz_fasta_rolls_back_empty_part_on_corrupt_gzip():
    # rollback_to == 0 edge case: an empty .part must stay empty after a failed append.
    _assert_rollback_on_truncated_gzip(b"")


def test_download_targets_upgrade_unlinks_stale_cds_only_sibling():
    # Migration case: a pre-existing CDS-only file for an accession must be
    # removed once a re-run successfully merges RNA, so exactly one file per
    # accession survives (README/spec: "a re-run replaces the earlier file").
    d = Path(tempfile.mkdtemp())
    stale = d / "GCF_1.1_ASM1_cds_from_genomic.fna"
    stale.write_text(">stale_cds\nAAAA\n")
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nCCCC\n"),
                     _rna_gz_bytes(">lcl|X_ncrna_1\nGGGG\n"))
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0, force=True)
    assert status == {"t": "merged"}
    assert Path(paths["t"]).name == "GCF_1.1_ASM1_cds_rna_from_genomic.fna"
    remaining = sorted(p.name for p in d.glob("GCF_1.1_ASM1_*"))
    assert remaining == ["GCF_1.1_ASM1_cds_rna_from_genomic.fna"]  # stale sibling gone
    assert not stale.exists()


def test_download_targets_forced_degrade_never_deletes_merged_sibling():
    # Asymmetric-safety case: a good pre-existing merged file must survive a
    # forced re-run that degrades to CDS-only (e.g. a transient FTP 404) —
    # stale merged data beats fresh degraded data and self-heals later.
    d = Path(tempfile.mkdtemp())
    good = d / "GCF_1.1_ASM1_cds_rna_from_genomic.fna"
    good.write_text(">good_merged\nAAAA\nGGGG\n")
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nCCCC\n"), rna_gz=None)  # RNA 404s
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0, force=True)
    assert status == {"t": "missing"}
    assert good.exists()                                   # merged file NOT deleted
    assert good.read_text() == ">good_merged\nAAAA\nGGGG\n"  # and untouched


def test_download_targets_force_refetches_despite_merged_cache():
    d = Path(tempfile.mkdtemp())
    (d / "GCF_1.1_ASM1_cds_rna_from_genomic.fna").write_text(">stale\nAAAA\n")
    s = _StubSession(_cds_zip("GCF_1.1", ">lcl|X_cds_1\nCCCC\n"),
                     _rna_gz_bytes(">lcl|X_ncrna_1\nGGGG\n"))
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0, force=True)
    assert status == {"t": "merged"}
    assert s.posts != [] and s.gets != []                 # network happened despite the cached file
    assert Path(paths["t"]).read_text() == (">lcl|X_cds_1\nCCCC\n"
                                            ">lcl|X_ncrna_1\nGGGG\n")


def test_download_targets_skips_rna_when_cds_missing_from_zip():
    d = Path(tempfile.mkdtemp())
    # The ZIP payload has no cds_from_genomic.fna entry for GCF_1.1 at all
    # (e.g. NCBI's package omitted it) — RNA must never become the sole content.
    s = _StubSession(_cds_zip("GCF_OTHER", ">lcl|Z_cds_1\nAAAA\n"),
                     _rna_gz_bytes(">lcl|X_ncrna_1\nGGGG\n"))
    paths, status = download_targets({"t": ("GCF_1.1", "ASM1")}, d, s, delay=0)
    assert status == {"t": "-"}
    assert paths == {"t": None}
    assert s.gets == []                                   # RNA fetch skipped entirely
    assert list(d.glob("*")) == []                         # nothing written for this accession


#----- download failures: recorded per accession, no temp ZIP outlives the call -----#

class _Abort(BaseException):
    """Anything that is not a download failure — e.g. Ctrl-C mid-run."""


class _BatchSession(_StubSession):
    """One POST outcome per batch, in order: ZIP bytes, an HTTP status, or an exception."""
    def __init__(self, batches, rna_gz=None):
        super().__init__(rna_gz=rna_gz)
        self.batches = list(batches)

    def post(self, url, **kw):
        self.posts.append(url)
        outcome = self.batches.pop(0)
        if isinstance(outcome, BaseException):
            raise outcome
        if isinstance(outcome, _Resp):
            return outcome
        return _Resp(outcome) if isinstance(outcome, int) else _Resp(200, outcome)


class _CutStream(_Resp):
    """A 200 whose body breaks off with `exc` after the first chunk."""
    def __init__(self, exc):
        super().__init__(200, b"PK\x03\x04partial")
        self.exc = exc

    def iter_content(self, chunk_size=1):
        yield self._body
        raise self.exc


_REAL_NAMED_TEMPFILE = tempfile.NamedTemporaryFile


class _FullDisk:
    """A temp file whose last flush fails on close, as it does on a full disk."""
    def __init__(self, **kw):
        self._fh = _REAL_NAMED_TEMPFILE(**kw)
        self.name, self.write = self._fh.name, self._fh.write

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._fh.close()
        raise OSError(errno.ENOSPC, "No space left on device")


@contextlib.contextmanager
def _one_accession_per_batch():
    """Batch size 1, and every temp file redirected into a fresh dir the test can inspect."""
    tmp = Path(tempfile.mkdtemp())
    with mock.patch.object(tempfile, "tempdir", str(tmp)), \
            mock.patch.object(fetch_genomes, "_BATCH_SIZE", 1):
        yield tmp


_TWO_HITS = {"a": ("GCF_1.1", "ASM1"), "b": ("GCF_2.2", "ASM2")}


def test_a_failed_batch_is_recorded_and_the_others_still_land():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"), 400])
    with _one_accession_per_batch() as tmp:
        paths, status = download_targets(_TWO_HITS, d, s, delay=0)
    assert status == {"a": "missing", "b": "-"}
    assert Path(paths["a"]).read_text() == ">lcl|X_cds_1\nAAAA\n"
    assert paths["b"] is None
    assert list(tmp.iterdir()) == []                      # batch 1's ZIP did not leak


def test_a_failed_first_batch_does_not_stop_the_later_ones():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([requests.ConnectionError("reset"),
                       _cds_zip("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n")])
    with _one_accession_per_batch() as tmp:
        paths, status = download_targets(_TWO_HITS, d, s, delay=0)
    assert status == {"a": "-", "b": "missing"}
    assert paths["a"] is None and Path(paths["b"]).exists()
    assert list(tmp.iterdir()) == []


def test_an_interrupted_download_deletes_the_zips_already_streamed():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"), _Abort()])
    with _one_accession_per_batch() as tmp:
        try:
            download_targets(_TWO_HITS, d, s, delay=0)
        except _Abort:
            pass
        else:
            raise AssertionError("a non-download error must still propagate")
        assert list(tmp.iterdir()) == []


def test_an_abort_mid_stream_deletes_the_zip_being_written():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"), _CutStream(_Abort())])
    with _one_accession_per_batch() as tmp:
        try:
            download_targets(_TWO_HITS, d, s, delay=0)
        except _Abort:
            pass
        else:
            raise AssertionError("a non-download error must still propagate")
        assert list(tmp.iterdir()) == []                  # neither batch 1's nor batch 2's


def test_a_zip_that_fails_to_close_is_recorded_and_deleted():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"),
                       _cds_zip("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n")])
    with _one_accession_per_batch() as tmp, \
            mock.patch.object(tempfile, "NamedTemporaryFile", _FullDisk):
        paths, status = download_targets(_TWO_HITS, d, s, delay=0)
    assert status == {"a": "-", "b": "-"}
    assert list(tmp.iterdir()) == []


def test_an_interrupted_extraction_deletes_every_zip():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([_cds_zip("GCF_1.1", ">lcl|X_cds_1\nAAAA\n"),
                       _cds_zip("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n")])
    with _one_accession_per_batch() as tmp, \
            mock.patch.object(fetch_genomes.shutil, "copyfileobj", side_effect=_Abort):
        try:
            download_targets(_TWO_HITS, d, s, delay=0)
        except _Abort:
            pass
        assert list(tmp.iterdir()) == []                  # batch 2's ZIP was never reached


def test_a_corrupt_package_is_recorded_and_leaks_nothing():
    d = Path(tempfile.mkdtemp())
    s = _BatchSession([b"not a zip", _cds_zip("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n")])
    with _one_accession_per_batch() as tmp:
        paths, status = download_targets(_TWO_HITS, d, s, delay=0)
    assert status == {"a": "-", "b": "missing"}
    assert list(tmp.iterdir()) == []


def _crc_broken_zip(pairs):
    """A stored (uncompressed) ZIP whose member holding "BADSEQ" fails its CRC on read."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as zf:
        for acc, text in pairs:
            zf.writestr(f"ncbi_dataset/data/{acc}/cds_from_genomic.fna", text)
    data = bytearray(buf.getvalue())
    i = data.index(b"BADSEQ")
    data[i] = ord("b")
    return bytes(data)


def test_a_member_failing_mid_read_is_recorded_and_leaves_no_part():
    d = Path(tempfile.mkdtemp())
    zip_bytes = _crc_broken_zip([("GCF_1.1", ">lcl|X_cds_1\nBADSEQ\n"),
                                 ("GCF_2.2", ">lcl|Y_cds_1\nCCCC\n")])
    s = _StubSession(zip_bytes, rna_gz=None)
    paths, status = download_targets(_TWO_HITS, d, s, delay=0)
    assert status == {"a": "-", "b": "missing"}
    assert paths["a"] is None
    assert list(d.glob("*.part")) == []


#----- _resolve_one: taxon-name search vs accession pin -----#

class _ReportResp:
    def __init__(self, payload, status_code=200):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        return self._payload


class _ReportSession:
    """Serves canned dataset_report JSON and records every (url, params) call.

    A `reports` value may be a plain list, or a dict with "reference"/"page" keys
    to serve the taxon path's two tiers separately (tier 1 is the call carrying
    filters.reference_only). An unknown url yields no reports, which is what NCBI
    returns for an accession that does not exist.
    """
    def __init__(self, reports):
        self.reports = reports
        self.calls = []

    def get(self, url, params=None, timeout=None):
        params = dict(params or {})
        self.calls.append((url, params))
        entry = self.reports.get(url, [])
        if isinstance(entry, dict):
            entry = entry["reference" if "filters.reference_only" in params else "page"]
        return _ReportResp({"reports": entry})


def _report(acc, asm, category="", level="Complete Genome"):
    return {"accession": acc,
            "assembly_info": {"assembly_name": asm, "refseq_category": category,
                              "assembly_level": level}}


def _url(kind, ident):
    return f"{_API_BASE}/genome/{kind}/{ident}/dataset_report"


class _BoomSession:
    """Every lookup raises, as a DNS failure or a dropped connection would."""
    def get(self, *a, **k):
        raise requests.ConnectionError("network down")


def test_resolve_assemblies_reports_which_lookups_errored():
    # A transient network failure is not "NCBI has no such genome". Collapsing the
    # two makes prepare_inputs guess the reason back from the taxon key's shape,
    # which reads a dropped connection on a pinned line as a typo.
    hits, errored = resolve_assemblies(["escherichia coli", "GCF_000284435.1"],
                                       _BoomSession(), 0)
    assert hits == {"escherichia coli": None, "GCF_000284435.1": None}
    assert errored == {"escherichia coli", "GCF_000284435.1"}


def test_resolve_assemblies_reports_no_errors_on_a_clean_miss():
    s = _ReportSession({})
    hits, errored = resolve_assemblies(["nosuchbug"], s, 0)
    assert hits == {"nosuchbug": None}
    assert errored == set()


def test_resolve_one_pins_an_accession_via_the_accession_endpoint():
    s = _ReportSession({_url("accession", "GCF_000284435.1"):
                        [_report("GCF_000284435.1", "ASM28443v1")]})
    assert _resolve_one("GCF_000284435.1", s) == ("GCF_000284435.1", "ASM28443v1")
    assert len(s.calls) == 1
    assert "/genome/accession/" in s.calls[0][0]


def test_resolve_one_accession_pin_sends_no_selection_filters():
    # reference_only would reject any non-reference assembly, and assembly_source
    # would reject a GenBank pin outright. A pin selects nothing: it names one row.
    s = _ReportSession({_url("accession", "GCA_000284435.1"):
                        [_report("GCA_000284435.1", "ASM28443v1")]})
    assert _resolve_one("GCA_000284435.1", s) == ("GCA_000284435.1", "ASM28443v1")
    params = s.calls[0][1]
    assert "filters.reference_only" not in params
    assert "filters.assembly_source" not in params


def test_resolve_one_bad_accession_never_falls_back_to_a_taxon_search():
    # A typo must fail closed. Falling through to the name path would resolve the
    # pin to some other assembly and write it into the pinned output tree.
    s = _ReportSession({})
    assert _resolve_one("GCF_999999999.9", s) is None
    assert len(s.calls) == 1
    assert "/genome/accession/" in s.calls[0][0]


def test_resolve_one_still_prefers_the_reference_assembly_for_a_name():
    s = _ReportSession({_url("taxon", "escherichia%20coli"):
                        {"reference": [_report("GCF_000005845.2", "ASM584v2", "reference genome")],
                         "page": []}})
    assert _resolve_one("escherichia coli", s) == ("GCF_000005845.2", "ASM584v2")
    assert len(s.calls) == 1
    assert "/genome/taxon/" in s.calls[0][0]
    assert s.calls[0][1]["filters.reference_only"] == "true"


def test_resolve_one_name_falls_back_to_ranking_a_candidate_page():
    s = _ReportSession({_url("taxon", "candidatus%20arthromitus"):
                        {"reference": [],
                         "page": [_report("GCF_1.1", "ASM1", "", "Contig"),
                                  _report("GCF_2.2", "ASM2", "representative genome")]}})
    assert _resolve_one("candidatus arthromitus", s) == ("GCF_2.2", "ASM2")
    assert len(s.calls) == 2


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
