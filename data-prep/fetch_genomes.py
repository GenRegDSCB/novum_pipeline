"""
fetch_genomes.py — resolve taxa to NCBI assemblies and download each one's target
FASTA. An importable library for prepare_inputs.py; it has no CLI.

A taxon field is either an organism name, ranked to an assembly by select_best,
or an assembly accession, which pins that exact one (see pairs.is_accession).

Each target is assembled from two different NCBI services, because Datasets v2
publishes no rna_from_genomic.fna for prokaryotes:

  CDS   POST /datasets/v2/genome/download, batched, extracted from the ZIP
  RNA   the FTP mirror, at a URL derived from (accession, assembly name)

Both stream into one <accession>_<assembly>_cds_rna_from_genomic.fna through a
staged .part file, so an interrupted run cannot leave a half-merged target that
the next run mistakes for complete. Only that merged name counts as a cache hit;
a CDS-only file is a degraded result and is retried.

Environment:
  NCBI_API_KEY     Optional. Raises rate limit from ~3 to ~10 req/s.
  NCBI_API_EMAIL   Optional. Appended to the User-Agent header (NCBI etiquette).
"""
from __future__ import annotations

import gzip
import os
import shutil
import tempfile
import time
import zipfile
import zlib
from pathlib import Path
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from pairs import is_accession, taxon_to_query

# NCBI Datasets v2 REST API. Switch this if NCBI moves the endpoint.
_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"

# NCBI etiquette: identify the client. Email comes from NCBI_API_EMAIL if set.
_DEFAULT_USER_AGENT = "ncbi-cds-downloader/1.0"

# The CDS set NCBI ships inside each ncbi_dataset/data/<accession>/ directory.
_CDS_FILENAME = "cds_from_genomic.fna"

# A corrupt gzip/deflate stream transfers cleanly and fails only on read, raising EOFError or
# zlib.error, neither of which subclasses OSError.
_CORRUPT_STREAM_ERRORS = (OSError, EOFError, zlib.error)

# Accessions per POST /genome/download. Bacterial CDS sets are small; this bounds
# how much NCBI packages server-side in one request.
_BATCH_SIZE = 200

# Datasets v2 serves no rna_from_genomic.fna for prokaryotes, so the RNA half comes
# from the FTP mirror, at a URL derived from (acc, asm) rather than queried.
_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
_RNA_SUFFIX = "rna_from_genomic.fna.gz"

# NCBI substitutes these characters in the assembly-name part of its paths.
_NAME_UNSAFE = str.maketrans({c: "_" for c in " #()/"})

# HTTP status codes worth retrying with backoff (rate limits + transient server errors).
_RETRIABLE_STATUS = (429, 500, 502, 503, 504)

# Per-request delay between taxon resolution calls. NCBI: 3 req/s without a key, 10 with one.
_RATE_DELAY_NO_KEY = 0.35
_RATE_DELAY_WITH_KEY = 0.12

# Assembly-selection ranking: lower tuple sorts first.
_CATEGORY_RANK = {"reference genome": 0, "representative genome": 1, "": 2}
_LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


#----- Picks the best assembly: category > completeness > newest release date -----#
def select_best(candidates: list[dict]) -> tuple[str, str] | None:
    if not candidates:
        return None

    def info(rep):
        return rep.get("assembly_info", {}) or {}

    # Newest-first pre-sort, then pick the best (category, level) with a stable
    # min — ties keep the newest because the pre-sort already ordered by date desc.
    ranked = sorted(candidates, key=lambda r: info(r).get("release_date", ""), reverse=True)
    best = min(
        ranked,
        key=lambda r: (
            _CATEGORY_RANK.get(info(r).get("refseq_category", "") or "", 3),
            _LEVEL_RANK.get(info(r).get("assembly_level", "") or "", 4),
        ),
    )
    acc = best.get("accession") or best.get("current_accession")
    if not acc:
        return None
    return acc, info(best).get("assembly_name", "")


#----- Session plus its throttle, from the environment; single owner of NCBI etiquette policy -----#
def session_from_env() -> tuple[requests.Session, float]:
    api_key = os.environ.get("NCBI_API_KEY")
    session = make_session(api_key=api_key, contact_email=os.environ.get("NCBI_API_EMAIL"))
    return session, (_RATE_DELAY_WITH_KEY if api_key else _RATE_DELAY_NO_KEY)


#----- Resolves one taxon name (reference first, else best of a page) or pins an accession -----#
def _resolve_one(name: str, session: requests.Session) -> tuple[str, str] | None:
    pinned = is_accession(name)
    kind = "accession" if pinned else "taxon"
    query = quote(taxon_to_query(name), safe="")
    base = f"{_API_BASE}/genome/{kind}/{query}/dataset_report"

    def reports_for(params):
        r = session.get(base, params=params, timeout=30)
        return r.json().get("reports", []) if r.status_code == 200 else None

    if pinned:
        # No fallback by design: a typo must fail closed rather than resolve to
        # some other assembly and be written into the pinned output tree.
        return select_best(reports_for({}) or [])
    # Tier 1: a designated reference assembly (clean filter, one row).
    rep = reports_for({"filters.reference_only": "true",
                       "filters.assembly_source": "refseq", "page_size": "1"})
    if rep:
        return select_best(rep)
    # Tier 2: no reference — rank a page of RefSeq candidates client-side.
    rep = reports_for({"filters.assembly_source": "refseq", "page_size": "50"})
    return select_best(rep) if rep is not None else None


#----- Resolves each taxon to (accession, assembly_name); returns the hits and which lookups errored -----#
def resolve_assemblies(names: list[str], session: requests.Session,
                       delay: float) -> tuple[dict[str, tuple[str, str] | None], set[str]]:
    result: dict[str, tuple[str, str] | None] = {}
    errored: set[str] = set()
    for i, name in enumerate(names, 1):
        print(f"  [{i:>3}/{len(names)}] {name!r} ... ", end="", flush=True)
        try:
            hit = _resolve_one(name, session)
        except requests.RequestException as exc:
            print(f"error ({exc.__class__.__name__})")
            hit = None
            errored.add(name)
        else:
            print(f"{hit[0]} ({hit[1]})" if hit else "no_assembly")
        result[name] = hit
        time.sleep(delay)
    return result, errored


#----- "<accession>_<assembly name>", the stem NCBI uses in both paths and filenames -----#
def assembly_stem(acc: str, asm: str) -> str:
    return f"{acc}_{asm.translate(_NAME_UNSAFE)}"


#----- FTP URL of an assembly's rna_from_genomic.fna.gz, derived from (acc, asm) -----#
def rna_ftp_url(acc: str, asm: str) -> str:
    # GCF_000005845.2 -> prefix "GCF", digits "000005845" -> .../GCF/000/005/845/
    prefix, _, rest = acc.partition("_")
    digits = rest.split(".", 1)[0]
    stem = assembly_stem(acc, asm)
    return (f"{_FTP_BASE}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}"
            f"/{stem}/{stem}_{_RNA_SUFFIX}")


#----- On-disk name for a target FASTA; `with_rna` picks the merged vs CDS-only name -----#
def target_dest(out_dir, acc: str, asm: str, with_rna: bool = True) -> Path:
    kind = "cds_rna" if with_rna else "cds"
    return Path(out_dir) / f"{assembly_stem(acc, asm)}_{kind}_from_genomic.fna"


#----- Appends a gzipped FASTA onto `part_path`, healing a missing newline; returns records added -----#
def append_gz_fasta(part_path, gz_path) -> int:
    part_path = Path(part_path)
    records = 0
    with part_path.open("rb+") as dst:
        dst.seek(0, os.SEEK_END)
        rollback_to = dst.tell()  # pristine CDS-only length; restore here on any failure
        # A CDS chunk missing its trailing newline would glue onto the first RNA header.
        if rollback_to:
            dst.seek(-1, os.SEEK_END)
            if dst.read(1) != b"\n":
                dst.write(b"\n")
        try:
            with gzip.open(gz_path, "rb") as src:
                for line in src:
                    if line.startswith(b">"):
                        records += 1
                    dst.write(line)
        except Exception:
            # A gzip failing mid-stream has already written partial RNA; roll back so
            # the caller's *_cds_from_genomic.fna name stays honest.
            dst.seek(rollback_to)
            dst.truncate()
            raise
    return records


#----- Yields (accession, open member) for each matching file in the ZIPs -----#
def iter_package_members(zip_paths, filename):
    for zip_path in zip_paths:
        try:
            with zipfile.ZipFile(zip_path) as zf:
                for member in zf.namelist():
                    parts = member.split("/")
                    # NCBI Datasets ZIP layout: ncbi_dataset/data/<accession>/<filename>
                    if (len(parts) == 4 and parts[0] == "ncbi_dataset"
                            and parts[1] == "data" and parts[3] == filename):
                        with zf.open(member) as src:
                            yield parts[2], src
        except zipfile.BadZipFile as exc:
            # A corrupt package drops only the accessions it holds.
            print(f"    unreadable CDS package skipped: {exc}")


#----- Extracts each accession's cds_from_genomic.fna from the batched ZIPs into a <stem>.part
#      file. The ZIPs live in one temp dir, so every exit path removes them -----#
def _download_cds_parts(accessions, asm_by_acc, out_dir, session):
    parts: dict[str, Path] = {}
    with tempfile.TemporaryDirectory() as zip_dir:
        zip_paths = download_batch(accessions, session, zip_dir)
        for acc, src in iter_package_members(zip_paths, _CDS_FILENAME):
            asm = asm_by_acc.get(acc)
            if asm is None:
                continue
            # Staged as .part so an interrupted run can't leave a half-merged
            # file the next run mistakes for complete.
            part = out_dir / f"{assembly_stem(acc, asm)}.part"
            try:
                with part.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
            except (*_CORRUPT_STREAM_ERRORS, zipfile.BadZipFile) as exc:
                # A member failing mid-read (bad CRC, truncated deflate) drops only its accession.
                part.unlink(missing_ok=True)
                print(f"    CDS extract failed for {acc}: {exc.__class__.__name__}")
                continue
            parts[acc] = part
    return parts


#----- Streams a response body to a temp file the caller owns; self-cleans a failed transfer -----#
def _stream_to_tempfile(response, suffix: str, dir=None) -> Path:
    fh = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, dir=dir)
    path = Path(fh.name)
    try:
        # The close is inside the try: its final flush can fail too (a full disk).
        with fh:
            for chunk in response.iter_content(chunk_size=1 << 20):  # 1 MiB
                fh.write(chunk)
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return path


#----- Fetches the assembly's RNA FASTA from FTP and appends it to `part`; returns merged/missing/error -----#
def _append_rna(part, acc, asm, session):
    url = rna_ftp_url(acc, asm)
    tmp = None
    try:
        # The api-key header is scoped to the API host; strip it for the FTP mirror.
        r = session.get(url, headers={"api-key": None}, stream=True, timeout=300)
        if r.status_code == 404:
            print(f"    no rna_from_genomic for {acc}")
            return "missing"
        r.raise_for_status()
        tmp = _stream_to_tempfile(r, ".gz")
        append_gz_fasta(part, tmp)
        return "merged"
    except (requests.RequestException, *_CORRUPT_STREAM_ERRORS) as exc:
        # The URL is printed so a derivation bug can't look like "NCBI has no RNA".
        print(f"    RNA fetch failed for {acc}: {exc.__class__.__name__} — {url}")
        return "error"
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


#----- Merges RNA onto a staged .part, renames it to its final name, drops a stale sibling -----#
def _finalize_target(part, acc, asm, out_dir, session, delay):
    status = _append_rna(part, acc, asm, session)
    merged = status == "merged"
    dest = target_dest(out_dir, acc, asm, with_rna=merged)
    part.replace(dest)
    if merged:
        # Drop any pre-migration CDS-only file so one per accession survives.
        # One-directional by design: a degraded result never deletes a merged
        # sibling, since stale merged data beats it and self-heals via the cache.
        target_dest(out_dir, acc, asm, with_rna=False).unlink(missing_ok=True)
        print(f"    {dest.name}")
    time.sleep(delay)
    return status, dest


#----- Downloads CDS (API) + rna_from_genomic (FTP) per accession, merged into one target FASTA -----#
def download_targets(hits: dict[str, tuple[str, str] | None], out_dir, session,
                     delay: float, force: bool = False) -> tuple[dict[str, str | None], dict[str, str]]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Dedup by accession: several taxa can resolve to the same assembly, so we
    # download each accession once and map the file back to every taxon below.
    asm_by_acc: dict[str, str] = {}
    for hit in hits.values():
        if hit is not None:
            asm_by_acc.setdefault(hit[0], hit[1])

    # Only the merged name counts as a cache hit. A CDS-only file is a degraded
    # result, so we retry it — which also migrates pre-merge trees automatically.
    final_by_acc: dict[str, Path] = {}
    status_by_acc: dict[str, str] = {}
    todo: list[str] = []
    for acc, asm in asm_by_acc.items():
        merged = target_dest(out_dir, acc, asm, with_rna=True)
        if not force and merged.exists():
            final_by_acc[acc], status_by_acc[acc] = merged, "cached"
        else:
            todo.append(acc)

    if todo:
        parts = _download_cds_parts(todo, asm_by_acc, out_dir, session)
        for acc in todo:
            part = parts.get(acc)
            if part is None:
                # CDS never landed; there is nothing to merge onto and RNA is
                # never the sole content of a target FASTA.
                status_by_acc[acc] = "-"
                continue
            status_by_acc[acc], final_by_acc[acc] = _finalize_target(
                part, acc, asm_by_acc[acc], out_dir, session, delay)

    # Map every resolved taxon to its accession's file (None if it never landed),
    # so taxa sharing one accession all receive the same path. final_by_acc holds
    # only confirmed successes — a cache hit or a completed rename.
    paths: dict[str, str | None] = {}
    rna_status: dict[str, str] = {}
    for taxon, hit in hits.items():
        if hit is None:
            continue
        dest = final_by_acc.get(hit[0])
        paths[taxon] = str(dest) if dest is not None else None
        rna_status[taxon] = status_by_acc[hit[0]]
    return paths, rna_status


#----- Builds a requests session with retry-on-429/5xx and NCBI-etiquette headers -----#
def make_session(api_key: str | None = None, contact_email: str | None = None) -> requests.Session:
    s = requests.Session()
    user_agent = _DEFAULT_USER_AGENT
    if contact_email:
        user_agent = f"{user_agent} ({contact_email})"
    s.headers["User-Agent"] = user_agent
    if api_key:
        s.headers["api-key"] = api_key
    # urllib3 handles exponential backoff and respects Retry-After headers, so we don't roll our own.
    # POST is not in urllib3's default allowed_methods, which would leave
    # /genome/download — the run's single most expensive request — as the only
    # one with no retry. Replay is safe: the endpoint is a read.
    retry = Retry(total=3, backoff_factor=2.0, status_forcelist=_RETRIABLE_STATUS,
                  allowed_methods=Retry.DEFAULT_ALLOWED_METHODS | {"POST"})
    s.mount("https://", HTTPAdapter(max_retries=retry))
    return s


#----- POSTs accessions to /genome/download in batches, each ZIP into `zip_dir`; skips a failed batch -----#
def download_batch(accessions: list[str], session: requests.Session, zip_dir) -> list[Path]:
    zip_paths: list[Path] = []
    n_batches = (len(accessions) + _BATCH_SIZE - 1) // _BATCH_SIZE
    for b, start in enumerate(range(0, len(accessions), _BATCH_SIZE), 1):
        chunk = accessions[start : start + _BATCH_SIZE]
        print(f"  Batch {b}/{n_batches}: {len(chunk)} accessions ...", flush=True)
        try:
            r = session.post(
                f"{_API_BASE}/genome/download",
                json={"accessions": chunk, "include_annotation_type": ["CDS_FASTA"]},
                stream=True,
                timeout=600,
            )
            r.raise_for_status()
            # Stream to disk rather than buffering: one batch is up to 200 CDS sets.
            zip_paths.append(_stream_to_tempfile(r, ".zip", dir=zip_dir))
        except (requests.RequestException, OSError) as exc:
            # Its accessions never land, so the report records them as download_failed.
            print(f"    batch {b} failed, {len(chunk)} accessions skipped: "
                  f"{exc.__class__.__name__} — {exc}")
    return zip_paths


