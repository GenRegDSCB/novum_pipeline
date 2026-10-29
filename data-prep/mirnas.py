"""mirnas.py — miRBase mature.fa: cache/download, index, name resolution, FASTA out."""
from __future__ import annotations

import gzip
from pathlib import Path

import requests

from config_writer import atomic_write

# miRBase current mature set. gz is auto-detected by magic bytes, so either the
# plain or the .gz URL works; adjust here if miRBase moves the path.
MIRBASE_URL = "https://www.mirbase.org/download/mature.fa"

# Species fallback order for a bare name. One spelling: prepare_inputs' CLI default
# is built from this, so changing it here actually changes the behaviour.
DEFAULT_PREFIX_CHAIN = ("hsa", "mmu")


class MatureUnavailable(RuntimeError):
    """miRBase unreachable and no cached mature.fa present — an external outage, not a bug."""

# A name is "bare" (no species prefix) when it starts with a miRNA stem.
_STEMS = ("mir-", "let-", "lin-")


#----- Parses mature.fa into {lowercased name: (canonical name, sequence)} -----#
def index_mature(path) -> dict[str, tuple[str, str]]:
    index: dict[str, tuple[str, str]] = {}
    name: str | None = None
    seq_parts: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                index[name.lower()] = (name, "".join(seq_parts))
            name = line[1:].split()[0]
            seq_parts = []
        elif (s := line.strip()):
            seq_parts.append(s)
    if name is not None:
        index[name.lower()] = (name, "".join(seq_parts))
    return index


#----- Exact name first, else prefix a bare mir-/let-/lin- stem; 3rd value flags the fallback -----#
def resolve_mirna(raw: str, index: dict[str, tuple[str, str]],
                  prefix_chain: tuple[str, ...] = DEFAULT_PREFIX_CHAIN):
    key = raw.strip().lower()
    if key in index:
        canonical, seq = index[key]
        return canonical, seq, False
    if key.startswith(_STEMS):
        for prefix in prefix_chain:
            hit = index.get(f"{prefix}-{key}")
            if hit:
                return hit[0], hit[1], True
    return None


#----- The query-FASTA path for a taxon: <rnas_dir>/<taxon>.fa (single owner) -----#
def query_fasta_path(rnas_dir, taxon_key: str) -> Path:
    return Path(rnas_dir) / f"{taxon_key}.fa"


#----- Writes one FASTA record per resolved miRNA into the taxon's query file -----#
def write_query_fasta(taxon_key: str, resolved, rnas_dir) -> Path:
    out = query_fasta_path(rnas_dir, taxon_key)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for canonical, seq in resolved:
            fh.write(f">{canonical}\n{seq}\n")
    return out


#----- Cached miRBase fetch (gzip auto-detected); an outage raises MatureUnavailable, not a traceback -----#
def ensure_mature_fa(mirnas_dir, session, force: bool = False, url: str = MIRBASE_URL) -> Path:
    out = Path(mirnas_dir) / "mature.fa"
    if out.exists() and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        resp = session.get(url, timeout=(10, 120))
        resp.raise_for_status()
    except requests.RequestException as e:
        raise MatureUnavailable(
            f"could not download miRBase mature.fa from {url}\n"
            f"  reason: {e}\n"
            f"miRBase is a single un-mirrored server and is periodically down.\n"
            f"Workaround: fetch mature.fa from any working network and place it at\n"
            f"  {out}\n"
            f"then re-run — the download is skipped whenever that file already exists."
        ) from e
    data = resp.content
    # Magic bytes only, not the URL suffix: with Content-Encoding: gzip requests
    # has already decompressed, and a .gz URL would then send plain FASTA here.
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    atomic_write(out, data)             # a partial file would read as a cache hit
    return out
