"""Dependency-light helpers shared across the Workflow/Scripts arms.

Kept deliberately free of third-party imports (no pandas): the RNAhybrid conda
env (`Workflow/Envs/rnahybrid.yaml`) ships no pandas, yet its scripts still need
`which_required` and `ensure_parent`, so importing this module must never pull a
heavy dependency. (The IntaRNA and postprocess envs do ship pandas.) Snakemake
puts each script's own directory on sys.path, so the sibling `from _common import
...` resolves when scripts run under the `script:` directive.
"""

import shutil
from pathlib import Path


#----- Resolves a required executable from PATH, trying each candidate name in order; the first candidate names the tool in the error -----#
def which_required(*candidates):
    for name in candidates:
        path = shutil.which(name)
        if path is not None:
            return path
    raise RuntimeError(f"Required executable not found in PATH: {candidates[0]}")


#----- Yields (header, sequence) per FASTA record: header is the text after '>', sequence is the residues with whitespace stripped -----#
def iter_fasta_records(path):
    header = None
    seq_parts = []

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        yield header, "".join(seq_parts)


#----- Yields each FASTA header (the text after '>') without accumulating sequences, for the
#      callers that only need record identity or a record count -----#
def iter_fasta_headers(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


#----- The FASTA join key: the first whitespace-delimited token of the header, which is
#      what RNAhybrid, IntaRNA and RNAcalibrate all report for a record -----#
def parse_header_id(header):
    return header.split()[0]


#----- The per-miRNA key. The same token, so rnahybrid.py splits query files to match
#      RNAcalibrate's column 1 -----#
query_key = parse_header_id


#----- Ensures the parent directory of `path` exists; returns `path` as a Path for chaining -----#
def ensure_parent(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
