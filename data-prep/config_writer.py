"""config_writer.py — surgically rewrite queries:/targets: blocks in config.yaml,
preserving every other line (comments, params, formatting)."""
from __future__ import annotations

import os
import re
import secrets
import stat
from pathlib import Path


#----- Swaps one top-level key's entries for `mapping`, carrying the block's comments through -----#
def replace_block(lines: list[str], key: str, mapping: dict[str, str]) -> list[str]:
    key_re = re.compile(rf"^{re.escape(key)}:\s*$")          # top-level, anchored
    i = next((idx for idx, ln in enumerate(lines) if key_re.match(ln)), None)
    if i is None:
        raise ValueError(f"expected top-level '{key}:' in config")

    # The block runs to the next top-level KEY. A column-0 comment must not end it, or the
    # old entries below it survive as silent duplicate mapping keys (last wins).
    end = i + 1
    for j in range(i + 1, len(lines)):
        ln = lines[j]
        if ln.strip() == "" or ln[:1] in (" ", "\t") or ln.lstrip().startswith("#"):
            end = j + 1
            continue
        break

    # Entries are replaced wholesale; the block's blanks and comments are re-emitted below them.
    kept = [ln for ln in lines[i + 1 : end]
            if ln.strip() == "" or ln.lstrip().startswith("#")]
    new_entries = [f"  {k}: {v}" for k, v in mapping.items()]
    return lines[: i + 1] + new_entries + kept + lines[end:]


#----- Writes text or bytes via a same-directory temp file + os.replace, so a crash never leaves a
#      partial file. 0o666 lets umask set a new file's mode; an existing file keeps its own -----#
def atomic_write(path, data: str | bytes) -> None:
    path = Path(path)
    tmp = path.with_name(f"{path.name}.{secrets.token_hex(8)}.tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
    try:
        with os.fdopen(fd, "wb" if isinstance(data, bytes) else "w") as fh:
            if path.exists():
                os.chmod(tmp, stat.S_IMODE(path.stat().st_mode))
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


#----- Rewrites queries: and targets: atomically -----#
def update_path_blocks(config_path, queries: dict[str, str], targets: dict[str, str]) -> None:
    path = Path(config_path)
    lines = path.read_text().splitlines()
    lines = replace_block(lines, "queries", queries)
    lines = replace_block(lines, "targets", targets)
    atomic_write(path, "\n".join(lines) + "\n")
