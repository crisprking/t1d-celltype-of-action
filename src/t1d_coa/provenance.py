"""Append-only provenance log for every artifact the pipeline writes.

One entry per file: relative path, source URL or step description, UTC timestamp,
first 16 hex of sha256, byte size, free-text notes. The file is a markdown table
so it stays human-readable as the pipeline grows.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

from .config import PROJECT_ROOT, PROVENANCE

_HEADER = (
    "# Data Provenance\n\n"
    "| Path | Source | Date (UTC) | sha256 | Bytes | Notes |\n"
    "|------|--------|------------|--------|-------|-------|\n"
)


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    """Stream a file through sha256 in 1 MiB chunks."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda: f.read(chunk), b""):
            h.update(buf)
    return h.hexdigest()


def log_artifact(path: Path, source: str, notes: str = "") -> None:
    """Append one row to PROVENANCE.md. Initializes the file if absent."""
    if not PROVENANCE.exists():
        PROVENANCE.parent.mkdir(parents=True, exist_ok=True)
        PROVENANCE.write_text(_HEADER)

    digest = sha256(path)
    size = path.stat().st_size
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    rel = path.relative_to(PROJECT_ROOT) if PROJECT_ROOT in path.parents else path

    src_safe = source.replace("|", "\\|")
    notes_safe = notes.replace("|", "\\|")
    with PROVENANCE.open("a") as f:
        f.write(
            f"| `{rel}` | {src_safe} | {now} | `{digest[:16]}…` | "
            f"{size:,} | {notes_safe} |\n"
        )
