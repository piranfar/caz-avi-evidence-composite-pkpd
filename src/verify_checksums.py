"""Compute SHA-256 checksums for the tracked repository files.

Build artefacts are skipped so the manifest is stable across runs; the
reproducibility workflow compares it against the committed copy.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE = {"SHA256SUMS.txt"}
EXCLUDE_DIRS = {".git", "__pycache__", ".venv", ".ipynb_checkpoints"}
EXCLUDE_SUFFIXES = {".pyc", ".pyo"}

for path in sorted(p for p in ROOT.rglob("*") if p.is_file()
                    and p.name not in EXCLUDE
                    and p.suffix not in EXCLUDE_SUFFIXES
                    and not EXCLUDE_DIRS & set(p.parts)):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    print(f"{h.hexdigest()}  {path.relative_to(ROOT).as_posix()}")
