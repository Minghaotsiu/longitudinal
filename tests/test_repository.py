#!/usr/bin/env python3
"""Fail if the public repository accidentally contains large/private files."""

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALLOWED_LARGE_SUFFIXES = {".png"}
FORBIDDEN_SUFFIXES = {
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".bam",
    ".sam",
    ".bai",
    ".tiff",
    ".tif",
    ".pptx",
    ".docx",
    ".zip",
}
FORBIDDEN_NAME_PARTS = (
    "password",
    "credential",
    "private_key",
    "id_rsa",
    "id_ed25519",
)
SECRET_PATTERNS = [
    re.compile(r"BEGIN (?:OPENSSH|RSA|EC) PRIVATE KEY"),
    re.compile(r"AccessKeySecret\s*[:=]"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
]


def main() -> None:
    problems: list[str] = []
    for path in ROOT.rglob("*"):
        if (
            not path.is_file()
            or ".git" in path.parts
            or "__pycache__" in path.parts
        ):
            continue
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            problems.append(f"Forbidden file type: {path.relative_to(ROOT)}")
        if any(part in path.name.lower() for part in FORBIDDEN_NAME_PARTS):
            problems.append(f"Sensitive filename: {path.relative_to(ROOT)}")
        if path.suffix.lower() in {".png", ".jpg", ".jpeg"}:
            continue
        if path.stat().st_size > 5_000_000:
            problems.append(f"Unexpected file larger than 5 MB: {path.relative_to(ROOT)}")
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for pattern in SECRET_PATTERNS:
            if pattern.search(text):
                problems.append(f"Possible secret in {path.relative_to(ROOT)}")

    with (ROOT / "metadata" / "scientific_data_sample_crosswalk.csv").open(
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 34:
        problems.append(f"Expected 34 metadata rows, found {len(rows)}")

    if problems:
        raise SystemExit("\n".join(problems))
    print("Repository safety and metadata checks passed.")


if __name__ == "__main__":
    main()

