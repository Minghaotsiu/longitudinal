#!/usr/bin/env python3
"""Validate the minimum inputs needed by the downstream analysis."""

from __future__ import annotations

import csv
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("PROCESSED_DATA_DIR", ROOT / "data" / "doudou_all34"))

REQUIRED_DATA = [
    "ALL_34_raw_counts_matrix.csv",
    "ALL_34_1M_counts_matrix.csv",
    "ALL_34_sample_metadata.csv",
    "ALL_34_QC_metrics.csv",
    "soybean_gene_symbol_map.csv",
]
REQUIRED_METADATA = [
    ROOT / "metadata" / "scientific_data_sample_crosswalk.csv",
    ROOT / "metadata" / "sample_metadata_34.tsv",
    ROOT / "metadata" / "qc_metrics_corrected_long.csv",
]


def check_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        raise ValueError(f"Empty file: {path}")


def main() -> None:
    for name in REQUIRED_DATA:
        check_file(DATA / name)
    for path in REQUIRED_METADATA:
        check_file(path)

    with (ROOT / "metadata" / "scientific_data_sample_crosswalk.csv").open(
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        crosswalk = list(csv.DictReader(handle))
    if len(crosswalk) != 34:
        raise ValueError(f"Expected 34 crosswalk rows, found {len(crosswalk)}")
    sample_ids = [row["CorrectedSampleID"] for row in crosswalk]
    if len(set(sample_ids)) != 34:
        raise ValueError("CorrectedSampleID is not unique")

    with (ROOT / "metadata" / "sample_metadata_34.tsv").open(
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    if len(metadata) != 34:
        raise ValueError(f"Expected 34 metadata rows, found {len(metadata)}")

    counts = [row["SamplingMethod"] for row in crosswalk]
    if counts.count("MN") != 20 or counts.count("Bulk") != 14:
        raise ValueError("Expected 20 MN and 14 Bulk samples")

    print(f"Data directory: {DATA}")
    print("Input validation passed.")


if __name__ == "__main__":
    main()

