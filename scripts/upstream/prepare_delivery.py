#!/usr/bin/env python3
"""Collect paired FASTQ files into one batch without duplicating raw data."""

from __future__ import annotations

import argparse
import glob
import os
import re


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", help="Directory containing delivered FASTQ files.")
    parser.add_argument("destination", help="Batch working directory.")
    parser.add_argument("expected", type=int, help="Expected number of paired runs.")
    parser.add_argument(
        "--pattern",
        default=r"(d[0-3]-\d+-1)",
        help="Regular expression used to derive a sample ID from each R1 file.",
    )
    parser.add_argument(
        "--r1-suffix",
        default="_combined_R1.fastq.gz",
        help="Suffix used to identify and pair R1 files.",
    )
    args = parser.parse_args()

    samples: dict[str, tuple[str, str]] = {}
    for r1 in glob.glob(
        os.path.join(args.source, "**", f"*{args.r1_suffix}"),
        recursive=True,
    ):
        match = re.search(args.pattern, os.path.basename(r1))
        if not match:
            continue
        sample = match.group(1)
        r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
        if os.path.exists(r2):
            samples[sample] = (r1, r2)

    if len(samples) != args.expected:
        raise SystemExit(
            f"Expected {args.expected} samples in {args.source}, found {len(samples)}"
        )

    os.makedirs(args.destination, exist_ok=True)
    for sample, (r1, r2) in sorted(samples.items()):
        for source, read in ((r1, "R1"), (r2, "R2")):
            target = os.path.join(args.destination, f"{sample}_{read}.fastq.gz")
            if not os.path.exists(target):
                try:
                    os.link(source, target)
                except OSError:
                    os.symlink(source, target)

    with open(
        os.path.join(args.destination, "samples.tsv"),
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write("sample\tR1\tR2\n")
        for sample in sorted(samples):
            handle.write(
                f"{sample}\t{sample}_R1.fastq.gz\t{sample}_R2.fastq.gz\n"
            )

    print("batch prepared:", args.destination)
    print(", ".join(sorted(samples)))


if __name__ == "__main__":
    main()
