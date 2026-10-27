from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


def attr_value(attributes: str, key: str) -> str | None:
    match = re.search(rf'{re.escape(key)} "([^"]+)"', attributes)
    return match.group(1) if match else None


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: parse_gtf_symbols.py GTF OUTPUT")
    gtf_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    rows: dict[str, tuple[str, str]] = {}
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_id = attr_value(fields[8], "gene_id")
            if not gene_id:
                continue
            gene_name = attr_value(fields[8], "gene_name")
            biotype = attr_value(fields[8], "gene_biotype") or ""
            if gene_name:
                symbol = gene_name
            elif gene_id.startswith("GLYMA_"):
                symbol = gene_id.replace("GLYMA_", "Glyma.", 1)
            else:
                symbol = gene_id
            rows[gene_id] = (symbol, biotype)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gene_ID", "Symbol", "Gene_has_GTF_name", "Gene_biotype"])
        for gene_id in sorted(rows):
            symbol, biotype = rows[gene_id]
            writer.writerow([
                gene_id,
                symbol,
                "yes" if symbol != gene_id.replace("GLYMA_", "Glyma.", 1) else "no",
                biotype,
            ])
    print(f"genes={len(rows)}")
    print(f"output={output_path}")


if __name__ == "__main__":
    main()
