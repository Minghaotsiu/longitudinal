import glob
import os
import re
import pandas as pd


ROOT = os.environ.get("UPSTREAM_DIR")
GTF = os.environ.get("SOYBEAN_GTF")
if not ROOT or not GTF:
    raise SystemExit("Set UPSTREAM_DIR and SOYBEAN_GTF before merging batches.")
batches = sorted(
    os.path.basename(path)
    for path in glob.glob(os.path.join(ROOT, "batch*"))
    if os.path.isdir(path)
)

lengths = {}
names = {}
biotypes = {}
with open(GTF) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        gene_id = re.search(r'gene_id "([^"]+)"', fields[8])
        gene_name = re.search(r'gene_name "([^"]+)"', fields[8])
        gene_biotype = re.search(r'gene_biotype "([^"]+)"', fields[8])
        if not gene_id:
            continue
        gid = gene_id.group(1)
        lengths[gid] = int(fields[4]) - int(fields[3]) + 1
        names[gid] = gene_name.group(1) if gene_name else gid.replace("GLYMA_", "Glyma.", 1)
        biotypes[gid] = gene_biotype.group(1) if gene_biotype else "unknown"


def merge_counts(pattern, output_path):
    frames = []
    for batch in batches:
        for file_path in sorted(glob.glob(os.path.join(ROOT, batch, "fenxi", pattern))):
            sample = (
                os.path.basename(file_path)
                .replace("_union_counts.txt", "")
                .replace("_1M", "")
            )
            frame = pd.read_csv(
                file_path,
                sep="\t",
                names=["Gene_ID", sample],
            )
            frame = frame[
                ~frame["Gene_ID"].astype(str).str.startswith("__")
            ]
            frames.append(frame.set_index("Gene_ID"))
    matrix = pd.concat(frames, axis=1).fillna(0).astype(int)
    matrix = matrix[sorted(matrix.columns)]
    matrix.index.name = "Gene_ID"
    matrix.to_csv(output_path, sep=",")
    return matrix


counts = merge_counts("*_union_counts.txt", os.path.join(ROOT, "ALL_8_raw_counts_matrix.csv"))
counts_1m = merge_counts("1M/*_union_counts.txt", os.path.join(ROOT, "ALL_8_1M_counts_matrix.csv"))

totals = counts.sum(axis=0)
fpkm = counts.astype(float).copy()
for gene_id in fpkm.index:
    fpkm.loc[gene_id] = (
        counts.loc[gene_id]
        / max(lengths.get(gene_id, 1), 1)
        / (totals / 1e9)
    )
fpkm.index.name = "Gene_ID"
fpkm.to_csv(os.path.join(ROOT, "ALL_8_FPKM.txt"), sep="\t")

symbol_map = pd.DataFrame(
    {
        "Gene_ID": list(counts.index),
        "Symbol": [names.get(gene_id, gene_id) for gene_id in counts.index],
        "Gene_biotype": [biotypes.get(gene_id, "unknown") for gene_id in counts.index],
    }
)
symbol_map.to_csv(os.path.join(ROOT, "ALL_8_gene_symbol_map.csv"), index=False)

protein_ids = symbol_map.loc[
    symbol_map["Gene_biotype"] == "protein_coding",
    "Gene_ID",
]
counts.loc[protein_ids].to_csv(
    os.path.join(ROOT, "ALL_8_protein_coding_counts_matrix.csv")
)
counts_1m.loc[protein_ids].to_csv(
    os.path.join(ROOT, "ALL_8_protein_coding_1M_counts_matrix.csv")
)

maps = []
for batch in batches:
    map_file = os.path.join(ROOT, batch, "fenxi", f"{batch}_FPKM_RNA_MAP.csv")
    if os.path.exists(map_file):
        maps.append(pd.read_csv(map_file))
if maps:
    sample_map = pd.concat(maps, ignore_index=True).dropna(subset=["Sample"])
    sample_map.to_csv(os.path.join(ROOT, "ALL_8_RNA_MAP.csv"), index=False)

print("batches:", ", ".join(batches))
print("raw counts:", counts.shape)
print("1M counts:", counts_1m.shape)
print("protein-coding counts:", counts.loc[protein_ids].shape)
