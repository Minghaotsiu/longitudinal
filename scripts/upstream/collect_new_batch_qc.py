import glob
import os
import pandas as pd


ROOT = os.environ.get("UPSTREAM_DIR")
if not ROOT:
    raise SystemExit("Set UPSTREAM_DIR before collecting QC metrics.")
batches = sorted(
    os.path.basename(path)
    for path in glob.glob(os.path.join(ROOT, "batch*"))
    if os.path.isdir(path)
)


def parse_log(file_path):
    values = {}
    with open(file_path) as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            values[key.strip()] = value.strip()
    return values


rows = []
for batch in batches:
    pattern = os.path.join(
        ROOT,
        batch,
        "trim_fastqc_*",
        "star_out_*",
        "log_final_out",
        "*Log.final.out",
    )
    for file_path in glob.glob(pattern):
        sample = os.path.basename(file_path).replace("Log.final.out", "")
        values = parse_log(file_path)

        def number(key):
            value = values.get(key, "")
            value = value.replace("%", "").replace(",", "")
            try:
                return float(value)
            except ValueError:
                return float("nan")

        rows.append(
            {
                "Sample": sample,
                "InputReads": int(number("Number of input reads")),
                "AvgLen": number("Average input read length"),
                "Unique%": number("Uniquely mapped reads %"),
                "Multi%": number("% of reads mapped to multiple loci"),
                "TooShort%": number("% of reads unmapped: too short"),
                "Mismatch%": number("Mismatch rate per base, %"),
                "MappedLen": number("Average mapped length"),
            }
        )

quality = pd.DataFrame(rows).set_index("Sample")

raw_frames = {}
one_m_frames = {}
for batch in batches:
    for file_path in glob.glob(
        os.path.join(ROOT, batch, "fenxi", "*_union_counts.txt")
    ):
        sample = os.path.basename(file_path).replace("_union_counts.txt", "")
        raw_frames[sample] = pd.read_csv(
            file_path, sep="\t", names=["g", "c"]
        )
    for file_path in glob.glob(
        os.path.join(ROOT, batch, "fenxi", "1M", "*_union_counts.txt")
    ):
        sample = os.path.basename(file_path).replace("_1M_union_counts.txt", "")
        one_m_frames[sample] = pd.read_csv(
            file_path, sep="\t", names=["g", "c"]
        )


def count_stats(frame):
    frame = frame[~frame["g"].astype(str).str.startswith("__")]
    values = frame["c"].values
    total = values.sum()
    top = sorted(values, reverse=True)[:100]
    return {
        "counts": int(total),
        "genes1": int((values >= 1).sum()),
        "genes5": int((values >= 5).sum()),
        "genes10": int((values >= 10).sum()),
        "top100pct": round(100.0 * sum(top) / total, 1) if total else 0.0,
    }


for sample, frame in raw_frames.items():
    for key, value in count_stats(frame).items():
        quality.loc[sample, key] = value
for sample, frame in one_m_frames.items():
    frame = frame[~frame["g"].astype(str).str.startswith("__")]
    quality.loc[sample, "genes1_1M"] = int((frame["c"].values >= 1).sum())

quality = quality.sort_values("Unique%", ascending=False)
columns = [
    "InputReads",
    "Unique%",
    "Multi%",
    "TooShort%",
    "Mismatch%",
    "counts",
    "genes1",
    "genes5",
    "genes10",
    "genes1_1M",
    "top100pct",
]
quality[columns].to_csv(
    os.path.join(ROOT, "QC_all8.tsv"),
    sep="\t",
)
print(quality[columns].to_string())
