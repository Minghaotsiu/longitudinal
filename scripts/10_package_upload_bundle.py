from __future__ import annotations

import hashlib
import os
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(
    os.environ.get(
        "PROJECT_ROOT",
        Path(__file__).resolve().parents[1],
    )
).resolve()
OUT = Path(
    os.environ.get(
        "SCIENTIFIC_DATA_OUTPUT_DIR",
        ROOT / "outputs" / "doudou_scientific_data",
    )
).resolve()
UPLOAD = Path(
    os.environ.get(
        "UPLOAD_BUNDLE_DIR",
        OUT / "submission" / "upload_ready",
    )
).resolve()
SOURCE = OUT / "source_data"
FIGURES = OUT / "figures"
UPLOAD.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


crosswalk = pd.read_csv(OUT / "tables" / "scientific_data_sample_crosswalk.csv")
delivery_alias_map = {
    "RJYY-20260923-L-02": "Delivery01_20260923",
    "RJYY-20261015-L-01": "Delivery02_20261015",
}
crosswalk["DeliveryAlias"] = crosswalk["SourceDelivery"].map(delivery_alias_map)
if crosswalk["DeliveryAlias"].isna().any():
    missing = sorted(crosswalk.loc[crosswalk["DeliveryAlias"].isna(), "SourceDelivery"].unique())
    raise ValueError(f"Unmapped source delivery values: {missing}")

crosswalk["RunAlias"] = [
    f"RUN{index:02d}_{sample_id}"
    for index, sample_id in enumerate(crosswalk["CorrectedSampleID"], start=1)
]

if crosswalk["RunAlias"].duplicated().any():
    raise ValueError("RunAlias values are not unique")

# The source deliveries contain identically named FASTQ files. A public
# run identifier and delivery alias are therefore required to keep all 34
# runs distinguishable without altering the original filenames.
filename_collisions = crosswalk.groupby("OriginalR1Filename")["DeliveryAlias"].nunique()
if (filename_collisions > 1).any() and crosswalk["DeliveryAlias"].isna().any():
    raise ValueError("Colliding FASTQ filenames cannot be disambiguated")

sample_metadata = pd.DataFrame(
    {
        "sample_name": crosswalk["CorrectedSampleID"],
        "run_alias": crosswalk["RunAlias"],
        "title": crosswalk["CorrectedSampleID"],
        "organism": "Glycine max",
        "cultivar": "Williams 82",
        "tissue": "seed",
        "developmental_stage": crosswalk["Day"],
        "sampling_method": crosswalk["SamplingMethod"],
        "seed_id": crosswalk["SeedID"],
        "library_strategy": "RNA-Seq",
        "library_layout": "PAIRED",
        "instrument_model": "DNBSEQ-T7",
        "read_length": 150,
        "target_output_gb": 2,
        "GSA_sample_accession": "",
        "GSA_experiment_accession": "",
        "GSA_run_accession": "",
    }
)
sample_metadata.to_csv(UPLOAD / "GSA_sample_metadata_template.tsv", sep="\t", index=False)

raw_template = crosswalk[
    [
        "RunAlias",
        "CorrectedSampleID",
        "SamplingMethod",
        "SeedID",
        "Day",
        "SourceDelivery",
        "DeliveryAlias",
        "OriginalR1Filename",
        "OriginalR2Filename",
    ]
].copy()
raw_template["Upload_FASTQ_R1_Filename"] = (
    raw_template["RunAlias"] + "_R1.fastq.gz"
)
raw_template["Upload_FASTQ_R2_Filename"] = (
    raw_template["RunAlias"] + "_R2.fastq.gz"
)
raw_template["Local_FASTQ_R1_Path"] = ""
raw_template["Local_FASTQ_R2_Path"] = ""
raw_template["MD5_R1"] = ""
raw_template["MD5_R2"] = ""
raw_template["Uploaded_Path"] = ""
raw_template.columns = [
    "RunAlias",
    "SampleID",
    "SamplingMethod",
    "SeedID",
    "Day",
    "SourceDelivery",
    "DeliveryAlias",
    "Local_FASTQ_R1_Filename",
    "Local_FASTQ_R2_Filename",
    "Upload_FASTQ_R1_Filename",
    "Upload_FASTQ_R2_Filename",
    "Local_FASTQ_R1_Path",
    "Local_FASTQ_R2_Path",
    "MD5_R1",
    "MD5_R2",
    "Repository_Path",
]
raw_template.to_csv(UPLOAD / "FASTQ_upload_and_MD5_template.tsv", sep="\t", index=False)

processed_files = sorted(SOURCE.rglob("*"))
processed_files = [path for path in processed_files if path.is_file()]
manifest_rows = []
for source_path in processed_files:
    relative_path = source_path.relative_to(OUT)
    target_path = UPLOAD / "processed" / relative_path
    target_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, target_path)
    manifest_rows.append(
        {
            "File": relative_path.as_posix(),
            "Description": (
                "Source data for main figures"
                if relative_path.parts[0] == "source_data"
                else "Legacy selected-5 analysis table"
            ),
            "SHA256": sha256(source_path),
            "Bytes": source_path.stat().st_size,
        }
    )

for stem in [
    "Figure1_sampling_workflow",
    "Figure2_technical_QC",
    "Figure3_expression_structure",
    "Figure4_MN_Bulk_comparison",
    "Figure5_pseudotime_comparison",
    "Figure6_dataset_supported_analyses",
]:
    for extension in ("pdf", "svg", "tiff", "png"):
        source_path = FIGURES / f"{stem}.{extension}"
        relative_path = source_path.relative_to(OUT)
        target_path = UPLOAD / "figures" / source_path.name
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target_path)
        manifest_rows.append(
            {
                "File": relative_path.as_posix(),
                "Description": "Main figure output",
                "SHA256": sha256(source_path),
                "Bytes": source_path.stat().st_size,
            }
        )

pd.DataFrame(manifest_rows).to_csv(
    UPLOAD / "processed_file_manifest_sha256.csv",
    index=False,
)

(UPLOAD / "DATASET_README.md").write_text(
    """# Longitudinal soybean seed-germination transcriptome dataset

## Summary

Bulk RNA-seq data from longitudinal active-microneedle sampling of five soybean seeds at Day0-Day3 and from 14 independent destructive bulk samples.

## Files

- `GSA_sample_metadata_template.tsv`: sample metadata for GSA/BioSample registration.
- `FASTQ_upload_and_MD5_template.tsv`: unique run aliases, delivery aliases,
  proposed upload filenames and MD5 fields. The two source deliveries contain
  identically named files; keep `SourceDelivery`/`DeliveryAlias` visible when
  resolving local paths and never collapse 34 runs into 17 filenames.
- `processed/`: processed expression matrices, QC tables, figure source data and legacy analysis tables.
- `figures/`: main figures in PDF, SVG, TIFF and PNG.
- `processed_file_manifest_sha256.csv`: file-level SHA256 checksums.

## Variables and units

- Day: Day0, Day1, Day2 or Day3 after the start of imbibition.
- SamplingMethod: MN (active microneedle) or Bulk (independent destructive sampling).
- RunAlias: one of 34 unique run identifiers. Source deliveries contain
  identically named FASTQ files, so the delivery alias and run alias must be
  retained when resolving local files or uploading to a repository.
- Count data: integer gene-level counts.
- Normalized data: DESeq2 normalized counts and variance-stabilizing transformed expression.
- Read length: 150 bp.

## Methods and provenance

Samples were collected as described in the manuscript. Reads were processed with fastp, STAR 2.7.11b and htseq-count 2.0.5 against Glycine max v2.1 and the Glycine_max_v2.1.63.chr.gtf annotation.

## Access and licence

The public repository accession and release date should be inserted after GSA/OMIX registration. The planned data licence is CC BY 4.0 and the planned code licence is MIT, subject to final author confirmation.
""",
    encoding="utf-8",
)

print(UPLOAD)
