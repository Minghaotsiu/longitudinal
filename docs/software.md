# Software and reference resources

## Command-line tools

| Tool | Version | Purpose |
|---|---:|---|
| FastQC | record in run log | raw and trimmed read QC |
| cutadapt | record in run log | adapter and quality trimming |
| BBTools `repair.sh` | record in run log | paired-read repair |
| STAR | 2.7.11b | read alignment |
| samtools | record in run log | SAM/BAM processing |
| Python | 3 | upstream and report scripts |
| R | 4.3.1 | statistical analysis and figures |
| DESeq2 | 1.40.2 | differential expression |

## Reference resources

The repository does not redistribute the soybean reference. Set:

- `SOYBEAN_GTF`
- `STAR_INDEX`
- `SOYBEAN_EXON_BED`
- `SOYBEAN_GENE_BED`
- `MT_GENE_FILE`
- `RP_GENE_FILE`

The analysis used `Glycine_max_v2.1.63.chr.gtf`.

