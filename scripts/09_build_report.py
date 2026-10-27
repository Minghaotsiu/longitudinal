from __future__ import annotations

import base64
import os
from html import escape
from io import BytesIO
from pathlib import Path

import pandas as pd
from PIL import Image


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
FIG = OUT / "figures"
SUP = FIG / "supplementary"
SOURCE = OUT / "source_data"


def image_data_uri(path: Path, thumbnail: bool = False) -> str:
    if not thumbnail:
        payload = path.read_bytes()
        mime = "image/png" if path.suffix.lower() == ".png" else "image/jpeg"
        return f"data:{mime};base64,{base64.b64encode(payload).decode('ascii')}"

    with Image.open(path) as image:
        image.thumbnail((1800, 1800))
        if image.mode not in ("RGB", "L"):
            image = image.convert("RGB")
        buffer = BytesIO()
        image.save(buffer, format="JPEG", quality=86, optimize=True)
    return f"data:image/jpeg;base64,{base64.b64encode(buffer.getvalue()).decode('ascii')}"


def table(df: pd.DataFrame, classes: str = "data-table") -> str:
    return df.to_html(
        index=False,
        escape=True,
        classes=classes,
        border=0,
        na_rep="NA",
    )


def figure_card(
    title: str,
    image: str,
    pdf: str,
    caption: str,
    figure_id: str,
) -> str:
    image_path = OUT / image
    return f"""
    <figure class="main-figure" id="{figure_id}">
      <div class="figure-head">
        <div>
          <div class="eyebrow">Main figure</div>
          <h3>{escape(title)}</h3>
        </div>
        <a class="button" href="{escape(pdf)}" target="_blank">Open PDF</a>
      </div>
      <a href="{escape(image)}" target="_blank">
        <img src="{image_data_uri(image_path)}" alt="{escape(title)}">
      </a>
      <figcaption>{caption}</figcaption>
    </figure>
    """


de_summary = pd.read_csv(OUT.parent / "doudou_selected5_v2" / "tables" / "DE_summary_selected_n5.csv")
de_summary["total"] = de_summary["up"] + de_summary["down"]
de_summary_display = de_summary.rename(
    columns={
        "comparison": "Comparison",
        "up": "Up",
        "down": "Down",
        "tested": "Tested genes",
        "total": "Significant total",
    }
)

candidate = pd.read_csv(
    OUT.parent / "doudou_selected5_v2" / "tables" / "candidate_gene_one_way_ANOVA.csv"
)
candidate_display = candidate.rename(
    columns={
        "Symbol": "Gene",
        "Category": "Annotation",
        "ANOVA_F": "ANOVA F",
        "ANOVA_p": "ANOVA P",
        "n_per_group": "n per Day",
    }
)

mfuzz = pd.read_csv(
    OUT.parent / "doudou_selected5_v2" / "tables" / "Mfuzz_cluster_interpretation.csv"
)
mfuzz_display = mfuzz.rename(
    columns={
        "Cluster": "Cluster",
        "n_genes": "Genes",
        "Peak_day": "Peak Day",
        "Trough_day": "Trough Day",
        "Relative_pattern": "Expression pattern",
        "Top_genes": "Representative genes",
        "Top_pathways": "Representative terms",
        "Biological_interpretation": "Description",
    }
)[
    [
        "Cluster",
        "Genes",
        "Peak Day",
        "Trough Day",
        "Expression pattern",
        "Representative genes",
        "Representative terms",
        "Description",
    ]
]

pseudotime_genes = pd.read_csv(
    OUT.parent / "doudou_selected5_v2" / "tables" / "pseudotime_gene_interpretation.csv"
)
pseudotime_gene_display = pseudotime_genes.rename(
    columns={
        "Symbol": "Gene",
        "Function": "Function",
        "Spearman_rho": "Spearman rho",
        "Direction": "Direction",
        "Interpretation": "Interpretation",
    }
)[["Gene", "Function", "Spearman rho", "Direction", "Interpretation"]]

qc_stats = pd.read_csv(SOURCE / "Figure2_ANOVA_statistics.csv")
qc_display = qc_stats.rename(
    columns={
        "Variable": "Metric",
        "PanelLabel": "Display label",
        "PValue": "ANOVA P",
        "Displayed": "Shown in Figure 2",
    }
)[["Metric", "Display label", "ANOVA P", "Shown in Figure 2"]]
qc_display["ANOVA P"] = qc_display["ANOVA P"].map(lambda value: f"{value:.3e}")

pt_stats = pd.read_csv(SOURCE / "Figure5_pseudotime_vs_day_source_data.csv")
pt_display = pt_stats.rename(
    columns={"Method": "Method", "SpearmanRho": "Spearman rho versus recorded Day"}
)

gallery_groups = [
    (
        "Differential expression across the germination time course",
        [
            ("S03_DE_Day1_vs_Day0.png", "Day1 versus Day0"),
            ("S04_DE_Day2_vs_Day0.png", "Day2 versus Day0"),
            ("S05_DE_Day3_vs_Day0.png", "Day3 versus Day0"),
            ("S06_DE_Day2_vs_Day1.png", "Day2 versus Day1"),
            ("S07_DE_Day3_vs_Day1.png", "Day3 versus Day1"),
            ("S08_DE_Day3_vs_Day2.png", "Day3 versus Day2"),
        ],
    ),
    (
        "Candidate developmental genes",
        [
            ("S09_candidate_gene_ANOVA.png", "Candidate-gene ANOVA"),
            ("S10_candidate_gene_Tukey.png", "Pairwise Tukey comparisons"),
        ],
    ),
    (
        "Temporal clustering and functional enrichment",
        [
            ("S11_Mfuzz_clusters.png", "Mfuzz clusters"),
            ("S12_Mfuzz_trajectories.png", "Mfuzz trajectories"),
            ("S13_Mfuzz_heatmap.png", "Mfuzz heatmap"),
            ("S15_Mfuzz_complete.png", "Mfuzz combined figure"),
            ("S14_Mfuzz_enrichment.png", "Mfuzz cluster enrichment"),
            ("S24_GO_KEGG_bubble.png", "GO/KEGG enrichment"),
            ("S25_GO_KEGG_gene_links.png", "Enriched terms and genes"),
            ("S26_GO_KEGG_sankey.png", "GO/KEGG Sankey"),
        ],
    ),
    (
        "Pseudotime and method sensitivity",
        [
            ("S16_pseudotime_UMAP.png", "Pseudotime in UMAP space"),
            ("S17_pseudotime_by_Day.png", "Pseudotime distribution by Day"),
            ("S18_pseudotime_gene_trends.png", "Pseudotime-associated gene trends"),
            ("S19_TSCAN_state_tree.png", "TSCAN state tree"),
            ("S20_TSCAN_state_UMAP.png", "TSCAN states in UMAP"),
            ("S21_SCORPIUS_trajectory.png", "SCORPIUS trajectory"),
            ("S22_pseudotime_method_comparison.png", "Comparison of pseudotime methods"),
            ("S23_TSCAN_pseudotime_by_Day.png", "TSCAN pseudotime by Day"),
        ],
    ),
]

gallery_html = ""
for heading, items in gallery_groups:
    cards = []
    for filename, title in items:
        path = SUP / filename
        if path.exists():
            cards.append(
                f"""
                <a class="gallery-card" href="{escape(path.relative_to(OUT).as_posix())}" target="_blank">
                  <img src="{image_data_uri(path, thumbnail=True)}" alt="{escape(title)}">
                  <span>{escape(title)}</span>
                </a>
                """
            )
    gallery_html += f"""
    <section class="gallery-group">
      <h4>{escape(heading)}</h4>
      <div class="gallery-grid">{''.join(cards)}</div>
    </section>
    """

feature_cards = [
    ("34", "sequencing runs"),
    ("20", "longitudinal MN samples"),
    ("14", "destructive Bulk samples"),
    ("6", "pairwise DE contrasts"),
    ("6", "temporal Mfuzz clusters"),
    ("3", "pseudotime methods"),
]
feature_html = "".join(
    f'<div class="metric"><strong>{value}</strong><span>{label}</span></div>'
    for value, label in feature_cards
)

main_figures = "".join(
    [
        figure_card(
            "Figure 1 | Data-acquisition workflow",
            "figures/Figure1_sampling_workflow.png",
            "figures/Figure1_sampling_workflow.pdf",
            "Five seeds are followed longitudinally with active microneedles from Day0 to Day3, while independent destructive bulk samples provide the cross-sectional comparator. The lower workflow traces material release, RNA extraction, library preparation, sequencing and construction of reusable expression and Seed-by-Day metadata tables.",
            "figure-1",
        ),
        figure_card(
            "Figure 2 | Technical metrics without significant MN–Bulk differences",
            "figures/Figure2_technical_QC.png",
            "figures/Figure2_technical_QC.pdf",
            "Displayed metrics have one-way ANOVA P ≥ 0.05 between sampling methods. The complete statistical table, including the input-read and clean-read metrics that did not meet this display criterion, is reported below.",
            "figure-2",
        ),
        figure_card(
            "Figure 3 | Global expression structure",
            "figures/Figure3_expression_structure.png",
            "figures/Figure3_expression_structure.pdf",
            "PCA, UMAP and t-SNE. Samples are circles, and colour encodes MN versus Bulk only. These panels are descriptive.",
            "figure-3",
        ),
        figure_card(
            "Figure 4 | Microneedle versus destructive bulk expression profiles",
            "figures/Figure4_MN_Bulk_comparison.png",
            "figures/Figure4_MN_Bulk_comparison.pdf",
            "Group-mean expression concordance, within-method sample similarity and the top 80 variable genes. The panels describe concordance and heterogeneity rather than formal equivalence.",
            "figure-4",
        ),
        figure_card(
            "Figure 5 | Data-derived ordering compared with recorded sampling Day",
            "figures/Figure5_pseudotime_comparison.png",
            "figures/Figure5_pseudotime_comparison.pdf",
            "Pseudotime is inferred from expression and is not equivalent to experimental time. The comparison quantifies agreement with recorded Day and disagreement among Slingshot, TSCAN and SCORPIUS.",
            "figure-5",
        ),
        figure_card(
            "Figure 6 | Developmental and reserve-mobilization candidate genes",
            "figures/Figure6_dataset_supported_analyses.png",
            "figures/Figure6_dataset_supported_analyses.pdf",
            "The finalized candidate-gene figure from the earlier selected-5 analysis is used without adding a differential-expression count summary. The panel supports reuse of the processed expression matrix; it is not used to establish causal regulation.",
            "figure-6",
        ),
    ]
)


def extract_figure(figure_id: str) -> str:
    marker = f'<figure class="main-figure" id="{figure_id}">'
    start = main_figures.find(marker)
    if start < 0:
        raise ValueError(figure_id)
    end = main_figures.find('<figure class="main-figure"', start + len(marker))
    return main_figures[start:] if end < 0 else main_figures[start:end]


html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Longitudinal soybean transcriptome data report</title>
<style>
:root {{
  --ink: #24313a;
  --muted: #66757f;
  --blue: #356c82;
  --teal: #78a58f;
  --sand: #c57c5a;
  --paper: #f7f8f5;
  --line: #dce2e2;
  --white: #ffffff;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  color: var(--ink);
  background: var(--paper);
  font-family: "Arial", "Helvetica Neue", sans-serif;
  line-height: 1.62;
}}
.hero {{
  padding: 64px 7vw 52px;
  color: white;
  background:
    radial-gradient(circle at 78% 20%, rgba(197,124,90,.40), transparent 25%),
    linear-gradient(135deg, #254653, #356c82 58%, #527d75);
}}
.hero h1 {{
  max-width: 1100px;
  margin: 10px 0 18px;
  font-size: clamp(32px, 5vw, 62px);
  line-height: 1.06;
  letter-spacing: -0.025em;
}}
.hero p {{ max-width: 900px; font-size: 18px; color: #e9f1ef; }}
.kicker {{ text-transform: uppercase; letter-spacing: .16em; font-size: 12px; font-weight: 700; }}
.metric-grid {{
  display: grid;
  grid-template-columns: repeat(6, minmax(110px, 1fr));
  gap: 12px;
  margin-top: 32px;
}}
.metric {{
  padding: 18px;
  border: 1px solid rgba(255,255,255,.22);
  border-radius: 14px;
  background: rgba(255,255,255,.10);
  backdrop-filter: blur(4px);
}}
.metric strong {{ display: block; font-size: 30px; line-height: 1; }}
.metric span {{ display: block; margin-top: 8px; color: #dfe9e7; font-size: 12px; }}
.layout {{ display: grid; grid-template-columns: 245px minmax(0, 1fr); gap: 34px; max-width: 1500px; margin: auto; padding: 38px 4vw 80px; }}
nav {{
  position: sticky;
  top: 20px;
  align-self: start;
  padding: 20px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: var(--white);
}}
nav a {{ display: block; padding: 7px 0; color: var(--blue); text-decoration: none; font-size: 14px; }}
main {{ min-width: 0; }}
section {{
  margin-bottom: 30px;
  padding: 30px;
  border: 1px solid var(--line);
  border-radius: 18px;
  background: var(--white);
  box-shadow: 0 8px 24px rgba(36,49,58,.045);
}}
h2 {{ margin: 0 0 14px; font-size: 28px; line-height: 1.2; }}
h3 {{ margin: 4px 0 12px; font-size: 21px; line-height: 1.25; }}
h4 {{ margin: 18px 0 10px; font-size: 17px; }}
p {{ margin: 9px 0; }}
.note {{
  padding: 15px 17px;
  border-left: 4px solid var(--sand);
  background: #fbf4f0;
  color: #654c41;
}}
.ok {{ border-left-color: var(--teal); background: #eef6f1; color: #345b48; }}
.main-figure {{ margin: 24px 0 34px; }}
.figure-head {{ display: flex; align-items: end; justify-content: space-between; gap: 20px; }}
.eyebrow {{ color: var(--sand); text-transform: uppercase; letter-spacing: .12em; font-size: 11px; font-weight: 700; }}
.main-figure img {{
  display: block;
  width: 100%;
  margin-top: 15px;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: white;
}}
figcaption {{ margin-top: 10px; color: var(--muted); font-size: 13px; }}
.button {{ padding: 8px 12px; border-radius: 9px; color: white; background: var(--blue); text-decoration: none; font-size: 12px; font-weight: 700; white-space: nowrap; }}
.table-wrap {{ width: 100%; overflow-x: auto; margin: 14px 0; }}
.data-table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
.data-table th {{ position: sticky; top: 0; padding: 9px; text-align: left; color: white; background: var(--blue); }}
.data-table td {{ padding: 8px 9px; border-bottom: 1px solid var(--line); vertical-align: top; }}
.data-table tr:nth-child(even) td {{ background: #f8faf9; }}
.callout-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 14px; }}
.callout {{ padding: 16px; border-radius: 12px; background: #f0f5f3; }}
.callout strong {{ color: var(--blue); }}
.gallery-group {{ margin-top: 22px; }}
.gallery-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }}
.gallery-card {{ padding: 12px; border: 1px solid var(--line); border-radius: 12px; color: var(--ink); background: #fbfcfb; text-decoration: none; }}
.gallery-card img {{ display: block; width: 100%; border-radius: 7px; background: white; }}
.gallery-card span {{ display: block; margin-top: 9px; font-size: 12px; color: var(--muted); }}
details {{ margin: 12px 0; }}
summary {{ cursor: pointer; font-weight: 700; color: var(--blue); }}
footer {{ color: var(--muted); font-size: 12px; text-align: center; padding: 20px; }}
@media (max-width: 1050px) {{
  .layout {{ grid-template-columns: 1fr; }}
  nav {{ position: static; }}
  .metric-grid {{ grid-template-columns: repeat(3, 1fr); }}
}}
@media (max-width: 680px) {{
  .metric-grid, .gallery-grid, .callout-grid {{ grid-template-columns: 1fr 1fr; }}
  .hero {{ padding: 42px 5vw; }}
  section {{ padding: 20px; }}
}}
</style>
</head>
<body>
<header class="hero">
  <div class="kicker">Scientific Data · data descriptor working report</div>
  <h1>Longitudinal soybean transcriptomes acquired by active microneedle sampling</h1>
  <p>A dataset-centred account of the sampling workflow, technical characteristics, longitudinal expression structure and reusable biological analyses across Day0–Day3.</p>
  <div class="metric-grid">{feature_html}</div>
</header>
<div class="layout">
  <nav>
    <strong>Contents</strong>
    <a href="#overview">Dataset overview</a>
    <a href="#workflow">Figure 1 · Workflow</a>
    <a href="#quality">Figure 2 · Technical comparability</a>
    <a href="#structure">Figure 3 · Expression structure</a>
    <a href="#concordance">Figure 4 · MN versus Bulk</a>
    <a href="#pseudotime">Figure 5 · Pseudotime comparison</a>
    <a href="#biology">Figure 6 · Dataset-supported analyses</a>
    <a href="#gallery">Supplementary analyses</a>
    <a href="#gallery">Supplementary figures</a>
    <a href="#completeness">Data completeness</a>
    <a href="manuscript/legacy_selected5_report.html" target="_blank">Legacy selected-5 report</a>
  </nav>
  <main>
    <section id="overview">
      <h2>Dataset overview</h2>
      <p>The dataset combines 20 longitudinal microneedle samples, representing five seeds followed from Day0 to Day3, with 14 independent destructive bulk samples. The primary reusable products are FASTQ-derived gene counts, normalized expression matrices, sample metadata linking sampling method, Seed ID and Day, quality-control summaries and scripted downstream analyses.</p>
      <div class="callout-grid">
        <div class="callout"><strong>Core dataset</strong><br>34 sequencing runs; four calendar days; one microneedle sample per seed per day; destructive bulk samples kept as a separate comparison set.</div>
        <div class="callout"><strong>Analysis boundary</strong><br>PCA and related ordinations are descriptive. Expression concordance describes the relationship between sampling routes but does not establish formal equivalence.</div>
        <div class="callout"><strong>Pseudotime</strong><br>Pseudotime is treated as a data-derived ordering and compared with recorded Day. It is not used as evidence of sample identity or as a substitute for experimental time.</div>
        <div class="callout"><strong>Traceability</strong><br>Original sequencing-run identifiers remain in the machine-readable crosswalk but are not used as colours or biological group labels in the main figures.</div>
      </div>
    </section>

    <section id="workflow">
      <h2>Data acquisition</h2>
      {extract_figure("figure-1")}
    </section>

    <section id="quality">
      <h2>Technical comparability</h2>
      {extract_figure("figure-2")}
      <div class="note">The figure intentionally displays only metrics with P ≥ 0.05. Input reads (P = 0.023) and clean reads (P = 0.025) were not displayed; the full test table is retained below so the selection is transparent.</div>
      <div class="table-wrap">{table(qc_display)}</div>
    </section>

    <section id="structure">
      <h2>Expression structure</h2>
      {extract_figure("figure-3")}
    </section>

    <section id="concordance">
      <h2>Microneedle versus destructive bulk profiles</h2>
      {extract_figure("figure-4")}
    </section>

    <section id="pseudotime">
      <h2>Pseudotime as a sensitivity analysis</h2>
      {extract_figure("figure-5")}
      <div class="table-wrap">{table(pt_display)}</div>
    </section>

    <section id="biology">
      <h2>Biological signals retained from the earlier selected-5 analysis</h2>
      {extract_figure("figure-6")}
      <p>These analyses are reported as exploratory characterizations of the dataset rather than as independent validation of the sampling platform. They increase the reuse value of the dataset by showing that the matrices support conventional differential-expression, temporal-clustering, functional-enrichment and trajectory analyses.</p>
      <h3>Differential expression across the time course</h3>
      <p>The largest significant DE sets occur for Day2 versus Day1 and Day3 versus Day2, whereas Day3 versus Day1 has the smallest set. This pattern is consistent with a non-linear transcriptional transition rather than a single monotonic trend.</p>
      <div class="table-wrap">{table(de_summary_display)}</div>

      <h3>Candidate developmental genes</h3>
      <p>Four developmental or reserve-mobilization candidates passed one-way ANOVA after the previously requested gene list was applied. The Tukey results indicate a transient Day2 peak for ICL1, NAC2, AS1 and TOP6B relative to Day0 and Day1, followed by lower values at Day3.</p>
      <div class="table-wrap">{table(candidate_display)}</div>

      <h3>Temporal expression modules</h3>
      <p>Mfuzz identified six expression modules with distinct Day1, Day2 and Day3 peaks. The modules include ribosome- and translation-related programmes, chloroplast and electron-transport terms, and reserve-mobilization processes. These are descriptive temporal associations and do not establish causal regulation.</p>
      <div class="table-wrap">{table(mfuzz_display)}</div>

      <h3>Pseudotime-associated genes retained from the earlier analysis</h3>
      <p>The named genes provide examples of coordinated change along the inferred ordering, including WRKY22, ICL1 and AS1. These gene trends describe the data-derived trajectory and are not interpreted as proof of chronological age.</p>
      <div class="table-wrap">{table(pseudotime_gene_display)}</div>
    </section>

    <section id="gallery">
      <h2>Supplementary analysis gallery</h2>
      <p>The original figures are retained in the bundle. Their visual language is preserved, while the main report clarifies which analyses are descriptive, which are supervised and where the pseudotime interpretation stops.</p>
      {gallery_html}
    </section>

    <section id="completeness">
      <h2>Data completeness and reuse</h2>
      <div class="callout-grid">
        <div class="callout"><strong>Sample layer</strong><br>Corrected sample identifiers, sampling method, Seed ID, Day, GTF annotation and run-level provenance.</div>
        <div class="callout"><strong>Expression layer</strong><br>Raw counts, DESeq2 normalized counts, VST expression and the filtered analysis matrix.</div>
        <div class="callout"><strong>Quality layer</strong><br>Read counts, Q20/Q30, retention, mapping fractions, assigned counts, gene-detection metrics and complete ANOVA tables.</div>
        <div class="callout"><strong>Analysis layer</strong><br>Six pairwise DE contrasts, candidate-gene ANOVA/Tukey, Mfuzz modules, enrichment tables and three pseudotime methods with source data.</div>
      </div>
      <div class="note ok">The complete analysis record is retained in a versioned bundle. The main figures preferentially present non-significant technical metrics where requested, while significant results remain available in Source Data and the supplementary gallery.</div>
    </section>
  </main>
</div>
<footer>Generated 2026-10-28 · R figure backend · source data and QA retained alongside the HTML</footer>
</body>
</html>
"""

output_path = OUT / "doudou_scientific_data_report_v3.html"
output_path.write_text(html, encoding="utf-8")
print(output_path)
