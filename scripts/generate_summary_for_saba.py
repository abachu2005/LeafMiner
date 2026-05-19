#!/usr/bin/env python3
"""Generate SUMMARY_FOR_SABA.md from aggregated run results.

Reads the per-run summary.json and pe_developmental_candidates.tsv files
from outputs/for_saba/<run>/ and produces a formatted Markdown report.

Usage:
  python3 scripts/generate_summary_for_saba.py
"""

import json
import csv
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional

OUTPUT_ROOT = Path("outputs/for_saba")

RUNS = {
    "run1_grcz12_stock_baseline": {
        "label": "Run 1: GRCz12tu stock (baseline)",
        "assembly": "GRCz12tu",
        "gtf_desc": "RefSeq stock GTF (GRCz12tu.ucsc.gtf)",
        "notes": "Reuses prior run artifacts (job 583631a5).",
    },
    "run2_grcz12_longorf": {
        "label": "Run 2: GRCz12tu + LongORF supplement",
        "assembly": "GRCz12tu",
        "gtf_desc": "Stock GTF + 24,694 LongORF NMD transcripts (GRCz12tu.ucsc.longorf.gtf)",
        "notes": "Reuses Stage 1 SJ files from Run 1. Only Stage 2 (LC2) and Stage 3 (PE inclusion) re-run with supplemented GTF.",
    },
    "run3_grcz11_stock": {
        "label": "Run 3: GRCz11 (danRer11) Ensembl 115",
        "assembly": "GRCz11 (danRer11)",
        "gtf_desc": "Ensembl 115 GTF (Danio_rerio.GRCz11.115.ucsc.gtf); LongORF liftover yield was 16.5% (<70% threshold), so stock GTF used.",
        "notes": "Fresh alignment of all 276 samples to GRCz11. Slurm array parallelization (20 concurrent tasks).",
    },
}


def load_summary(run_dir: Path) -> dict:
    for p in [
        run_dir / "out" / "pe_inclusion" / "summary.json",
        run_dir / "pe_inclusion" / "summary.json",
        run_dir / "summary.json",
    ]:
        if p.exists():
            with open(p) as fh:
                return json.load(fh)
    return {}


def load_candidates(run_dir: Path) -> list:
    for p in [
        run_dir / "out" / "pe_inclusion" / "pe_developmental_candidates.tsv",
        run_dir / "pe_inclusion" / "pe_developmental_candidates.tsv",
        run_dir / "pe_developmental_candidates.tsv",
    ]:
        if p.exists():
            with open(p) as fh:
                return list(csv.DictReader(fh, delimiter="\t"))
    return []


def safe_int(v, default="N/A"):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def summary_metric(summary: dict, metric: str):
    stage_c = summary.get("stage_c", {}) if isinstance(summary, dict) else {}
    aliases = {
        "n_pe_events_tested": ("n_events_tested_glm", "n_events_tested"),
        "n_events_passing_all_filters": ("n_events_passing_all_filters",),
        "n_significant_spearman": ("n_significant_spearman_q05",),
        "n_significant_dpsi": ("n_significant_combined",),
    }
    if metric in summary:
        return summary.get(metric)
    for key in aliases.get(metric, ()):
        if key in stage_c:
            return stage_c.get(key)
    return "N/A"


def generate():
    summaries: Dict[str, dict] = {}
    candidates: Dict[str, list] = {}
    gene_sets: Dict[str, set] = {}

    for run_name in RUNS:
        run_dir = OUTPUT_ROOT / run_name
        summaries[run_name] = load_summary(run_dir)
        candidates[run_name] = load_candidates(run_dir)
        genes = set()
        for row in candidates[run_name]:
            g = row.get("gene_symbol") or row.get("gene_name") or row.get("gene_id", "")
            if g:
                genes.add(g)
        gene_sets[run_name] = genes

    r1_genes = gene_sets.get("run1_grcz12_stock_baseline", set())
    r2_genes = gene_sets.get("run2_grcz12_longorf", set())
    r3_genes = gene_sets.get("run3_grcz11_stock", set())
    novel_r2 = r2_genes - r1_genes
    shared_r2_r3 = r2_genes & r3_genes

    md = []
    md.append("# Zebrafish Poison-Exon Developmental Analysis — Triple-Run Comparison\n")
    md.append(f"*Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n")

    md.append("## 1. Headline\n")
    md.append("Three parallel pipeline runs were executed on **276 zebrafish developmental RNA-seq")
    md.append("libraries** (6 timepoints: 24hpf → 5dpf) to assess the impact of:\n")
    md.append("1. **GTF annotation completeness** — adding 24,694 LongORF PTC+ NMD transcripts to the stock GTF")
    md.append("2. **Assembly version** — GRCz12tu (2025) vs GRCz11/danRer11 (2017, Ensembl 115)\n")

    md.append("### Key metrics\n")
    md.append("| Metric | Run 1 (GRCz12tu stock) | Run 2 (GRCz12tu + LongORF) | Run 3 (GRCz11 stock) |")
    md.append("|---|---|---|---|")
    metrics = [
        ("Events in PSI matrix", "n_events_in_psi"),
        ("Events passing all coverage filters", "n_events_passing_all_filters"),
        ("GLM-tested PE events", "n_pe_events_tested"),
        ("Candidate rows reported", "n_candidates"),
    ]
    for label, key in metrics:
        vals = []
        for run_name in RUNS:
            v = summary_metric(summaries.get(run_name, {}), key)
            vals.append(str(safe_int(v)))
        md.append(f"| {label} | {' | '.join(vals)} |")
    md.append(f"| Unique candidate genes | {len(r1_genes)} | {len(r2_genes)} | {len(r3_genes)} |")
    md.append(f"| Novel genes from LongORF | — | {len(novel_r2)} | — |")
    md.append(f"| Genes shared (Run 2 ∩ Run 3) | — | — | {len(shared_r2_r3)} |\n")

    md.append("## 2. Run Descriptions\n")
    for run_name, info in RUNS.items():
        md.append(f"### {info['label']}\n")
        md.append(f"- **Assembly:** {info['assembly']}")
        md.append(f"- **GTF:** {info['gtf_desc']}")
        md.append(f"- **Notes:** {info['notes']}\n")

    md.append("## 3. LongORF Supplement Details\n")
    md.append("Talha's `LongORF_PTC+_fromClustered.tsv` contains **24,757 PTC+ transcripts** across")
    md.append("**3,967 genes**. After filtering (63 rows with no genomic coordinates), **24,694 transcripts**")
    md.append("were converted to a proper GTF supplement with `transcript_biotype \"nonsense_mediated_decay\"`")
    md.append("so LeafCutter2 classifies their junctions as **UP** (unproductive).\n")
    md.append("**Liftover to GRCz11:** Attempted using UCSC chain `danRer11ToGCF_049306965.1.over.chain`")
    md.append("with a properly swapped chain + CrossMap. Yield was only **16.5%** (4,080 / 24,694 transcripts),")
    md.append("well below the 70% threshold. The massive coordinate divergence between GRCz12tu (2025 Verkko")
    md.append("assembly) and GRCz11 (2017) makes direct liftover impractical. Run 3 therefore used the")
    md.append("**stock Ensembl 115 GRCz11 GTF** without the LongORF supplement.\n")

    run3_stage_c = summaries.get("run3_grcz11_stock", {}).get("stage_c", {})
    if safe_int(run3_stage_c.get("n_events_tested_glm"), 0) == 0:
        md.append("### Run 3 PE-coordinate limitation\n")
        md.append("Run 3 LC2 completed cleanly on GRCz11 with real gene names, but the POISEN PE list")
        md.append("is still defined in GRCz12tu genomic coordinates. When those PE intervals were")
        md.append("queried against GRCz11 STAR junctions, no PE events passed coverage filters and")
        md.append("no GRCz11 PE figures were generated. Treat Run 3 as a validated GRCz11 LC2/assembly")
        md.append("control, not as a fully comparable GRCz11 PE-inclusion test, unless the PE list is")
        md.append("lifted or regenerated on GRCz11 coordinates.\n")

    if novel_r2:
        md.append("### Novel PE genes uncovered by LongORF supplement\n")
        md.append("These genes became testable in Run 2 but were not testable in Run 1 (baseline):\n")
        for g in sorted(novel_r2)[:25]:
            md.append(f"- {g}")
        if len(novel_r2) > 25:
            md.append(f"- *(... {len(novel_r2) - 25} more in `novel_in_longorf.tsv`)*")
        md.append("")

    md.append("## 4. Pipeline Parallelization\n")
    md.append("The pipeline was refactored to use **Slurm job arrays** for per-sample parallelization:")
    md.append("- QC (fastp trimming): `--array=0-275%20` (20 concurrent downloads)")
    md.append("- Stage 1 (STAR alignment): `--array=0-275%20` (20 concurrent alignments)")
    md.append("- Dependency chaining via `--dependency=afterany:JID` (tolerates individual task failures)")
    md.append("- Email notifications on completion/failure (no busy-polling)\n")
    md.append("This reduced the estimated Run 3 wall-clock time from **~6 days (serial)** to")
    md.append("**~4-8 hours (parallel)**.\n")

    md.append("## 5. Figures\n")
    md.append("### Per-run visualizations\n")
    for run_name, info in RUNS.items():
        md.append(f"**{info['label']}**\n")
        pe_dir = f"{run_name}/pe_inclusion" if (OUTPUT_ROOT / run_name / "pe_inclusion").exists() else run_name
        for fig, desc in [
            ("pe_psi_heatmap.png", "PSI heatmap"),
            ("pe_psi_vs_time_top.png", "PSI vs developmental time (top candidates)"),
            ("pe_dpsi_vs_pvalue.png", "ΔPSI volcano"),
            ("pe_spearman_distribution.png", "Spearman ρ distribution"),
            ("pe_lc2_concordance.png", "LC2 concordance"),
        ]:
            fig_path = OUTPUT_ROOT / run_name / "pe_inclusion" / fig
            if not fig_path.exists():
                fig_path = OUTPUT_ROOT / run_name / fig
            if fig_path.exists():
                rel = fig_path.relative_to(OUTPUT_ROOT)
                md.append(f"![{desc}]({rel})\n")
        md.append("")

    md.append("## 6. Methods\n")
    md.append("Same pipeline as previously documented in `outputs/stage3_v2/RESULTS_FOR_SABA.md`.")
    md.append("Key additions:\n")
    md.append("- **LongORF GTF supplement:** `scripts/longorf_to_gtf.py` converts the LongORF TSV")
    md.append("  exon chains to GTF format with NMD biotype annotation.")
    md.append("- **GTF merge:** `scripts/build_supplemented_gtf.sh` concatenates stock + supplement,")
    md.append("  sorts, and validates.")
    md.append("- **Parallelization:** `scripts/run_zebrafish_full_chain.py` orchestrates the full Slurm")
    md.append("  dependency chain with per-sample array tasks.")
    md.append("- **Liftover attempt:** `scripts/setup_grcz11_longorf_chain.sh` properly swaps UCSC chains")
    md.append("  and uses CrossMap for coordinate conversion.\n")

    md.append("## 7. Folder Map\n")
    md.append("```")
    md.append("outputs/for_saba/")
    md.append("├── run1_grcz12_stock_baseline/   # GRCz12tu stock baseline")
    md.append("│   ├── pe_inclusion/              # PE analysis results")
    md.append("│   └── lc2/                       # LeafCutter2 outputs")
    md.append("├── run2_grcz12_longorf/           # GRCz12tu + LongORF supplement")
    md.append("│   ├── pe_inclusion/")
    md.append("│   └── lc2/")
    md.append("├── run3_grcz11_stock/             # GRCz11 (danRer11)")
    md.append("│   ├── pe_inclusion/")
    md.append("│   └── lc2/")
    md.append("├── comparison_summary.tsv          # 3-way metric comparison")
    md.append("├── novel_in_longorf.tsv            # Genes novel in Run 2")
    md.append("├── combined_pe_candidates.xlsx     # All candidates (3 sheets)")
    md.append("└── SUMMARY_FOR_SABA.md             # This document")
    md.append("```\n")

    text = "\n".join(md)
    out_path = OUTPUT_ROOT / "SUMMARY_FOR_SABA.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(text)
    print(f"Wrote {out_path} ({len(text)} chars)")


if __name__ == "__main__":
    generate()
