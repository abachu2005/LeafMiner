#!/usr/bin/env python3
"""Aggregate zebrafish pipeline results from Quest into outputs/for_saba/.

For each completed run, SCP the relevant artifacts, then build a 3-way
comparison table and Excel workbook.

Usage:
  python3 scripts/aggregate_for_saba.py [--dry_run]
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

QUEST_HOST = "login.quest.northwestern.edu"
QUEST_USER = "iis1026"
QUEST_REPO = "/projects/p52853/iis1026/Leaf_Cutter"

OUTPUT_ROOT = Path("outputs/for_saba")

RUNS = {
    "run1_grcz12_stock_baseline": {
        "label": "GRCz12tu stock (baseline)",
        "assembly": "GRCz12tu",
        "gtf": "GRCz12tu.ucsc.gtf (stock RefSeq)",
        "remote_dir": None,
        "stage3_job": "583631a5-5149-4b91-af10-69fafe05224d",
        "stage1_job": "1c0556db-a2af-44e9-86d4-45e9c40bf49b",
        "local_source": "outputs/stage3_v2",
    },
    "run2_grcz12_longorf": {
        "label": "GRCz12tu + LongORF supplement",
        "assembly": "GRCz12tu",
        "gtf": "GRCz12tu.ucsc.longorf.gtf (stock + 24,694 NMD transcripts)",
        "chain_json": "outputs/for_saba/_runs/run2_grcz12_longorf.json",
    },
    "run3_grcz11_stock": {
        "label": "GRCz11 (danRer11) stock Ensembl 115",
        "assembly": "GRCz11 (danRer11)",
        "gtf": "Danio_rerio.GRCz11.115.ucsc.gtf (Ensembl 115; LongORF liftover yield <70%, used stock)",
        "chain_json": "outputs/for_saba/_runs/run3_grcz11_stock.json",
    },
}

ARTIFACTS = [
    "out/summary.json",
    "out/stage3_pe_inclusion.json",
    "out/pe_inclusion/pe_developmental_candidates.tsv",
    "out/pe_inclusion/pe_developmental_candidates.json",
    "out/pe_inclusion/pe_events.tsv",
    "out/pe_inclusion/summary.json",
    "out/pe_inclusion/pe_psi_heatmap.png",
    "out/pe_inclusion/pe_psi_vs_time_top.png",
    "out/pe_inclusion/pe_dpsi_vs_pvalue.png",
    "out/pe_inclusion/pe_spearman_distribution.png",
    "out/pe_inclusion/pe_lc2_concordance.png",
    "out/pe_inclusion/pe_low_coverage.tsv",
    "lc2/*_lc2.cluster_ratios.gz",
    "lc2/*_lc2.junction_counts.gz",
    "lc2/clustering/*_junction_classifications.txt",
]


def scp_file(remote_path: str, local_path: Path, dry_run: bool = False) -> bool:
    local_path.parent.mkdir(parents=True, exist_ok=True)
    src = f"{QUEST_USER}@{QUEST_HOST}:{remote_path}"
    if dry_run:
        print(f"  [DRY] scp {src} -> {local_path}")
        return True
    try:
        subprocess.run(
            ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", src, str(local_path)],
            check=True, capture_output=True, text=True,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def scp_glob(remote_pattern: str, local_dir: Path, dry_run: bool = False) -> int:
    """SCP files matching a remote glob pattern."""
    local_dir.mkdir(parents=True, exist_ok=True)
    try:
        result = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
             f"{QUEST_USER}@{QUEST_HOST}", f"ls {remote_pattern} 2>/dev/null"],
            capture_output=True, text=True, check=False,
        )
        files = [f.strip() for f in result.stdout.strip().splitlines() if f.strip()]
        count = 0
        for f in files:
            name = Path(f).name
            if scp_file(f, local_dir / name, dry_run):
                count += 1
        return count
    except Exception:
        return 0


def fetch_run(run_name: str, run_info: dict, dry_run: bool = False) -> Path:
    """Fetch all artifacts for a single run."""
    local_dir = OUTPUT_ROOT / run_name
    local_dir.mkdir(parents=True, exist_ok=True)

    if run_info.get("local_source"):
        local_src = Path(run_info["local_source"])
        if local_src.exists():
            print(f"  Copying from local source: {local_src}")
            if not dry_run:
                import shutil
                for item in local_src.iterdir():
                    dest = local_dir / item.name
                    if item.is_file():
                        shutil.copy2(item, dest)
                    elif item.is_dir():
                        if dest.exists():
                            shutil.rmtree(dest)
                        shutil.copytree(item, dest)

        if run_info.get("stage3_job"):
            remote_base = f"{QUEST_REPO}/jobs/{run_info['stage3_job']}"
            for art in ["lc2/*_lc2.cluster_ratios.gz", "lc2/*_lc2.junction_counts.gz",
                        "lc2/clustering/*_junction_classifications.txt"]:
                scp_glob(f"{remote_base}/{art}", local_dir / Path(art).parent, dry_run)
        return local_dir

    chain_json = run_info.get("chain_json", "")
    if chain_json and Path(chain_json).exists():
        with open(chain_json) as fh:
            chain = json.load(fh)
        remote_dir = chain["remote_run_dir"]
    else:
        print(f"  WARNING: No chain JSON for {run_name}")
        return local_dir

    print(f"  Remote dir: {remote_dir}")
    for art in ARTIFACTS:
        if "*" in art:
            scp_glob(f"{remote_dir}/{art}", local_dir / Path(art).parent, dry_run)
        else:
            scp_file(f"{remote_dir}/{art}", local_dir / art, dry_run)

    return local_dir


def load_summary(run_dir: Path) -> dict:
    """Load the PE inclusion summary.json for a run."""
    for candidate in [
        run_dir / "out" / "pe_inclusion" / "summary.json",
        run_dir / "pe_inclusion" / "summary.json",
        run_dir / "summary.json",
    ]:
        if candidate.exists():
            with open(candidate) as fh:
                return json.load(fh)
    return {}


def load_candidates(run_dir: Path) -> list:
    """Load PE developmental candidates TSV."""
    for candidate in [
        run_dir / "out" / "pe_inclusion" / "pe_developmental_candidates.tsv",
        run_dir / "pe_inclusion" / "pe_developmental_candidates.tsv",
        run_dir / "pe_developmental_candidates.tsv",
    ]:
        if candidate.exists():
            import csv
            with open(candidate) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                return list(reader)
    return []


def summary_metric(summary: dict, metric: str):
    """Return a report metric from current or older summary schemas."""
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


def build_comparison(dry_run: bool = False):
    """Build 3-way comparison TSV and optional Excel."""
    all_candidates: Dict[str, list] = {}
    all_summaries: Dict[str, dict] = {}

    for run_name in RUNS:
        run_dir = OUTPUT_ROOT / run_name
        all_summaries[run_name] = load_summary(run_dir)
        all_candidates[run_name] = load_candidates(run_dir)

    comparison_path = OUTPUT_ROOT / "comparison_summary.tsv"
    with open(comparison_path, "w") as fh:
        fh.write("metric\t" + "\t".join(RUNS.keys()) + "\n")

        metrics = [
            "n_events_in_psi", "n_events_passing_all_filters",
            "n_pe_events_tested", "n_candidates",
            "n_significant_spearman", "n_significant_dpsi",
        ]
        for metric in metrics:
            vals = []
            for run_name in RUNS:
                s = all_summaries.get(run_name, {})
                v = summary_metric(s, metric)
                vals.append(str(v))
            fh.write(f"{metric}\t" + "\t".join(vals) + "\n")

    gene_sets: Dict[str, set] = {}
    for run_name, candidates in all_candidates.items():
        genes = set()
        for row in candidates:
            gene = row.get("gene_symbol") or row.get("gene_name") or row.get("gene_id", "")
            if gene:
                genes.add(gene)
        gene_sets[run_name] = genes

    novel_path = OUTPUT_ROOT / "novel_in_longorf.tsv"
    run1_genes = gene_sets.get("run1_grcz12_stock_baseline", set())
    run2_genes = gene_sets.get("run2_grcz12_longorf", set())
    novel_genes = run2_genes - run1_genes
    with open(novel_path, "w") as fh:
        fh.write("gene\tnote\n")
        for g in sorted(novel_genes):
            fh.write(f"{g}\tnovel in Run 2 (LongORF supplement)\n")

    print(f"\nComparison summary: {comparison_path}")
    print(f"Novel LongORF genes: {len(novel_genes)} -> {novel_path}")

    try:
        import openpyxl
        wb = openpyxl.Workbook()
        for i, (run_name, candidates) in enumerate(all_candidates.items()):
            if i == 0:
                ws = wb.active
                ws.title = run_name[:31]
            else:
                ws = wb.create_sheet(title=run_name[:31])
            if candidates:
                headers = list(candidates[0].keys())
                ws.append(headers)
                for row in candidates:
                    ws.append([row.get(h, "") for h in headers])
        xlsx_path = OUTPUT_ROOT / "combined_pe_candidates.xlsx"
        wb.save(xlsx_path)
        print(f"Excel workbook: {xlsx_path}")
    except ImportError:
        print("openpyxl not installed — skipping Excel output")

    return all_summaries, gene_sets


def main():
    p = argparse.ArgumentParser(description="Aggregate zebrafish results for Saba")
    p.add_argument("--dry_run", action="store_true")
    args = p.parse_args()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    for run_name, run_info in RUNS.items():
        print(f"\n=== Fetching {run_name}: {run_info['label']} ===")
        fetch_run(run_name, run_info, args.dry_run)

    print("\n=== Building comparison ===")
    summaries, gene_sets = build_comparison(args.dry_run)

    print(f"\n=== Done. Results in {OUTPUT_ROOT} ===")
    for run_name, s in summaries.items():
        n_cand = summary_metric(s, "n_candidates")
        n_tested = summary_metric(s, "n_pe_events_tested")
        print(f"  {run_name}: {n_cand} candidate rows; {n_tested} GLM-tested events")


if __name__ == "__main__":
    main()
