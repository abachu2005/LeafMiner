#!/usr/bin/env python3
"""
merge_tissue_results.py

Combine per-tissue LeafCutter2 pipeline outputs into a single top-level
summary.json and merged artifacts, matching the output contract expected
by the webapp UI.

Memory-efficient: computes per-tissue metrics directly from individual
tissue files instead of loading everything into one giant DataFrame.

Usage:
    python merge_tissue_results.py \
        --rundir /path/to/job_root \
        --tissues "Brain - Cortex,Liver" \
        --prefix web_run
"""

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


VALID_CLASSIFICATIONS = {"UP", "PR", "NE", "IN"}

_STAGE_PATTERN = re.compile(
    r"\d+(hpf|dpf|cell|somite)"
    r"|dome|shield|bud|sphere|oblong|gastrula|pharyngula|segmentation"
    r"|blastula|cleavage|epiboly|hatching|larval",
    re.IGNORECASE,
)


def _is_developmental_stage_labels(labels) -> bool:
    return any(_STAGE_PATTERN.fullmatch(str(l).strip()) for l in labels)


def _group_noun(labels) -> str:
    return "developmental stage" if _is_developmental_stage_labels(labels) else "tissue"


def _safe_tissue_dir(tissue: str) -> str:
    out = tissue.replace(" - ", "_").replace(" ", "_")
    for ch in "()/'\"":
        out = out.replace(ch, "")
    return out


def _find_tissue_file(tdir: Path, pattern: str) -> Optional[Path]:
    candidates = list(tdir.glob(pattern))
    return candidates[0] if candidates else None


def collect_junction_filelists(tissue_dirs: Dict[str, Path], out_path: Path) -> Path:
    """Merge per-tissue junction_files.txt into one combined list."""
    all_paths: List[str] = []
    for tissue, tdir in sorted(tissue_dirs.items()):
        jf = tdir / "junctions_bed" / "junction_files.txt"
        if not jf.exists():
            bed_dir_files = list((tdir / "junctions_bed").glob("*/*.juncs.bed"))
            all_paths.extend(str(p) for p in bed_dir_files)
            continue
        all_paths.extend(l.strip() for l in jf.read_text().splitlines() if l.strip())
    out_path.write_text("\n".join(all_paths) + "\n")
    return out_path


def _compute_tissue_unproductive(
    tissue: str, jc_path: Path
) -> Optional[Dict[str, Any]]:
    """Compute unproductive-read stats for a single tissue's junction_counts.gz.

    Returns a dict with per-sample stats, or None if no classifications found.
    """
    with gzip.open(str(jc_path), "rt") as f:
        header_tokens = f.readline().strip().split()
        sample_cols = header_tokens[1:]
        n = len(sample_cols)
        if n == 0:
            return None

        total_reads = np.zeros(n, dtype=np.float64)
        up_reads = np.zeros(n, dtype=np.float64)
        has_classification = False

        for line in f:
            fields = line.strip().split()
            if len(fields) < 2:
                continue
            junc_id = fields[0]
            reads = np.array([float(x) for x in fields[1 : n + 1]])
            label = junc_id.rsplit(":", 1)[-1]
            if label in VALID_CLASSIFICATIONS:
                has_classification = True
            total_reads += reads
            if label == "UP":
                up_reads += reads

    if not has_classification:
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(total_reads > 0, up_reads / total_reads * 100.0, np.nan)

    valid = ~np.isnan(pct)
    if not valid.any():
        return None

    vals = pct[valid]
    return {
        "tissue": tissue,
        "n_samples": int(valid.sum()),
        "mean_pct": round(float(np.mean(vals)), 4),
        "median_pct": round(float(np.median(vals)), 4),
        "std_pct": round(float(np.std(vals, ddof=1)), 4) if len(vals) > 1 else 0.0,
        "total_reads": int(total_reads[valid].sum()),
        "unproductive_reads": int(up_reads[valid].sum()),
    }


def compute_all_tissue_metrics(
    tissue_dirs: Dict[str, Path], outdir: Path
) -> Optional[Dict[str, Any]]:
    """Compute unproductive-read metrics across all tissues without merging files."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        has_mpl = True
    except ImportError:
        has_mpl = False

    rows: List[Dict[str, Any]] = []
    for tissue, tdir in sorted(tissue_dirs.items()):
        jc = _find_tissue_file(tdir, "lc2/*_lc2.junction_counts.gz")
        if jc is None:
            print(f"  [WARN] No junction_counts for {tissue}, skipping metrics")
            continue
        print(f"  Computing metrics for {tissue}...")
        result = _compute_tissue_unproductive(tissue, jc)
        if result:
            rows.append(result)

    if not rows:
        return None

    import pandas as pd
    result_df = pd.DataFrame(rows).sort_values("mean_pct").reset_index(drop=True)

    all_pcts = []
    total_up = 0
    total_all = 0
    for r in rows:
        total_up += r["unproductive_reads"]
        total_all += r["total_reads"]

    global_weighted_mean = (total_up / total_all * 100.0) if total_all > 0 else 0.0
    all_means = [r["mean_pct"] for r in rows]
    global_median = float(np.median(all_means))
    n_tissues = len(rows)

    tsv_path = outdir / "unproductive_by_tissue.tsv"
    result_df.to_csv(str(tsv_path), sep="\t", index=False)

    json_path = outdir / "unproductive_by_tissue.json"
    payload = {
        "global_weighted_mean_pct": round(global_weighted_mean, 4),
        "global_median_pct": round(global_median, 4),
        "n_tissues": n_tissues,
        "per_tissue": rows,
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    png_path = outdir / "unproductive_by_tissue.png"
    if has_mpl and n_tissues > 0:
        fig, ax = plt.subplots(figsize=(max(8, n_tissues * 0.35), 5))
        tissues_sorted = result_df["tissue"].tolist()
        means = result_df["mean_pct"].tolist()
        stds = result_df["std_pct"].tolist()
        group = _group_noun(tissues_sorted)
        ax.barh(range(n_tissues), means, xerr=stds, color="#4C72B0", alpha=0.85)
        ax.set_yticks(range(n_tissues))
        ax.set_yticklabels(tissues_sorted, fontsize=7)
        ax.set_xlabel("Unproductive splice reads (%)")
        ax.set_title(f"Unproductive splicing by {group}")
        ax.axvline(global_weighted_mean, color="red", ls="--", lw=0.8, label=f"weighted mean {global_weighted_mean:.1f}%")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(str(png_path), dpi=150)
        plt.close(fig)

    return {
        "tsv": str(tsv_path),
        "json": str(json_path),
        "png": str(png_path) if png_path.exists() else None,
        "global_weighted_mean_pct": round(global_weighted_mean, 4),
        "global_median_pct": round(global_median, 4),
        "n_tissues": n_tissues,
    }


def _collect_per_tissue_lc2_paths(
    tissue_dirs: Dict[str, Path],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return dicts mapping tissue -> path for junction_counts and cluster_ratios."""
    jc_paths: Dict[str, str] = {}
    cr_paths: Dict[str, str] = {}
    for tissue, tdir in sorted(tissue_dirs.items()):
        jc = _find_tissue_file(tdir, "lc2/*_lc2.junction_counts.gz")
        if jc:
            jc_paths[tissue] = str(jc)
        cr = _find_tissue_file(tdir, "lc2/*_lc2.cluster_ratios.gz")
        if cr:
            cr_paths[tissue] = str(cr)
    return jc_paths, cr_paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-tissue LC2 results")
    parser.add_argument("--rundir", required=True, help="Top-level job run directory")
    parser.add_argument("--tissues", required=True, help="Comma-separated tissue names (same order as array)")
    parser.add_argument("--prefix", default="web_run")
    args = parser.parse_args()

    rundir = Path(args.rundir)
    tissues = [t.strip() for t in args.tissues.split(",") if t.strip()]
    prefix = args.prefix

    print(f"[MERGE] Merging results for {len(tissues)} tissues in {rundir}")

    tissue_dirs: Dict[str, Path] = {}
    failed_tissues: List[str] = []
    for tissue in tissues:
        tdir = rundir / f"tissue_{_safe_tissue_dir(tissue)}"
        summary = tdir / "out" / "summary.json"
        if tdir.exists() and summary.exists():
            tissue_dirs[tissue] = tdir
        else:
            failed_tissues.append(tissue)
            print(f"  [WARN] Tissue {tissue!r} missing or incomplete at {tdir}")

    if not tissue_dirs:
        print("[ERROR] No tissue results found to merge", file=sys.stderr)
        sys.exit(1)

    if failed_tissues:
        print(f"[WARN] {len(failed_tissues)} tissue(s) failed: {failed_tissues}")

    outdir = rundir / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    lc2_outdir = rundir / "lc2"
    lc2_outdir.mkdir(parents=True, exist_ok=True)

    print("[MERGE] Collecting junction file lists...")
    merged_filelist = rundir / "clusters" / f"{prefix}_junction_files.txt"
    merged_filelist.parent.mkdir(parents=True, exist_ok=True)
    collect_junction_filelists(tissue_dirs, merged_filelist)

    print("[MERGE] Collecting per-tissue LC2 output paths...")
    jc_paths, cr_paths = _collect_per_tissue_lc2_paths(tissue_dirs)

    print(f"[MERGE] Computing tissue metrics from {len(jc_paths)} tissues (streaming, one at a time)...")
    tissue_metrics = compute_all_tissue_metrics(tissue_dirs, outdir)

    per_tissue_summaries = {}
    total_junction_inputs = 0
    for tissue, tdir in sorted(tissue_dirs.items()):
        ts = tdir / "out" / "summary.json"
        if ts.exists():
            with open(ts) as f:
                per_tissue_summaries[tissue] = json.load(f)
                total_junction_inputs += per_tissue_summaries[tissue].get("n_junction_inputs", 0)

    summary: Dict[str, Any] = {
        "n_junction_inputs": total_junction_inputs,
        "junction_file_list": str(merged_filelist),
        "parallel_tissues": True,
        "n_tissues_requested": len(tissues),
        "n_tissues_succeeded": len(tissue_dirs),
        "failed_tissues": failed_tissues,
        "lc2_outputs": {
            "per_tissue_junction_counts": jc_paths,
            "per_tissue_cluster_ratios": cr_paths,
            "run_dir": str(lc2_outdir),
        },
        "cluster_files": {},
        "metrics": {},
        "per_tissue_summaries": per_tissue_summaries,
        "notes": "Parallel per-tissue pipeline run. Metrics computed per-tissue to avoid OOM.",
    }
    if tissue_metrics:
        summary["metrics"]["unproductive_by_tissue"] = tissue_metrics

    summary_path = outdir / "summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n[MERGE] Done.")
    print(f"  - Summary: {summary_path}")
    print(f"  - Tissues succeeded: {len(tissue_dirs)}/{len(tissues)}")
    if tissue_metrics:
        print(f"  - Tissue metrics: {tissue_metrics.get('n_tissues', 0)} tissues")
    if failed_tissues:
        print(f"  - Skipped tissues (no results): {failed_tissues}")


if __name__ == "__main__":
    main()
