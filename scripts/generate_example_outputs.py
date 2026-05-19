#!/usr/bin/env python3
"""
Generate realistic example pipeline outputs for inspection.
Simulates a Quest run with 2 tissues (Brain_Cortex, Liver) x 3 samples each.
"""

import gzip
import json
import os
import random
import sys
from pathlib import Path

import numpy as np

random.seed(42)
np.random.seed(42)

BASE = Path(__file__).resolve().parent.parent / "example_output"

TISSUES = ["Brain_Cortex", "Liver"]
SAMPLES_PER_TISSUE = 3

SAMPLE_IDS = {
    "Brain_Cortex": ["GTEX-1117F", "GTEX-111CU", "GTEX-111FC"],
    "Liver": ["GTEX-1122O", "GTEX-117XS", "GTEX-117YW"],
}

CHROMS = ["chr1", "chr2", "chr3", "chr10", "chr17"]
CLASSIFICATIONS = ["PR", "UP", "NE", "IN"]
CLASS_WEIGHTS = [0.55, 0.08, 0.30, 0.07]

PREFIX = "gtex_demo"
LC2_PREFIX = f"{PREFIX}_lc2"

N_JUNCTIONS_PER_BED = 800
N_CLUSTERS = 120
JUNCTIONS_PER_CLUSTER_RANGE = (2, 8)


def ensure(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p


def rand_intron(chrom):
    start = random.randint(10000, 50000000)
    length = random.choice([500, 1200, 5000, 15000, 45000, 120000])
    end = start + length
    strand = random.choice(["+", "-"])
    return chrom, start, end, strand


def generate_junction_beds(jdir: Path):
    all_beds = []
    for tissue in TISSUES:
        tdir = ensure(jdir / tissue)
        for sid in SAMPLE_IDS[tissue]:
            bed_path = tdir / f"{sid}.juncs.bed"
            with open(bed_path, "w") as f:
                for _ in range(N_JUNCTIONS_PER_BED):
                    chrom = random.choice(CHROMS)
                    ch, s, e, st = rand_intron(chrom)
                    score = max(1, int(np.random.lognormal(3.5, 1.5)))
                    f.write(f"{ch}\t{s}\t{e}\tJUNC\t{score}\t{st}\n")
            all_beds.append(bed_path)
    return all_beds


def generate_junction_filelist(bed_paths, out_path):
    with open(out_path, "w") as f:
        for p in bed_paths:
            f.write(str(p) + "\n")
    return out_path


def generate_perind_counts(sample_names, out_path, n_introns=400):
    clusters = []
    for c in range(1, N_CLUSTERS + 1):
        chrom = random.choice(CHROMS)
        strand = random.choice(["+", "-"])
        n_junc = random.randint(*JUNCTIONS_PER_CLUSTER_RANGE)
        base = random.randint(10000, 50000000)
        for j in range(n_junc):
            s = base + j * random.randint(200, 2000)
            e = s + random.randint(500, 50000)
            clusters.append((chrom, s, e, f"clu_{c}_{strand}"))

    with gzip.open(str(out_path), "wt") as f:
        f.write(" ".join(sample_names) + "\n")
        for chrom, s, e, clu in clusters:
            intron_id = f"{chrom}:{s}:{e}:{clu}"
            fracs = []
            for _ in sample_names:
                numer = max(0, int(np.random.lognormal(2.5, 1.2)))
                denom = numer + max(1, int(np.random.lognormal(3, 1.0)))
                fracs.append(f"{numer}/{denom}")
            f.write(intron_id + " " + " ".join(fracs) + "\n")


def generate_cluster_ratios(sample_names, out_path):
    with gzip.open(str(out_path), "wt") as f:
        f.write(" ".join(sample_names) + "\n")
        for c in range(1, N_CLUSTERS + 1):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            n_junc = random.randint(*JUNCTIONS_PER_CLUSTER_RANGE)
            base = random.randint(10000, 50000000)
            for j in range(n_junc):
                s = base + j * random.randint(200, 2000)
                e = s + random.randint(500, 50000)
                cls = random.choices(CLASSIFICATIONS, weights=CLASS_WEIGHTS, k=1)[0]
                intron_id = f"{chrom}:{s}:{e}:clu_{c}_{strand}:{cls}"
                fracs = []
                for _ in sample_names:
                    numer = max(0, int(np.random.lognormal(2.5, 1.2)))
                    denom = numer + max(1, int(np.random.lognormal(3, 1.0)))
                    fracs.append(f"{numer}/{denom}")
                f.write(intron_id + " " + " ".join(fracs) + "\n")


def generate_junction_counts(sample_names, out_path):
    """junction_counts.gz — raw integer read counts with classification labels."""
    with gzip.open(str(out_path), "wt") as f:
        f.write(" ".join(sample_names) + "\n")
        for c in range(1, N_CLUSTERS + 1):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            n_junc = random.randint(*JUNCTIONS_PER_CLUSTER_RANGE)
            base = random.randint(10000, 50000000)
            for j in range(n_junc):
                s = base + j * random.randint(200, 2000)
                e = s + random.randint(500, 50000)
                cls = random.choices(CLASSIFICATIONS, weights=CLASS_WEIGHTS, k=1)[0]
                intron_id = f"{chrom}:{s}:{e}:clu_{c}_{strand}:{cls}"
                reads = [max(0, int(np.random.lognormal(2.5, 1.2))) for _ in sample_names]
                f.write(intron_id + " " + " ".join(str(r) for r in reads) + "\n")


def generate_clustering_intermediates(clustering_dir: Path, prefix: str):
    """Generate the LC2 clustering intermediate files (clusters, exon_stats,
    junction_classifications, lowusage_introns) that live alongside the
    distance files."""
    gene_names = [
        "BRCA1", "TP53", "GAPDH", "ACTB", "MYC", "RB1", "PTEN", "EGFR",
        "KRAS", "VEGFA", "FOS", "JUN", "CDK2", "BCL2", "AKT1",
    ]

    # clusters_clusters — cluster membership
    with open(clustering_dir / f"{prefix}_clusters", "w") as f:
        for c in range(1, N_CLUSTERS + 1):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            n_junc = random.randint(*JUNCTIONS_PER_CLUSTER_RANGE)
            base = random.randint(10000, 50000000)
            introns = []
            for j in range(n_junc):
                s = base + j * random.randint(200, 2000)
                e = s + random.randint(500, 50000)
                introns.append(f"{chrom}:{s}:{e}")
            f.write(f"clu_{c}_{strand}\t" + "\t".join(introns) + "\n")

    # clusters_exon_stats.txt
    with open(clustering_dir / f"{prefix}_exon_stats.txt", "w") as f:
        f.write("junction_id\texon_length\tframe_status\tgene_id\tgene_name\n")
        for c in range(1, 80):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            s = random.randint(10000, 50000000)
            e = s + random.randint(500, 50000)
            jid = f"{chrom}:{s}:{e}:clu_{c}_{strand}"
            exon_len = random.randint(50, 3000)
            frame = random.choice(["in_frame", "out_of_frame", "unknown"])
            gene = random.choice(gene_names)
            gid = f"ENSG{random.randint(100000, 999999):06d}.{random.randint(1,20)}"
            f.write(f"{jid}\t{exon_len}\t{frame}\t{gid}\t{gene}\n")

    # clusters_junction_classifications.txt
    with open(clustering_dir / f"{prefix}_junction_classifications.txt", "w") as f:
        f.write("junction_id\tclassification\tconfidence\tgene_id\tgene_name\ttranscript_id\n")
        for c in range(1, 100):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            s = random.randint(10000, 50000000)
            e = s + random.randint(500, 50000)
            jid = f"{chrom}:{s}:{e}:clu_{c}_{strand}"
            cls = random.choices(CLASSIFICATIONS, weights=CLASS_WEIGHTS, k=1)[0]
            conf = round(random.uniform(0.6, 1.0), 3)
            gene = random.choice(gene_names)
            gid = f"ENSG{random.randint(100000, 999999):06d}.{random.randint(1,20)}"
            tid = f"ENST{random.randint(100000, 999999):06d}.{random.randint(1,10)}"
            f.write(f"{jid}\t{cls}\t{conf}\t{gid}\t{gene}\t{tid}\n")

    # clusters_lowusage_introns
    with open(clustering_dir / f"{prefix}_lowusage_introns", "w") as f:
        for c in range(1, 25):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            s = random.randint(10000, 50000000)
            e = s + random.randint(500, 50000)
            f.write(f"{chrom}:{s}:{e}:clu_{c}_{strand}\n")


def generate_distance_file(out_path, label):
    """long_exon_distances.txt / nuc_rule_distances.txt"""
    with open(out_path, "w") as f:
        f.write(f"junction_id\t{label}_distance\tgene_id\tgene_name\ttranscript_id\n")
        gene_names = [
            "BRCA1", "TP53", "GAPDH", "ACTB", "MYC", "RB1", "PTEN", "EGFR",
            "KRAS", "VEGFA", "FOS", "JUN", "CDK2", "BCL2", "AKT1", "RAF1",
            "SRC", "HIF1A", "NOTCH1", "WNT5A",
        ]
        for c in range(1, 60):
            chrom = random.choice(CHROMS)
            strand = random.choice(["+", "-"])
            base = random.randint(10000, 50000000)
            s = base
            e = s + random.randint(500, 50000)
            jid = f"{chrom}:{s}:{e}:clu_{c}_{strand}"
            dist = random.randint(-500, 2000)
            gene = random.choice(gene_names)
            gid = f"ENSG{random.randint(100000, 999999):06d}.{random.randint(1,20)}"
            tid = f"ENST{random.randint(100000, 999999):06d}.{random.randint(1,10)}"
            f.write(f"{jid}\t{dist}\t{gid}\t{gene}\t{tid}\n")


def generate_tissue_analysis(outdir, sample_names, junction_counts_path, filelist_path):
    """Produce unproductive_by_tissue.tsv, .json, and .png."""
    tissue_data = []
    for tissue in TISSUES:
        n = SAMPLES_PER_TISSUE
        mean_pct = np.clip(np.random.normal(4.5 if tissue == "Brain_Cortex" else 3.2, 0.8), 1, 10)
        std_pct = abs(np.random.normal(0.5, 0.2))
        total_r = random.randint(800000, 2000000)
        up_r = int(total_r * mean_pct / 100)
        tissue_data.append({
            "tissue": tissue,
            "n_samples": n,
            "mean_pct": round(float(mean_pct), 4),
            "median_pct": round(float(mean_pct - 0.1), 4),
            "std_pct": round(float(std_pct), 4),
            "total_reads": total_r,
            "unproductive_reads": up_r,
        })

    import csv
    tsv_path = outdir / "unproductive_by_tissue.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tissue_data[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(tissue_data)

    all_means = [t["mean_pct"] for t in tissue_data]
    total_up = sum(t["unproductive_reads"] for t in tissue_data)
    total_all = sum(t["total_reads"] for t in tissue_data)
    global_wmean = round(total_up / total_all * 100, 4)

    json_payload = {
        "n_tissues": len(TISSUES),
        "n_samples": len(sample_names),
        "global_weighted_mean_pct": global_wmean,
        "global_median_pct": round(float(np.median(all_means)), 4),
        "tissues": tissue_data,
    }
    json_path = outdir / "unproductive_by_tissue.json"
    with open(json_path, "w") as f:
        json.dump(json_payload, f, indent=2)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4.5), dpi=150)
        x = np.arange(len(TISSUES))
        means = [t["mean_pct"] for t in tissue_data]
        stds = [t["std_pct"] for t in tissue_data]
        colors = ["#4e79a7", "#e15759"]
        ax.bar(x, means, yerr=stds, capsize=4, color=colors, edgecolor="white",
               linewidth=0.6, alpha=0.85, zorder=2)
        ax.axhline(global_wmean, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.text(len(TISSUES) - 0.5, global_wmean, f"  mean = {global_wmean:.2f}%",
                va="bottom", ha="right", fontsize=7, color="grey")
        ax.set_xticks(x)
        ax.set_xticklabels(TISSUES, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("% unproductive junction reads")
        ax.set_title("Unproductive splicing by tissue", fontsize=11, fontweight="bold")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", linestyle="--", alpha=0.3, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        png_path = outdir / "unproductive_by_tissue.png"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  PNG: {png_path}")
    except ImportError:
        print("  [WARN] matplotlib not available; skipping PNG generation")

    return {
        "png": "unproductive_by_tissue.png",
        "tsv": "unproductive_by_tissue.tsv",
        "json": "unproductive_by_tissue.json",
        "n_tissues": len(TISSUES),
        "n_samples": len(sample_names),
        "global_weighted_mean_pct": global_wmean,
        "global_median_pct": json_payload["global_median_pct"],
    }


def generate_up_transcript_counts(outdir):
    """Produce synthetic up_transcript_counts.tsv and .png."""
    import csv
    timepoints = TISSUES
    rows = []
    for tissue in timepoints:
        n = SAMPLES_PER_TISSUE
        mean_juncs = np.clip(np.random.normal(45 if tissue == "Brain_Cortex" else 32, 8), 10, 100)
        mean_genes = np.clip(mean_juncs * np.random.uniform(0.55, 0.75), 5, 80)
        rows.append({
            "timepoint": tissue,
            "n_samples": n,
            "mean_up_junctions": round(float(mean_juncs), 1),
            "std_up_junctions": round(abs(np.random.normal(5, 2)), 1),
            "mean_up_genes": round(float(mean_genes), 1),
            "std_up_genes": round(abs(np.random.normal(3, 1.5)), 1),
        })

    tsv_path = outdir / "up_transcript_counts.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n_tp = len(rows)
        x = np.arange(n_tp)
        bar_w = 0.35
        fig, ax = plt.subplots(figsize=(max(6, n_tp * 1.4 + 2), 5), dpi=150)

        ax.bar(x - bar_w / 2, [r["mean_up_junctions"] for r in rows], bar_w,
               yerr=[r["std_up_junctions"] for r in rows], capsize=3,
               color="#4e79a7", edgecolor="white", linewidth=0.5, alpha=0.88,
               label="UP junctions")
        ax.bar(x + bar_w / 2, [r["mean_up_genes"] for r in rows], bar_w,
               yerr=[r["std_up_genes"] for r in rows], capsize=3,
               color="#e15759", edgecolor="white", linewidth=0.5, alpha=0.88,
               label="Unique genes")

        ax.set_xticks(x)
        ax.set_xticklabels([r["timepoint"] for r in rows], rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Count (mean per sample)", fontsize=9.5)
        ax.set_title("UP transcripts detected per tissue", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle="--", alpha=0.25, zorder=0)
        fig.tight_layout()
        png_path = outdir / "up_transcript_counts.png"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  PNG: {png_path}")
    except ImportError:
        print("  [WARN] matplotlib not available; skipping PNG generation")

    return {
        "tsv": "up_transcript_counts.tsv",
        "png": "up_transcript_counts.png",
        "timepoints": [r["timepoint"] for r in rows],
        "mean_up_junctions": [r["mean_up_junctions"] for r in rows],
        "mean_up_genes": [r["mean_up_genes"] for r in rows],
    }


def generate_pipeline_log(log_path, sample_names):
    lines = [
        ">> python3 scripts/lc2_pipeline.py --workdir example_output --prefix gtex_demo ...",
        "Using 6 pre-made junction BED files",
        "",
        ">> python3 tools/leafcutter2/scripts/leafcutter_make_clusters.py -j clusters/gtex_demo_junction_files.txt -m 50 -o gtex_demo -r clusters/ -l 500000",
        "Sorting...",
        "Merging...",
        "Refining clusters...",
        "Wrote 120 clusters",
        "",
        ">> python3 tools/leafcutter2/scripts/leafcutter2.py -j clusters/gtex_demo_junction_files.txt -r lc2 -o gtex_demo_lc2 -A refs/gencode.v46.annotation.gtf -G refs/GRCh38.fa -m 50 -l 500000 --keepannot",
        "Loading genome refs/GRCh38.fa ...done!",
        "Classifying splice junctions...",
        f"Loading annotation from refs/gencode.v46.annotation.gtf ...",
        "2026-03-18 14:22:01 pool_junc_reads ...",
        "2026-03-18 14:22:03 Refine clusters ...",
        "Sorting lc2/clustering/gtex_demo_lc2_sortedlibs ...",
        "2026-03-18 14:22:15 Merging ...",
        "Extracting numerators ...",
        " ... 10000 introns annotated.",
        " ... 20000 introns annotated.",
        "2026-03-18 14:25:33 Annotation done.",
        "Annotated 24317 introns.",
        "Filtered out 1823 introns with stdev(reads) < 0.5 or zero usage.",
        "",
        "2026-03-18 14:25:34 Done.",
        "",
        f"[TISSUE] 2 tissues, weighted mean = 3.85%",
        f"[TISSUE] TSV:  example_output/out/unproductive_by_tissue.tsv",
        f"[TISSUE] PNG:  example_output/out/unproductive_by_tissue.png",
        f"[TISSUE] JSON: example_output/out/unproductive_by_tissue.json",
        "",
        "Done.",
        f"- LC2 cluster ratios: example_output/lc2/{LC2_PREFIX}.cluster_ratios.gz",
        f"- LC2 junction counts: example_output/lc2/{LC2_PREFIX}.junction_counts.gz",
        f"- Summary: example_output/out/summary.json",
    ]
    with open(log_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    print(f"Generating example outputs in {BASE}/\n")

    jdir = ensure(BASE / "junctions")
    cdir = ensure(BASE / "clusters")
    lc2dir = ensure(BASE / "lc2")
    lc2_clustering = ensure(lc2dir / "clustering")
    outdir = ensure(BASE / "out")

    # 1. Junction BEDs
    print("[1/7] Junction BED files...")
    bed_paths = generate_junction_beds(jdir)

    # 2. Junction filelist
    print("[2/7] Junction file list...")
    filelist = generate_junction_filelist(bed_paths, cdir / f"{PREFIX}_junction_files.txt")

    # 3. Clustering outputs (perind)
    all_sample_names = []
    for tissue in TISSUES:
        for sid in SAMPLE_IDS[tissue]:
            all_sample_names.append(sid)

    print("[3/7] Clustering perind counts...")
    generate_perind_counts(all_sample_names, cdir / f"{PREFIX}_perind.counts.gz")
    generate_perind_counts(all_sample_names, cdir / f"{PREFIX}_perind_numers.counts.gz")

    # 4. LC2 classification outputs
    print("[4/7] LC2 cluster_ratios.gz...")
    generate_cluster_ratios(all_sample_names, lc2dir / f"{LC2_PREFIX}.cluster_ratios.gz")

    print("[5/7] LC2 junction_counts.gz...")
    jc_path = lc2dir / f"{LC2_PREFIX}.junction_counts.gz"
    generate_junction_counts(all_sample_names, jc_path)

    print("[6/8] Distance files + clustering intermediates...")
    generate_distance_file(lc2_clustering / f"{LC2_PREFIX}_long_exon_distances.txt", "long_exon")
    generate_distance_file(lc2_clustering / f"{LC2_PREFIX}_nuc_rule_distances.txt", "nuc_rule")
    generate_clustering_intermediates(lc2_clustering, LC2_PREFIX)

    # 5. Tissue analysis
    print("[7/9] Tissue unproductive analysis...")
    tissue_metrics = generate_tissue_analysis(outdir, all_sample_names, jc_path, filelist)

    # 5b. UP transcript counts
    print("[8/9] UP transcript counts...")
    up_count_metrics = generate_up_transcript_counts(outdir)

    # 6. Summary JSON
    summary = {
        "n_junction_inputs": len(bed_paths),
        "junction_file_list": str(filelist),
        "lc2_outputs": {
            "cluster_ratios": str(lc2dir / f"{LC2_PREFIX}.cluster_ratios.gz"),
            "junction_counts": str(jc_path),
            "long_exon_distances": str(lc2_clustering / f"{LC2_PREFIX}_long_exon_distances.txt"),
            "nuc_rule_distances": str(lc2_clustering / f"{LC2_PREFIX}_nuc_rule_distances.txt"),
            "run_dir": str(lc2dir),
        },
        "cluster_files": {
            "perind_counts": str(cdir / f"{PREFIX}_perind.counts.gz"),
            "perind_numers": str(cdir / f"{PREFIX}_perind_numers.counts.gz"),
        },
        "metrics": {},
        "notes": "Core pipeline run complete. Use long_exon/nuc_rule outputs for NMD-focused downstream analysis.",
    }
    if tissue_metrics:
        summary["metrics"]["unproductive_by_tissue"] = tissue_metrics
    if up_count_metrics:
        summary["metrics"]["up_transcript_counts"] = up_count_metrics

    with open(outdir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # 7. Pipeline log
    print("[9/9] Pipeline log...")
    generate_pipeline_log(BASE / "pipeline.log", all_sample_names)

    print(f"\nDone! Example outputs written to: {BASE}/")
    print(f"\nDirectory structure:")
    for root, dirs, files in os.walk(BASE):
        level = root.replace(str(BASE), "").count(os.sep)
        indent = "  " * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = "  " * (level + 1)
        for file in sorted(files):
            size = os.path.getsize(os.path.join(root, file))
            if size > 1024:
                size_str = f"{size/1024:.1f} KB"
            else:
                size_str = f"{size} B"
            print(f"{subindent}{file}  ({size_str})")


if __name__ == "__main__":
    main()
