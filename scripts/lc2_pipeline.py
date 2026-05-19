#!/usr/bin/env python3
"""
LeafCutter2 end-to-end pipeline (Python-only launcher)
- Works on macOS / Linux / Windows (Anaconda Prompt)
- No bash scripting required
"""

import argparse, gzip, os, re, sys, shutil, subprocess, json, tempfile
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from glob import glob

import numpy as np
import pandas as pd
from tqdm import tqdm

# ---------------------------
# GTF biotype normalization
# ---------------------------

# NCBI GTFs use different biotype values than Ensembl.
# LeafCutter2 expects Ensembl-style values (e.g. "protein_coding",
# "nonsense_mediated_decay").  This map covers NCBI→Ensembl for the
# biotypes that affect classification.
_NCBI_TO_ENSEMBL_BIOTYPE = {
    "mRNA": "protein_coding",
    "lnc_RNA": "lncRNA",
    "primary_transcript": "processed_transcript",
}


def normalize_gtf_for_lc2(gtf: Path, workdir: Path) -> Path:
    """Return a GTF whose transcript_biotype values match Ensembl conventions.

    If the GTF already uses Ensembl-style biotypes, return it unchanged.
    Otherwise create a normalised copy under *workdir* and return that path.
    Works for any organism / genome build (GRCz11, GRCz12, etc.).
    """
    needs_fix = False
    with open(gtf, "r") as fh:
        for i, line in enumerate(fh):
            if i >= 500:
                break
            if 'transcript_biotype "mRNA"' in line:
                needs_fix = True
                break

    if not needs_fix:
        return gtf

    out = workdir / (gtf.stem + ".lc2_compat.gtf")
    if out.exists() and out.stat().st_size > 0:
        print(f"[GTF] Reusing normalised GTF: {out}", flush=True)
        return out

    print(f"[GTF] NCBI-style biotypes detected — creating LeafCutter2-compatible GTF ...", flush=True)
    tmp = out.with_suffix(".tmp")
    with open(gtf, "r") as fin, open(tmp, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            for ncbi_val, ensembl_val in _NCBI_TO_ENSEMBL_BIOTYPE.items():
                line = line.replace(
                    f'transcript_biotype "{ncbi_val}"',
                    f'transcript_biotype "{ensembl_val}"',
                )
                line = line.replace(
                    f'gene_biotype "{ncbi_val}"',
                    f'gene_biotype "{ensembl_val}"',
                )
            fout.write(line)
    tmp.rename(out)
    print(f"[GTF] Wrote {out}", flush=True)
    return out


# ---------------------------
# Helpers
# ---------------------------

def check_exe(name: str):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required executable '{name}' not found on PATH. Install it or activate the right conda env.")
    return path

def run(cmd: List[str], workdir: Path = None):
    print(">>", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(workdir) if workdir else None, check=True)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    return p

def write_junction_filelist(paths: List[Path], out_path: Path) -> Path:
    with open(out_path, "w") as fh:
        for p in paths:
            fh.write(str(p) + "\n")
    return out_path

# ---------------------------
# Junction extraction
# ---------------------------

def extract_junctions_from_bams(bams: List[Path], out_dir: Path) -> List[Path]:
    """regtools junctions extract for each BAM -> BED"""
    check_exe("regtools")
    ensure_dir(out_dir)
    out_beds = []
    for bam in tqdm(bams, desc="Extracting junctions with regtools"):
        if not bam.exists():
            raise FileNotFoundError(bam)
        bed = out_dir / (bam.stem + ".juncs.bed")
        cmd = ["regtools", "junctions", "extract", "-o", str(bed), str(bam)]
        run(cmd)
        out_beds.append(bed)
    return out_beds

def convert_star_sj_to_bed(sj_tabs: List[Path], out_dir: Path) -> List[Path]:
    """
    Convert STAR SJ.out.tab -> BED12 (regtools-compatible) for leafcutter_cluster_regtools.py.

    The clustering script unpacks 12 columns and derives intron coords via:
        Aoff, Boff = blockSize.split(",")
        A, B = chromStart + Aoff, chromEnd - Boff + 1

    We use 1-bp exonic anchor blocks flanking each intron so the math recovers the
    correct 0-based intron start (A) and 1-based intron end (B).
    """
    ensure_dir(out_dir)
    beds = []
    for sj in tqdm(sj_tabs, desc="Converting STAR SJ.out.tab to BED12"):
        if not sj.exists():
            raise FileNotFoundError(sj)
        df = pd.read_csv(sj, sep="\t", header=None, comment="#")
        if df.shape[1] < 7:
            raise ValueError(f"{sj} appears malformed; expected >=7 columns.")
        intron_start_1 = df[1].astype(int)  # 1-based intron start
        intron_end_1   = df[2].astype(int)  # 1-based intron end
        strand = df[3].map({0: ".", 1: "+", 2: "-"}).fillna(".")
        score  = df[6].astype(int)

        chrom_start = intron_start_1 - 2       # 0-based; 1-bp anchor before intron
        chrom_end   = intron_end_1              # half-open end encompassing 1-bp anchor after intron
        block_start_2 = intron_end_1 - intron_start_1 + 1

        bed = pd.DataFrame({
            "chrom":       df[0],
            "chromStart":  chrom_start,
            "chromEnd":    chrom_end,
            "name":        "JUNC",
            "score":       score,
            "strand":      strand,
            "thickStart":  chrom_start,
            "thickEnd":    chrom_end,
            "rgb":         "255,0,0",
            "blockCount":  2,
            "blockSizes":  "1,1",
            "blockStarts": "0," + block_start_2.astype(str),
        })
        out = out_dir / (sj.parent.name + "_" + sj.stem + ".bed")
        bed.to_csv(out, sep="\t", header=False, index=False)
        beds.append(out)
    return beds

# ---------------------------
# LeafCutter clustering
# ---------------------------

def _find_cluster_script(leafcutter_repo: Path, leafcutter2_repo: Optional[Path] = None) -> Optional[Path]:
    """Search for a compatible LeafCutter clustering script."""
    candidates = [
        leafcutter_repo / "clustering" / "leafcutter_cluster_regtools.py",
        leafcutter_repo / "scripts" / "leafcutter_cluster_regtools.py",
    ]
    if leafcutter2_repo:
        candidates.append(leafcutter2_repo / "scripts" / "leafcutter_make_clusters.py")
    for c in candidates:
        if c.exists():
            return c
    return None


def leafcutter_cluster(
    leafcutter_repo: Path,
    bed_paths: List[Path],
    out_prefix: Path,
    min_reads=50,
    max_intron_len=500000,
    leafcutter2_repo: Optional[Path] = None,
):
    """Call LeafCutter clustering script (Python), using regtools BEDs list.

    Searches for the clustering script in the LC1 repo first, then falls back
    to the LC2 repo's ``leafcutter_make_clusters.py`` which has a compatible
    CLI.  Returns ``None`` if no clustering script is found.
    """
    cluster_script = _find_cluster_script(leafcutter_repo, leafcutter2_repo)
    if cluster_script is None:
        return None

    filelist = write_junction_filelist(bed_paths, out_prefix.parent / (out_prefix.name + "_junction_files.txt"))

    cmd = [
        sys.executable, str(cluster_script),
        "-j", str(filelist),
        "-m", str(min_reads),
        "-o", out_prefix.name,
        "-r", str(out_prefix.parent) + "/",
        "-l", str(max_intron_len),
    ]
    run(cmd)

    perind = out_prefix.with_name(out_prefix.name + "_perind.counts.gz")
    perind_num = out_prefix.with_name(out_prefix.name + "_perind_numers.counts.gz")
    if not perind.exists() or not perind_num.exists():
        print("[WARN] LeafCutter clustering did not produce expected perind files — continuing without them")
        return filelist, None, None
    return filelist, perind, perind_num

# ---------------------------
# LeafCutter2 classification
# ---------------------------

def _cache_annotation_pickles(gtf: Path, run_dir: Path, *, restore: bool):
    """Copy cached annotation pickles between a shared location (next to the GTF)
    and the per-job run_dir so repeated runs skip the slow GTF parse."""
    gtf_stem = gtf.name.split(".gtf")[0]
    pickle_names = [
        f"{gtf_stem}_SJC_annotations.pckle",
        f"txn2gene.{gtf_stem}_SJC_annotations.pckle",
        f"nmd_txn2gene.{gtf_stem}_SJC_annotations.pckle",
    ]
    cache_dir = gtf.parent
    for name in pickle_names:
        shared = cache_dir / name
        local = run_dir / name
        if restore:
            if shared.exists() and not local.exists():
                shutil.copy2(shared, local)
        else:
            if local.exists() and not shared.exists():
                shutil.copy2(local, shared)


def leafcutter2_classify(
    leafcutter2_repo: Path,
    fasta: Path,
    gtf: Path,
    junction_filelist: Path,
    run_dir: Path,
    out_prefix: str,
    min_cluster_reads: int = 30,
    max_intron_len: int = 100000,
    gene_type_tag: str = "gene_type",
    transcript_type_tag: str = "transcript_type",
    gene_name_tag: str = "gene_name",
    transcript_name_tag: str = "transcript_name",
) -> Dict[str, str]:
    lc2 = leafcutter2_repo / "leafcutter2.py"
    if not lc2.exists():
        # sometimes the script is in scripts/
        lc2 = leafcutter2_repo / "scripts" / "leafcutter2.py"
    if not lc2.exists():
        raise FileNotFoundError("leafcutter2.py not found in the provided repo path.")

    _cache_annotation_pickles(gtf, run_dir, restore=True)

    cmd = [
        sys.executable, str(lc2),
        "-j", str(junction_filelist),
        "-r", str(run_dir),
        "-o", out_prefix,
        "-A", str(gtf),
        "-G", str(fasta),
        "-m", str(min_cluster_reads),
        "-l", str(max_intron_len),
        "--keepannot",
        "--gene_type", str(gene_type_tag),
        "--transcript_type", str(transcript_type_tag),
        "--gene_name", str(gene_name_tag),
        "--transcript_name", str(transcript_name_tag),
    ]
    run(cmd)

    _cache_annotation_pickles(gtf, run_dir, restore=False)
    expected = [
        run_dir / f"{out_prefix}.cluster_ratios.gz",
        run_dir / f"{out_prefix}.junction_counts.gz",
        run_dir / "clustering" / f"{out_prefix}_long_exon_distances.txt",
        run_dir / "clustering" / f"{out_prefix}_nuc_rule_distances.txt",
    ]
    missing = [str(p) for p in expected if not p.exists()]
    if missing:
        raise RuntimeError("LeafCutter2 completed but expected outputs are missing:\n  - " + "\n  - ".join(missing))
    return {
        "cluster_ratios": str(expected[0]),
        "junction_counts": str(expected[1]),
        "long_exon_distances": str(expected[2]),
        "nuc_rule_distances": str(expected[3]),
        "run_dir": str(run_dir),
    }

# ---------------------------
# Defaults + parser (Option 1)
# ---------------------------

VALID_CLASSIFICATIONS = {"UP", "PR", "NE", "IN"}

_STAGE_PATTERN = re.compile(
    r"\d+(hpf|dpf|cell|somite)"
    r"|dome|shield|bud|sphere|oblong|gastrula|pharyngula|segmentation"
    r"|blastula|cleavage|epiboly|hatching|larval",
    re.IGNORECASE,
)


def _is_developmental_stage_labels(labels) -> bool:
    """Return True if any group label looks like a developmental timepoint."""
    return any(_STAGE_PATTERN.fullmatch(str(l).strip()) for l in labels)


def _group_noun(labels) -> str:
    """Return 'developmental stage' or 'tissue' depending on the label style."""
    return "developmental stage" if _is_developmental_stage_labels(labels) else "tissue"


def _timepoint_to_hours(label: str) -> float:
    """Convert a developmental-stage label to hours for chronological sorting.

    Handles 'hpf' (hours post-fertilisation) and 'dpf' (days) labels.
    Returns ``float('inf')`` for unrecognised strings so they sort last.
    """
    s = str(label).strip().lower()
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(hpf|dpf)", s)
    if m:
        val = float(m.group(1))
        return val if m.group(2) == "hpf" else val * 24.0
    return float("inf")


def _spearman_rho(x, y):
    """Compute Spearman rank correlation and approximate two-sided p-value.

    Pure-numpy implementation so scipy is not required on Quest.
    Returns ``(rho, p_value)`` or ``(None, None)`` if fewer than 3 points.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 3:
        return None, None

    def _rankdata(a):
        order = np.argsort(a)
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, n + 1, dtype=float)
        return ranks

    rx, ry = _rankdata(x), _rankdata(y)
    mx, my = rx.mean(), ry.mean()
    dx, dy = rx - mx, ry - my
    denom = np.sqrt((dx ** 2).sum() * (dy ** 2).sum())
    if denom == 0:
        return 0.0, 1.0
    rho = float((dx * dy).sum() / denom)
    t_stat = rho * np.sqrt((n - 2) / max(1.0 - rho ** 2, 1e-15))
    # Two-tailed p via regularised incomplete-beta (Abramowitz & Stegun 26.5.27)
    df = n - 2
    x_beta = df / (df + t_stat ** 2)
    p = _betai(df / 2.0, 0.5, x_beta)
    return rho, p


def _betai(a, b, x):
    """Regularised incomplete beta function via continued-fraction expansion."""
    if x < 0 or x > 1:
        return 1.0
    if x == 0 or x == 1:
        return 1.0 - x
    from math import lgamma, exp, log
    lbeta = lgamma(a) + lgamma(b) - lgamma(a + b)
    front = exp(a * log(x) + b * log(1.0 - x) - lbeta) / a
    # Lentz continued-fraction (modified)
    if x < (a + 1) / (a + b + 2):
        return front * _betacf(a, b, x)
    else:
        return 1.0 - exp(a * log(x) + b * log(1.0 - x) - lbeta) / b * _betacf(b, a, 1.0 - x)


def _betacf(a, b, x, max_iter=200, eps=3e-12):
    """Continued-fraction evaluation for the incomplete beta function."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = max(1.0 - qab * x / qap, 1e-30)
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = max(1.0 + aa * d, 1e-30)
        c = max(1.0 + aa / c, 1e-30)
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = max(1.0 + aa * d, 1e-30)
        c = max(1.0 + aa / c, 1e-30)
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _bed_path_to_col_name(bed_path: str) -> str:
    """Replicate the column-naming logic from leafcutter_cluster_regtools.py."""
    return Path(bed_path).name.split(".junc")[0]


def _resolve_tissue_labels(
    sample_cols: List[str],
    junction_filelist_path: Path,
    samples_tsv_path: Optional[Path] = None,
) -> List[str]:
    """Map junction_counts column names to tissue/group labels.

    Strategy (in priority order):
      1. BED parent directory names from the junction filelist (GTEx-style).
      2. ``condition`` column from samples_tsv (non-GTEx with metadata).
      3. Use each sample name as its own group (fallback).
    """
    col_to_tissue: Dict[str, str] = {}

    if junction_filelist_path.exists():
        bed_paths = [
            l.strip()
            for l in junction_filelist_path.read_text().splitlines()
            if l.strip()
        ]
        for bp in bed_paths:
            p = Path(bp)
            tissue = p.parent.name
            col_name = _bed_path_to_col_name(bp)
            col_to_tissue[col_name] = tissue
            col_to_tissue[p.name] = tissue
            col_to_tissue[p.stem] = tissue

    tissues = [col_to_tissue.get(c, col_to_tissue.get(Path(c).name, None)) for c in sample_cols]

    unique = set(t for t in tissues if t is not None)
    all_resolved = all(t is not None for t in tissues)
    single_group = len(unique) <= 1

    if (not all_resolved or single_group) and samples_tsv_path and samples_tsv_path.exists():
        meta = pd.read_csv(samples_tsv_path, sep="\t", dtype=str)
        if {"sample", "condition"}.issubset(meta.columns):
            cond_map: Dict[str, str] = {}
            sample_ids: List[str] = []
            for _, row in meta.iterrows():
                s = str(row["sample"])
                c = str(row["condition"])
                sample_ids.append(s)
                cond_map[s] = c
                cond_map[Path(s).name] = c
                cond_map[Path(s).stem] = c
                for suf in [".juncs.bed", ".juncs", ".bed", ".sorted.gz"]:
                    if s.endswith(suf):
                        cond_map[s[: -len(suf)]] = c

            def _lookup(col: str) -> str:
                hit = cond_map.get(col) or cond_map.get(Path(col).name) or cond_map.get(Path(col).stem)
                if hit:
                    return hit
                stripped = col
                for suf in [".SJ.out.bed", ".SJ.out", ".bed", ".sorted.gz"]:
                    if stripped.endswith(suf):
                        stripped = stripped[: -len(suf)]
                        break
                hit = cond_map.get(stripped)
                if hit:
                    return hit
                for sid in sample_ids:
                    if sid in stripped:
                        return cond_map[sid]
                return col

            tissues = [_lookup(c) for c in sample_cols]
            return tissues

    # Fallback: extract developmental stage from column name prefix
    # (e.g. "24hpf_ERR738214.SJ.out.bed" -> "24hpf")
    if single_group or not all_resolved:
        _run_id_pat = re.compile(r"[_-](ERR|SRR|DRR|SAMN|SAME)\d+", re.IGNORECASE)
        extracted = []
        for c in sample_cols:
            m = _run_id_pat.search(c)
            prefix = c[:m.start()] if m else c.split("_")[0]
            extracted.append(prefix.strip("_.- "))
        if len(set(extracted)) > 1 and _is_developmental_stage_labels(set(extracted)):
            return extracted

    return [t if t is not None else c for t, c in zip(tissues, sample_cols)]


def compute_unproductive_by_tissue(
    junction_counts_path: Path,
    junction_filelist_path: Path,
    outdir: Path,
    samples_tsv_path: Optional[Path] = None,
) -> Optional[Dict]:
    """Compute per-tissue unproductive junction read percentages and generate
    a bar-chart PNG, a TSV, and a JSON summary.

    Uses the LC2 classification labels (UP/PR/NE/IN) embedded in
    ``junction_counts.gz`` junction IDs.  Returns metadata dict for
    inclusion in ``summary.json``, or ``None`` on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib not available — skipping unproductive-by-tissue plot")
        return None

    if not junction_counts_path.exists():
        print(f"[WARN] junction_counts not found at {junction_counts_path}")
        return None

    # ---- 1. Parse junction_counts.gz ----
    with gzip.open(str(junction_counts_path), "rt") as f:
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
        print("[WARN] junction_counts has no UP/PR/NE/IN labels — classification may not have run")
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_per_sample = np.where(
            total_reads > 0, up_reads / total_reads * 100.0, np.nan
        )

    # ---- 2. Resolve tissue labels ----
    tissues = _resolve_tissue_labels(
        sample_cols, junction_filelist_path, samples_tsv_path
    )

    # ---- 3. Aggregate by tissue ----
    tissue_samples: Dict[str, List[float]] = defaultdict(list)
    for i, tissue in enumerate(tissues):
        if not np.isnan(pct_per_sample[i]):
            tissue_samples[tissue].append(pct_per_sample[i])

    if not tissue_samples:
        return None

    rows = []
    for tissue in sorted(tissue_samples):
        vals = np.array(tissue_samples[tissue])
        idxs = [i for i, t in enumerate(tissues) if t == tissue]
        rows.append(
            {
                "tissue": tissue,
                "n_samples": len(vals),
                "mean_pct": round(float(np.mean(vals)), 4),
                "median_pct": round(float(np.median(vals)), 4),
                "std_pct": round(float(np.std(vals, ddof=1)), 4) if len(vals) > 1 else 0.0,
                "total_reads": int(sum(total_reads[i] for i in idxs)),
                "unproductive_reads": int(sum(up_reads[i] for i in idxs)),
            }
        )

    result_df = pd.DataFrame(rows)

    # Sort chronologically when labels are developmental stages, else by mean_pct
    is_dev = _is_developmental_stage_labels(result_df["tissue"].values)
    if is_dev:
        result_df["_hours"] = result_df["tissue"].apply(_timepoint_to_hours)
        result_df = result_df.sort_values("_hours").reset_index(drop=True)
        result_df.drop(columns="_hours", inplace=True)
    else:
        result_df = result_df.sort_values("mean_pct").reset_index(drop=True)

    # Global statistics
    all_pcts = np.array([v for vs in tissue_samples.values() for v in vs])
    global_weighted_mean = float(np.sum(up_reads[~np.isnan(pct_per_sample)]) /
                                 np.sum(total_reads[~np.isnan(pct_per_sample)]) * 100.0)
    global_median = float(np.median(all_pcts))
    n_tissues = len(tissue_samples)

    # Spearman correlation of mean % vs developmental time
    spearman_rho, spearman_p = None, None
    if is_dev and n_tissues >= 3:
        hours = np.array([_timepoint_to_hours(t) for t in result_df["tissue"]])
        finite = np.isfinite(hours)
        if finite.sum() >= 3:
            spearman_rho, spearman_p = _spearman_rho(
                hours[finite], result_df["mean_pct"].values[finite]
            )

    # ---- 4. Write TSV ----
    tsv_path = outdir / "unproductive_by_tissue.tsv"
    result_df.to_csv(tsv_path, sep="\t", index=False)

    # ---- 5. Write JSON ----
    json_path = outdir / "unproductive_by_tissue.json"
    json_payload = {
        "n_tissues": n_tissues,
        "n_samples": int(np.sum(~np.isnan(pct_per_sample))),
        "global_weighted_mean_pct": round(global_weighted_mean, 4),
        "global_median_pct": round(global_median, 4),
        "tissues": rows,
    }
    if spearman_rho is not None:
        json_payload["spearman_rho"] = round(spearman_rho, 4)
        json_payload["spearman_p"] = round(spearman_p, 4)
    with open(json_path, "w") as fh:
        json.dump(json_payload, fh, indent=2)

    # ---- 6. Plot ----
    png_path = outdir / "unproductive_by_tissue.png"
    _plot_unproductive_chart(
        result_df, tissue_samples, global_weighted_mean, png_path,
        spearman_rho=spearman_rho, spearman_p=spearman_p,
    )

    group = _group_noun(result_df["tissue"].values)
    print(f"[GROUP] {n_tissues} {group}s, weighted mean = {global_weighted_mean:.3f}%")
    if spearman_rho is not None:
        print(f"[GROUP] Spearman rho = {spearman_rho:.4f}, p = {spearman_p:.4f}")
    print(f"[GROUP] TSV:  {tsv_path}")
    print(f"[GROUP] PNG:  {png_path}")
    print(f"[GROUP] JSON: {json_path}")

    ret: Dict = {
        "png": "unproductive_by_tissue.png",
        "tsv": "unproductive_by_tissue.tsv",
        "json": "unproductive_by_tissue.json",
        "n_tissues": n_tissues,
        "n_samples": int(np.sum(~np.isnan(pct_per_sample))),
        "global_weighted_mean_pct": round(global_weighted_mean, 4),
        "global_median_pct": round(global_median, 4),
    }
    if spearman_rho is not None:
        ret["spearman_rho"] = round(spearman_rho, 4)
        ret["spearman_p"] = round(spearman_p, 4)
    return ret


def _plot_unproductive_chart(
    result_df: pd.DataFrame,
    tissue_samples: Dict[str, List[float]],
    global_mean: float,
    out_path: Path,
    *,
    spearman_rho: Optional[float] = None,
    spearman_p: Optional[float] = None,
) -> None:
    """Render a sorted bar chart of per-tissue unproductive junction read %."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(result_df)
    bar_width = 0.6 if n <= 3 else 0.8
    fig_width = max(5, min(n * 0.7 + 2.5, 18))
    fig_height = 5
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=150)

    x = np.arange(n)
    palette = ["#4e79a7", "#e15759", "#59a14f", "#f28e2b", "#b07aa1",
               "#76b7b2", "#ff9da7", "#9c755f", "#edc948", "#bab0ac"]
    if n <= len(palette):
        colors = palette[:n]
    else:
        cmap = plt.get_cmap("tab20" if n <= 20 else "gist_rainbow")
        colors = [cmap(i / max(n - 1, 1)) for i in range(n)]

    means = result_df["mean_pct"].values
    stds = result_df["std_pct"].values
    labels = result_df["tissue"].values

    ax.bar(x, means, width=bar_width, yerr=stds, capsize=4, color=colors,
           edgecolor="white", linewidth=0.6, alpha=0.88, zorder=2,
           error_kw=dict(lw=1.2, capthick=1.0))

    jitter_half = bar_width * 0.35
    for i, tissue in enumerate(labels):
        vals = tissue_samples.get(tissue, [])
        if vals and len(vals) > 1:
            jitter = np.random.default_rng(42).uniform(
                -jitter_half, jitter_half, size=len(vals),
            )
            ax.scatter(
                np.full(len(vals), i) + jitter,
                vals,
                s=14, alpha=0.55, color="#333333", zorder=3, linewidths=0,
            )

    ax.axhline(global_mean, color="#888888", linestyle="--", linewidth=0.9,
               alpha=0.6, zorder=1)
    ax.annotate(
        f"mean = {global_mean:.2f}%",
        xy=(n - 0.5, global_mean), xytext=(6, 4),
        textcoords="offset points", fontsize=7.5, color="#666666",
        ha="right", va="bottom",
    )

    if spearman_rho is not None:
        sig = "*" if spearman_p is not None and spearman_p < 0.05 else ""
        p_str = f"{spearman_p:.3f}" if spearman_p is not None else "n/a"
        ax.text(
            0.98, 0.95,
            f"Spearman \u03c1 = {spearman_rho:.3f}{sig}\np = {p_str}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="#444444",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f0f0f0", ec="#cccccc",
                      alpha=0.85),
        )

    rotation = 0 if n <= 4 else (35 if n <= 10 else 55)
    ha = "center" if rotation == 0 else "right"
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=rotation, ha=ha, fontsize=8.5)
    ax.set_ylabel("% unproductive junction reads", fontsize=9.5)

    group = _group_noun(labels)
    subtitle = (f"{n} {group}{'s' if n != 1 else ''}  \u00b7  "
                f"weighted mean = {global_mean:.2f}%  \u00b7  "
                f"median = {np.median(means):.2f}%")
    ax.set_title(
        f"Unproductive splicing by {group}\n",
        fontsize=12, fontweight="bold", pad=4,
    )
    ax.text(
        0.5, 1.02, subtitle,
        transform=ax.transAxes, ha="center", fontsize=8, color="#888888",
    )

    margin = bar_width / 2 + 0.3
    ax.set_xlim(-margin, n - 1 + margin)

    all_vals = [v for vs in tissue_samples.values() for v in vs]
    if all_vals:
        data_min = min(min(all_vals), min(means - stds))
        data_max = max(max(all_vals), max(means + stds))
        pad = (data_max - data_min) * 0.15 or 0.2
        ax.set_ylim(max(0, data_min - pad), data_max + pad)
    else:
        ax.set_ylim(bottom=0)

    ax.grid(axis="y", linestyle="--", alpha=0.25, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=8.5)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------
# UP transcript / gene count per timepoint
# ----------------------------------------------------------

def compute_up_transcript_counts(
    junction_counts_path: Path,
    classification_dir: Path,
    lc2_prefix: str,
    junction_filelist_path: Path,
    outdir: Path,
    samples_tsv_path: Optional[Path] = None,
) -> Optional[Dict]:
    """Count unique UP junctions and genes detected per developmental timepoint.

    For each sample, a UP junction is "detected" if its read count > 0.
    Results are aggregated as mean +/- std across samples within each timepoint.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("[WARN] matplotlib not available — skipping UP transcript count plot")
        return None

    if not junction_counts_path.exists():
        print(f"[WARN] junction_counts not found at {junction_counts_path}")
        return None

    classif_path = classification_dir / f"{lc2_prefix}_junction_classifications.txt"

    # ---- 1. Build coord -> gene mapping ----
    coord_to_gene: Dict[str, str] = {}
    if classif_path.exists():
        with open(classif_path) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    coord_to_gene[parts[1]] = parts[0]

    # ---- 2. Parse junction_counts.gz — collect UP junction read vectors ----
    up_junc_reads: Dict[str, np.ndarray] = {}
    with gzip.open(str(junction_counts_path), "rt") as f:
        header_tokens = f.readline().strip().split()
        sample_cols = header_tokens[1:]
        n = len(sample_cols)
        if n == 0:
            return None

        for line in f:
            fields = line.strip().split()
            if len(fields) < 2:
                continue
            junc_id = fields[0]
            label = junc_id.rsplit(":", 1)[-1]
            if label != "UP":
                continue
            reads = np.array([float(x) for x in fields[1 : n + 1]])
            toks = junc_id.split(":")
            coord = f"{toks[0]}:{toks[1]}-{toks[2]}"
            up_junc_reads[coord] = reads

    if not up_junc_reads:
        print("[WARN] No UP junctions found — skipping transcript count analysis")
        return None

    # ---- 3. Resolve timepoint labels ----
    tissues = _resolve_tissue_labels(sample_cols, junction_filelist_path, samples_tsv_path)
    unique_tissues = sorted(set(tissues), key=lambda t: _timepoint_to_hours(t))
    tissue_idx: Dict[str, List[int]] = defaultdict(list)
    for i, t in enumerate(tissues):
        tissue_idx[t].append(i)

    # ---- 4. Per-sample counts of detected UP junctions and genes ----
    rows = []
    for tp in unique_tissues:
        idxs = tissue_idx[tp]
        sample_junc_counts = []
        sample_gene_counts = []
        for si in idxs:
            detected_juncs = set()
            detected_genes = set()
            for coord, reads in up_junc_reads.items():
                if reads[si] > 0:
                    detected_juncs.add(coord)
                    gene = coord_to_gene.get(coord, "unknown")
                    detected_genes.add(gene)
            sample_junc_counts.append(len(detected_juncs))
            sample_gene_counts.append(len(detected_genes))

        jc = np.array(sample_junc_counts, dtype=float)
        gc = np.array(sample_gene_counts, dtype=float)
        rows.append({
            "timepoint": tp,
            "n_samples": len(idxs),
            "mean_up_junctions": round(float(np.mean(jc)), 1),
            "std_up_junctions": round(float(np.std(jc, ddof=1)), 1) if len(jc) > 1 else 0.0,
            "mean_up_genes": round(float(np.mean(gc)), 1),
            "std_up_genes": round(float(np.std(gc, ddof=1)), 1) if len(gc) > 1 else 0.0,
        })

    df = pd.DataFrame(rows)

    # ---- 5. Write TSV ----
    tsv_path = outdir / "up_transcript_counts.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[UP-COUNTS] Wrote {tsv_path}")

    # ---- 6. Plot ----
    png_path = outdir / "up_transcript_counts.png"
    _plot_up_transcript_counts(df, png_path)

    for _, r in df.iterrows():
        print(f"[UP-COUNTS]   {r['timepoint']}: {r['mean_up_junctions']:.0f} junctions "
              f"from {r['mean_up_genes']:.0f} genes (n={r['n_samples']})")

    return {
        "tsv": "up_transcript_counts.tsv",
        "png": "up_transcript_counts.png",
        "timepoints": [r["timepoint"] for r in rows],
        "mean_up_junctions": [r["mean_up_junctions"] for r in rows],
        "mean_up_genes": [r["mean_up_genes"] for r in rows],
    }


def _plot_up_transcript_counts(df: pd.DataFrame, out_path: Path) -> None:
    """Grouped bar chart: unique UP junctions and genes per developmental timepoint."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(df)
    x = np.arange(n)
    bar_w = 0.35
    fig_width = max(6, n * 1.4 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, 5), dpi=150)

    bars_j = ax.bar(
        x - bar_w / 2, df["mean_up_junctions"], bar_w,
        yerr=df["std_up_junctions"], capsize=3,
        color="#4e79a7", edgecolor="white", linewidth=0.5, alpha=0.88,
        label="UP junctions", error_kw=dict(lw=1.0, capthick=0.8),
    )
    bars_g = ax.bar(
        x + bar_w / 2, df["mean_up_genes"], bar_w,
        yerr=df["std_up_genes"], capsize=3,
        color="#e15759", edgecolor="white", linewidth=0.5, alpha=0.88,
        label="Unique genes", error_kw=dict(lw=1.0, capthick=0.8),
    )

    for bar_set in (bars_j, bars_g):
        for bar in bar_set:
            h = bar.get_height()
            if h > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, h,
                    f"{h:.0f}", ha="center", va="bottom",
                    fontsize=7, color="#333333",
                )

    ax.set_xticks(x)
    rotation = 0 if n <= 4 else 35
    ha = "center" if rotation == 0 else "right"
    ax.set_xticklabels(df["timepoint"], rotation=rotation, ha=ha, fontsize=9)
    ax.set_ylabel("Count (mean per sample)", fontsize=9.5)
    ax.set_title("UP transcripts detected per developmental stage\n",
                 fontsize=12, fontweight="bold", pad=4)
    ax.text(0.5, 1.02,
            f"{n} timepoints  \u00b7  bars = mean \u00b1 s.d. across samples",
            transform=ax.transAxes, ha="center", fontsize=8, color="#888888")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.85)

    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.25, zorder=0)
    ax.tick_params(axis="both", which="major", labelsize=8.5)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[UP-COUNTS] Plot: {out_path}")


# ----------------------------------------------------------
# Poison-exon / developmental-UP analysis
# ----------------------------------------------------------

def compute_poison_exon_analysis(
    junction_counts_path: Path,
    cluster_ratios_path: Path,
    classification_dir: Path,
    lc2_prefix: str,
    junction_filelist_path: Path,
    outdir: Path,
    samples_tsv_path: Optional[Path] = None,
) -> Optional[Dict]:
    """Identify poison-exon candidates by correlating UP junction PSI with
    developmental time.  Generates a ranked TSV, heatmap, correlation scatter,
    and NMD-feature distribution plots.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError:
        print("[WARN] matplotlib not available — skipping poison-exon analysis")
        return None

    # --- paths for LC2 auxiliary files ---
    classif_path = classification_dir / f"{lc2_prefix}_junction_classifications.txt"
    long_exon_path = classification_dir / f"{lc2_prefix}_long_exon_distances.txt"
    nuc_rule_path = classification_dir / f"{lc2_prefix}_nuc_rule_distances.txt"
    exon_stats_path = classification_dir / f"{lc2_prefix}_exon_stats.txt"

    for fp, label in [
        (cluster_ratios_path, "cluster_ratios"),
        (classif_path, "junction_classifications"),
    ]:
        if not fp.exists():
            print(f"[WARN] {label} not found at {fp} — skipping poison-exon analysis")
            return None

    # ---- 1. Parse cluster_ratios.gz for UP junctions ----
    with gzip.open(str(cluster_ratios_path), "rt") as f:
        header = f.readline().strip().split()
        sample_cols = header[1:]
        n_samples = len(sample_cols)

        up_junctions: Dict[str, np.ndarray] = {}      # junc_id -> PSI array
        up_denominators: Dict[str, np.ndarray] = {}
        for line in f:
            fields = line.strip().split()
            junc_id = fields[0]
            if not junc_id.endswith(":UP"):
                continue
            psi = np.full(n_samples, np.nan)
            denom = np.zeros(n_samples)
            for j in range(n_samples):
                token = fields[j + 1] if j + 1 < len(fields) else "0/0"
                if "/" in token:
                    num, den = token.split("/")
                    num, den = float(num), float(den)
                    denom[j] = den
                    psi[j] = num / den if den > 0 else np.nan
                else:
                    psi[j] = np.nan
            up_junctions[junc_id] = psi
            up_denominators[junc_id] = denom

    if not up_junctions:
        print("[WARN] No UP junctions in cluster_ratios — skipping poison-exon analysis")
        return None

    # ---- 2. Build coord -> gene name mapping from classifications ----
    coord_to_gene: Dict[str, str] = {}
    if classif_path.exists():
        with open(classif_path) as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    coord_to_gene[parts[1]] = parts[0]

    def _junc_id_to_coord(jid: str) -> str:
        """chr1:2425958:2525850:clu_5_+:UP -> chr1:2425958-2525850"""
        toks = jid.split(":")
        return f"{toks[0]}:{toks[1]}-{toks[2]}"

    def _junc_id_to_cluster(jid: str) -> str:
        toks = jid.split(":")
        return toks[3] if len(toks) > 3 else ""

    # ---- 3. Resolve timepoint labels ----
    tissues = _resolve_tissue_labels(sample_cols, junction_filelist_path, samples_tsv_path)
    unique_tissues = sorted(set(tissues), key=lambda t: _timepoint_to_hours(t))
    is_dev = _is_developmental_stage_labels(unique_tissues)
    if not is_dev:
        print("[WARN] Labels are not developmental stages — skipping poison-exon analysis")
        return None

    tissue_hours = np.array([_timepoint_to_hours(t) for t in unique_tissues])
    tissue_idx: Dict[str, List[int]] = defaultdict(list)
    for i, t in enumerate(tissues):
        tissue_idx[t].append(i)

    # ---- 4. Per-junction per-timepoint mean PSI + Spearman ----
    rows = []
    for jid, psi_arr in up_junctions.items():
        coord = _junc_id_to_coord(jid)
        gene = coord_to_gene.get(coord, "unknown")
        cluster = _junc_id_to_cluster(jid)

        tp_means = []
        for tp in unique_tissues:
            idxs = tissue_idx[tp]
            vals = psi_arr[idxs]
            valid = vals[~np.isnan(vals)]
            tp_means.append(float(np.mean(valid)) if len(valid) > 0 else np.nan)

        tp_means_arr = np.array(tp_means)
        finite = np.isfinite(tp_means_arr) & np.isfinite(tissue_hours)
        rho, p = (None, None)
        if finite.sum() >= 3:
            rho, p = _spearman_rho(tissue_hours[finite], tp_means_arr[finite])

        total_reads = float(up_denominators[jid].sum())
        row = {
            "gene": gene,
            "junction_coord": coord,
            "cluster_id": cluster,
            "total_reads": round(total_reads),
        }
        for k, tp in enumerate(unique_tissues):
            row[f"psi_{tp}"] = round(tp_means[k], 6) if not np.isnan(tp_means[k]) else None
        row["spearman_rho"] = round(rho, 4) if rho is not None else None
        row["spearman_p"] = round(p, 4) if p is not None else None
        rows.append(row)

    # ---- 5. Merge NMD features ----
    long_exon_map: Dict[str, Dict] = {}
    if long_exon_path.exists():
        with open(long_exon_path) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    try:
                        long_exon_map[parts[1]] = {
                            "ptc_distance": float(parts[2]),
                            "exon_length": float(parts[3]),
                        }
                    except ValueError:
                        pass

    nuc_rule_map: Dict[str, float] = {}
    if nuc_rule_path.exists():
        with open(nuc_rule_path) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    try:
                        nuc_rule_map[parts[1]] = float(parts[2])
                    except ValueError:
                        pass

    exon_stats_map: Dict[str, Dict] = {}
    if exon_stats_path.exists():
        with open(exon_stats_path) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    try:
                        exon_stats_map[parts[1]] = {
                            "exons_before": float(parts[2]),
                            "exons_after": float(parts[3]),
                        }
                    except ValueError:
                        pass

    for row in rows:
        coord = row["junction_coord"]
        le = long_exon_map.get(coord, {})
        row["ptc_distance"] = le.get("ptc_distance")
        row["exon_length"] = le.get("exon_length")
        row["ejc_distance"] = nuc_rule_map.get(coord)
        es = exon_stats_map.get(coord, {})
        row["exons_before"] = es.get("exons_before")
        row["exons_after"] = es.get("exons_after")
        row["nmd_eligible"] = (
            (row["ejc_distance"] is not None and row["ejc_distance"] > 50)
            or (row["exon_length"] is not None and row["exon_length"] > 50)
        )

    # Sort by |spearman_rho| descending (best candidates first)
    rows.sort(key=lambda r: abs(r["spearman_rho"]) if r["spearman_rho"] is not None else -1, reverse=True)

    # ---- 6. Write TSV ----
    tsv_path = outdir / "poison_exon_candidates.tsv"
    df = pd.DataFrame(rows)
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[POISON] Wrote {len(rows)} candidates to {tsv_path}")

    # ---- 7. Positive-correlation candidates ----
    pos_rows = [
        r for r in rows
        if r.get("spearman_rho") is not None and r["spearman_rho"] > 0
        and r.get("spearman_p") is not None and r["spearman_p"] < 0.05
    ]
    pos_rows.sort(key=lambda r: r["spearman_rho"], reverse=True)
    n_pos = len(pos_rows)

    if pos_rows:
        pos_tsv_path = outdir / "positive_correlation_candidates.tsv"
        pd.DataFrame(pos_rows).to_csv(pos_tsv_path, sep="\t", index=False)
        print(f"[POISON] {n_pos} candidates with POSITIVE PSI-time correlation (p<0.05):")
        for pr in pos_rows:
            nmd_tag = " [NMD-eligible]" if pr.get("nmd_eligible") else ""
            print(f"[POISON]   {pr['gene']:20s}  rho={pr['spearman_rho']:+.4f}  "
                  f"p={pr['spearman_p']:.4f}{nmd_tag}")
        print(f"[POISON] Wrote {pos_tsv_path}")
    else:
        print("[POISON] No candidates with significant positive PSI-time correlation")

    # ---- 8. Write JSON ----
    n_nmd = sum(1 for r in rows if r.get("nmd_eligible"))
    n_sig = sum(1 for r in rows if r.get("spearman_p") is not None and r["spearman_p"] < 0.05)
    json_payload = {
        "n_candidates": len(rows),
        "n_nmd_eligible": n_nmd,
        "n_significant": n_sig,
        "n_positive_correlation": n_pos,
        "timepoints": unique_tissues,
        "positive_correlation": pos_rows,
        "candidates": rows,
    }
    json_path = outdir / "poison_exon_candidates.json"
    with open(json_path, "w") as fh:
        json.dump(json_payload, fh, indent=2)

    # ---- 9. Plots ----
    _plot_poison_exon_heatmap(rows, unique_tissues, outdir / "poison_exon_heatmap.png")
    _plot_up_psi_correlation(rows, outdir / "up_psi_vs_time.png")
    _plot_nmd_features(rows, long_exon_map, nuc_rule_map, outdir / "nmd_features.png")

    print(f"[POISON] {len(rows)} UP junctions, {n_nmd} NMD-eligible, "
          f"{n_sig} with p<0.05, {n_pos} positive correlation")

    return {
        "tsv": "poison_exon_candidates.tsv",
        "json": "poison_exon_candidates.json",
        "positive_tsv": "positive_correlation_candidates.tsv" if pos_rows else None,
        "heatmap_png": "poison_exon_heatmap.png",
        "correlation_png": "up_psi_vs_time.png",
        "nmd_png": "nmd_features.png",
        "n_candidates": len(rows),
        "n_nmd_eligible": n_nmd,
        "n_significant": n_sig,
        "n_positive_correlation": n_pos,
    }


def _plot_poison_exon_heatmap(rows, timepoints, out_path, max_genes=30):
    """Heatmap of top UP junction PSI values across developmental timepoints."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    top = [r for r in rows if r.get("spearman_rho") is not None][:max_genes]
    if not top:
        return

    psi_cols = [f"psi_{tp}" for tp in timepoints]
    matrix = []
    ylabels = []
    for r in top:
        vals = [r.get(c, np.nan) for c in psi_cols]
        vals = [v if v is not None else np.nan for v in vals]
        matrix.append(vals)
        sig = "*" if r.get("spearman_p") is not None and r["spearman_p"] < 0.05 else ""
        nmd = " [NMD]" if r.get("nmd_eligible") else ""
        ylabels.append(f"{r['gene']} ({r['junction_coord']}){sig}{nmd}")

    matrix = np.array(matrix, dtype=float)

    n_rows = len(top)
    fig_h = max(4, n_rows * 0.35 + 1.5)
    fig, ax = plt.subplots(figsize=(max(6, len(timepoints) * 1.2 + 3), fig_h), dpi=150)

    cmap = LinearSegmentedColormap.from_list("psi", ["#f0f0f0", "#4e79a7", "#e15759"])
    vmax = max(np.nanmax(matrix), 0.01) if matrix.size > 0 else 0.1
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")

    ax.set_xticks(range(len(timepoints)))
    ax.set_xticklabels(timepoints, fontsize=9)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.set_xlabel("Developmental stage", fontsize=10)
    ax.set_title("Poison-exon candidates: UP junction PSI across timepoints\n"
                 "(ranked by |Spearman \u03c1|; * = p<0.05)", fontsize=11, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Mean PSI (junction usage within cluster)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    for i in range(n_rows):
        for j in range(len(timepoints)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                        fontsize=6, color="white" if v > vmax * 0.6 else "#333333")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[POISON] Heatmap: {out_path}")


def _plot_up_psi_correlation(rows, out_path):
    """Lollipop plot of per-gene Spearman rho values, color-coded by direction and significance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    valid = [r for r in rows if r.get("spearman_rho") is not None]
    if not valid:
        return

    valid.sort(key=lambda r: r["spearman_rho"])
    genes = [r["gene"] for r in valid]
    rhos = [r["spearman_rho"] for r in valid]
    sig = [r.get("spearman_p", 1.0) is not None and r.get("spearman_p", 1.0) < 0.05 for r in valid]
    nmd = [r.get("nmd_eligible", False) for r in valid]
    positive = [rho > 0 for rho in rhos]

    n = len(valid)
    fig_h = max(4, n * 0.3 + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_h), dpi=150)

    y = np.arange(n)
    # Negative rho = blues, positive rho = reds/oranges
    # Significance and NMD status control shade/intensity
    colors = []
    for s, nm, pos in zip(sig, nmd, positive):
        if pos:
            if s and nm:
                colors.append("#c72e2e")     # dark red: sig + NMD + positive
            elif s:
                colors.append("#e15759")     # red: sig + positive
            else:
                colors.append("#f5a6a7")     # light red: not sig + positive
        else:
            if s and nm:
                colors.append("#2d5f8a")     # dark blue: sig + NMD + negative
            elif s:
                colors.append("#4e79a7")     # blue: sig + negative
            else:
                colors.append("#a8c4db")     # light blue: not sig + negative

    ax.barh(y, rhos, height=0.6, color=colors, edgecolor="white", linewidth=0.4, alpha=0.88)
    ax.axvline(0, color="#888888", linewidth=0.7, zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xlabel("Spearman \u03c1 (PSI vs developmental time)", fontsize=9.5)
    ax.set_title("UP junction PSI correlation with developmental time\n", fontsize=11, fontweight="bold")

    legend_handles = [
        Patch(facecolor="#2d5f8a", label="Sig + NMD (neg)"),
        Patch(facecolor="#4e79a7", label="Significant (neg)"),
        Patch(facecolor="#a8c4db", label="Not sig (neg)"),
        Patch(facecolor="#c72e2e", label="Sig + NMD (pos)"),
        Patch(facecolor="#e15759", label="Significant (pos)"),
        Patch(facecolor="#f5a6a7", label="Not sig (pos)"),
    ]
    ax.legend(handles=legend_handles, fontsize=6.5, loc="lower right",
              framealpha=0.85, ncol=2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle="--", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[POISON] Correlation plot: {out_path}")


def _plot_nmd_features(rows, long_exon_map, nuc_rule_map, out_path):
    """Two-panel histogram of PTC distances and EJC distances for UP junctions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ptc_dists = [r["ptc_distance"] for r in rows if r.get("ptc_distance") is not None]
    ejc_dists = [r["ejc_distance"] for r in rows if r.get("ejc_distance") is not None]

    # Also include ALL entries from the NMD files (not just UP junctions)
    all_exon_lengths = [v["exon_length"] for v in long_exon_map.values()]
    all_ejc = list(nuc_rule_map.values())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=150)

    # Panel A: Exon lengths from long_exon_distances
    ax = axes[0]
    if all_exon_lengths:
        ax.hist(all_exon_lengths, bins=30, color="#4e79a7", alpha=0.7, edgecolor="white", label="All UP junctions")
        if ptc_dists:
            ax.hist(ptc_dists, bins=30, color="#e15759", alpha=0.6, edgecolor="white", label="Clustered UP")
        ax.axvline(50, color="#333333", linestyle="--", linewidth=1.0, label="50 nt threshold")
        ax.set_xlabel("PTC distance from exon start (nt)", fontsize=9)
        ax.set_ylabel("Count", fontsize=9)
        ax.set_title("Premature termination codon position", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No PTC distance data", transform=ax.transAxes, ha="center", fontsize=10, color="#aaa")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: EJC distances from nuc_rule_distances
    ax = axes[1]
    if all_ejc:
        ax.hist(all_ejc, bins=30, color="#59a14f", alpha=0.7, edgecolor="white", label="All UP junctions")
        if ejc_dists:
            ax.hist(ejc_dists, bins=30, color="#e15759", alpha=0.6, edgecolor="white", label="Clustered UP")
        ax.axvline(50, color="#333333", linestyle="--", linewidth=1.0, label="50 nt rule")
        ax.set_xlabel("Distance from PTC to last exon-exon junction (nt)", fontsize=9)
        ax.set_ylabel("Count", fontsize=9)
        ax.set_title("EJC distance (NMD eligibility)", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "No EJC distance data", transform=ax.transAxes, ha="center", fontsize=10, color="#aaa")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("NMD feature distributions for unproductive junctions", fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[POISON] NMD features: {out_path}")


def _expand_star_files(patterns: List[str]) -> List[str]:
    """Expand globs and deduplicate while preserving order."""
    paths: List[str] = []
    for pat in patterns:
        expanded = sorted(glob(pat))
        if not expanded and ("*" not in pat and "?" not in pat and "[" not in pat):
            expanded = [pat]
        paths.extend(expanded)
    seen = set()
    uniq = []
    for p in paths:
        if p not in seen:
            uniq.append(p); seen.add(p)
    return uniq

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Python-only LeafCutter + LeafCutter2 pipeline (with sane defaults).")

    # Repos / refs
    p.add_argument("--leafcutter_repo", default="tools/leafcutter",
                   help="Path to LeafCutter repo (default: tools/leafcutter)")
    p.add_argument("--leafcutter2_repo", default="tools/leafcutter2",
                   help="Path to LeafCutter2 repo (default: tools/leafcutter2)")
    p.add_argument("--genome_fasta", default="refs/GRCh38.fa",
                   help="Reference genome FASTA (default: refs/GRCh38.fa)")
    p.add_argument("--gencode_gtf", default="refs/gencode.v46.annotation.gtf",
                   help="GENCODE GTF with CDS/start/stop (default: refs/gencode.v46.annotation.gtf)")

    # Samples
    p.add_argument("--samples_tsv", default="out/samples.tsv",
                   help="Optional sample sheet for downstream analysis; not required for core LC2 run.")

    # Junction inputs (STAR SJ files) — default to Brain_Cortex + Liver from recount3 conversion
    default_sj = (
        sorted(glob("star_sj/Brain_Cortex/*.SJ.out.tab")) +
        sorted(glob("star_sj/Liver/*.SJ.out.tab"))
    )
    p.add_argument("--star_sj", nargs="+", default=default_sj,
                   help="One or more STAR SJ.out.tab files (supports globs)")

    # Alternative BAM route (kept optional)
    p.add_argument("--bams", nargs="*", type=Path, default=None,
                   help="Alternative: BAMs to extract junctions via regtools (if not providing --star_sj)")

    # Pre-made BED files (e.g. from gtex_gct_to_bed.py or recount3_to_bed.py)
    p.add_argument("--junction_beds", nargs="+", default=None,
                   help="Pre-made junction BED files — skip SJ/BAM conversion and go straight to clustering")

    # Run options
    p.add_argument("--workdir", default=".",
                   help="Working/output directory (default: .)")
    p.add_argument("--prefix", default="gtex_demo",
                   help="Run prefix used for output file names (default: gtex_demo)")

    # LeafCutter / LC2 knobs
    p.add_argument("--min_reads", type=int, default=50,
                   help="Minimum reads per intron for clustering (default: 50)")
    p.add_argument("--max_intron_len", type=int, default=500000,
                   help="Maximum intron length for clustering (default: 500000)")
    p.add_argument("--nmd_distance", type=int, default=55,
                   help="Reserved for downstream analysis extensions; unused in core LC2 run.")
    p.add_argument("--gene_type_tag", default="gene_biotype",
                   help="GTF attribute tag for gene type (default: gene_biotype)")
    p.add_argument("--transcript_type_tag", default="transcript_biotype",
                   help="GTF attribute tag for transcript type (default: transcript_biotype)")
    p.add_argument("--gene_name_tag", default="gene_id",
                   help="GTF attribute tag for gene label (default: gene_id)")
    p.add_argument("--transcript_name_tag", default="transcript_id",
                   help="GTF attribute tag for transcript label (default: transcript_id)")

    return p

# ---------------------------
# Main
# ---------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()

    # Expand any globs passed to --star_sj (or used in defaults)
    sj_files = _expand_star_files(args.star_sj)

    # --- sanity checks ---
    problems = []
    for path in [args.leafcutter_repo, args.leafcutter2_repo, args.genome_fasta, args.gencode_gtf]:
        if not os.path.exists(path):
            problems.append(f"Missing: {path}")
    if not sj_files and not args.bams and not args.junction_beds:
        problems.append("No junction inputs found. Provide --star_sj, --bams, or --junction_beds.")
    if problems:
        msg = "Input validation failed:\n  - " + "\n  - ".join(problems)
        raise SystemExit(msg)

    # Prepare dirs
    wd = ensure_dir(Path(args.workdir))
    jdir = ensure_dir(wd / "junctions")
    cdir = ensure_dir(wd / "clusters")
    lc2dir = ensure_dir(wd / "lc2")
    outdir = ensure_dir(wd / "out")

    # Step 1: Junction files
    if args.junction_beds:
        bed_paths = [Path(p) for p in args.junction_beds]
        missing_beds = [str(p) for p in bed_paths if not p.exists()]
        if missing_beds:
            raise SystemExit(f"Missing BED files:\n  " + "\n  ".join(missing_beds))
        print(f"Using {len(bed_paths)} pre-made junction BED files")
    elif sj_files:
        bed_paths = convert_star_sj_to_bed([Path(p) for p in sj_files], jdir)
    else:
        check_exe("regtools")
        bed_paths = extract_junctions_from_bams([Path(p) for p in args.bams], jdir)

    # Step 2: Cluster (LeafCutter — optional)
    # LC2 does its own clustering internally, so this step is optional.
    # If the LC1 clustering script is available we run it for the perind
    # outputs; otherwise we just write the junction filelist and let LC2
    # handle clustering.
    prefix_path = cdir / args.prefix
    junction_filelist = write_junction_filelist(
        bed_paths, cdir / (args.prefix + "_junction_files.txt"),
    )
    perind = perind_num = None
    cluster_result = leafcutter_cluster(
        Path(args.leafcutter_repo),
        bed_paths,
        prefix_path,
        min_reads=args.min_reads,
        max_intron_len=args.max_intron_len,
        leafcutter2_repo=Path(args.leafcutter2_repo),
    )
    if cluster_result is not None:
        junction_filelist, perind, perind_num = cluster_result
    else:
        print("[INFO] No clustering script found — LeafCutter2 will cluster internally")

    # Step 3: Classify (LeafCutter2 — does clustering + classification)
    gtf_path = normalize_gtf_for_lc2(Path(args.gencode_gtf), wd)
    lc2_prefix = args.prefix + "_lc2"
    lc2_outputs = leafcutter2_classify(
        Path(args.leafcutter2_repo),
        Path(args.genome_fasta),
        gtf_path,
        junction_filelist,
        lc2dir,
        lc2_prefix,
        min_cluster_reads=args.min_reads,
        max_intron_len=args.max_intron_len,
        gene_type_tag=args.gene_type_tag,
        transcript_type_tag=args.transcript_type_tag,
        gene_name_tag=args.gene_name_tag,
        transcript_name_tag=args.transcript_name_tag,
    )

    # Step 4: Per-tissue unproductive junction read analysis
    samples_tsv = Path(args.samples_tsv) if args.samples_tsv else None
    tissue_metrics = None
    try:
        tissue_metrics = compute_unproductive_by_tissue(
            junction_counts_path=Path(lc2_outputs["junction_counts"]),
            junction_filelist_path=junction_filelist,
            outdir=outdir,
            samples_tsv_path=samples_tsv,
        )
    except Exception as exc:
        print(f"[WARN] Tissue unproductive analysis failed: {exc}")

    # Step 4b: UP transcript / gene counts per timepoint
    up_count_metrics = None
    try:
        up_count_metrics = compute_up_transcript_counts(
            junction_counts_path=Path(lc2_outputs["junction_counts"]),
            classification_dir=Path(lc2_outputs["run_dir"]) / "clustering",
            lc2_prefix=lc2_prefix,
            junction_filelist_path=junction_filelist,
            outdir=outdir,
            samples_tsv_path=samples_tsv,
        )
    except Exception as exc:
        print(f"[WARN] UP transcript count analysis failed: {exc}")

    # Step 5: Poison-exon / developmental-UP analysis
    poison_metrics = None
    try:
        poison_metrics = compute_poison_exon_analysis(
            junction_counts_path=Path(lc2_outputs["junction_counts"]),
            cluster_ratios_path=Path(lc2_outputs["cluster_ratios"]),
            classification_dir=Path(lc2_outputs["run_dir"]) / "clustering",
            lc2_prefix=lc2_prefix,
            junction_filelist_path=junction_filelist,
            outdir=outdir,
            samples_tsv_path=samples_tsv,
        )
    except Exception as exc:
        import traceback
        print(f"[WARN] Poison-exon analysis failed: {exc}")
        traceback.print_exc()

    # Save JSON summary
    summary = {
        "n_junction_inputs": len(bed_paths),
        "junction_file_list": str(junction_filelist),
        "lc2_outputs": lc2_outputs,
        "cluster_files": {},
        "metrics": {},
        "notes": "Core pipeline run complete. Use long_exon/nuc_rule outputs for NMD-focused downstream analysis.",
    }
    if perind:
        summary["cluster_files"]["perind_counts"] = str(perind)
    if perind_num:
        summary["cluster_files"]["perind_numers"] = str(perind_num)
    if tissue_metrics:
        summary["metrics"]["unproductive_by_tissue"] = tissue_metrics
    if up_count_metrics:
        summary["metrics"]["up_transcript_counts"] = up_count_metrics
    if poison_metrics:
        summary["metrics"]["poison_exon_analysis"] = poison_metrics

    with open(outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\nDone.", flush=True)
    print(f"- LC2 cluster ratios: {lc2_outputs['cluster_ratios']}")
    print(f"- LC2 junction counts: {lc2_outputs['junction_counts']}")
    print(f"- Summary: {outdir / 'summary.json'}", flush=True)

if __name__ == "__main__":
    main()
