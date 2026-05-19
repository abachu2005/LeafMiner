#!/usr/bin/env python3
"""
Quest-oriented zebrafish workflow with explicit stages:

- qc:             download FASTQs + fastp trim/QC + FastQC + MultiQC reports (serial)
- qc_array:       single-sample QC (one per SLURM_ARRAY_TASK_ID)
- qc_finalize:    aggregate per-sample QC results + MultiQC
- stage1:         STAR alignment (serial, all samples in one job)
- stage1_setup:   build STAR index (run once before array)
- stage1_array:   single-sample STAR alignment (one per SLURM_ARRAY_TASK_ID)
- stage1_finalize: aggregate SJ files into stage1 metadata
- stage2:         run LeafCutter2 analysis from previously persisted stage1 artifacts
- stage3:         POISEN PE differential inclusion analysis
"""


import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


def run(cmd: List[str]) -> None:
    print(">>", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def check_exe(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def safe_name(name: str) -> str:
    out = name.replace(" - ", "_").replace(" ", "_")
    for ch in "()/'\"":
        out = out.replace(ch, "")
    return out


def parse_manifest(path: Path, max_runs: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            run_id = (row.get("run_id") or "").strip()
            if not run_id:
                continue
            rows.append(
                {
                    "run_id": run_id,
                    "sample": (row.get("sample") or run_id).strip() or run_id,
                    "condition": (row.get("condition") or "zebrafish").strip() or "zebrafish",
                    "library_layout": (row.get("library_layout") or "SINGLE").strip().upper() or "SINGLE",
                    "fastq_ftp": (row.get("fastq_ftp") or "").strip(),
                }
            )
    if max_runs > 0:
        rows = rows[:max_runs]
    return rows


def download_fastqs(run_row: Dict[str, str], out_dir: Path) -> List[Path]:
    ftp = run_row.get("fastq_ftp", "")
    if not ftp:
        raise RuntimeError(f"Run {run_row['run_id']} has no fastq_ftp entry in manifest.")
    paths = [x.strip() for x in ftp.split(";") if x.strip()]
    if not paths:
        raise RuntimeError(f"Run {run_row['run_id']} has empty fastq_ftp paths.")

    local_paths: List[Path] = []
    for p in paths:
        url = p if p.startswith("http://") or p.startswith("https://") else f"https://{p}"
        fname = Path(p).name
        dest = out_dir / fname
        if dest.exists() and dest.stat().st_size > 0:
            print(f"[ZEBRAFISH] Reusing existing {fname}", flush=True)
        else:
            tmp = out_dir / f".{fname}.partial"
            if tmp.exists():
                tmp.unlink()
            run([
                "wget",
                "--tries=3",
                "--timeout=300",
                "--waitretry=30",
                "--progress=dot:mega",
                "-O", str(tmp),
                url,
            ])
            if not tmp.exists() or tmp.stat().st_size == 0:
                raise RuntimeError(f"Download failed — {fname} is empty or missing")
            tmp.rename(dest)
        local_paths.append(dest)
    return local_paths


def run_fastp(
    run_id: str,
    fastqs: List[Path],
    qc_dir: Path,
    threads: int = 4,
    min_length: int = 36,
    min_surviving_pct: float = 50.0,
) -> tuple:
    """Run fastp QC/trimming on raw FASTQs. Returns (trimmed_paths, qc_stats)."""
    qc_dir.mkdir(parents=True, exist_ok=True)
    report_dir = qc_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    trimmed_dir = qc_dir / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    is_paired = len(fastqs) >= 2
    trimmed: List[Path] = []
    json_report = report_dir / f"{run_id}.fastp.json"
    html_report = report_dir / f"{run_id}.fastp.html"

    cmd = [
        "fastp",
        "--thread", str(threads),
        "--length_required", str(min_length),
        "--detect_adapter_for_pe",
        "--json", str(json_report),
        "--html", str(html_report),
    ]

    if is_paired:
        out_r1 = trimmed_dir / f"{run_id}_1.trimmed.fastq.gz"
        out_r2 = trimmed_dir / f"{run_id}_2.trimmed.fastq.gz"
        cmd.extend(["-i", str(fastqs[0]), "-I", str(fastqs[1])])
        cmd.extend(["-o", str(out_r1), "-O", str(out_r2)])
        trimmed = [out_r1, out_r2]
    else:
        out_r1 = trimmed_dir / f"{run_id}.trimmed.fastq.gz"
        cmd.extend(["-i", str(fastqs[0])])
        cmd.extend(["-o", str(out_r1)])
        trimmed = [out_r1]

    print(f"[ZEBRAFISH][QC] Running fastp on {run_id} ({'PE' if is_paired else 'SE'}) ...", flush=True)
    run(cmd)

    stats: Dict = {}
    if json_report.exists():
        with open(json_report, "r") as fh:
            report = json.load(fh)
        summary = report.get("summary", {})
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})
        total_reads = before.get("total_reads", 0)
        passed_reads = after.get("total_reads", 0)
        pct_pass = (passed_reads / total_reads * 100) if total_reads > 0 else 0
        adapter_pct = report.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
        if total_reads > 0:
            adapter_pct = adapter_pct / total_reads * 100

        stats = {
            "total_reads_before": total_reads,
            "total_reads_after": passed_reads,
            "pct_passing": round(pct_pass, 2),
            "adapter_trimmed_pct": round(adapter_pct, 2),
            "q30_rate_before": before.get("q30_rate", 0),
            "q30_rate_after": after.get("q30_rate", 0),
        }
        print(
            f"[ZEBRAFISH][QC] {run_id}: {total_reads:,} reads -> {passed_reads:,} passed "
            f"({pct_pass:.1f}%), adapters trimmed: {adapter_pct:.1f}%",
            flush=True,
        )

        if pct_pass < min_surviving_pct:
            print(
                f"[ZEBRAFISH][QC] WARNING: {run_id} has only {pct_pass:.1f}% reads surviving "
                f"(threshold: {min_surviving_pct}%). Skipping this sample.",
                flush=True,
            )
            return [], stats

    return trimmed, stats


def run_fastqc(fastqs: List[Path], outdir: Path, threads: int = 4) -> None:
    """Run FastQC on a list of FASTQ files. Skips gracefully if fastqc is not available."""
    if shutil.which("fastqc") is None:
        print("[ZEBRAFISH][QC] FastQC not found on PATH — skipping.", flush=True)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["fastqc", "--outdir", str(outdir), "--threads", str(threads), "--quiet"]
    cmd.extend(str(p) for p in fastqs)
    print(f"[ZEBRAFISH][QC] Running FastQC on {len(fastqs)} file(s) ...", flush=True)
    try:
        run(cmd)
    except subprocess.CalledProcessError:
        print("[ZEBRAFISH][QC] FastQC failed — continuing without it.", flush=True)


def run_multiqc(qc_dir: Path, outdir: Path) -> None:
    """Run MultiQC on all reports in qc_dir. Skips gracefully if multiqc is not available."""
    if shutil.which("multiqc") is None:
        print("[ZEBRAFISH][QC] MultiQC not found on PATH — skipping.", flush=True)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "multiqc",
        str(qc_dir),
        "--outdir", str(outdir),
        "--filename", "multiqc_report.html",
        "--force",
        "--quiet",
    ]
    print("[ZEBRAFISH][QC] Running MultiQC ...", flush=True)
    try:
        run(cmd)
    except subprocess.CalledProcessError:
        print("[ZEBRAFISH][QC] MultiQC failed — continuing without it.", flush=True)


def write_qc_metadata(
    *,
    out_dir: Path,
    manifest: Path,
    all_runs: List[Dict[str, str]],
    passed_runs: List[Dict[str, str]],
    qc_stats: Dict[str, Dict],
    qc_dir: Path,
) -> None:
    """Write QC stage output metadata to out/qc_report.json and out/summary.json."""
    out_dir.mkdir(parents=True, exist_ok=True)

    trimmed_dir = qc_dir / "trimmed"
    trimmed_map: Dict[str, List[str]] = {}
    for row in passed_runs:
        rid = row["run_id"]
        candidates = sorted(trimmed_dir.glob(f"{rid}*trimmed*"))
        trimmed_map[rid] = [str(p) for p in candidates if p.exists()]

    stages = sorted({r.get("condition", "Unknown_Stage") for r in passed_runs})
    payload = {
        "stage": "qc_report",
        "manifest": str(manifest),
        "n_total": len(all_runs),
        "n_passed": len(passed_runs),
        "n_failed": len(all_runs) - len(passed_runs),
        "stages": stages,
        "qc_summary": qc_stats,
        "trimmed_fastqs": trimmed_map,
        "passed_runs": [
            {"run_id": r["run_id"], "condition": r["condition"]} for r in passed_runs
        ],
    }
    with open(out_dir / "qc_report.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(
            {
                "stage": "qc_report",
                "n_total": len(all_runs),
                "n_passed": len(passed_runs),
                "n_failed": len(all_runs) - len(passed_runs),
                "stages": stages,
                "qc_summary": qc_stats,
                "artifacts": {"qc_report": "qc_report.json"},
                "notes": "QC complete. Review reports then submit alignment (stage1).",
            },
            fh,
            indent=2,
        )


def run_qc(args: argparse.Namespace) -> None:
    """Stage QC: download FASTQs, run fastp + optional FastQC, aggregate with MultiQC."""
    workdir = Path(args.workdir).resolve()
    manifest = Path(args.manifest).resolve()
    fastq_dir = workdir / "inputs" / "fastq"
    qc_dir = workdir / "qc"
    outdir = workdir / "out"

    if not manifest.exists():
        raise SystemExit(f"Manifest not found: {manifest}")

    fastq_dir.mkdir(parents=True, exist_ok=True)

    runs = parse_manifest(manifest, args.max_runs)
    if not runs:
        raise SystemExit("No runs found in zebrafish manifest.")
    print(f"[ZEBRAFISH][qc] Runs to process: {len(runs)}", flush=True)

    skip_fastqc = getattr(args, "skip_fastqc", False)
    fastqc_dir = qc_dir / "fastqc"
    qc_stats_all: Dict[str, Dict] = {}
    passed_runs: List[Dict[str, str]] = []

    for row in runs:
        run_id = row["run_id"]
        print(f"[ZEBRAFISH][qc] Downloading FASTQ for {run_id}", flush=True)
        fastqs = download_fastqs(row, fastq_dir)
        if len(fastqs) > 2:
            fastqs = fastqs[:2]

        trimmed, qc_stats = run_fastp(
            run_id, fastqs, qc_dir,
            threads=min(4, args.star_threads),
        )
        qc_stats_all[run_id] = qc_stats

        if not trimmed:
            print(f"[ZEBRAFISH][qc] {run_id} failed QC — skipping", flush=True)
            for fq in fastqs:
                try:
                    fq.unlink()
                except OSError:
                    pass
            continue

        passed_runs.append(row)

        if not skip_fastqc:
            run_fastqc(trimmed, fastqc_dir, threads=min(4, args.star_threads))

        for fq in fastqs:
            try:
                fq.unlink()
                print(f"[ZEBRAFISH][qc] Cleaned up raw {fq.name}", flush=True)
            except OSError:
                pass

    if not passed_runs:
        raise SystemExit("No samples passed QC — nothing to write.")

    run_multiqc(qc_dir, outdir)

    report_dir = qc_dir / "reports"
    for html in sorted(report_dir.glob("*.fastp.html")):
        shutil.copy2(html, outdir / html.name)

    write_qc_metadata(
        out_dir=outdir, manifest=manifest, all_runs=runs,
        passed_runs=passed_runs, qc_stats=qc_stats_all, qc_dir=qc_dir,
    )
    print(f"[ZEBRAFISH][qc] Completed. Reports in {outdir}", flush=True)


# ---------------------------------------------------------------------------
# Array-compatible QC substages
# ---------------------------------------------------------------------------


def run_qc_array(args: argparse.Namespace) -> None:
    """Process a single sample for QC (download + fastp + optional FastQC).

    Designed to be called from a Slurm array task where --task_index ==
    $SLURM_ARRAY_TASK_ID.  Writes a per-sample result JSON so qc_finalize
    can aggregate later.
    """
    workdir = Path(args.workdir).resolve()
    manifest = Path(args.manifest).resolve()
    if not manifest.exists():
        raise SystemExit(f"Manifest not found: {manifest}")

    runs = parse_manifest(manifest, args.max_runs)
    if not runs:
        raise SystemExit("No runs found in zebrafish manifest.")
    idx = args.task_index
    if idx < 0 or idx >= len(runs):
        raise SystemExit(f"task_index {idx} out of range (0..{len(runs)-1})")

    row = runs[idx]
    run_id = row["run_id"]
    fastq_dir = workdir / "inputs" / "fastq"
    qc_dir = workdir / "qc"
    results_dir = qc_dir / "results_per_sample"
    results_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir.mkdir(parents=True, exist_ok=True)

    result_json = results_dir / f"{run_id}.json"
    if result_json.exists():
        print(f"[ZEBRAFISH][qc_array] {run_id} already processed — skipping", flush=True)
        return

    skip_fastqc = getattr(args, "skip_fastqc", False)
    fastqc_dir = qc_dir / "fastqc"

    print(f"[ZEBRAFISH][qc_array] Processing sample {idx}: {run_id}", flush=True)
    fastqs = download_fastqs(row, fastq_dir)
    if len(fastqs) > 2:
        fastqs = fastqs[:2]

    trimmed, qc_stats = run_fastp(
        run_id, fastqs, qc_dir,
        threads=min(4, args.star_threads),
    )

    passed = bool(trimmed)
    if passed and not skip_fastqc:
        run_fastqc(trimmed, fastqc_dir, threads=min(4, args.star_threads))

    for fq in fastqs:
        try:
            fq.unlink()
        except OSError:
            pass

    result = {
        "run_id": run_id,
        "condition": row["condition"],
        "library_layout": row.get("library_layout", "SINGLE"),
        "passed": passed,
        "qc_stats": qc_stats,
        "trimmed_paths": [str(p) for p in trimmed] if trimmed else [],
    }
    with open(result_json, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"[ZEBRAFISH][qc_array] {run_id}: {'PASSED' if passed else 'FAILED'}", flush=True)


def run_qc_finalize(args: argparse.Namespace) -> None:
    """Aggregate per-sample QC results and run MultiQC."""
    workdir = Path(args.workdir).resolve()
    manifest = Path(args.manifest).resolve()
    qc_dir = workdir / "qc"
    outdir = workdir / "out"
    results_dir = qc_dir / "results_per_sample"

    runs = parse_manifest(manifest, args.max_runs)
    if not runs:
        raise SystemExit("No runs found in zebrafish manifest.")

    qc_stats_all: Dict[str, Dict] = {}
    passed_runs: List[Dict[str, str]] = []

    for row in runs:
        rid = row["run_id"]
        rj = results_dir / f"{rid}.json"
        if not rj.exists():
            print(f"[ZEBRAFISH][qc_finalize] No result for {rid} — treating as failed", flush=True)
            continue
        with open(rj) as fh:
            result = json.load(fh)
        qc_stats_all[rid] = result.get("qc_stats", {})
        if result.get("passed"):
            passed_runs.append(row)

    if not passed_runs:
        raise SystemExit("No samples passed QC — nothing to write.")

    print(f"[ZEBRAFISH][qc_finalize] {len(passed_runs)}/{len(runs)} samples passed QC", flush=True)

    outdir.mkdir(parents=True, exist_ok=True)

    run_multiqc(qc_dir, outdir)

    report_dir = qc_dir / "reports"
    if report_dir.exists():
        for html in sorted(report_dir.glob("*.fastp.html")):
            try:
                shutil.copy2(html, outdir / html.name)
            except OSError:
                pass

    write_qc_metadata(
        out_dir=outdir, manifest=manifest, all_runs=runs,
        passed_runs=passed_runs, qc_stats=qc_stats_all, qc_dir=qc_dir,
    )
    print(f"[ZEBRAFISH][qc_finalize] Completed. Reports in {outdir}", flush=True)


# ---------------------------------------------------------------------------
# Array-compatible Stage1 substages
# ---------------------------------------------------------------------------


def run_stage1_setup(args: argparse.Namespace) -> None:
    """Build the STAR index if needed and write the manifest to the workdir."""
    workdir = Path(args.workdir).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf = Path(args.gencode_gtf).resolve()
    star_index = Path(args.star_index).resolve() if args.star_index else (workdir / "refs" / "star_index")

    if not genome_fasta.exists():
        raise SystemExit(f"Genome FASTA not found: {genome_fasta}")
    if not gtf.exists():
        raise SystemExit(f"GTF not found: {gtf}")

    ensure_star_index(genome_fasta, gtf, star_index, args.star_threads)
    print(f"[ZEBRAFISH][stage1_setup] STAR index ready at {star_index}", flush=True)


def run_stage1_array(args: argparse.Namespace) -> None:
    """Align a single sample via STAR. Designed for Slurm array tasks."""
    workdir = Path(args.workdir).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf = Path(args.gencode_gtf).resolve()
    star_index = Path(args.star_index).resolve() if args.star_index else (workdir / "refs" / "star_index")
    align_dir = workdir / "star_align"
    star_sj_root = workdir / "star_sj"

    qc_workdir = Path(args.qc_workdir).resolve() if getattr(args, "qc_workdir", "") else None
    use_qc = qc_workdir is not None

    if use_qc:
        qc_meta_path = qc_workdir / "out" / "qc_report.json"
        if not qc_meta_path.exists():
            raise SystemExit(f"QC metadata not found: {qc_meta_path}")
        with open(qc_meta_path, "r") as fh:
            qc_meta = json.load(fh)
        passed_info = qc_meta.get("passed_runs", [])
        trimmed_map = qc_meta.get("trimmed_fastqs", {})
        runs = [{"run_id": r["run_id"], "condition": r["condition"]} for r in passed_info]
    else:
        manifest = Path(args.manifest).resolve()
        if not manifest.exists():
            raise SystemExit(f"Manifest not found: {manifest}")
        runs = parse_manifest(manifest, args.max_runs)

    if not runs:
        raise SystemExit("No runs found.")
    idx = args.task_index
    if idx < 0 or idx >= len(runs):
        raise SystemExit(f"task_index {idx} out of range (0..{len(runs)-1})")

    row = runs[idx]
    run_id = row["run_id"]
    cond = safe_name(row.get("condition", "zebrafish"))
    cond_dir = star_sj_root / cond
    cond_dir.mkdir(parents=True, exist_ok=True)
    align_dir.mkdir(parents=True, exist_ok=True)

    out_sj = cond_dir / f"{run_id}.SJ.out.tab"
    if out_sj.exists() and out_sj.stat().st_size > 0:
        print(f"[ZEBRAFISH][stage1_array] {run_id} SJ already exists — skipping", flush=True)
        return

    fastq_dir = workdir / "inputs" / "fastq"
    qc_dir = workdir / "qc"
    skip_qc = getattr(args, "skip_qc", False)

    if use_qc:
        trimmed_paths = [Path(p) for p in trimmed_map.get(run_id, [])]
        missing = [p for p in trimmed_paths if not p.exists()]
        if not trimmed_paths or missing:
            print(f"[ZEBRAFISH][stage1_array] Trimmed FASTQs missing for {run_id} — skipping", flush=True)
            return
        align_inputs = trimmed_paths
    else:
        fastq_dir.mkdir(parents=True, exist_ok=True)
        print(f"[ZEBRAFISH][stage1_array] Downloading FASTQ for {run_id}", flush=True)
        fastqs = download_fastqs(row, fastq_dir)
        if len(fastqs) > 2:
            fastqs = fastqs[:2]
        align_inputs = fastqs

        if not skip_qc:
            trimmed, _ = run_fastp(
                run_id, fastqs, qc_dir,
                threads=min(4, args.star_threads),
            )
            if not trimmed:
                print(f"[ZEBRAFISH][stage1_array] {run_id} failed inline QC — skipping", flush=True)
                for fq in fastqs:
                    try:
                        fq.unlink()
                    except OSError:
                        pass
                return
            align_inputs = trimmed

    print(f"[ZEBRAFISH][stage1_array] Aligning {run_id}", flush=True)
    prefix = align_dir / f"{run_id}."
    sj = align_run(row, align_inputs, star_index, prefix, args.star_threads)
    shutil.copy2(sj, out_sj)

    for fq in align_inputs:
        try:
            fq.unlink()
        except OSError:
            pass
    if not use_qc and not skip_qc:
        for fq in fastqs:
            try:
                fq.unlink()
            except OSError:
                pass

    print(f"[ZEBRAFISH][stage1_array] {run_id} aligned -> {out_sj}", flush=True)


def run_stage1_finalize(args: argparse.Namespace) -> None:
    """Aggregate SJ files produced by stage1_array tasks into stage1 metadata."""
    workdir = Path(args.workdir).resolve()
    star_sj_root = workdir / "star_sj"
    samples_tsv = workdir / "inputs" / "samples.tsv"
    outdir = workdir / "out"

    qc_workdir = Path(args.qc_workdir).resolve() if getattr(args, "qc_workdir", "") else None
    use_qc = qc_workdir is not None

    if use_qc:
        qc_meta_path = qc_workdir / "out" / "qc_report.json"
        with open(qc_meta_path, "r") as fh:
            qc_meta = json.load(fh)
        manifest_path = Path(qc_meta.get("manifest", ""))
        runs = [{"run_id": r["run_id"], "condition": r["condition"]}
                for r in qc_meta.get("passed_runs", [])]
    else:
        manifest = Path(args.manifest).resolve()
        manifest_path = manifest
        runs = parse_manifest(manifest, args.max_runs)

    sj_files: List[Path] = []
    aligned_runs: List[Dict[str, str]] = []

    for row in runs:
        rid = row["run_id"]
        cond = safe_name(row.get("condition", "zebrafish"))
        sj_path = star_sj_root / cond / f"{rid}.SJ.out.tab"
        if sj_path.exists() and sj_path.stat().st_size > 0:
            sj_files.append(sj_path)
            aligned_runs.append(row)
        else:
            print(f"[ZEBRAFISH][stage1_finalize] Missing SJ for {rid} — skipping", flush=True)

    if not sj_files:
        raise SystemExit("No SJ files found — nothing to finalize.")

    print(f"[ZEBRAFISH][stage1_finalize] {len(sj_files)}/{len(runs)} samples have SJ files", flush=True)

    write_samples_tsv(aligned_runs, samples_tsv)
    write_stage1_metadata(
        out_dir=outdir, manifest=manifest_path, samples_tsv=samples_tsv,
        sj_files=sj_files, runs=aligned_runs,
    )
    print(f"[ZEBRAFISH][stage1_finalize] Completed. Metadata: {outdir / 'stage1_alignment.json'}", flush=True)


def ensure_star_index(genome_fasta: Path, gtf: Path, star_index: Path, threads: int) -> None:
    if (star_index / "Genome").exists():
        print(f"[ZEBRAFISH] Reusing existing STAR index at {star_index}", flush=True)
        return
    star_index.mkdir(parents=True, exist_ok=True)
    print(f"[ZEBRAFISH] Building STAR index at {star_index} ...", flush=True)
    run(
        [
            "STAR",
            "--runThreadN",
            str(threads),
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(star_index),
            "--genomeFastaFiles",
            str(genome_fasta),
            "--sjdbGTFfile",
            str(gtf),
            "--limitGenomeGenerateRAM",
            "150000000000",
        ]
    )


def align_run(run_row: Dict[str, str], fastqs: List[Path], star_index: Path, out_prefix: Path, threads: int) -> Path:
    cmd = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(star_index),
        "--readFilesIn",
    ]
    cmd.extend(str(p) for p in fastqs)
    cmd.extend(
        [
            "--readFilesCommand",
            "zcat",
            "--outFileNamePrefix",
            str(out_prefix),
            "--outSAMtype",
            "None",
        ]
    )
    run(cmd)
    sj = Path(str(out_prefix) + "SJ.out.tab")
    if not sj.exists():
        raise RuntimeError(f"STAR did not produce SJ.out.tab for {run_row['run_id']}")
    return sj


def write_samples_tsv(samples: List[Dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample", "condition"])
        for s in samples:
            writer.writerow([s["run_id"], s["condition"]])


def write_stage1_metadata(
    *,
    out_dir: Path,
    manifest: Path,
    samples_tsv: Path,
    sj_files: List[Path],
    runs: List[Dict[str, str]],
    qc_stats: Dict[str, Dict] = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sj_list_path = out_dir / "stage1_sj_files.txt"
    with open(sj_list_path, "w") as fh:
        for p in sj_files:
            fh.write(str(p) + "\n")

    stages = sorted({r.get("condition", "Unknown_Stage") for r in runs})
    payload: Dict = {
        "stage": "stage1_alignment",
        "manifest": str(manifest),
        "samples_tsv": str(samples_tsv),
        "sj_file_list": str(sj_list_path),
        "n_runs": len(runs),
        "n_stages": len(stages),
        "stages": stages,
    }
    if qc_stats:
        payload["qc_summary"] = qc_stats
    with open(out_dir / "stage1_alignment.json", "w") as fh:
        json.dump(payload, fh, indent=2)
    with open(out_dir / "summary.json", "w") as fh:
        summary: Dict = {
            "stage": "stage1_alignment",
            "n_runs": len(runs),
            "n_stages": len(stages),
            "stages": stages,
            "artifacts": {
                "stage1_alignment": "stage1_alignment.json",
                "sj_file_list": "stage1_sj_files.txt",
            },
            "notes": "Stage1 complete. Submit stage2 to run LeafCutter2 from stored STAR junctions.",
        }
        if qc_stats:
            summary["qc_summary"] = qc_stats
        json.dump(summary, fh, indent=2)


def validate_gtf_for_stage2(args: argparse.Namespace, gtf: Path, genome_fasta: Path) -> None:
    validator = Path(__file__).resolve().parent / "validate_gtf_for_lc2.py"
    cmd = [
        sys.executable,
        str(validator),
        str(gtf),
        "--fasta",
        str(genome_fasta),
        "--gene-name-tag",
        args.gene_name_tag,
        "--transcript-name-tag",
        args.transcript_name_tag,
        "--gene-type-tag",
        args.gene_type_tag,
        "--transcript-type-tag",
        args.transcript_type_tag,
        "--min-tag-frac",
        str(args.gtf_min_tag_frac),
    ]
    print("[ZEBRAFISH][stage2] Validating GTF for LC2", flush=True)
    run(cmd)


def count_sj_rows(paths: List[Path]) -> int:
    total = 0
    for path in paths:
        try:
            with open(path, "r") as fh:
                for _ in fh:
                    total += 1
        except OSError:
            pass
    return total


def validate_lc2_classification(lc2_dir: Path, prefix: str, max_unknown_frac: float) -> Dict[str, object]:
    path = lc2_dir / "clustering" / f"{prefix}_lc2_junction_classifications.txt"
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"LC2 classification file missing or empty: {path}")

    total = unknown = 0
    top: Dict[str, int] = {}
    with open(path, "r") as fh:
        header = next(fh, "")
        for line in fh:
            if not line.strip():
                continue
            total += 1
            gene = line.split("\t", 1)[0]
            if gene == "unknown":
                unknown += 1
            top[gene] = top.get(gene, 0) + 1
    if total == 0:
        raise SystemExit(f"LC2 classification file has no data rows: {path}")
    unknown_frac = unknown / total
    top_genes = sorted(top.items(), key=lambda kv: kv[1], reverse=True)[:10]
    print(
        f"[ZEBRAFISH][stage2] LC2 classifications: rows={total} unknown={unknown} "
        f"unknown_frac={unknown_frac:.4f} top_genes={top_genes}",
        flush=True,
    )
    if unknown_frac > max_unknown_frac:
        raise SystemExit(
            f"LC2 classification unknown fraction {unknown_frac:.3f} exceeds threshold {max_unknown_frac:.3f}"
        )
    return {
        "classification_rows": total,
        "unknown_rows": unknown,
        "unknown_fraction": unknown_frac,
        "top_genes": top_genes,
        "classification_file": str(path),
    }


def run_stage1(args: argparse.Namespace) -> None:
    """Stage 1 (alignment). When --qc_workdir is given, reads trimmed FASTQs from QC output;
    otherwise downloads fresh and optionally runs inline fastp."""
    workdir = Path(args.workdir).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf = Path(args.gencode_gtf).resolve()
    star_index = Path(args.star_index).resolve() if args.star_index else (workdir / "refs" / "star_index")
    align_dir = workdir / "star_align"
    star_sj_root = workdir / "star_sj"
    samples_tsv = workdir / "inputs" / "samples.tsv"
    outdir = workdir / "out"

    qc_workdir = Path(args.qc_workdir).resolve() if getattr(args, "qc_workdir", "") else None
    use_qc = qc_workdir is not None

    if not genome_fasta.exists():
        raise SystemExit(f"Genome FASTA not found: {genome_fasta}")
    if not gtf.exists():
        raise SystemExit(f"GTF not found: {gtf}")

    align_dir.mkdir(parents=True, exist_ok=True)
    star_sj_root.mkdir(parents=True, exist_ok=True)

    if use_qc:
        qc_meta_path = qc_workdir / "out" / "qc_report.json"
        if not qc_meta_path.exists():
            raise SystemExit(f"QC metadata not found: {qc_meta_path}")
        with open(qc_meta_path, "r") as fh:
            qc_meta = json.load(fh)
        passed_info = qc_meta.get("passed_runs", [])
        trimmed_map = qc_meta.get("trimmed_fastqs", {})
        manifest_path = Path(qc_meta.get("manifest", ""))
        runs = [{"run_id": r["run_id"], "condition": r["condition"]} for r in passed_info]
        if not runs:
            raise SystemExit("QC metadata has no passed runs.")
        print(f"[ZEBRAFISH][stage1] Using {len(runs)} QC-passed samples from {qc_meta_path}", flush=True)
    else:
        manifest = Path(args.manifest).resolve()
        if not manifest.exists():
            raise SystemExit(f"Manifest not found: {manifest}")
        manifest_path = manifest
        runs = parse_manifest(manifest, args.max_runs)
        if not runs:
            raise SystemExit("No runs found in zebrafish manifest.")

    print(f"[ZEBRAFISH][stage1] Runs to align: {len(runs)}", flush=True)
    ensure_star_index(genome_fasta, gtf, star_index, args.star_threads)

    fastq_dir = workdir / "inputs" / "fastq"
    qc_dir = workdir / "qc"
    skip_qc = getattr(args, "skip_qc", False)

    sj_files: List[Path] = []
    qc_stats_all: Dict[str, Dict] = {}
    aligned_runs: List[Dict[str, str]] = []

    for row in runs:
        run_id = row["run_id"]
        cond = safe_name(row["condition"])
        cond_dir = star_sj_root / cond
        cond_dir.mkdir(parents=True, exist_ok=True)

        if use_qc:
            trimmed_paths = [Path(p) for p in trimmed_map.get(run_id, [])]
            missing = [p for p in trimmed_paths if not p.exists()]
            if not trimmed_paths or missing:
                print(f"[ZEBRAFISH][stage1] Trimmed FASTQs missing for {run_id} — skipping", flush=True)
                continue
            align_inputs = trimmed_paths
        else:
            fastq_dir.mkdir(parents=True, exist_ok=True)
            print(f"[ZEBRAFISH][stage1] Downloading FASTQ for {run_id}", flush=True)
            fastqs = download_fastqs(row, fastq_dir)
            if len(fastqs) > 2:
                fastqs = fastqs[:2]
            align_inputs = fastqs

            if not skip_qc:
                trimmed, qc_stats = run_fastp(
                    run_id, fastqs, qc_dir,
                    threads=min(4, args.star_threads),
                )
                qc_stats_all[run_id] = qc_stats
                if not trimmed:
                    print(f"[ZEBRAFISH][stage1] Skipping {run_id} (failed inline QC)", flush=True)
                    for fq in fastqs:
                        try:
                            fq.unlink()
                        except OSError:
                            pass
                    continue
                align_inputs = trimmed

        print(f"[ZEBRAFISH][stage1] Aligning {run_id}", flush=True)
        prefix = align_dir / f"{run_id}."
        sj = align_run(row, align_inputs, star_index, prefix, args.star_threads)

        out_sj = cond_dir / f"{run_id}.SJ.out.tab"
        shutil.copy2(sj, out_sj)
        sj_files.append(out_sj)
        aligned_runs.append(row)

        for fq in align_inputs:
            try:
                fq.unlink()
                print(f"[ZEBRAFISH][stage1] Cleaned up {fq.name}", flush=True)
            except OSError:
                pass
        if not use_qc and not skip_qc:
            for fq in fastqs:
                try:
                    fq.unlink()
                except OSError:
                    pass

    if not sj_files:
        raise SystemExit("No samples aligned — nothing to write.")

    write_samples_tsv(aligned_runs, samples_tsv)
    write_stage1_metadata(
        out_dir=outdir, manifest=manifest_path, samples_tsv=samples_tsv,
        sj_files=sj_files, runs=aligned_runs, qc_stats=qc_stats_all if qc_stats_all else None,
    )
    print(f"[ZEBRAFISH][stage1] Completed. Metadata: {outdir / 'stage1_alignment.json'}", flush=True)


def run_stage2(args: argparse.Namespace) -> None:
    base_workdir = Path(args.workdir).resolve()
    is_smoke = args.stage == "stage2_smoke"
    workdir = (base_workdir / "smoke_lc2") if is_smoke else base_workdir
    stage1_workdir = Path(args.stage1_workdir).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()
    gtf = Path(args.gencode_gtf).resolve()

    stage1_meta = stage1_workdir / "out" / "stage1_alignment.json"
    if not stage1_meta.exists():
        raise SystemExit(f"Stage1 metadata not found: {stage1_meta}")
    with open(stage1_meta, "r") as fh:
        meta = json.load(fh)

    sj_list_path = Path(meta.get("sj_file_list", ""))
    samples_tsv = Path(meta.get("samples_tsv", ""))
    if not sj_list_path.exists():
        raise SystemExit(f"Stage1 SJ list missing: {sj_list_path}")
    if not samples_tsv.exists():
        raise SystemExit(f"Stage1 samples.tsv missing: {samples_tsv}")
    if not genome_fasta.exists():
        raise SystemExit(f"Genome FASTA not found: {genome_fasta}")
    if not gtf.exists():
        raise SystemExit(f"GTF not found: {gtf}")

    validate_gtf_for_stage2(args, gtf, genome_fasta)

    sj_files = [Path(x.strip()) for x in sj_list_path.read_text().splitlines() if x.strip()]
    if not sj_files:
        raise SystemExit("Stage1 SJ list is empty.")
    if is_smoke:
        n = max(1, args.smoke_n_sj_files)
        sj_files = sj_files[:n]
        print(f"[ZEBRAFISH][stage2_smoke] Using first {len(sj_files)} SJ files", flush=True)

    selected_stages = [t.strip() for t in (args.selected_stages or "").split(",") if t.strip()]
    if selected_stages:
        keep_runs = set()
        selected_set = set(selected_stages)
        filtered_rows = []
        with open(samples_tsv, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cond = (row.get("condition") or "").strip()
                sample = (row.get("sample") or "").strip()
                if cond in selected_set:
                    filtered_rows.append({"sample": sample, "condition": cond})
                    keep_runs.add(sample)
        if not keep_runs:
            raise SystemExit(
                f"No runs matched selected stages: {','.join(selected_stages)}"
            )
        filtered_sj = []
        for p in sj_files:
            rid = p.name.split(".SJ.out.tab")[0]
            if rid in keep_runs:
                filtered_sj.append(p)
        if not filtered_sj:
            raise SystemExit("No SJ files matched selected stage runs.")
        sj_files = filtered_sj
        local_samples_tsv = workdir / "inputs" / "samples.tsv"
        local_samples_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(local_samples_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["sample", "condition"], delimiter="\t")
            writer.writeheader()
            writer.writerows(filtered_rows)
        samples_tsv = local_samples_tsv
        print(
            f"[ZEBRAFISH][stage2] Stage filter applied: {','.join(selected_stages)} ({len(keep_runs)} runs)",
            flush=True,
        )

    total_sj_rows = count_sj_rows(sj_files)
    print(
        f"[ZEBRAFISH][stage2] Running LC2 on {len(sj_files)} SJ files "
        f"({total_sj_rows} STAR SJ rows)",
        flush=True,
    )
    # Persist the resolved SJ list path for stage3 consumption (matches the path
    # stage1 wrote, but re-emit here in case stage2 was invoked with a stage filter).
    stage2_sj_list = workdir / "out" / "stage2_sj_files.txt"
    stage2_sj_list.parent.mkdir(parents=True, exist_ok=True)
    with open(stage2_sj_list, "w") as fh:
        for p in sj_files:
            fh.write(str(p) + "\n")

    lc2_cmd = [
        sys.executable,
        str(Path(__file__).resolve().parent / "lc2_pipeline.py"),
        "--workdir",
        str(workdir),
        "--prefix",
        args.prefix,
        "--leafcutter_repo",
        args.leafcutter_repo,
        "--leafcutter2_repo",
        args.leafcutter2_repo,
        "--genome_fasta",
        str(genome_fasta),
        "--gencode_gtf",
        str(gtf),
        "--min_reads",
        str(args.min_reads),
        "--max_intron_len",
        str(args.max_intron_len),
        "--samples_tsv",
        str(samples_tsv),
        "--gene_type_tag",
        args.gene_type_tag,
        "--transcript_type_tag",
        args.transcript_type_tag,
        "--gene_name_tag",
        args.gene_name_tag,
        "--transcript_name_tag",
        args.transcript_name_tag,
        "--star_sj",
    ]
    lc2_cmd.extend(str(p) for p in sj_files)
    run(lc2_cmd)
    metrics = validate_lc2_classification(
        workdir / "lc2", args.prefix, args.max_unknown_gene_fraction
    )
    outdir = workdir / "out"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / ("stage2_smoke_qc.json" if is_smoke else "stage2_qc.json"), "w") as fh:
        json.dump(
            {
                "stage": "stage2_smoke" if is_smoke else "stage2",
                "n_sj_files": len(sj_files),
                "n_star_sj_rows": total_sj_rows,
                "gtf": str(gtf),
                "genome_fasta": str(genome_fasta),
                **metrics,
            },
            fh,
            indent=2,
        )


# ---------------------------------------------------------------------------
# Stage 3 — POISEN PE differential inclusion
# ---------------------------------------------------------------------------


def _resolve_pe_list(workdir: Path, override: str) -> Path:
    """Return the POISEN PE TSV path.

    The per-project artifact lives at ``<workdir>/refs/poisen_pe.tsv`` (analogous
    to ``<workdir>/inputs/samples.tsv``).  An ``--pe_list`` CLI override wins.
    """
    if override:
        p = Path(override).resolve()
        if not p.exists():
            raise SystemExit(f"--pe_list given but file not found: {p}")
        return p
    default = workdir / "refs" / "poisen_pe.tsv"
    if not default.exists():
        raise SystemExit(
            f"POISEN PE list not found at {default}. "
            "Place Talha's POISEN TSV there or pass --pe_list <path> explicitly."
        )
    return default


def _resolve_lc2_outputs(stage2_workdir: Optional[Path], lc2_prefix: str) -> Dict[str, Optional[Path]]:
    """Locate the stage2 LC2 cross-reference files when available.

    Tries ``<lc2_prefix>_*`` first, then falls back to a glob across known
    locations so we can still pick up cross-reference data when stage3 is being
    pointed at a stage2 workdir whose prefix differs from ours (e.g. webapp's
    ``web_run`` vs the pipeline's default ``zebrafish_run``).
    """
    if stage2_workdir is None:
        return {"classifications": None, "cluster_ratios": None, "samples_tsv": None}

    candidate_dirs = [
        stage2_workdir / "lc2" / "clustering",
        stage2_workdir / "lc2",
    ]

    def _first_existing(name_for_prefix, glob_pattern):
        for d in candidate_dirs:
            if not d.exists():
                continue
            p = d / name_for_prefix
            if p.exists():
                return p
        for d in candidate_dirs:
            if not d.exists():
                continue
            hits = sorted(d.glob(glob_pattern))
            if hits:
                return hits[0]
        return None

    classif = _first_existing(
        f"{lc2_prefix}_lc2_junction_classifications.txt",
        "*_lc2_junction_classifications.txt",
    )
    cluster_ratios = _first_existing(
        f"{lc2_prefix}_lc2.cluster_ratios.gz",
        "*_lc2.cluster_ratios.gz",
    )
    if cluster_ratios is None:
        cluster_ratios = _first_existing(
            f"{lc2_prefix}_lc2_cluster_ratios.gz",
            "*_lc2_cluster_ratios.gz",
        )

    samples_tsv = stage2_workdir / "inputs" / "samples.tsv"
    return {
        "classifications": classif,
        "cluster_ratios": cluster_ratios,
        "samples_tsv": samples_tsv if samples_tsv.exists() else None,
    }


def run_stage3(args: argparse.Namespace) -> None:
    """Stage 3 — POISEN PE differential inclusion analysis.

    Resolves stage1 SJ files (and optionally stage2 LC2 outputs) for the
    cross-reference columns, then delegates to ``pe_inclusion_analysis.py``.
    Stage3 is intentionally **opt-in** and can be re-run cheaply with an updated
    PE list without touching stage1/stage2 artifacts.
    """
    workdir = Path(args.workdir).resolve()
    stage1_workdir = Path(args.stage1_workdir).resolve()
    stage2_workdir = Path(args.stage2_workdir).resolve() if args.stage2_workdir else None
    outdir = workdir / "out"
    pe_outdir = outdir / "pe_inclusion"
    pe_outdir.mkdir(parents=True, exist_ok=True)

    stage1_meta = stage1_workdir / "out" / "stage1_alignment.json"
    if not stage1_meta.exists():
        raise SystemExit(f"Stage1 metadata not found: {stage1_meta}")
    with open(stage1_meta, "r") as fh:
        meta = json.load(fh)
    sj_list_path = Path(meta.get("sj_file_list", ""))
    samples_tsv = Path(meta.get("samples_tsv", ""))
    if not sj_list_path.exists():
        raise SystemExit(f"Stage1 SJ list missing: {sj_list_path}")
    if not samples_tsv.exists():
        raise SystemExit(f"Stage1 samples.tsv missing: {samples_tsv}")

    pe_list = _resolve_pe_list(workdir, args.pe_list)

    longorf_tsv = Path(args.longorf_tsv).resolve() if args.longorf_tsv else (workdir / "refs" / "longorf_ptc.tsv")
    longorf_arg: Optional[Path] = longorf_tsv if longorf_tsv.exists() else None

    lc2 = _resolve_lc2_outputs(stage2_workdir, args.prefix)
    lc2_samples_tsv = lc2["samples_tsv"] or samples_tsv

    chrom_map_arg: Optional[Path] = None
    if getattr(args, "chrom_map", ""):
        chrom_map_arg = Path(args.chrom_map).resolve()
    else:
        for guess in (
            workdir / "refs" / "assembly_report.txt",
            Path(args.genome_fasta).resolve().parent / "GRCz12tu_assembly_report.txt"
            if args.genome_fasta else None,
            Path(args.genome_fasta).resolve().with_suffix("").parent / (
                Path(args.genome_fasta).stem + "_assembly_report.txt"
            ) if args.genome_fasta else None,
        ):
            if guess is not None and guess.exists():
                chrom_map_arg = guess
                break

    cmd: List[str] = [
        "python3",
        str(Path(__file__).resolve().parent / "pe_inclusion_analysis.py"),
        "--pe-list", str(pe_list),
        "--sj-list", str(sj_list_path),
        "--samples-tsv", str(samples_tsv),
        "--outdir", str(pe_outdir),
        "--min-junction-reads", str(args.pe_min_junction_reads),
        "--min-locus-reads", str(args.pe_min_locus_reads),
        "--min-observability-frac", str(args.pe_min_observability_frac),
        "--min-timepoints-observed", str(args.pe_min_timepoints_observed),
        "--min-total-locus-reads", str(args.pe_min_total_locus_reads),
        "--min-dpsi", str(args.pe_min_dpsi),
        "--fdr-threshold", str(args.pe_fdr_threshold),
        "--max-intron-len", str(args.max_intron_len),
    ]
    if chrom_map_arg is not None:
        cmd += ["--chrom-map", str(chrom_map_arg)]
    if longorf_arg is not None:
        cmd += ["--longorf-tsv", str(longorf_arg)]
    if lc2["classifications"] is not None:
        cmd += ["--lc2-classifications", str(lc2["classifications"])]
    if lc2["cluster_ratios"] is not None:
        cmd += ["--lc2-cluster-ratios", str(lc2["cluster_ratios"])]
    if args.ds:
        cmd.append("--ds")
    # Auto-enable DM validation track when leafcutter repo is supplied.
    if args.leafcutter_repo:
        cmd += ["--leafcutter-repo", args.leafcutter_repo]

    print(f"[ZEBRAFISH][stage3] PE list: {pe_list}", flush=True)
    print(f"[ZEBRAFISH][stage3] LongORF: {longorf_arg or '(none)'}", flush=True)
    print(f"[ZEBRAFISH][stage3] LC2 classifications: {lc2['classifications'] or '(none)'}", flush=True)
    print(f"[ZEBRAFISH][stage3] LC2 cluster ratios: {lc2['cluster_ratios'] or '(none)'}", flush=True)
    run(cmd)

    pe_summary = pe_outdir / "summary.json"
    if pe_summary.exists():
        with open(pe_summary, "r") as fh:
            pe_payload = json.load(fh)
    else:
        pe_payload = {}

    out_meta = outdir / "stage3_pe_inclusion.json"
    with open(out_meta, "w") as fh:
        json.dump({
            "stage": "stage3_pe_inclusion",
            "stage1_workdir": str(stage1_workdir),
            "stage2_workdir": str(stage2_workdir) if stage2_workdir else None,
            "pe_list": str(pe_list),
            "longorf_tsv": str(longorf_arg) if longorf_arg else None,
            "lc2_classifications": str(lc2["classifications"]) if lc2["classifications"] else None,
            "lc2_cluster_ratios": str(lc2["cluster_ratios"]) if lc2["cluster_ratios"] else None,
            "pe_inclusion": pe_payload,
        }, fh, indent=2)

    summary_path = outdir / "summary.json"
    summary: Dict = {}
    if summary_path.exists():
        try:
            with open(summary_path, "r") as fh:
                summary = json.load(fh)
        except Exception:
            summary = {}
    summary["stage"] = "stage3_pe_inclusion"
    summary.setdefault("metrics", {})["pe_inclusion_analysis"] = {
        "summary_json": "pe_inclusion/summary.json",
        "candidates_tsv": "pe_inclusion/pe_developmental_candidates.tsv",
        "candidates_json": "pe_inclusion/pe_developmental_candidates.json",
        "events_tsv": "pe_inclusion/pe_events.tsv",
        "psi_tsv_gz": "pe_inclusion/pe_inclusion_psi.tsv.gz",
        "low_coverage_tsv": "pe_inclusion/pe_low_coverage.tsv",
        "heatmap_png": "pe_inclusion/pe_psi_heatmap.png",
        "psi_vs_time_png": "pe_inclusion/pe_psi_vs_time_top.png",
        "dpsi_volcano_png": "pe_inclusion/pe_dpsi_vs_pvalue.png",
        "spearman_distribution_png": "pe_inclusion/pe_spearman_distribution.png",
        "lc2_concordance_png": "pe_inclusion/pe_lc2_concordance.png",
        "n_candidates": pe_payload.get("n_candidates"),
        "n_events_in_psi": pe_payload.get("n_events_in_psi"),
        "ordered_timepoints": pe_payload.get("ordered_timepoints"),
        "stage_a": pe_payload.get("stage_a"),
        "stage_c": pe_payload.get("stage_c"),
        "top_hits": pe_payload.get("top_hits"),
    }
    summary.setdefault("artifacts", {})["pe_inclusion"] = "stage3_pe_inclusion.json"
    summary["notes"] = "Stage3 complete. POISEN PE differential inclusion candidates ready."
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[ZEBRAFISH][stage3] Completed. Metadata: {out_meta}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description="Run zebrafish end-to-end LC2 workflow.")
    ALL_STAGES = [
        "qc", "qc_array", "qc_finalize",
        "stage1", "stage1_setup", "stage1_array", "stage1_finalize",
        "stage2", "stage2_smoke", "stage3", "all", "all+pe",
    ]
    p.add_argument("--stage", choices=ALL_STAGES, required=True,
                   help=("Pipeline stage. 'all' runs stage1->stage2 only. "
                         "'all+pe' adds stage3 (requires --pe_list or "
                         "<workdir>/refs/poisen_pe.tsv). stage3 alone is opt-in "
                         "and re-runnable on its own. Array substages "
                         "(qc_array, stage1_array) process one sample via "
                         "--task_index (= $SLURM_ARRAY_TASK_ID)."))
    p.add_argument("--workdir", required=True)
    p.add_argument("--manifest", default="", help="Run manifest TSV from zebrafish_project_runs.py")
    p.add_argument("--qc_workdir", default="", help="Path to prior QC workdir for stage1 (reads trimmed FASTQs).")
    p.add_argument("--stage1_workdir", default="", help="Path to prior stage1 workdir for stage2/stage3.")
    p.add_argument("--stage2_workdir", default="", help="Path to prior stage2 workdir for stage3 (provides LC2 cross-reference columns).")
    p.add_argument("--prefix", default="zebrafish_run")
    p.add_argument("--leafcutter_repo", default="")
    p.add_argument("--leafcutter2_repo", default="")
    p.add_argument("--genome_fasta", default="", help="Genome FASTA (required for stage1/stage2, not for qc).")
    p.add_argument("--gencode_gtf", default="", help="GTF annotation (required for stage1/stage2, not for qc).")
    p.add_argument("--min_reads", type=int, default=50)
    p.add_argument("--max_intron_len", type=int, default=500000)
    p.add_argument("--star_threads", type=int, default=8)
    p.add_argument("--max_runs", type=int, default=0)
    p.add_argument("--star_index", default="", help="Optional STAR index dir; defaults to <workdir>/refs/star_index")
    p.add_argument("--gene_type_tag", default="gene_biotype")
    p.add_argument("--transcript_type_tag", default="transcript_biotype")
    p.add_argument("--gene_name_tag", default="gene_id")
    p.add_argument("--transcript_name_tag", default="transcript_id")
    p.add_argument("--selected_stages", default="", help="Comma-separated developmental stages to include in stage2 (default all).")
    p.add_argument("--smoke_n_sj_files", type=int, default=5,
                   help="For stage2_smoke, number of Stage1 SJ files to test.")
    p.add_argument("--max_unknown_gene_fraction", type=float, default=0.20,
                   help="Fail stage2/stage2_smoke when LC2 classifications exceed this unknown-gene fraction.")
    p.add_argument("--gtf_min_tag_frac", type=float, default=0.01,
                   help="Minimum fraction of GTF rows that must contain each requested LC2 tag.")
    p.add_argument("--skip_qc", action="store_true", help="Skip fastp QC/trimming in stage1 (align raw reads directly).")
    p.add_argument("--skip_fastqc", action="store_true", help="Skip FastQC in QC stage (fastp + MultiQC still run).")
    # ----- Stage 3 (POISEN PE differential inclusion) -----
    p.add_argument("--pe_list", default="",
                   help=("Override path to POISEN PE TSV (default: "
                         "<workdir>/refs/poisen_pe.tsv)."))
    p.add_argument("--longorf_tsv", default="",
                   help=("Override path to LongORF_PTC+_fromClustered.tsv "
                         "(default: <workdir>/refs/longorf_ptc.tsv if present)."))
    p.add_argument("--pe_min_junction_reads", type=int, default=3,
                   help="Per-junction read floor for stage3 PSI hygiene (mirrors LC2 pool_junc_reads).")
    p.add_argument("--pe_min_locus_reads", type=int, default=10,
                   help="Per-PE per-sample inclusion+skip read floor for stage3 "
                        "(LeafCutter2 default 10).")
    p.add_argument("--pe_min_observability_frac", type=float, default=0.60,
                   help="Stage3: event must be observed in >=this fraction of "
                        "developmental samples to be tested (LC2 default 0.60).")
    p.add_argument("--pe_min_timepoints_observed", type=int, default=3,
                   help="Stage3: event must have data in >=this many distinct "
                        "developmental timepoints to be tested.")
    p.add_argument("--pe_min_total_locus_reads", type=int, default=50,
                   help="Stage3: event-aggregate inclusion+skip reads floor; "
                        "LC2 'appreciable abundance' default 50.")
    p.add_argument("--pe_min_dpsi", type=float, default=0.10,
                   help="Stage3: |dPSI| effect-size bar between earliest and "
                        "latest observed timepoints (default 0.10).")
    p.add_argument("--pe_fdr_threshold", type=float, default=0.05,
                   help="Stage3: FDR threshold for combined-significance reporting.")
    p.add_argument("--ds", action="store_true",
                   help=("Stage3: also run leafcutter Dirichlet-multinomial test on "
                         "first vs last developmental timepoint as a validation track "
                         "(requires --leafcutter_repo). Auto-enabled when "
                         "--leafcutter_repo is provided."))
    p.add_argument("--chrom_map", default="",
                   help=("Stage3: NCBI assembly_report.txt (or 2-col TSV) mapping "
                         "POISEN's RefSeq accessions (e.g. NC_133177.1) to STAR "
                         "chrom names (e.g. chr2). Auto-detected from the genome "
                         "FASTA's directory if not provided."))
    p.add_argument("--task_index", type=int, default=-1,
                   help=("Array task index (0-based). Required for qc_array and "
                         "stage1_array. Typically set to $SLURM_ARRAY_TASK_ID."))
    args = p.parse_args()

    if args.stage == "qc":
        required = ["wget", "python3", "fastp"]
        for exe in required:
            check_exe(exe)
        if not args.manifest:
            raise SystemExit("--manifest is required for qc stage")
        run_qc(args)
        return

    if args.stage == "qc_array":
        for exe in ["wget", "python3", "fastp"]:
            check_exe(exe)
        if not args.manifest:
            raise SystemExit("--manifest is required for qc_array")
        if args.task_index < 0:
            raise SystemExit("--task_index is required for qc_array")
        run_qc_array(args)
        return

    if args.stage == "qc_finalize":
        check_exe("python3")
        if not args.manifest:
            raise SystemExit("--manifest is required for qc_finalize")
        run_qc_finalize(args)
        return

    if args.stage == "stage1_setup":
        check_exe("STAR")
        if not args.genome_fasta or not args.gencode_gtf:
            raise SystemExit("--genome_fasta and --gencode_gtf are required for stage1_setup")
        run_stage1_setup(args)
        return

    if args.stage == "stage1_array":
        for exe in ["STAR", "python3"]:
            check_exe(exe)
        if not args.genome_fasta or not args.gencode_gtf:
            raise SystemExit("--genome_fasta and --gencode_gtf are required for stage1_array")
        if args.task_index < 0:
            raise SystemExit("--task_index is required for stage1_array")
        if not args.qc_workdir and not args.manifest:
            raise SystemExit("--manifest or --qc_workdir is required for stage1_array")
        run_stage1_array(args)
        return

    if args.stage == "stage1_finalize":
        check_exe("python3")
        if not args.qc_workdir and not args.manifest:
            raise SystemExit("--manifest or --qc_workdir is required for stage1_finalize")
        run_stage1_finalize(args)
        return

    if args.stage in ("stage1", "all", "all+pe"):
        required = ["STAR", "python3"]
        if not args.qc_workdir:
            required.append("wget")
            if not args.skip_qc:
                required.append("fastp")
        for exe in required:
            check_exe(exe)
        if not args.qc_workdir and not args.manifest:
            raise SystemExit("--manifest or --qc_workdir is required for stage1")
        if not args.genome_fasta or not args.gencode_gtf:
            raise SystemExit("--genome_fasta and --gencode_gtf are required for stage1")
        run_stage1(args)
        if args.stage == "stage1":
            return
        # 'all' / 'all+pe': stage2 reads back from this same workdir
        if not args.stage1_workdir:
            args.stage1_workdir = args.workdir

    if args.stage in ("stage2", "stage2_smoke", "all", "all+pe"):
        check_exe("python3")
        if not args.stage1_workdir:
            raise SystemExit("--stage1_workdir is required for stage2")
        if not args.leafcutter_repo or not args.leafcutter2_repo:
            raise SystemExit("--leafcutter_repo and --leafcutter2_repo are required for stage2")
        if not args.genome_fasta or not args.gencode_gtf:
            raise SystemExit("--genome_fasta and --gencode_gtf are required for stage2")
        run_stage2(args)
        if args.stage in ("stage2", "stage2_smoke"):
            return
        if args.stage == "all":
            return
        # 'all+pe' falls through to stage3 with stage2_workdir defaulting to workdir
        if not args.stage2_workdir:
            args.stage2_workdir = args.workdir

    if args.stage in ("stage3", "all+pe"):
        check_exe("python3")
        if not args.stage1_workdir:
            raise SystemExit("--stage1_workdir is required for stage3")
        run_stage3(args)
        return


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
