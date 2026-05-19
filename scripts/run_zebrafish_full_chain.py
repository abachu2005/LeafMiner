#!/usr/bin/env python3
"""Submit a complete zebrafish pipeline as a Slurm dependency chain.

Stages:
  qc_setup (login download manifest)
  -> qc_array (--array=0-N per-sample fastp)
  -> qc_finalize (MultiQC)
  -> stage1_setup (ensure STAR index)
  -> stage1_array (--array=0-N per-sample STAR)
  -> stage1_finalize (aggregate SJ)
  -> stage2 (LeafCutter2)
  -> stage3 (POISEN PE inclusion)

With --reuse_stage1_workdir, the chain shortens to stage2 -> stage3 only
(reuses SJ files from a prior run).

Each sbatch sets --mail-type=END,FAIL so completion comes via email.
Slurm JIDs are recorded to a local JSON file.
"""

import argparse
import json
import shlex
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Dict, List, Optional


def ssh_cmd(host: str, user: str, cmd: str, check: bool = True) -> str:
    full = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
            f"{user}@{host}", f"bash --login -c {shlex.quote(cmd)}"]
    r = subprocess.run(full, capture_output=True, text=True, check=check)
    return r.stdout.strip()


def submit_sbatch(host: str, user: str, sbatch_body: str, dry_run: bool = False) -> str:
    """Submit sbatch script and return Slurm JID."""
    if dry_run:
        print("=== DRY RUN sbatch body ===")
        print(sbatch_body)
        print("===========================\n")
        return "DRY_RUN_JID"
    escaped = sbatch_body.replace("'", "'\\''")
    cmd = f"echo '{escaped}' | sbatch"
    result = ssh_cmd("", "", "")  # placeholder
    full = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
            f"{user}@{host}",
            f"bash --login -c 'echo {shlex.quote(sbatch_body)} | sbatch'"]
    r = subprocess.run(full, capture_output=True, text=True, check=True)
    for line in r.stdout.strip().splitlines():
        if "Submitted batch job" in line:
            return line.split()[-1]
    raise RuntimeError(f"Could not parse sbatch output: {r.stdout}")


def _submit_sbatch_via_file(host: str, user: str, sbatch_body: str,
                            remote_dir: str, script_name: str,
                            dry_run: bool = False) -> str:
    """Write sbatch to a remote file, then submit it."""
    if dry_run:
        print(f"=== DRY RUN: {script_name} ===")
        print(sbatch_body)
        print("===========================\n")
        return "DRY_RUN_JID"

    remote_path = f"{remote_dir}/{script_name}"
    escaped = sbatch_body.replace("\\", "\\\\").replace("'", "'\\''")
    write_cmd = f"mkdir -p {remote_dir} && cat > {remote_path} << 'SBATCH_EOF'\n{sbatch_body}\nSBATCH_EOF"
    subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
         f"{user}@{host}", write_cmd],
        check=True, capture_output=True, text=True,
    )
    r = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
         f"{user}@{host}", f"bash --login -c 'sbatch {remote_path}'"],
        capture_output=True, text=True, check=True,
    )
    for line in r.stdout.strip().splitlines():
        if "Submitted batch job" in line:
            return line.split()[-1]
    raise RuntimeError(f"Could not parse sbatch output: {r.stdout}\nstderr: {r.stderr}")


def build_preamble(
    account: str, partition: str, job_name: str,
    remote_run_dir: str, time: str = "04:00:00",
    mem: str = "16G", cpus: int = 1,
    mail_user: str = "", dependency: str = "",
    array_spec: str = "",
) -> List[str]:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --account={account}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --time={time}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={remote_run_dir}/slurm_%j_%a.out",
        f"#SBATCH --error={remote_run_dir}/slurm_%j_%a.err",
    ]
    if array_spec:
        lines.append(f"#SBATCH --array={array_spec}")
    if dependency:
        lines.append(f"#SBATCH --dependency={dependency}")
    if mail_user:
        lines.append(f"#SBATCH --mail-user={mail_user}")
        lines.append("#SBATCH --mail-type=END,FAIL")
    lines += [
        "",
        "set -eo pipefail",
        "module purge",
        "module load anaconda3/2022.05",
        "module load STAR/2.7.11b",
        "module load fastp/0.23.4",
        "module load samtools",
        "module load R/4.4.0",
        "",
    ]
    return lines


def main():
    p = argparse.ArgumentParser(description="Submit full zebrafish pipeline chain to Quest")
    p.add_argument("--run_name", required=True, help="Human-readable run name (e.g. run2_grcz12_longorf)")
    p.add_argument("--quest_host", default="login.quest.northwestern.edu",
                   help="Slurm login host (defaults to Northwestern Quest)")
    p.add_argument("--quest_user", required=True, help="SSH username on the Slurm host (e.g. your NetID)")
    p.add_argument("--quest_account", required=True, help="Slurm --account value")
    p.add_argument("--quest_partition", default="short")
    p.add_argument("--stage2_partition", default="",
                   help="Partition for Stage 2; defaults to --quest_partition.")
    p.add_argument("--stage3_partition", default="",
                   help="Partition for Stage 3; defaults to --quest_partition.")
    p.add_argument("--remote_repo_root", required=True,
                   help="Absolute path on the cluster where Leaf_Cutter is checked out (e.g. /projects/<account>/<user>/Leaf_Cutter)")
    p.add_argument("--mail_user", default="", help="Slurm --mail-user (optional)")

    p.add_argument("--projects", default="PRJEB7244", help="Comma-separated ENA project IDs")
    p.add_argument("--max_runs", type=int, default=0)
    p.add_argument("--array_concurrency", type=int, default=50, help="Max concurrent array tasks")

    p.add_argument("--genome_fasta", default="refs/GRCz12tu.fa")
    p.add_argument("--gencode_gtf", default="refs/GRCz12tu.ucsc.gtf")
    p.add_argument("--star_index", default="refs/star_index_grcz12tu")
    p.add_argument("--star_threads", type=int, default=8)

    p.add_argument("--prefix", default="zebrafish_run")
    p.add_argument("--leafcutter_repo", default="tools/leafcutter")
    p.add_argument("--leafcutter2_repo", default="tools/leafcutter2")
    p.add_argument("--gene_type_tag", default="gene_biotype")
    p.add_argument("--transcript_type_tag", default="transcript_biotype")
    p.add_argument("--gene_name_tag", default="gene_id")
    p.add_argument("--transcript_name_tag", default="transcript_id")
    p.add_argument("--stage2_time", default="04:00:00")
    p.add_argument("--stage3_time", default="04:00:00")
    p.add_argument("--max_unknown_gene_fraction", type=float, default=0.20)
    p.add_argument("--gtf_min_tag_frac", type=float, default=0.01)

    p.add_argument("--pe_list", default="UPDATED_STRAND_POISEN_HITS.tsv")
    p.add_argument("--longorf_tsv", default="LongORF_PTC+_fromClustered.tsv")

    p.add_argument("--reuse_stage1_workdir", default="",
                   help="Skip QC+alignment; reuse SJ files from this workdir. Chain starts at stage2.")

    p.add_argument("--dry_run", action="store_true", help="Print sbatch bodies without submitting")
    p.add_argument("--output_json", default="", help="Write chain JIDs to this JSON file")
    args = p.parse_args()

    host = args.quest_host
    user = args.quest_user
    rr = args.remote_repo_root.rstrip("/")
    stage2_partition = args.stage2_partition or args.quest_partition
    stage3_partition = args.stage3_partition or args.quest_partition
    run_id = str(uuid.uuid4())
    remote_run_dir = f"{rr}/jobs/{run_id}"
    reuse = args.reuse_stage1_workdir.strip()

    print(f"Run name: {args.run_name}")
    print(f"Run UUID: {run_id}")
    print(f"Remote dir: {remote_run_dir}")
    if reuse:
        print(f"Reusing stage1 from: {reuse}")

    if not args.dry_run:
        ssh_cmd(host, user, f"mkdir -p {remote_run_dir}/inputs")

    chain_jids: Dict[str, str] = {}

    if reuse:
        # ---------- REUSE PATH: stage2 -> stage3 only ----------
        stage2_body = "\n".join(build_preamble(
            account=args.quest_account, partition=stage2_partition,
            job_name=f"zf_s2_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time=args.stage2_time,
            mem="64G", cpus=args.star_threads, mail_user=args.mail_user,
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage2 \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --stage1_workdir {reuse} \\",
            f"  --prefix {args.prefix} \\",
            f"  --leafcutter_repo {rr}/{args.leafcutter_repo} \\",
            f"  --leafcutter2_repo {rr}/{args.leafcutter2_repo} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --gene_type_tag {args.gene_type_tag} \\",
            f"  --transcript_type_tag {args.transcript_type_tag} \\",
            f"  --gene_name_tag {args.gene_name_tag} \\",
            f"  --transcript_name_tag {args.transcript_name_tag} \\",
            f"  --max_unknown_gene_fraction {args.max_unknown_gene_fraction} \\",
            f"  --gtf_min_tag_frac {args.gtf_min_tag_frac} \\",
            f"  --star_threads {args.star_threads}",
        ]) + "\n"

        jid_s2 = _submit_sbatch_via_file(
            host, user, stage2_body, remote_run_dir, "stage2.sbatch", args.dry_run)
        chain_jids["stage2"] = jid_s2
        print(f"  stage2 JID: {jid_s2}")

        stage3_body = "\n".join(build_preamble(
            account=args.quest_account, partition=stage3_partition,
            job_name=f"zf_s3_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time=args.stage3_time,
            mem="64G", cpus=4, mail_user=args.mail_user,
            dependency=f"afterok:{jid_s2}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage3 \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --stage1_workdir {reuse} \\",
            f"  --stage2_workdir {remote_run_dir} \\",
            f"  --prefix {args.prefix} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --pe_list {rr}/{args.pe_list} \\",
            f"  --longorf_tsv {rr}/{args.longorf_tsv} \\",
            f"  --leafcutter_repo {rr}/{args.leafcutter_repo} \\",
            f"  --star_threads 4",
        ]) + "\n"

        jid_s3 = _submit_sbatch_via_file(
            host, user, stage3_body, remote_run_dir, "stage3.sbatch", args.dry_run)
        chain_jids["stage3"] = jid_s3
        print(f"  stage3 JID: {jid_s3}")

    else:
        # ---------- FULL PATH: generate manifest, then qc -> stage1 -> stage2 -> stage3 ----------

        # Step 0: Generate manifest on Quest
        projects = args.projects
        manifest_path = f"{remote_run_dir}/inputs/zebrafish_runs.tsv"
        if not args.dry_run:
            print("Generating manifest on Quest...")
            ssh_cmd(host, user,
                    f"cd {rr} && python3 {rr}/scripts/zebrafish_project_runs.py "
                    f"--project_ids \"{projects}\" --out \"{manifest_path}\" "
                    f"--max_runs {args.max_runs}")
            n_samples = ssh_cmd(host, user, f"tail -n +2 {manifest_path} | wc -l").strip()
        else:
            n_samples = "N"
        print(f"Manifest: {manifest_path} ({n_samples} samples)")
        n_int = int(n_samples) if n_samples.isdigit() else 276
        array_max = n_int - 1

        # Stage: qc_array
        qc_array_body = "\n".join(build_preamble(
            account=args.quest_account, partition=args.quest_partition,
            job_name=f"zf_qc_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time="02:00:00",
            mem="8G", cpus=4, mail_user=args.mail_user,
            array_spec=f"0-{array_max}%{args.array_concurrency}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage qc_array \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --manifest {manifest_path} \\",
            f"  --task_index $SLURM_ARRAY_TASK_ID \\",
            f"  --star_threads 4 \\",
            f"  --max_runs {args.max_runs}",
        ]) + "\n"

        jid_qc = _submit_sbatch_via_file(
            host, user, qc_array_body, remote_run_dir, "qc_array.sbatch", args.dry_run)
        chain_jids["qc_array"] = jid_qc
        print(f"  qc_array JID: {jid_qc}")

        # Stage: qc_finalize (afterany: tolerates partial QC failures)
        qc_fin_body = "\n".join(build_preamble(
            account=args.quest_account, partition=args.quest_partition,
            job_name=f"zf_qcf_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time="01:00:00",
            mem="16G", cpus=2, mail_user=args.mail_user,
            dependency=f"afterany:{jid_qc}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage qc_finalize \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --manifest {manifest_path} \\",
            f"  --max_runs {args.max_runs}",
        ]) + "\n"

        jid_qcf = _submit_sbatch_via_file(
            host, user, qc_fin_body, remote_run_dir, "qc_finalize.sbatch", args.dry_run)
        chain_jids["qc_finalize"] = jid_qcf
        print(f"  qc_finalize JID: {jid_qcf}")

        # Stage: stage1_setup
        s1_setup_body = "\n".join(build_preamble(
            account=args.quest_account, partition=args.quest_partition,
            job_name=f"zf_s1s_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time="02:00:00",
            mem="64G", cpus=args.star_threads, mail_user=args.mail_user,
            dependency=f"afterok:{jid_qcf}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage1_setup \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --star_index {rr}/{args.star_index} \\",
            f"  --star_threads {args.star_threads}",
        ]) + "\n"

        jid_s1s = _submit_sbatch_via_file(
            host, user, s1_setup_body, remote_run_dir, "stage1_setup.sbatch", args.dry_run)
        chain_jids["stage1_setup"] = jid_s1s
        print(f"  stage1_setup JID: {jid_s1s}")

        # Stage: stage1_array
        s1_array_body = "\n".join(build_preamble(
            account=args.quest_account, partition=args.quest_partition,
            job_name=f"zf_s1_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time="02:00:00",
            mem="48G", cpus=args.star_threads, mail_user=args.mail_user,
            dependency=f"afterok:{jid_s1s}",
            array_spec=f"0-{array_max}%{args.array_concurrency}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage1_array \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --qc_workdir {remote_run_dir} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --star_index {rr}/{args.star_index} \\",
            f"  --star_threads {args.star_threads} \\",
            f"  --task_index $SLURM_ARRAY_TASK_ID \\",
            f"  --max_runs {args.max_runs}",
        ]) + "\n"

        jid_s1a = _submit_sbatch_via_file(
            host, user, s1_array_body, remote_run_dir, "stage1_array.sbatch", args.dry_run)
        chain_jids["stage1_array"] = jid_s1a
        print(f"  stage1_array JID: {jid_s1a}")

        # Stage: stage1_finalize (afterany: tolerates partial alignment failures)
        s1_fin_body = "\n".join(build_preamble(
            account=args.quest_account, partition=args.quest_partition,
            job_name=f"zf_s1f_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time="00:30:00",
            mem="8G", cpus=1, mail_user=args.mail_user,
            dependency=f"afterany:{jid_s1a}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage1_finalize \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --qc_workdir {remote_run_dir} \\",
            f"  --manifest {manifest_path} \\",
            f"  --max_runs {args.max_runs}",
        ]) + "\n"

        jid_s1f = _submit_sbatch_via_file(
            host, user, s1_fin_body, remote_run_dir, "stage1_finalize.sbatch", args.dry_run)
        chain_jids["stage1_finalize"] = jid_s1f
        print(f"  stage1_finalize JID: {jid_s1f}")

        # Stage: stage2
        stage2_body = "\n".join(build_preamble(
            account=args.quest_account, partition=stage2_partition,
            job_name=f"zf_s2_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time=args.stage2_time,
            mem="64G", cpus=args.star_threads, mail_user=args.mail_user,
            dependency=f"afterok:{jid_s1f}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage2 \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --stage1_workdir {remote_run_dir} \\",
            f"  --prefix {args.prefix} \\",
            f"  --leafcutter_repo {rr}/{args.leafcutter_repo} \\",
            f"  --leafcutter2_repo {rr}/{args.leafcutter2_repo} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --gene_type_tag {args.gene_type_tag} \\",
            f"  --transcript_type_tag {args.transcript_type_tag} \\",
            f"  --gene_name_tag {args.gene_name_tag} \\",
            f"  --transcript_name_tag {args.transcript_name_tag} \\",
            f"  --max_unknown_gene_fraction {args.max_unknown_gene_fraction} \\",
            f"  --gtf_min_tag_frac {args.gtf_min_tag_frac} \\",
            f"  --star_threads {args.star_threads}",
        ]) + "\n"

        jid_s2 = _submit_sbatch_via_file(
            host, user, stage2_body, remote_run_dir, "stage2.sbatch", args.dry_run)
        chain_jids["stage2"] = jid_s2
        print(f"  stage2 JID: {jid_s2}")

        # Stage: stage3
        stage3_body = "\n".join(build_preamble(
            account=args.quest_account, partition=stage3_partition,
            job_name=f"zf_s3_{args.run_name[:10]}",
            remote_run_dir=remote_run_dir, time=args.stage3_time,
            mem="64G", cpus=4, mail_user=args.mail_user,
            dependency=f"afterok:{jid_s2}",
        ) + [
            f"cd {rr}",
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py \\",
            f"  --stage stage3 \\",
            f"  --workdir {remote_run_dir} \\",
            f"  --stage1_workdir {remote_run_dir} \\",
            f"  --stage2_workdir {remote_run_dir} \\",
            f"  --prefix {args.prefix} \\",
            f"  --genome_fasta {rr}/{args.genome_fasta} \\",
            f"  --gencode_gtf {rr}/{args.gencode_gtf} \\",
            f"  --pe_list {rr}/{args.pe_list} \\",
            f"  --longorf_tsv {rr}/{args.longorf_tsv} \\",
            f"  --leafcutter_repo {rr}/{args.leafcutter_repo} \\",
            f"  --star_threads 4",
        ]) + "\n"

        jid_s3 = _submit_sbatch_via_file(
            host, user, stage3_body, remote_run_dir, "stage3.sbatch", args.dry_run)
        chain_jids["stage3"] = jid_s3
        print(f"  stage3 JID: {jid_s3}")

    # Save chain info
    result = {
        "run_name": args.run_name,
        "run_id": run_id,
        "remote_run_dir": remote_run_dir,
        "reuse_stage1": reuse or None,
        "genome_fasta": args.genome_fasta,
        "gencode_gtf": args.gencode_gtf,
        "chain_jids": chain_jids,
        "dry_run": args.dry_run,
    }

    out_path = args.output_json or f"outputs/for_saba/_runs/{args.run_name}.json"
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nChain info saved to: {out_path}")
    print(f"Total stages in chain: {len(chain_jids)}")
    if not args.dry_run:
        print("Email notifications will arrive at completion/failure.")


if __name__ == "__main__":
    main()
