#!/usr/bin/env python3
"""
LeafCutter2 Web API — FastAPI backend for local and Quest/Slurm pipeline execution.
"""
from __future__ import annotations

import json
import os
import re
import signal
import shutil
import sqlite3
import subprocess
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

APP_ROOT = Path(__file__).resolve().parents[2]  # repo root
DATA_ROOT = APP_ROOT / "webapp" / "data"
JOBS_ROOT = DATA_ROOT / "jobs"
DB_PATH = DATA_ROOT / "jobs.db"
FRONTEND_DIR = APP_ROOT / "webapp" / "frontend"

GTEX_CACHE = DATA_ROOT / "gtex_cache"
GTEX_GCT_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
)
GTEX_GCT_FILENAME = "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
GTEX_ANNOT_URL = (
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/"
    "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
)
GTEX_ANNOT_FILENAME = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

# ---------------------------------------------------------------------------
# Thread-safe state
# ---------------------------------------------------------------------------

DB_LOCK = threading.Lock()
ACTIVE_WORKERS: Dict[str, threading.Thread] = {}
ACTIVE_PROCS: Dict[str, subprocess.Popen] = {}

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_paths() -> None:
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    JOBS_ROOT.mkdir(parents=True, exist_ok=True)
    GTEX_CACHE.mkdir(parents=True, exist_ok=True)


def download_if_missing(url: str, dest: Path) -> Path:
    """Download a file via curl/wget if it does not already exist locally."""
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("curl"):
        subprocess.run(
            ["curl", "-fSL", "--progress-bar", "-o", str(dest), url],
            check=True,
        )
    elif shutil.which("wget"):
        subprocess.run(
            ["wget", "-q", "--show-progress", "-O", str(dest), url],
            check=True,
        )
    else:
        import urllib.request
        urllib.request.urlretrieve(url, str(dest))
    return dest


def ensure_gtex_annotations() -> Path:
    """Download and cache the GTEx sample annotations; return local path."""
    return download_if_missing(GTEX_ANNOT_URL, GTEX_CACHE / GTEX_ANNOT_FILENAME)



def parse_gtex_tissues(annot_path: Path) -> List[Dict[str, Any]]:
    """Parse annotation TSV and return sorted list of {tissue, sample_count}.

    Only counts RNA-seq samples (SMAFRZE == "RNASEQ") so tissues that lack
    junction data (e.g. "Cells - Leukemia cell line (CML)" whose samples
    are all EXCLUDE) are filtered out automatically.
    """
    from collections import Counter

    counts: Counter = Counter()
    with open(annot_path, "r") as fh:
        header = fh.readline().strip().split("\t")
        try:
            tis_col = header.index("SMTSD")
            frz_col = header.index("SMAFRZE")
        except ValueError:
            return []
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) > max(tis_col, frz_col) and fields[tis_col]:
                if fields[frz_col] != "RNASEQ":
                    continue
                counts[fields[tis_col]] += 1

    return sorted(
        [{"tissue": t, "sample_count": c} for t, c in counts.items()],
        key=lambda x: x["tissue"],
    )


# ---------------------------------------------------------------------------
# Database helpers (SQLite)
# ---------------------------------------------------------------------------


def init_db() -> None:
    ensure_paths()
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                id              TEXT PRIMARY KEY,
                mode            TEXT NOT NULL,
                status          TEXT NOT NULL,
                created_at      TEXT NOT NULL,
                updated_at      TEXT NOT NULL,
                work_dir        TEXT NOT NULL,
                input_payload   TEXT NOT NULL,
                config_payload  TEXT NOT NULL,
                quest_job_id    TEXT,
                quest_account   TEXT,
                quest_partition TEXT,
                error           TEXT,
                summary_path    TEXT,
                artifacts_zip   TEXT,
                runner_ref      TEXT
            )
            """
        )
        for col, typedef in [
            ("progress_pct", "INTEGER DEFAULT 0"),
            ("progress_label", "TEXT DEFAULT ''"),
        ]:
            try:
                conn.execute(f"ALTER TABLE jobs ADD COLUMN {col} {typedef}")
            except sqlite3.OperationalError:
                pass
        conn.commit()


def db_execute(query: str, params: tuple = ()) -> None:
    with DB_LOCK:
        with sqlite3.connect(DB_PATH) as conn:
            conn.execute(query, params)
            conn.commit()


def db_fetch_one(query: str, params: tuple = ()) -> Optional[sqlite3.Row]:
    with DB_LOCK:
        with sqlite3.connect(DB_PATH) as conn:
            conn.row_factory = sqlite3.Row
            return conn.execute(query, params).fetchone()


def db_fetch_all(query: str, params: tuple = ()) -> List[sqlite3.Row]:
    with DB_LOCK:
        with sqlite3.connect(DB_PATH) as conn:
            conn.row_factory = sqlite3.Row
            return conn.execute(query, params).fetchall()


def update_job(job_id: str, **fields: Any) -> None:
    if not fields:
        return
    fields["updated_at"] = now_iso()
    keys = list(fields.keys())
    values = [fields[k] for k in keys]
    set_clause = ", ".join(f"{k}=?" for k in keys)
    db_execute(f"UPDATE jobs SET {set_clause} WHERE id=?", tuple(values + [job_id]))


def _row_to_dict(row: sqlite3.Row) -> Dict[str, Any]:
    d = dict(row)
    for k in ("input_payload", "config_payload"):
        if d.get(k):
            try:
                d[k] = json.loads(d[k])
            except (json.JSONDecodeError, TypeError):
                pass
    return d


def get_job_or_404(job_id: str) -> Dict[str, Any]:
    row = db_fetch_one("SELECT * FROM jobs WHERE id=?", (job_id,))
    if row is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return _row_to_dict(row)


# ---------------------------------------------------------------------------
# Shell / SSH helpers
# ---------------------------------------------------------------------------


def run_local_cmd(
    cmd: List[str], cwd: Optional[Path] = None
) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True, text=True, capture_output=True)


def validate_shell_script(path: Path) -> None:
    """Run ``bash -n`` on *path* and raise with diagnostics if syntax is invalid."""
    result = subprocess.run(
        ["bash", "-n", str(path)], capture_output=True, text=True,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or f"bash -n exited with code {result.returncode}"
        raise RuntimeError(
            f"Shell syntax error in {path.name}:\n{detail}"
        )


def _ssh_base_opts() -> List[str]:
    """Common SSH/SCP options: batch mode, auto-accept host keys, timeout, explicit key."""
    ssh_key = Path.home() / ".ssh" / "id_ed25519"
    opts = [
        "-o", "BatchMode=yes",
        "-o", "StrictHostKeyChecking=accept-new",
        "-o", "ConnectTimeout=30",
    ]
    if ssh_key.exists():
        opts += ["-i", str(ssh_key)]
    return opts


def run_ssh_cmd(host: str, user: str, remote_cmd: str) -> subprocess.CompletedProcess:
    cmd = ["ssh"] + _ssh_base_opts() + [f"{user}@{host}", remote_cmd]
    try:
        return run_local_cmd(cmd, cwd=APP_ROOT)
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip() if exc.stderr else f"exit status {exc.returncode}"
        raise subprocess.CalledProcessError(
            exc.returncode, exc.cmd, output=exc.stdout, stderr=exc.stderr
        ) from RuntimeError(f"SSH failed: {detail}")


def _fetch_remote_bytes(host: str, user: str, remote_path: str, byte_offset: int) -> str:
    """Fetch new content from a remote file starting at byte_offset (non-raising)."""
    cmd = (
        ["ssh"] + _ssh_base_opts()
        + [f"{user}@{host}", f"tail -c +{byte_offset + 1} {remote_path} 2>/dev/null || true"]
    )
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        return result.stdout or ""
    except Exception:
        return ""


def _sync_remote_log(
    host: str, user: str, remote_stdout: str, remote_stderr: str,
    local_log: Path, stdout_offset: int, stderr_offset: int,
) -> tuple:
    """Append new remote stdout/stderr to local pipeline.log. Returns updated offsets."""
    new_out = _fetch_remote_bytes(host, user, remote_stdout, stdout_offset)
    new_err = _fetch_remote_bytes(host, user, remote_stderr, stderr_offset)
    chunk = ""
    if new_out:
        chunk += new_out
    if new_err:
        for line in new_err.splitlines():
            if not line.startswith((" ", "\t")) and "......" not in line:
                chunk += line + "\n"
    if chunk:
        with open(local_log, "a") as f:
            f.write(chunk)
    return (
        stdout_offset + len(new_out.encode("utf-8", errors="replace")),
        stderr_offset + len(new_err.encode("utf-8", errors="replace")),
    )


def _list_remote_files(host: str, user: str, pattern: str) -> List[str]:
    """List concrete remote file paths matching a glob pattern (non-raising)."""
    cmd = (
        ["ssh"] + _ssh_base_opts()
        + [f"{user}@{host}", f"ls -1 {pattern} 2>/dev/null || true"]
    )
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        return [p for p in (result.stdout or "").splitlines() if p.strip()]
    except Exception:
        return []


def _sync_remote_log_multi(
    host: str, user: str, remote_patterns: List[str],
    local_log: Path, offsets: Dict[str, int],
) -> Dict[str, int]:
    """Incrementally sync multiple remote log files into a single local log.

    *offsets* maps concrete remote paths to the last-read byte offset.
    For patterns containing globs, files are listed first so each concrete
    path is tracked independently.  Returns the updated offsets dict.
    """
    concrete_paths: List[str] = []
    for pat in remote_patterns:
        if "*" in pat or "?" in pat:
            concrete_paths.extend(_list_remote_files(host, user, pat))
        else:
            concrete_paths.append(pat)

    chunk = ""
    for rpath in concrete_paths:
        offset = offsets.get(rpath, 0)
        new_data = _fetch_remote_bytes(host, user, rpath, offset)
        if not new_data:
            continue
        is_stderr = rpath.endswith(".err")
        if is_stderr:
            for line in new_data.splitlines():
                if not line.startswith((" ", "\t")) and "......" not in line:
                    chunk += line + "\n"
        else:
            chunk += new_data
        offsets[rpath] = offset + len(new_data.encode("utf-8", errors="replace"))

    if chunk:
        with open(local_log, "a") as f:
            f.write(chunk)
    return offsets


def run_scp_cmd(srcs: List[str], dest: str, recursive: bool = False) -> subprocess.CompletedProcess:
    cmd = ["scp"] + _ssh_base_opts()
    if recursive:
        cmd.append("-r")
    cmd += srcs + [dest]
    try:
        return run_local_cmd(cmd, cwd=APP_ROOT)
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip() if exc.stderr else f"exit status {exc.returncode}"
        raise subprocess.CalledProcessError(
            exc.returncode, exc.cmd, output=exc.stdout, stderr=exc.stderr
        ) from RuntimeError(f"SCP failed: {detail}")


def _scp_error_detail(exc: subprocess.CalledProcessError) -> str:
    """Extract a human-readable error from a failed SCP CalledProcessError."""
    parts = [f"exit status {exc.returncode}"]
    if exc.stderr and exc.stderr.strip():
        parts.append(exc.stderr.strip())
    return " — ".join(parts)


def _retrieve_outputs_with_retry(
    host: str, user: str, remote_out_dir: str, local_out: Path,
    max_retries: int = 3, backoff: float = 5.0,
) -> None:
    """SCP remote output directory to local, retrying on transient failures."""
    src = f"{user}@{host}:{remote_out_dir}/."
    last_exc: Optional[subprocess.CalledProcessError] = None
    for attempt in range(1, max_retries + 1):
        try:
            run_scp_cmd([src], str(local_out), recursive=True)
            return
        except subprocess.CalledProcessError as exc:
            last_exc = exc
            if attempt < max_retries:
                time.sleep(backoff * attempt)
    assert last_exc is not None
    detail = _scp_error_detail(last_exc)
    raise RuntimeError(
        f"Output retrieval failed after {max_retries} attempts: {detail}"
    )


# ---------------------------------------------------------------------------
# Pipeline command builder
# ---------------------------------------------------------------------------


def build_pipeline_cmd(
    run_dir: Path,
    star_sj_paths: List[Path],
    prefix: str,
    leafcutter_repo: str,
    leafcutter2_repo: str,
    genome_fasta: str,
    gencode_gtf: str,
    min_reads: int,
    max_intron_len: int,
    samples_tsv: Optional[Path] = None,
) -> List[str]:
    cmd = [
        "python3",
        str(APP_ROOT / "scripts" / "lc2_pipeline.py"),
        "--workdir", str(run_dir),
        "--prefix", prefix,
        "--leafcutter_repo", leafcutter_repo,
        "--leafcutter2_repo", leafcutter2_repo,
        "--genome_fasta", genome_fasta,
        "--gencode_gtf", gencode_gtf,
        "--min_reads", str(min_reads),
        "--max_intron_len", str(max_intron_len),
        "--star_sj",
    ]
    cmd.extend(str(p) for p in star_sj_paths)
    if samples_tsv is not None:
        cmd.extend(["--samples_tsv", str(samples_tsv)])
    return cmd


# ---------------------------------------------------------------------------
# Artifact packaging
# ---------------------------------------------------------------------------


def package_outputs(job_id: str, work_dir: Path) -> Optional[Path]:
    import zipfile

    out_dir = work_dir / "out"
    if not out_dir.exists():
        return None

    zip_path = work_dir / f"{job_id}_artifacts.zip"
    subdirs = ["out", "lc2", "clusters"]

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for subdir in subdirs:
            src = work_dir / subdir
            if not src.exists():
                continue
            for root, _dirs, files in os.walk(src):
                for fname in files:
                    full = Path(root) / fname
                    arcname = str(full.relative_to(work_dir))
                    zf.write(full, arcname)

    return zip_path


# ---------------------------------------------------------------------------
# Job submission helpers
# ---------------------------------------------------------------------------


def _insert_job_record(
    *,
    job_id: str,
    mode: str,
    work_dir: Path,
    input_payload: Dict[str, Any],
    config: Dict[str, Any],
    quest_account: Optional[str],
    quest_partition: Optional[str],
) -> None:
    db_execute(
        """
        INSERT INTO jobs
        (id, mode, status, created_at, updated_at, work_dir, input_payload, config_payload,
         quest_job_id, quest_account, quest_partition, error, summary_path, artifacts_zip, runner_ref)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, NULL, ?, ?, NULL, NULL, NULL, NULL)
        """,
        (
            job_id, mode, "queued", now_iso(), now_iso(), str(work_dir),
            json.dumps(input_payload), json.dumps(config),
            quest_account or None, quest_partition or None,
        ),
    )


def _submit_slurm_and_start_poller(
    *,
    job_id: str,
    work_dir: Path,
    remote_run_dir: str,
    sbatch_body: str,
    quest_host: str,
    quest_user: str,
) -> str:
    slurm_script = work_dir / "job.sbatch"
    slurm_script.write_text(sbatch_body)
    validate_shell_script(slurm_script)
    remote_script = f"{remote_run_dir}/job.sbatch"
    run_scp_cmd([str(slurm_script)], f"{quest_user}@{quest_host}:{remote_script}")

    submit = run_ssh_cmd(quest_host, quest_user, f"sbatch --parsable {remote_script}")
    quest_job_id = submit.stdout.strip().split(";")[0]
    if not quest_job_id:
        raise RuntimeError("sbatch returned no job ID")

    update_job(job_id, status="submitted", quest_job_id=quest_job_id)
    t = threading.Thread(
        target=slurm_worker,
        args=(job_id, quest_host, quest_user, remote_run_dir, f"{remote_run_dir}/out"),
        daemon=True,
    )
    ACTIVE_WORKERS[job_id] = t
    t.start()
    return quest_job_id


# ---------------------------------------------------------------------------
# Slurm state mapping
# ---------------------------------------------------------------------------

SLURM_MAP_SUBMITTED = {"PENDING", "CONFIGURING", "REQUEUED"}
SLURM_MAP_RUNNING = {"RUNNING", "COMPLETING", "SUSPENDED"}
SLURM_MAP_SUCCEEDED = {"COMPLETED"}
SLURM_MAP_CANCELLED = {"CANCELLED", "CANCELLED+"}
SLURM_MAP_FAILED = {"FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"}


def map_slurm_state(raw: str) -> str:
    state = (raw or "").upper().strip()
    if state in SLURM_MAP_SUBMITTED:
        return "submitted"
    if state in SLURM_MAP_RUNNING:
        return "running"
    if state in SLURM_MAP_SUCCEEDED:
        return "succeeded"
    if state in SLURM_MAP_CANCELLED or state.startswith("CANCELLED"):
        return "cancelled"
    if state in SLURM_MAP_FAILED:
        return "failed"
    return "submitted"


# ---------------------------------------------------------------------------
# Progress tracking helpers
# ---------------------------------------------------------------------------

_RE_GCT_ROWS = re.compile(r"\.\.\.\s+([\d,]+)/([\d,]+)\s+rows")
_RE_PIPELINE_CLUSTER = re.compile(r">>\s+.*leafcutter_cluster|Wrote \d+ clusters")
_RE_PIPELINE_LC2 = re.compile(r">>\s+.*leafcutter2|Refine clusters|Sorting .+\.bed|Extracting numerators")
_RE_PIPELINE_DONE = re.compile(r"^Done\.\s*$", re.MULTILINE)


def _read_log_tail(log_path: Path, tail_bytes: int = 8192) -> str:
    """Read the last tail_bytes of a log file."""
    if not log_path.exists():
        return ""
    try:
        size = log_path.stat().st_size
        with open(log_path, "rb") as f:
            f.seek(max(0, size - tail_bytes))
            return f.read().decode("utf-8", errors="replace")
    except Exception:
        return ""


def _gct_row_progress(log_tail: str, pct_min: int, pct_max: int) -> int:
    """Extract GCT conversion progress from row-count markers in log tail."""
    matches = list(_RE_GCT_ROWS.finditer(log_tail))
    if not matches:
        return pct_min
    last = matches[-1]
    current = int(last.group(1).replace(",", ""))
    total = int(last.group(2).replace(",", ""))
    if total <= 0:
        return pct_min
    frac = min(1.0, current / total)
    return pct_min + int(frac * (pct_max - pct_min))


def _upload_pipeline_progress(log_tail: str) -> tuple:
    """Progress for upload-source local pipeline (full 10-95 range)."""
    if _RE_PIPELINE_DONE.search(log_tail):
        return 95, "Pipeline complete"
    if _RE_PIPELINE_LC2.search(log_tail):
        return 70, "Classifying splicing events..."
    if _RE_PIPELINE_CLUSTER.search(log_tail):
        return 40, "Clustering junctions..."
    return 10, "Running pipeline..."


def _gtex_pipeline_progress(log_tail: str) -> tuple:
    """Progress for GTEx-source pipeline phase (60-92 range)."""
    if _RE_PIPELINE_DONE.search(log_tail):
        return 92, "Pipeline complete"
    if _RE_PIPELINE_LC2.search(log_tail):
        return 80, "Classifying splicing events..."
    if _RE_PIPELINE_CLUSTER.search(log_tail):
        return 65, "Clustering junctions..."
    return 60, "Running pipeline..."


def _slurm_running_progress(
    log_path: Path, source: str, pct_floor: int,
) -> tuple:
    """Derive finer progress from the pipeline log while Slurm state is RUNNING.

    Returns (pct, label).  *pct* is guaranteed >= *pct_floor* so progress
    never moves backwards.
    """
    tail = _read_log_tail(log_path)
    if not tail:
        return pct_floor, "Running on cluster"
    if source == "gtex":
        pct, label = _gtex_pipeline_progress(tail)
    else:
        pct, label = _upload_pipeline_progress(tail)
    return max(pct, pct_floor), label


def _poll_proc(
    job_id: str,
    proc: subprocess.Popen,
    log_path: Path,
    default_pct: int,
    default_label: str,
    progress_fn=None,
    interval: float = 3.0,
) -> int:
    """Poll subprocess until it exits, updating progress periodically.

    progress_fn(log_tail: str) -> (pct: int, label: str)
    Returns the subprocess exit code.
    """
    update_job(job_id, progress_pct=default_pct, progress_label=default_label)
    while proc.poll() is None:
        time.sleep(interval)
        if progress_fn:
            try:
                tail = _read_log_tail(log_path)
                pct, label = progress_fn(tail)
                update_job(job_id, progress_pct=pct, progress_label=label)
            except Exception:
                pass
    return proc.returncode


# ---------------------------------------------------------------------------
# Background workers
# ---------------------------------------------------------------------------


def local_worker(job_id: str, cmd: List[str], work_dir: Path) -> None:
    log_file = work_dir / "pipeline.log"
    update_job(job_id, status="running", progress_pct=5, progress_label="Starting pipeline...")
    try:
        with open(log_file, "w") as log:
            proc = subprocess.Popen(
                cmd, cwd=str(APP_ROOT), stdout=log, stderr=subprocess.STDOUT, text=True,
            )
            ACTIVE_PROCS[job_id] = proc
            update_job(job_id, runner_ref=str(proc.pid))
            rc = _poll_proc(
                job_id, proc, log_file, 10, "Running pipeline...",
                _upload_pipeline_progress,
            )

        if rc < 0:
            update_job(job_id, status="cancelled", progress_label="Terminated", error=f"Terminated by signal {-rc}")
            return
        if rc != 0:
            tail = ""
            try:
                lines = log_file.read_text().splitlines()
                tail = "\n".join(lines[-20:])
            except Exception:
                pass
            update_job(job_id, status="failed", progress_label="Failed", error=f"Exit code {rc}\n{tail}")
            return

        summary_path = work_dir / "out" / "summary.json"
        zip_path = package_outputs(job_id, work_dir)
        update_job(
            job_id,
            status="succeeded",
            progress_pct=100,
            progress_label="Complete",
            summary_path=str(summary_path) if summary_path.exists() else None,
            artifacts_zip=str(zip_path) if zip_path else None,
        )
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Error", error=f"{type(exc).__name__}: {exc}")
    finally:
        ACTIVE_PROCS.pop(job_id, None)
        ACTIVE_WORKERS.pop(job_id, None)


def slurm_worker(
    job_id: str,
    host: str,
    user: str,
    remote_run_dir: str,
    remote_out_dir: str,
) -> None:
    try:
        job = get_job_or_404(job_id)
        slurm_id = job["quest_job_id"]
        poll_interval = 15
        ssh_backoff = poll_interval
        max_backoff = 300
        update_job(job_id, progress_pct=5, progress_label="Submitted to cluster")

        inp = job.get("input_payload", {})
        source = inp.get("source", "upload")

        local_work = Path(job["work_dir"])
        local_log = local_work / "pipeline.log"
        remote_stdout = f"{remote_run_dir}/slurm_{slurm_id}.out"
        remote_stderr = f"{remote_run_dir}/slurm_{slurm_id}.err"
        out_offset = 0
        err_offset = 0
        last_pct = 5

        while True:
            ssh_ok = True
            try:
                sq = run_ssh_cmd(host, user, f"squeue -h -j {slurm_id} -o %T")
                sq_state = sq.stdout.strip()
            except subprocess.CalledProcessError:
                sq_state = ""
                ssh_ok = False

            out_offset, err_offset = _sync_remote_log(
                host, user, remote_stdout, remote_stderr,
                local_log, out_offset, err_offset,
            )

            if sq_state:
                ssh_backoff = poll_interval
                mapped = map_slurm_state(sq_state)
                if mapped == "running":
                    pct, label = _slurm_running_progress(local_log, source, last_pct)
                    last_pct = pct
                    update_job(job_id, status="running", progress_pct=pct, progress_label=label)
                else:
                    update_job(job_id, status=mapped, progress_pct=10, progress_label="Queued on cluster")
                time.sleep(poll_interval)
                continue

            try:
                sa = run_ssh_cmd(host, user, f"sacct -n -P -j {slurm_id} --format=State,ExitCode")
                first_line = sa.stdout.strip().splitlines()[0] if sa.stdout.strip() else ""
                final_state = first_line.split("|")[0].strip() if first_line else "UNKNOWN"
                ssh_ok = True
            except (subprocess.CalledProcessError, IndexError):
                final_state = "UNKNOWN"

            mapped = map_slurm_state(final_state)

            if mapped in ("running", "submitted"):
                ssh_backoff = poll_interval
                time.sleep(poll_interval)
                continue

            if not ssh_ok and mapped == "submitted":
                time.sleep(ssh_backoff)
                ssh_backoff = min(ssh_backoff * 2, max_backoff)
                continue

            # Final log sync
            _sync_remote_log(
                host, user, remote_stdout, remote_stderr,
                local_log, out_offset, err_offset,
            )

            if mapped == "cancelled":
                update_job(job_id, status="cancelled", progress_label="Cancelled", error=f"Slurm: {final_state}")
                return
            if mapped != "succeeded":
                # Fetch full stderr for error context
                full_err = _fetch_remote_bytes(host, user, remote_stderr, 0)
                err_tail = "\n".join(full_err.strip().splitlines()[-20:]) if full_err else ""
                update_job(job_id, status="failed", progress_label="Failed", error=f"Slurm: {final_state}\n{err_tail}")
                return

            # Succeeded — copy outputs back
            update_job(job_id, progress_pct=85, progress_label="Retrieving outputs...")
            local_work = Path(get_job_or_404(job_id)["work_dir"])
            local_out = local_work / "out"
            local_out.mkdir(parents=True, exist_ok=True)
            try:
                _retrieve_outputs_with_retry(host, user, remote_out_dir, local_out)
            except RuntimeError as e:
                update_job(job_id, status="failed", progress_label="Output retrieval failed", error=str(e))
                return

            summary_path = local_out / "summary.json"
            zip_path = package_outputs(job_id, local_work)
            update_job(
                job_id,
                status="succeeded",
                progress_pct=100,
                progress_label="Complete",
                summary_path=str(summary_path) if summary_path.exists() else None,
                artifacts_zip=str(zip_path) if zip_path else None,
            )
            return

    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Error", error=f"{type(exc).__name__}: {exc}")
    finally:
        ACTIVE_WORKERS.pop(job_id, None)


def slurm_array_worker(
    job_id: str,
    host: str,
    user: str,
    remote_run_dir: str,
    setup_jid: str,
    array_jid: str,
    merge_jid: str,
    n_tissues: int,
) -> None:
    """Poll a three-stage Slurm chain: setup -> array (per-tissue) -> merge."""
    try:
        poll_interval = 15
        ssh_backoff = poll_interval
        max_backoff = 300
        local_work = Path(get_job_or_404(job_id)["work_dir"])
        local_log = local_work / "pipeline.log"
        log_offsets: Dict[str, int] = {}

        stages = [
            (setup_jid, "Setup: downloading GTEx data", 5, 15),
            (array_jid, "Running per-tissue pipelines", 15, 80),
            (merge_jid, "Merging tissue results", 80, 92),
        ]

        for stage_idx, (jid, stage_label, pct_min, pct_max) in enumerate(stages):
            update_job(job_id, progress_pct=pct_min, progress_label=stage_label)
            is_array = (stage_idx == 1)

            log_patterns = [
                f"{remote_run_dir}/slurm_{jid}.out",
                f"{remote_run_dir}/slurm_{jid}.err",
            ]
            if is_array:
                log_patterns += [
                    f"{remote_run_dir}/slurm_{jid}_*.out",
                    f"{remote_run_dir}/slurm_{jid}_*.err",
                ]

            while True:
                ssh_ok = True
                try:
                    if is_array:
                        sq = run_ssh_cmd(host, user, f"squeue -h -j {jid} -o %T")
                        running_states = sq.stdout.strip().splitlines()
                    else:
                        sq = run_ssh_cmd(host, user, f"squeue -h -j {jid} -o %T")
                        running_states = [sq.stdout.strip()] if sq.stdout.strip() else []
                except subprocess.CalledProcessError:
                    running_states = []
                    ssh_ok = False

                log_offsets = _sync_remote_log_multi(
                    host, user, log_patterns, local_log, log_offsets,
                )

                if running_states and any(s for s in running_states):
                    ssh_backoff = poll_interval
                    n_running = len([s for s in running_states if s])
                    if is_array:
                        completed_approx = max(0, n_tissues - n_running)
                        frac = completed_approx / max(n_tissues, 1)
                        pct = pct_min + int(frac * (pct_max - pct_min))
                        label = f"{stage_label} ({completed_approx}/{n_tissues} done)"
                    else:
                        pct = pct_min + (pct_max - pct_min) // 2
                        label = stage_label

                    mapped = "running"
                    for s in running_states:
                        if s and map_slurm_state(s) in ("submitted",):
                            mapped = "submitted"
                            break

                    update_job(job_id, status="running" if mapped == "running" else "submitted",
                               progress_pct=pct, progress_label=label)
                    time.sleep(poll_interval)
                    continue

                # Stage left the queue — final log sync then check state
                log_offsets = _sync_remote_log_multi(
                    host, user, log_patterns, local_log, log_offsets,
                )

                try:
                    sa = run_ssh_cmd(host, user, f"sacct -n -P -j {jid} --format=State,ExitCode")
                    lines = [l for l in sa.stdout.strip().splitlines() if l.strip()]
                    ssh_ok = True
                except (subprocess.CalledProcessError, IndexError):
                    lines = []

                all_states = []
                for line in lines:
                    parts = line.split("|")
                    state = parts[0].strip() if parts else "UNKNOWN"
                    all_states.append(state)

                still_active = [s for s in all_states if map_slurm_state(s) in ("running", "submitted")]
                if still_active:
                    ssh_backoff = poll_interval
                    time.sleep(poll_interval)
                    continue

                failed_states = [s for s in all_states if map_slurm_state(s) == "failed"]
                cancelled_states = [s for s in all_states if map_slurm_state(s) == "cancelled"]

                if cancelled_states:
                    update_job(job_id, status="cancelled", progress_label="Cancelled",
                               error=f"Stage {stage_label} cancelled")
                    return

                if failed_states and not is_array:
                    err_patterns = [p for p in log_patterns if p.endswith(".err")]
                    err_offsets: Dict[str, int] = {}
                    _sync_remote_log_multi(host, user, err_patterns, local_log, err_offsets)
                    err_tail = _read_log_tail(local_log, 4096)
                    err_lines = err_tail.strip().splitlines()[-30:] if err_tail else []
                    update_job(
                        job_id, status="failed", progress_label=f"{stage_label} failed",
                        error=f"Stage {stage_label}: {len(failed_states)} task(s) failed\n" + "\n".join(err_lines),
                    )
                    return

                if not all_states:
                    if not ssh_ok:
                        time.sleep(ssh_backoff)
                        ssh_backoff = min(ssh_backoff * 2, max_backoff)
                    else:
                        time.sleep(poll_interval)
                    continue

                # Stage completed (possibly with some failed array tasks)
                ssh_backoff = poll_interval
                if failed_states and is_array:
                    n_failed = len([s for s in all_states if map_slurm_state(s) == "failed"])
                    update_job(job_id, progress_pct=pct_max,
                               progress_label=f"{stage_label} complete ({n_failed} task(s) had no data)")
                else:
                    update_job(job_id, progress_pct=pct_max, progress_label=f"{stage_label} complete")
                break

        # All stages done — retrieve outputs
        update_job(job_id, progress_pct=92, progress_label="Retrieving outputs...")
        local_work = Path(get_job_or_404(job_id)["work_dir"])
        local_out = local_work / "out"
        local_out.mkdir(parents=True, exist_ok=True)
        remote_out_dir = f"{remote_run_dir}/out"
        try:
            _retrieve_outputs_with_retry(host, user, remote_out_dir, local_out)
        except RuntimeError as e:
            update_job(job_id, status="failed", progress_label="Output retrieval failed",
                       error=str(e))
            return

        summary_path = local_out / "summary.json"
        zip_path = package_outputs(job_id, local_work)
        update_job(
            job_id,
            status="succeeded",
            progress_pct=100,
            progress_label="Complete",
            summary_path=str(summary_path) if summary_path.exists() else None,
            artifacts_zip=str(zip_path) if zip_path else None,
        )

    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Error",
                   error=f"{type(exc).__name__}: {exc}")
    finally:
        ACTIVE_WORKERS.pop(job_id, None)


def gtex_local_worker(job_id: str, tissues: str, work_dir: Path, config: Dict[str, Any]) -> None:
    """Background worker: download GTEx GCT, convert to BED, run pipeline."""
    log_file = work_dir / "pipeline.log"
    update_job(job_id, status="running", progress_pct=2, progress_label="Preparing...")
    try:
        with open(log_file, "w") as log:
            # Phase 1: ensure annotations + GCT are cached
            log.write("[GTEx] Ensuring annotations are cached ...\n")
            log.flush()
            update_job(job_id, progress_pct=5, progress_label="Downloading annotations...")

            def _dl(url: str, dest: Path, label: str) -> bool:
                if dest.exists():
                    return True
                log.write(f"[GTEx] Downloading {label} from {url}\n")
                log.flush()
                if shutil.which("curl"):
                    cmd = ["curl", "-fSL", "--progress-bar", "-o", str(dest), url]
                elif shutil.which("wget"):
                    cmd = ["wget", "-q", "--show-progress", "-O", str(dest), url]
                else:
                    log.write("[GTEx] ERROR: neither curl nor wget found\n")
                    return False
                rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True).returncode
                return rc == 0

            annot_path = GTEX_CACHE / GTEX_ANNOT_FILENAME
            if not _dl(GTEX_ANNOT_URL, annot_path, "annotations"):
                update_job(job_id, status="failed", progress_label="Download failed", error="Annotations download failed")
                return

            update_job(job_id, progress_pct=10, progress_label="Downloading junction data...")
            gct_path = GTEX_CACHE / GTEX_GCT_FILENAME
            if not gct_path.exists():
                log.write("[GTEx] This will take a while on the first run (~4 GB) ...\n")
                log.flush()
            if not _dl(GTEX_GCT_URL, gct_path, "junction GCT"):
                gct_path.unlink(missing_ok=True)
                update_job(job_id, status="failed", progress_label="Download failed", error="GCT download failed")
                return
            else:
                log.write("[GTEx] Using cached GCT\n")
                log.flush()

            update_job(job_id, progress_pct=20, progress_label="Data ready")

            # Phase 2: convert GCT -> per-sample BEDs
            bed_dir = work_dir / "junctions_bed"
            log.write(f"[GTEx] Converting GCT for tissues: {tissues}\n")
            log.flush()
            convert_cmd = [
                "python3",
                str(APP_ROOT / "scripts" / "gtex_gct_to_bed.py"),
                "--gct", str(gct_path),
                "--annotations", str(annot_path),
                "--tissues", tissues,
                "--outdir", str(bed_dir),
                "--max_samples_per_tissue", str(config.get("max_samples_per_tissue", 0)),
            ]
            proc = subprocess.Popen(
                convert_cmd, cwd=str(APP_ROOT), stdout=log, stderr=subprocess.STDOUT, text=True,
            )
            ACTIVE_PROCS[job_id] = proc
            rc = _poll_proc(
                job_id, proc, log_file, 25, "Converting GCT to BED...",
                lambda tail: (_gct_row_progress(tail, 25, 55), "Converting GCT to BED..."),
            )
            if rc != 0:
                update_job(job_id, status="failed", progress_label="Conversion failed", error=f"GCT-to-BED conversion failed (exit {rc})")
                return

            filelist = bed_dir / "junction_files.txt"
            if not filelist.exists():
                update_job(job_id, status="failed", progress_label="Conversion failed", error="Conversion produced no junction_files.txt")
                return

            bed_files = [l.strip() for l in filelist.read_text().splitlines() if l.strip()]
            log.write(f"[GTEx] Conversion done: {len(bed_files)} BED files\n")
            log.flush()
            update_job(job_id, progress_pct=55, progress_label="Conversion complete")

            # Phase 3: run the LeafCutter2 pipeline using --junction_beds
            pipeline_cmd = [
                "python3",
                str(APP_ROOT / "scripts" / "lc2_pipeline.py"),
                "--workdir", str(work_dir),
                "--prefix", config.get("prefix", "gtex_run"),
                "--leafcutter_repo", config.get("leafcutter_repo", "tools/leafcutter"),
                "--leafcutter2_repo", config.get("leafcutter2_repo", "tools/leafcutter2"),
                "--genome_fasta", config.get("genome_fasta", "refs/GRCh38.fa"),
                "--gencode_gtf", config.get("gencode_gtf", "refs/gencode.v46.annotation.gtf"),
                "--min_reads", str(config.get("min_reads", 50)),
                "--max_intron_len", str(config.get("max_intron_len", 500000)),
                "--junction_beds",
            ] + bed_files

            log.write("[Pipeline] Starting LeafCutter2 pipeline ...\n")
            log.flush()
            proc = subprocess.Popen(
                pipeline_cmd, cwd=str(APP_ROOT), stdout=log, stderr=subprocess.STDOUT, text=True,
            )
            ACTIVE_PROCS[job_id] = proc
            update_job(job_id, runner_ref=str(proc.pid))
            rc = _poll_proc(
                job_id, proc, log_file, 60, "Starting pipeline...",
                _gtex_pipeline_progress,
            )

        if rc != 0:
            tail = ""
            try:
                lines = log_file.read_text().splitlines()
                tail = "\n".join(lines[-20:])
            except Exception:
                pass
            update_job(job_id, status="failed", progress_label="Pipeline failed", error=f"Pipeline failed (exit {rc})\n{tail}")
            return

        update_job(job_id, progress_pct=95, progress_label="Packaging outputs...")
        summary_path = work_dir / "out" / "summary.json"
        zip_path = package_outputs(job_id, work_dir)
        update_job(
            job_id,
            status="succeeded",
            progress_pct=100,
            progress_label="Complete",
            summary_path=str(summary_path) if summary_path.exists() else None,
            artifacts_zip=str(zip_path) if zip_path else None,
        )
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Error", error=f"{type(exc).__name__}: {exc}")
    finally:
        ACTIVE_PROCS.pop(job_id, None)
        ACTIVE_WORKERS.pop(job_id, None)


# ---------------------------------------------------------------------------
# File upload helper
# ---------------------------------------------------------------------------


def write_upload(upload: UploadFile, target: Path) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "wb") as fh:
        while True:
            chunk = upload.file.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    return target


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------

app = FastAPI(title="LeafCutter2 Pipeline", version="0.2.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def startup() -> None:
    init_db()
    _recover_orphaned_jobs()


def _recover_orphaned_jobs() -> None:
    """Re-attach polling threads to slurm jobs orphaned by a server restart."""
    rows = db_fetch_all(
        "SELECT * FROM jobs WHERE status IN ('queued','submitted','running')"
    )
    for row in rows:
        job = _row_to_dict(row)
        job_id = job["id"]

        if job["mode"] == "slurm" and job.get("quest_job_id"):
            cfg = job.get("config_payload", {})
            host = cfg.get("quest_host", "login.quest.northwestern.edu")
            user = cfg.get("quest_user")
            remote_work_root = cfg.get("remote_work_root", "")
            if not remote_work_root.strip():
                rr = cfg.get("remote_repo_root", "").strip()
                if rr:
                    remote_work_root = f"{rr.rstrip('/')}/jobs"
            if not user or not remote_work_root:
                update_job(job_id, status="failed", progress_label="Recovery failed",
                           error="Missing quest_user or remote_work_root — cannot re-attach poller")
                continue
            remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"

            # Check if this is a parallel array job
            if cfg.get("slurm_array_jid"):
                inp = job.get("input_payload", {})
                n_tissues = len(inp.get("gtex_tissues", []))
                t = threading.Thread(
                    target=slurm_array_worker,
                    args=(job_id, host, user, remote_run_dir,
                          cfg["slurm_setup_jid"], cfg["slurm_array_jid"],
                          cfg["slurm_merge_jid"], n_tissues),
                    daemon=True,
                )
            else:
                remote_out_dir = f"{remote_run_dir}/out"
                t = threading.Thread(
                    target=slurm_worker,
                    args=(job_id, host, user, remote_run_dir, remote_out_dir),
                    daemon=True,
                )
            ACTIVE_WORKERS[job_id] = t
            t.start()
        else:
            update_job(job_id, status="failed", progress_label="Interrupted",
                       error="Server restarted while job was active")


# Serve frontend static files
app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR)), name="static")


@app.get("/")
def root() -> FileResponse:
    index = FRONTEND_DIR / "index.html"
    if not index.exists():
        raise HTTPException(status_code=404, detail="Frontend not found")
    return FileResponse(str(index), media_type="text/html")


# ---------------------------------------------------------------------------
# GET /gtex/tissues — list available GTEx tissues
# ---------------------------------------------------------------------------


@app.get("/gtex/tissues")
def gtex_tissues() -> JSONResponse:
    try:
        annot = ensure_gtex_annotations()
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Failed to fetch GTEx annotations: {exc}")
    return JSONResponse(parse_gtex_tissues(annot))


# ---------------------------------------------------------------------------
# POST /jobs/zebrafish/stage1 — submit alignment stage
# ---------------------------------------------------------------------------


@app.post("/jobs/zebrafish/qc")
async def create_zebrafish_qc_job(
    mode: str = Form("slurm"),
    zebrafish_projects_csv: str = Form(""),
    prefix: str = Form("zebrafish_run"),
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    zebrafish_max_runs: int = Form(0),
    zebrafish_star_threads: int = Form(8),
) -> JSONResponse:
    mode = mode.strip().lower()
    if mode != "slurm":
        raise HTTPException(status_code=400, detail="Zebrafish QC currently supports slurm mode only.")

    projects = [p.strip() for p in zebrafish_projects_csv.split(",") if p.strip()]
    if not projects:
        raise HTTPException(status_code=400, detail="Provide at least one project ID (PRJNA/ERP/DRP).")

    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        raise HTTPException(status_code=400, detail=f"Slurm mode requires: {', '.join(missing)}")

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "prefix": prefix,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "zebrafish_max_runs": zebrafish_max_runs,
        "zebrafish_star_threads": zebrafish_star_threads,
        "zebrafish_stage": "qc_report",
    }
    input_payload = {
        "source": "zebrafish",
        "zebrafish_stage": "qc_report",
        "zebrafish_projects": projects,
    }
    _insert_job_record(
        job_id=job_id,
        mode="slurm",
        work_dir=work_dir,
        input_payload=input_payload,
        config=config,
        quest_account=quest_account,
        quest_partition=quest_partition,
    )

    try:
        if not remote_work_root.strip():
            remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
        if remote_work_root.startswith("~"):
            home_result = run_ssh_cmd(quest_host, quest_user, "echo $HOME")
            remote_home = home_result.stdout.strip()
            remote_work_root = remote_work_root.replace("~", remote_home, 1)
        config["remote_work_root"] = remote_work_root
        db_execute("UPDATE jobs SET config_payload=? WHERE id=?", (json.dumps(config), job_id))

        remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
        run_ssh_cmd(quest_host, quest_user, f"mkdir -p {remote_run_dir}/inputs")

        rr = remote_repo_root.rstrip("/")
        slurm_time = {"short": "04:00:00", "normal": "48:00:00", "long": "168:00:00"}.get(
            quest_partition, "48:00:00"
        )
        sbatch_preamble = [
            "#!/bin/bash",
            f"#SBATCH --account={quest_account}",
            f"#SBATCH --partition={quest_partition}",
            f"#SBATCH --time={slurm_time}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=32G",
            f"#SBATCH --job-name=zf_qc_{job_id[:8]}",
            f"#SBATCH --output={remote_run_dir}/slurm_%j.out",
            f"#SBATCH --error={remote_run_dir}/slurm_%j.err",
            f"#SBATCH --cpus-per-task={min(4, int(zebrafish_star_threads))}",
            "",
            "set -eo pipefail",
            f"cd {rr}",
            "module purge",
            "module load anaconda3/2022.05",
            "module load fastp/0.23.4",
            "module load fastqc/0.12.0",
            "module load multiqc/1.11",
            "",
        ]
        zf_projects = ",".join(projects)
        zf_manifest = f"{remote_run_dir}/inputs/zebrafish_runs.tsv"
        run_ssh_cmd(
            quest_host,
            quest_user,
            (
                f"python3 {rr}/scripts/zebrafish_project_runs.py"
                f" --project_ids \"{zf_projects}\""
                f" --out \"{zf_manifest}\""
                f" --max_runs {int(zebrafish_max_runs)}"
            ),
        )
        sbatch_body = "\n".join(sbatch_preamble + [
            (
                f"python3 {rr}/scripts/zebrafish_quest_pipeline.py"
                f" --stage qc"
                f" --workdir {remote_run_dir}"
                f" --manifest \"{zf_manifest}\""
                f" --max_runs {int(zebrafish_max_runs)}"
                f" --star_threads {min(4, int(zebrafish_star_threads))}"
            ),
        ]) + "\n"
        quest_job_id = _submit_slurm_and_start_poller(
            job_id=job_id,
            work_dir=work_dir,
            remote_run_dir=remote_run_dir,
            sbatch_body=sbatch_body,
            quest_host=quest_host,
            quest_user=quest_user,
        )
        return JSONResponse({
            "job_id": job_id,
            "status": "submitted",
            "mode": "slurm",
            "source": "zebrafish",
            "stage": "qc_report",
            "quest_job_id": quest_job_id,
        })
    except HTTPException:
        raise
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"SSH/SCP error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"QC submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"QC submission failed: {exc}")


# ---------------------------------------------------------------------------
# POST /jobs/zebrafish/stage1 — submit alignment stage (step 2)
# ---------------------------------------------------------------------------


@app.post("/jobs/zebrafish/stage1")
async def create_zebrafish_stage1_job(
    mode: str = Form("slurm"),
    zebrafish_projects_csv: str = Form(""),
    qc_job_id: str = Form(""),
    prefix: str = Form("zebrafish_run"),
    # Quest/Slurm fields
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    # Zebrafish stage1 params
    zebrafish_max_runs: int = Form(0),
    zebrafish_star_threads: int = Form(8),
    zebrafish_genome_fasta: str = Form("refs/GRCz12tu.fa"),
    zebrafish_gencode_gtf: str = Form("refs/GRCz12tu.ucsc.gtf"),
) -> JSONResponse:
    mode = mode.strip().lower()
    if mode != "slurm":
        raise HTTPException(status_code=400, detail="Zebrafish stage1 currently supports slurm mode only.")

    qc_job_id = qc_job_id.strip()
    use_qc = bool(qc_job_id)

    if use_qc:
        qc_job = get_job_or_404(qc_job_id)
        qc_input = qc_job.get("input_payload", {}) or {}
        if qc_input.get("source") != "zebrafish" or qc_input.get("zebrafish_stage") != "qc_report":
            raise HTTPException(status_code=400, detail="Provided qc_job_id is not a zebrafish QC job.")
        if qc_job.get("status") != "succeeded":
            raise HTTPException(status_code=400, detail="QC job must have succeeded before running alignment.")
    else:
        projects = [p.strip() for p in zebrafish_projects_csv.split(",") if p.strip()]
        if not projects:
            raise HTTPException(status_code=400, detail="Provide a qc_job_id or at least one project ID.")

    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        raise HTTPException(status_code=400, detail=f"Slurm mode requires: {', '.join(missing)}")

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "prefix": prefix,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "zebrafish_max_runs": zebrafish_max_runs,
        "zebrafish_star_threads": zebrafish_star_threads,
        "zebrafish_genome_fasta": zebrafish_genome_fasta,
        "zebrafish_gencode_gtf": zebrafish_gencode_gtf,
        "zebrafish_stage": "stage1_alignment",
    }
    if use_qc:
        config["qc_job_id"] = qc_job_id
    input_payload = {
        "source": "zebrafish",
        "zebrafish_stage": "stage1_alignment",
    }
    if use_qc:
        input_payload["qc_job_id"] = qc_job_id
    else:
        input_payload["zebrafish_projects"] = projects
    _insert_job_record(
        job_id=job_id,
        mode="slurm",
        work_dir=work_dir,
        input_payload=input_payload,
        config=config,
        quest_account=quest_account,
        quest_partition=quest_partition,
    )

    try:
        if not remote_work_root.strip():
            remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
        if remote_work_root.startswith("~"):
            home_result = run_ssh_cmd(quest_host, quest_user, "echo $HOME")
            remote_home = home_result.stdout.strip()
            remote_work_root = remote_work_root.replace("~", remote_home, 1)
        config["remote_work_root"] = remote_work_root
        db_execute("UPDATE jobs SET config_payload=? WHERE id=?", (json.dumps(config), job_id))

        remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
        run_ssh_cmd(quest_host, quest_user, f"mkdir -p {remote_run_dir}/inputs")

        rr = remote_repo_root.rstrip("/")
        slurm_time = {"short": "04:00:00", "normal": "48:00:00", "long": "168:00:00"}.get(
            quest_partition, "48:00:00"
        )
        sbatch_preamble = [
            "#!/bin/bash",
            f"#SBATCH --account={quest_account}",
            f"#SBATCH --partition={quest_partition}",
            f"#SBATCH --time={slurm_time}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=180G",
            f"#SBATCH --job-name=zf_align_{job_id[:8]}",
            f"#SBATCH --output={remote_run_dir}/slurm_%j.out",
            f"#SBATCH --error={remote_run_dir}/slurm_%j.err",
            f"#SBATCH --cpus-per-task={int(zebrafish_star_threads)}",
            "",
            "set -eo pipefail",
            f"cd {rr}",
            "module purge",
            "module load anaconda3/2022.05",
            "module load STAR/2.7.11b",
            "module load samtools",
            "",
        ]
        shared_star_index = f"{rr}/refs/star_index_grcz12tu"

        if use_qc:
            qc_remote_dir = f"{remote_work_root.rstrip('/')}/{qc_job_id}"
            run_ssh_cmd(
                quest_host, quest_user,
                f"[ -f \"{qc_remote_dir}/out/qc_report.json\" ] || (echo \"missing QC artifacts\" && exit 1)",
            )
            pipeline_cmd = (
                f"python3 {rr}/scripts/zebrafish_quest_pipeline.py"
                f" --stage stage1"
                f" --workdir {remote_run_dir}"
                f" --qc_workdir {qc_remote_dir}"
                f" --prefix {prefix}"
                f" --genome_fasta {rr}/{zebrafish_genome_fasta}"
                f" --gencode_gtf {rr}/{zebrafish_gencode_gtf}"
                f" --star_threads {int(zebrafish_star_threads)}"
                f" --star_index {shared_star_index}"
            )
        else:
            sbatch_preamble.insert(-1, "module load fastp/0.23.4")
            zf_projects = ",".join(projects)
            zf_manifest = f"{remote_run_dir}/inputs/zebrafish_runs.tsv"
            run_ssh_cmd(
                quest_host, quest_user,
                (
                    f"python3 {rr}/scripts/zebrafish_project_runs.py"
                    f" --project_ids \"{zf_projects}\""
                    f" --out \"{zf_manifest}\""
                    f" --max_runs {int(zebrafish_max_runs)}"
                ),
            )
            pipeline_cmd = (
                f"python3 {rr}/scripts/zebrafish_quest_pipeline.py"
                f" --stage stage1"
                f" --workdir {remote_run_dir}"
                f" --manifest \"{zf_manifest}\""
                f" --prefix {prefix}"
                f" --genome_fasta {rr}/{zebrafish_genome_fasta}"
                f" --gencode_gtf {rr}/{zebrafish_gencode_gtf}"
                f" --star_threads {int(zebrafish_star_threads)}"
                f" --max_runs {int(zebrafish_max_runs)}"
                f" --star_index {shared_star_index}"
            )

        sbatch_body = "\n".join(sbatch_preamble + [pipeline_cmd]) + "\n"
        quest_job_id = _submit_slurm_and_start_poller(
            job_id=job_id,
            work_dir=work_dir,
            remote_run_dir=remote_run_dir,
            sbatch_body=sbatch_body,
            quest_host=quest_host,
            quest_user=quest_user,
        )
        return JSONResponse({
            "job_id": job_id,
            "status": "submitted",
            "mode": "slurm",
            "source": "zebrafish",
            "stage": "stage1_alignment",
            "quest_job_id": quest_job_id,
        })
    except HTTPException:
        raise
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"SSH/SCP error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"Stage1 submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"Stage1 submission failed: {exc}")


# ---------------------------------------------------------------------------
# POST /jobs/zebrafish/stage2 — submit LC2 stage from prior stage1
# ---------------------------------------------------------------------------


@app.post("/jobs/zebrafish/stage2")
async def create_zebrafish_stage2_job(
    mode: str = Form("slurm"),
    stage1_job_id: str = Form(""),
    selected_stages_csv: str = Form(""),
    prefix: str = Form("zebrafish_run"),
    min_reads: int = Form(50),
    max_intron_len: int = Form(500000),
    leafcutter_repo: str = Form("tools/leafcutter"),
    leafcutter2_repo: str = Form("tools/leafcutter2"),
    # Quest/Slurm fields
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    # Zebrafish stage2 params
    zebrafish_genome_fasta: str = Form("refs/GRCz12tu.fa"),
    zebrafish_gencode_gtf: str = Form("refs/GRCz12tu.ucsc.gtf"),
    zebrafish_gene_type_tag: str = Form("gene_biotype"),
    zebrafish_transcript_type_tag: str = Form("transcript_biotype"),
    zebrafish_gene_name_tag: str = Form("gene_id"),
    zebrafish_transcript_name_tag: str = Form("transcript_id"),
) -> JSONResponse:
    mode = mode.strip().lower()
    if mode != "slurm":
        raise HTTPException(status_code=400, detail="Zebrafish stage2 currently supports slurm mode only.")
    stage1_job_id = stage1_job_id.strip()
    if not stage1_job_id:
        raise HTTPException(status_code=400, detail="Provide stage1_job_id.")

    stage1 = get_job_or_404(stage1_job_id)
    stage1_input = stage1.get("input_payload", {}) or {}
    if stage1_input.get("source") != "zebrafish" or stage1_input.get("zebrafish_stage") != "stage1_alignment":
        raise HTTPException(status_code=400, detail="Provided stage1_job_id is not a zebrafish stage1 job.")
    if stage1.get("status") != "succeeded":
        raise HTTPException(status_code=400, detail="Stage1 job must be succeeded before running stage2.")

    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        raise HTTPException(status_code=400, detail=f"Slurm mode requires: {', '.join(missing)}")

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "prefix": prefix,
        "min_reads": min_reads,
        "max_intron_len": max_intron_len,
        "leafcutter_repo": leafcutter_repo,
        "leafcutter2_repo": leafcutter2_repo,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "zebrafish_genome_fasta": zebrafish_genome_fasta,
        "zebrafish_gencode_gtf": zebrafish_gencode_gtf,
        "zebrafish_gene_type_tag": zebrafish_gene_type_tag,
        "zebrafish_transcript_type_tag": zebrafish_transcript_type_tag,
        "zebrafish_gene_name_tag": zebrafish_gene_name_tag,
        "zebrafish_transcript_name_tag": zebrafish_transcript_name_tag,
        "zebrafish_stage": "stage2_analysis",
        "stage1_job_id": stage1_job_id,
        "selected_stages": [t.strip() for t in selected_stages_csv.split(",") if t.strip()],
    }
    input_payload = {
        "source": "zebrafish",
        "zebrafish_stage": "stage2_analysis",
        "stage1_job_id": stage1_job_id,
        "selected_stages": [t.strip() for t in selected_stages_csv.split(",") if t.strip()],
    }
    _insert_job_record(
        job_id=job_id,
        mode="slurm",
        work_dir=work_dir,
        input_payload=input_payload,
        config=config,
        quest_account=quest_account,
        quest_partition=quest_partition,
    )

    try:
        if not remote_work_root.strip():
            remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
        if remote_work_root.startswith("~"):
            home_result = run_ssh_cmd(quest_host, quest_user, "echo $HOME")
            remote_home = home_result.stdout.strip()
            remote_work_root = remote_work_root.replace("~", remote_home, 1)
        config["remote_work_root"] = remote_work_root
        db_execute("UPDATE jobs SET config_payload=? WHERE id=?", (json.dumps(config), job_id))

        remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
        stage1_remote_dir = f"{remote_work_root.rstrip('/')}/{stage1_job_id}"
        run_ssh_cmd(quest_host, quest_user, f"mkdir -p {remote_run_dir}/inputs")
        run_ssh_cmd(
            quest_host,
            quest_user,
            f"[ -f \"{stage1_remote_dir}/out/stage1_alignment.json\" ] || (echo \"missing stage1 artifacts\" && exit 1)",
        )

        rr = remote_repo_root.rstrip("/")
        slurm_time = {"short": "04:00:00", "normal": "48:00:00", "long": "168:00:00"}.get(
            quest_partition, "48:00:00"
        )
        sbatch_preamble = [
            "#!/bin/bash",
            f"#SBATCH --account={quest_account}",
            f"#SBATCH --partition={quest_partition}",
            f"#SBATCH --time={slurm_time}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=64G",
            f"#SBATCH --job-name=zf_s2_{job_id[:8]}",
            f"#SBATCH --output={remote_run_dir}/slurm_%j.out",
            f"#SBATCH --error={remote_run_dir}/slurm_%j.err",
            "",
            "set -eo pipefail",
            f"cd {rr}",
            "module purge",
            "module load anaconda3/2022.05",
            "",
        ]
        sbatch_body = "\n".join(sbatch_preamble + [
            (
                f"python3 {rr}/scripts/zebrafish_quest_pipeline.py"
                f" --stage stage2"
                f" --workdir {remote_run_dir}"
                f" --stage1_workdir {stage1_remote_dir}"
                f" --prefix {prefix}"
                f" --leafcutter_repo {rr}/{leafcutter_repo}"
                f" --leafcutter2_repo {rr}/{leafcutter2_repo}"
                f" --genome_fasta {rr}/{zebrafish_genome_fasta}"
                f" --gencode_gtf {rr}/{zebrafish_gencode_gtf}"
                f" --min_reads {min_reads}"
                f" --max_intron_len {max_intron_len}"
                f" --gene_type_tag {zebrafish_gene_type_tag}"
                f" --transcript_type_tag {zebrafish_transcript_type_tag}"
                f" --gene_name_tag {zebrafish_gene_name_tag}"
                f" --transcript_name_tag {zebrafish_transcript_name_tag}"
                f" --selected_stages \"{selected_stages_csv}\""
            ),
        ]) + "\n"
        quest_job_id = _submit_slurm_and_start_poller(
            job_id=job_id,
            work_dir=work_dir,
            remote_run_dir=remote_run_dir,
            sbatch_body=sbatch_body,
            quest_host=quest_host,
            quest_user=quest_user,
        )
        return JSONResponse({
            "job_id": job_id,
            "status": "submitted",
            "mode": "slurm",
            "source": "zebrafish",
            "stage": "stage2_analysis",
            "stage1_job_id": stage1_job_id,
            "quest_job_id": quest_job_id,
        })
    except HTTPException:
        raise
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"SSH/SCP error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"Stage2 submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"Stage2 submission failed: {exc}")


# ---------------------------------------------------------------------------
# POST /jobs/zebrafish/full — run all stages automatically (dependency chain)
# ---------------------------------------------------------------------------


@app.post("/jobs/zebrafish/full")
async def create_zebrafish_full_job(
    mode: str = Form("slurm"),
    zebrafish_projects_csv: str = Form("PRJEB7244"),
    prefix: str = Form("zebrafish_run"),
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    zebrafish_max_runs: int = Form(0),
    zebrafish_star_threads: int = Form(8),
    zebrafish_genome_fasta: str = Form("refs/GRCz12tu.fa"),
    zebrafish_gencode_gtf: str = Form("refs/GRCz12tu.ucsc.gtf"),
    zebrafish_star_index: str = Form("refs/star_index_grcz12tu"),
    zebrafish_pe_list: str = Form("UPDATED_STRAND_POISEN_HITS.tsv"),
    zebrafish_longorf_tsv: str = Form("LongORF_PTC+_fromClustered.tsv"),
    reuse_stage1_workdir: str = Form(""),
    array_concurrency: int = Form(50),
    gene_type_tag: str = Form("gene_type"),
    transcript_type_tag: str = Form("transcript_type"),
    gene_name_tag: str = Form("gene_name"),
    transcript_name_tag: str = Form("transcript_name"),
) -> JSONResponse:
    """Submit the complete zebrafish pipeline as a Slurm dependency chain."""
    mode = mode.strip().lower()
    if mode != "slurm":
        raise HTTPException(status_code=400, detail="Full pipeline currently supports slurm mode only.")

    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        raise HTTPException(status_code=400, detail=f"Slurm mode requires: {', '.join(missing)}")

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "prefix": prefix,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "zebrafish_stage": "full_pipeline",
        "zebrafish_max_runs": zebrafish_max_runs,
        "zebrafish_star_threads": zebrafish_star_threads,
        "zebrafish_genome_fasta": zebrafish_genome_fasta,
        "zebrafish_gencode_gtf": zebrafish_gencode_gtf,
        "zebrafish_star_index": zebrafish_star_index,
        "reuse_stage1_workdir": reuse_stage1_workdir,
    }
    input_payload = {
        "source": "zebrafish",
        "zebrafish_stage": "full_pipeline",
        "zebrafish_projects": [p.strip() for p in zebrafish_projects_csv.split(",") if p.strip()],
    }
    _insert_job_record(
        job_id=job_id,
        mode="slurm",
        work_dir=work_dir,
        input_payload=input_payload,
        config=config,
        quest_account=quest_account,
        quest_partition=quest_partition,
    )

    try:
        if not remote_work_root.strip():
            remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
        config["remote_work_root"] = remote_work_root

        # Build the orchestrator command and run it in a background thread
        rr = remote_repo_root.rstrip("/")
        orchestrator_args = [
            "python3", str(Path(__file__).resolve().parent.parent / "scripts" / "run_zebrafish_full_chain.py"),
            "--run_name", f"webapp_{job_id[:8]}",
            "--quest_host", quest_host,
            "--quest_user", quest_user,
            "--quest_account", quest_account,
            "--quest_partition", quest_partition,
            "--remote_repo_root", remote_repo_root,
            "--projects", zebrafish_projects_csv,
            "--max_runs", str(zebrafish_max_runs),
            "--array_concurrency", str(array_concurrency),
            "--genome_fasta", zebrafish_genome_fasta,
            "--gencode_gtf", zebrafish_gencode_gtf,
            "--star_index", zebrafish_star_index,
            "--star_threads", str(zebrafish_star_threads),
            "--prefix", prefix,
            "--pe_list", zebrafish_pe_list,
            "--longorf_tsv", zebrafish_longorf_tsv,
            "--gene_type_tag", gene_type_tag,
            "--transcript_type_tag", transcript_type_tag,
            "--gene_name_tag", gene_name_tag,
            "--transcript_name_tag", transcript_name_tag,
            "--output_json", str(work_dir / "chain_info.json"),
        ]
        if reuse_stage1_workdir.strip():
            orchestrator_args += ["--reuse_stage1_workdir", reuse_stage1_workdir.strip()]

        import subprocess as sp
        result = sp.run(orchestrator_args, capture_output=True, text=True, check=True)
        print(result.stdout, flush=True)

        chain_info = {}
        chain_path = work_dir / "chain_info.json"
        if chain_path.exists():
            with open(chain_path) as fh:
                chain_info = json.load(fh)

        update_job(job_id, status="running", progress_pct=5,
                   progress_label="Full chain submitted to Slurm")

        return JSONResponse({
            "job_id": job_id,
            "status": "submitted",
            "mode": "slurm",
            "source": "zebrafish",
            "stage": "full_pipeline",
            "chain_jids": chain_info.get("chain_jids", {}),
            "remote_run_dir": chain_info.get("remote_run_dir", ""),
        })
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"Orchestrator error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"Full pipeline submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"Full pipeline submission failed: {exc}")


# ---------------------------------------------------------------------------
# POST /jobs/zebrafish/stage3 — POISEN PE differential inclusion
# ---------------------------------------------------------------------------


@app.post("/jobs/zebrafish/stage3")
async def create_zebrafish_stage3_job(
    mode: str = Form("slurm"),
    stage1_job_id: str = Form(""),
    stage2_job_id: str = Form(""),
    prefix: str = Form("zebrafish_run"),
    max_intron_len: int = Form(500000),
    leafcutter_repo: str = Form("tools/leafcutter"),
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    zebrafish_pe_list: str = Form("UPDATED_STRAND_POISEN_HITS.tsv"),
    zebrafish_longorf_tsv: str = Form("LongORF_PTC+_fromClustered.tsv"),
    zebrafish_genome_fasta: str = Form("refs/GRCz12tu.fa"),
    zebrafish_chrom_map: str = Form("refs/GRCz12tu_assembly_report.txt"),
    pe_min_junction_reads: int = Form(3),
    pe_min_locus_reads: int = Form(10),
    ds_enabled: str = Form("false"),
) -> JSONResponse:
    mode = mode.strip().lower()
    if mode != "slurm":
        raise HTTPException(status_code=400, detail="Zebrafish stage3 currently supports slurm mode only.")
    stage1_job_id = stage1_job_id.strip()
    stage2_job_id = stage2_job_id.strip()
    if not stage1_job_id:
        raise HTTPException(status_code=400, detail="Provide stage1_job_id (alignment) for stage3 SJ files.")

    stage1 = get_job_or_404(stage1_job_id)
    stage1_input = stage1.get("input_payload", {}) or {}
    if stage1_input.get("source") != "zebrafish" or stage1_input.get("zebrafish_stage") != "stage1_alignment":
        raise HTTPException(status_code=400, detail="stage1_job_id is not a zebrafish stage1 job.")
    if stage1.get("status") != "succeeded":
        raise HTTPException(status_code=400, detail="Stage1 job must be succeeded before running stage3.")

    if stage2_job_id:
        stage2 = get_job_or_404(stage2_job_id)
        stage2_input = stage2.get("input_payload", {}) or {}
        if stage2_input.get("source") != "zebrafish" or stage2_input.get("zebrafish_stage") != "stage2_analysis":
            raise HTTPException(status_code=400, detail="stage2_job_id is not a zebrafish stage2 job.")

    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        raise HTTPException(status_code=400, detail=f"Slurm mode requires: {', '.join(missing)}")

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    work_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "prefix": prefix,
        "max_intron_len": max_intron_len,
        "leafcutter_repo": leafcutter_repo,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "zebrafish_pe_list": zebrafish_pe_list,
        "zebrafish_longorf_tsv": zebrafish_longorf_tsv,
        "pe_min_junction_reads": pe_min_junction_reads,
        "pe_min_locus_reads": pe_min_locus_reads,
        "ds_enabled": ds_enabled.strip().lower() == "true",
        "zebrafish_stage": "stage3_pe_inclusion",
        "stage1_job_id": stage1_job_id,
        "stage2_job_id": stage2_job_id or None,
    }
    input_payload = {
        "source": "zebrafish",
        "zebrafish_stage": "stage3_pe_inclusion",
        "stage1_job_id": stage1_job_id,
        "stage2_job_id": stage2_job_id or None,
    }
    _insert_job_record(
        job_id=job_id,
        mode="slurm",
        work_dir=work_dir,
        input_payload=input_payload,
        config=config,
        quest_account=quest_account,
        quest_partition=quest_partition,
    )

    try:
        if not remote_work_root.strip():
            remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
        if remote_work_root.startswith("~"):
            home_result = run_ssh_cmd(quest_host, quest_user, "echo $HOME")
            remote_home = home_result.stdout.strip()
            remote_work_root = remote_work_root.replace("~", remote_home, 1)
        config["remote_work_root"] = remote_work_root
        db_execute("UPDATE jobs SET config_payload=? WHERE id=?", (json.dumps(config), job_id))

        remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
        stage1_remote_dir = f"{remote_work_root.rstrip('/')}/{stage1_job_id}"
        stage2_remote_dir = f"{remote_work_root.rstrip('/')}/{stage2_job_id}" if stage2_job_id else ""
        run_ssh_cmd(quest_host, quest_user, f"mkdir -p {remote_run_dir}/inputs")
        run_ssh_cmd(
            quest_host,
            quest_user,
            f"[ -f \"{stage1_remote_dir}/out/stage1_alignment.json\" ] || (echo \"missing stage1 artifacts\" && exit 1)",
        )

        rr = remote_repo_root.rstrip("/")
        slurm_time = {"short": "04:00:00", "normal": "12:00:00", "long": "48:00:00"}.get(quest_partition, "12:00:00")
        sbatch_preamble = [
            "#!/bin/bash",
            f"#SBATCH --account={quest_account}",
            f"#SBATCH --partition={quest_partition}",
            f"#SBATCH --time={slurm_time}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=32G",
            f"#SBATCH --job-name=zf_s3_{job_id[:8]}",
            f"#SBATCH --output={remote_run_dir}/slurm_%j.out",
            f"#SBATCH --error={remote_run_dir}/slurm_%j.err",
            "",
            "set -eo pipefail",
            f"cd {rr}",
            "module purge",
            "module load anaconda3/2022.05",
            "",
        ]
        cmd_parts = [
            f"python3 {rr}/scripts/zebrafish_quest_pipeline.py",
            "--stage stage3",
            f"--workdir {remote_run_dir}",
            f"--stage1_workdir {stage1_remote_dir}",
            f"--prefix {prefix}",
            f"--max_intron_len {max_intron_len}",
            f"--pe_list {rr}/{zebrafish_pe_list}",
            f"--longorf_tsv {rr}/{zebrafish_longorf_tsv}",
            f"--pe_min_junction_reads {pe_min_junction_reads}",
            f"--pe_min_locus_reads {pe_min_locus_reads}",
            f"--genome_fasta {rr}/{zebrafish_genome_fasta}",
        ]
        if zebrafish_chrom_map.strip():
            cmd_parts.append(f"--chrom_map {rr}/{zebrafish_chrom_map}")
        if stage2_remote_dir:
            cmd_parts.append(f"--stage2_workdir {stage2_remote_dir}")
        if config["ds_enabled"]:
            cmd_parts.append("--ds")
            cmd_parts.append(f"--leafcutter_repo {rr}/{leafcutter_repo}")

        sbatch_body = "\n".join(sbatch_preamble + [" ".join(cmd_parts)]) + "\n"
        quest_job_id = _submit_slurm_and_start_poller(
            job_id=job_id,
            work_dir=work_dir,
            remote_run_dir=remote_run_dir,
            sbatch_body=sbatch_body,
            quest_host=quest_host,
            quest_user=quest_user,
        )
        return JSONResponse({
            "job_id": job_id,
            "status": "submitted",
            "mode": "slurm",
            "source": "zebrafish",
            "stage": "stage3_pe_inclusion",
            "stage1_job_id": stage1_job_id,
            "stage2_job_id": stage2_job_id or None,
            "quest_job_id": quest_job_id,
        })
    except HTTPException:
        raise
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"SSH/SCP error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"Stage3 submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"Stage3 submission failed: {exc}")


# ---------------------------------------------------------------------------
# POST /jobs — create and dispatch a pipeline job
# ---------------------------------------------------------------------------


@app.post("/jobs")
async def create_job(
    mode: str = Form("local"),
    source: str = Form("upload"),
    gtex_tissues_csv: str = Form(""),
    zebrafish_projects_csv: str = Form(""),
    prefix: str = Form("web_run"),
    min_reads: int = Form(50),
    max_intron_len: int = Form(500000),
    leafcutter_repo: str = Form("tools/leafcutter"),
    leafcutter2_repo: str = Form("tools/leafcutter2"),
    genome_fasta: str = Form("refs/GRCh38.fa"),
    gencode_gtf: str = Form("refs/gencode.v46.annotation.gtf"),
    files: List[UploadFile] = File(default=[]),
    samples_tsv: Optional[UploadFile] = File(default=None),
    # Quest/Slurm fields
    quest_host: str = Form("login.quest.northwestern.edu"),
    quest_user: str = Form(""),
    quest_account: str = Form(""),
    quest_partition: str = Form(""),
    remote_repo_root: str = Form(""),
    remote_work_root: str = Form(""),
    remote_juncfiles: str = Form(""),
    max_samples_per_tissue: int = Form(0),
    run_tissues_in_parallel: str = Form("false"),
    zebrafish_max_runs: int = Form(0),
    zebrafish_star_threads: int = Form(8),
    zebrafish_genome_fasta: str = Form("refs/GRCz12tu.fa"),
    zebrafish_gencode_gtf: str = Form("refs/GRCz12tu.ucsc.gtf"),
    zebrafish_gene_type_tag: str = Form("gene_biotype"),
    zebrafish_transcript_type_tag: str = Form("transcript_biotype"),
    zebrafish_gene_name_tag: str = Form("gene_id"),
    zebrafish_transcript_name_tag: str = Form("transcript_id"),
) -> JSONResponse:
    mode = mode.strip().lower()
    if mode not in {"local", "slurm"}:
        raise HTTPException(status_code=400, detail="mode must be 'local' or 'slurm'")
    source = source.strip().lower()
    if source not in {"upload", "gtex", "zebrafish"}:
        raise HTTPException(status_code=400, detail="source must be 'upload', 'gtex', or 'zebrafish'")
    if source == "zebrafish":
        raise HTTPException(
            status_code=400,
            detail="Use /jobs/zebrafish/stage1 and /jobs/zebrafish/stage2 endpoints for zebrafish workflows.",
        )

    job_id = str(uuid.uuid4())
    work_dir = JOBS_ROOT / job_id
    inputs_dir = work_dir / "inputs"
    sj_dir = inputs_dir / "star_sj"
    work_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded files
    local_sj_files: List[Path] = []
    for upload in files:
        if not upload.filename:
            continue
        local_sj_files.append(write_upload(upload, sj_dir / Path(upload.filename).name))

    local_samples: Optional[Path] = None
    if samples_tsv and samples_tsv.filename:
        local_samples = write_upload(samples_tsv, inputs_dir / "samples.tsv")

    config = {
        "prefix": prefix,
        "min_reads": min_reads,
        "max_intron_len": max_intron_len,
        "leafcutter_repo": leafcutter_repo,
        "leafcutter2_repo": leafcutter2_repo,
        "genome_fasta": genome_fasta,
        "gencode_gtf": gencode_gtf,
        "quest_host": quest_host,
        "quest_user": quest_user,
        "quest_account": quest_account,
        "quest_partition": quest_partition,
        "remote_repo_root": remote_repo_root,
        "remote_work_root": remote_work_root,
        "max_samples_per_tissue": max_samples_per_tissue,
        "run_tissues_in_parallel": run_tissues_in_parallel.strip().lower() == "true",
        "zebrafish_max_runs": zebrafish_max_runs,
        "zebrafish_star_threads": zebrafish_star_threads,
        "zebrafish_genome_fasta": zebrafish_genome_fasta,
        "zebrafish_gencode_gtf": zebrafish_gencode_gtf,
        "zebrafish_gene_type_tag": zebrafish_gene_type_tag,
        "zebrafish_transcript_type_tag": zebrafish_transcript_type_tag,
        "zebrafish_gene_name_tag": zebrafish_gene_name_tag,
        "zebrafish_transcript_name_tag": zebrafish_transcript_name_tag,
    }
    gtex_tissues_list = [t.strip() for t in gtex_tissues_csv.split(",") if t.strip()]
    zebrafish_projects_list = [p.strip() for p in zebrafish_projects_csv.split(",") if p.strip()]

    input_payload = {
        "source": source,
        "star_sj_files": [str(p) for p in local_sj_files],
        "samples_tsv": str(local_samples) if local_samples else None,
        "remote_juncfiles": remote_juncfiles or None,
        "gtex_tissues": gtex_tissues_list if source == "gtex" else None,
        "zebrafish_projects": zebrafish_projects_list if source == "zebrafish" else None,
    }

    db_execute(
        """
        INSERT INTO jobs
        (id, mode, status, created_at, updated_at, work_dir, input_payload, config_payload,
         quest_job_id, quest_account, quest_partition, error, summary_path, artifacts_zip, runner_ref)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, NULL, ?, ?, NULL, NULL, NULL, NULL)
        """,
        (
            job_id, mode, "queued", now_iso(), now_iso(), str(work_dir),
            json.dumps(input_payload), json.dumps(config),
            quest_account or None, quest_partition or None,
        ),
    )

    # ---- LOCAL MODE ----
    if mode == "local":
        if source == "gtex":
            if not gtex_tissues_list:
                update_job(job_id, status="failed", error="No GTEx tissues selected.")
                return JSONResponse(status_code=400, content={"job_id": job_id, "error": "No GTEx tissues selected."})
            t = threading.Thread(
                target=gtex_local_worker,
                args=(job_id, ",".join(gtex_tissues_list), work_dir, config),
                daemon=True,
            )
            ACTIVE_WORKERS[job_id] = t
            t.start()
            return JSONResponse({"job_id": job_id, "status": "queued", "mode": "local", "source": "gtex"})
        if source == "zebrafish":
            update_job(job_id, status="failed", error="Zebrafish source is currently supported only in Slurm mode.")
            return JSONResponse(
                status_code=400,
                content={"job_id": job_id, "error": "Zebrafish source is currently supported only in Slurm mode."},
            )

        # source == "upload" (default)
        if not local_sj_files:
            update_job(job_id, status="failed", error="No STAR SJ files uploaded.")
            return JSONResponse(status_code=400, content={"job_id": job_id, "error": "No STAR SJ files uploaded."})
        cmd = build_pipeline_cmd(
            run_dir=work_dir,
            star_sj_paths=local_sj_files,
            prefix=prefix,
            leafcutter_repo=leafcutter_repo,
            leafcutter2_repo=leafcutter2_repo,
            genome_fasta=genome_fasta,
            gencode_gtf=gencode_gtf,
            min_reads=min_reads,
            max_intron_len=max_intron_len,
            samples_tsv=local_samples,
        )
        t = threading.Thread(target=local_worker, args=(job_id, cmd, work_dir), daemon=True)
        ACTIVE_WORKERS[job_id] = t
        t.start()
        return JSONResponse({"job_id": job_id, "status": "queued", "mode": "local"})

    # ---- SLURM MODE ----
    missing = []
    if not quest_user:
        missing.append("quest_user")
    if not quest_account:
        missing.append("quest_account")
    if not quest_partition:
        missing.append("quest_partition")
    if not remote_repo_root:
        missing.append("remote_repo_root")
    if missing:
        msg = f"Slurm mode requires: {', '.join(missing)}"
        update_job(job_id, status="failed", error=msg)
        raise HTTPException(status_code=400, detail=msg)

    try:
        # Default work root to project space (avoids home directory quota).
        # Fall back to ~/leafcutter_jobs only if no repo root is provided.
        if not remote_work_root.strip():
            if remote_repo_root.strip():
                remote_work_root = f"{remote_repo_root.rstrip('/')}/jobs"
            else:
                remote_work_root = "~/leafcutter_jobs"

        # Resolve ~ to absolute home path on the remote host so that
        # SBATCH directives (which don't expand ~) and double-quoted
        # shell strings all work correctly.
        if remote_work_root.startswith("~"):
            home_result = run_ssh_cmd(quest_host, quest_user, "echo $HOME")
            remote_home = home_result.stdout.strip()
            remote_work_root = remote_work_root.replace("~", remote_home, 1)

        config["remote_work_root"] = remote_work_root
        db_execute(
            "UPDATE jobs SET config_payload=? WHERE id=?",
            (json.dumps(config), job_id),
        )

        remote_run_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
        run_ssh_cmd(quest_host, quest_user, f"mkdir -p {remote_run_dir}/inputs")

        rr = remote_repo_root.rstrip("/")
        slurm_time = {"short": "04:00:00", "normal": "48:00:00", "long": "168:00:00"}.get(
            quest_partition, "48:00:00"
        )
        sbatch_preamble = [
            "#!/bin/bash",
            f"#SBATCH --account={quest_account}",
            f"#SBATCH --partition={quest_partition}",
            f"#SBATCH --time={slurm_time}",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --mem=32G",
            f"#SBATCH --job-name=lc2_{job_id[:8]}",
            f"#SBATCH --output={remote_run_dir}/slurm_%j.out",
            f"#SBATCH --error={remote_run_dir}/slurm_%j.err",
            "",
            "set -eo pipefail",
            f"cd {rr}",
            "module purge",
            "module load anaconda3/2022.05",
            "",
        ]

        if source == "gtex":
            if not gtex_tissues_list:
                raise HTTPException(status_code=400, detail="No GTEx tissues selected.")
            tissues_escaped = ",".join(gtex_tissues_list)
            gtex_shared_cache = f"{rr}/gtex_cache"
            max_spt = config.get("max_samples_per_tissue", 0)
            parallel = config.get("run_tissues_in_parallel", False)

            # Use a shared cache so GTEx data is downloaded once across all jobs.
            # Download on the login node since compute nodes lack internet access.
            run_ssh_cmd(quest_host, quest_user, f"mkdir -p {gtex_shared_cache}")
            annot_remote = f"{gtex_shared_cache}/{GTEX_ANNOT_FILENAME}"
            gct_remote = f"{gtex_shared_cache}/{GTEX_GCT_FILENAME}"
            run_ssh_cmd(
                quest_host, quest_user,
                f'[ -f "{annot_remote}" ] || wget -q -O "{annot_remote}" "{GTEX_ANNOT_URL}"',
            )
            update_job(job_id, progress_pct=3, progress_label="Downloading GTEx junction data (may take a while)...")
            run_ssh_cmd(
                quest_host, quest_user,
                f'[ -f "{gct_remote}" ] || wget -q -O "{gct_remote}" "{GTEX_GCT_URL}"',
            )

            if parallel and len(gtex_tissues_list) > 1:
                # --- PARALLEL: Slurm job array per tissue ---
                # 1) Setup script: verify GTEx data is present
                setup_body = "\n".join(sbatch_preamble + [
                    f'ANNOT="{annot_remote}"',
                    f'GCT="{gct_remote}"',
                    '',
                    '[ -f "$ANNOT" ] || { echo "ERROR: annotations file missing" >&2; exit 1; }',
                    '[ -f "$GCT" ] || { echo "ERROR: GCT file missing" >&2; exit 1; }',
                    "",
                    'echo "GTEx data ready"',
                ]) + "\n"
                setup_script = work_dir / "setup.sbatch"
                setup_script.write_text(setup_body)

                # 2) Write the tissue list as a file for the array tasks to index
                tissue_list_local = work_dir / "tissues.txt"
                tissue_list_local.write_text("\n".join(gtex_tissues_list) + "\n")

                # 3) Array task script: each task processes one tissue
                n_tissues = len(gtex_tissues_list)
                array_preamble = [
                    "#!/bin/bash",
                    f"#SBATCH --account={quest_account}",
                    f"#SBATCH --partition={quest_partition}",
                    f"#SBATCH --time={slurm_time}",
                    "#SBATCH --nodes=1",
                    "#SBATCH --ntasks=1",
                    "#SBATCH --mem=32G",
                    f"#SBATCH --job-name=lc2_arr_{job_id[:8]}",
                    f"#SBATCH --output={remote_run_dir}/slurm_%A_%a.out",
                    f"#SBATCH --error={remote_run_dir}/slurm_%A_%a.err",
                    f"#SBATCH --array=0-{n_tissues - 1}",
                    "",
                    "set -eo pipefail",
                    f"cd {rr}",
                    "module purge",
                    "module load anaconda3/2022.05",
                    "",
                ]
                array_body = "\n".join(array_preamble + [
                    f'TISSUE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {remote_run_dir}/tissues.txt)',
                    'echo "Processing tissue: $TISSUE (task $SLURM_ARRAY_TASK_ID)"',
                    "",
                    '# Sanitize tissue name for directory',
                    "TDIR_NAME=$(echo \"$TISSUE\" | sed 's/ - /_/g; s/ /_/g; s/[()]//g')",
                    f'TISSUE_WORKDIR="{remote_run_dir}/tissue_${{TDIR_NAME}}"',
                    'mkdir -p "$TISSUE_WORKDIR"',
                    "",
                    f'ANNOT="{annot_remote}"',
                    f'GCT="{gct_remote}"',
                    "",
                    (
                        f'python3 {rr}/scripts/gtex_gct_to_bed.py'
                        f' --gct "$GCT"'
                        f' --annotations "$ANNOT"'
                        f' --tissues "$TISSUE"'
                        f' --outdir "$TISSUE_WORKDIR/junctions_bed"'
                        f' --max_samples_per_tissue {max_spt}'
                    ),
                    "",
                    (
                        f'python3 {rr}/scripts/lc2_pipeline.py'
                        f' --workdir "$TISSUE_WORKDIR"'
                        f' --prefix {prefix}'
                        f' --leafcutter_repo {rr}/{leafcutter_repo}'
                        f' --leafcutter2_repo {rr}/{leafcutter2_repo}'
                        f' --genome_fasta {rr}/{genome_fasta}'
                        f' --gencode_gtf {rr}/{gencode_gtf}'
                        f' --min_reads {min_reads}'
                        f' --max_intron_len {max_intron_len}'
                        f' --junction_beds $(cat "$TISSUE_WORKDIR/junctions_bed/junction_files.txt")'
                    ),
                ]) + "\n"
                array_script = work_dir / "array.sbatch"
                array_script.write_text(array_body)

                # 4) Merge script: runs after all array tasks finish
                merge_preamble = list(sbatch_preamble)  # copy
                merge_body = "\n".join(merge_preamble + [
                    (
                        f'python3 {rr}/scripts/merge_tissue_results.py'
                        f' --rundir {remote_run_dir}'
                        f' --tissues "{tissues_escaped}"'
                        f' --prefix {prefix}'
                    ),
                ]) + "\n"
                merge_script = work_dir / "merge.sbatch"
                merge_script.write_text(merge_body)

                # Preflight: validate all generated scripts
                for script in (setup_script, array_script, merge_script):
                    validate_shell_script(script)

                # Upload all scripts and tissues.txt
                remote_setup = f"{remote_run_dir}/setup.sbatch"
                remote_array = f"{remote_run_dir}/array.sbatch"
                remote_merge = f"{remote_run_dir}/merge.sbatch"
                remote_tissues = f"{remote_run_dir}/tissues.txt"
                run_scp_cmd(
                    [str(setup_script), str(array_script), str(merge_script), str(tissue_list_local)],
                    f"{quest_user}@{quest_host}:{remote_run_dir}/",
                )

                # Submit chain: setup -> array (--dependency) -> merge (--dependency)
                submit_setup = run_ssh_cmd(quest_host, quest_user, f"sbatch --parsable {remote_setup}")
                setup_jid = submit_setup.stdout.strip().split(";")[0]
                if not setup_jid:
                    raise RuntimeError("sbatch setup returned no job ID")

                submit_array = run_ssh_cmd(
                    quest_host, quest_user,
                    f"sbatch --parsable --dependency=afterok:{setup_jid} {remote_array}",
                )
                array_jid = submit_array.stdout.strip().split(";")[0]
                if not array_jid:
                    raise RuntimeError("sbatch array returned no job ID")

                submit_merge = run_ssh_cmd(
                    quest_host, quest_user,
                    f"sbatch --parsable --dependency=afterany:{array_jid} {remote_merge}",
                )
                merge_jid = submit_merge.stdout.strip().split(";")[0]
                if not merge_jid:
                    raise RuntimeError("sbatch merge returned no job ID")

                quest_job_id = merge_jid
                update_job(
                    job_id,
                    status="submitted",
                    quest_job_id=quest_job_id,
                    progress_label=f"Submitted: setup={setup_jid}, array={array_jid}, merge={merge_jid}",
                )

                # Store all Slurm IDs for tracking/cancellation
                config["slurm_setup_jid"] = setup_jid
                config["slurm_array_jid"] = array_jid
                config["slurm_merge_jid"] = merge_jid
                db_execute(
                    "UPDATE jobs SET config_payload=? WHERE id=?",
                    (json.dumps(config), job_id),
                )

                t = threading.Thread(
                    target=slurm_array_worker,
                    args=(job_id, quest_host, quest_user, remote_run_dir,
                          setup_jid, array_jid, merge_jid, n_tissues),
                    daemon=True,
                )
                ACTIVE_WORKERS[job_id] = t
                t.start()
                return JSONResponse({
                    "job_id": job_id, "status": "submitted", "mode": "slurm",
                    "quest_job_id": quest_job_id, "parallel_tissues": True,
                    "n_tissues": n_tissues,
                })

            # --- SEQUENTIAL (legacy): single sbatch for all tissues ---
            bed_dir_remote = f"{remote_run_dir}/junctions_bed"
            sbatch_body = "\n".join(sbatch_preamble + [
                f'ANNOT="{annot_remote}"',
                f'GCT="{gct_remote}"',
                "",
                '[ -f "$ANNOT" ] || { echo "ERROR: annotations file missing" >&2; exit 1; }',
                '[ -f "$GCT" ] || { echo "ERROR: GCT file missing" >&2; exit 1; }',
                "",
                (
                    f'python3 {rr}/scripts/gtex_gct_to_bed.py'
                    f' --gct "$GCT"'
                    f' --annotations "$ANNOT"'
                    f' --tissues "{tissues_escaped}"'
                    f' --outdir {bed_dir_remote}'
                    f' --max_samples_per_tissue {max_spt}'
                ),
                "",
                (
                    f"python3 {rr}/scripts/lc2_pipeline.py"
                    f" --workdir {remote_run_dir}"
                    f" --prefix {prefix}"
                    f" --leafcutter_repo {rr}/{leafcutter_repo}"
                    f" --leafcutter2_repo {rr}/{leafcutter2_repo}"
                    f" --genome_fasta {rr}/{genome_fasta}"
                    f" --gencode_gtf {rr}/{gencode_gtf}"
                    f" --min_reads {min_reads}"
                    f" --max_intron_len {max_intron_len}"
                    f" --junction_beds $(cat {bed_dir_remote}/junction_files.txt)"
                ),
            ]) + "\n"
        elif source == "zebrafish":
            if not zebrafish_projects_list:
                raise HTTPException(
                    status_code=400,
                    detail="No zebrafish project IDs provided. Add one or more PRJNA/ERP/DRP IDs.",
                )
            zf_projects = ",".join(zebrafish_projects_list)
            zf_threads = int(config.get("zebrafish_star_threads", 8))
            zf_max_runs = int(config.get("zebrafish_max_runs", 0))
            zf_fasta = str(config.get("zebrafish_genome_fasta", "refs/GRCz12tu.fa"))
            zf_gtf = str(config.get("zebrafish_gencode_gtf", "refs/GRCz12tu.ucsc.gtf"))
            zf_gene_type = str(config.get("zebrafish_gene_type_tag", "gene_biotype"))
            zf_tx_type = str(config.get("zebrafish_transcript_type_tag", "transcript_biotype"))
            zf_gene_name = str(config.get("zebrafish_gene_name_tag", "gene_id"))
            zf_tx_name = str(config.get("zebrafish_transcript_name_tag", "transcript_id"))
            zf_manifest = f"{remote_run_dir}/inputs/zebrafish_runs.tsv"

            # Discover public runs on the login node so compute nodes do not
            # need direct internet access for ENA metadata lookup.
            run_ssh_cmd(
                quest_host,
                quest_user,
                (
                    f"python3 {rr}/scripts/zebrafish_project_runs.py"
                    f" --project_ids \"{zf_projects}\""
                    f" --out \"{zf_manifest}\""
                    f" --max_runs {zf_max_runs}"
                ),
            )

            sbatch_body = "\n".join(sbatch_preamble + [
                (
                    f"python3 {rr}/scripts/zebrafish_quest_pipeline.py"
                    f" --workdir {remote_run_dir}"
                    f" --manifest \"{zf_manifest}\""
                    f" --prefix {prefix}"
                    f" --leafcutter_repo {rr}/{leafcutter_repo}"
                    f" --leafcutter2_repo {rr}/{leafcutter2_repo}"
                    f" --genome_fasta {rr}/{zf_fasta}"
                    f" --gencode_gtf {rr}/{zf_gtf}"
                    f" --min_reads {min_reads}"
                    f" --max_intron_len {max_intron_len}"
                    f" --star_threads {zf_threads}"
                    f" --max_runs {zf_max_runs}"
                    f" --gene_type_tag {zf_gene_type}"
                    f" --transcript_type_tag {zf_tx_type}"
                    f" --gene_name_tag {zf_gene_name}"
                    f" --transcript_name_tag {zf_tx_name}"
                ),
            ]) + "\n"
        else:
            # Upload source — resolve junction file list on Quest
            remote_junc_list = remote_juncfiles.strip()
            if not remote_junc_list:
                if not local_sj_files:
                    raise HTTPException(status_code=400, detail="Provide remote_juncfiles or upload STAR SJ files.")
                remote_lines: List[str] = []
                for p in local_sj_files:
                    rpath = f"{remote_run_dir}/inputs/{p.name}"
                    run_scp_cmd([str(p)], f"{quest_user}@{quest_host}:{rpath}")
                    remote_lines.append(rpath)
                local_list = work_dir / "remote_junction_files.txt"
                local_list.write_text("\n".join(remote_lines) + "\n")
                remote_junc_list = f"{remote_run_dir}/inputs/junction_files.txt"
                run_scp_cmd([str(local_list)], f"{quest_user}@{quest_host}:{remote_junc_list}")

            sbatch_body = "\n".join(sbatch_preamble + [
                (
                    f"python3 {rr}/scripts/lc2_pipeline.py"
                    f" --workdir {remote_run_dir}"
                    f" --prefix {prefix}"
                    f" --leafcutter_repo {rr}/{leafcutter_repo}"
                    f" --leafcutter2_repo {rr}/{leafcutter2_repo}"
                    f" --genome_fasta {rr}/{genome_fasta}"
                    f" --gencode_gtf {rr}/{gencode_gtf}"
                    f" --min_reads {min_reads}"
                    f" --max_intron_len {max_intron_len}"
                    f" --star_sj $(cat {remote_junc_list})"
                ),
            ]) + "\n"

        slurm_script = work_dir / "job.sbatch"
        slurm_script.write_text(sbatch_body)
        validate_shell_script(slurm_script)
        remote_script = f"{remote_run_dir}/job.sbatch"
        run_scp_cmd([str(slurm_script)], f"{quest_user}@{quest_host}:{remote_script}")

        submit = run_ssh_cmd(quest_host, quest_user, f"sbatch --parsable {remote_script}")
        quest_job_id = submit.stdout.strip().split(";")[0]
        if not quest_job_id:
            raise RuntimeError("sbatch returned no job ID")

        update_job(job_id, status="submitted", quest_job_id=quest_job_id)

        t = threading.Thread(
            target=slurm_worker,
            args=(job_id, quest_host, quest_user, remote_run_dir, f"{remote_run_dir}/out"),
            daemon=True,
        )
        ACTIVE_WORKERS[job_id] = t
        t.start()
        return JSONResponse({
            "job_id": job_id, "status": "submitted", "mode": "slurm", "quest_job_id": quest_job_id,
        })

    except HTTPException:
        raise
    except subprocess.CalledProcessError as exc:
        ssh_err = exc.stderr.strip() if exc.stderr else str(exc)
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"SSH/SCP error: {ssh_err}")
        raise HTTPException(status_code=500, detail=f"Slurm submission failed: {ssh_err}")
    except Exception as exc:
        update_job(job_id, status="failed", progress_label="Submission failed", error=f"{type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail=f"Slurm submission failed: {exc}")


# ---------------------------------------------------------------------------
# GET /jobs — list all jobs
# ---------------------------------------------------------------------------


@app.get("/jobs")
def list_jobs() -> JSONResponse:
    rows = db_fetch_all("SELECT * FROM jobs ORDER BY created_at DESC LIMIT 100")
    return JSONResponse([_row_to_dict(r) for r in rows])


# ---------------------------------------------------------------------------
# GET /jobs/{id}
# ---------------------------------------------------------------------------


@app.get("/jobs/{job_id}")
def get_job_status(job_id: str) -> JSONResponse:
    return JSONResponse(get_job_or_404(job_id))


# ---------------------------------------------------------------------------
# GET /jobs/{id}/results
# ---------------------------------------------------------------------------


@app.get("/jobs/{job_id}/results")
def get_job_results(job_id: str) -> JSONResponse:
    job = get_job_or_404(job_id)
    summary = None
    sp = job.get("summary_path")
    if sp and Path(sp).exists():
        try:
            summary = json.loads(Path(sp).read_text())
        except Exception:
            pass

    log_tail = None
    log_path = Path(job["work_dir"]) / "pipeline.log"
    if log_path.exists():
        try:
            lines = log_path.read_text().splitlines()
            log_tail = lines[-50:]
        except Exception:
            pass

    return JSONResponse({
        "job_id": job_id,
        "status": job["status"],
        "error": job.get("error"),
        "summary": summary,
        "log_tail": log_tail,
    })


# ---------------------------------------------------------------------------
# GET /jobs/{id}/logs — incremental log streaming
# ---------------------------------------------------------------------------


@app.get("/jobs/{job_id}/logs")
def get_job_logs(job_id: str, offset: int = 0, limit: int = 65536) -> JSONResponse:
    job = get_job_or_404(job_id)
    log_path = Path(job["work_dir"]) / "pipeline.log"
    terminal = job["status"] in {"succeeded", "failed", "cancelled"}

    if not log_path.exists():
        return JSONResponse({"lines": [], "next_offset": 0, "eof": terminal})

    try:
        size = log_path.stat().st_size
    except OSError:
        return JSONResponse({"lines": [], "next_offset": offset, "eof": True})

    if offset >= size:
        return JSONResponse({"lines": [], "next_offset": offset, "eof": terminal})

    with open(log_path, "rb") as f:
        f.seek(offset)
        raw = f.read(limit)

    new_offset = offset + len(raw)
    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    eof = new_offset >= size and terminal

    return JSONResponse({"lines": lines, "next_offset": new_offset, "eof": eof})


# ---------------------------------------------------------------------------
# GET /jobs/{id}/download
# ---------------------------------------------------------------------------


@app.get("/jobs/{job_id}/download")
def download_artifacts(job_id: str) -> FileResponse:
    job = get_job_or_404(job_id)
    zp = job.get("artifacts_zip")
    if not zp or not Path(zp).exists():
        raise HTTPException(status_code=404, detail="Artifacts not available yet.")
    return FileResponse(path=zp, filename=f"{job_id}_artifacts.zip", media_type="application/zip")


# ---------------------------------------------------------------------------
# GET /jobs/{id}/artifacts/{filename} — single artifact download
# ---------------------------------------------------------------------------

_ARTIFACT_MIME = {
    ".png": "image/png",
    ".tsv": "text/tab-separated-values",
    ".json": "application/json",
    ".txt": "text/plain",
    ".gz": "application/gzip",
    ".html": "text/html",
}


@app.get("/jobs/{job_id}/artifacts/{filename:path}")
def download_single_artifact(job_id: str, filename: str) -> FileResponse:
    job = get_job_or_404(job_id)
    if ".." in filename.split("/") or filename.startswith("/"):
        raise HTTPException(status_code=400, detail="Invalid filename")

    out_root = (Path(job["work_dir"]) / "out").resolve()
    artifact = (out_root / filename).resolve()
    try:
        artifact.relative_to(out_root)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid filename")
    if not artifact.exists() or not artifact.is_file():
        raise HTTPException(status_code=404, detail=f"Artifact '{filename}' not found")

    mime = _ARTIFACT_MIME.get(artifact.suffix, "application/octet-stream")
    return FileResponse(path=str(artifact), filename=artifact.name, media_type=mime)


# ---------------------------------------------------------------------------
# POST /jobs/{id}/cancel
# ---------------------------------------------------------------------------


@app.post("/jobs/{job_id}/cancel")
def cancel_job(job_id: str) -> JSONResponse:
    job = get_job_or_404(job_id)

    if job["status"] in {"succeeded", "failed", "cancelled"}:
        return JSONResponse({"job_id": job_id, "status": job["status"], "message": "Job already terminal."})

    # Local cancel
    if job["mode"] == "local":
        proc = ACTIVE_PROCS.get(job_id)
        if proc and proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
        else:
            pid_str = job.get("runner_ref")
            if pid_str:
                try:
                    os.kill(int(pid_str), signal.SIGTERM)
                except (ProcessLookupError, ValueError, PermissionError):
                    pass
        update_job(job_id, status="cancelled", error="Cancelled by user")
        ACTIVE_PROCS.pop(job_id, None)
        return JSONResponse({"job_id": job_id, "status": "cancelled"})

    # Slurm cancel — cancel all related job IDs (parallel chain or single)
    cfg = job.get("config_payload", {})
    slurm_id = job.get("quest_job_id")
    host = cfg.get("quest_host", "login.quest.northwestern.edu")
    user = cfg.get("quest_user")
    jids_to_cancel = set()
    if slurm_id:
        jids_to_cancel.add(slurm_id)
    for key in ("slurm_setup_jid", "slurm_array_jid", "slurm_merge_jid"):
        val = cfg.get(key)
        if val:
            jids_to_cancel.add(val)
    if jids_to_cancel and user:
        for jid in jids_to_cancel:
            try:
                run_ssh_cmd(host, user, f"scancel {jid}")
            except Exception:
                pass
    update_job(job_id, status="cancelled", error="Cancelled by user")
    return JSONResponse({"job_id": job_id, "status": "cancelled", "quest_job_id": slurm_id})


# ---------------------------------------------------------------------------
# DELETE /jobs/{id} — delete a single job (files + DB row)
# ---------------------------------------------------------------------------


def _delete_job_impl(job: Dict[str, Any]) -> Dict[str, Any]:
    """Shared logic: cancel if active, remove local/remote files, delete DB row."""
    job_id = job["id"]
    result: Dict[str, Any] = {"job_id": job_id, "local_cleaned": False, "remote_cleaned": False}

    # Cancel if still active
    if job["status"] in ACTIVE:
        if job["mode"] == "local":
            proc = ACTIVE_PROCS.get(job_id)
            if proc and proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
            ACTIVE_PROCS.pop(job_id, None)
        else:
            cfg = job.get("config_payload", {})
            host = cfg.get("quest_host", "login.quest.northwestern.edu")
            user = cfg.get("quest_user")
            jids = set()
            if job.get("quest_job_id"):
                jids.add(job["quest_job_id"])
            for key in ("slurm_setup_jid", "slurm_array_jid", "slurm_merge_jid"):
                val = cfg.get(key)
                if val:
                    jids.add(val)
            if jids and user:
                for jid in jids:
                    try:
                        run_ssh_cmd(host, user, f"scancel {jid}")
                    except Exception:
                        pass
        ACTIVE_WORKERS.pop(job_id, None)

    # Remove local work directory
    work_dir = Path(job["work_dir"])
    if work_dir.exists():
        shutil.rmtree(work_dir, ignore_errors=True)
        result["local_cleaned"] = True

    # Remove remote job directory (best-effort)
    if job["mode"] == "slurm":
        cfg = job.get("config_payload", {})
        host = cfg.get("quest_host", "login.quest.northwestern.edu")
        user = cfg.get("quest_user")
        remote_work_root = cfg.get("remote_work_root", "")
        if user and remote_work_root:
            remote_dir = f"{remote_work_root.rstrip('/')}/{job_id}"
            try:
                run_ssh_cmd(host, user, f"rm -rf {remote_dir}")
                result["remote_cleaned"] = True
            except Exception:
                pass

    # Delete DB row
    db_execute("DELETE FROM jobs WHERE id=?", (job_id,))
    return result


ACTIVE = {"queued", "submitted", "running"}


@app.delete("/jobs/{job_id}")
def delete_job(job_id: str) -> JSONResponse:
    job = get_job_or_404(job_id)
    result = _delete_job_impl(job)
    return JSONResponse(result)


# ---------------------------------------------------------------------------
# DELETE /jobs — delete all jobs
# ---------------------------------------------------------------------------


@app.delete("/jobs")
def delete_all_jobs() -> JSONResponse:
    rows = db_fetch_all("SELECT * FROM jobs")
    deleted = 0
    for row in rows:
        try:
            _delete_job_impl(_row_to_dict(row))
            deleted += 1
        except Exception:
            pass
    return JSONResponse({"deleted": deleted, "total": len(rows)})


def main() -> None:
    """Console entry point: launch the LeafCutter2 web UI via uvicorn."""
    import os
    import uvicorn

    host = os.environ.get("LEAFCUTTER2_HOST", "127.0.0.1")
    port = int(os.environ.get("LEAFCUTTER2_PORT", "8000"))
    uvicorn.run("webapp.backend.main:app", host=host, port=port, reload=False)


if __name__ == "__main__":
    main()
