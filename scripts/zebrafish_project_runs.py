#!/usr/bin/env python3
"""
Resolve public zebrafish project IDs to a normalized run manifest TSV.

This script queries ENA's filereport endpoint and writes a manifest used by
the Quest zebrafish workflow.
"""


import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.parse import urlencode
from urllib.request import urlopen


ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = [
    "run_accession",
    "study_accession",
    "study_title",
    "sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "description",
    "scientific_name",
    "library_layout",
    "fastq_ftp",
    "submitted_ftp",
    "sra_ftp",
]


STAGE_RULES: List[Tuple[str, str]] = [
    (r"(\d+)\s*hpf", "hpf"),
    (r"(\d+)\s*dpf", "dpf"),
    (r"(\d+)\s*cell", "cell"),
    (r"(\d+)\s*somite", "somite"),
]

NAMED_STAGES: List[Tuple[str, str]] = [
    (r"\bdome\b", "dome"),
    (r"\bshield\b", "shield"),
    (r"\bbud\b", "bud"),
    (r"\bsphere\b", "sphere"),
    (r"\boblong\b", "oblong"),
    (r"\bgastrula\b", "gastrula"),
    (r"\bpharyngula\b", "pharyngula"),
    (r"\bsegmentation\b", "segmentation"),
    (r"\bblastula\b", "blastula"),
    (r"\bcleavage\b", "cleavage"),
    (r"\bepiboly\b", "epiboly"),
    (r"\bhatching\b", "hatching"),
    (r"\blarva(l|e)?\b", "larval"),
]


def fetch_project_rows(project_id: str) -> List[Dict[str, str]]:
    params = {
        "accession": project_id,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
    }
    url = f"{ENA_BASE}?{urlencode(params)}"
    with urlopen(url, timeout=120) as resp:
        payload = resp.read().decode("utf-8", errors="replace")
    rows = list(csv.DictReader(payload.splitlines(), delimiter="\t"))
    return rows


def infer_stage(row: Dict[str, str]) -> Tuple[str, str, str]:
    """Extract a developmental stage label from ENA sample metadata.

    Returns (raw_text, normalized_label, confidence).
    """
    text_fields = [
        row.get("sample_alias", ""),
        row.get("sample_title", ""),
        row.get("sample_description", ""),
        row.get("description", ""),
        row.get("study_title", ""),
    ]
    raw = " | ".join(x.strip() for x in text_fields if x and x.strip())
    blob = " ".join(x.strip().lower() for x in text_fields if x and x.strip())
    blob = re.sub(r"[^a-z0-9\s_-]+", " ", blob)
    blob = re.sub(r"\s+", " ", blob).strip()

    for pattern, unit in STAGE_RULES:
        m = re.search(pattern, blob)
        if m:
            label = f"{m.group(1)}{unit}"
            return raw or label, label, "high"

    for pattern, label in NAMED_STAGES:
        if re.search(pattern, blob):
            return raw or label, label, "high"

    if raw:
        return raw, "Unknown_Stage", "low"
    return "NA", "Unknown_Stage", "none"


def normalize_rows(project_ids: List[str], max_runs: int) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    for project_id in project_ids:
        rows = fetch_project_rows(project_id)
        for row in rows:
            sci = (row.get("scientific_name") or "").strip().lower()
            if sci and "danio rerio" not in sci:
                continue
            run_id = (row.get("run_accession") or "").strip()
            if not run_id:
                continue
            study = (row.get("study_accession") or project_id).strip() or project_id
            sample = (row.get("sample_accession") or run_id).strip() or run_id
            stage_raw, stage_norm, stage_conf = infer_stage(row)
            rec = {
                "run_id": run_id,
                "project_id": project_id,
                "sample": sample,
                "condition": stage_norm,
                "stage_raw": stage_raw,
                "stage": stage_norm,
                "stage_confidence": stage_conf,
                "library_layout": (row.get("library_layout") or "").strip().upper() or "SINGLE",
                "fastq_ftp": (row.get("fastq_ftp") or "").strip(),
                "submitted_ftp": (row.get("submitted_ftp") or "").strip(),
                "sra_ftp": (row.get("sra_ftp") or "").strip(),
                "study_accession": study,
                "study_title": (row.get("study_title") or "").strip(),
            }
            records.append(rec)

    # De-duplicate by run_id while preserving first appearance.
    uniq: List[Dict[str, str]] = []
    seen = set()
    for rec in records:
        run_id = rec["run_id"]
        if run_id in seen:
            continue
        seen.add(run_id)
        uniq.append(rec)

    if max_runs > 0:
        return uniq[:max_runs]
    return uniq


def main() -> None:
    p = argparse.ArgumentParser(description="Fetch zebrafish runs from ENA projects.")
    p.add_argument("--project_ids", required=True, help="Comma-separated project IDs (PRJNA/ERP/DRP).")
    p.add_argument("--out", required=True, help="Output manifest TSV path.")
    p.add_argument("--max_runs", type=int, default=0, help="Optional cap on total runs (0=all).")
    args = p.parse_args()

    project_ids = [x.strip() for x in args.project_ids.split(",") if x.strip()]
    if not project_ids:
        raise SystemExit("No project IDs provided.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    records = normalize_rows(project_ids, args.max_runs)
    if not records:
        raise SystemExit("No zebrafish runs found for supplied project IDs.")

    fields = [
        "run_id",
        "project_id",
        "sample",
        "condition",
        "stage",
        "stage_raw",
        "stage_confidence",
        "study_accession",
        "study_title",
        "library_layout",
        "fastq_ftp",
        "submitted_ftp",
        "sra_ftp",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

    print(f"[ZEBRAFISH] Projects: {','.join(project_ids)}")
    print(f"[ZEBRAFISH] Runs discovered: {len(records)}")
    print(f"[ZEBRAFISH] Manifest: {out_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
