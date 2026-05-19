#!/usr/bin/env python3
"""Validate a GTF before using it with LeafCutter2 classification."""

import argparse
import gzip
import re
import subprocess
from collections import Counter
from pathlib import Path
from typing import Dict, Set


ATTR_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_.-]*)\s+"([^"]*)"')


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_attrs(attr: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(attr)}


def fasta_names(path: Path) -> Set[str]:
    fai = path.with_suffix(path.suffix + ".fai")
    if fai.exists():
        with fai.open() as fh:
            return {line.split("\t", 1)[0] for line in fh if line.strip()}
    if path.suffix == ".gz":
        opener = lambda: gzip.open(path, "rt")
    else:
        opener = path.open
    names = set()  # type: Set[str]
    with opener() as fh:
        for line in fh:
            if line.startswith(">"):
                names.add(line[1:].split()[0])
    return names


def ensure_fai(path: Path) -> None:
    if path.with_suffix(path.suffix + ".fai").exists():
        return
    try:
        subprocess.run(["samtools", "faidx", str(path)], check=True)
    except Exception:
        return


def validate(args: argparse.Namespace) -> int:
    required = [
        args.gene_name_tag,
        args.transcript_name_tag,
        args.gene_type_tag,
        args.transcript_type_tag,
    ]
    required = [x for x in required if x]

    rows = comments = parseable = malformed_attrs = 0
    feature_counts: Counter[str] = Counter()
    tag_counts: Counter[str] = Counter()
    missing_gene_id_rows = missing_tx_id_rows = 0
    seqnames: Counter[str] = Counter()
    samples = []

    with open_text(args.gtf) as fh:
        for lineno, line in enumerate(fh, start=1):
            if line.startswith("#"):
                comments += 1
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                print(f"ERROR: {args.gtf}:{lineno}: expected 9 tab-delimited fields, found {len(fields)}")
                return 2
            rows += 1
            feature = fields[2]
            feature_counts[feature] += 1
            seqnames[fields[0]] += 1
            attrs = parse_attrs(fields[8])
            if attrs:
                parseable += 1
            else:
                malformed_attrs += 1
                if malformed_attrs <= 5:
                    print(f"ERROR: {args.gtf}:{lineno}: no parseable key \"value\" attributes")

            for tag in required:
                if tag in attrs:
                    tag_counts[tag] += 1

            if feature in {"transcript", "exon", "CDS", "stop_codon", "start_codon"}:
                if args.gene_name_tag and args.gene_name_tag not in attrs:
                    missing_gene_id_rows += 1
                if args.transcript_name_tag and args.transcript_name_tag not in attrs:
                    missing_tx_id_rows += 1

            if len(samples) < args.sample_rows and feature in {"gene", "transcript", "exon", "CDS"}:
                samples.append((fields[0], lineno, feature, attrs))

    if rows == 0:
        print("ERROR: no non-comment GTF rows found")
        return 2

    print(f"GTF: {args.gtf}")
    print(f"Rows: {rows}; comments: {comments}; parseable_attrs: {parseable}")
    print(f"Feature counts: {dict(feature_counts.most_common(12))}")
    print(f"Top seqnames: {seqnames.most_common(10)}")
    print(f"Required tag counts: {dict(tag_counts)}")
    for seq, lineno, feature, attrs in samples:
        preview = {k: attrs.get(k) for k in required if k in attrs}
        print(f"Sample row {lineno} {seq} {feature}: {preview}")

    if parseable / rows < args.min_parseable_frac:
        print(f"ERROR: parseable attribute fraction {parseable / rows:.3f} below {args.min_parseable_frac}")
        return 2

    for tag in required:
        frac = tag_counts[tag] / rows
        if frac < args.min_tag_frac:
            print(f"ERROR: tag {tag!r} present in only {frac:.3f} of rows, below {args.min_tag_frac}")
            return 2

    relevant = sum(feature_counts[f] for f in ["transcript", "exon", "CDS", "stop_codon", "start_codon"])
    if relevant:
        gene_missing_frac = missing_gene_id_rows / relevant
        tx_missing_frac = missing_tx_id_rows / relevant
        if gene_missing_frac > args.max_missing_identifier_frac:
            print(f"ERROR: missing {args.gene_name_tag} fraction {gene_missing_frac:.3f} too high")
            return 2
        if tx_missing_frac > args.max_missing_identifier_frac:
            print(f"ERROR: missing {args.transcript_name_tag} fraction {tx_missing_frac:.3f} too high")
            return 2

    if args.fasta:
        ensure_fai(args.fasta)
        fasta = fasta_names(args.fasta)
        overlap = set(seqnames) & fasta
        print(f"FASTA seqnames: {len(fasta)}; GTF seqnames: {len(seqnames)}; overlap: {len(overlap)}")
        if not overlap:
            print("ERROR: no overlap between GTF seqnames and FASTA seqnames")
            return 2

    print("GTF validation passed")
    return 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("gtf", type=Path)
    p.add_argument("--fasta", type=Path)
    p.add_argument("--gene-name-tag", default="gene_name")
    p.add_argument("--transcript-name-tag", default="transcript_id")
    p.add_argument("--gene-type-tag", default="gene_biotype")
    p.add_argument("--transcript-type-tag", default="transcript_biotype")
    p.add_argument("--min-parseable-frac", type=float, default=0.95)
    p.add_argument("--min-tag-frac", type=float, default=0.05)
    p.add_argument("--max-missing-identifier-frac", type=float, default=0.05)
    p.add_argument("--sample-rows", type=int, default=5)
    raise SystemExit(validate(p.parse_args()))


if __name__ == "__main__":
    main()
