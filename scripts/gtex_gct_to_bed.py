#!/usr/bin/env python3
"""
gtex_gct_to_bed.py

Convert the GTEx V8 STAR junction GCT matrix (gzipped) into per-sample BED6
files filtered by tissue, ready for LeafCutter/LeafCutter2.

The GCT is streamed line-by-line so the full ~4 GB matrix is never held in RAM.
Only columns (samples) belonging to the requested tissue(s) are written out.

Usage:
    python gtex_gct_to_bed.py \\
        --gct  GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz \\
        --annotations GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt \\
        --tissues "Brain - Cortex,Liver" \\
        --outdir /path/to/output \\
        --min_count 1

Outputs:
    <outdir>/<TissueName>/<sample>.juncs.bed   (one per matching sample)
    <outdir>/junction_files.txt                (list of all BED paths)
"""


import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

STAR_STRAND_MAP = {"0": ".", "1": "+", "2": "-"}


def parse_junction_name(name: str) -> Optional[Tuple[str, int, int, str]]:
    """Parse a STAR junction row name into (chrom, bed_start, bed_end, strand).

    STAR GCT junction names use the pattern ``chr_intronStart_intronEnd_strandCode``
    where start/end are 1-based and strandCode is 0/1/2.  Chromosome names that
    themselves contain underscores (e.g. ``chrUn_gl000220``) are handled by
    splitting from the right.
    """
    parts = name.rsplit("_", 3)
    if len(parts) == 4:
        chrom, start_s, end_s, strand_code = parts
        strand = STAR_STRAND_MAP.get(strand_code, ".")
    elif len(parts) == 3:
        chrom, start_s, end_s = parts
        strand = "."
    else:
        return None

    try:
        start_1 = int(start_s)
        end_1 = int(end_s)
    except ValueError:
        return None

    return (chrom, start_1 - 1, end_1, strand)


def load_sample_tissue_map(annot_path: Path) -> Dict[str, str]:
    """Return ``{SAMPID: SMTSD}`` from the GTEx sample-attributes file."""
    mapping: Dict[str, str] = {}
    with open(annot_path, "r") as fh:
        header = fh.readline().strip().split("\t")
        try:
            sid_col = header.index("SAMPID")
            tis_col = header.index("SMTSD")
        except ValueError:
            raise RuntimeError(
                f"Annotations file missing SAMPID or SMTSD columns. Header: {header[:10]}"
            )
        for line in fh:
            fields = line.strip().split("\t")
            if len(fields) > max(sid_col, tis_col):
                mapping[fields[sid_col]] = fields[tis_col]
    return mapping


def tissue_to_dirname(tissue: str) -> str:
    """Sanitize a GTEx tissue subtype name for use as a directory name."""
    out = tissue.replace(" - ", "_").replace(" ", "_")
    for ch in "()/'\"":
        out = out.replace(ch, "")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert GTEx STAR junction GCT to per-sample BED6, filtered by tissue."
    )
    parser.add_argument("--gct", required=True, help="GTEx junction GCT file (may be gzipped)")
    parser.add_argument("--annotations", required=True, help="GTEx sample annotations TSV")
    parser.add_argument("--tissues", required=True,
                        help="Comma-separated tissue names (must match SMTSD column)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--min_count", type=int, default=1,
                        help="Minimum read count to emit a junction (default: 1)")
    parser.add_argument("--max_samples_per_tissue", type=int, default=0,
                        help="Cap number of samples per tissue (0 = unlimited)")
    args = parser.parse_args()

    gct_path = Path(args.gct)
    annot_path = Path(args.annotations)
    outdir = Path(args.outdir)
    requested: Set[str] = {t.strip() for t in args.tissues.split(",") if t.strip()}
    min_count = args.min_count
    max_per_tissue = args.max_samples_per_tissue

    if not gct_path.exists():
        raise SystemExit(f"GCT not found: {gct_path}")
    if not annot_path.exists():
        raise SystemExit(f"Annotations not found: {annot_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Sample-tissue mapping ----
    print(f"[ANNOT] Loading {annot_path.name} ...")
    sample_tissue = load_sample_tissue_map(annot_path)
    available_tissues = sorted(set(sample_tissue.values()))
    print(f"[ANNOT] {len(sample_tissue)} samples across {len(available_tissues)} tissues")

    bad = requested - set(available_tissues)
    if bad:
        print(f"[WARN] Requested tissues not found: {bad}", file=sys.stderr)
        print(f"[INFO] Available: {available_tissues}", file=sys.stderr)

    # ---- 2. Read GCT header ----
    print(f"[GCT]  Opening {gct_path.name} ...")
    opener = gzip.open if str(gct_path).endswith(".gz") else open

    with opener(gct_path, "rt") as gct_fh:
        version_line = gct_fh.readline().strip()
        if not version_line.startswith("#1."):
            print(f"[WARN] Unexpected GCT version: {version_line}")

        dims = gct_fh.readline().strip().split("\t")
        n_junctions, n_samples = int(dims[0]), int(dims[1])
        print(f"[GCT]  {n_junctions:,} junctions x {n_samples:,} samples")

        header_fields = gct_fh.readline().rstrip("\n").split("\t")
        gct_sample_ids = header_fields[2:]

        # ---- 3. Identify columns that match requested tissues ----
        keep_col_indices: List[int] = []
        keep_sample_ids: List[str] = []
        tissue_counts: Dict[str, int] = {}
        for i, sid in enumerate(gct_sample_ids):
            tissue = sample_tissue.get(sid)
            if tissue and tissue in requested:
                n = tissue_counts.get(tissue, 0)
                if max_per_tissue > 0 and n >= max_per_tissue:
                    continue
                tissue_counts[tissue] = n + 1
                keep_col_indices.append(i)
                keep_sample_ids.append(sid)

        if not keep_col_indices:
            print(f"[WARN] No samples matched tissues {requested} — writing empty filelist")
            filelist = outdir / "junction_files.txt"
            filelist.write_text("")
            return
        cap_msg = f" (capped at {max_per_tissue}/tissue)" if max_per_tissue > 0 else ""
        print(f"[FILTER] {len(keep_col_indices):,} samples match{cap_msg}")
        for t, c in sorted(tissue_counts.items()):
            print(f"  {t}: {c} samples")

        # ---- 4. Open per-sample BED file handles ----
        bed_paths: List[Path] = []
        handles: Dict[int, "IO"] = {}
        for pos, (col_idx, sid) in enumerate(zip(keep_col_indices, keep_sample_ids)):
            tissue = sample_tissue[sid]
            tdir = outdir / tissue_to_dirname(tissue)
            tdir.mkdir(parents=True, exist_ok=True)
            safe = sid.replace("/", "_").replace(" ", "_")
            bp = tdir / f"{safe}.juncs.bed"
            bed_paths.append(bp)
            handles[pos] = open(bp, "w")

        # ---- 5. Stream data rows ----
        print(f"[GCT]  Streaming {n_junctions:,} rows ...")
        written = 0
        skipped = 0
        DATA_OFFSET = 2  # Name + Description columns before sample data

        for row_i, line in enumerate(gct_fh):
            if row_i % 100_000 == 0 and row_i > 0:
                print(f"  ... {row_i:,}/{n_junctions:,} rows  ({written:,} entries)")

            fields = line.rstrip("\n").split("\t")
            coords = parse_junction_name(fields[0])
            if coords is None:
                skipped += 1
                continue

            chrom, bed_start, bed_end, strand = coords
            label = f"{chrom}:{bed_start}-{bed_end}:{strand}"

            for pos, col_idx in enumerate(keep_col_indices):
                try:
                    count = int(round(float(fields[col_idx + DATA_OFFSET])))
                except (ValueError, IndexError):
                    continue
                if count >= min_count:
                    handles[pos].write(
                        f"{chrom}\t{bed_start}\t{bed_end}\t{label}\t{count}\t{strand}\n"
                    )
                    written += 1

        for h in handles.values():
            h.close()

    print(f"[GCT]  Done: {written:,} entries across {len(bed_paths)} samples")
    if skipped:
        print(f"[GCT]  Skipped {skipped:,} unparseable junction names")

    # ---- 6. Write junction filelist ----
    filelist = outdir / "junction_files.txt"
    with open(filelist, "w") as fh:
        for p in bed_paths:
            fh.write(str(p) + "\n")

    print(f"[LIST] {filelist}")
    print("[DONE] Per-sample BEDs ready for LeafCutter clustering.")


if __name__ == "__main__":
    main()
