#!/usr/bin/env python3
"""
recount3_to_star_sj.py  (no-args edition, with tqdm + progress logging)

- Scans $LEAFCUTTER_ROOT/junctions for recount3 triplets
  (*.UNIQUE.MM.gz, *.UNIQUE.RR.gz, *.UNIQUE.ID.gz) for BRAIN and SPLEEN.
- Converts them into per-sample STAR SJ.out.tab files and writes to:
    $LEAFCUTTER_ROOT/star_sj/<TISSUE>/<sample>.SJ.out.tab

STAR SJ.out.tab columns written:
 1 chr, 2 intronStart(1-based), 3 intronEnd(1-based), 4 strand(0/1/2),
 5 intronMotif=0, 6 annotated=0, 7 uniquely_mapping_reads=count,
 8 multi_mapping_reads=0, 9 max_spliced_alignment_overhang=0
"""

from pathlib import Path
from glob import glob
from collections import defaultdict
import gzip
import pandas as pd
import sys
import re
from tqdm import tqdm

# ------------------ CONSTANT PATHS (edit if you move things) ------------------

import os as _os
_ROOT = Path(_os.environ.get("LEAFCUTTER_ROOT", _os.getcwd())).expanduser()
JUNCTIONS_DIR = _ROOT / "junctions"
OUT_ROOT      = _ROOT / "star_sj"
TISSUES       = ("BRAIN", "SPLEEN")      # auto-detects these in filenames
MIN_COUNT     = 1                        # drop per-sample junctions with count < 1
MAX_SAMPLES   = None                     # or small int for smoke tests

# ------------------ Helpers ------------------

def _ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def _open_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def _to_chr(x: str) -> str:
    x = str(x).strip()
    if not x.lower().startswith("chr"):
        if x and x[0].isdigit():
            return "chr" + x
        if x and (x.startswith("X") or x.startswith("Y") or x.startswith("M")):
            return "chr" + x
    return x

def _strand_to_star(s: str) -> int:
    s = str(s).strip()
    return 1 if s == "+" else 2 if s == "-" else 0

# ------------------ recount3 loaders ------------------

def load_sample_ids(id_path: Path):
    print(f"[IDs] Loading sample IDs from: {id_path}")
    ids = []
    with _open_maybe_gzip(id_path) as fh:
        for line in fh:
            t = line.strip()
            if t:
                ids.append(t)
    if not ids:
        raise RuntimeError(f"No sample IDs found in {id_path}")
    print(f"[IDs] Loaded {len(ids)} sample IDs")
    return ids

def load_row_ranges(rr_path: Path) -> pd.DataFrame:
    print(f"[RR] Loading rowRanges from: {rr_path}")
    # try standard read
    try:
        df = pd.read_csv(rr_path, sep="\t", dtype=str, low_memory=False, comment="#")
    except Exception:
        df = pd.read_csv(rr_path, sep="\t", dtype=str, low_memory=False, comment="#", skiprows=1)

    cols = {c.lower(): c for c in df.columns}
    def _pick(*names):
        for n in names:
            if n in cols: return cols[n]
        return None

    chr_col   = _pick("seqnames","seqname","chrom","chr")
    start_col = _pick("start")
    end_col   = _pick("end")
    strand_col= _pick("strand")

    if not all([chr_col, start_col, end_col, strand_col]):
        print("[RR] Header looked odd; trying headerless heuristic...")
        df2 = pd.read_csv(rr_path, sep="\t", header=None, dtype=str, low_memory=False, comment="#")
        strand_j = None
        for j in range(min(8, df2.shape[1])):
            vals = set(df2[j].dropna().astype(str).str.strip().unique())
            if vals.issubset({"+","-",".","*"}):
                strand_j = j; break
        if strand_j is None or df2.shape[1] < 4:
            raise RuntimeError(f"Couldn't parse rowRanges columns in {rr_path}")
        chr_j   = strand_j - 3
        start_j = strand_j - 2
        end_j   = strand_j - 1
        df = df2[[chr_j,start_j,end_j,strand_j]].copy()
        df.columns = ["seqnames","start","end","strand"]
    else:
        df = df[[chr_col,start_col,end_col,strand_col]].copy()
        df.columns = ["seqnames","start","end","strand"]

    # clean
    before = len(df)
    df["seqnames"] = df["seqnames"].astype(str).map(_to_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    df = df.dropna(subset=["start","end"])
    df = df[(df["start"] > 0) & (df["end"] > 0)]
    dropped = before - len(df)
    if dropped > 0:
        print(f"[RR] Dropped {dropped} malformed junction rows")

    df["star_strand"] = df["strand"].astype(str).map(_strand_to_star).astype(int)
    print(f"[RR] Final usable junction rows: {len(df)}")
    return df.reset_index(drop=True)

def stream_mm_with_meta(mm_path: Path):
    """
    Stream triplets (i, j, x) from Matrix Market with metadata.
    Returns (nrow, ncol, nnz, iterator).
    """
    print(f"[MM] Scanning header: {mm_path}")
    f = _open_maybe_gzip(mm_path)
    header = f.readline()
    if not header.startswith("%%MatrixMarket"):
        f.close()
        raise RuntimeError(f"{mm_path} is not Matrix Market")
    # Skip comments
    line = f.readline()
    while line.startswith("%") or not line.strip():
        line = f.readline()
    nrow, ncol, nnz = map(int, line.strip().split())

    def _iter():
        with _open_maybe_gzip(mm_path) as fh:
            # skip header + comments + dims again
            _ = fh.readline()
            line = fh.readline()
            while line.startswith("%") or not line.strip():
                line = fh.readline()
            # now the triplets
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                i_str, j_str, x_str = line.split()
                i = int(i_str) - 1
                j = int(j_str) - 1
                # recount3 sometimes stores ints as float strings
                x = int(float(x_str))
                if x > 0:
                    yield i, j, x
    return nrow, ncol, nnz, _iter()

# ------------------ Discovery & Conversion ------------------

def find_triplet(jdir: Path, tissue: str):
    pat = re.compile(rf"\.{tissue}\.", re.IGNORECASE)
    def pick(glob_pat):
        cands = [Path(p) for p in glob(str(jdir / glob_pat))]
        cands = [p for p in cands if pat.search(p.name)]
        if not cands:
            raise FileNotFoundError(f"No {glob_pat} for {tissue} in {jdir}")
        cands.sort()
        return cands[-1]   # most recent/lexicographically last
    mm  = pick("*UNIQUE.MM.gz")
    rr  = pick("*UNIQUE.RR.gz")
    idf = pick("*UNIQUE.ID.gz")
    print(f"[DISCOVER] {tissue}:")
    print(f"  MM: {mm.name}\n  RR: {rr.name}\n  ID: {idf.name}")
    return mm, rr, idf

def convert_one_tissue(mm_path: Path, rr_path: Path, id_path: Path, out_dir: Path,
                       min_count: int = MIN_COUNT, max_samples = MAX_SAMPLES):
    # Load RR & IDs
    rr = load_row_ranges(rr_path)
    sample_ids = load_sample_ids(id_path)
    if max_samples is not None:
        print(f"[IDs] Limiting to first {max_samples} samples (debug)")
        sample_ids = sample_ids[:max_samples]

    # Stream MM with a progress bar using nnz
    nrow, ncol, nnz, trip_iter = stream_mm_with_meta(mm_path)
    if max_samples is not None:
        print(f"[MM] Matrix dims: rows={nrow}, cols={ncol} -> limiting to cols={max_samples}")
    else:
        print(f"[MM] Matrix dims: rows={nrow}, cols={ncol}, nnz={nnz}")

    buckets = defaultdict(lambda: defaultdict(int))  # sample -> {row: count}
    skipped_low = 0
    skipped_cols = 0

    for i, j, x in tqdm(trip_iter, total=nnz, desc=f"[MM] Counting (nnz={nnz})", unit="nz"):
        if max_samples is not None and j >= max_samples:
            skipped_cols += 1
            continue
        if x >= min_count:
            buckets[j][i] += x
        else:
            skipped_low += 1

    print(f"[MM] Done streaming. Buckets: {len(buckets)} samples with >={min_count} counts "
          f"(skipped low-count={skipped_low}, skipped_out_of_range_cols={skipped_cols})")

    # Write STAR files with progress over samples
    print(f"[STAR] Writing per-sample SJ.out.tab to: {out_dir}")
    for sidx, sname in tqdm(list(enumerate(sample_ids)), desc="[STAR] Samples", unit="smp"):
        sdir = _ensure_dir(out_dir / sname)
        outp = sdir / f"{sname}.SJ.out.tab"
        rows = buckets.get(sidx, {})
        if not rows:
            # still write an empty file (valid but has no junctions)
            open(outp, "w").close()
            continue
        with open(outp, "w") as fw:
            for ridx, cnt in rows.items():
                if ridx >= len(rr):  # guard
                    continue
                chrom = rr.at[ridx, "seqnames"]
                start = int(rr.at[ridx, "start"])   # 1-based intron start
                end   = int(rr.at[ridx, "end"])     # 1-based intron end
                strand= int(rr.at[ridx, "star_strand"])
                fw.write(f"{chrom}\t{start}\t{end}\t{strand}\t0\t0\t{int(cnt)}\t0\t0\n")
    print(f"[STAR] Completed writing STAR files for {len(sample_ids)} samples.")

def main():
    print("[SETUP] Verifying directories...")
    if not JUNCTIONS_DIR.exists():
        raise SystemExit(f"junctions dir not found: {JUNCTIONS_DIR}")
    _ensure_dir(OUT_ROOT)
    print(f"[SETUP] junctions={JUNCTIONS_DIR}")
    print(f"[SETUP] out_root ={OUT_ROOT}")

    for tissue in TISSUES:
        print(f"\n=== TISSUE: {tissue} ===")
        mm, rr, idf = find_triplet(JUNCTIONS_DIR, tissue)
        t_out = _ensure_dir(OUT_ROOT / tissue.upper())
        convert_one_tissue(mm, rr, idf, t_out, MIN_COUNT, MAX_SAMPLES)

    print("\n[DONE] STAR SJ files are under:", OUT_ROOT)

if __name__ == "__main__":
    main()
