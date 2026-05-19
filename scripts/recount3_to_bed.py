#!/usr/bin/env python3
"""
recount3_to_bed.py  (no-args, with tqdm + robust RR parsing + progress logs)

Convert recount3 junction triplets (MM/RR/ID) directly into per-sample BED6 files
compatible with LeafCutter/LeafCutter2—no STAR step.

Inputs (auto-discovered by tissue in $LEAFCUTTER_ROOT/junctions):
  *.UNIQUE.MM.gz  (Matrix Market, rows=junctions, cols=samples)
  *.UNIQUE.RR.gz  (rowRanges: coordinates + strand)
  *.UNIQUE.ID.gz  (sample IDs, one per line)

Outputs:
  $LEAFCUTTER_ROOT/junctions_bed/<TISSUE>/<sample>.juncs.bed
  $LEAFCUTTER_ROOT/out/junction_files_<TISSUE>.txt  (list of BEDs)

BED6 columns written:
  chrom  start(0-based)  end  name(junction_id)  score(count)  strand(+/-/.)
Notes:
- recount3 RR are 1-based, closed; BED is 0-based, half-open → start_bed = start_rr - 1; end_bed = end_rr
- We stream the sparse MM and bucket counts per sample to avoid dense memory use.
"""

from pathlib import Path
from glob import glob
from collections import defaultdict
import gzip, re, sys
import pandas as pd
from tqdm import tqdm

# ---------- Paths (override with LEAFCUTTER_ROOT or --root) ----------
import os as _os
ROOT           = Path(_os.environ.get("LEAFCUTTER_ROOT", _os.getcwd())).expanduser()
JUNCTIONS_DIR  = ROOT / "junctions"
BED_ROOT       = ROOT / "junctions_bed"
OUT_DIR        = ROOT / "out"
TISSUES        = ("BRAIN", "SPLEEN")
MIN_COUNT      = 1         # drop per-sample junctions with count < 1
MAX_SAMPLES    = None      # set small int to smoke-test (e.g., 10)

# ---------- Helpers ----------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def open_maybe_gzip(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")

def to_chr(x: str) -> str:
    x = str(x).strip()
    if not x.lower().startswith("chr"):
        if x and (x[0].isdigit() or x in ("X","Y","M")):
            return "chr" + x
    return x

def strand_clean(s: str) -> str:
    s = str(s).strip()
    if s in {"+","-"}:
        return s
    return "."

# ---------- Loaders ----------
def load_sample_ids(id_path: Path):
    print(f"[IDs] {id_path.name}")
    ids = []
    with open_maybe_gzip(id_path) as fh:
        for line in fh:
            t = line.strip()
            if t:
                ids.append(t)
    if not ids:
        raise RuntimeError(f"No sample IDs found in {id_path}")
    print(f"[IDs] loaded {len(ids)} samples")
    return ids

def load_row_ranges(rr_path: Path) -> pd.DataFrame:
    """
    Robust RR reader.
    Accepts headers like: seqnames | start | end | strand
    ...or Bioconductor-style: seqnames | ranges.start | ranges.end | strand
    ...or other variants (chrom/chr/seqname). Falls back to headerless heuristic.
    """
    print(f"[RR] {rr_path.name}")

    def _try_headered_read(skiprows=None):
        # don't treat '#' as comments; some RR put useful headers there
        return pd.read_csv(
            rr_path,
            sep="\t",
            dtype=str,
            low_memory=False,
            comment=None,
            skiprows=skiprows
        )

    # 1) Try normal header read; if it fails, try skipping one line
    try_orders = [None, 1]
    last_err = None
    df = None
    for sk in try_orders:
        try:
            df = _try_headered_read(skiprows=sk)
            break
        except Exception as e:
            last_err = e
            df = None
    if df is None:
        raise RuntimeError(f"[RR] Could not read {rr_path}: {last_err}")

    # 2) Regex-based header matching (lowercased names)
    cols = {c.lower(): c for c in df.columns}
    lc_names = list(cols.keys())

    def pick(regexes, exclude_regexes=()):
        import re as _re
        for key in lc_names:
            if any(_re.search(rx, key) for rx in regexes) and not any(_re.search(erx, key) for erx in exclude_regexes):
                return cols[key]
        return None

    chr_col   = pick([r"^seqnames?$", r"^seqname$", r"^chrom$", r"^chr$", r"seqname", r"chrom"])
    start_col = pick([r"(^|[\.\:_])start($|[\.\:_])", r"^ranges\.start$"])
    end_col   = pick([r"(^|[\.\:_])end($|[\.\:_])",   r"^ranges\.end$"], exclude_regexes=[r"friend"])
    strand_col= pick([r"^strand$"])

    have_header = all([chr_col, start_col, end_col, strand_col])

    # 3) If header mapping failed, attempt headerless heuristic
    if not have_header:
        print("[RR] Header looked non-standard; trying headerless heuristic...")
        df2 = pd.read_csv(rr_path, sep="\t", header=None, dtype=str, low_memory=False, comment=None)
        strand_j = None
        for j in range(min(12, df2.shape[1])):
            vals = set(df2[j].dropna().astype(str).str.strip().unique())
            if vals.issubset({"+","-",".","*"}):
                strand_j = j
                break
        if strand_j is None or df2.shape[1] < 4:
            # tiny preview to help debug
            print("[RR] Preview of first 5 non-empty lines:")
            with open_maybe_gzip(rr_path) as fh:
                shown = 0
                for line in fh:
                    t = line.strip()
                    if not t:
                        continue
                    print("   ", t[:200])
                    shown += 1
                    if shown >= 5:
                        break
            raise RuntimeError(f"Could not parse {rr_path} (no obvious strand column found).")

        chr_j   = max(0, strand_j - 3)
        start_j = max(0, strand_j - 2)
        end_j   = max(0, strand_j - 1)
        df = df2[[chr_j, start_j, end_j, strand_j]].copy()
        df.columns = ["seqnames","start","end","strand"]
    else:
        df = df[[chr_col, start_col, end_col, strand_col]].copy()
        df.columns = ["seqnames","start","end","strand"]

    # 4) Clean / coerce
    before = len(df)
    df["seqnames"] = df["seqnames"].astype(str).map(to_chr)
    df["start"]    = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"]      = pd.to_numeric(df["end"],   errors="coerce").astype("Int64")
    df["strand"]   = df["strand"].astype(str).str.strip().map(strand_clean)
    df = df.dropna(subset=["start","end"])
    df = df[(df["start"] > 0) & (df["end"] > 0)]
    dropped = before - len(df)
    if dropped:
        print(f"[RR] Dropped {dropped} malformed junction rows")

    # BED conversion (recount3 is 1-based closed; BED is 0-based half-open)
    df["bed_start"] = df["start"].astype(int) - 1
    df["bed_end"]   = df["end"].astype(int)
    return df.reset_index(drop=True)

def stream_mm(mm_path: Path):
    # read header/dims
    with open_maybe_gzip(mm_path) as fh:
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise RuntimeError(f"{mm_path} not Matrix Market")
        line = fh.readline()
        while line.startswith("%") or not line.strip():
            line = fh.readline()
        nrow, ncol, nnz = map(int, line.strip().split())

    # iterator over triplets
    def it():
        with open_maybe_gzip(mm_path) as f2:
            _ = f2.readline()                               # header
            line = f2.readline()
            while line.startswith("%") or not line.strip(): # comments
                line = f2.readline()
            for line in f2:
                line = line.strip()
                if not line:
                    continue
                i_str, j_str, x_str = line.split()
                i = int(i_str) - 1
                j = int(j_str) - 1
                x = int(float(x_str))
                if x > 0:
                    yield i, j, x
    return nrow, ncol, nnz, it()

# ---------- Discovery ----------
def find_triplet(jdir: Path, tissue: str):
    pat = re.compile(rf"\.{tissue}\.", re.IGNORECASE)
    def pick(glob_pat):
        cands = [Path(p) for p in glob(str(jdir / glob_pat))]
        cands = [p for p in cands if pat.search(p.name)]
        if not cands:
            raise FileNotFoundError(f"No {glob_pat} for {tissue} in {jdir}")
        cands.sort()
        return cands[-1]
    mm  = pick("*UNIQUE.MM.gz")
    rr  = pick("*UNIQUE.RR.gz")
    idf = pick("*UNIQUE.ID.gz")
    print(f"[DISCOVER] {tissue}\n  MM: {mm.name}\n  RR: {rr.name}\n  ID: {idf.name}")
    return mm, rr, idf

# ---------- Conversion ----------
def write_beds_for_tissue(tissue: str):
    mm, rr, idf = find_triplet(JUNCTIONS_DIR, tissue)
    rr_df = load_row_ranges(rr)
    sample_ids = load_sample_ids(idf)
    if MAX_SAMPLES is not None:
        sample_ids = sample_ids[:MAX_SAMPLES]
        print(f"[IDs] limiting to first {len(sample_ids)} samples")

    nrow, ncol, nnz, nz_iter = stream_mm(mm)
    print(f"[MM] rows={nrow} cols={ncol} nnz={nnz}")

    # buckets: sample_index -> { row_index: count }
    buckets = defaultdict(lambda: defaultdict(int))
    skipped_low = 0
    skipped_col = 0

    for i, j, x in tqdm(nz_iter, total=nnz, unit="nz", desc=f"[MM] {tissue} nonzeros"):
        if MAX_SAMPLES is not None and j >= MAX_SAMPLES:
            skipped_col += 1
            continue
        if x >= MIN_COUNT:
            buckets[j][i] += x
        else:
            skipped_low += 1

    print(f"[MM] done. kept_samples={len(buckets)} skipped_low={skipped_low} skipped_cols_out_of_range={skipped_col}")

    # Write one BED per sample
    tissue_dir = ensure_dir(BED_ROOT / tissue.upper())
    bed_list = []

    for sidx, sname in tqdm(list(enumerate(sample_ids)), unit="smp", desc=f"[BED] {tissue} samples"):
        rows = buckets.get(sidx, {})
        outp = tissue_dir / f"{sname}.juncs.bed"
        bed_list.append(outp)
        with open(outp, "w") as fw:
            # BED6: chrom  start(0-based)  end  name  score  strand
            for ridx, cnt in rows.items():
                if ridx >= len(rr_df):  # guard
                    continue
                chrom  = rr_df.at[ridx, "seqnames"]
                start  = int(rr_df.at[ridx, "bed_start"])
                end    = int(rr_df.at[ridx, "bed_end"])
                strand = rr_df.at[ridx, "strand"]
                name   = f"{chrom}:{start}-{end}:{strand}"  # junction ID (human-readable)
                fw.write(f"{chrom}\t{start}\t{end}\t{name}\t{int(cnt)}\t{strand}\n")

    # Write a junction_files list for LeafCutter
    ensure_dir(OUT_DIR)
    list_path = OUT_DIR / f"junction_files_{tissue.upper()}.txt"
    with open(list_path, "w") as fh:
        for p in bed_list:
            fh.write(str(p) + "\n")
    print(f"[BED] wrote {len(bed_list)} BEDs  →  {tissue_dir}")
    print(f"[LIST] {list_path} (for leafcutter_cluster_regtools.py)")

# ---------- Main ----------
def main():
    print("[SETUP] ROOT:", ROOT)
    if not JUNCTIONS_DIR.exists():
        raise SystemExit(f"junctions dir not found: {JUNCTIONS_DIR}")
    ensure_dir(BED_ROOT); ensure_dir(OUT_DIR)

    for tissue in TISSUES:
        print(f"\n=== {tissue} ===")
        write_beds_for_tissue(tissue)

    print("\n[DONE] Per-sample BEDs are under:", BED_ROOT)
    print("[NEXT] Run LeafCutter clustering with the written junction_files_<TISSUE>.txt lists.")

if __name__ == "__main__":
    main()
