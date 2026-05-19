#!/usr/bin/env python3
"""Rename only column 1 of a GTF using a two-column mapping file.

This intentionally preserves columns 2-9 exactly, especially the attributes
column. It fails fast on malformed non-comment rows because splitting GTF
attributes into extra tab fields can silently corrupt downstream LC2 parsing.
"""

import argparse
from pathlib import Path
from typing import Dict


def load_mapping(path: Path) -> Dict[str, str]:
    mapping = {}  # type: Dict[str, str]
    with path.open() as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                raise SystemExit(f"{path}:{lineno}: expected at least two columns")
            old, new = parts[0], parts[1]
            mapping[old] = new
    if not mapping:
        raise SystemExit(f"No mapping rows found in {path}")
    return mapping


def rename_gtf(in_gtf: Path, mapping: Dict[str, str], out_gtf: Path) -> None:
    rows = comments = renamed = unchanged = missing = 0
    missing_names = {}  # type: Dict[str, int]
    out_gtf.parent.mkdir(parents=True, exist_ok=True)

    tmp = out_gtf.with_suffix(out_gtf.suffix + ".tmp")
    with in_gtf.open() as fin, tmp.open("w") as fout:
        for lineno, line in enumerate(fin, start=1):
            if line.startswith("#"):
                comments += 1
                fout.write(line)
                continue
            if not line.strip():
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SystemExit(
                    f"{in_gtf}:{lineno}: expected 9 tab-delimited GTF fields, found {len(fields)}"
                )

            rows += 1
            old = fields[0]
            new = mapping.get(old)
            if new is None:
                missing += 1
                missing_names[old] = missing_names.get(old, 0) + 1
            elif new != old:
                fields[0] = new
                renamed += 1
            else:
                unchanged += 1
            fout.write("\t".join(fields) + "\n")

    tmp.replace(out_gtf)
    print(f"Input: {in_gtf}")
    print(f"Output: {out_gtf}")
    print(f"Comments: {comments}")
    print(f"Rows: {rows}")
    print(f"Renamed: {renamed}")
    print(f"Unchanged: {unchanged}")
    print(f"Missing mapping rows: {missing}")
    if missing_names:
        top = sorted(missing_names.items(), key=lambda kv: kv[1], reverse=True)[:20]
        print("Top missing seqnames:")
        for name, count in top:
            print(f"  {name}\t{count}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_gtf", type=Path)
    p.add_argument("mapping_tsv", type=Path)
    p.add_argument("output_gtf", type=Path)
    args = p.parse_args()
    rename_gtf(args.input_gtf, load_mapping(args.mapping_tsv), args.output_gtf)


if __name__ == "__main__":
    main()
