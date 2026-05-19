#!/usr/bin/env python3
"""Convert LongORF_PTC+_fromClustered.tsv to a proper GTF supplement.

Parses exon chains from the ``genomic_coordinates`` column, translates
RefSeq accession chromosome names (NC_133*.1) to UCSC-style chr* using
an NCBI assembly report, and emits gene / transcript / exon / CDS lines.

PTC+ transcripts receive ``transcript_biotype "nonsense_mediated_decay"``
so LeafCutter2's classifier labels junctions as UP (matching the GTF
normalization in ``lc2_pipeline.py``).
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def load_chrom_map(assembly_report: Path) -> Dict[str, str]:
    """Build RefSeq-accession -> chr* map from an NCBI assembly_report.txt."""
    mapping: Dict[str, str] = {}
    with open(assembly_report) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            ucsc_name = cols[0].strip()
            refseq_accn = cols[6].strip()
            if refseq_accn and ucsc_name:
                mapping[refseq_accn] = ucsc_name
    return mapping


def parse_genomic_coordinates(coord_str: str) -> List[Tuple[str, int, int]]:
    """Parse '(NC_133188.1:42668650-42668808);(...)' into [(chrom, start, end), ...]."""
    exons: List[Tuple[str, int, int]] = []
    for m in re.finditer(r"\(([^:]+):(\d+)-(\d+)\)", coord_str):
        chrom = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        exons.append((chrom, start, end))
    return exons


def format_gtf_attributes(attrs: Dict[str, str]) -> str:
    parts = []
    for k, v in attrs.items():
        parts.append(f'{k} "{v}"')
    return "; ".join(parts) + ";"


def convert_longorf_to_gtf(
    longorf_tsv: Path,
    assembly_report: Path,
    output_gtf: Path,
    source_tag: str = "LongORF",
) -> dict:
    """Main conversion. Returns stats dict."""
    chrom_map = load_chrom_map(assembly_report)
    if not chrom_map:
        raise RuntimeError(f"No chrom mappings found in {assembly_report}")
    print(f"Loaded {len(chrom_map)} chrom mappings from {assembly_report}", file=sys.stderr)

    stats = {
        "total_rows": 0,
        "skipped_no_coords": 0,
        "skipped_unmapped_chrom": 0,
        "transcripts_written": 0,
        "exons_written": 0,
        "genes_written": 0,
        "ptc_count": 0,
    }

    gene_spans: Dict[str, Dict] = {}
    transcript_records: List[Dict] = []

    with open(longorf_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            stats["total_rows"] += 1
            tid = row["transcript_id"].strip()
            gene_sym = (row.get("gene_symbol") or tid).strip()
            chrom_raw = row["chrom"].strip()
            strand = row["strand"].strip()
            is_ptc = row.get("is_ptc", "").strip().lower() == "true"
            coord_str = row.get("genomic_coordinates", "").strip()

            if not coord_str:
                stats["skipped_no_coords"] += 1
                continue

            exons = parse_genomic_coordinates(coord_str)
            if not exons:
                stats["skipped_no_coords"] += 1
                continue

            chr_name = chrom_map.get(chrom_raw)
            if chr_name is None:
                stats["skipped_unmapped_chrom"] += 1
                continue

            if is_ptc:
                stats["ptc_count"] += 1

            safe_tid = tid.replace("/", "_")
            gene_id = f"LongORF_{gene_sym}"

            sorted_exons = sorted(exons, key=lambda e: e[1])
            tx_start = sorted_exons[0][1]
            tx_end = sorted_exons[-1][2]

            biotype = "nonsense_mediated_decay" if is_ptc else "protein_coding"

            if gene_id not in gene_spans:
                gene_spans[gene_id] = {
                    "chrom": chr_name, "start": tx_start, "end": tx_end,
                    "strand": strand, "gene_symbol": gene_sym,
                }
            else:
                g = gene_spans[gene_id]
                g["start"] = min(g["start"], tx_start)
                g["end"] = max(g["end"], tx_end)

            transcript_records.append({
                "gene_id": gene_id,
                "transcript_id": safe_tid,
                "chrom": chr_name,
                "strand": strand,
                "tx_start": tx_start,
                "tx_end": tx_end,
                "exons": [(chr_name, e[1], e[2]) for e in sorted_exons],
                "biotype": biotype,
                "gene_symbol": gene_sym,
            })

    with open(output_gtf, "w") as out:
        out.write(f"##gtf-version 2\n")
        out.write(f"##source: longorf_to_gtf.py from {longorf_tsv.name}\n")

        for gid, ginfo in sorted(gene_spans.items()):
            attrs = format_gtf_attributes({
                "gene_id": gid,
                "gene_name": ginfo["gene_symbol"],
                "gene_biotype": "protein_coding",
            })
            out.write(
                f"{ginfo['chrom']}\t{source_tag}\tgene\t{ginfo['start']}\t{ginfo['end']}\t"
                f".\t{ginfo['strand']}\t.\t{attrs}\n"
            )
            stats["genes_written"] += 1

        for rec in transcript_records:
            tx_attrs = format_gtf_attributes({
                "gene_id": rec["gene_id"],
                "transcript_id": rec["transcript_id"],
                "gene_name": rec["gene_symbol"],
                "gene_biotype": "protein_coding",
                "transcript_biotype": rec["biotype"],
            })
            out.write(
                f"{rec['chrom']}\t{source_tag}\ttranscript\t{rec['tx_start']}\t{rec['tx_end']}\t"
                f".\t{rec['strand']}\t.\t{tx_attrs}\n"
            )
            stats["transcripts_written"] += 1

            for i, (echrom, estart, eend) in enumerate(rec["exons"], 1):
                exon_attrs = format_gtf_attributes({
                    "gene_id": rec["gene_id"],
                    "transcript_id": rec["transcript_id"],
                    "exon_number": str(i),
                    "gene_name": rec["gene_symbol"],
                    "gene_biotype": "protein_coding",
                    "transcript_biotype": rec["biotype"],
                })
                out.write(
                    f"{echrom}\t{source_tag}\texon\t{estart}\t{eend}\t"
                    f".\t{rec['strand']}\t.\t{exon_attrs}\n"
                )
                stats["exons_written"] += 1

    print(f"\nConversion stats:", file=sys.stderr)
    for k, v in stats.items():
        print(f"  {k}: {v}", file=sys.stderr)
    print(f"\nWrote {output_gtf}", file=sys.stderr)
    return stats


def main():
    p = argparse.ArgumentParser(description="Convert LongORF TSV to GTF supplement")
    p.add_argument("--longorf-tsv", required=True, help="Path to LongORF_PTC+_fromClustered.tsv")
    p.add_argument("--assembly-report", required=True, help="NCBI assembly_report.txt with RefSeq->chr mapping")
    p.add_argument("--output-gtf", required=True, help="Output GTF file path")
    p.add_argument("--source-tag", default="LongORF", help="GTF source field value")
    args = p.parse_args()

    convert_longorf_to_gtf(
        longorf_tsv=Path(args.longorf_tsv),
        assembly_report=Path(args.assembly_report),
        output_gtf=Path(args.output_gtf),
        source_tag=args.source_tag,
    )


if __name__ == "__main__":
    main()
