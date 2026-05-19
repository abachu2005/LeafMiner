#!/usr/bin/env bash
# build_supplemented_gtf.sh — concatenate a stock GTF with a LongORF supplement,
# sort by chrom/start, and validate the result.
#
# Usage:
#   bash scripts/build_supplemented_gtf.sh \
#     refs/GRCz12tu.ucsc.gtf \
#     refs/longorf_supplement.grcz12.gtf \
#     refs/GRCz12tu.ucsc.longorf.gtf
set -euo pipefail

STOCK_GTF="${1:?Usage: $0 <stock.gtf> <supplement.gtf> <output.gtf>}"
SUPPLEMENT_GTF="${2:?Usage: $0 <stock.gtf> <supplement.gtf> <output.gtf>}"
OUTPUT_GTF="${3:?Usage: $0 <stock.gtf> <supplement.gtf> <output.gtf>}"

validate_gtf_9_fields() {
    local path="$1"
    local label="$2"
    awk -F'\t' -v label="$label" '
        /^#/ { next }
        NF != 9 {
            printf("ERROR: %s has malformed GTF row at line %d: expected 9 tab fields, found %d\n", label, NR, NF) > "/dev/stderr"
            exit 2
        }
        END {
            if (NR == 0) {
                printf("ERROR: %s is empty\n", label) > "/dev/stderr"
                exit 2
            }
        }
    ' "$path"
}

validate_gtf_attrs() {
    local path="$1"
    local label="$2"
    python3 - "$path" "$label" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
label = sys.argv[2]
attr_re = re.compile(r'([A-Za-z_][A-Za-z0-9_.-]*)\s+"([^"]*)"')
rows = parseable = transcript_rows = transcript_id_rows = 0
with path.open() as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        rows += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise SystemExit(f"ERROR: {label}: expected 9 fields, found {len(fields)}")
        attrs = dict(attr_re.findall(fields[8]))
        if attrs:
            parseable += 1
        if fields[2] in {"transcript", "exon", "CDS"}:
            transcript_rows += 1
            if "transcript_id" in attrs:
                transcript_id_rows += 1
if rows == 0:
    raise SystemExit(f"ERROR: {label}: no feature rows")
if parseable / rows < 0.95:
    raise SystemExit(f"ERROR: {label}: only {parseable / rows:.1%} rows have parseable attributes")
if transcript_rows and transcript_id_rows / transcript_rows < 0.80:
    raise SystemExit(
        f"ERROR: {label}: only {transcript_id_rows / transcript_rows:.1%} transcript/exon/CDS rows have transcript_id"
    )
print(f"[build_supplemented_gtf] {label}: rows={rows} parseable_attrs={parseable}")
PY
}

echo "[build_supplemented_gtf] Stock GTF:      $STOCK_GTF"
echo "[build_supplemented_gtf] Supplement GTF:  $SUPPLEMENT_GTF"
echo "[build_supplemented_gtf] Output GTF:      $OUTPUT_GTF"

if [ ! -s "$STOCK_GTF" ]; then
    echo "ERROR: Stock GTF not found or empty: $STOCK_GTF" >&2
    exit 1
fi
if [ ! -s "$SUPPLEMENT_GTF" ]; then
    echo "ERROR: Supplement GTF not found or empty: $SUPPLEMENT_GTF" >&2
    exit 1
fi

echo "[build_supplemented_gtf] Validating inputs..."
validate_gtf_9_fields "$STOCK_GTF" "stock"
validate_gtf_9_fields "$SUPPLEMENT_GTF" "supplement"
validate_gtf_attrs "$STOCK_GTF" "stock"
validate_gtf_attrs "$SUPPLEMENT_GTF" "supplement"

TMPFILE="${OUTPUT_GTF}.tmp"

echo "[build_supplemented_gtf] Concatenating..."
{
    grep '^#' "$STOCK_GTF" || true
    echo "## LongORF supplement appended by build_supplemented_gtf.sh"
    grep -v '^#' "$STOCK_GTF"
    grep -v '^#' "$SUPPLEMENT_GTF"
} > "$TMPFILE"

echo "[build_supplemented_gtf] Sorting by chrom + start..."
{
    grep '^#' "$TMPFILE"
    grep -v '^#' "$TMPFILE" | sort -k1,1 -k4,4n -k5,5n
} > "$OUTPUT_GTF"
rm -f "$TMPFILE"

echo "[build_supplemented_gtf] Validating output..."
validate_gtf_9_fields "$OUTPUT_GTF" "merged"
validate_gtf_attrs "$OUTPUT_GTF" "merged"

N_STOCK=$(awk '$3=="transcript"' "$STOCK_GTF" | wc -l | tr -d ' ')
N_SUPP=$(awk '$3=="transcript"' "$SUPPLEMENT_GTF" | wc -l | tr -d ' ')
N_MERGED=$(awk '$3=="transcript"' "$OUTPUT_GTF" | wc -l | tr -d ' ')

echo "[build_supplemented_gtf] Transcript counts: stock=$N_STOCK supplement=$N_SUPP merged=$N_MERGED"

EXPECTED=$((N_STOCK + N_SUPP))
if [ "$N_MERGED" -ne "$EXPECTED" ]; then
    echo "WARNING: merged transcript count ($N_MERGED) != stock+supplement ($EXPECTED)" >&2
fi

echo "[build_supplemented_gtf] Done: $OUTPUT_GTF"
