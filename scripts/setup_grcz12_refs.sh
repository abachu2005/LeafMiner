#!/usr/bin/env bash
#
# Download GRCz12tu (zebrafish) genome + GTF from NCBI and convert
# chromosome names from RefSeq accessions to UCSC-style chr* names.
#
# Usage:  bash scripts/setup_grcz12_refs.sh [DEST_DIR]
#         DEST_DIR defaults to ./refs
#
# Requires: ncbi-datasets-cli (module load ncbi-datasets-cli/v2 on Quest)

set -euo pipefail

ACCESSION="GCF_049306965.1"
DEST="${1:-$(dirname "$0")/../refs}"
DEST="$(cd "$DEST" && pwd)"
WORK="$DEST/.grcz12tu_tmp"

FASTA_OUT="$DEST/GRCz12tu.fa"
GTF_OUT="$DEST/GRCz12tu.ucsc.gtf"
REPORT_OUT="$DEST/GRCz12tu_assembly_report.txt"
MAP_FILE="$WORK/chr_map.tsv"

echo "=== GRCz12tu Reference Setup ==="
echo "Accession : $ACCESSION"
echo "Dest dir  : $DEST"
echo ""

if [ -f "$FASTA_OUT" ] && [ -f "$GTF_OUT" ]; then
    echo "Output files already exist:"
    echo "  $FASTA_OUT"
    echo "  $GTF_OUT"
    echo "Delete them first if you want to re-download."
    exit 0
fi

mkdir -p "$WORK"

# ── 1. Download via NCBI datasets CLI ────────────────────────────────
echo "[1/5] Downloading genome + GTF from NCBI ..."
if ! command -v datasets &>/dev/null; then
    echo "ERROR: 'datasets' CLI not found. On Quest run: module load ncbi-datasets-cli/v2"
    exit 1
fi

datasets download genome accession "$ACCESSION" \
    --include genome,gtf \
    --filename "$WORK/ncbi_dataset.zip"

echo "[1/5] Extracting ..."
unzip -o "$WORK/ncbi_dataset.zip" -d "$WORK/ncbi_extract"

RAW_FASTA="$(find "$WORK/ncbi_extract" -name '*.fna' | head -1)"
RAW_GTF="$(find "$WORK/ncbi_extract" -name '*.gtf' | head -1)"

if [ -z "$RAW_FASTA" ] || [ -z "$RAW_GTF" ]; then
    echo "ERROR: Could not find FASTA or GTF in downloaded package."
    ls -R "$WORK/ncbi_extract"
    exit 1
fi
echo "  FASTA: $RAW_FASTA"
echo "  GTF  : $RAW_GTF"

# ── 2. Download assembly report ──────────────────────────────────────
echo "[2/5] Downloading assembly report ..."
REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/306/965/GCF_049306965.1_GRCz12tu/GCF_049306965.1_GRCz12tu_assembly_report.txt"
wget -q -O "$WORK/assembly_report.txt" "$REPORT_URL" || {
    echo "Direct FTP download failed; trying datasets ..."
    datasets summary genome accession "$ACCESSION" --as-json-lines > "$WORK/summary.jsonl" 2>/dev/null || true
    REPORT_URL_ALT="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/306/965/${ACCESSION}_GRCz12tu/${ACCESSION}_GRCz12tu_assembly_report.txt"
    wget -q -O "$WORK/assembly_report.txt" "$REPORT_URL_ALT"
}
cp "$WORK/assembly_report.txt" "$REPORT_OUT"

# ── 3. Build chromosome name mapping ────────────────────────────────
echo "[3/5] Building chromosome name mapping ..."

# assembly_report.txt is tab-delimited with columns:
#   1: Sequence-Name   (e.g. "1", "2", ... "25", "MT")
#   2: Sequence-Role   (assembled-molecule, unplaced-scaffold, ...)
#   3: Assigned-Molecule
#   4: Assigned-Molecule-loc/type
#   5: GenBank-Accn
#   6: Relationship
#   7: RefSeq-Accn     (e.g. NC_007112.8)
#   8: Assembly-Unit
#   9: Sequence-Length
#   10: UCSC-style-name (may be "na")
#
# We map RefSeq-Accn -> chr{Sequence-Name}, with MT -> chrM.

awk -F'\t' '
    /^#/ { next }
    $2 == "assembled-molecule" {
        refseq = $7
        name = $1
        # Sequence-Name may already have chr prefix (GRCz12tu uses chr1..chr25)
        if (name == "MT") name = "chrM"
        else if (name !~ /^chr/) name = "chr" name
        print refseq "\t" name
    }
' "$WORK/assembly_report.txt" > "$MAP_FILE"

# Also include unplaced scaffolds (keep their RefSeq accession as-is,
# prefixed with "chrUn_" to follow UCSC convention)
awk -F'\t' '
    /^#/ { next }
    $2 == "unplaced-scaffold" {
        refseq = $7
        safe = refseq
        gsub(/\./, "v", safe)
        print refseq "\tchrUn_" safe
    }
' "$WORK/assembly_report.txt" >> "$MAP_FILE"

N_MAPPED=$(wc -l < "$MAP_FILE")
echo "  Mapped $N_MAPPED sequences"
echo "  First 5 entries:"
head -5 "$MAP_FILE" | sed 's/^/    /'

# ── 4. Rename chromosomes in FASTA and GTF ───────────────────────────
echo "[4/5] Renaming chromosomes in FASTA ..."

# Build a sed script from the mapping
SED_SCRIPT="$WORK/rename.sed"
> "$SED_SCRIPT"
while IFS=$'\t' read -r old new; do
    # Escape dots in the accession for sed regex
    old_esc=$(echo "$old" | sed 's/\./\\./g')
    echo "s/^>${old_esc}\\b/>${new}/" >> "$SED_SCRIPT"
    echo "s/^${old_esc}\\t/${new}\\t/" >> "$SED_SCRIPT"
done < "$MAP_FILE"

sed -f "$SED_SCRIPT" "$RAW_FASTA" > "$FASTA_OUT"

echo "  Renaming chromosomes in GTF ..."
# For GTF, only the first column (seqname) needs renaming.
# Use awk for precise first-column replacement.
awk -F'\t' -v OFS='\t' '
    BEGIN {
        while ((getline line < "'"$MAP_FILE"'") > 0) {
            split(line, a, "\t")
            map[a[1]] = a[2]
        }
    }
    /^#/ { print; next }
    {
        if ($1 in map) $1 = map[$1]
        print
    }
' "$RAW_GTF" > "$GTF_OUT"

# ── 5. Verify ────────────────────────────────────────────────────────
echo "[5/5] Verifying ..."

FASTA_CHRS=$(grep '^>' "$FASTA_OUT" | sed 's/^>//' | cut -d' ' -f1 | sort)
GTF_CHRS=$(awk -F'\t' '!/^#/{print $1}' "$GTF_OUT" | sort -u)

FASTA_N=$(echo "$FASTA_CHRS" | wc -l)
GTF_N=$(echo "$GTF_CHRS" | wc -l)

echo "  FASTA sequences: $FASTA_N"
echo "  GTF chromosomes: $GTF_N"

MISSING=$(comm -23 <(echo "$GTF_CHRS") <(echo "$FASTA_CHRS"))
if [ -z "$MISSING" ]; then
    echo "  OK: All GTF chromosomes found in FASTA."
else
    echo "  WARNING: GTF references chromosomes not in FASTA:"
    echo "$MISSING" | head -10 | sed 's/^/    /'
fi

echo ""
echo "  Assembled chromosomes in FASTA:"
grep '^>chr[0-9M]' "$FASTA_OUT" | head -30 | sed 's/^/    /'

# Clean up temp files
rm -rf "$WORK"

echo ""
echo "=== Done ==="
echo "  FASTA : $FASTA_OUT"
echo "  GTF   : $GTF_OUT"
echo "  Report: $REPORT_OUT"
echo ""
echo "Next: The pipeline will auto-build the STAR index on first run."
