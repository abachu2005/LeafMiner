#!/usr/bin/env bash
# setup_grcz11_longorf_chain.sh — properly swap the UCSC chain so CrossMap can
# lift LongORF supplement from GRCz12tu to GRCz11 (danRer11).
#
# The chain `danRer11ToGCF_049306965.1.over.chain` maps danRer11 -> GRCz12tu.
# We swap it correctly (handling minus-strand queries per UCSC spec) so
# CrossMap can convert FROM NC_* (GRCz12tu) TO chr* (danRer11).
#
# Usage (on Quest):
#   bash scripts/setup_grcz11_longorf_chain.sh \
#     refs/longorf_supplement.grcz12.gtf \
#     refs/longorf_supplement.grcz11.gtf
#
# If liftover yield < 70% the script exits with code 2 (fallback signal).
set -euo pipefail

INPUT_GTF="${1:?Usage: $0 <input_grcz12.gtf> <output_grcz11.gtf>}"
OUTPUT_GTF="${2:?Usage: $0 <input_grcz12.gtf> <output_grcz11.gtf>}"
REFS_DIR="$(dirname "$OUTPUT_GTF")"
CHAIN_DIR="${REFS_DIR}/liftover_chains"
MIN_YIELD="${3:-70}"

mkdir -p "$CHAIN_DIR"

CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/liftOver/danRer11ToGCF_049306965.1.over.chain.gz"
CHAIN_GZ="$CHAIN_DIR/danRer11ToGCF_049306965.1.over.chain.gz"
CHAIN_RAW="$CHAIN_DIR/danRer11ToGCF_049306965.1.over.chain"
SWAPPED_CHAIN="$CHAIN_DIR/GCF_049306965.1_to_danRer11.over.chain"

echo "[liftover] Step 1: Install CrossMap if needed..."
pip install --user CrossMap 2>&1 | tail -2
export PATH="$HOME/.local/bin:$PATH"
CROSSMAP=$(which CrossMap.py 2>/dev/null || which CrossMap)
echo "  CrossMap at: $CROSSMAP"

echo "[liftover] Step 2: Download + decompress chain file..."
if [ ! -s "$CHAIN_GZ" ]; then
    wget -q -O "$CHAIN_GZ" "$CHAIN_URL"
fi
if [ ! -s "$CHAIN_RAW" ]; then
    gunzip -k "$CHAIN_GZ"
fi

echo "[liftover] Step 3: Build chr* -> NC_* mapping..."
ASSEMBLY_REPORT="${REFS_DIR}/GRCz12tu_assembly_report.txt"
if [ ! -s "$ASSEMBLY_REPORT" ]; then
    echo "ERROR: Assembly report not found: $ASSEMBLY_REPORT" >&2
    exit 1
fi
CHR_TO_NC="$CHAIN_DIR/grcz12tu_chr_to_nc.tsv"
awk 'BEGIN{FS="\t"} !/^#/ && NF>=7 {print $1, $7}' "$ASSEMBLY_REPORT" > "$CHR_TO_NC"
echo "  Chrom mapping entries: $(wc -l < "$CHR_TO_NC")"

echo "[liftover] Step 4: Convert supplement GTF from chr* to NC_*..."
INPUT_NC_GTF="$CHAIN_DIR/longorf_supplement.grcz12.nc_names.gtf"
python3 << PYEOF
mapping = {}
with open("$CHR_TO_NC") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
n_mapped = n_skip = 0
with open("$INPUT_GTF") as fin, open("$INPUT_NC_GTF", "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        chrom = line.split("\t", 1)[0]
        if chrom in mapping:
            fout.write(mapping[chrom] + "\t" + line.split("\t", 1)[1])
            n_mapped += 1
        else:
            n_skip += 1
print(f"Mapped {n_mapped} lines, skipped {n_skip} unmapped")
PYEOF

echo "[liftover] Step 5: Swap chain direction with proper minus-strand handling..."
export CHAIN_RAW_PATH="$CHAIN_RAW"
export SWAPPED_PATH="$SWAPPED_CHAIN"
python3 << 'SWAP_EOF'
import sys, os

chain_raw = os.environ["CHAIN_RAW_PATH"]
swapped = os.environ["SWAPPED_PATH"]

def write_chain(fout, hdr, blocks):
    if len(hdr) < 13:
        return
    score = hdr[1]
    tN, tSz, tStr, tSt, tEn = hdr[2], int(hdr[3]), hdr[4], int(hdr[5]), int(hdr[6])
    qN, qSz, qStr, qSt, qEn = hdr[7], int(hdr[8]), hdr[9], int(hdr[10]), int(hdr[11])
    cid = hdr[12]

    if qStr == "+":
        ntN, ntSz, ntStr, ntSt, ntEn = qN, qSz, "+", qSt, qEn
        nqN, nqSz, nqStr, nqSt, nqEn = tN, tSz, "+", tSt, tEn
        nb = []
        for b in blocks:
            p = b.split()
            nb.append(f"{p[0]}\t{p[2]}\t{p[1]}" if len(p) == 3 else b)
    else:
        ntN, ntSz, ntStr = qN, qSz, "+"
        ntSt = qSz - qEn
        ntEn = qSz - qSt
        nqN, nqSz, nqStr = tN, tSz, "-"
        nqSt = tSz - tEn
        nqEn = tSz - tSt
        parsed = []
        for b in blocks:
            p = b.split()
            if len(p) == 3:
                parsed.append((int(p[0]), int(p[1]), int(p[2])))
            elif len(p) == 1:
                parsed.append((int(p[0]),))
        nb = []
        if parsed:
            last = parsed[-1]
            rest = parsed[:-1]
            rev = list(reversed(rest))
            for item in rev:
                if len(item) == 3:
                    nb.append(f"{item[0]}\t{item[2]}\t{item[1]}")
            if len(last) == 1:
                nb.append(str(last[0]))

    fout.write(f"chain {score} {ntN} {ntSz} {ntStr} {ntSt} {ntEn} "
               f"{nqN} {nqSz} {nqStr} {nqSt} {nqEn} {cid}\n")
    for b in nb:
        fout.write(b + "\n")
    fout.write("\n")

with open(chain_raw) as fin, open(swapped, "w") as fout:
    hdr = None
    blks = []
    for line in fin:
        line = line.rstrip("\n")
        if line.startswith("chain"):
            if hdr is not None:
                write_chain(fout, hdr, blks)
            hdr = line.split()
            blks = []
        elif line.strip() == "":
            if hdr is not None:
                write_chain(fout, hdr, blks)
                hdr = None
                blks = []
        else:
            blks.append(line.strip())
    if hdr is not None:
        write_chain(fout, hdr, blks)

print("Chain swap complete")
SWAP_EOF

echo "[liftover] Step 6: Run CrossMap liftover..."
"$CROSSMAP" gff "$SWAPPED_CHAIN" "$INPUT_NC_GTF" "$OUTPUT_GTF" 2>&1 | tail -5

echo "[liftover] Step 7: Validate and compute yield..."
N_IN=$(awk '$3=="transcript"' "$INPUT_GTF" | wc -l)
N_OUT=$(awk '$3=="transcript"' "$OUTPUT_GTF" 2>/dev/null | wc -l || echo 0)

if [ "$N_IN" -gt 0 ]; then
    YIELD=$(python3 -c "print(f'{$N_OUT/$N_IN*100:.1f}')")
else
    YIELD=0
fi

echo ""
echo "=== Liftover Summary ==="
echo "  Input transcripts:   $N_IN"
echo "  Output transcripts:  $N_OUT"
echo "  Yield:               ${YIELD}%"
echo ""

YIELD_INT=$(python3 -c "print(int($N_OUT/$N_IN*100) if $N_IN>0 else 0)")
if [ "$YIELD_INT" -lt "$MIN_YIELD" ]; then
    echo "WARNING: Liftover yield (${YIELD}%) is below ${MIN_YIELD}% threshold!" >&2
    echo "Falling back to stock GRCz11 GTF for Run 3." >&2
    exit 2
fi

echo "[liftover] Done: $OUTPUT_GTF"
