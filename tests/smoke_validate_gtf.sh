#!/usr/bin/env bash
set -euo pipefail

python3 scripts/validate_gtf_for_lc2.py tests/fixtures/valid_small.gtf \
  --gene-name-tag gene_name \
  --transcript-name-tag transcript_id \
  --gene-type-tag gene_biotype \
  --transcript-type-tag transcript_biotype \
  --min-tag-frac 0.05

if python3 scripts/validate_gtf_for_lc2.py tests/fixtures/malformed_split_attributes.gtf \
  --gene-name-tag gene_name \
  --transcript-name-tag transcript_id \
  --gene-type-tag gene_biotype \
  --transcript-type-tag transcript_biotype \
  --min-tag-frac 0.05 >/tmp/malformed_gtf.err 2>&1; then
  echo "Expected malformed_split_attributes.gtf to fail" >&2
  exit 1
fi
grep -q "expected 9" /tmp/malformed_gtf.err

if python3 scripts/validate_gtf_for_lc2.py tests/fixtures/missing_required_attrs.gtf \
  --gene-name-tag gene_name \
  --transcript-name-tag transcript_id \
  --gene-type-tag gene_biotype \
  --transcript-type-tag transcript_biotype \
  --min-tag-frac 0.05 >/tmp/missing_attrs.err 2>&1; then
  echo "Expected missing_required_attrs.gtf to fail" >&2
  exit 1
fi

echo "GTF validation smoke tests passed"
