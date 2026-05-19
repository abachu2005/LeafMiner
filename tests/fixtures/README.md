# Test fixtures

Tiny, synthetic inputs that let the pipeline + setup wizard run end-to-end in seconds without downloading real genome data.

| File | Format | Purpose |
|---|---|---|
| `sample.SJ.out.tab` | STAR splice-junction tab | Used by `bin/leafcutter-setup` for the smoke test and by `tests/smoke_validate_gtf.sh` examples |

STAR `SJ.out.tab` columns (1-based): chrom, intron_start, intron_end, strand (0=undefined,1=+,2=-), motif (0=non-canonical,1=GT/AG,2=CT/AC,3=GC/AG,4=CT/GC,5=AT/AC,6=GT/AT), annotated (0/1), unique-reads, multi-reads, max-overhang.
