# Run 3 Recovery Notes

Started: 2026-04-28 15:17:29 CDT

## Why the live Run 3 was stopped

The live Run 3 Stage 2 job (`6687898`) was using `refs/Danio_rerio.GRCz11.115.ucsc.longorf.gtf`.
That GTF was malformed: the attributes column had been split into extra tab-delimited fields.
LeafCutter2 therefore parsed field 9 as only `gene_id`, lost all gene/transcript attributes, and
classified all observed introns under gene name `unknown`.

Observed bad-run evidence before cancellation:

- Stage 2 job: `6687898`, `zf_s2_run3`, elapsed about 5h25m.
- Dependent Stage 3 job: `6687899`, `zf_s3_run3`.
- Partial classification file: `lc2/clustering/zebrafish_run_lc2_junction_classifications.txt`.
- Partial file had about 1.23M classification rows, all with gene name `unknown`, and had only reached early chromosomes.

Action taken:

- Issued `scancel 6687898 6687899` on Quest.
- Partial outputs are intentionally preserved in the original run directory for debugging, but must not be treated as final deliverables.

## Fixed rerun outcome

Completed: 2026-04-28 16:00:31 CDT

- Rebuilt `refs/Danio_rerio.GRCz11.115.ucsc.gtf` with a column-1-only seqname conversion and validated 9-field, parseable GTF attributes.
- Kept Run 3 as the conservative stock GRCz11 fallback because the LongORF liftover yield was only 16.5%.
- Ran LC2 smoke on 5 SJ files successfully: 3,837 classification rows, 0 `unknown` genes.
- Submitted fixed Run 3 Stage 2/3 rerun in `/projects/p52853/iis1026/Leaf_Cutter/jobs/8e56d39b-2e33-4cd1-943c-0b7c3b6bb1c8`.
- Stage 2 job `6717964` completed in 4m48s with 39,468 classification rows and 0 `unknown` genes.
- Stage 3 job `6717965` completed in 3m48s.

Important limitation:

- The POISEN PE list remains in GRCz12tu coordinates. Against GRCz11 STAR junctions, no PE events passed coverage filters, so Run 3 is a validated GRCz11 LC2/assembly control rather than a fully comparable GRCz11 PE-inclusion test.
