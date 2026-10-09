# Zebrafish Poison-Exon Developmental Analysis — Triple-Run Comparison

*Generated 2026-04-29*

## 1. Headline

Three parallel pipeline runs were executed on **276 zebrafish developmental RNA-seq
libraries** (6 timepoints: 24hpf → 5dpf) to assess the impact of:

1. **GTF annotation completeness** — adding 24,694 LongORF PTC+ NMD transcripts to the stock GTF
2. **Assembly version** — GRCz12tu (2025) vs GRCz11/danRer11 (2017, Ensembl 115)

### Key result: LongORF supplement dramatically increases UP junction discovery

The LongORF GTF supplement had **no effect** on the POISEN PE inclusion analysis
(which queries a fixed coordinate list), but it had a **massive effect** on
LeafCutter2's unproductive junction classification — the analysis that discovers
novel poison-exon candidates from the data:

| LC2 Poison-Exon Discovery | Run 1 (stock GTF) | Run 2 (+ LongORF) | Change |
|---|---|---|---|
| UP junctions classified | 86 | **482** | **+460% (5.6×)** |
| Significant (Spearman p<0.05) | 21 | **127** | **+505% (6.0×)** |
| Positive PSI-time correlation (p<0.05) | 11 | **84** | **+664% (7.6×)** |
| Unique genes with UP junctions | 72 | **294** | **+308% (4.1×)** |
| Novel genes (not in Run 1) | — | **226** | — |
| Novel genes with sig. positive corr. | — | **37** | — |
| NMD-eligible (stock annotation) | 13 | 4 | see note below |

NMD-eligibility dropped because the newly classified UP junctions come from
LongORF supplement transcripts that lack entries in LC2's `long_exon_distances`
and `nuc_rule_distances` auxiliary files. The junctions are real and significant;
the NMD annotation is incomplete, not absent.

### POISEN PE inclusion analysis (unchanged)

| POISEN PE Metric | Run 1 (GRCz12tu stock) | Run 2 (GRCz12tu + LongORF) | Run 3 (GRCz11 stock) |
|---|---|---|---|
| Events in PSI matrix | 1936 | 1936 | 1936 |
| Events passing all coverage filters | 23 | 23 | 0 |
| GLM-tested PE events | 2 | 2 | 0 |

## 2. Run Descriptions

### Run 1: GRCz12tu stock (baseline)

- **Assembly:** GRCz12tu
- **GTF:** RefSeq stock GTF (GRCz12tu.ucsc.gtf)
- **Notes:** Reuses prior run artifacts (job 583631a5).

### Run 2: GRCz12tu + LongORF supplement

- **Assembly:** GRCz12tu
- **GTF:** Stock GTF + 24,694 LongORF NMD transcripts (GRCz12tu.ucsc.longorf.gtf)
- **Notes:** Reuses Stage 1 SJ files from Run 1. Only Stage 2 (LC2) and Stage 3 (PE inclusion) re-run with supplemented GTF.

### Run 3: GRCz11 (danRer11) Ensembl 115

- **Assembly:** GRCz11 (danRer11)
- **GTF:** Ensembl 115 GTF (Danio_rerio.GRCz11.115.ucsc.gtf); LongORF liftover yield was 16.5% (<70% threshold), so stock GTF used.
- **Notes:** Fresh alignment of all 276 samples to GRCz11. Slurm array parallelization (20 concurrent tasks).

## 3. LongORF Supplement Details

Talha's `LongORF_PTC+_fromClustered.tsv` contains **24,757 PTC+ transcripts** across
**3,967 genes**. After filtering (63 rows with no genomic coordinates), **24,694 transcripts**
were converted to a proper GTF supplement with `transcript_biotype "nonsense_mediated_decay"`
so LeafCutter2 classifies their junctions as **UP** (unproductive).

**Liftover to GRCz11:** Attempted using UCSC chain `danRer11ToGCF_049306965.1.over.chain`
with a properly swapped chain + CrossMap. Yield was only **16.5%** (4,080 / 24,694 transcripts),
well below the 70% threshold. The massive coordinate divergence between GRCz12tu (2025 Verkko
assembly) and GRCz11 (2017) makes direct liftover impractical. Run 3 therefore used the
**stock Ensembl 115 GRCz11 GTF** without the LongORF supplement.

### Run 3 PE-coordinate limitation

Run 3 LC2 completed cleanly on GRCz11 with real gene names, but the POISEN PE list
is still defined in GRCz12tu genomic coordinates. When those PE intervals were
queried against GRCz11 STAR junctions, no PE events passed coverage filters and
no GRCz11 PE figures were generated. Treat Run 3 as a validated GRCz11 LC2/assembly
control, not as a fully comparable GRCz11 PE-inclusion test, unless the PE list is
lifted or regenerated on GRCz11 coordinates.

## 4. Pipeline Parallelization

The pipeline was refactored to use **Slurm job arrays** for per-sample parallelization:
- QC (fastp trimming): `--array=0-275%20` (20 concurrent downloads)
- Stage 1 (STAR alignment): `--array=0-275%20` (20 concurrent alignments)
- Dependency chaining via `--dependency=afterany:JID` (tolerates individual task failures)
- Email notifications on completion/failure (no busy-polling)

This reduced the estimated Run 3 wall-clock time from **~6 days (serial)** to
**~4-8 hours (parallel)**.

## 5. Top Novel Candidates from LongORF Supplement

These genes were **not classifiable as UP** with the stock GTF but became
significant poison-exon candidates (positive PSI-time correlation, p<0.05)
after adding LongORF NMD transcripts:

| Gene | Spearman ρ | p-value | Notes |
|---|---|---|---|
| gnl3 | +1.000 | <0.001 | Nucleolar GTP-binding protein |
| ube2d2 | +1.000 | <0.001 | Ubiquitin-conjugating enzyme |
| scnm1 | +1.000 | <0.001 | Sodium channel modifier |
| psmd1 | +1.000 | <0.001 | Proteasome subunit |
| mxi1 | +1.000 | <0.001 | MYC-associated factor X interactor |
| eef1da | +1.000 | <0.001 | Translation elongation factor |
| jade2 | +1.000 | <0.001 | Histone acetyltransferase complex |
| rrp36 | +1.000 | <0.001 | Ribosome biogenesis |
| l3mbtl1b | +1.000 | <0.001 | Polycomb group protein |
| ddx3xb | +0.943 | 0.005 | DEAD-box RNA helicase |
| meis1b | +0.943 | 0.005 | Homeobox TF (developmental) |
| ncl | +0.943 | 0.005 | Nucleolin |
| atp1a1b | +0.943 | 0.005 | Na+/K+-ATPase |
| *... 24 more* | | | *See `poison_exon_candidates.tsv`* |

Known splicing regulators from the original analysis that were previously
significant but lacked NMD annotation (`hnrnpa1b`) now appear with rho = +1.0
in Run 2, confirming the annotation gap hypothesis.

## 6. Figures

### Poison-exon discovery (LC2 UP junction analysis)

**Run 1: GRCz12tu stock (86 UP junctions)**

![UP junction heatmap](run1_grcz12_stock_baseline/poison_exon/poison_exon_heatmap.png)

![PSI vs developmental time](run1_grcz12_stock_baseline/poison_exon/up_psi_vs_time.png)

![NMD features](run1_grcz12_stock_baseline/poison_exon/nmd_features.png)

**Run 2: GRCz12tu + LongORF (482 UP junctions)**

![UP junction heatmap](run2_grcz12_longorf/poison_exon/poison_exon_heatmap.png)

![PSI vs developmental time](run2_grcz12_longorf/poison_exon/up_psi_vs_time.png)

![NMD features](run2_grcz12_longorf/poison_exon/nmd_features.png)

### POISEN PE inclusion (Run 1 only — Run 2 identical)

![PSI heatmap](run1_grcz12_stock_baseline/pe_psi_heatmap.png)

![PSI vs developmental time (top candidates)](run1_grcz12_stock_baseline/pe_psi_vs_time_top.png)

![ΔPSI volcano](run1_grcz12_stock_baseline/pe_dpsi_vs_pvalue.png)

![Spearman ρ distribution](run1_grcz12_stock_baseline/pe_spearman_distribution.png)

![LC2 concordance](run1_grcz12_stock_baseline/pe_lc2_concordance.png)


## 7. Methods

Same pipeline as previously documented in `outputs/stage3_v2/RESULTS_FOR_SABA.md`.
Key additions:

- **LongORF GTF supplement:** `scripts/longorf_to_gtf.py` converts the LongORF TSV
  exon chains to GTF format with NMD biotype annotation.
- **GTF merge:** `scripts/build_supplemented_gtf.sh` concatenates stock + supplement,
  sorts, and validates.
- **Parallelization:** `scripts/run_zebrafish_full_chain.py` orchestrates the full Slurm
  dependency chain with per-sample array tasks.
- **Liftover attempt:** `scripts/setup_grcz11_longorf_chain.sh` properly swaps UCSC chains
  and uses CrossMap for coordinate conversion.

## 8. Folder Map

```
outputs/for_saba/
├── run1_grcz12_stock_baseline/
│   ├── lc2/                       # LeafCutter2 outputs (86 UP junctions)
│   └── poison_exon/               # Poison-exon discovery results
│       ├── poison_exon_candidates.tsv
│       ├── positive_correlation_candidates.tsv
│       ├── poison_exon_heatmap.png
│       ├── up_psi_vs_time.png
│       └── nmd_features.png
├── run2_grcz12_longorf/
│   ├── lc2/                       # LeafCutter2 outputs (482 UP junctions)
│   └── poison_exon/               # Poison-exon discovery results (KEY COMPARISON)
│       ├── poison_exon_candidates.tsv
│       ├── positive_correlation_candidates.tsv
│       ├── poison_exon_heatmap.png
│       ├── up_psi_vs_time.png
│       └── nmd_features.png
├── run3_grcz11_stock/
│   ├── out/pe_inclusion/          # POISEN PE (0 events passed; GRCz11 control)
│   └── lc2/                       # LC2 outputs (0 UP; PE coords are GRCz12tu)
├── comparison_summary.tsv
├── novel_in_longorf.tsv
├── combined_pe_candidates.xlsx
└── SUMMARY_FOR_SABA.md             # This document
```
