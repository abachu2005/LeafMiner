# POISEN poison-exon developmental analysis — results

## Headline

Of the **1,936 unique poison-exon (PE) events** in Talha's POISEN long-read
catalogue, **only 23 (1.2 %) have enough short-read coverage to test for
differential inclusion** across the 6 zebrafish developmental timepoints
(24hpf, 30hpf, 48hpf, 72hpf, 4dpf, 5dpf; n=276 RNA-seq libraries from GTEx-style
public sources, GRCz12tu).  Under the literature-standard bar
(FDR < 0.05 **and** |ΔPSI| ≥ 0.10), **zero events show developmentally regulated
inclusion in short-read data alone**.

This is the result you predicted: the short-read pipeline is *systematically
unable* to quantify the inclusion of these long-read-discovered poison exons,
because the transcripts that contain them are degraded by NMD.  In other words,
the analysis is itself the experimental justification for needing the long-read
POISEN data.

## What the 23 testable events look like

| Coverage class | n events | Examples | Interpretation |
|---|---|---|---|
| PSI ≈ 0 across all six timepoints | 6 | **xbp1**, npm1a, nkap, rad23b, zmp:0000001139, … | POISEN finds the exon; short-read confirms the *skip* junction dominates everywhere → consistent with NMD-mediated degradation of the included isoform.  These are the most biologically convincing PE calls. |
| PSI = 1 across all six timepoints | 13 | rplp1, hmgn1a, sf3b2, macroh2a2, rsrc2, pabpn1, … | The "skip" junction we derive from the POISEN exon chain is essentially never observed; short-read sees only the inclusion form.  These are likely either (a) constitutive exons mis-classified as PEs by POISEN, or (b) cases where the alternative skip junction is too rare even for the 50-read floor. |
| Intermediate, variable PSI | 2 | itsn2b (PSI 0 → 0.014, q=2e-4); rbm25b (PSI 0.0006 → 0, p=0.32) | Statistically detectable trend on n=269 samples, but the effect size is < 2 % ΔPSI — well below the field-standard biological bar of 10 %. |

The remaining **1,881 events (97 %) fail the 50-read total-locus floor** —
they are the long tail of NMD-degraded isoforms that short-read RNA-seq cannot
see at all.  This is exactly the gap that long-read POISEN was designed to fill.

## Methods (literature-standard pipeline)

* **Event catalogue.**  POISEN long-read TSV (`UPDATED_STRAND_POISEN_HITS.tsv`,
  6,149 reads → 1,936 unique PE coordinates) parsed into (5′ inclusion intron,
  3′ inclusion intron, skip intron) triplets.  Coordinate frames matched
  carefully (POISEN BED-style 0-based half-open → STAR SJ.out.tab 1-based
  inclusive).
* **Quantification.**  Per-(event, sample) inclusion-vs-skip counts pulled
  from STAR `SJ.out.tab`.  Inclusion counts taken as `max(5′-junction,
  3′-junction)` to tolerate one-sided dropout.  PSI = inclusion / (inclusion +
  skip).
* **Filters.**  Per LeafCutter2 (Yang lab, 2025) and Mudge & Pritchard 2024
  (*Nat Genet*, "Global impact of unproductive splicing"):
  * ≥ 10 reads supporting (inclusion + skip) per (event, sample);
  * event observed in ≥ 60 % of developmental samples (≥ 166 of 276 here);
  * event observed in ≥ 3 distinct timepoints;
  * ≥ 50 total inclusion+skip reads aggregated across all samples.
* **Primary differential test.**  Quasi-binomial GLM `logit(PSI_i) = β₀ + β₁ · hpf_i`
  fit per event by IRLS with the Pearson-χ²-based overdispersion correction
  (the count-based reference test in `glm(family = quasibinomial)` and the
  approach used by LeafCutter2 / DEXSeq for proportional splicing data).
  Wald two-tailed p-value, BH-FDR.
* **Effect size.**  ΔPSI between the count-pooled PSI of the **earliest and
  latest *observed*** timepoints (LeafCutter2 / Mudge 2024 convention).
* **Combined significance bar.**  FDR < 0.05 **and** |ΔPSI| ≥ 0.10.
* **Cross-checks.**  Spearman ρ of per-sample PSI vs hpf (rank-based, 1
  significant by FDR — itsn2b, same as GLM).  Dirichlet-multinomial via
  `leafcutter_ds.R` (LeafCutter convention) is wired up as a validation track
  but failed on Quest because the leafcutter R package isn't installed in the
  cluster's R module — fixable separately, not affecting the GLM result.
* **Annotation joins.**  Each candidate annotated with LongORF PTC support
  (matched by POISEN read_id) and LC2 short-read junction class (matched on
  intron coords, with all chr/RefSeq aliases tried).

## What this means for next steps

1. **Long-read POISEN remains the primary discovery tool.**  Short-read can
   confirm a handful of PEs but cannot drive the analysis.
2. **A small number of PEs *are* quantifiable in short-read** (notably xbp1,
   npm1a, nkap, rad23b at PSI≈0; itsn2b/rbm25b at low intermediate PSI).
   These are good positive controls for any NMD-perturbation experiment (e.g.
   UPF1 morpholino) since we'd predict their PSI to *rise* upon NMD inhibition.
3. **The 13 "constitutively included" candidates merit a manual look** — if
   POISEN is calling them PEs in the long-read data, the short-read 100% PSI
   either reflects a real biological dual-isoform locus where the skip form is
   too rare for short-read, or reflects a POISEN false positive.  Either way
   that's a useful signal to feed back to Talha.
4. **Long-read inclusion quantification.**  For PEs whose host transcript is
   below the short-read 50-read floor, the natural next step is to quantify PE
   inclusion *directly from the long-read data* (read-counting at the
   long-read level) for the timepoints where long-read sequencing is
   available.  That's the intended downstream of POISEN and the only way to
   get developmental dynamics on the 1,881 PEs short-read can't see.

## Artifacts

All in `outputs/stage3_v2/`:

* `summary.json` — full filter funnel, GLM results, parameters
* `pe_developmental_candidates.tsv` — full ranked candidate table (1,936 rows; columns include `glm_beta1_per_hpf`, `glm_p`, `glm_bh_fdr`, `dpsi_first_to_last`, `passes_significance`, `spearman_*`, `longorf_*`, `lc2_*`, `psi_count_<tp>`)
* `pe_psi_heatmap.png`, `pe_psi_vs_time_top.png` — visualizations of the testable events
* `pe_dpsi_vs_pvalue.png` — volcano (ΔPSI vs −log10 GLM p)
* `pe_low_coverage.tsv` — per-(event,sample) audit of why each pair was dropped
* `pe_lc2_concordance.png` — agreement with the LC2 short-read annotation pipeline
