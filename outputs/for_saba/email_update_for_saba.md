Hi Saba,

Following up on my last email about the zebrafish poison-exon analysis — I have some significant results to share regarding the annotation gap I flagged.

**The problem (recap):** In the initial run, LeafCutter2 only identified 86 UP (unproductive) junctions, with 21 significant and just 2 final candidates (cald1a, thrap3a). I suspected the limited number was due to incomplete NMD transcript annotation in the new GRCz12tu assembly, since many of the significant hits (krt91, trdn, hnrnpc, hnrnpa1b) lacked proper NMD annotations in the stock GTF.

**What I did:** I ran three parallel comparisons on all 276 samples to test this:

1. **Run 1** — GRCz12tu with stock RefSeq GTF (baseline, same as before)
2. **Run 2** — GRCz12tu with Talha's LongORF PTC+ supplement added (24,694 NMD transcripts from 3,967 genes appended to the GTF)
3. **Run 3** — GRCz11 (danRer11) with Ensembl 115 GTF (cross-assembly validation)

**Key results — the LongORF supplement made a massive difference:**

|                                     | Run 1 (stock) | Run 2 (+ LongORF) | Run 3 (GRCz11) |
|-------------------------------------|---------------|---------------------|----------------|
| UP junctions classified             | 86            | **482 (5.6×)**      | 210 (2.4×)     |
| Significant (Spearman p<0.05)       | 21            | **127 (6×)**        | 53 (2.5×)      |
| Positive PSI-time correlation       | 11            | **84 (7.6×)**       | 31 (2.8×)      |
| Unique genes with UP junctions      | 72            | **294 (4.1×)**      | 134 (1.9×)     |
| NMD-eligible                        | 13            | 4*                  | 18             |

*NMD-eligibility appears to drop in Run 2 only because the newly discovered UP junctions from LongORF transcripts lack entries in LC2's auxiliary annotation files — the junctions themselves are real and significant.

This confirms the annotation gap hypothesis: the stock GRCz12tu GTF was missing most of the NMD transcript models needed for LC2 to classify junctions as unproductive. Adding the LongORF supplement fixed that.

**GRCz11 cross-check:** Interestingly, the older Ensembl 115 GRCz11 annotation also found substantially more UP junctions than stock GRCz12 (210 vs 86), consistent with GRCz11 having more mature NMD annotations. However, the gene overlap between assemblies is low — only 34 genes are shared between GRCz12 stock and GRCz11, and just **3 genes** show significant positive developmental correlation in both assemblies:

- **fam133b** (rho = +0.94 in both)
- **eif3jb** (rho = +0.94 in both)
- **dedd1** (rho = +0.94 / +0.83)

These are the highest-confidence cross-assembly validated hits.

**Notable new candidates from Run 2 (LongORF supplement):** 37 novel genes with significant positive PSI-developmental time correlation appeared only after adding the LongORF annotation. The top hits include:

- hnrnpa1b (rho = +1.0) — a known splicing regulator I flagged in my previous email as significant but lacking NMD annotation, now confirmed
- gnl3, ube2d2, eef1da, jade2, meis1b, ddx3xb, psmd1, mxi1 (all rho ≥ +0.94)

The full candidate lists are in the attached TSV files.

**What I think this means:** The original 2-candidate result (cald1a, thrap3a) was severely limited by annotation completeness, not by biology. With proper NMD annotation, we're seeing an order of magnitude more candidates. The next question is which of these 84 positive-correlation candidates warrant deeper validation — I think starting with the 3 cross-assembly confirmed hits (fam133b, eif3jb, dedd1) plus the known splicing regulators (hnrnpa1b, ddx3xb) would make the strongest case.

I also parallelized the pipeline using Slurm job arrays, which cut the full-run wall-clock time from ~6 days to ~4-8 hours.

Happy to walk through the detailed results. All output files are in `outputs/for_saba/` — the key comparison is between `run1_grcz12_stock_baseline/poison_exon/` and `run2_grcz12_longorf/poison_exon/`.

Thanks,
Abhi
