Hi Saba,

Weekly update on the zebrafish developmental splicing project.

**Summary:** I attempted to expand the UP junction count by supplementing the GRCz12 GTF with Talha's LongORF PTC+ NMD transcript models (24,694 transcripts from 3,967 genes). In the process I found and fixed a critical pipeline bug, but the supplement itself had no effect on the results. The original Run 1 numbers from my previous email remain the correct and final results for GRCz12. I also parallelized the full pipeline using Slurm job arrays, cutting wall-clock time from ~6 days to ~4-8 hours per full run.

**1. Bug found and fixed**

While validating the LongORF supplemented run, I noticed the unproductive read percentage jumped from 0.3% to 13% — which was clearly wrong. After investigating, I traced it to a GTF attribute tag mismatch in the pipeline's submission script: the NCBI GRCz12 GTF uses `transcript_biotype` as its attribute name, but the script was passing `transcript_type` to LeafCutter2's classifier. This caused LC2 to fail to recognize any protein-coding transcripts, misclassifying ~40,000 junctions as unproductive.

The fix was straightforward (correcting the default tag parameters in three scripts), and I've re-run the analysis with the correct settings. The bug only affected Run 2 (LongORF supplement); Run 1 (submitted through the webapp, which had the correct tags) and Run 3 (GRCz11 Ensembl GTF, which uses `transcript_type`) were unaffected.

**2. LongORF supplement: no effect**

After fixing the bug, the corrected Run 2 produces identical results to Run 1:

|                        | Run 1 (GRCz12 stock) | Run 2 (GRCz12 + LongORF) | Run 3 (GRCz11 Ensembl) |
|------------------------|----------------------|--------------------------|------------------------|
| UP junctions           | 86                   | 86                       | 210                    |
| UP read %              | 0.31%                | 0.31%                    | 1.17%                  |
| PR junctions           | 3,052                | 3,052                    | 2,234                  |

Zero of the 2,608 new junction coordinates from the LongORF NMD transcripts overlap with any observed splice junction in the RNA-seq data. The LongORF transcripts define splice sites at different genomic positions than what's actually observed in the White et al. reads — likely because NMD degrades these isoforms before junction-spanning reads accumulate to detectable levels.

**3. What this means**

The limited number of UP junctions (86) in GRCz12 is not due to an annotation gap that can be fixed by adding more transcript models. The bottleneck is that NMD-targeted isoforms are degraded, so their unique splice junctions simply aren't present in the sequencing data. LeafCutter2 can only classify junctions that have reads.

That said, GRCz11 with the more mature Ensembl 115 annotation does find 2.4x more UP junctions (210 vs 86), suggesting the older assembly's annotation captures some NMD-relevant transcript structures that the newer GRCz12 RefSeq annotation doesn't. This likely comes from Ensembl's inclusion of more alternative transcript models that share splice sites with productive isoforms.

**4. Validated results (unchanged from previous email)**

The Run 1 results remain the definitive findings:
- Global UP% trend: Spearman rho = -0.83, p = 0.042 (UP splicing decreases during development)
- 86 UP junctions, 21 statistically significant, 13 NMD-eligible
- 2 high-confidence candidates: cald1a, thrap3a
- Notable hits lacking full NMD annotation: krt91 (rho = -1.0), trdn (rho = -1.0), hnrnpc, hnrnpa1b (rho = -0.89)

**5. Pipeline improvements this week**
- Fixed GTF attribute tag defaults across all CLI scripts to match NCBI convention (`gene_biotype`/`transcript_biotype`)
- Parallelized pipeline via Slurm job arrays (full run: ~4-8 hours vs ~6 days serial)

I've attached the GRCz11 results since those are new:
- **poison_exon_candidates.tsv** — full ranked list of 210 UP junctions with PSI per timepoint and Spearman stats
- **poison_exon_heatmap.png** — PSI heatmap for the top candidates across developmental stages

Some notable GRCz11 hits: krt91 shows the same dramatic PSI drop (0.83 → 0.15) as in GRCz12, and genes like nsfa, nap1l1, med11, tcea3, col12a1b show strong positive correlations (rho = +1.0) that weren't captured in GRCz12.

Happy to discuss next steps — in particular whether we should focus on deeper validation of the existing candidates or explore the GRCz11 results further.

Thanks,
Abhi
