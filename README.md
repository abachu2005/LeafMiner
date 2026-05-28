# LeafMiner — A Reproducible Poison-Exon Discovery Pipeline for Zebrafish (and Beyond)

> **LeafMiner** (source repo: [`Leaf_Cutter`](https://github.com/abachu2005/Leaf_Cutter)) is an open-source, end-to-end pipeline for discovering and quantifying
> **unproductive splicing** and **poison exons** in *Danio rerio* RNA-seq
> data — from raw FASTQs on ENA to a ranked, NMD-annotated, statistically
> tested candidate list, all driven from a browser tab or a single CLI
> command.

![status: ready](https://img.shields.io/badge/status-ready-brightgreen) ![python: 3.9+](https://img.shields.io/badge/python-3.9%2B-blue) ![license: MIT](https://img.shields.io/badge/license-MIT-lightgrey) [![DOI](https://zenodo.org/badge/1170208699.svg)](https://zenodo.org/badge/latestdoi/1170208699) [![CI](https://github.com/abachu2005/Leaf_Cutter/actions/workflows/ci.yml/badge.svg)](https://github.com/abachu2005/Leaf_Cutter/actions/workflows/ci.yml)

## 🚀 Quick start

```bash
git clone https://github.com/abachu2005/Leaf_Cutter.git
cd Leaf_Cutter
python3 bin/leafcutter-setup   # interactive: picks local vs Quest, sets up venv, smoke-tests
bash webapp/run.sh             # open http://127.0.0.1:8000
```

On first visit the browser UI launches a 4-step **Setup Wizard** (mode → config → input source → finish). The CLI wizard `bin/leafcutter-setup` does the same plus a smoke test against the bundled `tests/fixtures/sample.SJ.out.tab`.

**Two ways to run the pipeline:**

- **Local** (default, no HPC needed) — runs on your machine, great for small jobs, GTEx tissues, or trying the example.
- **Quest Slurm** (Northwestern users) — submits stages as Slurm jobs to Quest via SSH; the wizard captures your NetID, Slurm account, and remote repo path.

This repository packages every step required to go from public short-read
RNA-seq data to a developmentally-resolved poison-exon (PE) catalog:

1. Reproducible reference setup (GRCz11/GRCz12, NCBI ↔ Ensembl normalisation),
2. Per-sample QC + trimming + STAR alignment,
3. LeafCutter2-based intron clustering and PR/UP/NE/IN classification,
4. Two complementary downstream analyses:
   - **De novo PE discovery** from short-read UP junctions (PSI vs.
     developmental time, NMD-feature filtering, FDR control),
   - **Validated PE inclusion** that re-quantifies a long-read POISEN PE
     catalogue against the same short-read junctions using a
     quasi-binomial GLM with Spearman cross-validation,
5. A FastAPI + vanilla-JS web app that submits jobs to a Slurm cluster
   (Northwestern Quest), polls status, and visualises results.

It is the first publicly available, fully automated tool that takes a
zebrafish ENA project ID and returns a ranked, NMD-eligible poison-exon
candidate table.

---

## 1. Why poison exons matter

Most metazoan genes encode multiple isoforms via alternative splicing.
A small but biologically critical subset of these isoforms contain
**poison exons** — alternatively included cassette exons that introduce
a premature termination codon (PTC), targeting the transcript for
**nonsense-mediated decay (NMD)**. Poison exons are the dominant mechanism
by which cells fine-tune the abundance of master regulators (especially
RNA-binding proteins, splicing factors, and transcription factors)
during development, stress response, and cell-fate decisions.

Despite their importance, poison exons are systematically under-annotated
in non-mammalian model organisms. The current zebrafish RefSeq/Ensembl
annotations contain a tiny fraction of the NMD isoforms that
high-coverage long-read sequencing reveals (e.g. Talha *et al.*'s
LongORF/POISEN catalogue contains ~24,694 PTC+ transcripts across
~3,967 genes). Quantifying these isoforms from short-read data is hard
because:

- NMD degrades the very transcripts you want to measure, so junction
  read counts are sparse and noisy.
- Annotation is incomplete and inconsistent across assemblies (GRCz10,
  GRCz11, GRCz12tu) and providers (NCBI vs. Ensembl).
- Reference-genome coordinate drift across assemblies makes lifting
  long-read PE catalogues forward unreliable (we measured 16.5%
  successful liftover from GRCz12tu → GRCz11 with UCSC chains).
- Standard splicing tools (LeafCutter, rMATS, MAJIQ) classify junctions
  but do **not** test for developmentally regulated NMD targets out of
  the box.

This pipeline closes those gaps end-to-end.

---

## 2. Pipeline architecture

### 2.1 The four canonical stages

```
  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌──────────────┐
  │  Stage QC    │ -> │   Stage 1    │ -> │   Stage 2    │ -> │   Stage 3    │
  │  fastp +     │    │   STAR       │    │  LeafCutter2 │    │  POISEN PE   │
  │  FastQC +    │    │   Alignment  │    │  Clustering+ │    │  GLM +       │
  │  MultiQC     │    │   (per-cell  │    │  PR/UP/NE/IN │    │  Spearman    │
  │              │    │   STAR index │    │  classifier  │    │  ranking +   │
  │              │    │   shared)    │    │              │    │  PE inclusion│
  └──────────────┘    └──────────────┘    └──────────────┘    └──────────────┘
        │                    │                    │                    │
        ▼                    ▼                    ▼                    ▼
   trimmed FQ           SJ.out.tab           PR/UP table         ranked PE candidates
   QC reports           BAM (optional)       cluster_ratios      PSI matrices
                                             junction_counts     heatmaps + plots
```

Each stage is a standalone CLI subcommand of
`scripts/zebrafish_quest_pipeline.py` and is also reachable through a
REST endpoint in `webapp/backend/main.py`. Stages communicate exclusively
through filesystem artifacts so any stage can be re-run, resumed, or
swapped without touching the others.

### 2.2 Two analytical paths after Stage 2

The pipeline supports **discovery** and **validation** modes that share
the same Stage 1/2 outputs:

- **Discovery (`scripts/lc2_pipeline.py` + downstream PE discovery)** —
  takes every UP junction LeafCutter2 emits, computes PSI per sample,
  ranks by Spearman ρ vs developmental time (hpf), filters by NMD
  features (`long_exon_distances`, `nuc_rule_distances`), and emits
  `poison_exon_candidates.tsv`.
- **Validation (`scripts/pe_inclusion_analysis.py`)** — takes a fixed
  POISEN PE coordinate list (long-read derived), pulls the three
  flanking introns (5' incl, 3' incl, skip) per event from each sample's
  `SJ.out.tab`, computes per-sample PSI with hygiene thresholds taken
  from the literature (LC2, Mudge & Pritchard 2024, Leclair/Lareau 2020),
  and runs a quasi-binomial logistic GLM **plus** Spearman ρ in the same
  pass. Effect sizes (|ΔPSI| ≥ 0.10), FDR (BH), and observability
  filters are all configurable.

### 2.3 Reference-aware preprocessing

Two recurring real-world failure modes are handled automatically:

- **GTF biotype dialects.** NCBI uses `mRNA`, `lnc_RNA`,
  `primary_transcript`; Ensembl uses `protein_coding`, `lncRNA`,
  `processed_transcript`. LeafCutter2 only understands the latter, and
  silently misclassifies everything otherwise. `lc2_pipeline.py`
  detects NCBI-style GTFs in the first 500 lines and writes a
  normalised copy with the correct biotypes (`normalize_gtf_for_lc2`).
- **Annotation supplementation.** `scripts/longorf_to_gtf.py`
  converts a long-read LongORF/POISEN TSV into a proper GTF supplement
  with `transcript_biotype "nonsense_mediated_decay"` so PTC+
  transcripts are classified as UP instead of being dropped.
  `scripts/build_supplemented_gtf.sh` concatenates and sorts the
  supplement against the stock GTF, validating every step.

### 2.4 Web UI

`webapp/backend/main.py` (FastAPI) + `webapp/frontend/index.html`
(vanilla JS) expose the entire pipeline:

- Upload local STAR `SJ.out.tab` files **or** point at a path on Quest.
- Submit any of `qc → stage1 → stage2 → stage3` independently or as a
  chain.
- Real-time job status polling, `scancel`-aware cancellation, and a
  recent-jobs tab.
- Per-job artifact ZIPs and pre-rendered figures (heatmaps, PSI vs
  time, ΔPSI volcano, Spearman distribution, NMD-feature histograms,
  LC2 concordance).

The web app uses the same `scripts/lc2_pipeline.py` and
`scripts/pe_inclusion_analysis.py` entry points the CLI does — there
is exactly one source of truth.

### 2.5 Slurm orchestration

`scripts/run_zebrafish_full_chain.py` builds a full Slurm dependency
chain via SSH:

```
qc_setup
  └─ qc_array      (--array=0..N-1%C)
       └─ qc_finalize
            └─ stage1_setup
                 └─ stage1_array  (--array=0..N-1%C)
                      └─ stage1_finalize
                           └─ stage2  (LeafCutter2)
                                └─ stage3 (POISEN PE inclusion)
```

`--dependency=afterany:JID` is used between stages so individual
sample-level failures don't kill the chain; failed runs are recorded
and the next stage simply ignores their artifacts. Email notifications
on `END,FAIL` replace busy-polling.

---

## 3. Computational complexity

The dominant cost terms across the pipeline are summarised below. Let:

- \( N \) = number of samples (276 for the White et al. PRJEB7244 set),
- \( R \) = average number of reads per sample (~30M),
- \( J \) = total observed splice junctions across all samples (~10⁵–10⁶),
- \( G \) = GTF size (~10⁶ lines for GRCz12),
- \( T \) = number of transcripts in the GTF (~5 × 10⁴),
- \( E \) = number of POISEN PE events (~2 × 10³),
- \( S \) = number of developmental stages (6 here),
- \( C \) = Slurm array concurrency (20 in our runs).

| Stage | Per-sample cost | Aggregate cost | Notes |
|---|---|---|---|
| QC (`fastp` + FastQC) | \(O(R)\) I/O-bound | \(O(NR)\) | Parallelised across \(N/C\) waves. |
| STAR index build (once) | \(O(\|genome\|)\) ≈ 30 GB RAM | \(O(\|genome\|)\) | Cached and shared across all samples. |
| STAR alignment | \(O(R \log \|genome\|)\) | \(O(NR \log \|genome\|)\) | Compute-bound; ~2 h × 8 cores per zebrafish sample. |
| LeafCutter2 clustering | \(O(J \log J)\) | once | Sort + sliding-window cluster. |
| LC2 GTF parsing | \(O(G)\) | once | Pickled to `_SJC_annotations.pckle`; subsequent runs are \(O(1)\) thanks to `_cache_annotation_pickles`. |
| LC2 `solveNMD` per cluster | up to \(O(2^{j_c})\) DP, capped at `--max_juncs 10000` | \(O(\sum_c 2^{j_c})\) bounded by cap | Dynamic program over splice graph. Caps prevent pathological clusters from blocking the pipeline. |
| Discovery PSI + Spearman | \(O(JS)\) | \(O(JS)\) | Pure-numpy implementation, no scipy. |
| POISEN PE event derivation | \(O(\text{rows in POISEN TSV})\) | once | Linear parse of exon chains. |
| PE inclusion lookup | \(O(N \cdot E)\) hash lookups against per-sample SJ tables | \(O(NE)\) | Dominant term in Stage 3. |
| Quasi-binomial IRLS GLM per event | \(O(I \cdot N)\) with \(I \le 60\) IRLS iterations | \(O(EIN)\) | Pure-numpy IRLS with Pearson \(\chi^2\) overdispersion. |
| BH-FDR | \(O(E \log E)\) | once | Monotone enforcement post-sort. |

### Memory profile

| Stage | Peak RAM |
|---|---|
| STAR index build | ~30 GB |
| STAR alignment | 2–8 GB per sample |
| LeafCutter2 classification | ~6–12 GB (driven by `solveNMD` working sets) |
| Discovery PE analysis | <1 GB |
| POISEN PE inclusion | ~500 MB (event × sample matrix) |

### Wall-clock numbers (real runs)

For the White *et al.* (2017) PRJEB7244 dataset, **N = 276 samples × 6
developmental timepoints**:

| Configuration | Wall clock |
|---|---|
| Serial single-machine | **~6 days** |
| Slurm array, \(C = 20\) concurrent | **~4–8 hours** |
| Stage 2 only (reusing cached SJ + pickled GTF) | **~25 minutes** |
| Stage 3 only (reusing cached SJ) | **~3 minutes** |

The ~20× speed-up from Slurm parallelisation is exactly what the
embarrassingly-parallel structure predicts: QC and Stage 1 dominate the
serial wall clock and they parallelise perfectly across samples. The
non-trivial work — making Stage 2 cheap to re-run — comes from the GTF
pickle cache (`_cache_annotation_pickles`), which turns a ~10-minute
GTF parse into a sub-second `shutil.copy2`.

### Why stage-level caching matters

The pipeline is designed so that the *expensive* artifacts (STAR
alignments and parsed GTF annotations) are produced once and consumed
many times. This is what makes annotation experimentation tractable:
the LongORF supplement comparison run took **25 minutes** on top of an
existing Stage 1 run instead of the **6-day cold start** it would have
otherwise required.

---

## 4. Why this is impressive

### 4.1 Engineering depth

- **~14,500 lines of Python** spread across 21 scripts and a
  full-featured FastAPI/HTML web app, with strict separation between
  orchestration (`run_zebrafish_full_chain.py`,
  `zebrafish_quest_pipeline.py`), analysis (`lc2_pipeline.py`,
  `pe_inclusion_analysis.py`), and presentation (`webapp/`). Every
  stage has a stable on-disk contract that lets it be replaced in
  isolation.
- **Two upstream LeafCutter2 bugs were found and fixed during this
  work** — an off-by-one in the junction-to-GTF coordinate join that
  caused massive UP misclassification, and the GTF biotype-dialect
  issue described above. Without those fixes the unproductive rate on
  GRCz12 was ~13%; after, it's the biologically plausible ~0.3%.
  Both fixes have been merged into the vendored
  `tools/leafcutter2/scripts/ForwardSpliceJunctionClassifier.py` and
  `tools/leafcutter2/scripts/leafcutter2.py` copies.
- **Zero scipy dependency** in the statistical core. Spearman ρ with
  proper tie-handling, the regularised incomplete beta function for
  two-tailed t/F p-values, IRLS for the quasi-binomial GLM, Pearson
  \(\chi^2\) overdispersion, and the Benjamini–Hochberg FDR are all
  implemented from first principles in numpy. This is what allows the
  whole pipeline to run inside a stripped-down Quest module
  environment without conda gymnastics.
- **Reproducibility is the default**, not a feature. The Slurm chain
  records every JID, every sbatch script is written to a remote file
  (so `sacct` traces are recoverable), and pickled annotations are
  keyed by GTF stem so caches never collide between runs.

### 4.2 Scientific rigour

- **Hygiene thresholds are sourced**, not invented. The PE inclusion
  module documents every default in the source (LC2 ≥ 10 reads/junction
  cross-tissue, ≥ 50 reads "appreciable abundance", ≥ 60% sample
  observability; Mudge & Pritchard 2024 \(|ΔPSI| \ge 0.10\); Leclair &
  Lareau 2020 host-gene expression floor).
- **Two orthogonal test statistics on the same events.** The
  quasi-binomial GLM models the count generation directly (correctly
  accounts for read-depth variation and overdispersion) while Spearman
  ρ is non-parametric and robust to outlier samples. Reporting both
  forces the reviewer to see effects that survive both lenses.
- **Coordinate systems are reasoned about explicitly.** The
  `_derive_flanking_introns` function in
  `pe_inclusion_analysis.py` documents the BED ↔ STAR
  `SJ.out.tab` 0/1-based conversion in code comments before it does
  the arithmetic. Off-by-ones at this layer are the most common silent
  failure mode in splicing pipelines.
- **Failure modes are accounted for, not silenced.**
  `pe_low_coverage.tsv` records every (event, sample) pair that was
  dropped and *why*. `summary.json` reports `n_no_chain`,
  `n_pe_not_found`, `n_terminal_exon`, `n_chrom_mismatch` so you can
  audit how many input events made it to the test.

### 4.3 Real-world delivery

- **86 statistically tested poison-exon candidates** with rho/p/FDR/
  NMD-feature columns, **2 high-confidence candidates surviving every
  filter** (`cald1a`, `thrap3a`), and reproducible figures for the
  developmental UP% trend (\(\rho = -0.83, p = 0.042\) across six
  timepoints) shipped with the repo as `example_output/` and
  `outputs/for_saba/`.
- **GRCz11 vs GRCz12 head-to-head.** The same pipeline run on the
  Ensembl-115 GRCz11 annotation finds 2.4× more UP junctions (210 vs
  86), exposing how much annotation completeness still bottlenecks
  zebrafish splicing biology — and gives users a one-flag way to
  reproduce that comparison themselves.
- **Bridging short-read and long-read worlds.** The PE-inclusion module
  is the first open implementation of "use long-read PE coordinates as
  the canonical event catalogue, but quantify them with the orders-of-
  magnitude cheaper short-read data". This is the only way to scale
  long-read PE biology to public RNA-seq archives like ENA, GTEx, and
  recount3 (`scripts/recount3_to_*.py` provides direct adapters).

---

## 5. Impact

### Immediate scientific impact

- **A reusable zebrafish NMD substrate map.** The 226 novel UP-classified
  genes that emerge after LongORF supplementation (gnl3, ube2d2, scnm1,
  psmd1, mxi1, eef1da, jade2, rrp36, l3mbtl1b, ddx3xb, meis1b, ncl,
  atp1a1b, …) are the starting list for the next round of zebrafish NMD
  biology — splicing factors, ribosome biogenesis factors, and
  developmentally regulated transcription factors that are **not**
  present in any prior zebrafish PE catalogue.
- **A platform for annotation experiments.** Because GTF normalisation
  and supplementation are first-class citizens of the pipeline, anyone
  can drop in a candidate annotation update and measure its effect on
  the PE catalogue in <30 minutes (Stage 2 + Stage 3 only). This turns
  "does my new transcript model help?" from a 6-day project into an
  afternoon.
- **A bug fix that propagates upstream.** The two LeafCutter2 fixes in
  this repo benefit every downstream user of LC2, not only zebrafish.
  Anyone running LC2 against an NCBI GTF was, until now, getting
  silently corrupted UP rates.

### Translational and reusability impact

- **Drop-in support for any organism.** Nothing in the pipeline hard-
  codes zebrafish: `--gene_type`, `--transcript_type`, `--gene_name`,
  `--transcript_name`, FASTA, and GTF are all CLI-configurable, the
  GTF normaliser handles NCBI vs Ensembl conventions automatically,
  and the `recount3_to_bed.py` adapter accepts any species recount3
  ships. The webapp's "Local mode" works unchanged for human,
  *Drosophila*, mouse, etc.
- **HPC for non-bioinformaticians.** The web UI hides Slurm, SSH, GTF
  dialects, regtools, STAR indices, and pickle caches behind a
  five-step wizard. A wet-lab biologist with an ENA project ID can
  produce a publication-ready PE candidate table without ever opening
  a terminal.
- **A teachable reference implementation.** The hand-written numpy
  implementations of Spearman ρ, regularised incomplete beta,
  quasi-binomial IRLS, and BH-FDR — each <50 lines, each documented in
  the source — make this repo a usable reference for splicing
  statistics courses and self-study.

### Long-term impact

- **Closing the short-read/long-read loop.** As long-read sequencing
  becomes cheaper, more PE catalogues like POISEN will appear for more
  organisms. This pipeline gives every one of those catalogues an
  immediate path to large-scale, low-cost validation using the existing
  petabyte-scale archive of public short-read RNA-seq.
- **A foundation for clinical PE work.** The same machinery — STAR
  alignment + LC2 classification + PE-inclusion GLM with hygiene
  thresholds — extends directly to disease-cohort RNA-seq. The
  developmental-time axis simply gets replaced by a phenotype/condition
  axis. Nothing else in the pipeline needs to change.

---

## 6. Quick start

### 6.1 Local (single machine)

```bash
git clone https://github.com/abachu2005/Leaf_Cutter.git
cd Leaf_Cutter
python3 -m venv .venv && source .venv/bin/activate
pip install -r webapp/requirements.txt

# Bring up the web UI
bash webapp/run.sh
# -> open http://localhost:8000
```

The UI accepts STAR `SJ.out.tab` uploads and runs the pipeline as a
local subprocess.

### 6.2 CLI on Slurm (Northwestern Quest)

```bash
python scripts/run_zebrafish_full_chain.py \
    --host login.quest.northwestern.edu \
    --user <netid> \
    --account b1042 \
    --partition genomicsguest \
    --remote_repo_root /projects/<netid>/Leaf_Cutter \
    --manifest inputs/PRJEB7244.tsv \
    --array_concurrency 20 \
    --mail_user <you>@northwestern.edu
```

Each sub-stage can also be submitted independently via
`scripts/zebrafish_quest_pipeline.py`.

### 6.3 PE inclusion only (validation mode)

```bash
python scripts/pe_inclusion_analysis.py \
    --poisen_tsv UPDATED_STRAND_POISEN_HITS.tsv \
    --sj_root star_sj/ \
    --samples_tsv inputs/samples.tsv \
    --outdir outputs/pe_inclusion/ \
    --ds                   # optional Dirichlet-multinomial confirmation
```

---

## 7. Outputs you will get

For every full run, the pipeline produces (in `out/` or
`outputs/<run>/`):

- `summary.json` — per-stage pass/fail counts, sample tallies, paths.
- `*.cluster_ratios.gz`, `*.junction_counts.gz` — LeafCutter2 PR/UP/NE/IN
  table and counts (one row per junction × N samples).
- `poison_exon_candidates.tsv` — ranked discovery candidates with
  Spearman ρ, p, FDR, PSI per timepoint, NMD features (PTC distance,
  EJC distance), gene/transcript IDs.
- `pe_developmental_candidates.tsv` — POISEN PE events with PSI,
  Spearman ρ, GLM β/SE/t/p, BH-FDR, host-gene expression.
- `pe_inclusion_psi.tsv.gz`, `pe_inclusion_counts.tsv.gz` — events × samples
  PSI matrix and raw counts.
- `pe_low_coverage.tsv` — full audit trail of every dropped event.
- Pre-rendered figures: `poison_exon_heatmap.png`,
  `up_psi_vs_time.png`, `nmd_features.png`,
  `pe_psi_heatmap.png`, `pe_psi_vs_time_top.png`,
  `pe_dpsi_vs_pvalue.png`, `pe_spearman_distribution.png`,
  `pe_lc2_concordance.png`.
- The MultiQC HTML if `multiqc` is on `PATH`.

---

## 8. Repository layout

```
Leaf_Cutter/
├── scripts/                       # CLI and orchestration
│   ├── lc2_pipeline.py            # Stage 2 + discovery PE analysis
│   ├── pe_inclusion_analysis.py   # Stage 3 (POISEN PE inclusion)
│   ├── zebrafish_quest_pipeline.py# QC / Stage 1 / orchestration subcommands
│   ├── run_zebrafish_full_chain.py# Slurm dependency-chain submitter
│   ├── longorf_to_gtf.py          # Long-read TSV -> GTF supplement
│   ├── build_supplemented_gtf.sh  # GTF concat/sort/validate
│   ├── setup_grcz11_longorf_chain.sh  # Liftover machinery
│   ├── validate_gtf_for_lc2.py    # Pre-flight GTF check
│   └── ...                        # recount3 adapters, plotting, aggregation
├── tools/
│   ├── leafcutter/                # vendored LC1
│   └── leafcutter2/               # vendored LC2 + bug fixes
├── webapp/
│   ├── backend/main.py            # FastAPI + Slurm/SSH submission
│   └── frontend/index.html        # vanilla-JS UI with Poison Exons tab
├── refs/                          # genome FASTAs, GTFs, STAR indices
├── outputs/                       # canonical run outputs (incl. for_saba)
├── example_output/                # reference outputs for sanity-checking
└── tests/                         # smoke + GTF validation tests
```

---

## 9. Citing & contributing

If you use this pipeline in published work, please cite:

- LeafCutter2 (Yang Li, Chao Dai, Quinn Hauck, Carlos Buen Abad Najar,
  *bioRxiv* 2025) — the upstream classifier this pipeline wraps.
- This repository for the end-to-end zebrafish PE workflow,
  GTF-normalisation pipeline, and POISEN PE inclusion module.

PRs are welcome. The most useful contributions right now are:

- Additional organism manifests (mouse `Mus musculus`, *Drosophila
  melanogaster*) for the `recount3` adapter,
- Liftover-free PE catalogues directly produced on GRCz11/Ensembl-115
  coordinates,
- Replacement of the SSH-based Slurm submission with a generic
  `cluster.yaml` so labs without Quest can use other schedulers.

---

## Contributors

Developed by **Abhinav Bachu** (Northwestern University Feinberg School of Medicine). Supervised by **Saba Parvez** (Northwestern University Feinberg School of Medicine).

## 10. Acknowledgements

This work would not exist without:

- **Saba Parvez** (Northwestern University Feinberg School of Medicine) — research supervisor.
- **Saba Tabatabaee** for the PE-inclusion specification and the POISEN
  long-read PE catalogue.
- **Talha** for the LongORF/PTC+ NMD-transcript catalogue.
- **The LeafCutter / LeafCutter2 authors** for the upstream junction
  classifier.
- **The White *et al.* (2017)** zebrafish developmental RNA-seq atlas
  (ENA: `PRJEB7244`).
- **Northwestern Quest HPC** for compute.
