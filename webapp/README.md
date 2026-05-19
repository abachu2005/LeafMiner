# LeafCutter2 Web App

Web interface for the LeafCutter2 splicing analysis pipeline. Supports local
execution and remote submission to Northwestern Quest via Slurm.

## Quick Start

```bash
cd /path/to/Leaf_Cutter
bash webapp/run.sh
```

Then open **http://localhost:8000** in your browser.

Or manually:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r webapp/requirements.txt
uvicorn webapp.backend.main:app --host 0.0.0.0 --port 8000 --reload
```

## Features

- Upload STAR `SJ.out.tab` junction files from browser
- Run pipeline locally or submit to Quest Slurm cluster
- Separate zebrafish Quest flow: project IDs -> run discovery -> STAR -> LeafCutter2 -> POISEN PE inclusion (Step 4)
- Long-read PE tab: ranks POISEN PE events by Spearman \u03c1(PSI vs developmental hpf)
  with FDR, LongORF/LC2 cross-reference, heatmap and per-event PSI vs time plots
- Real-time job status polling
- Cancel running jobs (local process kill or Slurm `scancel`)
- View pipeline summary and log output in-browser
- Download output artifacts as ZIP
- Recent jobs list

## API Endpoints

| Method | Path                    | Description                    |
|--------|-------------------------|--------------------------------|
| POST   | `/jobs`                 | Create and launch a job        |
| GET    | `/jobs`                 | List recent jobs               |
| GET    | `/jobs/{id}`            | Get job status and metadata    |
| GET    | `/jobs/{id}/results`    | Get summary JSON and log tail  |
| GET    | `/jobs/{id}/download`   | Download artifacts ZIP         |
| POST   | `/jobs/{id}/cancel`     | Cancel a running job           |
| POST   | `/jobs/zebrafish/qc`    | Step 1: QC + trim (Quest)      |
| POST   | `/jobs/zebrafish/stage1`| Step 2: STAR alignment (Quest) |
| POST   | `/jobs/zebrafish/stage2`| Step 3: LeafCutter2 (Quest)    |
| POST   | `/jobs/zebrafish/stage3`| Step 4: POISEN PE inclusion (Quest) |

## Execution Modes

### Local Mode

Runs the pipeline as a subprocess on the machine hosting this web app.
Suitable for testing and small datasets.

### Quest Slurm Mode

Submits the pipeline as a batch job to Northwestern Quest via SSH.
Requires:

- **SSH key access** from your machine to `login.quest.northwestern.edu`
- A valid **Slurm account** (e.g., `b1042` for GCC, or `pXXXXX` project)
- A valid **partition** (e.g., `genomicsguest`, `genomics`, `genomics-gpu`)
- The **Leaf_Cutter repository cloned on Quest** with refs available

You can either upload STAR files (they get `scp`'d to Quest) or provide a
path to a junction file list already present on Quest.

### Zebrafish Mode (Quest only, two-stage)

The UI includes a dedicated **Zebrafish (Public Projects)** source with two
separate actions:

1. **Stage 1 (Alignment)**
   - expands project IDs (`PRJNA`/`ERP`/`DRP`) to runs via ENA filereport,
   - extracts tissue labels from ENA sample metadata,
   - downloads FASTQ files and aligns with STAR,
   - persists STAR junction artifacts and stage1 metadata.
2. **Stage 2 (LeafCutter2)**
   - takes a completed Stage 1 job ID,
   - optionally selects a subset of discovered tissues in the UI,
   - reuses stored STAR junctions/samples metadata,
   - runs `scripts/lc2_pipeline.py` to produce final LC2 outputs.

Required references on Quest (under your repo path by default):

- `refs/danRer11.fa` (and `.fai`)
- `refs/Danio_rerio.GRCz11.115.gtf` (or equivalent GRCz11 annotation)

The zebrafish GTF tag defaults are Ensembl-style:

- `gene_biotype`
- `transcript_biotype`
- `gene_id`
- `transcript_id`

Stage 2 can be submitted hours or days later as long as Stage 1 artifacts are
still present on Quest.

## Job Lifecycle

```
queued -> submitted (Slurm accepted) -> running -> succeeded
                                                -> failed
                                                -> cancelled
```

## Storage

- Job metadata: `webapp/data/jobs.db` (SQLite)
- Job files: `webapp/data/jobs/<job_id>/`
- The `webapp/data/` directory is gitignored

## Pipeline Contract

The web app invokes `scripts/lc2_pipeline.py` which calls:

1. STAR SJ.out.tab -> BED conversion
2. LeafCutter clustering (`leafcutter_cluster_regtools.py`)
3. LeafCutter2 classification (`leafcutter2.py -j -r -o -A -G`)

Outputs are written to `<workdir>/out/summary.json` and include paths to:

- `*.cluster_ratios.gz`
- `*.junction_counts.gz`
- `clustering/*_long_exon_distances.txt`
- `clustering/*_nuc_rule_distances.txt`

For zebrafish jobs, the contract is intentionally the same, so the UI reuses
the existing summary, tissue chart, log, and artifact download components.

## Troubleshooting

**Slurm submission fails with "Invalid account"**
- Verify your account name with `groups` on Quest
- Check allocation status with `checkproject <account>`

**Job stays in "submitted" forever**
- Check Quest queue: `squeue -u <netid>`
- The partition may be full; try `short` for test runs

**Pipeline fails with missing refs**
- Ensure `refs/GRCh38.fa` and `refs/gencode.v46.annotation.gtf` exist
- On Quest, use absolute paths via `remote_repo_root`
