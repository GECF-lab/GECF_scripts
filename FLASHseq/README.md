# FLASH-seq (single-cell) RNA-seq pipeline

Run-to-result pipeline for **Alithea FLASH-seq** single-cell RNA-seq data
(FastQC → alignment → gene-level counts → merged matrix), plus standalone
per-sample RSeQC QC on the resulting STAR BAMs. Each run is re-runnable:
finished steps are skipped automatically (resume), and a failing sample never
blocks the others.

## Contents

| File | Purpose |
|---|---|
| `flashseq_pipeline.sh` | Core pipeline: FASTQ → FastQC → STAR → featureCounts → merged genes × samples matrix |
| `flashseq_rseqc.sh` | Per-sample RSeQC QC (gene-body coverage + read distribution) on existing BAMs |
| `gbcov_redraw.py` | Merged gene-body coverage plot: PDF/PNG/JPEG (via R) or self-contained interactive HTML (via plotly) |
| `flashseq-env.yaml` | Conda environment definition for all of the above |

`gbcov_redraw.py` must remain in the **same directory** as
`flashseq_rseqc.sh` — the QC script resolves it relative to its own location.

## Environment

Both scripts use the `flashseq` conda environment:

```bash
conda activate flashseq
# or, to rebuild it from this folder:
conda env create -f flashseq-env.yaml
```

Core tools provided by the environment: STAR, samtools, featureCounts,
fastqc, RSeQC, gffread, R, and Python (with plotly for HTML output). The
scripts hard-fail at startup, naming any missing tool, so nothing half-starts.

## Input requirements

A project directory whose sample directories each contain exactly one
paired-end set of FASTQ files, e.g.:

```
<project>/
├── avid12345/
│   ├── lung_001_S1_R1_001.fastq.gz
│   └── lung_001_S1_R2_001.fastq.gz
└── avid12346/
    └── ...
```

- Missing pairs → the sample is skipped and listed in `logs/_failed.samples`.
- The sample **label** is derived from the R1 FASTQ filename
  (`lung_001_S1_R1_001.fastq.gz` → `lung_001`) and names the BAMs, count
  tables, and matrix columns.
- The run label (AVT id) is inferred from the `--fastq-dir` path when
  possible; override with `--avt`.

---

## 1. `flashseq_pipeline.sh` — core pipeline

Per sample:

1. **FastQC** — raw-read QC (zip reports per FASTQ)
2. **STAR** — splice-aware alignment, coordinate-sorted BAM per sample
3. **QC gate** — uniquely-mapped % recorded to `logs/_gates.csv`
   (warning if < 70 %; informational, never filters)
4. **featureCounts** — gene-level counts
   (`-p -T 1 -t exon -g gene_id --fracOverlap 0.25`)

After all samples finish, an R merge step combines every per-sample count
table into one **genes × samples matrix** and writes a per-cell
genes-detected table (thresholds ≥ 1 and ≥ 5).

### Usage

```bash
bash flashseq_pipeline.sh --species <human|mouse|name> --fastq-dir <DIR> [options]
```

### Options

| Flag | Meaning | Default |
|---|---|---|
| `--species` | `human` / `mouse` use the built-in references; any other value is *custom* — then `--gtf` is required, plus a prebuilt `--star-index` or `--ref-fasta` (the STAR index is built on demand from the FASTA + GTF) | — |
| `--fastq-dir` | Root containing the `avid*/` sample directories (required) | — |
| `--avt` | Run label; inferred from `--fastq-dir` if omitted | inferred |
| `--samples a,b` | Process only these sample directories | all found |
| `--exclude a,b` | Skip these sample directories entirely | none |
| `--max-jobs N` | Samples in parallel (a concurrency cap, not a sample limit) | 4 |
| `--threads-per N` | Cores per STAR/featureCounts job | 8 |
| `--outdir PATH` | Output root override | `/mnt/data/FLASHseq_output` |
| `--gtf PATH` | GTF annotation (required for custom species; optional override for human/mouse) | built-in |
| `--star-index DIR` | Prebuilt STAR index (custom species) | — |
| `--ref-fasta PATH` | Genome FASTA; the STAR index is built from this + `--gtf` | — |
| `--ref-name LABEL` | Gene-set label shown in the log/summary | `custom` |
| `--sjdb-overhang N` | STAR `--sjdbOverhang` when building an index | 74 |
| `--skip-fastqc` | Omit the FastQC stage | — |
| `--skip-star` | Reuse an existing STAR BAM; skip alignment | — |
| `--count-only` | BAMs already exist — just featureCounts + merge | — |
| `--no-resume` | Never skip finished steps (full redo) | resume on |
| `--dry-run` | Print the planned commands, execute nothing | — |

> **Sizing:** total cores in flight = `--max-jobs × --threads-per`. On a
> 32-core machine, 4 × 8 is comfortable; 10 × 8 over-subscribes the box and
> can exhaust memory. The script warns when free memory drops below ~30 GB.
> Use `nohup … &` for unattended runs.

### Example

```bash
nohup bash flashseq_pipeline.sh \
  --species human \
  --fastq-dir /mnt/PTEGraw/…/AVT0551/…/Samples/<project>_Jun26/ \
  > pipeline.out 2>&1 &
```

### Outputs

Each invocation gets its own timestamped run directory:

```
<outdir>/scrnaseq/<AVT>_<species>_<YYYYMMDD_HHMM>/
├── fastqc/               # FastQC zip reports (one per R1/R2 FASTQ)
├── bam/                  # <label>_Aligned.sortedByCoord.out.bam + STAR logs
├── featurecounts/        # <label>.featureCounts.txt per sample
├── matrix/
│   ├── counts_matrix.tsv    # merged genes × samples count matrix
│   └── genes_per_cell.csv   # per-cell genes detected at >=1 and >=5
└── logs/
    ├── pipeline_run.log     # master run log + final SUMMARY block
    ├── _gates.csv           # per-sample unique-mapped % (QC gate)
    ├── _failed.samples      # samples that failed or were skipped
    └── _merge.R             # generated R merge script (reproducible)
```

A clean run ends with a `=== SUMMARY ===` block in `logs/pipeline_run.log`;
any shortfall is listed in `logs/_failed.samples`.

---

## 2. `flashseq_rseqc.sh` — per-sample RSeQC QC

Runs RSeQC QC on the STAR BAMs of an **existing** pipeline run — no FASTQ,
no alignment, no counting. For each sample it runs:

- **`geneBody_coverage.py`** — 5′→3′ gene-body coverage (tabular
  100-percentile matrix, ranked by skewness) — reference = the per-species
  *housekeeping-gene* BED (`--refbed` to override)
- **`read_distribution.py`** — reads per region (exon, 5′/3′ UTR, intron,
  intergenic; TSS/TES ±10/500 bp) — reference = the genome-wide GTF baked to
  a BED12 once per run (`--gtf` to override)

The script discovers every BAM under `<run dir>/bam/`, auto-creates a missing
`.bam.bai` with `samtools index`, is resume-safe (samples with existing
outputs are skipped), and is soft-fail per sample (a QC failure is listed in
`logs/_failed.samples` and never blocks the rest).

### Usage

```bash
bash flashseq_rseqc.sh --species <human|mouse|name> --outdir <RUN DIR> [options]
```

`--outdir` is the **exact** existing pipeline run directory (it must contain
`bam/`).

### Options

| Flag | Meaning | Default |
|---|---|---|
| `--species` | `human` / `mouse` (built-in references) or a custom name | — |
| `--outdir` | Existing pipeline run directory (required) | — |
| `--steps` | `both` \| `gbcov` \| `rdist` \| any comma combination | `gbcov,rdist` |
| `--refbed PATH` | Ready-made housekeeping-gene BED12 for gene-body coverage | per-species default |
| `--gtf PATH` | GTF for `read_distribution` (baked to BED12) | per-species default |
| `--samples a,b` | Process only these sample labels | all found |
| `--max-jobs N` | Samples in parallel; raise to ~100 on an idle machine to speed up the large merged plot | 4 |
| `--per-sample` | Per-sample gene-body coverage (table + plots per sample) instead of the merged combined plot | merged |
| `--gbcov-format FMT` | Merged plot format: `pdf` \| `png` \| `jpeg` \| `html`. `html` writes self-contained interactive browser pages (hover a curve to read its name, clickable legend — the legible choice at hundreds of samples) and needs no R | `pdf` |
| `--dry-run` | Print the planned commands, execute nothing | — |

> `--gbcov-format` applies to the **merged** plot only. In `--per-sample`
> mode the plots are produced by RSeQC's own R step (PDF).
> `geneBody_coverage.py` is the slow stage (several minutes per sample); the
> merged plot is built in parallel across samples, so `--max-jobs` controls
> the wall time.

### Example

```bash
# both tools, merged plot (default) — all samples of an existing run
nohup bash flashseq_rseqc.sh \
  --species human \
  --outdir /mnt/data/FLASHseq_output/scrnaseq/AVT0551_human_20260910_0930 \
  --max-jobs 10 > rseqc.out 2>&1 &

# interactive HTML merged plot instead of PDF
bash flashseq_rseqc.sh --species human --outdir <RUN DIR> --steps gbcov \
  --gbcov-format html

# per-sample tables + plots
bash flashseq_rseqc.sh --species human --outdir <RUN DIR> --per-sample
```

### Outputs (under the existing `<run dir>`)

| Path | Content |
|---|---|
| `rseqc/geneBodyCombined.geneBodyCoverage.txt` | Merged gene-body coverage table (all samples, 100 percentiles) |
| `rseqc/geneBodyCombined.heatMap.{png,pdf,html}` | Merged heatmap (≥ 3 samples), ranked by skewness |
| `rseqc/geneBodyCombined.curves.{png,pdf,html}` | Merged coverage curves |
| `rseqc/coverage/<label>/` | Per-sample mode: `<label>.geneBodyCoverage.txt` + plots per sample |
| `rseqc/<label>.readDistribution.txt` | `read_distribution.py` report per sample |
| `rseqc/<gtf basename>.bed12` | BED12 reference baked once per run (shared) |
| `logs/<label>.rseqc.log` | Full per-sample RSeQC stdout/stderr |
| `logs/rseqc_run.log` | Master RSeQC run log + SUMMARY |
| `bam/*.bam.bai` | Index auto-created for any BAM that lacked one |

---

## Resuming and re-running

- **Pipeline:** resume is on by default — FastQC/STAR/featureCounts steps are
  skipped when their outputs exist. `--no-resume` forces a full redo.
- **RSeQC:** resume is always on — a sample is skipped when its outputs
  already exist. Delete a sample's `rseqc/` outputs (or its table) to force a
  redo.
- **Partial re-count:** `--count-only` skips FastQC and STAR and goes straight
  to featureCounts + matrix merge (BAMs must already be in `<run dir>/bam/`).

## Common pitfalls

- `--outdir` means different things: **pipeline** = the root under which a
  timestamped run directory is created; **RSeQC** = the exact existing run
  directory (must already contain `bam/`).
- This pipeline counts by `gene_id` (`-g gene_id -t exon` with
  `--fracOverlap 0.25`, paired-end) — do not mix in a single-end or
  other counting convention for these matrices.
