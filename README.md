# GECF Scripts

Scripts developed at the **Gene Expression Core Facility (GECF)** for Oxford Nanopore Technology (ONT) data processing. These scripts cover the full workflow from raw POD5 files to basecalled, demultiplexed, aligned reads, plus custom 10x Genomics reference/template generation and qPCR plate validation.

---

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Scripts Overview](#scripts-overview)
- [dorado_basecaller_prod.sh](#dorado_basecaller_prodsh)
- [format_sample_sheet.sh](#format_sample_sheetsh)
- [demux.sh](#demuxsh)
- [get_stats.sh](#get_statssh)
- [adaptive_sampling_filter.sh](#filter_adaptive_samplingsh)
- [mapping_bams.sh](#mapping_bamssh)
- [update_dorado.sh](#update_doradosh)
- [cellranger_pipelines/ — 10x Genomics CellRanger](#cellranger_pipelines--10x-genomics-cellranger)
- [FLASHseq/ — FLASH-seq (single-cell) RNA-seq](#flashseq--flash-seq-single-cell-rna-seq)
- [fasta_to_gtf.py](#fasta_to_gtfpy)
- [qpcr-template-validator.py](#qpcr-template-validatorpy)
- [Test Fixtures](#test-fixtures)
- [Supporting Files](#supporting-files)
- [Typical Workflows](#typical-workflows)

---

## Requirements

- Linux system
- [Conda](https://docs.conda.io/en/latest/)
- [Dorado](https://github.com/nanoporetech/dorado) (installed separately, see below)

All other dependencies are managed through the provided conda environment.

---

## Installation

**1. Clone the repository**

```bash
git clone https://github.com/danielcav/GECF_scripts.git
cd GECF_scripts
```

**2. Create the conda environment**

```bash
conda env create -f basecaller_env.yml
conda activate basecaller
```

**3. Install Dorado**

Dorado is not available via conda and must be installed separately. You can use the provided update script:

```bash
bash update_dorado.sh
```

Or download it manually from the [Dorado releases page](https://github.com/nanoporetech/dorado/releases).

**4. Make scripts executable**

```bash
chmod +x *.sh
```

> All scripts must be kept in the same folder to work correctly together.

---

## Scripts Overview

| Script | Description |
|---|---|
| `dorado_basecaller_prod.sh` | Basecalls raw pod5 files using Dorado |
| `format_sample_sheet.sh` | Generates a sample sheet for demultiplexing |
| `demux.sh` | Demultiplexes barcoded basecalled BAM files |
| `get_stats.sh` | Computes sequencing statistics from a BAM file |
| `adaptive_sampling_filter.sh` | Filters pod5 files using adaptive sampling decisions |
| `mapping_bams.sh` | Aligns FASTQ/BAM reads to a reference genome with minimap2 |
| `update_dorado.sh` | Downloads and installs the latest stable Dorado release |
| `cellranger_pipelines/` | 10x Genomics CellRanger pipeline scripts (gex, atac, arc) |
| `cellranger_pipelines/common/` | Custom-gene (transgene) expression report — per-gene status, UMI/read stats, plots |
| `FLASHseq/` | FLASH-seq (single-cell) RNA-seq pipeline + RSeQC QC |
| `fasta_to_gtf.py` | Converts custom multi-FASTA to cellranger-arc-compatible GTF |
| `qpcr-template-validator.py` | Validates qPCR DNA/PRIMER Hamilton Excel templates before sending to robot |

---

## dorado_basecaller_prod.sh

Basecalls raw pod5 files using Dorado. Supports simplex and duplex modes, modified base calling, and both barcoded and non-barcoded runs. Automatically downloads the selected model and generates QC reports and sequencing statistics.

### Usage

```bash
bash dorado_basecaller_prod.sh [--qscore <int>] [--legacy]
```

### Options

| Option | Description | Default |
|---|---|---|
| `--qscore <int>` | Minimum qscore filter for pass/fail read splitting | 10 |
| `--legacy` | Show older model versions in model selection | false |
| `-h`, `--help` | Show help message | |

### The script will interactively ask for

- Project name (run name)
- Path to the pod5 directory
- Basecalling mode (simplex or duplex)
- Accuracy mode (fast, hac, sup)
- Model selection
- Whether to call modified bases (5mC, 5hmC, 6mA...)
- Whether reads are barcoded

### Output structure

```
SUP_basecalling_<run_name>/
├── ubam_files/          # Basecalled BAM files (pass and fail)
├── fastq_files/         # Compressed FASTQ files
├── summaries/           # Sequencing summary files
├── QC_reports/          # Sequali QC reports (HTML + JSON)
├── logs/                # Log file with parameters, versions and stats
└── checksums_<run>.txt  # MD5 checksums for BAM and FASTQ files
```

### Example

```bash
conda activate basecaller
bash dorado_basecaller_prod.sh --qscore 10
```

> If reads are barcoded, the script will stop after basecalling and ask you to run `demux.sh` next.

---
 
## format_sample_sheet.sh
 
Generates a properly formatted sample sheet CSV file required by `demux.sh`. It takes a simple text file mapping barcodes to sample names and enriches it with run metadata (kit name, experiment name, sequencer position), either by auto-detecting them from a pod5 file or by manual input.
 
> Always run this script before `demux.sh`.
 
### Usage
 
```bash
bash format_sample_sheet.sh <barcode_aliases.txt>
```
 
### Input file format
 
A plain text file with one entry per line, formatted as `barcode;sample_name`:
 
```
barcode01;sample1
barcode02;sample2
barcode03;HEK293
barcode04;HeLa
```
 
### The script will interactively ask for
 
- Whether to auto-detect metadata from a pod5 file (recommended)
- If auto-detect: path to a pod5 file or directory
- If manual: kit name, experiment name, and sequencer position (A or B)
### Output
 
A CSV sample sheet saved in the current directory:
 
```
sample_sheet_<experiment_name>.csv
```
 
With the following structure:
 
```
experiment_id,kit,position_id,barcode,alias
SOL0046,SQK-NBD114-24,P2S-00697-A,barcode01,sample1
SOL0046,SQK-NBD114-24,P2S-00697-A,barcode02,sample2
```
 
### Example
 
```bash
bash format_sample_sheet.sh barcodes_SOL0046.txt
# When prompted, choose auto-detect and provide the pod5 directory
```
 
---

## demux.sh

Demultiplexes barcoded BAM files produced by `dorado_basecaller_prod.sh`. Handles both simplex and duplex data, including false positive duplex read detection and removal. Produces submission-ready folder structure with QC reports and statistics for each sample.

### Usage

```bash
bash demux.sh --input <file.bam> --sample-sheet <file.csv> [options]
```

### Options

| Option | Description | Default |
|---|---|---|
| `--input <file.bam>` | Input BAM file (required) | |
| `--sample-sheet <file>` | Sample sheet CSV or TXT file (required) | |
| `--qscore <int>` | Minimum qscore filter | 10 |
| `--greedy` | Classify duplex reads using only one parent if the other is unclassified | false |
| `--barcode-both-ends` | Require barcodes on both ends of the read | false |
| `-h`, `--help` | Show help message | |

### Sample sheet format

The sample sheet must be a CSV or TXT file. It can be generated using the provided `format_sample_sheet.sh` script. It should contain at minimum the run name and kit name.

### Output structure

```
demultiplexed_data_<run_name>/
├── ubam_files/          # Per-sample BAM files (all, pass, fail)
├── fastq_files/         # Per-sample compressed FASTQ files
├── summaries/           # Per-sample sequencing summaries
├── QC_reports/          # Per-sample Sequali QC reports
└── logs/                # Log file with parameters and stats
```

### Example

```bash
bash demux.sh \
    --input /data/SUP_basecalling_p2sid0046/ubam_files/p2sid0046_untrimmed.all.bam \
    --sample-sheet /data/sample_sheet_p2sid0046.csv \
    --qscore 10
```

---

## get_stats.sh

Computes sequencing statistics directly from a BAM file. Called automatically by `dorado_basecaller_prod.sh` and `demux.sh`, but can also be used standalone on any BAM file.

### Usage

```bash
bash get_stats.sh [-o output_log_file] <your_file.bam>
```

### Options

| Option | Description |
|---|---|
| `-o <file>` | Append statistics to a log file (optional) |
| `-h`, `--help` | Show help message |

### Output metrics

- Total number of reads
- Total yield (Gb)
- Mean and median read length
- Read length N50
- Mean and median read quality
- Number of reads greater than 40kb and 100kb
- Number of active channels
- Top 5 longest reads with their quality scores

### Example

```bash
# Print stats to terminal only
bash get_stats.sh /data/my_sample.bam

# Print stats to terminal and append to a log file
bash get_stats.sh -o my_run.log /data/my_sample.bam
```

---

## filter_adaptive_sampling.sh

When sequencing with adaptive sampling, MinKNOW generates a decisions file (`AS_decisions_*.csv`) that records which reads were accepted or rejected in real time. This script extracts the IDs of accepted reads and creates a filtered pod5 file containing only those reads, significantly reducing basecalling time.

### Usage

```bash
bash filter_adaptive_sampling.sh \
    --as-csv <AS_decisions.csv> \
    --pod5-dir <pod5_directory> \
    --project <project_name> \
    [--output-dir <output_directory>]
```

### Options

| Option | Description | Default |
|---|---|---|
| `--as-csv <file>` | Path to the adaptive sampling decisions CSV file (required) | |
| `--pod5-dir <dir>` | Path to the directory containing pod5 files (required) | |
| `--project <name>` | Project name, e.g. p2sid0046 (required) | |
| `--output-dir <dir>` | Output directory | Parent of pod5 directory |
| `-h`, `--help` | Show help message | |

### Output

- `accepted_read_ids_<project>.txt` — list of accepted read IDs
- `accepted_reads_<project>.pod5` — filtered pod5 file containing only accepted reads

> Some read IDs may be missing from the pod5 files. This is normal behaviour with adaptive sampling and is handled automatically with the `--missing-ok` flag.

### Example

```bash
bash filter_adaptive_sampling.sh \
    --as-csv /data/adaptive_sampling/AS_decisions_SOL0046.csv \
    --pod5-dir /data/SOL0046/pod5/ \
    --project p2sid0046
```

Then basecall the filtered pod5 file:

```bash
bash dorado_basecaller_prod.sh
# When prompted, point to the filtered pod5 file
```

---

## mapping_bams.sh

Aligns FASTQ or BAM reads to a reference genome using **minimap2**, producing sorted and indexed BAM files. Processes single files or entire directories. The script also prints a `samtools flagstat` summary for each file.

### Usage

```bash
bash mapping_bams.sh -r <reference.fa> -i <input_file_or_folder> [-p <pattern>]
```

### Options

| Option | Description | Default |
|---|---|---|
| `-r <path>` | Reference genome in FASTA format (required) | |
| `-i <path>` | Path to a single file or directory (required) | |
| `-p <pattern>` | File pattern/extension to match (e.g. `*.pass.bam`, `*.fastq.gz`) | `*.fq, *.fastq, *.fq.gz, *.fastq.gz, *.bam` |
| `-h` | Show help message | |

### Supported inputs

- FASTQ / `.fq` (gzipped or uncompressed)
- BAM — streamed through `samtools fastq` first, then aligned

### Output

For each input file, a sorted and indexed BAM (`<sample>.sorted.bam`) plus its `.bai` index. A flagstat summary is printed to the terminal.

> Uses all available CPU cores via `nproc`.

---

## update_dorado.sh

Automatically downloads and installs the **latest stable** Dorado release from GitHub. Checks for pre-release builds and skips them by default. Detects your OS/architecture automatically.

### Prerequisites

Dorado binaries are platform-specific — this script detects your OS (macOS / Linux) and architecture (x64 / ARM) so you get the right build:

```bash
conda activate basecaller
bash update_dorado.sh
```

> Dorado is required before you can run basecalling (`dorado_basecaller_prod.sh`). The script is also part of the installation process; see [Installation](#installation).

---

## cellranger_pipelines/ — 10x Genomics CellRanger

All 10x Genomics pipeline scripts live in `cellranger_pipelines/`, one subfolder per assay. Every pipeline follows the same conventions:

| Pipeline | Folder | Tool | Archive destination |
|----------|--------|------|---------------------|
| atac (single-cell ATAC) | `cellranger_pipelines/atac/` | cellranger-atac count | `/AA/NGSruns/ZZ_10XG/atac/` |
| gex (single-cell RNA) | `cellranger_pipelines/gex/` | cellranger count | `/AA/NGSruns/ZZ_10XG/count/` |
| arc (Multiome RNA+ATAC) | `cellranger_pipelines/arc/` | cellranger-arc count | `/AA/NGSruns/ZZ_10XG/arc/` |

- **Standard archive rule:** results are rsync'd to `/AA/NGSruns/ZZ_10XG/<modality>/<run folder>/` — a shared, sequencer-agnostic root (created 2026-11-24; mirrors the former `AVT/AV240401/ZZ_10XG/` layout, so historical AVT results stay in place and new results from any sequencer land in the same tree).
- **Generator-style pipelines** (gex, atac, arc): a Python generator reads the sample Excel/datasheet, writes one ready-to-run `.sh` into the **current working directory**; review it, then run it — results land locally there and are deployed to the archive automatically.
- **Legacy drivers removed from this repo (2026-11-24):** the multi (TotalSeq/VDJ) driver and the old arc bash driver `CRarc_patched_full` used to live here; they are no longer maintained in this repo. Server-side copies of those driver workflows under `/home/gecf/NGS_pipelines/CR_pipeline/` are untouched and still available there.

### atac — generate_cellranger_template_atac.py

Single-cell ATAC pipeline. Given an AVT Excel tracking file, generates a ready-to-run shell script that executes `cellranger-atac count` for every sample and writes it out as `cellranger-atac_AVT0XXX.sh` into the current directory.

```bash
python3 cellranger_pipelines/atac/generate_cellranger_template_atac.py "AVT0571_Arianna_GECF Copy.xlsx"
```

The prompt flow (AVT run, sample selection / add sheet-less samples, reference genome) and the run-wide + per-sample settings work the same as GEX and ARC — tweak them in the generated script's editable blocks, then run it. By default each sample's output folder is numbered in sample order — `01_<sample>`, `02_<sample>`, … — so the folders stay in order even when the names are not alphabetical; add `--no-number-prefix` to the generator command to keep bare sample names instead. With several Excel inputs, the sample order is the common order of all files when they agree on it, otherwise the first file's order (samples appearing only in later files are appended). Full step-by-step: see [ATAC_pipeline_howto.html](cellranger_pipelines/atac/ATAC_pipeline_howto.html). The Excel sheet needs a sample-name column plus a per-sample fastq ID (or full FASTQ path); the AVT run can come from the sheet, the file name, or the prompt.

### gex — generate_cellranger_template_gex.py

Single-cell RNA (GEX) pipeline. Given an AVT Excel tracking file, generates a ready-to-run shell script that executes `cellranger count` for every sample and writes it out as `cellranger-count_AVT0XXX.sh` into the current directory.

```bash
python3 cellranger_pipelines/gex/generate_cellranger_template_gex.py "AVT0571_Arianna_GECF Copy.xlsx"
```

Tweak the run-wide and per-sample settings in the generated script, then run it. By default each sample's output folder is numbered in sample order — `01_<sample>`, `02_<sample>`, … — so the folders stay in order even when the names are not alphabetical; add `--no-number-prefix` to the generator command to keep bare sample names instead. With several Excel inputs, the sample order is the common order of all files when they agree on it, otherwise the first file's order (samples appearing only in later files are appended). The script also measures the TSO% library-quality metric on the Read-2 fastq of **every** avid run of every sample and writes one central `tso_report.txt` (one section per avid run) to the archive. Full step-by-step: see [GEX_pipeline_howto.html](cellranger_pipelines/gex/GEX_pipeline_howto.html).

### arc — generate_cellranger_template_arc.py

Multiome (RNA + ATAC) pipeline. Given an AVT Excel tracking file, generates a ready-to-run shell script that executes `cellranger-arc count` for every sample (one `--libraries` CSV with the GEX leg before the ATAC leg) and writes it out as `cellranger-arc_AVT0XXX.sh` into the current directory. Legs are paired by sample name: `Sample_GEX` + `Sample_ATAC` become one sample (one `cellranger-arc count`), and a name without a leg suffix counts as GEX — see the howto for the full naming rule and the sample-selection prompt.

```bash
python3 cellranger_pipelines/arc/generate_cellranger_template_arc.py "AVT0571_Arianna_GECF Copy.xlsx"
```

Tweak the run-wide and per-sample settings in the generated script, then run it. At the end of the run the script writes a `<avt>_stats.xlsx` into the run folder — one sheet per modality (GEX, ATAC), one row per sample per modality (official sample name, GECF library label, seq. run, fastq code, cellranger-arc version, genome assembly, intron setting) — and ships it to the archive with the rest. The TSO% report is collected the same way as in GEX: one section per GEX leg, i.e. **every** avid run with a GEX leg is measured (ATAC-only samples are skipped). Non-fatal fallbacks: without an Excel library it writes the same thing as two CSVs, and without python3 it warns and skips. Full step-by-step: see [ARC_pipeline_howto.html](cellranger_pipelines/arc/ARC_pipeline_howto.html).

### common — custom_gene_report.py

Custom-gene (transgene) expression report (v1.1.0). Given a cellranger `count`
output and a list of custom-gene names, it writes one **self-contained HTML report**
(no CDN, fully offline) with one stats table per sample, per-gene plots (UMAP + UMI
histogram + total-UMI box) and a read-start coverage track; `--project-dir` mode
adds an "Across samples" comparison block.

```bash
python3 cellranger_pipelines/common/custom_gene_report.py --sample-dir <sample_dir> --genes genes.txt [-o OUT_FOLDER]
python3 cellranger_pipelines/common/custom_gene_report.py --project-dir <project_dir> --genes genes.txt [--min-umi N] [--bam BAM] [--per-file N]
python3 cellranger_pipelines/common/custom_gene_report.py --help
```

It is run under a plain `python`: if the interpreter you launched it with lacks
`plotly`/`numpy`/`scipy`, the script transparently re-runs itself under the lab's
conda `brbseq` env (then `basecaller`) — nothing to activate. Set
`CGR_NO_ENV_SWITCH=1` to keep the current interpreter. A long gene list splits into
an index + `…_part01.html`, `…_part02.html`, … (`--per-file`, default 50; `0` = one
file). The main file is always `custom_genes_<input-dir-name>.html`.

**Where every value comes from.** The cell-count and UMI columns are read entirely
from `outs/filtered_feature_bc_matrix/`: `features.tsv` (the feature table),
`barcodes.tsv` (the set of cells called in the sample), and `matrix.mtx` (the sparse
gene×cell UMI matrix). Two **optional** inputs enrich the report; when absent the
affected columns/panels degrade to `n/a` (the matrix columns always work):
- **BAM** — `possorted_genome_bam.bam` / `possorted_bam.bam` /
  `gex_possorted_bam.bam` in the sample dir, or `--bam`. Source of the
  *reads (called cells)* column and the read-start coverage track (needs
  `samtools`; resolved from PATH then the known lab envs).
- **UMAP projection CSV** — a `barcode,UMAP-1,UMAP-2` file next to the analysis
  output, auto-detected. Source of the UMAP panel.

**Transgene row vs. endogenous row — the key rule.** For each gene you list, the
script resolves its row(s) in `features.tsv`:
- **Transgene row** = the feature whose **ID equals the gene name** (the custom
  contig you appended; cellranger keeps the contig name as the feature ID).
  *All* detection/UMI statistics below are read from this row.
- **Endogenous row** = any *other* feature with the **same name but a different ID**
  (e.g. an Ensembl ID carrying the gene's display name) — this is the
  *endogenous cells / endogenous UMIs* context (expression of the same-named
  endogenous gene).
- **Status** (judged on the transgene row): `EXPRESSED` = row present and ≥1 cell
  reaches the detection threshold; `SILENT` = row present but 0 cells reach it;
  `NO-ROW` = no feature with an ID equal to the gene name exists (the gene was not
  in the reference used for that run).

**How each reported value is computed.** "Detected" always means a cell's UMI for
that gene ≥ threshold (`--min-umi`, default 3). `cells_total` = the number of
barcodes in `barcodes.tsv`.

| Column | Computed as | Source |
|---|---|---|
| Detected cells (`cells ≥ N UMI`) | # cells whose transgene-row UMI ≥ threshold | matrix |
| `% cells` (detected) | detected cells ÷ `cells_total` | matrix + barcodes.tsv |
| `cells ≥ 1` / `cells ≥ 3` UMI | # cells at that fixed cutoff — a sensitivity check independent of the threshold; the column equal to `--min-umi` is dropped so a count is never shown twice | matrix |
| `total UMIs` | sum of the transgene row across all cells | matrix |
| `reads (called cells)` | # **primary** alignments on the gene's contig carrying the gene's `GN:` tag that come from cells in `barcodes.tsv` — **reads, not UMIs** | BAM (`samtools`); `n/a` if no BAM |
| `mean UMI` | `total UMIs` ÷ # expressing cells (cells with ≥1 UMI — expressing only, not all cells) | matrix |
| `median UMI` | median of per-cell UMI over expressing cells (≥1 UMI) | matrix |
| `endogenous cells` | # cells with UMI > 0 in an endogenous (same-name) row | matrix |
| `endogenous UMIs` | sum of UMI across the endogenous (same-name) rows | matrix |

**Plots (per sample).** One gene selector re-renders the three figures in place
(one figure of each type, never one per gene): the UMAP (colour = that cell's UMI
for the selected gene, log10 scale, coordinates from the UMAP CSV); a histogram of
per-cell UMI of expressing cells (≥1 UMI, log bins); and a box plot of **total
transcriptome** UMI per cell, split into detected vs not-detected for the selected
gene. The **read-start coverage** panel is `samtools depth` along the gene's
transgene contig (one per sample).

Full step-by-step: see
[custom_gene_report_howto.html](cellranger_pipelines/common/custom_gene_report_howto.html).

---

## FLASHseq/ — FLASH-seq (single-cell) RNA-seq

Run-to-result pipeline for **Alithea FLASH-seq** single-cell RNA-seq: FASTQ →
FastQC → STAR → featureCounts → one merged genes × samples count matrix, plus
standalone per-sample RSeQC QC (gene-body coverage + read distribution) on the
STAR BAMs. Resume-safe; a failing sample never blocks the others.

| File | Purpose |
|---|---|
| `FLASHseq/flashseq_pipeline.sh` | Core pipeline (FastQC → STAR → featureCounts → merged matrix) |
| `FLASHseq/flashseq_rseqc.sh` | Per-sample RSeQC QC on existing BAMs |
| `FLASHseq/gbcov_redraw.py` | Merged gene-body plot: PDF/PNG/JPEG (via R) or interactive HTML (via plotly) |
| `FLASHseq/flashseq-env.yaml` | Conda environment definition |

See [FLASHseq/README.md](FLASHseq/README.md) for full usage, options, outputs,
and worked examples.

### Quick start

```bash
conda env create -f FLASHseq/flashseq-env.yaml   # once, if you don't have it
conda activate flashseq

# 1) the full pipeline (FastQC -> STAR -> featureCounts -> merged matrix)
nohup bash FLASHseq/flashseq_pipeline.sh \
  --species human \
  --fastq-dir /mnt/PTEGraw/…/AVT0551/…/Samples/<project>_Jun26/ \
  > pipeline.out 2>&1 &

# 2) RSeQC QC on the finished run (interactive HTML plot at hundreds of samples)
bash FLASHseq/flashseq_rseqc.sh \
  --species human \
  --outdir /mnt/data/FLASHseq_output/scrnaseq/AVT0551_human_20260910_0930 \
  --steps gbcov --gbcov-format html
```

> `gbcov_redraw.py` must stay in the **same directory** as
> `flashseq_rseqc.sh`.

---

## fasta_to_gtf.py

Converts a multi-FASTA file (one contig per custom gene) into four GTF lines (gene + transcript + exon + CDS). The output pairs with your FASTA for building a **cellranger-arc** reference (`mkref`). Also writes a cleaned FASTA matching all names used in the GTF.

Contig names longer than **32 characters** are rejected with an error before any output is written. The cellranger-arc 2.x ATAC signal-track stage parses bedgraph rows with a 35-character field, so names over 35 are silently truncated and the run later fails with a misleading "invalid contig name" error. The 32-character cap keeps names safely under that limit.

### Usage

```bash
python3 fasta_to_gtf.py input.fasta output.gtf
```

### Output

| Path | Purpose |
|---|---|
| `<output>.gtf` | 4 lines per FASTA record (+ strand) |
| `<output>.clean.fasta` | Cleaned FASTA headers matching GTF contig names |

> Names are sanitized: underscores become hyphens, `|` in headers drops metadata after the first pipe. Names must be at most 32 characters after sanitization, or the script aborts with an error.

---

## qpcr-template-validator.py

Validates DNA and PRIMER Excel workbooks used as input templates for the **Hamilton qPCR workstation** (**v4.03 Hamilton Workstation**) before sending them to the robot. Catches errors early so plate preparation runs smoothly.

The validator checks both **DNA list** and **PRIMER list** sheets (a workbook must contain one or the other, never both). It verifies headers, sample/primer names, tasks, quantities, controls, plate capacity limits, duplicate entries, TaqMan reporter types, and more.

Two independent files work together:

- `qpcr_validation.py` — pure validation core (no UI imports)
- `qpcr-template-validator.py` — dashboard wrapper with **GUI** and **CLI** modes; loads the core from its sibling path

### Usage

**CLI mode** (headless or scripted):

```bash
# Single file
python qpcr-template-validator.py --cli /path/to/DNA_list.xls

# Multiple files
python qpcr-template-validator.py --cli file1.xls file2.xlsx

# Scan a directory recursively for Excel workbooks
python qpcr-template-validator.py --cli --dir /path/to/workbook_folder/
```

**GUI mode:**

```bash
# File picker popup
python qpcr-template-validator.py

# Pre-load files into the dashboard
python qpcr-template-validator.py file1.xls file2.xlsx
```

### Requirements

- `openpyxl` (for `.xlsx`/`.xlsm`)
- `xlrd` (for legacy `.xls`)

### Output

Exit code:

| Code | Meaning |
|---|---|
| 0 | All workbooks passed |
| 1 | One or more workbooks has failures |

Results are printed immediately. Each validation checks:

- Correct sheet name (`DNA` or `PRIMER`)
- Valid headers
- Task values (STND, NTC, UNKN)
- Sample/primer name format and uniqueness
- Quantity validation (NTC = 0, no decimal standards)
- Presence of controls (NTC + standard curve)
- Plate capacity limits (192 rows max for DNA, 64 unique primers)
- TaqMan reporter type

---

## Test Fixtures

A set of **20** pre-built test Excel workbooks live in `test_fixtures/`, covering every validation rule. To install them:

```bash
cd test_fixtures
python create_test_fixtures.py
```

Then verify the validator against all fixtures:

```bash
python run_fixture_checks.py
```

See `test_fixtures/README.md` for details on each individual fixture case. Fixtures are not tracked in git (see `.gitignore`).

---

## Supporting Files

| Path | Purpose |
|---|---|
| `basecaller_env.yml` | Conda environment definition for Dorado basecalling (pre-built models excluded) |
| `.gitignore` | Excludes large artifacts: `.xlsx` test fixtures, `.local/`, and generated BAMs |
| `test_fixtures/` | Test fixture generators and runner; see [Test Fixtures](#test-fixtures) above |

---

## Typical Workflows

### Standard simplex or duplex run (non-barcoded)

```bash
conda activate basecaller
bash dorado_basecaller_prod.sh
```

### Barcoded run (simplex or duplex)

```bash
conda activate basecaller
bash dorado_basecaller_prod.sh
# Script stops after basecalling for barcoded data
bash demux.sh --input <basecalled.bam> --sample-sheet <sheet.csv>
```

### Adaptive sampling run

```bash
conda activate basecaller
bash filter_adaptive_sampling.sh \
    --as-csv <AS_decisions.csv> \
    --pod5-dir <pod5/> \
    --project <run_name>
bash dorado_basecaller_prod.sh
# When prompted, point to the filtered pod5 file
```

### Compute stats on any BAM file

```bash
conda activate basecaller
bash get_stats.sh /path/to/any_file.bam
```

### FLASH-seq (single-cell) RNA-seq run

```bash
conda activate flashseq
bash FLASHseq/flashseq_pipeline.sh \
  --species human \
  --fastq-dir /mnt/PTEGraw/…/AVT0551/…/Samples/<project>_Jun26/

# then, once the pipeline has finished, RSeQC QC on the run:
bash FLASHseq/flashseq_rseqc.sh \
  --species human \
  --outdir /mnt/data/FLASHseq_output/scrnaseq/AVT0551_human_20260910_0930 \
  --gbcov-format html
```