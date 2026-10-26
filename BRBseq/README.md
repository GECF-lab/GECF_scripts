# BRB-seq (MERCURIUS DRUG-seq) single-pool RNA-seq pipeline

`brbseq_pipeline.sh` (v1.1.0) — end-to-end pre-processing for
**one pool/run** of MERCURIUS 3' RNA-seq (BRB-seq / DRUG-seq) data:
barcode prep → STARsolo alignment + gene/UMI counts → mandatory duplication
statistics → sample-named matrices → FastReadCounter read-distribution report
→ optional per-sample demux.

The STAR request always requests **both** UMI dedup methods
(`--soloUMIdedup NoDedup 1MM_Directional` in one pass) so the mandatory
duplication step (and the two delivered matrices) can be produced from a
single alignment.

## Contents

| File | Purpose |
|---|---|
| `brbseq_pipeline.sh` | Core single-pool pipeline (this README documents v1.1.0) |
| `build_star_index.sh` | One-time STAR genomeDir build (Mercurius guide 1.1; no Ensembl fetch) |
| `brbseq-env.yaml` | Conda environment (`brbseq`) — STAR, samtools, fastqc, picard, fqtk, RSeQC tools, python/openpyxl |
| `utils/prepare_barcodes.sh` | Builds `barcode_whitelist.txt` + demux maps from the manifest (TSV `Sample_ID` → `bar`) |
| `utils/parse_libprep_xlsx.py` | Reads the Mercurius lib-prep workbook (samples, plate wells, Kit ref → `bar`) |
| `utils/rename_matrix.py` | Renames STARsolo's `Gene/raw` barcode-matrix to `Sample_ID` columns (into `star/samples/` or `star/samples_all_wells/`) |
| `utils/prepare_umi_matrix.py` | Flattens `star/samples/umiDedup-*.mtx` into the two matrix CSVs (`umi_matrix`, `reads_matrix`) |
| `utils/get_duplication_statistics.py` | Per-well + pool duplication stats from NoDedup vs 1MM_Directional matrices |
| `utils/demultiplex_bam.sh` | Vendored picard-based BAM demux helper (one `FilterSamReads` per sample, `CR` tag) |

## Environment

The pipeline locates its tools **from the `brb`/`brbseq` Conda environment
first** — even when a working copy of STAR already sits on the server PATH —
so a stray old system Python can no longer hijack a run. Discovery order:

1. `BRBSEQ_ENV=<path>` (explicit override, highest priority);
2. the currently active Conda env, if it is named `brbseq` or `brb`;
3. a `brbseq` then `brb` env under the standard production anaconda dir.

In practice **no `conda activate` step is needed** — just run it:

```bash
cd /path/to/BRBseq
./brbseq_pipeline.sh --run-name Alireza ...
```

The env is **created** once if it doesn't yet exist (or point the pipeline at
an existing one with `BRBSEQ_ENV=<path>`):

```bash
conda env create -f brbseq-env.yaml   # creates env 'brbseq'
```

| Tool (from the env) | Why it's needed |
|---|---|
| `STAR` | aligns reads + CB/UMI extraction (BRB-seq `STAR_Solo`) |
| `fastqc` | read-quality QC |
| `samtools` | index/BAM manipulation |
| `read_distribution.py` | RSeQC read-length distribution |
| `picard` | `FilterSamReads` — per-sample BAM demux (writes `demux_bam/`) |
| `fqtk` | per-sample FASTQ demux (writes `demux_fastq/`; on by default, disable with `--skip-dmux-fq`) |
| `python` + `openpyxl` | workbook mode + UMI dedup (needs Python **3.7+**) |
| `java` (OpenJDK, shipped with the `brbseq` Conda env) | runs **FastReadCounter** (the jar is outside Conda — a production-system file at `/home/gecf/softwares/FastReadCounter/`) |

> **Why env-first?** Previously, with no `conda activate`, the script picked up
> the server's old system `python3` (< 3.7) for the workbook parser and crashed
> on `from __future__ import annotations`. Tool resolution is now env-first to
> prevent this class of breakage.

The pipeline hard-fails naming any missing **required** tool before it starts.
A missing **optional**-step tool (fastqc / picard / fqtk / RSeQC) logs a warning
and skips that step rather than aborting.

## Input contract

- **R1 barcode** read: 14 bp CB + 14 bp UMI (`--soloCBlen 14 --soloUMIlen 14`).
- **Pool label** is **inferred from the R1/R2 FASTQ filenames**, not from a
  flag: `basename` → cut at the first `_S<digit>` occurrence. E.g.
  `pool_poolA_S1_R1_001.fastq.gz` → pool `pool_poolA`. (No `--pool` flag, no
  `AVT` env override.)
- **Manifest** — exactly one of:
  - `--barcode-tsv <tsv>` (TSV of `Sample_ID<TAB>barcode`; a header row is
    auto-detected and skipped automatically when present), **or**
  - `--libprep-xlsx <xlsx>` — Mercurius lib-prep workbook; the barcode plate
    is joined by well, plate family inferred from the `Kit ref` row
    (V5D / V5A10) unless `--plate <file>` is given.
- `--species human|mouse|<custom>`; custom species require
  `--star-genome <genomeDir>` built with `build_star_index.sh`.

## Workbook mode (`--libprep-xlsx`) — always on, no flag needed

In lib-prep workbook mode:

* Every **unoccupied plate well** is kept in the barcode whitelist and named
  `empty_well_<WELL>` (reserved prefix: a real sample name may not start with
  `empty_well`).  This is the debugging/troubleshooting default so empty
  reads are never silently dropped by the STAR whitelist.
* The sample-named count matrix is therefore written **twice** from the SAME
  single STAR output (one alignment, no second STAR call):
  - `star/samples_all_wells/` — samples + `empty_well_*` (full plate)
  - `star/samples/` — samples only; **the final UMI/reads CSVs are built from
    THIS directory**
* Per-sample demux output (`demux_bam/`, `demux_fastq/`) contains **real
  samples only**: `empty_well_*` rows are filtered out of the demux maps
  before picard/fqtk see them (each picard `FilterSamReads` is a full BAM
  scan, so one per empty well would be pure waste).  Empty wells remain in
  the STAR whitelist and in `samples_all_wells/` for debugging.

## Read-distribution report (FastReadCounter)

After the matrices are built (optional demux still to come), the pipeline runs
**FastReadCounter** over the same STAR BAM and the same 96-barcode whitelist,
then derives two reports for the CLC / QC side:

- `files_for_report_and_CLC/<pool>_report_summary.csv` — a whole-plate grid.
  One row per well (A01…H12) plus a trailing `Unknown` row; columns
  `depth_mio, mapped_reads_mio, no_feature_mio, not_aligned_mio` (millions).
  This is the sequencing-quality table that reports how much of each sample's
  reads landed on the expected target vs. not.
- `files_for_report_and_CLC/<pool>_all_stats.csv` — the same per-sample
  numbers as raw integer counts (unrounded), including all FRC buckets
  (`__no_feature`, `__ambiguous`, `__too_low_aQual`, `__not_aligned`,
  `__alignment_not_unique`, total mapped).

This is a **port of the legacy BRB-seq pipeline's FRC + `report_stats` chain**
(FastReadCounter → rename columns to Sample_ID → two awk passes), so the
output is directly comparable with historical runs.  The report is
**mandatory**; it always uses the STAR BAM already on disk (no second
alignment), the gene annotation registered in the STAR index
(`sjdbGTFfile` from `genomeParameters.txt`), and the same 96-barcode
whitelist STAR used.  The jar itself is **outside** the Conda env — a
pre-installed Java program shipped on the production server; its path is the
single hardcoded line `FRC_JAR=` near the top of the FRC section of
`brbseq_pipeline.sh` (update this one line if the jar ever moves).

> **Note on cross-GTF comparability.**  `depth` and `not_aligned` are
> annotation-independent and reproduce historical values exactly.  The split
> between `mapped_reads` and `no_feature` is annotation-dependent: a pool
> re-run against a newer genome annotation (e.g. GRCm39.114 → GRCm39.116)
> will shift some reads between those two columns while leaving the two
> annotation-independent columns identical.  That is expected, not a bug.

## Output layout

Single-level leaf under `--out-dir`:

```
OUT_DIR/BRBseq_<run-name>_<YYYYMMDD_HHMMSS>/
├── barcodes/             barcode_whitelist.txt, manifest*.tsv,
│                         fastq_sample_metadata.tsv, bam_barcodes.tsv
├── star/
│   ├── STAR.log          STAR invocation + runtime output
│   ├── <pool>_Aligned.sortedByCoord.out.bam
│   ├── <pool>_Solo.out/  STARsolo's own output — Gene/raw/ holds the raw
│   │                     UMI matrix (its ONLY location, no mirror copies)
│   ├── samples/          Sample_ID-named matrix (SAMPLES ONLY — the input
│   │                     to the two matrix CSVs below)
│   └── samples_all_wells/  Sample_ID-named matrix (SAMPLES + empty_well_*)
│                         [lib-prep workbook mode only]
│                         (both also hold sample_barcode_map.tsv and the
│                          umiDedup-{NoDedup,1MM_Directional}.mtx matrices)
├── rseqc/                read_distribution (optional)
├── fastreadcounter/      FastReadCounter output (counts.detailed.txt, etc.)
├── files_for_report_and_CLC/
│   ├── <pool>_report_summary.csv   per-well read-distribution grid
│   │                               (depth/mapped/no-feature/no-alignment, in
│   │                               millions) — one row per plate well + Unknown
│   └── <pool>_all_stats.csv        same read-distribution, raw totals per sample
├── demux_bam/            one BAM per sample (optional)
├── demux_fastq/          demuxed FASTQs (optional)
├── duplication_statistics.csv
├── matrices/             the two deliverable matrices (canonical PCA/DE input
│   ├── <pool>_umi_matrix.csv    UMI-deduped — CANONICAL
│   └── <pool>_reads_matrix.csv  raw-read matrix (duplication context only)
├── pipeline.log          supervisor log (step status + exact commands)
├── logs/                 one output log per step
└── README.html           auto-generated folder index — describes every file
                          above; written at the end of the run, open in a
                          browser; safe to share with users as the orientation
                          point
```

## Changelog

Bumped in `SCRIPT_VERSION` in `brbseq_pipeline.sh` (currently **v1.1.0** —
the current release; the latest entry is the v1.1.0 FastReadCounter report).
v1.0.1 was reserved but never shipped — the version went straight from the
1.0.0 baseline to 1.0.2.

### v1.0.0 — baseline single-pool pipeline

First end-to-end version of the MERCURIUS (BRB-seq / DRUG-seq) single-pool
workflow. Everything below is the shipped baseline:

- **One STARsolo pass, both dedup methods** — `--soloUMIdedup NoDedup
  1MM_Directional` in a single alignment, so the duplication statistics and
  the two deliverable matrices all come out of the same run.
- **Matrix CSVs (mandatory) are the two deliverables** — both in the
  `matrices/` folder at the run root:
  - `matrices/<pool>_umi_matrix.csv` — **CANONICAL**, from
    `umiDedup-1MM_Directional.mtx` (UMI molecule counts; the correct input
    for PCA / differential expression).
  - `matrices/<pool>_reads_matrix.csv` — from `umiDedup-NoDedup.mtx` (raw
    read counts; kept for duplication context, **not** a PCA/DE input).
  - Built by `utils/prepare_umi_matrix.py` from `star/samples/*.mtx`
    (rows = genes in `features.tsv` order, cols = `Sample_ID`s from the
    manifest).  `Gene_id` is col 1 (unique per row), `gene_symbol` is col 2
    (annotation).  No folding, no value change — a 1:1 value-verbatim
    projection of the matrix rows.
- **Empty wells are the default in lib-prep workbook mode** (no flag, no
  per-run choice) — details in the **Workbook mode** section.  In
  `--barcode-tsv` mode the manifest is samples-only (no plate, no empty-well
  concept).
- **Per-sample demux runs over samples only** — `empty_well_*` rows are
  filtered out of the demux maps before picard/fqtk run.
- **Two-layer logging** — the run dir has:
  - `pipeline.log` — supervisor log only: the pipeline's own messages
    (RUN/SKIP/done, die/warn) + a pointer to each step's log.  No raw tool
    output.
  - `logs/<step>.log` — one file per pipeline step; the step's tool stdout
    + stderr land here (fastqc, samtools, picard, fqtk, python helpers).
    Step names are sanitised (non-alphanumeric → `_`).
  - STAR's own `star/STAR.log` is unchanged — STAR runtime output + the full
    invocation is written there (the invocation is also echoed as the first
    line of that file for provenance).
  - `run_or_skip` and the pipeline still emit `RUN <step>: <cmd>` to stdout
    (and therefore to `pipeline.log`) so the exact command line of every step
    is recoverable from the supervisor log.
- **Env-first tool resolution** (brbseq env located automatically, no
  `conda activate` needed) and a hard fail on any missing required tool.
- **Pool label inferred from the FASTQ names** (no `--pool` flag) and the
  Mercurius **lib-prep workbook mode** (`--libprep-xlsx`).

### v1.0.2 — resume + robustness fixes

- **`--resume-dir`** — continue an earlier run that didn't finish, in its own
  folder.  Steps with complete output are skipped; only the incomplete step(s)
  re-run, using the sample→barcode list already stored in that folder.  See
  the **Resuming an incomplete run** section.
- **Fixed the FASTQ demux step** so an interrupted/failed attempt is
  detected and re-run instead of being wrongly skipped.
- **Sentinel marker logic** (`run_or_skip` "sentinel" mode) — a step only
  records its done-marker once it truly finished, and a failed optional step
  (e.g. FastQC) no longer leaves a stale marker; it correctly re-runs on the
  next pass.
- **Run-folder index is now `README.html`** — the auto-generated orientation
  page written at the end of each run is HTML (open in a browser) rather than
  plain-text Markdown.

### v1.1.0 — FastReadCounter read-distribution report (CLC)

- **New mandatory `FRC` step** — runs **FastReadCounter** over the STAR
  BAM (already on disk; no second alignment) and the 96-barcode whitelist,
  then derives the legacy-style report:
  - `files_for_report_and_CLC/<pool>_report_summary.csv` — the whole-plate
    read-distribution grid (one row per well A01…H12 + `Unknown`):
    `depth_mio, mapped_reads_mio, no_feature_mio, not_aligned_mio`.
  - `files_for_report_and_CLC/<pool>_all_stats.csv` — the same per-sample
    numbers as raw integer counts (all FRC buckets).
- **Port of the legacy pipeline's FRC + `report_stats` chain**
  (FastReadCounter → rename barcode columns to Sample_ID → two awk passes),
  reusing the pipeline's own `sample_barcode_map.tsv` (all wells) and the
  lib-prep provenance for the well column — output is directly comparable
  with historical runs (validated against AVT0615: identical structure,
  wells, ordering, `depth`/`not_aligned`; `mapped`/`no_feature` only shift
  with a different genome annotation, as expected).
- **No new CLI flags, no new Conda dependency.**  The FastReadCounter jar is
  read from a single hardcoded path at the top of the FRC section
  (`FRC_JAR=/home/gecf/softwares/FastReadCounter/FastReadCounter-1.2.0.jar`
  on the production server); the gene annotation is taken automatically from
  the STAR index (`sjdbGTFfile` in `genomeParameters.txt`).  Java ships with
  the `brbseq` Conda env (OpenJDK) — nothing to install.

## Quick start

```bash
# no conda activate needed — the brb/brbseq env is located automatically

# Barcode TSV manifest
./brbseq_pipeline.sh \
    --species mouse \
    --r1  /data/POOL/S1_R1_001.fastq.gz \
    --r2  /data/POOL/S1_R2_001.fastq.gz \
    --barcode-tsv /data/POOL/sample_barcode.tsv \
    --run-name    Alireza \
    --out-dir /data/BRBseq_output

# Or use the Mercurius lib-prep workbook (barcodes read from it;
# plate inferred from the 'Kit ref' row unless --plate is given)
./brbseq_pipeline.sh \
    --species human \
    --r1  /data/POOL/S1_R1_001.fastq.gz \
    --r2  /data/POOL/S1_R2_001.fastq.gz \
    --libprep-xlsx "Mercurius HS BRB-seq RNA lib prep GECF.xlsx" \
    --plate /data/plates/Lib_01.tsv \
    --run-name Daniel

# Dry run (no side effects) — emits the planned commands + exits 0
./brbseq_pipeline.sh ... --dry-run
```

## Resuming an incomplete run

Every run writes into its own freshly-named folder
(`OUT_DIR/BRBseq_<run-name>_<YYYYMMDD_HHMMSS>/`).  If a run stops before
it finishes — a step failed, the machine ran out of resources, or someone
stopped it — the folder it started is not lost and you do *not* have to
start a brand new run from scratch: re-run with `--resume-dir` pointing at
that folder:

```bash
./brbseq_pipeline.sh \
    --species mouse \
    --r1  /data/POOL/S1_R1_001.fastq.gz \
    --r2  /data/POOL/S1_R2_001.fastq.gz \
    --resume-dir /data/BRBseq_output/BRBseq_Alireza_20261020_134512
```

- No new timestamped folder is created; the same folder the original run
  used is continued in place.
- Each step checks whether its own output files already exist before doing
  anything, so steps the original run finished are simply skipped — you
  don't need to know which step failed or where to look for it.
- The sample→barcode list is taken from the one already stored inside that
  same folder (`barcodes/`), so there is no `--barcode-tsv` /
  `--libprep-xlsx` / `--plate` on this command (passing one is an error —
  if you actually need a *different* sample→barcode list, that's a fresh
  run: drop `--resume-dir` and use `--run-name` like normal).
- An older practice of manually deleting a leftover `_done` marker file
  (e.g. `demux_fastq/_done`) before re-running is no longer necessary —
  the script now detects an incomplete step's output on its own.

## Status

`SCRIPT_VERSION="1.1.0"` (in `brbseq_pipeline.sh`) — the current, working
version of the pipeline.  The first real-pool deployment (against an
AVT0608/AVT0615-style pool) is built on top of this baseline.  See **Changelog**
above for what each release shipped; reference changes by feature or date (not
by the old intermediate version numbers).

## Reference (authoritative)

- Alithea **MERCURIUS 3' RNA-seq Data Pre-Processing User Guide** (July 2026)
  — §1.4 STAR request, §1.6 duplication step, §1.6.x matrix deliverable.
- `brbseq_pipeline.sh --help` — full flag reference + worked examples.
- The **brb-seq** Hermes skill (research category) — full design notes,
  provenance, and historical decisions.
