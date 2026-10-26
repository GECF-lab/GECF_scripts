#!/usr/bin/env bash
#=============================================================================
# brbseq_pipeline.sh — MERCURIUS DRUG-seq (BRB-seq) single-pool pipeline
#
# MERCURIUS 3' RNA-seq (BRB-seq/DRUG-seq) Data Pre-Processing pipeline
# for ONE pool/run, in step order:
#
#   STAR genome indexing            (one-time, NOT a pipeline step)
#   Creating barcode files          [MANDATORY ]  prepare_barcodes.sh
#   FASTQ quality (fastqc)          [ON by default]  off: --skip-fastq
#   STARsolo alignment + counts     [MANDATORY ]  STAR
#   Read distribution (RSeQC)       [ON by default]  off: --skip-rseqc
#   Duplication statistics          [MANDATORY ]  get_duplication_statistics.py
#   BAM demultiplex (picard)        [ON by default]  off: --skip-dmux-bam
#   FASTQ demultiplex (fqtk)        [ON by default]  off: --skip-dmux-fq
#   Merge multi-lane FASTQ          [note       ]  pre-merge or comma-list
#
# Logging (automatic, no flag)
#   * pipeline.log — supervisor log: every pipeline message (RUN/SKIP/done,
#     die/warn) + a pointer to each step's own log.  No raw tool output.
#   * logs/<step>.log — one file per pipeline step; tool stdout + stderr go
#     here (fastqc, samtools, picard, fqtk, python helpers).  Step name is
#     sanitised (non-alphanumeric chars → "_").
#   * star/STAR.log — STAR's own runtime output; the full STAR invocation is
#     also written as its first line for provenance.  The exact command line
#     of every other step is logged the same way (run_or_skip prints
#     "RUN  <step>: <cmd>").
#
# Note on the STAR invocation
#   The STAR request below asks for BOTH dedup methods (NoDedup and
#   1MM_Directional) in one alignment pass, so both result matrices are
#   produced at once. This is required by the mandatory duplication step,
#   which consumes both umiDedup-NoDedup.mtx (before dedup) and
#   umiDedup-1MM_Directional.mtx (after dedup). It matches the Alithea
#   production run:
#     * --soloUMIdedup NoDedup 1MM_Directional
#
# Usage notes
#   * All inputs are supplied via flags; no paths are hard-coded.
#   * Portable across hosts: activate your environment before running
#     (e.g. 'conda activate brbseq') so STAR, samtools, RSeQC, fastqc,
#     picard, fqtk and python+openpyxl are on PATH.  If the env does not
#     exist yet on the target machine, create it first from the env file
#     shipped next to this script:
#     conda env create -f brbseq-env.yaml            (creates env 'brbseq')
#     conda activate brbseq
#     (helper scripts are resolved from ./utils, shipped with this directory).
#   * Re-runnable: a step is skipped when its output already exists
#     (re-run the same command to continue; finished steps are skipped).
#   * The optional steps (FASTQ QC, RSeQC, BAM demux, FASTQ demux) are ON by
#     default (disable each with its --skip-* flag); a missing tool causes
#     that step to log a warning and skip rather than abort the run.
#=============================================================================
set -euo pipefail

#--- constants ---------------------------------------------------------------
SCRIPT_VERSION="1.1.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="${UTILS_DIR:-${SCRIPT_DIR}/utils}"
PREPARE_BARCODES="${UTILS_DIR}/prepare_barcodes.sh"
DEMUX_BAM="${UTILS_DIR}/demultiplex_bam.sh"
DUP_STATS="${UTILS_DIR}/get_duplication_statistics.py"
PARSE_LIBPREP="${UTILS_DIR}/parse_libprep_xlsx.py"
RENAME_MATRIX="${UTILS_DIR}/rename_matrix.py"
UMI_MATRIX_PREP="${UTILS_DIR}/prepare_umi_matrix.py"

#--- defaults (all overridable via flags) -----------------------------------
RUN_NAME="${RUN_NAME:-}"         # run label for the output leaf: BRBseq_<run-name>_<timestamp>
AVT_ID="${AVT:-}"                # DEPRECATED: no longer a pipeline parameter; value (if set) is used only as a fallback run name
POOL=""                           # pool label — inferred from the FASTQ names
THREADS="${THREADS:-8}"
R1_BARCODE=""                    # R1: 14 bp CB + 14 bp UMI  (barcode read)
R2_GENOMIC=""                    # R2: genomic/cDNA read
BARCODE_TSV=""                   # manifest: Sample_ID<TAB>Barcode
# Lib-prep workbook mode (alternative to --barcode-tsv):
LIBPREP_XLSX=""                  # Mercurius lib-prep .xlsx (sample code + Oligo-dT primer well)
BARCODES_PLATE=""                # well-matched barcode plate (well + 14-nt CB per line); empty = infer from the workbook's 'Kit ref' row
SPECIES=""                       # human | mouse | <custom name>
STAR_INDEX=""                    # STAR --genomeDir; resolved from --species if unset

#--- built-in references -----------------------------------------------------
# One line per built-in species: <species>|<STAR index dir>.
# Server-independent defaults for the production server; override per call
# with --star-genome or by exporting REF_HUMAN / REF_MOUSE.
REFS=(
  "human|${REF_HUMAN:-/home/gecf/references/STAR_references/BRB-seq/current_human}"
  "mouse|${REF_MOUSE:-/home/gecf/references/STAR_references/BRB-seq/current_mouse}"
)
BED_FILE=""                      # RSeQC BED12
OUT_DIR="${OUT_DIR:-/mnt/data/BRB-seq_output}"
MANIFEST_HEADER=1                # resolved to 1 (header row) or 0 (none) from auto-detection
# Optional analyses are ON by default: fastqc, read_distribution
# (RSeQC), BAM demux (picard), FASTQ demux (fqtk).  Disable any of them
# with --skip-fastq / --skip-rseqc / --skip-dmux-bam / --skip-dmux-fq.
# Disable any of them with --skip-fastq / --skip-rseqc / --skip-dmux-bam /
# --skip-dmux-fq.
RUN_OPT_FASTQ="${RUN_OPT_FASTQ:-1}"
RUN_OPT_RSEQC="${RUN_OPT_RSEQC:-1}"
RUN_OPT_DEMUX_BAM="${RUN_OPT_DEMUX_BAM:-1}"
RUN_OPT_DEMUX_FQ="${RUN_OPT_DEMUX_FQ:-1}"
# did the user pass the --skip-* opt-out for each?
S_FASTQ=0; S_RSEQC=0; S_DMBAM=0; S_DMFQ=0
# Alithea BRB-seq production SAM tag set (fixed, not user-configurable).
SAM_ATTRS="NH HI nM AS CR UR CB UB GX GN sS sQ sM"
RESUME_DIR="${RESUME_DIR:-}"     # --resume-dir: target an existing run dir in place
DRY_RUN="${DRY_RUN:-0}"
PYTHON="${PYTHON:-python3}"

#--- log helpers -------------------------------------------------------------
_ts()  { date +%H:%M:%S; }
log()  { printf '%s [INFO ] %s\n' "$(_ts)" "$*"; }
warn() { printf '%s [WARN ] %s\n' "$(_ts)" "$*" >&2; }
die()  { printf '%s [ERROR] %s\n' "$(_ts)" "$*" >&2; exit 1; }

usage() {
  printf 'brbseq_pipeline.sh — MERCURIUS DRUG-seq (BRB-seq) single-pool pipeline (v%s).\n\n' "$SCRIPT_VERSION"
  cat <<'EOF'
MERCURIUS 3' RNA-seq – Data Pre-Processing (BRB-seq / DRUG-seq) pipeline.

REQUIRED FLAGS
  --species <sp>        human | mouse | <custom name>
  --r1 <fq.gz>          R1 barcode read (14 bp CB + 14 bp UMI)
  --r2 <fq.gz>          R2 genomic / cDNA read
  --run-name <name>     run label for the output leaf
                        BRBseq_<run-name>_<timestamp>; in lib-prep mode
                        taken from the workbook 'Researcher' row (first
                        name) unless given here. Required in --barcode-tsv
                        mode (no workbook to read).

  MANIFEST SOURCE — exactly one of these two:
    (a)  --barcode-tsv <tsv>
    (b)  --libprep-xlsx <xlsx>   (+ optional --plate <file>)
         Mercurius lib-prep workbook (sample code + Oligo-dT primer well). The
         barcode plate is joined by well. If --plate is omitted, it is
         inferred from the workbook's 'Kit ref' row (V5A/V5D ->
         barcodes_96_<family>_*.txt in the plate directory
         /home/gecf/NGS_pipelines/BRB-seq/BRB_plates, override: $BRB_PLATE_DIR);
         if given, it is authoritative — a family disagreement with the
         workbook's 'Kit ref' only produces a loud warning (the named file is
         used as-is).

OPTIONAL FLAGS
  --star-genome <dir>   STAR genomeDir (built with build_star_index.sh);
                        required for custom --species, overrides built-in
  --bed-file <bed12>    RSeQC BED12 gene model (override). If omitted, bakes
                        one at run time from the GTF shipped inside the STAR
                        index dir (build_star_index.sh moves it there)
  --out-dir <dir>       root output dir (default /mnt/data/BRB-seq_output)
  --threads <n>         CPU threads (default 8); number of parallel threads
                        used by STAR (alignment) and fqtk (per-sample FASTQ
                        demux)
  --resume-dir <dir>    continue an earlier run that didn't finish, in its
                        own folder. Steps that already have complete output
                        are skipped; only the incomplete step(s) re-run,
                        using the sample→barcode list already stored in that
                        folder (barcodes/) — so pass --r1 --r2 --species
                        again, not a new --libprep-xlsx/--barcode-tsv/--plate.
                        (See the "Resuming an incomplete run" section of
                        README.html.)
  --skip-fastq          skip fastqc  (on by default)
  --skip-rseqc          skip read_distribution  (on by default)
  --skip-dmux-bam       skip picard BAM demux  (on by default)
  --skip-dmux-fq        skip fqtk FASTQ demux  (on by default)
  --dry-run             emit planned commands and exit (no side effects)

  -h, --help            this help
  --version             print version and exit

EXAMPLES
  # TSV manifest:
  ./brbseq_pipeline.sh --species mouse \
      --r1 /data/POOL/S1_R1_001.fastq.gz \
      --r2 /data/POOL/S1_R2_001.fastq.gz \
      --barcode-tsv /data/POOL/sample_barcode.tsv \
      --run-name Daniel --threads 16 --out-dir /mnt/data/BRB-seq_output

  # Lib-prep workbook instead of a TSV (plate inferred from the workbook):
  ./brbseq_pipeline.sh --species mouse \
      --r1 /data/POOL/S1_R1_001.fastq.gz \
      --r2 /data/POOL/S1_R2_001.fastq.gz \
      --libprep-xlsx "Mercurius HS BRB-seq RNA lib prep GECF.xlsx"

  # Custom species: build a STAR index once, then pass it:
  ./build_star_index.sh --fasta genome.fa --gtf genes.gtf --outdir star_index --sjdbOverhang 89
  ./brbseq_pipeline.sh --species mouse-custom --star-genome star_index ...

  # Dry run (no side effects):
  ./brbseq_pipeline.sh ... --dry-run

  # See the README for environment setup (conda), full flags, and output layout.
EOF
}

#--- argument parsing --------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r1)   R1_BARCODE="${2:?missing --r1 value}";   shift 2;;
    --r2)   R2_GENOMIC="${2:?missing --r2 value}";   shift 2;;
    --barcode-tsv)  BARCODE_TSV="${2:?missing --barcode-tsv value}"; shift 2;;
    --libprep-xlsx) LIBPREP_XLSX="${2:?missing --libprep-xlsx value}"; shift 2;;
    --plate) BARCODES_PLATE="${2:?missing --plate value}"; shift 2;;
    --species)      SPECIES="${2:?missing --species value}";          shift 2;;
    --star-genome)   STAR_INDEX="${2:?missing --star-genome value}";   shift 2;;
    --run-name)   RUN_NAME="${2:?missing --run-name value}";           shift 2;;
    --bed-file)     BED_FILE="${2:?missing --bed-file value}";      shift 2;;
    --out-dir)      OUT_DIR="${2:?missing --out-dir value}";        shift 2;;
    --threads)      THREADS="${2:?missing --threads value}";        shift 2;;
    --skip-dmux-fq) S_DMFQ=1;                                         shift;;
    --skip-fastq)   S_FASTQ=1;                                        shift;;
    --skip-rseqc)   S_RSEQC=1;                                        shift;;
    --skip-dmux-bam)S_DMBAM=1;                                        shift;;
    --resume-dir)   RESUME_DIR="${2:?missing --resume-dir value}";    shift 2;;
    --dry-run)      DRY_RUN=1;                                        shift;;
    -h|--help)      usage; exit 0;;
    --version)      echo "brbseq_pipeline.sh ${SCRIPT_VERSION}"; exit 0;;
    -*)             die "I don't recognise the option '$1'. Please check the spelling (options start with '--' or '-'). The full list of options is printed by:  brbseq_pipeline.sh --help";;
    *)              die "I don't recognise '$1' as a valid option or value. The full list of options is printed by:  brbseq_pipeline.sh --help";;
  esac
done

#--- apply --skip-* opt-outs (steps run by default otherwise) ---------------
[[ "$S_FASTQ"  -eq 1 ]] && RUN_OPT_FASTQ=0
[[ "$S_RSEQC"  -eq 1 ]] && RUN_OPT_RSEQC=0
[[ "$S_DMBAM"  -eq 1 ]] && RUN_OPT_DEMUX_BAM=0
[[ "$S_DMFQ"   -eq 1 ]] && RUN_OPT_DEMUX_FQ=0

#--- validate required inputs ------------------------------------------------
# Exactly one manifest source must be provided:
#   --barcode-tsv  (2-col Sample_ID/Barcode TSV),  OR
#   --libprep-xlsx  (+ optional --plate).
# The plate is optional in workbook mode: parse_libprep_xlsx.py infers it from
# the workbook's 'Kit ref' row (barcodes_96_<family>_*.txt next to the
# workbook); an explicitly given plate is authoritative (a family mismatch with
# the Kit ref is a loud warning, not an abort).
if [[ -n "$BARCODE_TSV" && -n "$LIBPREP_XLSX" ]]; then
  die "You specified both --barcode-tsv and --libprep-xlsx, but only one of them can be used (both define the sample ↔ barcode mapping).  Remove one of the two flags and re-run. (With --resume-dir the manifest inside the run dir is used — no new manifest source is needed.)"
fi
LIBPREP_MODE=0
if [[ -n "$LIBPREP_XLSX" ]]; then
  LIBPREP_MODE=1
fi
# --resume-dir: the run was already STARTED — the sample→barcode list must
# exist inside that run dir (barcodes/manifest_from_libprep.tsv for workbook
# mode, or the user TSV for --barcode-tsv mode).  No new manifest source is
# generated or required for this re-run.
if [[ -n "$RESUME_DIR" ]]; then
  LIBPREP_MODE=1
  if [[ -n "$LIBPREP_XLSX" || -n "$BARCODE_TSV" || -n "$BARCODES_PLATE" ]]; then
    die "--resume-dir uses the sample→barcode list already stored inside ${RESUME_DIR} — do not also pass --libprep-xlsx / --barcode-tsv / --plate. (If you WANT to regenerate the list from a new workbook that is a fresh run: drop --resume-dir and give --run-name for a new folder.)"
  fi
  [[ -d "$RESUME_DIR" ]] || die "--resume-dir folder not found: ${RESUME_DIR} (this should be an existing pipeline output folder named BRBseq_<name>_<timestamp>). Use --help for the full list of options."
  _rd_manifest=""
  for _mf in "${RESUME_DIR%/}/barcodes/manifest_from_libprep.tsv" "${RESUME_DIR%/}/barcodes/barcode_manifest.tsv"; do
    if [[ -s "$_mf" ]]; then _rd_manifest="$_mf"; break; fi
  done
  [[ -n "$_rd_manifest" ]] || die "--resume-dir needs the sample→barcode list (${RESUME_DIR}/barcodes/manifest_from_libprep.tsv) which the original run produced, but it is missing or empty there. Either that run never got as far as the barcode step, or the folder was pruned. Check the folder's pipeline.log to see how far it got, then fix that step first (see --help) and re-run into the same folder."
  RESUME_MANIFEST="$_rd_manifest"
  [[ -s "${RESUME_DIR%/}/barcodes/barcode_whitelist.txt" ]] || die "--resume-dir needs the barcode whitelist (${RESUME_DIR}/barcodes/barcode_whitelist.txt) which the original run produced, but it is missing or empty — the barcode-prep step did not complete in that run. See its entry in that folder's pipeline.log, fix that step first (see --help), then re-run into the same folder."
fi
MISSING=()
[[ -n "$R1_BARCODE"  ]] || MISSING+=( "--r1" )
[[ -n "$R2_GENOMIC"  ]] || MISSING+=( "--r2" )
# manifest source: exactly one of (a) --barcode-tsv  or  (b) --libprep-xlsx.
if [[ -z "$LIBPREP_XLSX" && -n "$BARCODES_PLATE" && -z "$RESUME_DIR" ]]; then
  die "--plate only makes sense together with the lib-prep workbook, but --libprep-xlsx was not given. Either add --libprep-xlsx <workbook.xlsx>, or remove --plate and let the workbook's 'Kit ref' row choose the plate automatically."
elif [[ "$LIBPREP_MODE" -ne 1 && -z "$BARCODE_TSV" ]]; then
  MISSING+=( "manifest source: --barcode-tsv <tsv>   or   --libprep-xlsx <xlsx>  (+ optional --plate <file>)" )
fi
[[ -n "$SPECIES"     ]] || MISSING+=( "--species" )
[[ -n "$RUN_NAME" || -n "$AVT_ID" || "$LIBPREP_MODE" -eq 1 ]]   || MISSING+=( "--run-name <name>" )
[[ "${#MISSING[@]}" -eq 0 ]] || die "The command is missing required option(s): ${MISSING[*]}. Add the listed option(s) to the command and re-run. (Full list of options: brbseq_pipeline.sh --help)"

#--- pool label ---------------------------------------------------------------
# The pool label is INFERRED from the FASTQ names (no user parameter).
# Rule: take the basename, cut at the lane marker "_S<n>" — everything
# before "_S" is kept verbatim.
#   pool_poolA_S1_R1_001.fastq.gz  -> pool_poolA
#   pool_poolB_S2_R2_001.fastq.gz  -> pool_poolB
#   sample_S1_R1_001.fastq.gz      -> sample
# R1 and R2 must agree; a disagreement is a WARNING (R1 wins).
infer_pool() {
  local base; base="$(basename "$1")"
  base="${base%_S[0-9]*}"     # cut at the lane marker "_S<n>"
  if [[ -n "$base" && ! "$base" =~ [./] && "$base" != *S[0-9]* ]]; then
    printf '%s' "$base"
  else
    printf ''
  fi
}
POOL_R1=""; POOL_R2=""
if [[ -n "$R1_BARCODE" ]]; then POOL_R1="$(infer_pool "$R1_BARCODE")"; fi
if [[ -n "$R2_GENOMIC" ]]; then POOL_R2="$(infer_pool "$R2_GENOMIC")"; fi
if [[ -n "$POOL_R1" && -n "$POOL_R2" && "$POOL_R1" == "$POOL_R2" ]]; then
  log "pool       = ${POOL_R1}  (inferred from FASTQ names)"
  POOL="$POOL_R1"
elif [[ -n "$POOL_R1" && -n "$POOL_R2" ]]; then
  warn "CAUTION: the two FASTQ files you gave disagree about which pooled experiment this is. The R1 (barcode) file suggests the pool is named  '$POOL_R1'  and the R2 (genomic) file suggests  '$POOL_R2'. This is only used to LABEL the output folders (e.g.  ${POOL}_aligned.bam), it does not change the biological analysis. The pipeline is using '$POOL_R1' (taken from the R1/barcode file). If the two really are from the same pooled run, just rename the FASTQ files so both spell the pool name the same way, and it will stop warning."
  POOL="$POOL_R1"
elif [[ -n "$POOL_R1" ]]; then
  warn "CAUTION: the file name '$(basename "$R2_GENOMIC")' (R2/genomic) does not let the pipeline read a 'pool' name out of it — it expected something like <pool>_<label>_S<lane>_...fastq.gz. The pool is only used to label output folders, so the pipeline is using '$POOL_R1' (read from the R1/barcode file) instead. If that label looks wrong, rename this FASTQ file so it follows the expected pattern, and it will stop warning."
  POOL="$POOL_R1"
elif [[ -n "$POOL_R2" ]]; then
  warn "CAUTION: the file name '$(basename "$R1_BARCODE")' (R1/barcode) does not let the pipeline read a 'pool' name out of it — it expected something like <pool>_<label>_S<lane>_...fastq.gz. The pool is only used to label output folders, so the pipeline is using '$POOL_R2' (read from the R2/genomic file) instead. If that label looks wrong, rename this FASTQ file so it follows the expected pattern, and it will stop warning."
  POOL="$POOL_R2"
else
  die "The pipeline reads the 'pool' name from the FASTQ file names, but neither of your two files uses a recognisable name. A valid name has the form <pool>_<label>_S<lane>_R<read>_...fastq.gz — for example pool_poolA_S1_R1_001.fastq.gz (everything before _S1 is the pool name). Please rename the FASTQ files so they follow that pattern, then re-run."
fi

#--- reference resolution ------------------------------------------------------
# --species selects the STAR genome:
#   human|mouse : the index is taken from the built-in table above; an
#                 explicit --star-genome wins over it.
#   <other>     : custom species — a prebuilt index must be passed via
#                 --star-genome (build it once with build_star_index.sh).
if [[ -z "$STAR_INDEX" ]]; then
  for line in "${REFS[@]}"; do
    IFS='|' read -r g idx _ <<<"$line"
    if [[ "$g" == "$SPECIES" ]]; then STAR_INDEX="$idx"; break; fi
  done
fi
# Built-in species point at a "current" wrapper dir that holds exactly one
# inner folder — the STAR genomeDir (plus the .fa/.gtf it was built from).
# Genome rebuilds swap that inner folder while the wrapper path stays valid.
# Resolve: if $STAR_INDEX itself is not a finished index, look inside for the
# single subfolder that is.
if [[ -n "$STAR_INDEX" && ! -f "$STAR_INDEX/Genome" && ! -f "$STAR_INDEX/SA" ]]; then
  inner=""
  for d in "$STAR_INDEX"/*/; do
    [[ -f "${d}Genome" && -f "${d}SA" ]] && inner="${d%/}"
  done
  [[ -n "$inner" ]] || die "The folder you gave with --star-genome ($STAR_INDEX) does not look like a STAR genome index — a STAR index is a folder that contains one or more subfolders, each holding two small files called 'Genome' and 'SA' (built by build_star_index.sh). This folder has no such subfolder, so the pipeline can't align the reads. Either check that the path you passed is really the index folder (and not, e.g., the folder one level above it), or build a fresh index with:  ./build_star_index.sh --fasta <genome.fa> --gtf <annotation.gtf> --outdir <index-dir> --sjdbOverhang 89  and re-run the pipeline with that <index-dir> as --star-genome."
  STAR_INDEX="$inner"
fi
if [[ -z "$STAR_INDEX" ]]; then
  die "--species '$SPECIES' has no built-in reference and no --star-genome was given. Either choose a supported species (human or mouse), or point the pipeline at a pre-built STAR index with --star-genome <index_dir>. If no such index exists yet, build one once with: ./build_star_index.sh --fasta <genome.fa> --gtf <annotation.gtf> --outdir <index-dir> --sjdbOverhang 89, then re-run the pipeline with that <index-dir>."
fi
#--- RSeQC reference resolution (bake BED12 at run time if no --bed-file) ----
# RSeQC read_distribution.py needs a 12-column gene-model BED (BED12). If the
# user did not pass --bed-file, we bake one at run time from the GTF that
# build_star_index.sh moved into the STAR index dir — that GTF is the exact
# gene model the index was built from, so provenance is self-contained in the
# index folder. The baked file lands in <OUT>/rseqc/<gtf>.bed12.
#
# A clean index dir contains exactly ONE GTF (plain build: the original
# <annotation>.gtf; ERCC build: annotation_with_ERCC.gtf).  We do NOT guess
# between several: 0 or >=2 GTFs in the index dir is a reference-packaging
# problem and is reported HERE, at validation time (start of run), not at
# RSeQC:
#   1 GTF          -> RSEQC_GTF set, BED12 baked at run time (unless --bed-file given)
#   >1 GTFs        -> die (explicit list of all of them) — or warn and
#                     continue if --bed-file is given (baking never happens)
#   0 GTFs         -> read_distribution auto-skips (warn), rest of pipeline proceeds
RSEQC_GTF=""
if [[ -d "$STAR_INDEX" ]]; then
  RSEQC_GTF_CANDIDATES=()
  for ext in .gtf .gtf.gz; do
    for f in "$STAR_INDEX"/*"$ext"; do
      [[ -f "$f" ]] && RSEQC_GTF_CANDIDATES+=("$f")
    done
  done
  # Prefer the annotation this index was actually built from: STAR records
  # the exact path in genomeParameters.txt (sjdbGTFfile). This keeps the
  # right choice even when the index dir legitimately holds several
  # GTFs (e.g. the original genome GTF next to the merged ..._with_ERCC.gtf).
  GP_GTF="$(awk '$1=="sjdbGTFfile"{print $2; exit}' "$STAR_INDEX/genomeParameters.txt" 2>/dev/null || true)"
  if [[ -n "$GP_GTF" && ${#RSEQC_GTF_CANDIDATES[@]} -gt 1 ]]; then
    for f in "${RSEQC_GTF_CANDIDATES[@]}"; do
      if [[ "$f" == "$GP_GTF" || "$(basename "$f")" == "$(basename "$GP_GTF")" ]]; then
        RSEQC_GTF_CANDIDATES=("$f")
        break
      fi
    done
  fi
  n=${#RSEQC_GTF_CANDIDATES[@]}
  if [[ "$n" -eq 1 ]]; then
    RSEQC_GTF="${RSEQC_GTF_CANDIDATES[0]}"
  elif [[ "$n" -gt 1 ]]; then
    # Several annotations in one index dir is a packaging bug (multiple
    # genome builds mixed in the same folder).  If RSeQC is requested we
    # must not guess — die now, at the start of the run, listing every
    # offender.  If RSeQC is not in this run (or --bed-file makes the
    # choice moot) we still say so explicitly, but do not block the run.
    if [[ "$RUN_OPT_RSEQC" -eq 1 && -z "$BED_FILE" ]]; then
      # read_distribution is ON by default, so a broken gene model must not
      # abort the whole run: skip the read-distribution step only and continue.
      warn "WARNING (read-distribution step): the STAR index folder $STAR_INDEX contains more than one GTF gene annotation ($n files found), and the pipeline does not know which one is the correct gene model for this genome, so it can't build the RSeQC BED12 file automatically this run. The read-distribution report will NOT be produced; the rest of the pipeline continues normally. To turn the report back on: pass the correct gene model already converted to BED12 format with --bed-file <file.bed12>, or remove the extra .gtf files from that folder so only the correct one remains."
      RUN_OPT_RSEQC=0
    else
      local_list="$(printf '%s, ' "${RSEQC_GTF_CANDIDATES[@]}")"
      warn "WARNING (read-distribution step): the STAR index folder $STAR_INDEX contains more than one GTF gene annotation (${local_list%, }). The pipeline does not know which one of these is the correct gene model for this genome, so it will not try to build the RSeQC BED12 file automatically this run — the rest of the pipeline continues normally, but the read-distribution report will not be produced. To turn on the read-distribution report, supply the correct gene model directly with  --bed-file <file.bed12>."
    fi
  fi
fi
if [[ -z "$BED_FILE" && -n "$RSEQC_GTF" ]]; then
  log "RSeQC: no --bed-file given; will bake BED12 from $RSEQC_GTF (the single GTF in the STAR index dir — the gene model the index was built from)"
elif [[ -z "$BED_FILE" && -z "$RSEQC_GTF" ]]; then
  # Nothing to bake from; the check below skips read_distribution (warn).
  :
fi
log "reference = ${SPECIES}  ->  STAR index $STAR_INDEX"
#--- (end RSeQC reference resolution) -------------------------------------------------

#--- validate file/directory inputs ------------------------------------------
[[ -f "$R1_BARCODE"  ]] || die "The R1 (barcode-channel) FASTQ file was not found: $R1_BARCODE.  Check that the path in --r1 is correct and the file is still there, then re-run."
[[ -f "$R2_GENOMIC"  ]] || die "The R2 (genomic-channel) FASTQ file was not found: $R2_GENOMIC.  Check that the path in --r2 is correct and the file is still there, then re-run."
# In lib-prep mode the manifest is generated later; only validate the workbook
# here.  The plate is optional: when omitted, the
# parser infers it from the workbook's 'Kit ref' row.
if [[ "$LIBPREP_MODE" -eq 1 && -z "$RESUME_DIR" ]]; then
  [[ -f "$LIBPREP_XLSX" ]] || die "The lib-prep workbook (--libprep-xlsx) file could not be found: $LIBPREP_XLSX. Check the path and re-run."
  [[ -z "$BARCODES_PLATE" || -f "$BARCODES_PLATE" ]] || die "The barcode plate file (--plate) could not be found: $BARCODES_PLATE. Check the path and re-run."
elif [[ "$LIBPREP_MODE" -eq 1 && -n "$RESUME_DIR" ]]; then
  :   # resume mode: the manifest lives inside the run dir (validated above)
else
  [[ -f "$BARCODE_TSV" ]] || die "The sample manifest (--barcode-tsv) file could not be found: $BARCODE_TSV. It should be a plain-text file with exactly two tab-separated columns: sample_id  [TAB]  barcode (one sample per line). Check the path and re-run."
fi
[[ -d "$STAR_INDEX"  ]] || die "The STAR genome index directory could not be found: $STAR_INDEX. Check the --star-genome value (or the built-in path for your --species) and re-run. If the index was deleted or moved, you will need to rebuild it with build_star_index.sh."
if [[ "$RUN_OPT_RSEQC" -eq 1 ]]; then
  if [[ -z "$BED_FILE" && -z "$RSEQC_GTF" ]]; then
    # read_distribution is ON by default: a missing gene model skips only that
    # step (the rest of the pipeline — alignment, counts, demux — still runs).
    warn "WARNING (read-distribution step): no gene-model file is available for it — neither --bed-file nor a GTF inside the STAR index directory was found, so the pipeline cannot build the BED12 it needs. The read-distribution report will NOT be produced; the rest of the pipeline continues normally. To turn the report back on: give the file directly with --bed-file <file.bed12>, or make sure the GTF annotation is present inside the STAR index directory (build_star_index.sh normally copies it there)."
    RUN_OPT_RSEQC=0
  fi
  if [[ -n "$BED_FILE" && ! -f "$BED_FILE" ]]; then
    die "The BED12 gene-model file passed with --bed-file was not found: $BED_FILE. Check the path and re-run. If you do not have a BED12 file, leave --bed-file out and instead make sure the original GTF annotation is inside the STAR index directory, so the pipeline can convert it to BED12 automatically."
  fi
fi

#--- validate shipped helpers ------------------------------------------------
for req in "$PREPARE_BARCODES" "$DEMUX_BAM" "$DUP_STATS" "$RENAME_MATRIX" "$UMI_MATRIX_PREP"; do
  [[ -f "$req" ]] || die "The pipeline's own supporting helper file is missing from this installation: $req — these small helper scripts should live in a 'utils' folder right next to brbseq_pipeline.sh. This means the pipeline installation itself is incomplete (for example the 'utils' folder was not copied when the pipeline was being set up on this machine). Fix: make sure the whole BRBseq folder (including its 'utils' subfolder) is intact in one place, point the pipeline at it again, and re-run. If you just copied files between machines, this is usually where a step got dropped."
done
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  [[ -f "$PARSE_LIBPREP" ]] || die "The pipeline's own small helper program that reads the lib-prep Excel workbook is missing from this installation: $PARSE_LIBPREP — it should sit in a 'utils' folder right next to brbseq_pipeline.sh. This means the pipeline installation itself is incomplete (for example the 'utils' folder was not copied when the pipeline was being set up on this machine). Fix: make sure the whole BRBseq folder (including its 'utils' subfolder) is intact in one place, point the pipeline at it again, and re-run."
fi

#--- run name (resolution, pre-layout) ---------------------------------------
# Priority: --run-name flag  >  AVT env (deprecated)  >  workbook 'Researcher' row
# (lib-prep mode only).  In lib-prep mode we ask the parser module for the
# researcher name now, before the output dir is known, so the leaf can be named;
# the manifest is generated into <run>/barcodes/ later, after the layout.
if [[ -z "$RUN_NAME" && -n "$RESUME_DIR" ]]; then
  _rd_base="$(basename "${RESUME_DIR%/}")"
  [[ "$_rd_base" =~ ^BRBseq_(.+)_.{15}$ ]] || die "--resume-dir points at '${RESUME_DIR}' which is not named like an existing run folder (expected BRBseq_<name>_<timestamp>). Re-check the path, or pass --run-name explicitly."
  RUN_NAME="${BASH_REMATCH[1]}"
  log "run name   = ${RUN_NAME}  (from the --resume-dir folder name; pass --run-name to override)"
fi
if [[ -z "$RUN_NAME" ]]; then
  if [[ -n "$AVT_ID" ]]; then
    RUN_NAME="$AVT_ID"
    log "run name   = ${RUN_NAME}  (from AVT env; use --run-name to override)"
  elif [[ "$LIBPREP_MODE" -eq 1 ]]; then
    RUN_NAME="$("$PYTHON" - "$PARSE_LIBPREP" "$LIBPREP_XLSX" <<'PYEOF'
import sys, importlib.util
spec = importlib.util.spec_from_file_location("_pl", sys.argv[1])
mod  = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
import openpyxl
wb = openpyxl.load_workbook(sys.argv[2], data_only=True)
for ws in wb.worksheets:
    nm = mod.extract_researcher_first_name(ws)
    if nm:
        print(nm); break
PYEOF
)"
    [[ -n "$RUN_NAME" ]] || die "could not infer a run name from the workbook (no 'Researcher' row) — pass --run-name <name>"
    log "run name   = ${RUN_NAME}  (inferred from the workbook 'Researcher' row)"
  else
    die "A run name is required so the output folder can be labelled, but no name was found. Pass it explicitly with --run-name <name> (or use --libprep-xlsx, whose workbook records the researcher's name and supplies it automatically)."
  fi
fi
RUN_NAME="$(tr -dc 'A-Za-z0-9._-' <<<"${RUN_NAME// /_}" )"   # filesystem-safe
[[ -n "$RUN_NAME" ]] || die "The run name you supplied (--run-name) is made up only of characters that are not allowed in file names, so after cleaning it became an empty name and the pipeline has no way to label the output folder. Please use a run name made of letters, numbers, dots, dashes or underscores (spaces will be turned into underscores) — for example 'PHX537_poolA' — and re-run."

#--- output layout -----------------------------------------------------------
if [[ -n "$RESUME_DIR" ]]; then
  # Re-run an EXISTING run folder in place (no new timestamped dir).
  # Every step re-checks its own completion state (see the resume
  # logic per step below), so finished steps are skipped and the one that
  # did not finish (e.g. the failed fastq demultiplex) is re-run.
  RUN_DIR="${RESUME_DIR%/}"
else
  TS="$(date +%Y%m%d_%H%M%S)"
  RUN_DIR="${OUT_DIR%/}/BRBseq_${RUN_NAME}_${TS}"
fi
BARCODES_DIR="${RUN_DIR}/barcodes"
FASTQC_DIR="${RUN_DIR}/fastqc"
STAR_DIR="${RUN_DIR}/star"
RSEQC_DIR="${RUN_DIR}/rseqc"
DMUX_BAM_DIR="${RUN_DIR}/demux_bam"
DMUX_FQ_DIR="${RUN_DIR}/demux_fastq"
RENAME_DIR="${STAR_DIR}/samples"
RENAME_ALL_DIR="${STAR_DIR}/samples_all_wells"
DUP_CSV="${RUN_DIR}/duplication_statistics.csv"
MATRIX_DIR="${RUN_DIR}/matrices"
UMI_MATRIX_CSV="${MATRIX_DIR}/${POOL}_umi_matrix.csv"
READS_MATRIX_CSV="${MATRIX_DIR}/${POOL}_reads_matrix.csv"
# In --barcode-tsv mode the manifest lives where the user put it (no generation);
# in lib-prep mode the manifest is generated into <run>/barcodes/ and survives re-runs.
MANIFEST="$BARCODE_TSV"
if [[ -n "${RESUME_MANIFEST:-}" ]]; then
  # --resume-dir: reuse exactly the sample→barcode list the run was started with.
  MANIFEST="$RESUME_MANIFEST"
  LIBPREP_PROV=""
elif [[ "$LIBPREP_MODE" -eq 1 ]]; then
  MANIFEST="${BARCODES_DIR}/manifest_from_libprep.tsv"
  LIBPREP_PROV="${MANIFEST}_provenance.tsv"
else
  LIBPREP_PROV=""
fi

#--- manifest header decision (auto-detected, no flag) ----------------------
# parse_libprep_xlsx.py (GECF utility parser) always emits a header row on
# line 1, so in lib-prep mode the answer is known without inspecting the file.
# For user-supplied --barcode-tsv, auto-detect from line 1: a valid line-1
# barcode cell (14-nt ACGTN, matches the Alithea BRB-seq kit constraint the
# STAR call hardcodes via --soloCBlen 14) means "no header"; anything else
# (a word, empty, or wrong length) means "has header".
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  MANIFEST_HEADER=1
else
  # Read line 1 as a pure string; guard with `|| true` so a missing file
  # doesn't trip `set -e` — the file-existence check happens earlier and the
  # real validation happens in Alithea's prepare_barcodes.sh downstream.
  line1="$(head -n1 -- "$MANIFEST" 2>/dev/null || true)"
  c2="${line1#*$'\t'}"; [[ "$c2" == "$line1" ]] && c2=""
  c2_clean="${c2%$'\r'}"; c2_up="$(printf '%s' "$c2_clean" | tr '[:lower:]' '[:upper:]')"
  if [[ ${#c2_up} -eq 14 && "$c2_up" =~ ^[ACGTN]+$ ]]; then
    MANIFEST_HEADER=0
    log "manifest : line 1 barcode cell '$c2_up' is a valid 14-nt ACGTN cell barcode - treating line 1 as DATA (no header)"
  else
    MANIFEST_HEADER=1
    log "manifest : line 1 barcode cell '${c2_clean:-<empty>}' is not a valid 14-nt ACGTN cell barcode - treating line 1 as HEADER (skipped)"
  fi
fi
if [[ "$MANIFEST_HEADER" -eq 0 ]]; then
  log "header   = NONE (manifest treated as headerless)"
else
  log "header   = YES   (first row of $MANIFEST is a header and is skipped)"
fi

#--- tool discovery (env-first) ----------------------------------------------
# Tool resolution policy (2026-10-22, user direction): use the pipeline tools
# from the brbseq Conda environment, not from copies that merely happen to sit
# on the server/PATH.  Previously the env was only consulted when STAR was NOT
# already in PATH, so a stray older system python3 won and crashed
# parse_libprep_xlsx.py (`from __future__ import annotations` needs Python 3.7+).
# Now the env is tried FIRST even when STAR is already on PATH.
have() { command -v "$1" >/dev/null 2>&1; }
_py_has_openpyxl() { [[ -n "$1" ]] && "$1" -c "import openpyxl" >/dev/null 2>&1; }
_py_modern() {  # True for a Python that supports `from __future__ import annotations`
  [[ -n "$1" ]] && "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 7) else 1)' >/dev/null 2>&1; }
_detect_tools() {
  STAR_BIN=""
  for cand in STAR-avx2 STAR-avx STAR; do
    if have "$cand"; then STAR_BIN="$cand"; break; fi
  done
  HAVE_STAR=0;      [[ -n "$STAR_BIN" ]] && have "$STAR_BIN" && HAVE_STAR=1
  HAVE_FASTQC=0;    have fastqc               && HAVE_FASTQC=1
  HAVE_SAMTOOLS=0;  have samtools             && HAVE_SAMTOOLS=1
  HAVE_RSEQC=0;     have read_distribution.py && HAVE_RSEQC=1
  HAVE_PICARD=0;    have picard               && HAVE_PICARD=1
  HAVE_FQTK=0;      have fqtk                 && HAVE_FQTK=1
  HAVE_JAVA=0;      have java                 && HAVE_JAVA=1
}

# Locate the brbseq conda env; first candidate with a bin/ wins:
#   1) explicit BRBSEQ_ENV (caller override, highest priority)
#   2) the active conda env, if it is named brbseq or brb
#   3) a brb*/brbseq env in the production anaconda envs/ directory — the
#      known ptegsrv3 path contains LITERAL square brackets, so it must stay
#      quoted (an unquoted '[' is a glob metacharacter).
_cands=( "${BRBSEQ_ENV:-}" )
if [[ ( "${CONDA_DEFAULT_ENV:-}" == "brbseq" || "${CONDA_DEFAULT_ENV:-}" == "brb" ) && -n "${CONDA_PREFIX:-}" ]]; then
  _cands+=( "$CONDA_PREFIX" )
fi
_cands+=( "/home/gecf/softwares/anaconda/[/home/gecf/softwares/anaconda]/envs/brbseq" )
BRBSEQ_ENV=""
for _c in "${_cands[@]}"; do
  [[ -n "$_c" && -d "${_c}/bin" ]] && { BRBSEQ_ENV="$_c"; break; }
done
if [[ -z "$BRBSEQ_ENV" ]]; then
  _ana="/home/gecf/softwares/anaconda/[/home/gecf/softwares/anaconda]"
  if [[ -d "${_ana}/envs" ]]; then
    for _base in brbseq brb; do
      if [[ -d "${_ana}/envs/${_base}/bin" ]]; then BRBSEQ_ENV="${_ana}/envs/${_base}"; break; fi
    done
  fi
fi
unset _c _cands _ana _base _env 2>/dev/null || true

if [[ -n "$BRBSEQ_ENV" ]]; then
  export PATH="${BRBSEQ_ENV}/bin:${PATH}"
  [[ -x "${BRBSEQ_ENV}/bin/python" ]] && PYTHON="${BRBSEQ_ENV}/bin/python"
  if [[ -n "$PYTHON" && -x "$PYTHON" ]] && ! _py_modern "$PYTHON"; then
    die "The 'brbseq' Conda environment at ${BRBSEQ_ENV} ships a Python that is too old for this pipeline (needs Python 3.7+, found: $PYTHON). Rebuild the environment with a modern Python:  conda env create -f ${SCRIPT_DIR}/brbseq-env.yaml  then  conda activate brbseq  and re-run; or point the pipeline at a working environment with BRBSEQ_ENV=<path>."
  fi
  _detect_tools
  log "env      = ${BRBSEQ_ENV}  (all pipeline tools resolved from this environment; the caller's PATH is not used for tools)"
else
  PYTHON="${PYTHON:-python3}"
  _detect_tools
  log "env      = (brbseq Conda environment not found — tools resolved from the current PATH)"
fi

if [[ "$HAVE_STAR" -eq 0 ]]; then
  die "STAR — the program that maps reads to the genome — was not found either in the brbseq Conda environment or on the current PATH. It normally ships inside the 'brbseq' Conda environment. Two ways to fix:  1) create that environment once from this script's folder:  conda env create -f ${SCRIPT_DIR}/brbseq-env.yaml  then  conda activate brbseq  and re-run;  2) or, if the environment already lives elsewhere on this machine, point the pipeline at it:  BRBSEQ_ENV=<path-to-env> brbseq_pipeline.sh ..."
fi

if [[ "$LIBPREP_MODE" -eq 1 ]] && ! _py_has_openpyxl "$PYTHON"; then
  die "Reading the lib-prep Excel workbook needs a Python that has the 'openpyxl' package, but the Python currently in use ($PYTHON) does not. The 'brbseq' Conda environment ships it — that environment is normally auto-detected, so this usually means the brbseq environment is missing or incomplete on this machine. Fix: create it once from this script's folder —  conda env create -f ${SCRIPT_DIR}/brbseq-env.yaml  then  conda activate brbseq  and re-run; or point the pipeline at it:  BRBSEQ_ENV=<path-to-env>."
fi

#--- generate manifest from lib-prep workbook (workbook mode) ----------------
# parse_libprep_xlsx.py needs python+openpyxl (verified above); the caller
# activates their own environment before launching.
if [[ "$LIBPREP_MODE" -eq 1 && ! -s "$MANIFEST" && "$DRY_RUN" -eq 0 ]]; then
  mkdir -p "$BARCODES_DIR"
  log "lib-prep workbook mode: generating manifest from $LIBPREP_XLSX"
  parse_args=( --xlsx "$LIBPREP_XLSX" --out "$MANIFEST" --provenance "$LIBPREP_PROV" )
  if [[ -n "$BARCODES_PLATE" ]]; then
    parse_args+=( --plate "$BARCODES_PLATE" )
  fi
  # Capture the parser's exit code explicitly (a bare "cmd; rc=$?" is killed
  # by set -e before the assignment; "if ! cmd" hides $? behind the inverted
  # status).  "cmd || rc=$?" is the set-e-safe idiom.
  parse_rc=0
  "$PYTHON" "$PARSE_LIBPREP" "${parse_args[@]}" || parse_rc=$?
  if [[ "$parse_rc" -ne 0 ]]; then
    die "The step that reads the lib-prep Excel workbook and builds the sample→barcode list failed (error code ${parse_rc}; the specific reason is printed in the lines above this message — start there). In this 'workbook' mode the pipeline normally finds the matching barcode plate file automatically using the 'Kit ref' row of the workbook, and expects a file like  barcodes_96_<family>_*.txt  in the plate directory /home/gecf/NGS_pipelines/BRB-seq/BRB_plates (override: \$BRB_PLATE_DIR). If that file is missing or named differently, either fix the plate directory, or pass the file explicitly with  --plate <file>. Then re-run."
  fi
  [[ -s "$MANIFEST" ]] || die "The sample→barcode list that should have been produced from the Excel workbook ended up completely empty (file: $MANIFEST) — meaning none of the samples in the workbook were successfully matched to a barcode. See the messages printed above this one for the specific mismatches (e.g. a well position in the workbook that does not exist on the plate file). Fix those in the workbook or the plate file, then delete this same output folder (or the stale file at $MANIFEST) and re-run so the list is generated afresh."
  log "manifest   = ${MANIFEST}  (generated from lib-prep workbook)"
elif [[ "$LIBPREP_MODE" -eq 1 && ! -s "$MANIFEST" && "$DRY_RUN" -eq 1 ]]; then
  log "manifest   = (workbook mode; would generate from lib-prep workbook — skipped in dry-run)"
elif [[ "$LIBPREP_MODE" -eq 1 ]]; then
  # manifest already on disk (resume).  Workbook mode ALWAYS emits
  # empty_well_<WELL> rows (unconditionally since 2026-10-22), so a manifest
  # without them predates that default.  FAIL rather than regenerate in
  # place: downstream resume steps may already depend on the
  # existing whitelist.  Start a new run (--run-name) to regenerate with
  # empty wells included.
  ew_in_manifest=0
  grep -q $'^empty_well_[A-H][0-9][0-9]	' "$MANIFEST" && ew_in_manifest=1
  if [[ "$ew_in_manifest" -eq 0 ]]; then
    die "The sample→barcode list already saved from a previous run of this experiment (file: $MANIFEST) does not include any empty-well entries — it was generated before including unoccupied plate wells became the default behaviour of this pipeline. The current pipeline would rebuild it with empty_well_<WELL> rows and a different barcode whitelist than the one this run directory was already computed with, so mixing would give an inconsistent result. To fix: start a fresh run with a new --run-name (the list will be regenerated including the empty wells), or, if you prefer to keep this run exactly as it was, remove and re-add the old manifest yourself — but the simpler, intended path is the fresh run."
  fi
  log "manifest   = ${MANIFEST}  (already present, reusing — includes $(grep -c $'^empty_well_' "$MANIFEST" 2>/dev/null || echo '?') empty-well row(s))"
fi

#--- banner -------------------------------------------------------------------
log "==== MERCURIUS DRUG-seq pipeline v${SCRIPT_VERSION} (dry_run=${DRY_RUN}) ===="
log "run name   = ${RUN_NAME}"
log "pool       = ${POOL}"
log "threads    = ${THREADS}"
log "R1 (CB+UMI)= ${R1_BARCODE}"
log "R2 (cDNA)  = ${R2_GENOMIC}"
log "manifest   = ${MANIFEST}  (header-mode=${MANIFEST_HEADER})"
log "STAR index = ${STAR_INDEX}   (star=${STAR_BIN:-NOT_FOUND})"
log "output     = ${RUN_DIR}"
log "optional:  fastqc=${RUN_OPT_FASTQ}  rseqc=${RUN_OPT_RSEQC}  demux_bam=${RUN_OPT_DEMUX_BAM}  demux_fq=${RUN_OPT_DEMUX_FQ}"
log "tools:     STAR=${HAVE_STAR}  fastqc=${HAVE_FASTQC}  samtools=${HAVE_SAMTOOLS}  rseqc=${HAVE_RSEQC}  picard=${HAVE_PICARD}  fqtk=${HAVE_FQTK}  java=${HAVE_JAVA}"
log "python     = ${PYTHON}"
if [[ -n "$LIBPREP_PROV" ]]; then
  log "lib-prep   = ${LIBPREP_XLSX}"
  if [[ -n "$BARCODES_PLATE" ]]; then
    log "plate      = ${BARCODES_PLATE}"
  else
    log "plate      = (not given — inferred from the workbook's 'Kit ref' row)"
  fi
  log "empty wells = ON (default: unoccupied plate wells are always emitted as empty_well_<WELL> rows)"
fi

#--- execute (or dry-run) ----------------------------------------------------
if [[ "$DRY_RUN" -eq 1 ]]; then
  log "---- DRY RUN: planned steps ----"
  hdr_flag=""
  [[ "$MANIFEST_HEADER" -eq 1 ]] && hdr_flag=" -H"
  log "  prepare_barcodes.sh -i <manifest> -w <OUT>/barcodes/barcode_whitelist.txt -f <OUT>/barcodes/fastq_sample_metadata.tsv -b <OUT>/barcodes/bam_barcodes.tsv${hdr_flag}"
  if [[ "$RUN_OPT_FASTQ" -eq 1 ]]; then
    log "  fastqc  --threads 32 --outdir <OUT>/fastqc <R1> <R2>"
  fi
  log "  STAR    (one STARsolo call, both dedup methods in list form)"
  log "      ${STAR_BIN} --runMode alignReads --runThreadN ${THREADS} --genomeDir ${STAR_INDEX} --readFilesIn <R2> <R1> --soloCBwhitelist <OUT>/barcodes/barcode_whitelist.txt --soloCBlen 14 --soloUMIlen 14 --soloUMIstart 15 --soloUMIdedup NoDedup 1MM_Directional --outSAMattributes '${SAM_ATTRS}' --quantMode GeneCounts --outFileNamePrefix <OUT>/star/${POOL}_ ...   ->  ${POOL}_Solo.out/Gene/raw/ (raw UMI matrix, kept in place)"
  if [[ "$RUN_OPT_RSEQC" -eq 1 ]]; then
    if [[ -n "$BED_FILE" ]]; then
      log "  RSeQC   reference (given --bed-file): ${BED_FILE}"
    else
      log "  RSeQC   bake:   gtfToGenePred ${RSEQC_GTF} tmp.gp ; genePredToBed tmp.gp <OUT>/rseqc/$(basename "${RSEQC_GTF%.gz}").bed12   (UCSC chain)"
    fi
    log "      samtools index <OUT>/star/${POOL}_Aligned.sortedByCoord.out.bam   [-@ 16]"
    log "      read_distribution.py -i <OUT>/star/${POOL}_Aligned.sortedByCoord.out.bam -r <bed12> > <OUT>/rseqc/${POOL}_read_distribution.txt"
  fi
  log "  dups    ${PYTHON} <UTILS>/get_duplication_statistics.py --raw_dir <OUT>/star/${POOL}_Solo.out/Gene/raw --dedup_method 1MM_Directional --output_csv <OUT>/duplication_statistics.csv"
  log "  rename  ${PYTHON} <UTILS>/rename_matrix.py --raw_dir <OUT>/star/${POOL}_Solo.out/Gene/raw --mapping <manifest> --out-dir <OUT>/star/samples -> barcodes named by Sample_ID (raw barcode matrix kept untouched in ${POOL}_Solo.out/Gene/raw)"
  if [[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]]; then
    log "  demuxBAM picard FilterSamReads per sample in <OUT>/barcodes/bam_barcodes.tsv MINUS empty_well_* rows (samples only, CR tag) -> <OUT>/demux_bam/<sample>.bam"
  fi
  if [[ "$RUN_OPT_DEMUX_FQ" -eq 1 ]]; then
    log "  demuxFQ fqtk demux <R1> <R2> --samples fastq_sample_metadata.tsv MINUS empty_well_* rows -> <OUT>/demux_fastq"
  fi
  log "---- DRY RUN: no files modified ----"
  exit 0
fi

# --- actually mkdir (skipped in dry-run) -----------------------------------
mkdir -p "$BARCODES_DIR" "$STAR_DIR" "$MATRIX_DIR"
[[ "$RUN_OPT_FASTQ" -eq 1 ]] && mkdir -p "$FASTQC_DIR"
[[ "$RUN_OPT_RSEQC" -eq 1 ]] && mkdir -p "$RSEQC_DIR"
[[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]] && mkdir -p "$DMUX_BAM_DIR"
[[ "$RUN_OPT_DEMUX_FQ" -eq 1 ]] && mkdir -p "$DMUX_FQ_DIR"
LOGS_DIR="${RUN_DIR}/logs"
mkdir -p "$LOGS_DIR"

#--- logging layout -------------------------------------------------------------
# pipeline.log      : supervisor log — RUN/SKIP/done lines, step-log pointers,
#                     die()/warn() messages.  Does NOT contain raw tool output.
# logs/<step>.log   : one file per pipeline step — tool stdout + stderr.
# star/STAR.log     : STAR's own log (already co-located with STAR output; unchanged).
# rseqc/*.txt etc.  : data products, not logs.
LOG_FILE="${RUN_DIR}/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2> >(tee -a "$LOG_FILE" >&2)
# Drain the tee pipes before the shell dies, or the final buffered lines can
# be lost on exit (process-substitution pipe flush is async).
trap 'sleep 0.5' EXIT
log "pipeline.log = ${LOG_FILE}  (supervisor — no raw tool output); step logs in ${LOGS_DIR}/"

# --- run_or_skip -------------------------------------------------------------
# Execute the command unless the marker file already exists (auto mode).
# Failure policy: a failing OPTIONAL step logs a warning and does not abort;
# a failing MANDATORY step aborts the pipeline.
#   $1=step name  $2=marker file (empty => always run)
#   $3=tool name  (empty => no tool gate, e.g. a Python helper)
#   $4=mandatory|optional
#   [$4.5]  optionally the literal word  stdout_to  followed by the product
#           file path — tells run_or_skip to redirect the command's stdout
#           to that file (the step's data product) and stderr to the step log.
#           Without this, both stdout and stderr go to the step log.
#   $5...=command
run_or_skip() {
  local name="$1"; shift
  local marker="${1:-}"; shift
  local tool="${1:-}"; shift
  local mode="${1:-mandatory}"; shift
  # Optional "sentinel" flag: the marker is a pure completion sentinel (empty
  # file, created by touch) and is tested with -e (exists) rather than -s
  # (non-empty).  Use it for steps whose marker has no data content.
  local sentinel=0
  if [[ "${1:-}" == "sentinel" ]]; then sentinel=1; shift; fi
  # Optional 2-token flag:  stdout_to <path>  → command's stdout is the
  # step's data product; run_or_skip does both redirects explicitly so the
  # log calls above it still land in pipeline.log, not the product.
  local product_path=""
  if [[ "${1:-}" == "stdout_to" ]]; then
    product_path="${2:?run_or_skip: stdout_to requires a file path}"
    shift 2
  fi
  local cmd=( "$@" )
  local marker_ok=0
  if [[ -n "$marker" ]]; then
    if [[ "$sentinel" -eq 1 ]]; then
      [[ -e "$marker" ]] && marker_ok=1     # empty sentinel: existence is the signal
    else
      [[ -s "$marker" ]] && marker_ok=1     # data marker: must be non-empty
    fi
  fi
  if [[ "$marker_ok" -eq 1 ]]; then
    log "SKIP ${name} (already done: ${marker})"
    return 0
  fi
  if [[ -n "$tool" ]] && ! have "$tool"; then
    if [[ "$mode" == "mandatory" ]]; then
      die "The step '${name}' needs a program called '${tool}' to run, but that program is not installed in the current environment. It normally ships inside the 'brbseq' Conda environment — activate it with 'conda activate brbseq' (create it first with 'conda env create -f brbseq-env.yaml' if you have never made it), then re-run the pipeline. Because this is a required step, the pipeline has stopped here."
    else
      warn "Skipping optional step '${name}': the program it needs ('${tool}') is not installed in the current environment. This does not stop the pipeline — everything else will still run — but the output this step would have produced will be missing. If you want it, install it into the 'brbseq' environment ('conda activate brbseq'; create it with 'conda env create -f brbseq-env.yaml' if needed) and re-run."
      return 0
    fi
  fi
  local step_name_safe
  step_name_safe="$(printf '%s' "$name" | sed 's/[^[:alnum:]]/_/g')"
  local step_log="${LOGS_DIR}/${step_name_safe}.log"
  log "RUN  ${name}: ${cmd[*]}"
  log "  step log -> ${step_log}"
  local rc=0
  if [[ -n "$product_path" ]]; then
    # stdout is the step's data product; stderr -> step log
    "${cmd[@]}" >"$product_path" 2>"$step_log" || rc=$?
  else
    # all output -> step log (removed from pipeline.log)
    "${cmd[@]}" >"$step_log" 2>&1 || rc=$?
  fi
  if [[ "$rc" -ne 0 ]]; then
    log "  step log tail (last 30 lines of ${step_log}):"
    tail -30 "$step_log" >&2 2>/dev/null || true
    if [[ "$mode" == "mandatory" ]]; then
      die "The required step '${name}' ran into an error and could not be completed, so the pipeline has stopped. Full output is in the step log: ${step_log}. Re-run the pipeline with --resume-dir pointing at this same run folder so the steps that already finished are not redone."
    else
      warn "The optional step '${name}' ran into an error, so it did not finish — but this does not stop the pipeline; everything else will still run. The reason is in the step log: ${step_log}. If you want this step's output, fix the underlying problem and re-run (it is safe to re-run just this one)."
    fi
    # A FAILED step must not keep (or gain) its completion marker, otherwise a
    # later resume would skip a step that never actually finished.
    if [[ -n "$marker" && -e "$marker" ]]; then
      rm -f -- "$marker"
      log "  (completion marker removed on failure: ${marker} — this step will be re-run on the next resume)"
    fi
  else
    [[ -n "$marker" ]] && touch "$marker"
  fi
}

#=== Creating barcode files  [MANDATORY] ==================================
log ""; log "Creating barcode files (mandatory)"
WHITELIST="${BARCODES_DIR}/barcode_whitelist.txt"
FQ_META="${BARCODES_DIR}/fastq_sample_metadata.tsv"
BAM_BARCODES="${BARCODES_DIR}/bam_barcodes.tsv"
HDR_ARGS=()
[[ "$MANIFEST_HEADER" -eq 1 ]] && HDR_ARGS+=("-H")

run_or_skip "prepare_barcodes" "$WHITELIST" "" "mandatory" \
  bash "$PREPARE_BARCODES" \
    -i "$MANIFEST" \
    -w "$WHITELIST" \
    -f "$FQ_META" \
    -b "$BAM_BARCODES" \
    "${HDR_ARGS[@]}"
# post-check
[[ -s "$WHITELIST"    ]] || die "The barcode file needed to align the reads ended up empty (file: $WHITELIST) — meaning none of the samples in your list were successfully turned into a usable barcode. See the message just above for the reason (usually a problem with the sample→barcode list itself). Fix the source list and re-run."
[[ -f "$FQ_META"      ]] || die "The per-sample barcode list used later to split the raw sequence files (file: $FQ_META) was not produced. See the message just above for the reason, fix it, and re-run."
[[ -f "$BAM_BARCODES" ]] || die "The per-sample barcode list used later to split the aligned results file (file: $BAM_BARCODES) was not produced. See the message just above for the reason, fix it, and re-run."
log "barcode prep done: whitelist=$(wc -l < "$WHITELIST") barcodes; fastq_meta+bam_map ready"

#=== FASTQ Quality Assessment  [OPTIONAL] ==================================
if [[ "$RUN_OPT_FASTQ" -eq 1 ]]; then
  log ""; log "fastqc (optional)"
  FQ_MARKER="${FASTQC_DIR}/QC_done"
  run_or_skip "fastqc" "$FQ_MARKER" "fastqc" "optional" sentinel \
    fastqc --threads 32 --outdir "$FASTQC_DIR" "$R1_BARCODE" "$R2_GENOMIC"
  log "fastqc outputs in: $FASTQC_DIR"
else
  log ""; log "fastqc — skipped (--skip-fastq)"
fi

#=== STARsolo  [MANDATORY] =================================================
log ""; log "STARsolo alignment + UMI count matrix (mandatory)"
[[ "$HAVE_STAR" -eq 1 ]] || die "The STAR aligner (the program that maps reads to the genome) was not found on this machine. It ships inside the 'brbseq' Conda environment. To fix: first make sure that environment exists and is active — run 'conda activate brbseq' — then try again. If you have never created it, create it once with: conda env create -f brbseq-env.yaml  (that file sits next to this pipeline script), then 'conda activate brbseq', then re-run."
[[ -d "$STAR_INDEX" ]] || die "The STAR genome index directory could not be found: $STAR_INDEX. Check the --star-genome value (or the built-in path for your --species) and re-run. If the index was deleted or moved, rebuild it with build_star_index.sh."
[[ -f "${STAR_INDEX}/Genome" ]] || warn "The STAR index folder $STAR_INDEX appears incomplete — the 'Genome' file is missing, and STAR will not be able to start if the index is unfinished. If you are unsure this index is finished, rebuild it with build_star_index.sh, then re-run."
[[ -f "${STAR_INDEX}/SA"     ]] || warn "The STAR index folder $STAR_INDEX appears incomplete — the 'SA' (sequence-to-genome) file is missing, and STAR will not be able to start if the index is unfinished. If you are unsure this index is finished, rebuild it with build_star_index.sh, then re-run."

# Fixed STAR build (brbseq conda env) => deterministic output names: STAR
# appends its fixed suffixes to --outFileNamePrefix, and the prefix carries a
# trailing underscore (<POOL>_) so the names stay legible:
#   ${STAR_DIR}/${POOL}_Aligned.sortedByCoord.out.bam
#   ${STAR_DIR}/${POOL}_Solo.out/Gene/raw/   <- the raw UMI matrix, its only location
# Both deduplication methods are requested in a single invocation via the list
# form of --soloUMIdedup (NoDedup 1MM_Directional). One alignment pass yields a
# per-method umiDedup-<method>.mtx in ${POOL}_Solo.out/Gene/raw/, avoiding a
# second, nondeterministic alignment of the same reads.
STAR_BASE="${STAR_DIR}/${POOL}_"
STAR_BAM="${STAR_DIR}/${POOL}_Aligned.sortedByCoord.out.bam"
RAW_DIR="${STAR_DIR}/${POOL}_Solo.out/Gene/raw"

#--- helpers: one STAR run, both matrices -------------------------------------
# The list of --soloUMIdedup methods requested in the single STAR call.
DEDUP_METHODS=( "NoDedup" "1MM_Directional" )

# Run STAR once, requesting both dedup methods via the list form of
# --soloUMIdedup.  Outputs a BAM + a SOLOout/Gene/raw dir containing
# umiDedup-<method>.mtx for every method in DEDUP_METHODS.
# Args:  $1=out_prefix_base  (bare --outFileNamePrefix; STAR appends its own suffixes)
star_solo() {
  local prefix="$1"
  local log_file="${STAR_DIR}/STAR.log"
  # --outSAMattributes emits the Alithea production tag set
  # (NH HI nM AS CR UR CB UB GX GN sS sQ sM) -- fixed.
  # STAR 2.7.11b requires each tag as a SEPARATE argument; passing the
  # whole list as one quoted string is a FATAL INPUT ERROR
  # ("unknown/unimplemented SAM attribute (tag): NH HI nM ...").
  local -a sam_attr_tokens=()
  read -r -a sam_attr_tokens <<< "${SAM_ATTRS}"
  # Build the invocation as one array: logged verbatim (pipeline.log) and
  # written as the first line of STAR.log for provenance.
  local -a STAR_CMD=(
    "${STAR_BIN}"
    --runMode alignReads
    --runThreadN "${THREADS}"
    --genomeDir "${STAR_INDEX}"
    --readFilesIn "${R2_GENOMIC}" "${R1_BARCODE}"
    --readFilesCommand zcat
    --soloCBwhitelist "${WHITELIST}"
    --soloCBstart 1
    --soloCBlen 14
    --soloUMIstart 15
    --soloUMIlen 14
    --soloUMIdedup "${DEDUP_METHODS[@]}"
    --soloType CB_UMI_Simple
    --soloCBmatchWLtype 1MM
    --soloStrand Forward
    --soloCellFilter None
    --soloBarcodeReadLength 0
    --outFilterMultimapNmax 1
    --outSAMattributes "${sam_attr_tokens[@]}"
    --outSAMmapqUnique 60
    --outSAMunmapped Within
    --outSAMtype BAM SortedByCoordinate
    --outBAMsortingThreadN 1
    --outSAMmultNmax 1
    --outFilterIntronMotifs None
    --outFileNamePrefix "${prefix}"
    --quantMode GeneCounts
  )
  log "RUN  STARsolo: ${STAR_CMD[*]}"
  printf '%s\n' "${STAR_CMD[*]}" > "${log_file}"
  if "${STAR_CMD[@]}" >> "${log_file}" 2>&1; then
    return 0
  else
    log "The STAR alignment/counting step (STARsolo, UMI counting methods: ${DEDUP_METHODS[*]}) reported that it failed. The last lines of its own log file are printed directly above this line — read them for the specific reason (for example too few reads, an index it cannot read, or a run out of memory). Fix that, then re-run the pipeline with --resume-dir pointing at this same run folder so the steps that already finished are not redone."
    tail -20 "${log_file}" >&2 || true
    return 1
  fi
}

S14_MARKER="${STAR_DIR}/_s14_done"
S14_DONE=0
if [[ -f "$S14_MARKER" && -s "$S14_MARKER" ]]; then S14_DONE=1; fi

if [[ "$S14_DONE" -eq 0 ]]; then
  # One STAR invocation, --soloUMIdedup in list form (both methods).
  star_solo "${STAR_BASE}"
  [[ -s "$STAR_BAM" ]] || die "The alignment step (STAR) ran to completion, but the aligned reads file it should have produced ($STAR_BAM) was not found or is empty. This usually means STAR itself stopped early with an error — see the messages printed just above this line, and the STAR log ($STAR_DIR/STAR.log), for the specific reason. Fix that, then re-run the pipeline (use --resume-dir pointing at this same folder so earlier steps are not redone)."
  echo "$(date +%FT%T)" > "${STAR_DIR}/_s14_done"
fi
log "STAR BAM: $STAR_BAM"
log "raw UMI matrix dir: $RAW_DIR"

#--- bake helper: GTF -> BED12 (UCSC genePred chain) -------------------------
# Same conversion as flashseq_rseqc.sh so both pipelines produce identical
# BED12 from the same GTF.  Writes $RSEQC_DIR/<gtf-basename>.bed12
# (.gz GTFs are decompressed on the fly).  Sets global RSEQC_BED to the
# output path (NOT via stdout capture — log() writes to stdout, so the value
# is returned in a variable to avoid mixing prose into the path).
#
# Chain (UCSC genePred, exactly what UCSC docs recommend and what RSeQC
# expects):
#          gtfToGenePred GTF   gp   (transcript-level genePred, 8 cols)
#          genePredToBed gp    bed  (12-col BED12)
#
# WHY UCSC CHAIN, NOT gffread --bed: for transcripts whose GTF has no explicit
# CDS feature (many lncRNAs / pseudogenes / many BRB-seq annotations), gffread
# fills thickStart/thickEnd over the WHOLE transcription-unit span — RSeQC
# read_distribution.py then classifies every exon as "CDS_exons". The UCSC
# chain leaves thickStart=thickEnd (empty CDS), which is the biologically
# faithful reading of "no CDS in this transcript". That is the correct,
# RSeQC-recommended partition for CDS/5'UTR/3'UTR/intron bucketing.
#
# Two live bugs to guard against (see flashseq_rseqc.sh bake_bed12):
#   1. genePredToBed takes <in.gp> <out.bed> positionally — it DOES NOT print
#      to stdout, so we must call it as `genePredToBed "$gp" "$bed"` (not
#      `genePredToBed "$gp" > "$bed"`).
#   2. genePredToBed prints usage to stderr and RETURNS 0 in the usage-error
#      case, so `if ! genePredToBed ...` will not catch failure — we MUST
#      check the output is non-empty (`[[ -s "$bed" ]]`) after every step.
RSEQC_BED=""
bake_rseqc_bed12() {
  local gtf="$1" out
  out="${RSEQC_DIR}/$(basename "${gtf%.gz}").bed12"
  if [[ -s "$out" ]]; then
    log "BED12 already baked (resume): $out  (from $gtf)"
    RSEQC_BED="$out"
    return 0
  fi
  # Tool gate: BOTH UCSC chain binaries must be reachable. (UCSC tools have
  # no --version flag — presence on PATH is what we check here; provenance is
  # logged by the caller via the no-arg usage string.)
  local _missing=""
  command -v gtfToGenePred >/dev/null 2>&1 || _missing="$_missing gtfToGenePred"
  command -v genePredToBed >/dev/null 2>&1 || _missing="$_missing genePredToBed"
  [[ -z "$_missing" ]] || die "The read-distribution step needs one or two UCSC conversion tools (missing: ${_missing// /, }) to turn the GTF gene model (file: $gtf) into the BED12 format that RSeQC reads, and they are currently not installed in this environment. Fix: install them into the 'brbseq' Conda environment with 'conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed', then run 'conda activate brbseq' again and re-run the pipeline. Alternatively, skip the automatic conversion entirely by pointing the pipeline at a finished BED12 file: --bed-file <file.bed12>."
  local tmpdir rc=1
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/brbseq_bed12.XXXXXX")"
  local gp="$tmpdir/genes.gp" bedtmp="$tmpdir/genes.bed"
  # -- GTF -> genePred -------------------------------------------------------
  # pipefail/errexit are off inside if-conditions: capture rc explicitly.
  if   [[ "$gtf" == *.gz  ]]; then gzip -dc "$gtf"  | gtfToGenePred -  "$gp"; rc=$?
  elif [[ "$gtf" == *.bz2 ]]; then bzcat -d "$gtf"  | gtfToGenePred -  "$gp"; rc=$?
  else                             gtfToGenePred     "$gtf"           "$gp"; rc=$?; fi
  (( rc == 0 ))      || { rm -rf "$tmpdir"; die "Could not convert the GTF gene-model file into an intermediate genePred format (this is step 1 of 2 before RSeQC) — the conversion tool gtfToGenePred stopped with an error on file: $gtf. This usually means the GTF file itself is malformed (for example, missing or broken 'gene'/'exon' lines, or unusual characters in chromosome names). Fix or re-download the GTF file, or supply a finished BED12 instead with --bed-file <file.bed12>."; }
  [[ -s "$gp" ]]     || { rm -rf "$tmpdir"; die "The conversion of the GTF gene-model file produced an empty intermediate result (step 1 of 2 before RSeQC) for file: $gtf — this usually means the GTF does not contain any recognisable 'gene'/'exon' entries, so there was nothing to convert. Double-check that this is really a GTF annotation file (not, say, a FASTA genome file) and re-run, or supply a finished BED12 instead with --bed-file <file.bed12>."; }
  # -- genePred -> BED12 (genePredToBed takes <in.gp> <out.bed>; output file is
  # NOT stdout, and it returns 0 on a usage error, so check the file is
  # non-empty, not just the exit code) --------------------------------------
  if ! genePredToBed "$gp" "$bedtmp"; then
    rm -rf "$tmpdir"
    die "Could not finish converting the gene model into the BED12 format RSeQC needs (this is step 2 of 2) for file: $gtf. The conversion tool itself reported an error. If you have a BED12 file already prepared for this same gene model, you can point the pipeline at it directly with --bed-file <file.bed12> and skip this conversion entirely.";
  fi
  [[ -s "$bedtmp" ]] || { rm -rf "$tmpdir"; die "The BED12 conversion finished, but it produced an empty result for file: $gtf — meaning the gene model did not get converted into the format RSeQC needs. If you have a BED12 file already prepared for this same gene model, you can point the pipeline at it directly with --bed-file <file.bed12> instead."; }
  local cols
  cols="$(awk -F'\t' '!/^#/{print NF; exit}' "$bedtmp")"
  if [[ "${cols:-0}" -lt 12 ]]; then
    rm -rf "$tmpdir"
    die "The finished BED12 file was produced, but its format is wrong (it has ${cols:-0} columns per line; a valid BED12 gene-model file needs exactly 12). This points to a problem in the original GTF gene-model file: $gtf. If you have a known-good BED12 for this gene model, point the pipeline at it with --bed-file <file.bed12> instead."
  fi
  mv -f "$bedtmp" "$out"
  rm -rf "$tmpdir"
  local rows; rows="$(wc -l < "$out" | tr -d '[:space:]')"
  log "baked BED12 (RSeQC): $out   [from GTF: $gtf]  (${rows} rows x ${cols} cols)"
  log "provenance:  gtfToGenePred $gtf tmp.gp ; genePredToBed tmp.gp $out   (UCSC chain)"
  RSEQC_BED="$out"
}

#=== Read distribution (RSeQC)  [OPTIONAL] ==================================
if [[ "$RUN_OPT_RSEQC" -eq 1 ]]; then
  log ""; log "read_distribution (optional)"
  RSEQC_OUT="${RSEQC_DIR}/${POOL}_read_distribution.txt"
  # Resolve the BED12: explicit --bed-file wins; otherwise bake from the GTF
  # shipped inside the STAR index dir (single source of truth = that GTF).
  if [[ -n "$BED_FILE" ]]; then
    log "RSeQC reference in use (as given --bed-file): $BED_FILE"
  else
    [[ -n "$RSEQC_GTF" ]] || die "The read-distribution step needs a gene-model file in a specific format (BED12), but the pipeline could not find the source file (a GTF annotation) to build one from, and no finished BED12 was given either. Normally build_star_index.sh saves the GTF inside the STAR index folder at index-build time, so check that $STAR_INDEX still contains a .gtf (or .gtf.gz) file; if it does not, either rebuild the index so it is saved there again, or supply a finished BED12 directly with --bed-file <file.bed12>."
    # Bake only when RSeQC actually runs; the bake is idempotent
    # (bake_rseqc_bed12 reuses an existing .bed12).  Sets RSEQC_BED.
    bake_rseqc_bed12 "$RSEQC_GTF"
    [[ -n "$RSEQC_BED" && -s "$RSEQC_BED" ]] || die "The pipeline tried to build the BED12 gene-model file that RSeQC needs, but the result was missing or zero-sized. This is the same class of problem as a broken/malformed GTF; if you have a known-good BED12 for this gene model, point the pipeline at it with --bed-file <file.bed12> instead."
    BED_FILE="$RSEQC_BED"
  fi
  log "RSeQC reference in use: $BED_FILE"
  log "RSeQC provenance: baked from GTF $RSEQC_GTF   (recipe: gtfToGenePred <gtf> tmp.gp ; genePredToBed tmp.gp <gtf>.bed12   [UCSC chain])"
  run_or_skip "samtools index" "${STAR_BAM}.bai" "samtools" "optional" \
    samtools index -@ 16 "$STAR_BAM"
  run_or_skip "read_distribution" "$RSEQC_OUT" "read_distribution.py" "optional" \
    stdout_to "$RSEQC_OUT" \
    read_distribution.py -i "$STAR_BAM" -r "$BED_FILE"
  if [[ -s "$RSEQC_OUT" ]]; then
    log "distribution report: $RSEQC_OUT"
  fi
  # Post-use cleanup: remove the baked BED12 (a genome-wide intermediate, not a
  # per-sample result) now that read_distribution has consumed it. Only ours —
  # a user-supplied --bed-file leaves RSEQC_BED empty and is left untouched.
  if [[ -n "$RSEQC_BED" && -f "$RSEQC_BED" ]]; then
    rm -f "$RSEQC_BED"
    log "removed baked BED12 intermediate (cleaned up after use): $RSEQC_BED"
  fi
else
  log ""; log "RSeQC — skipped (--skip-rseqc, or auto-skipped: no gene model / read-distribution tools)"
fi

#--- dup-stats helpers (kept here to anchor the section) --------------------
#=== Duplication statistics  [MANDATORY] ==================================
log ""; log "duplication statistics (mandatory)"
for req in "$RAW_DIR/barcodes.tsv" \
           "$RAW_DIR/umiDedup-NoDedup.mtx" \
           "$RAW_DIR/umiDedup-1MM_Directional.mtx"; do
  [[ -s "$req" ]] || die "The step that measures how duplicated / how 'empty' each sample's reads are needs three input files (the raw UMI count tables STAR produces), and one of them is missing or zero-sized: $req. All three must exist in $RAW_DIR before this step can run — this points to a problem earlier in the pipeline (the STAR / UMI-counting step). See the messages printed above for the likely reason; fix that, then re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
done
run_or_skip "duplication stats" "$DUP_CSV" "" "mandatory" \
  "$PYTHON" "$DUP_STATS" \
    --raw_dir "$RAW_DIR" \
    --dedup_method "1MM_Directional" \
    --manifest "$MANIFEST" \
    --output_csv "$DUP_CSV"
[[ -s "$DUP_CSV" ]] || die "The step that measures how duplicated / how 'empty' each sample's reads are did not produce its result file ($DUP_CSV). See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
log "duplication stats: $DUP_CSV"
if [[ -f "${DUP_CSV%.*}.summary.csv" ]]; then
  log "summary           : ${DUP_CSV%.*}.summary.csv"
fi

#=== Sample-named count matrix  [MANDATORY] =================================
# Rename the barcode-labelled matrix columns (barcodes.tsv in
# ${POOL}_Solo.out/Gene/raw) to the Sample_IDs from the run manifest:
#   workbook mode (always): two outputs from the SAME single STAR matrix
#     (one STAR call — the rename is a cheap post-processing step, not a
#     second alignment):
#     ${STAR_DIR}/samples_all_wells/  samples + empty_well_<WELL> (debug:
#                                     where the unoccupied-well reads went)
#     ${STAR_DIR}/samples/            samples only (empty_well_* excluded);
#                                     the final UMI/reads CSVs are built
#                                     from THIS directory
#   --barcode-tsv mode:
#     ${STAR_DIR}/samples/            samples only (there is no plate and no
#                                     empty-well concept); samples_all_wells
#                                     is simply absent
# The raw barcode-labelled matrix in ${POOL}_Solo.out/Gene/raw is left untouched —
# in workbook mode it is the full-plate view (up to 96 barcode columns).
log ""
log "sample-named count matrix (mandatory)"
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  run_or_skip "rename matrix (all wells: samples + empty wells)" "${RENAME_ALL_DIR}/barcodes.tsv" "" "mandatory" \
    "$PYTHON" "$RENAME_MATRIX" \
      --raw-dir   "$RAW_DIR" \
      --mapping   "$MANIFEST" \
      --out-dir   "$RENAME_ALL_DIR" \
      --dedup-methods "${DEDUP_METHODS[*]}"
  [[ -s "${RENAME_ALL_DIR}/barcodes.tsv" ]] || die "The full-plate (all-wells) sample-named matrix (for debugging empty-well reads) did not finish: ${RENAME_ALL_DIR}/barcodes.tsv is missing or empty. See the message(s) printed just above for the specific reason — often a sample in your list does not match any barcode found in the data. Fix and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
  log "all-wells matrix dir: ${RENAME_ALL_DIR}  (barcode->sample map: ${RENAME_ALL_DIR}/sample_barcode_map.tsv)"
fi
run_or_skip "rename matrix (samples only)" "${RENAME_DIR}/barcodes.tsv" "" "mandatory" \
  "$PYTHON" "$RENAME_MATRIX" \
    --raw-dir   "$RAW_DIR" \
    --mapping   "$MANIFEST" \
    --out-dir   "$RENAME_DIR" \
    --dedup-methods "${DEDUP_METHODS[*]}" \
    --exclude-prefix empty_well_
[[ -s "${RENAME_DIR}/barcodes.tsv" ]] || die "The step that gives every sample its proper name (matching the Sample_ID from your list) did not finish: the renamed barcode list ($RENAME_DIR/barcodes.tsv) is missing or empty. See the message(s) printed just above for the specific reason — often a sample in your list does not match any barcode found in the data, or vice versa. Fix the mismatch in your sample list and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
[[ -s "${RENAME_DIR}/umiDedup-1MM_Directional.mtx" ]] || die "The step that gives every sample its proper name (matching the Sample_ID from your list) did not finish: the renamed UMI count table ($RENAME_DIR/umiDedup-1MM_Directional.mtx) is missing or empty. See the message(s) printed just above for the specific reason — often a sample in your list does not match any barcode found in the data, or vice versa. Fix the mismatch in your sample list and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
log "sample-named matrix dir: ${RENAME_DIR}"
log "  (barcode->sample map: ${RENAME_DIR}/sample_barcode_map.tsv; raw matrix unchanged in ${RAW_DIR})"

#=== FastReadCounter (read-distribution report)  [MANDATORY] ==============
# FRC + report_stats — port of the legacy BRB-seq pipeline's chain:
#   FastReadCounter → rename_count_table (barcode → Sample_ID) → report_stats
# (old pipeline lines 592-618 + the report_stats function).  The invariant:
# FRC reads the same STAR BAM, the same GTF (the one STAR used), and the same
# 96-barcode whitelist as the legacy run, so the well grid in
# <pool>_report_summary.csv stays comparable with legacy runs.
#
# Writes:  <run>/fastreadcounter/                    FRC outputs (counts.*)
#          <run>/files_for_report_and_CLC/<pool>_all_stats.csv
#          <run>/files_for_report_and_CLC/<pool>_report_summary.csv
#
# The jar path is the production location (same as the legacy pipeline);
# update FRC_JAR here when the jar is upgraded.
FRC_JAR="/home/gecf/softwares/FastReadCounter/FastReadCounter-1.2.0.jar"
[[ -f "$FRC_JAR" ]] || die "The FastReadCounter program (the jar counted in the legacy pipeline) could not be found at $FRC_JAR. It normally ships on the production server at that path. If it has moved, update the FRC_JAR line near the 'FastReadCounter' section of the pipeline script, then re-run."

# GTF = the gene annotation STAR itself used (registered in the index);
# STAR writes it VERBATIM — often RELATIVE (e.g. "Mus_musculus.GRCm39.116/...gtf"),
# unresolvable from our CWD. Fallback: the single .gtf shipped IN the index dir.
FRC_GTF_RAW="$(awk '$1=="sjdbGTFfile"{print $2; exit}' "${STAR_INDEX}/genomeParameters.txt" 2>/dev/null || true)"
if [[ -n "$FRC_GTF_RAW" && -f "$FRC_GTF_RAW" ]]; then
    FRC_GTF="$FRC_GTF_RAW"
else
    FRC_GTF="$(ls "${STAR_INDEX}"/*.gtf 2>/dev/null | head -1 || true)"
fi
[[ -n "$FRC_GTF" && -f "$FRC_GTF" ]] || die "FastReadCounter needs the gene-annotation model (a GTF file) that STAR aligned against: sjdbGTFfile in ${STAR_INDEX}/genomeParameters.txt is '${FRC_GTF_RAW:-<absent>}' (not resolvable) and no .gtf was found in ${STAR_INDEX}. Put the GTF in the index dir or fix the path in genomeParameters.txt."
log "FastReadCounter gene model (GTF): $FRC_GTF"

FRC_DIR="${RUN_DIR}/fastreadcounter"
REPORTS_DIR="${RUN_DIR}/files_for_report_and_CLC"
FRC_DETAILED="${FRC_DIR}/counts.detailed.txt"
FRC_RENAMED="${FRC_DIR}/counts.detailed_withsamplesID.txt"
ALL_STATS="${REPORTS_DIR}/${POOL}_all_stats.csv"
REPORT_CSV="${REPORTS_DIR}/${POOL}_report_summary.csv"
mkdir -p "$FRC_DIR" "$REPORTS_DIR"

log ""
log "FastReadCounter + report (mandatory)"

#--- FRC call (one line, as in the legacy pipeline) -------------------------
run_or_skip "FRC FastReadCounter" "$FRC_DETAILED" "java" "mandatory" \
  java -jar "$FRC_JAR" \
    --bam    "$STAR_BAM" \
    --gtf    "$FRC_GTF" \
    --umi-dedup none \
    --barcodeFile "$WHITELIST" \
    -o "$FRC_DIR"
[[ -s "$FRC_DETAILED" ]] || die "FastReadCounter did not finish cleanly (its result file ${FRC_DETAILED} is missing or empty). See the FRC step log in ${LOGS_DIR} for the specific reason, fix it, and re-run (use --resume-dir pointing at this same folder so the earlier steps are not redone)."

#--- rename_count_table (ported): barcode columns → Sample_ID ---------------
# The all-wells barcode→name table from the rename-matrix step already covers
# all 96 plate barcodes (real samples + empty_well_<well>).  Falls back to
# the samples-only map in --barcode-tsv mode (no plate / empty-well concept).
SAMPLE_MAP="${RENAME_ALL_DIR}/sample_barcode_map.tsv"
if [[ ! -s "$SAMPLE_MAP" ]]; then
  SAMPLE_MAP="${STAR_DIR}/samples/sample_barcode_map.tsv"
fi
[[ -s "$SAMPLE_MAP" ]] || die "The barcode→sample-name table needed to name the FastReadCounter columns is missing: ${SAMPLE_MAP}. It is written by the sample-renaming step earlier in the pipeline — check the messages printed above for why that step did not finish."

if [[ -s "$FRC_RENAMED" ]]; then
  log "SKIP rename FRC count table (already done: $FRC_RENAMED)"
else
  log "RUN  rename FRC count table (barcode -> Sample_ID)"
  awk -v lookup="$SAMPLE_MAP" '
  BEGIN {
    FS = OFS = "\t"
    while ((getline line < lookup) > 0) {
      split(line, a, "\t")
      map[a[1]] = a[2]
    }
    close(lookup)
  }
  NR == 1 { for (i = 1; i <= NF; i++) if ($i in map) $i = map[$i] }
  { print }
  ' "$FRC_DETAILED" > "$FRC_RENAMED"
  [[ -s "$FRC_RENAMED" ]] || die "The step that names the FastReadCounter columns by sample (barcode → Sample_ID) did not produce a result file ($FRC_RENAMED). See the messages printed above for the reason, fix it, and re-run."
fi

#--- report_stats (ported), step 1: <pool>_all_stats.csv ---------------------
# One row per sample column (samples + empty wells + Unknown):
#   sample_name, __no_feature, __ambiguous, __too_low_aQual,
#   __not_aligned, __alignment_not_unique, total_mapped_reads
if [[ -s "$ALL_STATS" ]]; then
  log "SKIP all_stats (already done: $ALL_STATS)"
else
  log "RUN  report_stats all_stats"
  awk '
  BEGIN { FS = "\t"; OFS = "," }
  NR == 1 {
    n = NF - 2
    for (i = 3; i <= NF; i++) sample[i-2] = $i
    next
  }
  $1 !~ /^__/ {
    for (i = 3; i <= NF; i++) total[i-2] += $i
    next
  }
  /^__no_feature/           { for (i=3;i<=NF;i++) nf[i-2]   = $i; next }
  /^__ambiguous/            { for (i=3;i<=NF;i++) amb[i-2]  = $i; next }
  /^__too_low_aQual/        { for (i=3;i<=NF;i++) tl[i-2]   = $i; next }
  /^__not_aligned/          { for (i=3;i<=NF;i++) na[i-2]   = $i; next }
  /^__alignment_not_unique/ { for (i=3;i<=NF;i++) nu[i-2]   = $i; next }
  END {
    print "sample_name","__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique","total_mapped_reads"
    for (i = 1; i <= n; i++)
      print sample[i], nf[i]+0, amb[i]+0, tl[i]+0, na[i]+0, nu[i]+0, total[i]+0
  }
  ' "$FRC_RENAMED" > "$ALL_STATS"
  [[ -s "$ALL_STATS" ]] || die "report_stats did not produce ${ALL_STATS} (or it is empty). Check that $FRC_RENAMED looks right (a header plus per-gene rows), fix the cause, and re-run."
  log "  all_stats: $ALL_STATS"
fi

#--- report_stats (ported), step 2: <pool>_report_summary.csv ----------------
# Join each sample row with its plate well.  Well source (in order):
#   * real samples: the lib-prep provenance TSV (Sample_ID <TAB> ... <TAB> Well)
#   * empty wells  : the well is already embedded in 'empty_well_<well>'
#   * Unknown / no well available: NA  (same fallback as the legacy pipeline)
#
# In --barcode-tsv mode there is no plate or lib-prep provenance, so every
# sample gets well 'NA' — the legacy pipeline's 'well not found' fallback.
WLOOKUP="$(mktemp)"
if [[ -n "${LIBPREP_PROV:-}" && -s "${LIBPREP_PROV}" ]]; then
  awk '
  # First file: lib-prep provenance — columns: Sample_ID, Sample_raw, Well, Barcode
  FNR == NR {
    # skip header; keep rows whose Well column matches A-H + 2 digits
    if (FNR > 1) {
      well = $3
      if (well ~ /^[A-H][0-9][0-9]$/) well_of[$1] = well
    }
    next
  }
  # Second file: sample_barcode_map.tsv — barcode<TAB>name
  {
    if ($2 ~ /^empty_well_/) well_of[$2] = substr($2, 12)
    if ($2 == "Unknown")     well_of["Unknown"] = "NA"
  }
  END { for (n in well_of) printf "%s\t%s\n", n, well_of[n] }
  ' FS="\t" OFS="\t" "${LIBPREP_PROV}" "$SAMPLE_MAP" > "$WLOOKUP"
else
  awk 'BEGIN{FS="\t"} {print $2 "\tNA"}' "$SAMPLE_MAP" > "$WLOOKUP"
  echo -e "Unknown\tNA" >> "$WLOOKUP"
fi
[[ -s "$WLOOKUP" ]] || { rm -f "$WLOOKUP"; die "The well-lookup table for the report grid could not be built. Inputs checked: lib-prep provenance (\"${LIBPREP_PROV:-<none>}\") and sample map ($SAMPLE_MAP). See earlier messages for why neither was usable. (In --barcode-tsv mode with no plate there is no provenance — the pipeline still writes 'Unknown,NA' + all samples as 'NA'.)"; }

if [[ -s "$REPORT_CSV" ]]; then
  log "SKIP report_summary (already done: $REPORT_CSV)"
else
  log "RUN  report_stats report_summary"
  awk '
  BEGIN { FS = OFS = "," }
  NR == FNR { well_of[$1] = $2; next }
  FNR == 1 { next }
  {
    report_sample_name   = $1
    no_feature           = $2
    ambiguous            = $3
    too_low_aQual        = $4
    not_aligned          = $5
    alignment_not_unique = $6
    total_mapped         = $7

    depth = total_mapped + no_feature + ambiguous + too_low_aQual + not_aligned + alignment_not_unique

    depth_mio       = depth        / 1000000
    mapped_mio      = total_mapped / 1000000
    no_feature_mio  = no_feature   / 1000000
    not_aligned_mio = not_aligned  / 1000000

    report_well = (report_sample_name in well_of) ? well_of[report_sample_name] : "NA"

    printf "%s,%s,%.2f,%.2f,%.2f,%.2f\n", \
      report_well, report_sample_name, depth_mio, mapped_mio, no_feature_mio, not_aligned_mio
  }
  ' FS="	" "$WLOOKUP" FS="," "$ALL_STATS" | \
  sort -t',' -k1.2,1.3n -k1.1,1.1 | \
  awk 'BEGIN { print "plate_well,sample_name,depth_mio,mapped_reads_mio,no_feature_mio,not_aligned_mio" } { print }' \
    > "$REPORT_CSV"
  rm -f "$WLOOKUP"
  [[ -s "$REPORT_CSV" ]] || die "report_stats did not produce a usable ${REPORT_CSV}. See the messages printed above for the reason, fix it, and re-run."
  log "  report_csv: $REPORT_CSV"
fi

log ""
log "  FRC:        ${FRC_DIR}  (counts.detailed.txt, counts.txt, …)"
log "  all_stats:  $ALL_STATS"
log "  report_csv: $REPORT_CSV"

#=== matrix CSVs  [MANDATORY] ================================================
# Flatten the sample-named umiDedup matrices into single gene x sample CSVs.
# Row key: Gene_id (column 1, unique) then gene_symbol (column 2, annotation)
# — same layout as the old FastReadCounter deliveries (Gene_id, Gene_name).
# Cell values are written exactly as they are in the MTX.
# One column per Sample_ID from barcodes.tsv.
# This is the "ready for PCA" table — the equivalent of the legacy
# <pool>_counts.detailed_forPCA.csv — and the input for any downstream
# clustering/differential-expression step.
log ""
log "UMI matrix CSV (mandatory)"
run_or_skip "UMI matrix (deduplicated UMI counts)" "$UMI_MATRIX_CSV" "" "mandatory" \
  "$PYTHON" "$UMI_MATRIX_PREP" \
    --samples-dir "$RENAME_DIR" \
    --out-csv   "$UMI_MATRIX_CSV" \
    --dedup     "1MM_Directional"
[[ -s "$UMI_MATRIX_CSV" ]] || die "The final UMI-based gene count table (one row per gene, one column per sample — the primary result table) was not produced at $UMI_MATRIX_CSV. See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
log "UMI matrix CSV: $UMI_MATRIX_CSV"

log ""
log "Reads matrix CSV (mandatory — NoDedup / all PCR copies)"
run_or_skip "Reads matrix (NoDedup, all PCR copies)" "$READS_MATRIX_CSV" "" "mandatory" \
  "$PYTHON" "$UMI_MATRIX_PREP" \
    --samples-dir "$RENAME_DIR" \
    --out-csv   "$READS_MATRIX_CSV" \
    --dedup     "NoDedup"
[[ -s "$READS_MATRIX_CSV" ]] || die "The final gene count table based on total read counts rather than unique UMIs (one row per gene, one column per sample) was not produced at $READS_MATRIX_CSV. See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
log "Reads matrix CSV: $READS_MATRIX_CSV"

#=== Sample demultiplex (BAM, picard)  [OPTIONAL] ==========================
if [[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]]; then
  log ""; log "BAM demultiplex via picard (optional)"
  [[ -f "$BAM_BARCODES" ]] || die "This optional step (splitting the aligned results file into one file per sample) needs the per-sample barcode list ($BAM_BARCODES), which is normally produced by the earlier 'barcode prep' step. That earlier step apparently did not produce it in this run — see its messages in the run's pipeline.log for the reason. Fix that earlier step's problem and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
  [[ -s "$STAR_BAM"     ]] || die "This optional step (splitting the aligned results file into one file per sample) needs the aligned reads file produced by the STAR step ($STAR_BAM), but it is missing or empty. That points to a problem in the STAR alignment step itself — see its messages in the run's pipeline.log for the reason. Fix that and re-run (use --resume-dir pointing at this same folder so those steps are not redone)."
  done_bam="${DMUX_BAM_DIR}/_done"
  # Skip only if the step truly finished: per-sample BAMs present (a _done
  # marker from an OLD failed attempt does not count — see the FASTQ step
  # below for the same fix).
  n_bam=$(find "$DMUX_BAM_DIR" -maxdepth 1 -name 'unmatched.bam' -o -name '*.bam' 2>/dev/null | wc -l)
  if [[ "$n_bam" -gt 0 ]]; then
    [[ -e "$done_bam" ]] && log "SKIP BAM demux (per-sample BAMs already present: $DMUX_BAM_DIR)" \
                        || log "SKIP BAM demux (per-sample BAMs already present: $DMUX_BAM_DIR; restoring _done marker)"
    [[ -e "$done_bam" ]] || touch "$done_bam"
  else
    if [[ "$HAVE_PICARD" -eq 1 ]]; then
      # The vendor map (bam_barcodes.tsv) also lists empty_well_<WELL> rows —
      # and the vendor demultiplex_bam.sh is FROZEN (no skip flag, blindly
      # demuxes every row; each FilterSamReads is a full BAM scan).  So filter
      # the map HERE, at our call site: only real samples get a per-sample BAM.
      # (STAR whitelist + samples_all_wells/ keep the empty wells intact.)
      local_bam_map="${DMUX_BAM_DIR}/barcodes_samples.tsv"
      awk -F'\t' 'NR==1 && $1=="sample_id" {print; next} $1 !~ /^empty_well_/ {print}' "$BAM_BARCODES" > "$local_bam_map"
      ew_bam=$(( $(wc -l < "$BAM_BARCODES") - $(wc -l < "$local_bam_map") ))
      picard_log="${LOGS_DIR}/demux_bam.log"
      log "RUN  picard FilterSamReads per sample (tag CR)"
      log "  demux map: $local_bam_map  ($(wc -l < "$local_bam_map") line(s); ${ew_bam} empty_well row(s) excluded)"
      log "  step log -> ${picard_log}"
      bash "$DEMUX_BAM" "$DMUX_BAM_DIR" "$STAR_BAM" "$local_bam_map" "CR" >"$picard_log" 2>&1 \
        || { tail -20 "$picard_log" >&2 2>/dev/null || true; \
             warn "picard demux reported an error — see step log: ${picard_log}"; }
      touch "$done_bam"
      log "per-sample BAMs in: $DMUX_BAM_DIR"
    else
      warn "This optional step (splitting the aligned-results file into one file per sample) cannot run because 'picard' is not installed in the current environment. It ships inside the 'brbseq' Conda environment — activate it with 'conda activate brbseq' (create it first with 'conda env create -f brbseq-env.yaml' if you have never made it). The rest of the pipeline will continue as normal; only this per-sample split will be missing from this run's output."
    fi
  fi
else
  log ""; log "BAM demux — skipped (--skip-dmux-bam)"
fi

#=== Sample demultiplex (FASTQ, fqtk)  [OPTIONAL] ==========================
if [[ "$RUN_OPT_DEMUX_FQ" -eq 1 ]]; then
  log ""; log "FASTQ demultiplex via fqtk (optional)"
  # This step needs BOTH reads: fqtk demuxes off the 14-bp CB in R1 and
  # emits the cDNA read from R2. If R2 wasn't supplied, fail out clearly now.
  n_fq=$(find "$DMUX_FQ_DIR" -maxdepth 1 -name 'unmatched.fq.gz' -o -name '*.fq.gz' 2>/dev/null | wc -l)
  if [[ -z "$R2_GENOMIC" ]]; then
    warn "Skipping FASTQ demux: a cDNA read (R2) is required to demultiplex. Re-run with --r2 <R2.fastq.gz>."
  # Skip only if the step TRULY finished: per-sample FASTQs present (a _done
  # marker left behind by an OLD failed attempt does not count — that was the
  # very bug that left this run's fastq demux failed).
  elif [[ "$n_fq" -gt 0 ]]; then
    [[ -e "${DMUX_FQ_DIR}/_done" ]] && log "SKIP FASTQ demux (per-sample FASTQs already present: $DMUX_FQ_DIR)" \
                                   || log "SKIP FASTQ demux (per-sample FASTQs already present: $DMUX_FQ_DIR; restoring _done marker)"
    [[ -e "${DMUX_FQ_DIR}/_done" ]] || touch "${DMUX_FQ_DIR}/_done"
  else
    # A failed attempt (old code) may have left a marker behind: clear it so
    # this re-run is not confused with a completed one.
    [[ -e "${DMUX_FQ_DIR}/_done" ]] && { rm -f -- "${DMUX_FQ_DIR}/_done"; log "  (stale _done marker from an incomplete attempt removed: ${DMUX_FQ_DIR}/_done)"; }
    if [[ "$HAVE_FQTK" -eq 1 ]]; then
      # Same empty-well exclusion as the BAM demux above: the vendor map
      # (fastq_sample_metadata.tsv) lists empty_well_* rows; fqtk demuxes
      # every row, so filter the map at our call site (header row kept).
      local_fq_map="${DMUX_FQ_DIR}/fastq_sample_metadata.tsv"
      awk -F'\t' 'NR==1 && $1=="sample_id" {print; next} $1 !~ /^empty_well_/ {print}' "$FQ_META" > "$local_fq_map"
      ew_fq=$(( $(wc -l < "$FQ_META") - $(wc -l < "$local_fq_map") ))
      fqtk_log="${LOGS_DIR}/demux_fq.log"
      log "RUN  fqtk demux"
      log "  demux map: $local_fq_map  ($(wc -l < "$local_fq_map") line(s) incl. header; ${ew_fq} empty_well row(s) excluded)"
      log "  R1 barcode = 14 BP (CB) + 14 BP (UMI)   |   R2 = 90 BP (cDNA)"
      log "  step log -> ${fqtk_log}"
      # Reads are demultiplexed off the 14-bp CB in R1 (R1 = 14B + 14M; R2 = 90T cDNA),
      # per the Mercurius BRB-seq guide §1.8. Output: one cDNA FASTQ per sample (named
      # *_R1.fq.gz by fqtk, despite holding the R2/cDNA bases) + unmatched + demux-metrics.
      fqtk demux \
        -i "$R1_BARCODE" "$R2_GENOMIC" \
        -r "14B14M" "90T" \
        -s "$local_fq_map" \
        -o "$DMUX_FQ_DIR" \
        -u "unmatched" \
        -t "$THREADS" \
        >"$fqtk_log" 2>&1 \
        || { tail -20 "$fqtk_log" >&2 2>/dev/null || true; \
             warn "'fqtk' (the tool that splits the combined FASTQ into one file per sample) reported a problem. See step log: ${fqtk_log}. The per-sample FASTQ files it managed to write are still in $DMUX_FQ_DIR and can be used, but some samples may be missing or incomplete. If you expected every sample to be present, check that folder and compare against your sample list; if something is missing, you can re-run just this step (it is safe to re-run)."; }
      touch "${DMUX_FQ_DIR}/_done"
      log "per-sample FASTQs in: $DMUX_FQ_DIR"
    else
      warn "This optional step (splitting the raw sequence files into one per sample) cannot run because 'fqtk' is not installed in the current environment. It ships inside the 'brbseq' Conda environment — activate it with 'conda activate brbseq' (create it first with 'conda env create -f brbseq-env.yaml' if you have never made it). The rest of the pipeline will continue as normal; only this per-sample split will be missing from this run's output."
    fi
  fi
else
  log ""; log "FASTQ demux — skipped (--skip-dmux-fq)"
fi

#=== Run README (folder guide)  [ALWAYS] ===================================
# A short, user-readable index of what this folder holds — written at the end
# of a run so it always matches what is actually there.  Self-contained HTML
# (no external CSS/JS) so it opens in a browser on any machine by double-click.
# Pure documentation: it inspects finished artifacts (done markers / product
# files) and writes README.html; it changes no data and cannot affect results.
write_run_readme() {
  local README="${RUN_DIR}/README.html"
  local rel_bam rel_raw rel_smp rel_smpall rel_dups rel_umi rel_reads
  rel_bam="${STAR_BAM#"${RUN_DIR}"/}"
  rel_raw="${RAW_DIR#"${RUN_DIR}"/}"
  rel_smp="${RENAME_DIR#"${RUN_DIR}"/}"
  rel_smpall="${RENAME_ALL_DIR#"${RUN_DIR}"/}"
  rel_dups="${DUP_CSV#"${RUN_DIR}"/}"
  rel_umi="${UMI_MATRIX_CSV#"${RUN_DIR}"/}"
  rel_reads="${READS_MATRIX_CSV#"${RUN_DIR}"/}"

  # Presence probes: a step is listed only when its marker / product survived
  # on disk (covers both fresh runs and resumed ones).
  local _b="" _f="" _rseqc=0 _fastqc=0
  # NOTE: the demux _done markers are zero-byte (plain `touch`), so probe -e, not -s.
  if [[ "${RUN_OPT_DEMUX_BAM}" -eq 1 && -e "${DMUX_BAM_DIR}/_done" ]]; then _b="${DMUX_BAM_DIR#"${RUN_DIR}"/}"; fi
  if [[ "${RUN_OPT_DEMUX_FQ}"  -eq 1 && -e "${DMUX_FQ_DIR}/_done" ]];  then _f="${DMUX_FQ_DIR#"${RUN_DIR}"/}";  fi
  if [[ -s "${RSEQC_DIR}/${POOL}_read_distribution.txt" ]]; then _rseqc=1; fi
  if [[ -d "${FASTQC_DIR}" ]]; then _fastqc=1; fi
  if [[ -s "${REPORTS_DIR}/${POOL}_report_summary.csv" ]]; then _frc=1; fi

  {
    echo '<!DOCTYPE html>'
    echo '<html lang="en"><head><meta charset="utf-8">'
    echo '<meta name="viewport" content="width=device-width, initial-scale=1">'
    echo "<title>BRB-seq run — ${RUN_NAME} (${POOL})</title>"
    echo '<style>'
    echo '  body { font-family: system-ui, sans-serif; margin: 2rem auto; max-width: 52rem; line-height: 1.55; color: #1c1e21; padding: 0 1rem; }'
    echo '  h1 { font-size: 1.5rem; border-bottom: 2px solid #ddd; padding-bottom: .4rem; }'
    echo '  h2 { font-size: 1.15rem; margin-top: 1.6rem; }'
    echo '  code { background: #f2f3f5; padding: .1em .35em; border-radius: 4px; font-size: .92em; }'
    echo '  table { border-collapse: collapse; width: 100%; margin: .8rem 0; font-size: .95em; }'
    echo '  th, td { border: 1px solid #ddd; padding: .45rem .6rem; text-align: left; vertical-align: top; }'
    echo '  th { background: #f2f3f5; }'
    echo '  .meta { color: #555; font-size: .9em; margin-bottom: 1rem; }'
    echo '  ul { padding-left: 1.3rem; }'
    echo '  li { margin: .35rem 0; }'
    echo '  .primary { font-weight: 600; }'
    echo '</style></head><body>'
    echo "<h1>BRB-seq (MERCURIUS DRUG-seq) — run output folder</h1>"
    echo "<p class=\"meta\">Generated by brbseq_pipeline.sh v${SCRIPT_VERSION} · <code>${RUN_NAME}</code> · pool <code>${POOL}</code> · species <code>${SPECIES:-not specified}</code> · ${DATE_GEN:-$(date +%Y-%m-%d)}</p>"

    echo '<h2>What to use (the two deliverable tables)</h2>'
    echo '<table>'
    echo '<tr><th>File</th><th>Meaning</th><th>Use for</th></tr>'
    echo '<tr><td><code>'"${rel_umi}"'</code></td><td>Gene x sample UMI counts (deduplicated molecules)</td><td class="primary">PCA / differential expression (Seurat, edgeR) — the primary result table</td></tr>'
    echo '<tr><td><code>'"${rel_reads}"'</code></td><td>Gene x sample total read counts (all PCR copies)</td><td>Duplication context only — <strong>not</strong> a PCA/DE input</td></tr>'
    if [[ $_frc -eq 1 ]]; then
      echo '<tr><td><code>'"$(echo "${REPORTS_DIR}/${POOL}_report_summary.csv" | sed "s|${RUN_DIR}/||")"'</code></td><td>Per-well sequencing depth / mapped / no-feature / no-alignment read distribution (whole-plate grid + Unknown)</td><td class="primary">Sequencing QC &amp; CLC report — one row per plate well</td></tr>'
      echo '<tr><td><code>'"$(echo "${REPORTS_DIR}/${POOL}_all_stats.csv" | sed "s|${RUN_DIR}/||")"'</code></td><td>Same read-distribution counts, raw totals per sample (no mio rounding)</td><td>Underlying numbers behind the grid above</td></tr>'
    fi
    echo '</table>'
    echo '<p>Both tables have one row per gene (columns <code>Gene_id</code> then <code>gene_symbol</code>) and one column per sample (your <code>Sample_ID</code>). <code>empty_well_*</code> columns are unoccupied plate wells — drop them before downstream analysis.</p>'

    echo '<h2>Folder contents</h2>'
    echo '<ul>'
    if [[ -s "$UMI_MATRIX_CSV" ]]; then
      echo '<li class="primary"><code>'"${rel_umi}"'</code> — UMI count table (primary deliverable): gene x sample, deduplicated UMI (molecule) counts from STARsolo.</li>'
    fi
    if [[ -s "$READS_MATRIX_CSV" ]]; then
      echo '<li><code>'"${rel_reads}"'</code> — read count table (all PCR copies): gene x sample; duplication context only, not a PCA/DE input.</li>'
    fi
    if [[ -s "$DUP_CSV" ]]; then
      echo '<li><code>'"${rel_dups}"'</code> — per-sample UMI duplication statistics (UMI count vs read count, duplication rate).</li>'
    fi
    if [[ -s "$STAR_BAM" ]]; then
      echo '<li><code>'"${rel_bam}"'</code> — aligned, coordinate-sorted STAR BAM (all samples in one file; per-well barcode in the CB/CR tags) plus its <code>.bai</code> index.</li>'
    fi
    if [[ -d "$RAW_DIR" ]]; then
      echo '<li><code>'"${rel_raw}"'/</code> — STAR raw count tables (sparse MTX, one column per barcode) — the source the two CSV tables above are built from; kept untouched.</li>'
    fi
    if [[ -d "$RENAME_DIR" ]]; then
      echo '<li><code>'"${rel_smp}"'/</code> — the STAR count tables renamed by <code>Sample_ID</code> (real samples only).</li>'
    fi
    if [[ -d "$RENAME_ALL_DIR" ]]; then
      echo '<li><code>'"${rel_smpall}"'/</code> — the STAR count tables with all wells kept (real samples + <code>empty_well_*</code>), for debugging well assignment.</li>'
    fi
    if [[ -s "${BARCODES_DIR}/barcode_whitelist.txt" ]]; then
      echo '<li><code>barcodes/</code> — STAR barcode whitelist plus sample/barcode maps (<code>manifest*.tsv</code>, <code>fastq_sample_metadata.tsv</code>, <code>bam_barcodes.tsv</code>).</li>'
    fi
    if [[ $_fastqc -eq 1 ]]; then
      echo '<li><code>fastqc/</code> — FastQC quality reports for R1 and R2 (open the <code>*_fastqc.html</code> files in a browser).</li>'
    fi
    if [[ $_rseqc -eq 1 ]]; then
      echo '<li><code>rseqc/'"${POOL}"'_read_distribution.txt</code> — read-coverage distribution across exons (RSeQC).</li>'
    fi
    if [[ $_frc -eq 1 ]]; then
      echo '<li><code>files_for_report_and_CLC/'"${POOL}"'_report_summary.csv</code> — the whole-plate read-distribution grid (depth / mapped / no-feature / no-alignment per well, in millions — the CLC/sequencing report).</li>'
      echo '<li><code>files_for_report_and_CLC/'"${POOL}"'_all_stats.csv</code> — the same read-distribution counts as raw totals per sample (unrounded).</li>'
      echo '<li><code>fastreadcounter/</code> — FastReadCounter intermediates (<code>counts.detailed.txt</code> et al.) — the source the two report CSVs above are built from.</li>'
    fi
    if [[ -n "$_b" ]]; then
      echo '<li><code>'"${_b}"'/</code> — one aligned BAM per sample (split from the pool BAM; real samples only).</li>'
    fi
    if [[ -n "$_f" ]]; then
      echo '<li><code>'"${_f}"'/</code> — one demultiplexed cDNA FASTQ per sample (real samples only). Files are named <code>Sample_ID_R1.fq.gz</code> by fqtk convention but hold the 90-bp cDNA read; the 14-bp cell barcode and 14-bp UMI are carried in each read name. Plus <code>unmatched.R1.fq.gz</code> (reads with no matching barcode) and <code>demux-metrics.txt</code> (per-assignment stats).</li>'
    fi
    echo "<li><code>pipeline.log</code> — supervisor log: every step's status and exact command.</li>"
    echo '<li><code>logs/</code> — one full output log per pipeline step.</li>'
    echo '</ul>'

    echo '<h2>Notes</h2>'
    echo '<ul>'
    echo '<li><code>empty_well_*</code> columns are plate wells with no sample loaded — reads that landed there; keep them out of Seurat/edgeR analysis.</li>'
    echo '<li>The UMI table (<code>*umi_matrix.csv</code>) is the PCA/DE input; the reads table (<code>*reads_matrix.csv</code>) is duplication context only.</li>'
    echo '</ul>'
    echo '</body></html>'
  } > "$README"

  log "run README: $README  (written)"
}

write_run_readme

#=== Merging multi-lane FASTQ =============================================
log ""; log "multi-lane FASTQ merge not auto-run (it is optional;"
log "              pre-merge per-lane files or pass a comma list as --r1/--r2)"

log ""; log "══════ PIPELINE COMPLETE (v${SCRIPT_VERSION}) ══════"
log "  run directory         = ${RUN_DIR}"
log "  barcodes               = ${WHITELIST}  ($(wc -l < "$WHITELIST" 2>/dev/null || echo '?') barcodes)"
[[ -s "$STAR_BAM"     ]] && log "  STAR BAM             = ${STAR_BAM}"
[[ -d "$RAW_DIR"      ]] && log "  Gene/raw mtx         = ${RAW_DIR}"
[[ -s "$DUP_CSV"      ]] && log "  dup stats            = ${DUP_CSV}"
[[ -d "$RENAME_DIR"   ]] && log "  sample-named mtx     = ${RENAME_DIR}  (samples only)"
[[ -d "$RENAME_ALL_DIR" ]] && log "  sample-named mtx     = ${RENAME_ALL_DIR}  (all wells, debug)"
[[ -s "$UMI_MATRIX_CSV"  ]] && log "  UMI matrix CSV       = ${UMI_MATRIX_CSV}"
[[ -s "$READS_MATRIX_CSV" ]] && log "  Reads matrix CSV     = ${READS_MATRIX_CSV}"
[[ -s "$REPORT_CSV"       ]] && log "  report grid (CLC)    = $REPORT_CSV          (per-well depth/mapped/read distribution)"
[[ -s "$ALL_STATS"        ]] && log "  all_stats (report)   = $ALL_STATS          (raw per-sample read totals)"
[[ -s "${DMUX_DIR:-/nonexistent}"     ]] && :
[[ -d "${DMUX_BAM_DIR}/_done" && -s "${DMUX_BAM_DIR}/_done" ]] && log "  demux BAM              = ${DMUX_BAM_DIR}"
[[ -d "${DMUX_FQ_DIR}/_done"  && -s "${DMUX_FQ_DIR}/_done"  ]] && log "  demux FASTQ            = ${DMUX_FQ_DIR}"
[[ -s "${RUN_DIR}/README.html" ]] && log "  run README           = ${RUN_DIR}/README.html  (folder guide for users; open in a browser)"
log ""
log "Next step (outside this pipeline): downstream analysis (e.g. Seurat/Scanpy)"
log "clustering) on the produced matrices; this script stops at the"
log "pre-processing boundary."
exit 0
