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
#   * Every pipeline message AND the stdout/stderr of the tools it launches
#     (fastqc, samtools, picard, fqtk, python helpers) is mirrored to
#     <run-dir>/pipeline.log — no need to wrap the call in `tee`.
#   * The exact command line of each step is logged (run_or_skip prints
#     "RUN  <step>: <cmd>"; STAR's full invocation is logged the same way and
#     is also written as the first line of <run-dir>/star/STAR.log for
#     provenance). STAR's own runtime output stays in star/STAR.log.
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
#     (disable skipping with --resume force).
#   * The optional steps (FASTQ QC, RSeQC, BAM demux, FASTQ demux) are off by
#     default; a missing tool causes that step to log a warning and skip rather
#     than abort the run.
#=============================================================================
set -euo pipefail

#--- constants ---------------------------------------------------------------
SCRIPT_VERSION="1.0.0"
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
EMPTY_WELLS=0                    # pass --empty-wells to the parser (emit unoccupied plate wells)
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
RESUME="${RESUME:-auto}"         # auto | force
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
  --r1-barcode <fq.gz>  R1 barcode read (14 bp CB + 14 bp UMI)
  --r2-genomic <fq.gz>  R2 genomic / cDNA read
  --run-name <name>     run label for the output leaf
                        BRBseq_<run-name>_<timestamp>; in lib-prep mode
                        taken from the workbook 'Researcher' row (first
                        name) unless given here. Required in --barcode-tsv
                        mode (no workbook to read).

  MANIFEST SOURCE — exactly one of these two:
    (a)  --barcode-tsv <tsv>
    (b)  --libprep-xlsx <xlsx>   (+ optional --barcodes-plate <file>)
         Mercurius lib-prep workbook (sample code + Oligo-dT primer well). The
         barcode plate is joined by well. If --barcodes-plate is omitted, it is
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
  --threads <n>         threads (default 8)
  --resume <auto|force> auto = skip completed steps (default); force = rerun
  --skip-fastq          skip fastqc  (on by default)
  --skip-rseqc          skip read_distribution  (on by default)
  --skip-dmux-bam       skip picard BAM demux  (on by default)
  --skip-dmux-fq        skip fqtk FASTQ demux  (on by default)
  --empty-wells         workbook mode: also emit every unoccupied plate well as
                        an empty_well_<WELL> column (old-pipeline convention)
                        so all plate wells show up in the STARsolo matrix.
                        Reserved prefix: sample names must not start with
                        'empty_well'
  --dry-run             emit planned commands and exit (no side effects)
  -h, --help            this help
  --version             print version and exit

EXAMPLES
  # TSV manifest:
  ./brbseq_pipeline.sh --species mouse \
      --r1-barcode /data/POOL/S1_R1_001.fastq.gz \
      --r2-genomic /data/POOL/S1_R2_001.fastq.gz \
      --barcode-tsv /data/POOL/sample_barcode.tsv \
      --run-name Daniel --threads 16 --out-dir /mnt/data/BRB-seq_output

  # Lib-prep workbook instead of a TSV (plate inferred from the workbook):
  ./brbseq_pipeline.sh --species mouse \
      --r1-barcode /data/POOL/S1_R1_001.fastq.gz \
      --r2-genomic /data/POOL/S1_R2_001.fastq.gz \
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
    --r1-barcode)   R1_BARCODE="${2:?missing --r1-barcode value}";   shift 2;;
    --r2-genomic)   R2_GENOMIC="${2:?missing --r2-genomic value}";   shift 2;;
    --barcode-tsv)  BARCODE_TSV="${2:?missing --barcode-tsv value}"; shift 2;;
    --libprep-xlsx) LIBPREP_XLSX="${2:?missing --libprep-xlsx value}"; shift 2;;
    --barcodes-plate) BARCODES_PLATE="${2:?missing --barcodes-plate value}"; shift 2;;
    --empty-wells)    EMPTY_WELLS=1;                                           shift;;
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
    --resume)       RESUME="${2:?missing --resume value}";            shift 2;;
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
#   --libprep-xlsx  (+ optional --barcodes-plate).
# The plate is optional in workbook mode: parse_libprep_xlsx.py infers it from
# the workbook's 'Kit ref' row (barcodes_96_<family>_*.txt next to the
# workbook); an explicitly given plate is authoritative (a family mismatch with
# the Kit ref is a loud warning, not an abort).
if [[ -n "$BARCODE_TSV" && -n "$LIBPREP_XLSX" ]]; then
  die "You specified both --barcode-tsv and --libprep-xlsx, but only one of them can be used (both define the sample ↔ barcode mapping).  Remove one of the two flags and re-run."
fi
LIBPREP_MODE=0
if [[ -n "$LIBPREP_XLSX" ]]; then
  LIBPREP_MODE=1
fi
MISSING=()
[[ -n "$R1_BARCODE"  ]] || MISSING+=( "--r1-barcode" )
[[ -n "$R2_GENOMIC"  ]] || MISSING+=( "--r2-genomic" )
# manifest source: exactly one of (a) --barcode-tsv  or  (b) --libprep-xlsx.
if [[ -z "$LIBPREP_XLSX" && -n "$BARCODES_PLATE" ]]; then
  die "--barcodes-plate only makes sense together with the lib-prep workbook, but --libprep-xlsx was not given. Either add --libprep-xlsx <workbook.xlsx>, or remove --barcodes-plate and let the workbook's 'Kit ref' row choose the plate automatically."
elif [[ "$LIBPREP_MODE" -ne 1 && -z "$BARCODE_TSV" ]]; then
  MISSING+=( "manifest source: --barcode-tsv <tsv>   or   --libprep-xlsx <xlsx>  (+ optional --barcodes-plate <file>)" )
fi
[[ -n "$SPECIES"     ]] || MISSING+=( "--species" )
[[ -n "$RUN_NAME" || -n "$AVT_ID" || "$LIBPREP_MODE" -eq 1 ]]   || MISSING+=( "--run-name <name>" )
[[ "$RESUME" == "auto" || "$RESUME" == "force" ]] || die "--resume can only be 'auto' (skip steps that already finished) or 'force' (run every step again). Use one of those two and re-run."
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
[[ -f "$R1_BARCODE"  ]] || die "The R1 (barcode-channel) FASTQ file was not found: $R1_BARCODE.  Check that the path in --r1-barcode is correct and the file is still there, then re-run."
[[ -f "$R2_GENOMIC"  ]] || die "The R2 (genomic-channel) FASTQ file was not found: $R2_GENOMIC.  Check that the path in --r2-genomic is correct and the file is still there, then re-run."
# In lib-prep mode the manifest is generated later; only validate the workbook
# here.  The plate is optional: when omitted, the
# parser infers it from the workbook's 'Kit ref' row.
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  [[ -f "$LIBPREP_XLSX" ]] || die "The lib-prep workbook (--libprep-xlsx) file could not be found: $LIBPREP_XLSX. Check the path and re-run."
  [[ -z "$BARCODES_PLATE" || -f "$BARCODES_PLATE" ]] || die "The barcode plate file (--barcodes-plate) could not be found: $BARCODES_PLATE. Check the path and re-run."
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
TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_DIR%/}/BRBseq_${RUN_NAME}_${TS}"
BARCODES_DIR="${RUN_DIR}/barcodes"
FASTQC_DIR="${RUN_DIR}/fastqc"
STAR_DIR="${RUN_DIR}/star"
RSEQC_DIR="${RUN_DIR}/rseqc"
DMUX_BAM_DIR="${RUN_DIR}/demux_bam"
DMUX_FQ_DIR="${RUN_DIR}/demux_fastq"
RENAME_DIR="${STAR_DIR}/samples"
DUP_CSV="${RUN_DIR}/duplication_statistics.csv"
UMI_MATRIX_CSV="${RUN_DIR}/${POOL}_umi_matrix.csv"
READS_MATRIX_CSV="${RUN_DIR}/${POOL}_reads_matrix.csv"
# In --barcode-tsv mode the manifest lives where the user put it (no generation);
# in lib-prep mode the manifest is generated into <run>/barcodes/ and survives re-runs.
MANIFEST="$BARCODE_TSV"
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
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

#--- tool discovery ----------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }
# STAR ships in the brbseq conda env (bioconda); autodetect the binary.
STAR_BIN=""
for cand in STAR-avx2 STAR-avx STAR; do
  if have "$cand"; then STAR_BIN="$cand"; break; fi
done
HAVE_STAR=0;      [[ -n "$STAR_BIN" ]] && have "$STAR_BIN" && HAVE_STAR=1
HAVE_FASTQC=0;    have fastqc              && HAVE_FASTQC=1
HAVE_SAMTOOLS=0;  have samtools            && HAVE_SAMTOOLS=1
HAVE_RSEQC=0;     have read_distribution.py && HAVE_RSEQC=1
HAVE_PICARD=0;    have picard              && HAVE_PICARD=1
HAVE_FQTK=0;      have fqtk                && HAVE_FQTK=1
have "$PYTHON"    || PYTHON="$(have python3 && echo python3 || echo python)"
#--- auto-activate the known production env (ptegsrv3) ----------------------
# The brbseq env on ptegsrv3 lives at this LITERAL prefix — the square
# brackets are real characters in the path, so it must stay quoted (an
# unquoted '[' would glob). Override for other hosts via the environment:
#   BRBSEQ_ENV=<prefix>/envs/brbseq ./brbseq_pipeline.sh ...
BRBSEQ_ENV="${BRBSEQ_ENV:-/home/gecf/softwares/anaconda/[/home/gecf/softwares/anaconda]/envs/brbseq}"
if [[ "$HAVE_STAR" -eq 0 && -d "${BRBSEQ_ENV}/bin" ]]; then
  export PATH="${BRBSEQ_ENV}/bin:${PATH}"
  [[ -x "${BRBSEQ_ENV}/bin/python" ]] && PYTHON="${BRBSEQ_ENV}/bin/python"
  # PATH changed: re-run tool discovery so HAVE_* reflects the env.
  STAR_BIN=""
  for cand in STAR-avx2 STAR-avx STAR; do
    if have "$cand"; then STAR_BIN="$cand"; break; fi
  done
  HAVE_STAR=0;      [[ -n "$STAR_BIN" ]] && have "$STAR_BIN" && HAVE_STAR=1
  HAVE_FASTQC=0;    have fastqc              && HAVE_FASTQC=1
  HAVE_SAMTOOLS=0;  have samtools            && HAVE_SAMTOOLS=1
  HAVE_RSEQC=0;     have read_distribution.py && HAVE_RSEQC=1
  HAVE_PICARD=0;    have picard              && HAVE_PICARD=1
  HAVE_FQTK=0;      have fqtk                && HAVE_FQTK=1
  log "env      = ${BRBSEQ_ENV}"
fi
if [[ "$HAVE_STAR" -eq 0 ]]; then
  die "The Python environment that runs this pipeline cannot be found. This is the set of programs (STAR, samtools, etc.) the pipeline needs. Two ways to fix it: 1) create it once, from this script's folder: conda env create -f ${SCRIPT_DIR}/brbseq-env.yaml, then activate it: conda activate brbseq — then re-run;  2) or, if you already have such an environment, point to it:  BRBSEQ_ENV=<folder>/envs/<name> brbseq_pipeline.sh ..."
fi

#--- python with openpyxl (workbook mode) -----------------------------------
# parse_libprep_xlsx.py needs 'import openpyxl'.  The caller activates their
# own environment (e.g. 'conda activate brbseq') before launching; this script
# only verifies the active python has it and fails with a clear message if not.
py_has_openpyxl() { [[ "$1" != "" ]] && "$1" -c "import openpyxl" >/dev/null 2>&1; }
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  py_has_openpyxl "$PYTHON" || {
    have python3 || PYTHON="python3"
  }
  py_has_openpyxl "$PYTHON" || die "Reading the lib-prep Excel workbook needs a Python setup that includes the 'openpyxl' package, but the Python currently in use does not have it. Usually this is fixed by activating the 'brbseq' Conda environment:  conda activate brbseq  — then re-run the pipeline. If you have never created that environment, create it once from the 'brbseq-env.yaml' file that sits next to this pipeline script:  conda env create -f brbseq-env.yaml  and then 'conda activate brbseq', and re-run."
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
  if [[ "$EMPTY_WELLS" -eq 1 ]]; then
    parse_args+=( --empty-wells )
  fi
  # Capture the parser's exit code explicitly (a bare "cmd; rc=$?" is killed
  # by set -e before the assignment; "if ! cmd" hides $? behind the inverted
  # status).  "cmd || rc=$?" is the set-e-safe idiom.
  parse_rc=0
  "$PYTHON" "$PARSE_LIBPREP" "${parse_args[@]}" || parse_rc=$?
  if [[ "$parse_rc" -ne 0 ]]; then
    die "The step that reads the lib-prep Excel workbook and builds the sample→barcode list failed (error code ${parse_rc}; the specific reason is printed in the lines above this message — start there). In this 'workbook' mode the pipeline normally finds the matching barcode plate file automatically using the 'Kit ref' row of the workbook, and expects a file like  barcodes_96_<family>_*.txt  saved in the same folder as the workbook. If that file is missing or named differently, either rename/move it into place, or pass it explicitly with  --barcodes-plate <file>. Then re-run."
  fi
  [[ -s "$MANIFEST" ]] || die "The sample→barcode list that should have been produced from the Excel workbook ended up completely empty (file: $MANIFEST) — meaning none of the samples in the workbook were successfully matched to a barcode. See the messages printed above this one for the specific mismatches (e.g. a well position in the workbook that does not exist on the plate file). Fix those in the workbook or the plate file, then delete this same output folder (or the stale file at $MANIFEST) and re-run so the list is generated afresh."
  log "manifest   = ${MANIFEST}  (generated from lib-prep workbook)"
elif [[ "$LIBPREP_MODE" -eq 1 && ! -s "$MANIFEST" && "$DRY_RUN" -eq 1 ]]; then
  log "manifest   = (workbook mode; would generate from lib-prep workbook — skipped in dry-run)"
elif [[ "$LIBPREP_MODE" -eq 1 ]]; then
  # manifest already on disk (resume).  It must match this call:
  #  - flag set, rows missing  -> stale artifact (resume created without the
  #    flag).  FAIL rather than regenerate in place: downstream steps of
  #    --resume may already depend on the existing whitelist.  Re-run with
  #    --resume force (fresh manifest + full re-run) or start a new run.
  #  - flag unset, rows present -> same reason, reverse direction.
  ew_in_manifest=0
  grep -q $'^empty_well_[A-H][0-9][0-9]	' "$MANIFEST" && ew_in_manifest=1
  if [[ "$EMPTY_WELLS" -eq 1 && "$ew_in_manifest" -eq 0 ]]; then
    die "You asked to include empty wells (--empty-wells), but the sample→barcode list already saved from a previous run of this experiment (file: $MANIFEST) does not include any empty-well entries — it was generated earlier without that option. Mixing the two would give an inconsistent list, so the pipeline is stopping rather than silently proceeding. To fix: delete the stale file (or the whole previous run's output folder) and re-run, and it will regenerate the list including the empty wells."
  fi
  if [[ "$EMPTY_WELLS" -eq 0 && "$ew_in_manifest" -eq 1 ]]; then
    die "The existing sample→barcode list from a previous run of this experiment (file: $MANIFEST) does include empty-well entries, but this call does not ask for them (no --empty-wells flag). Mixing the two would give an inconsistent list, so the pipeline is stopping rather than silently proceeding. To fix: either add --empty-wells to keep the list as-is, or delete the stale file (or the whole previous run's output folder) and re-run without that flag, and the list will be regenerated without empty wells."
  fi
  log "manifest   = ${MANIFEST}  (already present, reusing)"
fi

#--- banner -------------------------------------------------------------------
log "==== MERCURIUS DRUG-seq pipeline v${SCRIPT_VERSION} (resume=${RESUME} dry_run=${DRY_RUN}) ===="
log "run name   = ${RUN_NAME}"
log "pool       = ${POOL}"
log "threads    = ${THREADS}"
log "R1 (CB+UMI)= ${R1_BARCODE}"
log "R2 (cDNA)  = ${R2_GENOMIC}"
log "manifest   = ${MANIFEST}  (header-mode=${MANIFEST_HEADER})"
log "STAR index = ${STAR_INDEX}   (star=${STAR_BIN:-NOT_FOUND})"
log "output     = ${RUN_DIR}"
log "optional:  fastqc=${RUN_OPT_FASTQ}  rseqc=${RUN_OPT_RSEQC}  demux_bam=${RUN_OPT_DEMUX_BAM}  demux_fq=${RUN_OPT_DEMUX_FQ}"
log "tools:     STAR=${HAVE_STAR}  fastqc=${HAVE_FASTQC}  samtools=${HAVE_SAMTOOLS}  rseqc=${HAVE_RSEQC}  picard=${HAVE_PICARD}  fqtk=${HAVE_FQTK}"
log "python     = ${PYTHON}"
if [[ -n "$LIBPREP_PROV" ]]; then
  log "lib-prep   = ${LIBPREP_XLSX}"
  if [[ -n "$BARCODES_PLATE" ]]; then
    log "plate      = ${BARCODES_PLATE}"
  else
    log "plate      = (not given — inferred from the workbook's 'Kit ref' row)"
  fi
  [[ "$EMPTY_WELLS" -eq 1 ]] && log "empty wells = ON (unoccupied plate wells emitted as empty_well_<WELL>)"
fi

#--- execute (or dry-run) ----------------------------------------------------
if [[ "$DRY_RUN" -eq 1 ]]; then
  log "---- DRY RUN: planned steps ----"
  hdr_flag=""
  [[ "$MANIFEST_HEADER" -eq 1 ]] && hdr_flag=" -H"
  log "  prepare_barcodes.sh -i <manifest> -w <OUT>/barcodes/barcode_whitelist.txt -f <OUT>/barcodes/fastq_sample_metadata.tsv -b <OUT>/barcodes/bam_barcodes.tsv${hdr_flag}"
  if [[ "$RUN_OPT_FASTQ" -eq 1 ]]; then
    log "  fastqc  --threads 16 --outdir <OUT>/fastqc <R1> <R2>"
  fi
  log "  STAR    (one STARsolo call, both dedup methods in list form)"
  log "      ${STAR_BIN} --runMode alignReads --runThreadN ${THREADS} --genomeDir ${STAR_INDEX} --readFilesIn <R2> <R1> --soloCBwhitelist <OUT>/barcodes/barcode_whitelist.txt --soloCBlen 14 --soloUMIlen 14 --soloUMIstart 15 --soloUMIdedup NoDedup 1MM_Directional --outSAMattributes '${SAM_ATTRS}' --quantMode GeneCounts ..."
  log "      consolidate  <prefix>Solo.out/Gene/raw/umiDedup-<method>.mtx -> <OUT>/star/Gene/raw/"
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
  log "  rename  ${PYTHON} <UTILS>/rename_matrix.py --raw_dir <OUT>/star/${POOL}_Solo.out/Gene/raw --mapping <manifest> --out-dir <OUT>/star/samples -> barcodes named by Sample_ID (raw barcode matrix kept untouched in Gene/raw)"
  if [[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]]; then
    log "  demuxBAM picard FilterSamReads per sample in <OUT>/barcodes/bam_barcodes.tsv (CR tag) -> <OUT>/demux_bam/<sample>.bam"
  fi
  if [[ "$RUN_OPT_DEMUX_FQ" -eq 1 ]]; then
    log "  demuxFQ fqtk demux <R1> <R2> --samples <OUT>/barcodes/fastq_sample_metadata.tsv --output <OUT>/demux_fastq"
  fi
  log "---- DRY RUN: no files modified ----"
  exit 0
fi

# --- actually mkdir (skipped in dry-run) -----------------------------------
mkdir -p "$BARCODES_DIR" "$STAR_DIR"
[[ "$RUN_OPT_FASTQ" -eq 1 ]] && mkdir -p "$FASTQC_DIR"
[[ "$RUN_OPT_RSEQC" -eq 1 ]] && mkdir -p "$RSEQC_DIR"
[[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]] && mkdir -p "$DMUX_BAM_DIR"
[[ "$RUN_OPT_DEMUX_FQ" -eq 1 ]] && mkdir -p "$DMUX_FQ_DIR"

#--- automatic pipeline.log ----------------------------------------------------
# Mirror everything (pipeline messages + tool stdout/stderr) to
# <RUN_DIR>/pipeline.log in addition to the terminal.  The per-step exact
# commands are already logged by run_or_skip ("RUN  <step>: <cmd>") and by
# star_solo (into star/STAR.log); tool chatter now lands in the file too.
LOG_FILE="${RUN_DIR}/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2> >(tee -a "$LOG_FILE" >&2)
# Drain the tee pipes before the shell dies, or the final buffered lines can
# be lost on exit (process-substitution pipe flush is async).
trap 'sleep 0.5' EXIT
log "log file = ${LOG_FILE}  (automatic; exact commands + tool output mirrored here)"

# --- run_or_skip -------------------------------------------------------------
# Execute the command ($5...) unless the marker file ($2) already exists in
# auto mode.  Failure policy: a failing OPTIONAL step logs a warning and does
# not abort the run; a failing MANDATORY step aborts the pipeline.
#   $1=step name  $2=marker file (empty => always run)
#   $3=tool name  (empty => no tool gate, e.g. a Python helper)
#   $4=mandatory|optional  $5...=command
run_or_skip() {
  local name="$1"; shift
  local marker="${1:-}"; shift
  local tool="${1:-}"; shift
  local mode="${1:-mandatory}"; shift
  local cmd=( "$@" )
  if [[ -n "$marker" && "$RESUME" == "auto" && -s "$marker" ]]; then
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
  log "RUN  ${name}: ${cmd[*]}"
  local rc=0
  "${cmd[@]}" || rc=$?
  if [[ "$rc" -ne 0 ]]; then
    if [[ "$mode" == "mandatory" ]]; then
      die "The required step '${name}' ran into an error and could not be completed, so the pipeline has stopped. The specific reason is in the messages printed in the lines just above this one (also saved in the run's pipeline.log). Fix that underlying problem, then re-run the pipeline — you can use --resume auto so steps that already finished are not redone."
    else
      warn "The optional step '${name}' ran into an error, so it did not finish — but this does not stop the pipeline; everything else will still run. The likely reason is in the messages printed just above this one (and in the run's pipeline.log). If you want this step's output, fix the underlying problem and re-run (it is safe to re-run just this one)."
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
  run_or_skip "fastqc" "$FQ_MARKER" "fastqc" "optional" \
    fastqc --threads 16 --outdir "$FASTQC_DIR" "$R1_BARCODE" "$R2_GENOMIC"
  touch "$FQ_MARKER"  # idempotent completion marker
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

# The suffixes appended to --outFileNamePrefix differ across STAR builds:
#   * older Alithea STAR : <prefix>_Aligned...  /  <prefix>Solo.out
#   * STAR 2.7.11b       : <prefix>Aligned...   /  <prefix>SOLOout
# A bare prefix (no trailing underscore) is therefore used, and the resulting
# BAM and SOLO directories are located by glob after the run.
#
# Both deduplication methods are requested in a single invocation via the list
# form of --soloUMIdedup (NoDedup 1MM_Directional). One alignment pass yields a
# per-method umiDedup-<method>.mtx under <prefix>SOLOout/Gene/raw/, matching the
# expected output and avoiding a second, nondeterministic alignment of the same reads.
STAR_BASE="${STAR_DIR}/${POOL}"

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
    log "The STAR alignment/counting step (STARsolo, UMI counting methods: ${DEDUP_METHODS[*]}) reported that it failed. The last lines of its own log file are printed directly above this line — read them for the specific reason (for example too few reads, an index it cannot read, or a run out of memory). Fix that, then re-run the pipeline with --resume auto so the steps that already finished are not redone."
    tail -20 "${log_file}" >&2 || true
    return 1
  fi
}

# Consolidate the single STAR run's SOLOout/Gene/raw (already per-method
# named) into the canonical ${STAR_DIR}/Gene/raw layout the dup-stats step
# expects:  bars.tsv, features.tsv, umiDedup-NoDedup.mtx,
# umiDedup-1MM_Directional.mtx.  Symlink when possible, copy as fallback.
consolidate_solo_output() {
  local dest="${STAR_DIR}/Gene/raw"
  mkdir -p "$dest"
  # consolidate_solo_output(): STAR appends a fixed suffix to --outFileNamePrefix
  # that is NOT stable across builds
  # (we have seen <prefix>Solo.out, <prefix>_Solo.out, <prefix>SOLOout).
  # We glob for the first dir under STAR_DIR matching any of those suffixes
  # rather than hard-coding one form.
  local solo_dir
  solo_dir="$(find "$STAR_DIR" -maxdepth 2 -type d \( -iname 'runSolo.out' -o -iname '*Solo.out' -o -iname '*SOLOout' \) 2>/dev/null | head -1)"
  [[ -n "$solo_dir" && -d "$solo_dir" ]] || die "STAR finished, but the pipeline did not find its output folder (normally a folder called 'Solo.out' or similar inside $STAR_DIR/). This usually means STAR itself stopped early with an error — check the STAR log files inside $STAR_DIR for the actual reason, then re-run the pipeline (set --resume force to redo every step)."
  log "consolidation: solo_dir=$solo_dir"
  local method kind
  for method in "${DEDUP_METHODS[@]}"; do
    local src
    for src in "${solo_dir}/Gene/raw/umiDedup-${method}.mtx" \
               "${solo_dir}/Gene/umiDedup-${method}.mtx"; do
      [[ -f "$src" ]] && { ln -f "$src" "${dest}/umiDedup-${method}.mtx" 2>/dev/null || cp -f "$src" "${dest}/umiDedup-${method}.mtx"; break; }
    done
    [[ -f "${dest}/umiDedup-${method}.mtx" ]] || die "STAR finished, but the raw UMI count matrix for the '${method}' counting mode was not found inside $solo_dir. This usually means STAR stopped early with an error — check the STAR log files inside $STAR_DIR for the actual reason, then re-run the pipeline (set --resume force to redo every step)."
    log "consolidate: $src  ->  ${dest}/umiDedup-${method}.mtx"
  done
  for kind in barcodes.tsv features.tsv; do
    for src in "${solo_dir}/Gene/raw/${kind}" "${solo_dir}/Gene/${kind}"; do
      [[ -f "$src" ]] && { ln -f "$src" "${dest}/${kind}" 2>/dev/null || cp -f "$src" "${dest}/${kind}"; break; }
    done
  done
  log "consolidated both matrices into $dest"
}

# Helper: find STAR's actual output BAM (any prefix+Aligned+sortedByCoord layout)
find_star_bam() {
  local d="$1"
  local f
  for f in "$(find "$d" -maxdepth 2 -type f -name '*Aligned.sortedByCoord.out.bam' 2>/dev/null)"; do
    printf '%s\n' "$f"; return 0
  done
  return 1
}

S14_MARKER="${STAR_DIR}/_s14_done"
S14_DONE=0
if [[ -f "$S14_MARKER" && -s "$S14_MARKER" ]]; then S14_DONE=1; fi
STAR_BAM=""

if [[ "$RESUME" == "auto" && "$S14_DONE" -eq 1 ]]; then
  log "SKIP STARsolo (already done per marker: $S14_MARKER)"
  STAR_BAM="$(find_star_bam "$STAR_DIR" | head -n 1)"
fi
if [[ "$RESUME" != "auto" || "$S14_DONE" -eq 0 || -z "$STAR_BAM" ]]; then
  # One STAR invocation, --soloUMIdedup in list form (both methods).
  star_solo "${STAR_BASE}"
  consolidate_solo_output
  STAR_BAM="$(find_star_bam "$STAR_DIR" | head -n 1)"
  [[ -n "$STAR_BAM" && -s "$STAR_BAM" ]] || die "The alignment step (STAR) ran to completion, but the aligned reads file it should have produced (Aligned.sortedByCoord.out.bam) was not found or is empty under $STAR_DIR. This usually means STAR itself stopped early with an error — see the STAR log files inside $STAR_DIR (e.g. STAR.log.out, STAR.log) for the specific reason. Fix that, then re-run the pipeline (use --resume auto so earlier steps are not redone)."
  echo "$(date +%FT%T)" > "${STAR_DIR}/_s14_done"
fi
log "STAR BAM: $STAR_BAM"
log "consolidated mtx dir: ${STAR_DIR}/Gene/raw"
SOLO_OUT="${STAR_DIR}/Gene/raw"   # alias used by the dup-stats step

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
    read_distribution.py -i "$STAR_BAM" -r "$BED_FILE" > "$RSEQC_OUT"
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
RAW_DIR="${STAR_DIR}/Gene/raw"
for req in "$RAW_DIR/barcodes.tsv" \
           "$RAW_DIR/umiDedup-NoDedup.mtx" \
           "$RAW_DIR/umiDedup-1MM_Directional.mtx"; do
  [[ -s "$req" ]] || die "The step that measures how duplicated / how 'empty' each sample's reads are needs three input files (the raw UMI count tables STAR produces), and one of them is missing or zero-sized: $req. All three must exist in $RAW_DIR before this step can run — this points to a problem earlier in the pipeline (the STAR / UMI-counting step). See the messages printed above for the likely reason; fix that, then re-run (use --resume auto so steps that already finished are not redone)."
done
run_or_skip "duplication stats" "$DUP_CSV" "" "mandatory" \
  "$PYTHON" "$DUP_STATS" \
    --raw_dir "$RAW_DIR" \
    --dedup_method "1MM_Directional" \
    --manifest "$MANIFEST" \
    --output_csv "$DUP_CSV"
[[ -s "$DUP_CSV" ]] || die "The step that measures how duplicated / how 'empty' each sample's reads are did not produce its result file ($DUP_CSV). See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume auto so steps that already finished are not redone)."
log "duplication stats: $DUP_CSV"
if [[ -f "${DUP_CSV%.*}.summary.csv" ]]; then
  log "summary           : ${DUP_CSV%.*}.summary.csv"
fi

#=== Sample-named count matrix  [MANDATORY] =================================
# Rename the barcode-labelled matrix columns (barcodes.tsv in Gene/raw) to the
# Sample_IDs from the run manifest, writing copies to ${STAR_DIR}/samples/.
# The raw barcode-labelled matrix in Gene/raw is left untouched.
log ""
log "sample-named count matrix (mandatory)"
run_or_skip "rename matrix" "${RENAME_DIR}/barcodes.tsv" "" "mandatory" \
  "$PYTHON" "$RENAME_MATRIX" \
    --raw-dir   "$RAW_DIR" \
    --mapping   "$MANIFEST" \
    --out-dir   "$RENAME_DIR" \
    --dedup-methods "${DEDUP_METHODS[*]}"
[[ -s "${RENAME_DIR}/barcodes.tsv" ]] || die "The step that gives every sample its proper name (matching the Sample_ID from your list) did not finish: the renamed barcode list ($RENAME_DIR/barcodes.tsv) is missing or empty. See the message(s) printed just above for the specific reason — often a sample in your list does not match any barcode found in the data, or vice versa. Fix the mismatch in your sample list and re-run (use --resume auto so steps that already finished are not redone)."
[[ -s "${RENAME_DIR}/umiDedup-1MM_Directional.mtx" ]] || die "The step that gives every sample its proper name (matching the Sample_ID from your list) did not finish: the renamed UMI count table ($RENAME_DIR/umiDedup-1MM_Directional.mtx) is missing or empty. See the message(s) printed just above for the specific reason — often a sample in your list does not match any barcode found in the data, or vice versa. Fix the mismatch in your sample list and re-run (use --resume auto so steps that already finished are not redone)."
log "sample-named matrix dir: ${RENAME_DIR}"
log "  (barcode->sample map: ${RENAME_DIR}/sample_barcode_map.tsv; raw matrix unchanged in ${RAW_DIR})"

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
[[ -s "$UMI_MATRIX_CSV" ]] || die "The final UMI-based gene count table (one row per gene, one column per sample — the primary result table) was not produced at $UMI_MATRIX_CSV. See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume auto so steps that already finished are not redone)."
log "UMI matrix CSV: $UMI_MATRIX_CSV"

log ""
log "Reads matrix CSV (mandatory — NoDedup / all PCR copies)"
run_or_skip "Reads matrix (NoDedup, all PCR copies)" "$READS_MATRIX_CSV" "" "mandatory" \
  "$PYTHON" "$UMI_MATRIX_PREP" \
    --samples-dir "$RENAME_DIR" \
    --out-csv   "$READS_MATRIX_CSV" \
    --dedup     "NoDedup"
[[ -s "$READS_MATRIX_CSV" ]] || die "The final gene count table based on total read counts rather than unique UMIs (one row per gene, one column per sample) was not produced at $READS_MATRIX_CSV. See the message(s) printed just above for the specific reason, fix it, and re-run (use --resume auto so steps that already finished are not redone)."
log "Reads matrix CSV: $READS_MATRIX_CSV"

#=== Sample demultiplex (BAM, picard)  [OPTIONAL] ==========================
if [[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]]; then
  log ""; log "BAM demultiplex via picard (optional)"
  [[ -f "$BAM_BARCODES" ]] || die "This optional step (splitting the aligned results file into one file per sample) needs the per-sample barcode list ($BAM_BARCODES), which is normally produced by the earlier 'barcode prep' step. That earlier step apparently did not produce it in this run — see its messages in the run's pipeline.log for the reason. Fix that earlier step's problem and re-run (use --resume auto so steps that already finished are not redone)."
  [[ -s "$STAR_BAM"     ]] || die "This optional step (splitting the aligned results file into one file per sample) needs the aligned reads file produced by the STAR step ($STAR_BAM), but it is missing or empty. That points to a problem in the STAR alignment step itself — see its messages in the run's pipeline.log for the reason. Fix that and re-run (use --resume auto so steps that already finished are not redone)."
  done_bam="${DMUX_BAM_DIR}/_done"
  if [[ "$RESUME" == "auto" && -s "$done_bam" ]]; then
    log "SKIP BAM demux (done marker present: $done_bam)"
  else
    if [[ "$HAVE_PICARD" -eq 1 ]]; then
      log "RUN  picard FilterSamReads per sample (tag CR)"
      bash "$DEMUX_BAM" "$DMUX_BAM_DIR" "$STAR_BAM" "$BAM_BARCODES" "CR"
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
  DONE_FQ="${DMUX_FQ_DIR}/_done"
  if [[ "$RESUME" == "auto" && -s "$DONE_FQ" ]]; then
    log "SKIP FASTQ demux (done marker present: $DONE_FQ)"
  else
    if [[ "$HAVE_FQTK" -eq 1 ]]; then
      log "RUN  fqtk demux"
      fqtk demux \
        -i "$R1_BARCODE" \
        -r "14I14I" \
        -s "$FQ_META" \
        -o "$DMUX_FQ_DIR" \
        -u "unmatched" \
        -t "$THREADS" \
        || warn "'fqtk' (the tool that splits the combined FASTQ into one file per sample) finished, but it reported a problem along the way. In plain terms: the per-sample FASTQ files it managed to write are still in $DMUX_FQ_DIR and can be used, but some samples may be missing or incomplete. If you expected every sample to be present, check that folder and compare against your sample list; if something is missing, you can re-run just this step (it is safe to re-run)."
      touch "$DONE_FQ"
      log "per-sample FASTQs in: $DMUX_FQ_DIR"
    else
      warn "This optional step (splitting the raw sequence files into one per sample) cannot run because 'fqtk' is not installed in the current environment. It ships inside the 'brbseq' Conda environment — activate it with 'conda activate brbseq' (create it first with 'conda env create -f brbseq-env.yaml' if you have never made it). The rest of the pipeline will continue as normal; only this per-sample split will be missing from this run's output."
    fi
  fi
else
  log ""; log "FASTQ demux — skipped (--skip-dmux-fq)"
fi

#=== Merging multi-lane FASTQ ==============================================
log ""; log "multi-lane FASTQ merge not auto-run (it is optional;"
log "              pre-merge per-lane files or pass a comma list as --r1-barcode/--r2-genomic)"

#=== summary =================================================================
log ""; log "══════ PIPELINE COMPLETE (v${SCRIPT_VERSION}) ══════"
log "  run directory         = ${RUN_DIR}"
log "  barcodes               = ${WHITELIST}  ($(wc -l < "$WHITELIST" 2>/dev/null || echo '?') barcodes)"
[[ -s "$STAR_BAM"     ]] && log "  STAR BAM             = ${STAR_BAM}"
[[ -d "$RAW_DIR"      ]] && log "  Gene/raw mtx         = ${RAW_DIR}"
[[ -s "$DUP_CSV"      ]] && log "  dup stats            = ${DUP_CSV}"
[[ -s "$UMI_MATRIX_CSV"  ]] && log "  UMI matrix CSV       = ${UMI_MATRIX_CSV}"
[[ -s "$READS_MATRIX_CSV" ]] && log "  Reads matrix CSV     = ${READS_MATRIX_CSV}"
[[ -s "${DMUX_DIR:-/nonexistent}"     ]] && :
[[ -d "${DMUX_BAM_DIR}/_done" && -s "${DMUX_BAM_DIR}/_done" ]] && log "  demux BAM              = ${DMUX_BAM_DIR}"
[[ -d "${DMUX_FQ_DIR}/_done"  && -s "${DMUX_FQ_DIR}/_done"  ]] && log "  demux FASTQ            = ${DMUX_FQ_DIR}"
log ""
log "Next step (outside this pipeline): downstream analysis (e.g. Seurat/Scanpy)"
log "clustering) on the produced matrices; this script stops at the"
log "pre-processing boundary."
exit 0
