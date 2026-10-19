#!/usr/bin/env bash
#=============================================================================
# brbseq_pipeline.sh — MERCURIUS DRUG-seq (BRB-seq) single-pool pipeline
#
# MERCURIUS 3' RNA-seq (BRB-seq/DRUG-seq) Data Pre-Processing pipeline
# for ONE pool/run, in step order:
#
#   STAR genome indexing            (one-time, NOT a pipeline step)
#   Creating barcode files          [MANDATORY ]  prepare_barcodes.sh
#   FASTQ quality (fastqc)          [optional   ]  --run-fastq
#   STARsolo alignment + counts     [MANDATORY ]  STAR
#   Read distribution (RSeQC)       [optional   ]  --run-rseqc
#   Duplication statistics          [MANDATORY ]  get_duplication_statistics.py
#   BAM demultiplex (picard)        [optional   ]  --run-dmux-bam
#   FASTQ demultiplex (fqtk)        [optional   ]  --run-dmux-fq
#   Merge multi-lane FASTQ          [note       ]  pre-merge or comma-list
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
#   * Portable across hosts: --conda-env activates a package environment
#     before running, --star selects the STAR binary, and helper scripts are
#     resolved from ./utils (shipped with this directory).
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
PCA_CSV="${UTILS_DIR}/prepare_for_pca.py"

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
FORCE_PLATE=0                    # pass --force-plate to the parser (deliberate kit/plate family override)
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
MANIFEST_HEADER=1                # 1 = has header row (default)
RUN_OPT_FASTQ="${RUN_OPT_FASTQ:-0}"
RUN_OPT_RSEQC="${RUN_OPT_RSEQC:-0}"
RUN_OPT_DEMUX_BAM="${RUN_OPT_DEMUX_BAM:-0}"
RUN_OPT_DEMUX_FQ="${RUN_OPT_DEMUX_FQ:-0}"
# Alithea BRB-seq production SAM tag set (fixed, not user-configurable).
SAM_ATTRS="NH HI nM AS CR UR CB UB GX GN sS sQ sM"
CONDA_ENV="${CONDA_ENV:-}"
STAR_BIN="${STAR:-}"             # default autodetected
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
  STAR genome indexing is a one-time setup step and is NOT executed here:
  built-in human/mouse genomes are resolved from --species; for any other
  species build an index once with build_star_index.sh (same directory) and
  pass it via --star-genome.

  MANDATORY:   barcode prep | STARsolo alignment + counts | duplication stats
               | PCA-ready count CSV (genes x samples)
  OPTIONAL:    fastqc (--run-fastq) | RSeQC (--run-rseqc) |
               BAM demux (--run-dmux-bam) | FASTQ demux (--run-dmux-fq)

  REFERENCES (built-in):  human, mouse  (--star-genome may override)
  REFERENCES (custom):    any other --species name + --star-genome

REQUIRED FLAGS
  --species     <sp>      human | mouse | <custom name>
  --r1-barcode  <fq.gz>   R1 barcode read (14 bp CB + 14 bp UMI)
  --r2-genomic  <fq.gz>   R2 genomic / cDNA read
  --run-name   <name>     run label for the output leaf BRBseq_<name>_<timestamp>;
                          in lib-prep mode taken from the workbook 'Researcher'
                          row (first name) unless given here.  Required in
                          --barcode-tsv mode (no workbook to read).
  pool label is INFERRED from the FASTQ names — everything before the lane
  marker "_S<n>" (kept verbatim).  e.g. pool_poolA_S1_R1_001.fastq.gz ->
  pool_poolA.  R1 and R2 must agree; a disagreement is a warning (R1 wins).

  MANIFEST SOURCE — exactly one of these two:
  (a)  --barcode-tsv <tsv>
        ready-made manifest: Sample_ID<TAB>Barcode (2 cols)
  (b)  --libprep-xlsx <xlsx>  (+ optional --barcodes-plate <file>)
        Mercurius lib-prep workbook (sample code + Oligo-dT primer well).
        The barcode plate is joined by well.  If --barcodes-plate is omitted,
        it is inferred from the workbook's 'Kit ref' row (V5A/V5D →
        barcodes_96_<family>_*.txt next to the xlsx); if given, the file's
        family is cross-checked against the Kit ref and a mismatch aborts the
        run (pass --force-plate to override a deliberate re-prep).
        Needs python+openpyxl — detected from the active python, the brbseq
        conda env, or a --conda-env you pin.


OPTIONAL FLAGS
  --star-genome <dir>     STAR genomeDir (built with build_star_index.sh);
                          required for custom --species, overrides built-in
  --bed-file    <bed12>   RSeQC BED12 gene model (override). If omitted,
                          bakes one at run time from the GTF shipped
                          inside the STAR index dir (build_star_index.sh
                          moves it there) — needs gffread in PATH
  --out-dir     <dir>     root output dir (default /mnt/data/BRB-seq_output)
  --threads     <n>       threads (default 8)
  --run-fastq               fastqc (needs fastqc in PATH)
  --run-rseqc               read_distribution (needs samtools,
                            read_distribution.py, and a BED12: --bed-file or
                            a GTF inside the STAR index dir)
  --run-dmux-bam            picard BAM demux (needs picard)
  --run-dmux-fq             fqtk FASTQ demux (needs fqtk; 14 bp CB + 14 bp UMI,
                            Alithea BRB-seq standard — hardcoded)
  --no-header               manifest has NO header row
                          (default: has a header)
  --force-plate             suppress the kit/plate family cross-check in
                          workbook mode (override a deliberate re-prep where
                          the --barcodes-plate family differs from the work-
                          book's 'Kit ref' row)
  --empty-wells             workbook mode: also emit every unoccupied plate
                          well as an empty_well_<WELL> column (old-pipeline
                          convention) so all plate wells show up in the
                          STARsolo matrix.  Reserved prefix: sample names
                          must not start with 'empty_well'.
  --conda-env  <name>       activate a conda env before running (portable hosts)
  --star       <bin>        STAR binary name (default autodetect: STAR-avx2,
                          STAR-avx, STAR)
  --resume     auto|force   auto = skip completed steps (default); force = rerun
  --dry-run                 emit planned commands and exit (no side effects)
  -h, --help                this help
  --version                 print version and exit

ENVIRONMENT OVERRIDES (advanced)
  STAR, UTILS_DIR, PYTHON, OUT_DIR, THREADS, AVT (legacy: treated as --run-name),
  SPECIES, RUN_NAME, RUN_OPT_FASTQ, RUN_OPT_RSEQC, RUN_OPT_DEMUX_BAM,
  RUN_OPT_DEMUX_FQ.

OUTPUT LAYOUT (under OUT_DIR):
  OUT_DIR/BRBseq_<run-name>_<timestamp>/
      barcodes/          barcode_whitelist.txt, fastq_sample_metadata
                          .tsv, bam_barcodes.tsv
      star/              <pool>/*.out.{bam,bam.bai,tar}, <pool>_Solo.out/
                         .../Gene/raw/{barcodes.tsv,features.tsv,
                         umiDedup-NoDedup.mtx, umiDedup-1MM_Directional.mtx}
                         [raw, barcode-labelled — never modified]
      star/samples/      the same matrix with barcode columns renamed to
                         Sample_IDs: barcodes.tsv, features.tsv,
                         umiDedup-<method>.mtx, sample_barcode_map.tsv
      rseqc/             <pool>_read_distribution.txt  (if --run-rseqc)
                          (+ <gtf>.bed12 auto-baked from the index-dir GTF
                          when --bed-file is not given)
      duplication_statistics.csv   barcodes x reads/UMIs/umi_%
      <pool>_counts_for_pca.csv    gene rows x sample cols (gene_symbol header);
                          the 1MM_Directional umiDedup MTX flattened for
                          downstream PCA / DE — replaces the legacy
                          _counts.detailed_for_pca.csv
      demux_bam/         per-sample BAMs                     (if --run-dmux-bam)
      demux_fastq/       per-sample FASTQs                   (if --run-dmux-fq)

EXAMPLES
  # Full mandatory pipeline on a mouse pool (built-in reference):
  ./brbseq_pipeline.sh \
      --species mouse \
      --r1-barcode  /data/AVT0608/pool_poolA_S1_R1_001.fastq.gz \
      --r2-genomic  /data/AVT0608/pool_poolA_S1_R2_001.fastq.gz \
      --barcode-tsv /data/AVT0608/sample_barcode.tsv \
      --run-name Eleni --threads 16 \
      --out-dir     /mnt/data/BRB-seq_output

  # Same, but feed the Mercurius lib-prep workbook directly (no hand-made
  # TSV; manifest goes to <run>/barcodes/, run name inferred from the
  # workbook 'Researcher' row).  The barcode plate is optional here: when
  # omitted it is inferred from the workbook's 'Kit ref' row — here V5A,
  # i.e. barcodes_96_V5A*.txt that sits next to the xlsx file:
  ./brbseq_pipeline.sh \
      --species mouse \
      --r1-barcode  /data/PHX537/pool_Eleni_S1_R1_001.fastq.gz \
      --r2-genomic  /data/PHX537/pool_Eleni_S1_R2_001.fastq.gz \
      --libprep-xlsx  "Mercurius HS BRB-seq RNA lib prep GECF Eleni_vA.02_PHX537.xlsx" \
      --out-dir /mnt/data/BRB-seq_output

  # Custom species: build a STAR index once with build_star_index.sh, then pass it:
  ./build_star_index.sh \
      --fasta Mus_musculus.GRCm39.dna.primary_assembly.fa \
      --gtf  Mus_musculus.GRCm39.114.gtf \
      --outdir star_index/Mus_musculus.GRCm39.114 \
      --sjdbOverhang 89
  ./brbseq_pipeline.sh \      --species mouse-custom --star-genome star_index/Mus_musculus.GRCm39.114 \
      ...

  # On a host that needs a conda env to see STAR / samtools / ...:
  conda activate brbseq      # or
  ./brbseq_pipeline.sh ... --conda-env brbseq

  # Dry run (print what WOULD be done, do not execute):
  ./brbseq_pipeline.sh ... --dry-run
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
    --force-plate)    FORCE_PLATE=1;                                           shift;;
    --empty-wells)    EMPTY_WELLS=1;                                           shift;;
    --species)      SPECIES="${2:?missing --species value}";          shift 2;;
    --star-genome)   STAR_INDEX="${2:?missing --star-genome value}";   shift 2;;
    --run-name)   RUN_NAME="${2:?missing --run-name value}";           shift 2;;
    --bed-file)     BED_FILE="${2:?missing --bed-file value}";      shift 2;;
    --out-dir)      OUT_DIR="${2:?missing --out-dir value}";        shift 2;;
    --threads)      THREADS="${2:?missing --threads value}";        shift 2;;
    --no-header)    MANIFEST_HEADER=0;                               shift;;
    --run-fastq)    RUN_OPT_FASTQ=1;                                  shift;;
    --run-rseqc)    RUN_OPT_RSEQC=1;                                  shift;;
    --run-dmux-bam) RUN_OPT_DEMUX_BAM=1;                              shift;;
    --run-dmux-fq)  RUN_OPT_DEMUX_FQ=1;                               shift;;
    --conda-env)    CONDA_ENV="${2:?missing --conda-env value}";      shift 2;;
    --star)         STAR_BIN="${2:?missing --star value}";           shift 2;;
    --resume)       RESUME="${2:?missing --resume value}";          shift 2;;
    --dry-run)      DRY_RUN=1;                                        shift;;
    -h|--help)      usage; exit 0;;
    --version)      echo "brbseq_pipeline.sh ${SCRIPT_VERSION}"; exit 0;;
    -*)             die "unknown option: $1 (see --help)";;
    *)              die "unexpected argument: $1 (see --help)";;
  esac
done

#--- validate required inputs ------------------------------------------------
# Exactly one manifest source must be provided:
#   --barcode-tsv  (2-col Sample_ID/Barcode TSV),  OR
#   --libprep-xlsx  (+ optional --barcodes-plate).
# The plate is optional in workbook mode: parse_libprep_xlsx.py infers it from
# the workbook's 'Kit ref' row (barcodes_96_<family>_*.txt next to the
# workbook) — or cross-checks an explicitly given plate against it (exit 4 on
# a family mismatch, --force-plate to override).
if [[ -n "$BARCODE_TSV" && -n "$LIBPREP_XLSX" ]]; then
  die "manifest conflict: pass either --barcode-tsv OR --libprep-xlsx, not both"
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
  die "--barcodes-plate is a workbook-mode flag: pass --libprep-xlsx as well — or drop the plate and let the workbook's 'Kit ref' row select it"
elif [[ "$LIBPREP_MODE" -ne 1 && -z "$BARCODE_TSV" ]]; then
  MISSING+=( "manifest source: --barcode-tsv <tsv>   or   --libprep-xlsx <xlsx>  (+ optional --barcodes-plate <file>)" )
fi
[[ -n "$SPECIES"     ]] || MISSING+=( "--species" )
[[ -n "$RUN_NAME" || -n "$AVT_ID" || "$LIBPREP_MODE" -eq 1 ]]   || MISSING+=( "--run-name <name>" )
[[ "$RESUME" == "auto" || "$RESUME" == "force" ]] || die "--resume must be auto or force"
[[ "${#MISSING[@]}" -eq 0 ]] || die "missing required flag(s): ${MISSING[*]} (see --help)"

#--- pool label ---------------------------------------------------------------
# The pool label is INFERRED from the FASTQ names (no user parameter).
# Rule: take the basename, cut at the lane marker "_S<n>" — everything
# before "_S" is kept verbatim.
#   pool_poolA_S1_R1_001.fastq.gz  -> pool_poolA
#   pool_Eleni_S2_R2_001.fastq.gz  -> pool_Eleni
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
  warn "pool label disagrees between FASTQ names: R1 suggests '$POOL_R1', R2 suggests '$POOL_R2' — using R1 ('$POOL_R1'). Rename the FASTQs to fix this."
  POOL="$POOL_R1"
elif [[ -n "$POOL_R1" ]]; then
  warn "could not infer pool from '$(basename "$R2_GENOMIC")' (no <pool>_<label>_S<lane> pattern); using '$POOL_R1' from R1. Rename the FASTQs to fix this."
  POOL="$POOL_R1"
elif [[ -n "$POOL_R2" ]]; then
  warn "could not infer pool from '$(basename "$R1_BARCODE")' (no <pool>_<label>_S<lane> pattern); using '$POOL_R2' from R2. Rename the FASTQs to fix this."
  POOL="$POOL_R2"
else
  die "could not infer the pool label from either FASTQ name (expected <pool>_<label>_S<lane>_R<read>_...fastq, e.g. pool_poolA_S1_R1_001.fastq.gz)"
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
  [[ -n "$inner" ]] || die "no STAR index found in $STAR_INDEX (no subfolder contains Genome + SA)"
  STAR_INDEX="$inner"
fi
if [[ -z "$STAR_INDEX" ]]; then
  die "--species '$SPECIES' has no built-in reference and no --star-genome was given.
  Build a STAR index once with build_star_index.sh and pass it:
    ./build_star_index.sh --fasta <genome.fa> --gtf <annotation.gtf> \
        --outdir <index-dir> --sjdbOverhang 89
    ./brbseq_pipeline.sh ... --species $SPECIES --star-genome <index-dir>"
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
#   0 GTFs         -> --run-rseqc without --bed-file dies at the gate below
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
      {
        echo "ERROR: $STAR_INDEX contains $n GTF annotations — cannot tell which one is the gene model to bake the RSeQC BED12 from:"
        printf '%s\n' "${RSEQC_GTF_CANDIDATES[@]}" | sed 's/^/    /' >&2
        echo "Fix: pass --bed-file <bed12> explicitly, or clean the index dir so it holds exactly one .gtf/.gtf.gz (one genome build per dir)." >&2
      }
      exit 1
    else
      local_list="$(printf '%s, ' "${RSEQC_GTF_CANDIDATES[@]}")"
      warn "RSeQC: $STAR_INDEX holds multiple GTF annotations (${local_list%, }) — no RSeQC reference can be selected from the index dir; use --bed-file explicitly"
    fi
  fi
fi
if [[ -z "$BED_FILE" && -n "$RSEQC_GTF" ]]; then
  log "RSeQC: no --bed-file given; will bake BED12 from $RSEQC_GTF (the single GTF in the STAR index dir — the gene model the index was built from)"
elif [[ -z "$BED_FILE" && -z "$RSEQC_GTF" ]]; then
  # Nothing to bake from; the gate below dies when --run-rseqc is enabled.
  :
fi
log "reference = ${SPECIES}  ->  STAR index $STAR_INDEX"
#--- (end RSeQC reference resolution) -------------------------------------------------

#--- validate file/directory inputs ------------------------------------------
[[ -f "$R1_BARCODE"  ]] || die "R1 barcode  read not found: $R1_BARCODE"
[[ -f "$R2_GENOMIC"  ]] || die "R2 genomic  read not found: $R2_GENOMIC"
# In lib-prep mode the manifest is generated later (needs the conda-env python);
# only validate the workbook here.  The plate is optional: when omitted, the
# parser infers it from the workbook's 'Kit ref' row.
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  [[ -f "$LIBPREP_XLSX" ]] || die "lib-prep workbook not found: $LIBPREP_XLSX"
  [[ -z "$BARCODES_PLATE" || -f "$BARCODES_PLATE" ]] || die "barcode plate not found: $BARCODES_PLATE"
else
  [[ -f "$BARCODE_TSV" ]] || die "manifest (barcode TSV) not found: $BARCODE_TSV"
fi
[[ -d "$STAR_INDEX"  ]] || die "STAR index dir  not found: $STAR_INDEX"
if [[ "$RUN_OPT_RSEQC" -eq 1 ]]; then
  if [[ -z "$BED_FILE" && -z "$RSEQC_GTF" ]]; then
    die "--run-rseqc needs a BED12 gene model: pass --bed-file <bed12>, or have a GTF inside the STAR index dir (build_star_index.sh moves it there) so we can bake one."
  fi
  if [[ -n "$BED_FILE" && ! -f "$BED_FILE" ]]; then
    die "BED12 file not found: $BED_FILE"
  fi
fi

#--- validate shipped helpers ------------------------------------------------
for req in "$PREPARE_BARCODES" "$DEMUX_BAM" "$DUP_STATS" "$RENAME_MATRIX" "$PCA_CSV"; do
  [[ -f "$req" ]] || die "missing helper script: $req  (must live in BRBseq/utils/)"
done
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  [[ -f "$PARSE_LIBPREP" ]] || die "missing helper script: $PARSE_LIBPREP  (BRBseq/utils/)"
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
    die "no run name — pass --run-name <name> (or use --libprep-xlsx which names the researcher)"
  fi
fi
RUN_NAME="$(tr -dc 'A-Za-z0-9._-' <<<"${RUN_NAME// /_}" )"   # filesystem-safe
[[ -n "$RUN_NAME" ]] || die "run name sanitised to empty"

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
PCA_OUT_CSV="${RUN_DIR}/${POOL}_counts_for_pca.csv"
# In --barcode-tsv mode the manifest lives where the user put it (no generation);
# in lib-prep mode the manifest is generated into <run>/barcodes/ and survives re-runs.
MANIFEST="$BARCODE_TSV"
if [[ "$LIBPREP_MODE" -eq 1 ]]; then
  MANIFEST="${BARCODES_DIR}/manifest_from_libprep.tsv"
  LIBPREP_PROV="${MANIFEST}_provenance.tsv"
else
  LIBPREP_PROV=""
fi

#--- tool discovery ----------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }
if [[ -z "$STAR_BIN" ]]; then
  for cand in STAR-avx2 STAR-avx STAR; do
    if have "$cand"; then STAR_BIN="$cand"; break; fi
  done
fi
HAVE_STAR=0;      [[ -n "$STAR_BIN" ]] && have "$STAR_BIN" && HAVE_STAR=1
HAVE_FASTQC=0;    have fastqc              && HAVE_FASTQC=1
HAVE_SAMTOOLS=0;  have samtools            && HAVE_SAMTOOLS=1
HAVE_RSEQC=0;     have read_distribution.py && HAVE_RSEQC=1
HAVE_PICARD=0;    have picard              && HAVE_PICARD=1
HAVE_FQTK=0;      have fqtk                && HAVE_FQTK=1
have "$PYTHON"    || PYTHON="$(have python3 && echo python3 || echo python)"

#--- conda env (portable hosts) ---------------------------------------------
# In lib-prep mode the manifest parser needs openpyxl.  If the active python
# can already import it, nothing extra is required.  Otherwise, prefer the
# brbseq env — via --conda-env activation when conda is available, or by
# pointing directly at the env's python when the env dir is on this host
# (no conda needed).  Only if neither works do we fail, with a message that
# names the missing piece.
py_has_openpyxl() { [[ "$1" != "" ]] && "$1" -c "import openpyxl" >/dev/null 2>&1; }
py_ok=0
py_has_openpyxl "$PYTHON" && py_ok=1

if [[ "$LIBPREP_MODE" -eq 1 && "$py_ok" -ne 1 ]]; then
  resolved_openpyxl=0
  if command -v conda >/dev/null 2>&1; then
    # try the brbseq env (or the pinned one) via conda
    want="${CONDA_ENV:-brbseq}"
    if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$want"; then
      log "activating conda env for openpyxl: ${want}"
      # shellcheck disable=SC1091
      eval "$(conda shell.bash hook)" 2>/dev/null || true
      if conda activate "$want" 2>/dev/null; then
        have "$PYTHON" || PYTHON="python3"
        py_has_openpyxl "$PYTHON" && resolved_openpyxl=1
      fi
    fi
  fi
  if [[ "$resolved_openpyxl" -ne 1 ]]; then
    # conda not usable; fall back to the env's python directly (if on host)
    for cand in \
      "${HOME:-/nonexistent}/.conda/envs/brbseq/bin/python3" \
      "${HOME:-/nonexistent}/miniconda3/envs/brbseq/bin/python3" \
      "/opt/conda/envs/brbseq/bin/python3"; do
      if [[ -x "$cand" ]] && py_has_openpyxl "$cand"; then
        PYTHON="$cand"
        resolved_openpyxl=1
        break
      fi
    done
  fi
  if [[ "$resolved_openpyxl" -ne 1 ]]; then
    die "lib-prep mode needs python+openpyxl; none available.  Fix: pass --conda-env <name with openpyxl> AND activate conda, or run with a python that can 'import openpyxl' (e.g. the brbseq conda env's python3)."
  fi
fi

if [[ -n "$CONDA_ENV" ]]; then
  have conda || die "--conda-env requested but 'conda' not in PATH"
  log "activating conda env: ${CONDA_ENV}"
  # shellcheck disable=SC1091
  eval "$(conda shell.bash hook)"
  conda activate "$CONDA_ENV" || die "conda activate ${CONDA_ENV} failed"
  # re-detect tools after activation
  HAVE_STAR=0;      [[ -n "$STAR_BIN" ]] && have "$STAR_BIN" && HAVE_STAR=1
  HAVE_FASTQC=0;    have fastqc              && HAVE_FASTQC=1
  HAVE_SAMTOOLS=0;  have samtools            && HAVE_SAMTOOLS=1
  HAVE_RSEQC=0;     have read_distribution.py && HAVE_RSEQC=1
  HAVE_PICARD=0;    have picard              && HAVE_PICARD=1
  HAVE_FQTK=0;      have fqtk                && HAVE_FQTK=1
fi

#--- generate manifest from lib-prep workbook (workbook mode) ----------------
# Must run after conda activation: parse_libprep_xlsx.py needs openpyxl, which
# ships in the brbseq env (or any --conda-env the caller supplies).
if [[ "$LIBPREP_MODE" -eq 1 && ! -s "$MANIFEST" ]]; then
  mkdir -p "$BARCODES_DIR"
  log "lib-prep workbook mode: generating manifest from $LIBPREP_XLSX"
  parse_args=( --xlsx "$LIBPREP_XLSX" --out "$MANIFEST" --provenance "$LIBPREP_PROV" )
  if [[ -n "$BARCODES_PLATE" ]]; then
    parse_args+=( --plate "$BARCODES_PLATE" )
  fi
  if [[ "$FORCE_PLATE" -eq 1 ]]; then
    parse_args+=( --force-plate )
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
    if [[ "$parse_rc" -eq 4 ]]; then
      die "kit/plate barcode family mismatch (see error above) — re-run with a --barcodes-plate whose family matches the workbook's 'Kit ref' row, or override deliberately with --force-plate"
    fi
    die "parse_libprep_xlsx.py failed (exit ${parse_rc}; see error above) — in workbook mode the plate is inferred from the workbook's 'Kit ref' row (needs a barcodes_96_<family>_*.txt next to the xlsx), or pass --barcodes-plate <file>"
  fi
  [[ -s "$MANIFEST" ]] || die "generated manifest is empty: $MANIFEST"
  log "manifest   = ${MANIFEST}  (generated from lib-prep workbook)"
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
    die "manifest ${MANIFEST} has no empty_well_* rows but this call passes --empty-wells (stale artifact from a run without the flag).  Remove the stale manifest (or delete the run dir) so it regenerates, then re-run — the barcodes step will rebuild the whitelist accordingly"
  fi
  if [[ "$EMPTY_WELLS" -eq 0 && "$ew_in_manifest" -eq 1 ]]; then
    die "the existing manifest ${MANIFEST} contains empty_well_* rows, but this call does NOT pass --empty-wells.  Remove the stale manifest (or delete the run dir) and re-run, or start a fresh run"
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
log "manifest   = ${MANIFEST}  (header=${MANIFEST_HEADER})"
log "STAR index = ${STAR_INDEX}   (star=${STAR_BIN:-NOT_FOUND})"
log "output     = ${RUN_DIR}"
log "optional:  fastqc=${RUN_OPT_FASTQ}  rseqc=${RUN_OPT_RSEQC}  demux_bam=${RUN_OPT_DEMUX_BAM}  demux_fq=${RUN_OPT_DEMUX_FQ}"
log "tools:     STAR=${HAVE_STAR}  fastqc=${HAVE_FASTQC}  samtools=${HAVE_SAMTOOLS}  rseqc=${HAVE_RSEQC}  picard=${HAVE_PICARD}  fqtk=${HAVE_FQTK}"
log "python     = ${PYTHON}"
if [[ -n "$LIBPREP_PROV" ]]; then
  log "lib-prep   = ${LIBPREP_XLSX}"
  if [[ -n "$BARCODES_PLATE" ]]; then
    plate_note=""
    [[ "$FORCE_PLATE" -eq 1 ]] && plate_note="  [--force-plate]"
    log "plate      = ${BARCODES_PLATE}${plate_note}"
  else
    log "plate      = (not given — inferred from the workbook's 'Kit ref' row)"
  fi
  [[ "$EMPTY_WELLS" -eq 1 ]] && log "empty wells = ON (unoccupied plate wells emitted as empty_well_<WELL>)"
fi

#--- execute (or dry-run) ----------------------------------------------------
if [[ "$DRY_RUN" -eq 1 ]]; then
  log "---- DRY RUN: planned steps ----"
  log "  prepare_barcodes.sh -i <manifest> -w <OUT>/barcodes/barcode_whitelist.txt -f <OUT>/barcodes/fastq_sample_metadata.tsv -b <OUT>/barcodes/bam_barcodes.tsv -$([[ "$MANIFEST_HEADER" -eq 1 ]] && echo H)"
  if [[ "$RUN_OPT_FASTQ" -eq 1 ]]; then
    log "  fastqc  --threads 2 --outdir <OUT>/fastqc <R1> <R2>"
  fi
  log "  STAR    (one STARsolo call, both dedup methods in list form)"
  log "      ${STAR_BIN} --runMode alignReads --runThreadN ${THREADS} --genomeDir ${STAR_INDEX} --readFilesIn <R2> <R1> --soloCBwhitelist <OUT>/barcodes/barcode_whitelist.txt --soloCBlen 14 --soloUMIlen 14 --soloUMIstart 15 --soloUMIdedup NoDedup 1MM_Directional --outSAMattributes '${SAM_ATTRS}' --quantMode GeneCounts ..."
  log "      consolidate  <prefix>Solo.out/Gene/raw/umiDedup-<method>.mtx -> <OUT>/star/Gene/raw/"
  if [[ "$RUN_OPT_RSEQC" -eq 1 ]]; then
    if [[ -n "$BED_FILE" ]]; then
      log "  RSeQC   reference (given --bed-file): ${BED_FILE}"
    else
      log "  RSeQC   bake:   gffread ${RSEQC_GTF} --bed > <OUT>/rseqc/$(basename "${RSEQC_GTF%.gz}").bed12"
    fi
    log "      samtools index <OUT>/star/${POOL}_Aligned.sortedByCoord.out.bam"
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
      die "${name} requires tool '${tool}' which is not in PATH."
    else
      warn "${name}: optional step skipped (tool '${tool}' not in PATH)"
      return 0
    fi
  fi
  log "RUN  ${name}: ${cmd[*]}"
  local rc=0
  "${cmd[@]}" || rc=$?
  if [[ "$rc" -ne 0 ]]; then
    if [[ "$mode" == "mandatory" ]]; then
      die "${name} failed (exit ${rc}) — mandatory step"
    else
      warn "${name} failed (exit ${rc}) — optional step, continuing"
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
[[ -s "$WHITELIST"    ]] || die "barcode prep failed: whitelist is empty"
[[ -f "$FQ_META"      ]] || die "barcode prep failed: fastq_sample_metadata.tsv missing"
[[ -f "$BAM_BARCODES" ]] || die "barcode prep failed: bam_barcodes.tsv missing"
log "barcode prep done: whitelist=$(wc -l < "$WHITELIST") barcodes; fastq_meta+bam_map ready"

#=== FASTQ Quality Assessment  [OPTIONAL] ==================================
if [[ "$RUN_OPT_FASTQ" -eq 1 ]]; then
  log ""; log "fastqc (optional)"
  FQ_MARKER="${FASTQC_DIR}/QC_done"
  run_or_skip "fastqc" "$FQ_MARKER" "fastqc" "optional" \
    fastqc --threads 2 --outdir "$FASTQC_DIR" "$R1_BARCODE" "$R2_GENOMIC"
  touch "$FQ_MARKER"  # idempotent completion marker
  log "fastqc outputs in: $FASTQC_DIR"
else
  log ""; log "fastqc — not requested (use --run-fastq to enable)"
fi

#=== STARsolo  [MANDATORY] =================================================
log ""; log "STARsolo alignment + UMI count matrix (mandatory)"
[[ "$HAVE_STAR" -eq 1 ]] || die "STAR binary not found (looked for STAR-avx2/STAR-avx/STAR); use --star or --conda-env."
[[ -d "$STAR_INDEX" ]] || die "STAR index dir not found: $STAR_INDEX"
[[ -f "${STAR_INDEX}/Genome" ]] || warn "STAR index: no 'Genome' file in $STAR_INDEX — STAR will fail if index is incomplete."
[[ -f "${STAR_INDEX}/SA"     ]] || warn "STAR index: no 'SA' file in $STAR_INDEX — STAR will fail if index is incomplete."

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
  local -a attr_tokens=()
  read -r -a attr_tokens <<< "${SAM_ATTRS}"
  local STAR_ATTR_ARGS=( --outSAMattributes "${attr_tokens[@]}" )
  log "RUN  STARsolo (--soloUMIdedup ${DEDUP_METHODS[*]}, list form) -> ${prefix}"
  if "${STAR_BIN}" \
      --runMode alignReads \
      --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" \
      --readFilesIn "${R2_GENOMIC}" "${R1_BARCODE}" \
      --readFilesCommand zcat \
      --soloCBwhitelist "${WHITELIST}" \
      --soloCBstart 1 \
      --soloCBlen 14 \
      --soloUMIstart 15 \
      --soloUMIlen 14 \
      --soloUMIdedup "${DEDUP_METHODS[@]}" \
      --soloType CB_UMI_Simple \
      --soloCBmatchWLtype 1MM \
      --soloStrand Forward \
      --soloCellFilter None \
      --soloBarcodeReadLength 0 \
      --outFilterMultimapNmax 1 \
      "${STAR_ATTR_ARGS[@]}" \
      --outSAMmapqUnique 60 \
      --outSAMunmapped Within \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN 1 \
      --outSAMmultNmax 1 \
      --outFilterIntronMotifs None \
      --outFileNamePrefix "${prefix}" \
      --quantMode GeneCounts  > "${log_file}" 2>&1
  then
    return 0
  else
    log "ERROR: STARsolo (--soloUMIdedup ${DEDUP_METHODS[*]}) failed; tail of ${log_file}:"
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
  [[ -n "$solo_dir" && -d "$solo_dir" ]] || die "consolidation: no Solo.out dir under $STAR_DIR/"
  log "consolidation: solo_dir=$solo_dir"
  local method kind
  for method in "${DEDUP_METHODS[@]}"; do
    local src
    for src in "${solo_dir}/Gene/raw/umiDedup-${method}.mtx" \
               "${solo_dir}/Gene/umiDedup-${method}.mtx"; do
      [[ -f "$src" ]] && { ln -f "$src" "${dest}/umiDedup-${method}.mtx" 2>/dev/null || cp -f "$src" "${dest}/umiDedup-${method}.mtx"; break; }
    done
    [[ -f "${dest}/umiDedup-${method}.mtx" ]] || die "consolidation: no umiDedup-${method}.mtx anywhere under $solo_dir"
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
  [[ -n "$STAR_BAM" && -s "$STAR_BAM" ]] || die "finished but no Aligned.sortedByCoord.out.bam found under $STAR_DIR"
  echo "$(date +%FT%T)" > "${STAR_DIR}/_s14_done"
fi
log "STAR BAM: $STAR_BAM"
log "consolidated mtx dir: ${STAR_DIR}/Gene/raw"
SOLO_OUT="${STAR_DIR}/Gene/raw"   # alias used by the dup-stats step

#--- bake helper: GTF -> BED12 (house recipe: gffread --bed) -----------------
# Same conversion as flashseq_rseqc.sh so both pipelines produce identical
# BED12 from the same GTF.  Writes $RSEQC_DIR/<gtf-basename>.bed12
# (.gz GTFs are decompressed on the fly).  Sets global RSEQC_BED to the
# output path (NOT via stdout capture — log() writes to stdout, so the value
# is returned in a variable to avoid mixing prose into the path).
RSEQC_BED=""
bake_rseqc_bed12() {
  local gtf="$1" out
  out="${RSEQC_DIR}/$(basename "${gtf%.gz}").bed12"
  if [[ -s "$out" ]]; then
    log "BED12 already baked (resume): $out  (from $gtf)"
    RSEQC_BED="$out"
    return 0
  fi
  command -v gffread >/dev/null 2>&1 || die "gffread not in PATH — needed to bake the BED12 from $gtf (add the env with cufflinks, pass --bed-file, or activate with --conda-env)"
  local tmp rc=1
  tmp="$(mktemp "${TMPDIR:-/tmp}/brbseq_bed12.XXXXXX")"
  # pipefail/errexit are off inside if-conditions: capture rc explicitly.
  if   [[ "$gtf" == *.gz  ]]; then gzip -dc "$gtf" | gffread --bed -           > "$tmp"; rc=$?
  elif [[ "$gtf" == *.bz2 ]]; then bzcat -d "$gtf" | gffread --bed -           > "$tmp"; rc=$?
  else                             gffread       "$gtf" --bed                 > "$tmp"; rc=$?; fi
  if [[ $rc -ne 0 ]] || [[ ! -s "$tmp" ]]; then
    rm -f "$tmp"
    die "bake failed: gffread --bed produced no usable output for $gtf (rc=$rc) — check the GTF (gene/exon entries, quotes, contig names)"
  fi
  local cols
  cols="$(awk -F'\t' '!/^#/{print NF; exit}' "$tmp")"
  if [[ "${cols:-0}" -lt 12 ]]; then
    rm -f "$tmp"
    die "bake failed: output for $gtf has ${cols:-0} column(s) (< 12) — not a valid BED12 gene-model file"
  fi
  mv -f "$tmp" "$out"
  local rows; rows="$(wc -l < "$out" | tr -d '[:space:]')"
  log "baked BED12 (RSeQC): $out   [from GTF: $gtf]  (${rows} rows x ${cols} cols)"
  log "provenance:  gffread $gtf --bed > $out"
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
    [[ -n "$RSEQC_GTF" ]] || die "cannot bake BED12: no GTF found in $STAR_INDEX (build_star_index.sh moves it there) and no --bed-file was given."
    # Bake only when RSeQC actually runs; the bake is idempotent
    # (bake_rseqc_bed12 reuses an existing .bed12).  Sets RSEQC_BED.
    bake_rseqc_bed12 "$RSEQC_GTF"
    [[ -n "$RSEQC_BED" && -s "$RSEQC_BED" ]] || die "bake failed: no BED12 produced."
    BED_FILE="$RSEQC_BED"
  fi
  log "RSeQC reference in use: $BED_FILE"
  log "RSeQC provenance: baked from GTF $RSEQC_GTF   (recipe: gffread <gtf> --bed > <gtf>.bed12)"
  run_or_skip "samtools index" "${STAR_BAM}.bai" "samtools" "optional" \
    samtools index -@ 8 "$STAR_BAM"
  run_or_skip "read_distribution" "$RSEQC_OUT" "read_distribution.py" "optional" \
    read_distribution.py -i "$STAR_BAM" -r "$BED_FILE" > "$RSEQC_OUT"
  if [[ -s "$RSEQC_OUT" ]]; then
    log "distribution report: $RSEQC_OUT"
  fi
else
  log ""; log "RSeQC — not requested (use --run-rseqc to enable)"
fi

#--- dup-stats helpers (kept here to anchor the section) --------------------
#=== Duplication statistics  [MANDATORY] ==================================
log ""; log "duplication statistics (mandatory)"
RAW_DIR="${STAR_DIR}/Gene/raw"
for req in "$RAW_DIR/barcodes.tsv" \
           "$RAW_DIR/umiDedup-NoDedup.mtx" \
           "$RAW_DIR/umiDedup-1MM_Directional.mtx"; do
  [[ -s "$req" ]] || die "required dup-stats matrix input missing: $req"
done
run_or_skip "duplication stats" "$DUP_CSV" "" "mandatory" \
  "$PYTHON" "$DUP_STATS" \
    --raw_dir "$RAW_DIR" \
    --dedup_method "1MM_Directional" \
    --output_csv "$DUP_CSV"
[[ -s "$DUP_CSV" ]] || die "dup stats failed: duplication_statistics.csv not produced"
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
[[ -s "${RENAME_DIR}/barcodes.tsv" ]] || die "rename failed: no remapped barcodes.tsv"
[[ -s "${RENAME_DIR}/umiDedup-1MM_Directional.mtx" ]] || die "rename failed: no remapped umiDedup-1MM_Directional.mtx"
log "sample-named matrix dir: ${RENAME_DIR}"
log "  (barcode->sample map: ${RENAME_DIR}/sample_barcode_map.tsv; raw matrix unchanged in ${RAW_DIR})"

#=== PCA-ready count CSV  [MANDATORY] ========================================
# Flatten the 1MM_Directional sample-named umiDedup matrix into a single
# gene x sample CSV.  One row per gene (duplicate symbols folded by summing,
# matching HTSeq 'union' semantics), one column per Sample_ID from barcodes.tsv.
# This is the "ready for PCA" table — the equivalent of the legacy
# <pool>_counts.detailed_forPCA.csv — and the input for any downstream
# clustering/differential-expression step.
log ""
log "PCA-ready count CSV (mandatory)"
run_or_skip "PCA count CSV" "$PCA_OUT_CSV" "" "mandatory" \
  "$PYTHON" "$PCA_CSV" \
    --samples-dir "$RENAME_DIR" \
    --out-csv   "$PCA_OUT_CSV" \
    --dedup     "1MM_Directional"
[[ -s "$PCA_OUT_CSV" ]] || die "PCA CSV step failed: $PCA_OUT_CSV not produced"
log "PCA-ready CSV: $PCA_OUT_CSV"

#=== Sample demultiplex (BAM, picard)  [OPTIONAL] ==========================
if [[ "$RUN_OPT_DEMUX_BAM" -eq 1 ]]; then
  log ""; log "BAM demultiplex via picard (optional)"
  [[ -f "$BAM_BARCODES" ]] || die "needs bam_barcodes.tsv (produced in barcode prep)"
  [[ -s "$STAR_BAM"     ]] || die "needs the STAR BAM from the STAR step"
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
      warn "optional step skipped: 'picard' not in PATH (see --conda-env or install)" 
    fi
  fi
else
  log ""; log "BAM demux — not requested (use --run-dmux-bam + picard)"
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
        || warn "fqtk reported non-zero exit (check log)" 
      touch "$DONE_FQ"
      log "per-sample FASTQs in: $DMUX_FQ_DIR"
    else
      warn "optional step skipped: 'fqtk' not in PATH (see --conda-env or install)" 
    fi
  fi
else
  log ""; log "FASTQ demux — not requested (use --run-dmux-fq + fqtk)"
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
[[ -s "$PCA_OUT_CSV"  ]] && log "  PCA-ready CSV        = ${PCA_OUT_CSV}"
[[ -s "${DMUX_DIR:-/nonexistent}"     ]] && :
[[ -d "${DMUX_BAM_DIR}/_done" && -s "${DMUX_BAM_DIR}/_done" ]] && log "  demux BAM              = ${DMUX_BAM_DIR}"
[[ -d "${DMUX_FQ_DIR}/_done"  && -s "${DMUX_FQ_DIR}/_done"  ]] && log "  demux FASTQ            = ${DMUX_FQ_DIR}"
log ""
log "Next step (outside this pipeline): downstream analysis (e.g. Seurat/Scanpy)"
log "clustering) on the produced matrices; this script stops at the"
log "pre-processing boundary."
exit 0
