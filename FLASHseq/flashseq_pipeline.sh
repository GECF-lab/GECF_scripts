#!/usr/bin/env bash
#=============================================================================
# flashseq_pipeline.sh  — single-cell FLASH-seq read counting
#
# FASTQ -> FastQC -> STAR alignment -> featureCounts
# -> merged genes x samples count matrix (gene_id). Per-sample uniquely
# mapped % is reported to _gates.csv (warning if <70%, never a stopper).
# CLI form of Part 4 of the Alithea "Mercurius Single-Cell FLASH-seq"
# guide (May 2026).
#
# Per sample, in order:
#   1. FastQC        raw-read QC
#   2. STAR          splice-aware alignment to the reference genome
#   3. QC report     warn if <70% of reads uniquely mapped (no filtering)
#   4. featureCounts gene-level read counts
# All samples:
#   6. R merge       combine into one genes x samples matrix; per-cell genes
#                     detected (counts >=1 and >=5) -> matrix/genes_per_cell.csv
#   6b. metadata    best-effort sample->group CSV from the matrix header
#                     -> matrix/counts_matrix_metadata.csv
#                     (failure is a warning only, never stops the pipeline)
#
# Re-runnable: a step is skipped when its output already exists
# (--no-resume forces a full redo). One sample failing does not stop
# the others — it is flagged and skipped. STAR/featureCounts parameters
# follow the Alithea guide.
#=============================================================================
set -uo pipefail


# --- defaults ---------------------------------------------------------------
SCRIPT_VERSION="1.2.0"   # Bump on any change that alters outputs/behaviour.
AVT=""
SPECIES=""
STAR_BIN="STAR-avx2"
# fastq root (root of the avid* sample dirs) — user MUST supply via --fastq-dir
FASTQ_DIR=""
STAR_INDEX=""
OUT=""
MAX_JOBS=4
THREADS_PER=8
DO_FASTQC=1
SKIP_STAR=0
COUNT_ONLY=0
SAMPLES=""
EXCLUDE=""
RESUME=1
DRY_RUN=0
REF_FASTA=""       # reference FASTA to build a STAR index (optional for custom species)
GTF_ARG=""         # annotation override (--gtf)
REF_NAME=""        # gene-set label for the matrix (custom species)
SJDB_OVERHANG=74   # STAR splice-database overhang used when building an index

# Fixed output root. NOT checked for existence/writability — the mount lives
# on the storage host, which is not necessarily this machine. --outdir overrides.
OUT_BASE="/mnt/data/FLASHseq_output"

# --- Built-in references ------------------------------------------------------
# One line per built-in species, pipe-delimited:
#     <species>|<STAR index dir>|<GTF file>|<gene-set label>
#   <species>        value of --species that selects this line (human/mouse)
#   <STAR index dir> prebuilt STAR genome dir (must contain Genome + SA)
#   <GTF file>       annotation for featureCounts (-a)
#   <gene-set label> cosmetic name shown in logs/summary (e.g. GRCh38.116)
#
# To use custom built-in references: edit the two lines below. To use a
# species NOT in this table, pass any name to --species plus --gtf and a
# prebuilt --star-index (see usage).
REFS=(
  "human|/home/gecf/references/STAR_references/FLASHseq/current_human_star/Homo_sapiens.GRCh38.116|/home/gecf/references/STAR_references/FLASHseq/current_human_star/Homo_sapiens.GRCh38.116.gtf|GRCh38.116"
  "mouse|/home/gecf/references/STAR_references/FLASHseq/current_mouse_star/Mus_musculus.GRCm39.116|/home/gecf/references/STAR_references/FLASHseq/current_mouse_star/Mus_musculus.GRCm39.116.gtf|GRCm39.116"
)

GTF=""; GENESET=""
declare -A SAMPLES_LIST
ID=""; R1=""; R2=""; STAR_BAM=""; GATE_PCT="NA"
LBL=""                                  # sample label for matrix columns (from FASTQ name)
declare -A LBL_SEEN=()                  # collision guard: label -> owning sample dir
RUN_LOG=""; PIDS=()

# ---------------------------------------------------------------------------
# Log helpers. All output is timestamped and tagged:
#   [INFO]/[WARN] = progress / advisory, [ERROR] + die() = aborted.
# ---------------------------------------------------------------------------
info() { echo "$(date +%H:%M:%S) [INFO]   $1"; }
warn() { echo "$(date +%H:%M:%S) [WARN]   $1"; }
err()  { echo "$(date +%H:%M:%S) [ERROR]  $1"; }
die()  { err "$1"; exit 1; }

usage() {
cat <<'USAGE'
Single-cell FLASH-seq: FASTQ -> FastQC -> STAR -> featureCounts ->
one merged genes x samples count matrix
(gene_id), plus per-cell genes detected at two thresholds (>=1 and >=5)
-> matrix/genes_per_cell.csv (informational).  Per-sample uniquely
mapped % is reported to _gates.csv (warning, not a filter).  CLI form of
Part 4 of the Alithea "Mercurius Single-Cell FLASH-seq" guide (May 2026).
Re-runnable: finished steps are skipped (--no-resume forces a full redo).
One sample failing does not stop the others.

REQUIRED (all species)
  --species  human | mouse | <name>
             human/mouse = built-in Ensembl STAR references under
             /home/gecf/references/STAR_references/FLASHseq/;
             anything else = custom: --gtf required, plus either --star-index
             (prebuilt) or --ref-fasta (index is built from FASTA + your GTF).
  --fastq-dir  DIR   root containing the avid*/ sample dirs.

OPTIONAL
  --avt AVT0551          run label; inferred from --fastq-dir if omitted
  --samples a,b          process only these sample dirs (default: all)
  --exclude a,b          skip these avidNNNNN sample dirs (default: none)
  --max-jobs N           samples in parallel (default 4)
  --threads-per N        cores per STAR/featureCounts job (default 8)
  --outdir PATH          output root (default /mnt/data/FLASHseq_output)
  Custom species:
  --gtf PATH             GTF annotation (REQUIRED for custom species)
  --star-index DIR       prebuilt STAR index
  --ref-fasta PATH       genome FASTA; STAR index is built from this + your --gtf
  --ref-name LABEL       gene-set label in the matrix (default custom)
  --sjdb-overhang N      STAR --sjdbOverhang for index building (default 74)
  Stages / resuming:
  --skip-fastqc          omit FastQC
  --skip-star            reuse existing STAR BAM (skip alignment)
  --count-only           STAR+BAM already exist; just count + merge
  --no-resume            never skip finished steps
  --dry-run              print planned commands, execute nothing
  --version              print the script version and exit

EXAMPLES
  # typical: human, AVT and samples inferred/auto
  flashseq_pipeline.sh --species human --fastq-dir /data/AVTs/AVT0551/Samples
  # single sample, 32 cores, one job at a time
  flashseq_pipeline.sh --species human --fastq-dir /data/AVTs/AVT0551/Samples \
      --samples avid13442 --max-jobs 1 --threads-per 32
  # custom species with a prebuilt STAR index
  flashseq_pipeline.sh --species zebrafish --fastq-dir /data/AVTs/AVT0600/zf/Samples \
      --star-index /data/references/zebrafish_star --gtf /data/references/zebrafish.gtf \
      --ref-name GRCz11
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --avt)          AVT="$2"; shift 2;;
    --species)      SPECIES="$2"; shift 2;;
    --fastq-dir)    FASTQ_DIR="$2"; shift 2;;
    --star-index)   STAR_INDEX="$2"; shift 2;;
    --ref-fasta)    REF_FASTA="$2"; shift 2;;
    --gtf)          GTF_ARG="$2"; shift 2;;
    --ref-name)     REF_NAME="$2"; shift 2;;
    --sjdb-overhang) SJDB_OVERHANG="$2"; shift 2;;
    --outdir)       OUT="$2"; shift 2;;
    --max-jobs)     MAX_JOBS="$2"; shift 2;;
    --threads-per)  THREADS_PER="$2"; shift 2;;
    --skip-fastqc)  DO_FASTQC=0; shift;;
    --skip-star)    SKIP_STAR=1; shift;;
    --count-only)   COUNT_ONLY=1; DO_FASTQC=0; SKIP_STAR=1; shift;;
    --samples)      SAMPLES="$2"; shift 2;;
    --exclude)      EXCLUDE="$2"; shift 2;;
    --no-resume)    RESUME=0; shift;;
    --dry-run)      DRY_RUN=1; shift;;
    --version)      echo "flashseq_pipeline.sh $SCRIPT_VERSION"; exit 0;;
    -h|--help)      usage; exit 0;;
    *) die "Unknown option '$1' (see --help)";;
  esac
done

# --- tool resolution --------------------------------------------------------
# Prefer tools on PATH from an active conda env; otherwise fall back to the
# known flashseq env locations. Hard fail if a core tool is still missing.
_prepend_if_present() {
  local d="$1"
  [[ -d "$d" ]] || return 0
  case ":$PATH:" in *":$d:"*) : ;; *) export PATH="$d:$PATH" ;; esac
}
if ! command -v "$STAR_BIN"     >/dev/null 2>&1 \
   || ! command -v featureCounts >/dev/null 2>&1 \
   || ! command -v Rscript        >/dev/null 2>&1; then
  [[ -n "${CONDA_PREFIX:-}" ]] && _prepend_if_present "${CONDA_PREFIX}/bin"
  _prepend_if_present "/home/gecf/softwares/anaconda/envs/flashseq/bin"
  _prepend_if_present "/home/hermes/.conda/envs/flashseq/bin"
fi

_missing=""
command -v "$STAR_BIN"      >/dev/null 2>&1 || _missing="$_missing $STAR_BIN"
command -v featureCounts    >/dev/null 2>&1 || _missing="$_missing featureCounts"
command -v Rscript          >/dev/null 2>&1 || _missing="$_missing Rscript"
if [[ -n "$_missing" ]]; then
  die "core tool(s) not found on PATH:$_missing
Fix: activate the flashseq conda env and re-run."
fi

# The merge step requires these R packages; verify at launch.
_rmissing="$(Rscript -e 'm <- character(0)
    for (p in c("dplyr","readr","stringr","tibble","tidyr"))
        if (!requireNamespace(p, quietly = TRUE)) m <- c(m, p)
    if (length(m)) cat(paste(m, collapse = " "))
    ' 2>/dev/null)" || _rmissing="Rscript error"
[[ -n "$_rmissing" ]] && die "R packages missing for the merge step: $_rmissing"

for _t in "$STAR_BIN" featureCounts Rscript; do
  info "${_t} -> $(command -v "$_t" 2>/dev/null || echo MISSING)"
done

[[ -n "$SPECIES" ]] || { err "--species is required (human|mouse|<custom>)"; echo; usage; exit 1; }
[[ -n "$FASTQ_DIR" ]] || { err "--fastq-dir is required (root of the avid* sample dirs)"; echo; usage; exit 1; }

# --- Conflict guard: built-in species + --ref-fasta + no --star-index ---------
# --ref-fasta is only consumed when building an index (STAR --runMode
# genomeGenerate). For a built-in species (human|mouse) without an explicit
# --star-index, the index defaults to the built-in dir, which already has
# Genome+SA, so the build path is never reached and --ref-fasta is silently
# ignored. That is the same silent-ignoring class as the old custom-species
# misrouting, so hard-stop it here rather than proceeding with the built-in
# genome.
#   Allowed:   --species human --gtf X            (annotation override)
#   Allowed:   --species human --star-index I     (user prebuilt index)
#   Allowed:   --species human --star-index J --ref-fasta R --gtf X
#               (J empty -> build into J from R)
#   Rejected:  --species human --ref-fasta R      (no target dir to build into)
case "${SPECIES:-}" in
  human|mouse)
    if [[ -n "$REF_FASTA" && -z "$STAR_INDEX" ]]; then
      err "--ref-fasta (build an index) was given with the built-in species '$SPECIES' but no --star-index."
      err "That combination has no place to build the index, so --ref-fasta would be silently"
      err "ignored and the built-in genome would be used instead."
      err "Options:"
      err "  - drop --ref-fasta (use the built-in '$SPECIES' index + GTF), or"
      err "  - add --star-index <fresh-dir> to build a STAR index in that dir from --ref-fasta + --gtf."
      exit 1
    fi
    ;;
esac

# --- AVT resolution -----------------------------------------------------------
# Priority: explicit --avt > inferred from --fastq-dir (first AVT[0-9]{3,5}
# token anywhere in the path, leftmost wins) > error.
if [[ -z "$AVT" ]]; then
  AVT="$(printf '%s' "$FASTQ_DIR" | grep -oE 'AVT[0-9]{3,5}' | head -n1)"
  if [[ -n "$AVT" ]]; then
    info "AVT inferred from --fastq-dir: $AVT (override with --avt if this is wrong)"
  fi
fi
if [[ -z "$AVT" ]]; then
  err "could not infer an AVT id from --fastq-dir ('$FASTQ_DIR')"
  err "no path segment matched AVT[0-9]{3,5} (e.g. /.../AVT0551/...)"
  err "pass it explicitly with: --avt AVTXXXX"
  echo; usage; exit 1
fi
AVT="${AVT^^}"

# --- reference resolution ---------------------------------------------------
# Two tiers, selected by --species:
#   human|mouse : STAR index + GTF + gene-set label come from the table below;
#                 --star-index / --gtf optional overrides.
#   other       : fully explicit. --gtf is REQUIRED. The STAR genome is either
#                 --star-index (prebuilt) or built on demand by STAR
#                 --runMode genomeGenerate from --ref-fasta (needs write access).
for line in "${REFS[@]}"; do
  IFS='|' read -r g s gt ref <<<"$line"
  if [[ "$g" == "$SPECIES" ]]; then GTF="$gt"; GENESET="$ref"; break; fi
done

if [[ -z "$GTF" ]]; then
  # custom species: the user supplies the annotation and the gene-set label
  GTF="$GTF_ARG"
  [[ -n "$GTF" ]] || { err "--gtf is required for a non-built-in species"; echo; usage; exit 1; }
  GENESET="${REF_NAME:-custom}"
else
  # built-in species: --gtf (if given) wins over the table entry
  GTF="${GTF_ARG:-$GTF}"
fi

# --- output dir -------------------------------------------------------------
# Runs are filed under a fixed root with an explicit scRNAseq modality level
# and a per-invocation timestamp on the leaf so re-runs do not collide with a
# prior run's outputs:
#   /mnt/data/FLASHseq_output/scrnaseq/<AVT>_<species>_<YYYYMMDD_HHMM>
# The root is NOT checked for existence or writability (different host); an
# explicit --outdir overrides the whole path.
RUN_TS="$(date +%Y%m%d_%H%M)"
if [[ -z "$OUT" ]]; then
  OUT="${OUT_BASE}/scrnaseq/${AVT}_${SPECIES}_${RUN_TS}"
fi
OUT="${OUT%/}"
# The output root may not exist / be writable on the host we're invoking from
# (different server). In dry-run we never touch it; in real runs we mkdir the
# stage sub-dirs best-effort and route the run log through a real path.
if [[ "$DRY_RUN" == "1" ]]; then
  RUN_LOG="/dev/null"
else
  for d in fastqc bam featurecounts matrix logs; do
    mkdir -p "$OUT/$d" 2>/dev/null || true
  done
  RUN_LOG="$OUT/logs/pipeline_run.log"
  : > "$OUT/logs/_failed.samples"
fi
exec 3>>"$RUN_LOG"

# --- checks -----------------------------------------------------------------
# Validate the STAR genomeDir: flat Genome + SA files present and readable
# (the two non-index files every STAR genomeDir must expose).
check_genome_readable() {
  [[ -d "$STAR_INDEX" && -f "$STAR_INDEX/Genome" && -f "$STAR_INDEX/SA" ]] \
    || { err "STAR genome dir invalid: $STAR_INDEX (need flat Genome + SA files)"; return 1; }
  [[ -r "$STAR_INDEX/Genome" && -r "$STAR_INDEX/SA" ]] \
    || { err "STAR genome not readable: $STAR_INDEX"
         err "  Fix with: chmod -R o+rX $STAR_INDEX"; return 1; }
  if [[ -f "$STAR_INDEX/genomeParameters.txt" ]]; then
    local wg ws gg gs
    wg=$(awk '/genomeFileSizes/{print $2}' "$STAR_INDEX/genomeParameters.txt")
    ws=$(awk '/genomeFileSizes/{print $3}' "$STAR_INDEX/genomeParameters.txt")
    gg=$(stat -c%s "$STAR_INDEX/Genome" 2>/dev/null); gs=$(stat -c%s "$STAR_INDEX/SA" 2>/dev/null)
    if [[ -n "${wg:-}" && -n "${ws:-}" ]]; then
      if [[ "$gg" == "$wg" && "$gs" == "$ws" ]]; then
        info "STAR genome integrity: Genome+SA sizes match genomeParameters.txt"
      else
        warn "STAR genome size mismatch (want ${wg}/${ws} got ${gg}/${gs}) — proceeding"
      fi
    fi
  fi
  return 0
}

# Ensure $STAR_INDEX is a valid STAR genome, building it on demand when the
# user supplies a reference FASTA and the index dir does not (yet) contain
# the flat Genome + SA pair that STAR alignment expects.
resolve_star_index() {
  # Not needed when the STAR alignment stage is bypassed in this invocation.
  if [[ "$SKIP_STAR" -eq 1 || "$COUNT_ONLY" -eq 1 ]]; then
    info "STAR index: not needed (skip-star=$SKIP_STAR count-only=$COUNT_ONLY)"
    return 0
  fi

  # Default location for a species' index (Ensembl-116 FLASHseq tree).
  # A custom species (anything other than human|mouse) has no default index:
  # it must be mapped to a prebuilt --star-index or built from --ref-fasta, and
  # never falls back to a built-in genome (which would silently re-target the run).
  if [[ -z "$STAR_INDEX" ]]; then
    case "${SPECIES:-}" in
      human) STAR_INDEX="/home/gecf/references/STAR_references/FLASHseq/current_human_star/Homo_sapiens.GRCh38.116" ;;
      mouse) STAR_INDEX="/home/gecf/references/STAR_references/FLASHseq/current_mouse_star/Mus_musculus.GRCm39.116" ;;
      "")
        die "no species specified: pass --species <human|mouse|name>, or --star-index <dir>"
        ;;
      *)
        # custom species: only buildable if a FASTA was supplied
        [[ -n "$REF_FASTA" ]] \
          || die "species '$SPECIES' has no default STAR index: pass --star-index, or --ref-fasta (with a GTF) to build one"
        STAR_INDEX="/home/gecf/references/STAR_references/FLASHseq/custom/${SPECIES}_star"
        ;;
    esac
  fi

  # Prebuilt? Done.
  if [[ -f "$STAR_INDEX/Genome" && -f "$STAR_INDEX/SA" ]]; then
    info "STAR index: using $STAR_INDEX"
    return 0
  fi

  # Otherwise build a splice-aware STAR index from --ref-fasta + $GTF
  # (GTF already required/validated for custom species above); the build is
  # idempotent because of the presence check above.
  [[ -n "$REF_FASTA" ]] || die "STAR index not present at $STAR_INDEX — provide a prebuilt one via --star-index, or generate one by adding --ref-fasta <genome.fa>"
  [[ -f "$REF_FASTA" && -r "$REF_FASTA" ]] || die "--ref-fasta not readable: $REF_FASTA"

  local cmd=("$STAR_BIN"
      "--runMode"            genomeGenerate
      "--runThreadN"         "$THREADS_PER"
      "--genomeFastaFiles"   "$REF_FASTA"
      "--sjdbGTFfile"        "$GTF"
      "--genomeDir"          "$STAR_INDEX"
      "--sjdbOverhang"       "$SJDB_OVERHANG")

  if [[ "$DRY_RUN" == "1" ]]; then
    info "[dry-run] STAR genomeGenerate: ${cmd[*]}"
    return 0
  fi
  mkdir -p "$STAR_INDEX" 2>/dev/null || die "cannot create STAR index dir: $STAR_INDEX (check permissions)"
  info "STAR index absent — building at $STAR_INDEX (--sjdbOverhang $SJDB_OVERHANG)"
  "${cmd[@]}" >>"$RUN_LOG" 2>&1 || die "STAR genomeGenerate failed rc=$? — see $RUN_LOG"
  [[ -f "$STAR_INDEX/Genome" && -f "$STAR_INDEX/SA" ]] \
    || die "STAR genomeGenerate finished but $STAR_INDEX has no Genome/SA"
  info "STAR index built: $STAR_INDEX"
}

# Sample selection: --samples list wins (warn on any not present in
# --fastq-dir); otherwise all avid*/ subdirs of --fastq-dir.
# --exclude drops avidNNNNN sample dirs by name before either path, so they
# never enter the count matrix.
EXCL=()
if [[ -n "$EXCLUDE" ]]; then
  IFS=',' read -ra EXCL <<<"$EXCLUDE"
  for i in "${!EXCL[@]}"; do EXCL[$i]="${EXCL[$i]//[[:space:]]/}"; done
fi

# True (0) if the sample dir name is in the --exclude list.
is_excluded() {
  local n="$1" cand
  for cand in "${EXCL[@]}"; do
    [[ -n "$cand" && "$cand" == "$n" ]] && return 0
  done
  return 1
}

list_samples() {
  local n d
  local -a EXCLUDED_LOG=()
  if [[ -n "$SAMPLES" ]]; then
    local arr; IFS=',' read -ra arr <<<"$SAMPLES"
    for n in "${arr[@]}"; do
      n="${n// /}"; [[ -z "$n" ]] && continue
      if [[ -d "$FASTQ_DIR/$n" ]]; then
        if is_excluded "$n"; then EXCLUDED_LOG+=("$n"); continue; fi
        SAMPLES_LIST["$n"]=1
      else warn "sample '$n' not found under $FASTQ_DIR (skipped)"; fi
    done
  else
    for d in "$FASTQ_DIR"/avid*/; do
      [[ -d "$d" ]] || continue
      n=$(basename "$d"); [[ "$n" == avid* ]] || continue
      if is_excluded "$n"; then EXCLUDED_LOG+=("$n"); continue; fi
      SAMPLES_LIST["$n"]=1
    done
  fi
  if [[ ${#EXCLUDED_LOG[@]} -gt 0 ]]; then
    info "--exclude dropped ${#EXCLUDED_LOG[@]} sample(s): ${EXCLUDED_LOG[*]}"
    warn "excluded samples are skipped entirely (no FastQC/STAR/counts), so they never enter the count matrix"
  fi
}

# Locate the paired reads for a sample; sets R1/R2, requires both.
find_pairs() {
  local dir="$FASTQ_DIR/$ID" f
  R1=""; R2=""
  for f in "$dir"/*_R1_001.fastq.gz; do R1="$f"; break; done
  for f in "$dir"/*_R2_001.fastq.gz; do R2="$f"; break; done
  [[ -n "$R1" && -n "$R2" ]]
}

# Derive the sample label used for matrix columns / count-table files:
# the R1 FASTQ basename with everything from _S<n> onward removed.
#   lung_001_S1_R1_001.fastq.gz -> lung_001
#   SampleA_S3_R1_28_R2_90.fastq.gz -> SampleA
# Unmatched names pass through (minus the read-type suffix) so nothing
# is ever silently renamed to a bare read number.
derive_label() {
  local base="${R1##*/}"
  local name="$base"
  name="${name%.fastq.gz}"; name="${name%.fastq}"
  if [[ "$name" =~ ^(.*)_S[0-9]+_ ]]; then REPLY="${BASH_REMATCH[1]}"
  elif [[ "$name" =~ ^(.*)_[Rr][12]_ ]]; then REPLY="${BASH_REMATCH[1]}"
  else REPLY="$name"; fi
  # Collision guard: two sample dirs claiming the same label -> fall back to
  # the sample dir name (unique by construction) so no column overwrites.
  if [[ -n "${LBL_SEEN[$REPLY]:-}" && "${LBL_SEEN[$REPLY]}" != "$ID" ]]; then
    warn "$ID: label '$REPLY' already used by ${LBL_SEEN[$REPLY]} — using '$ID'"
    REPLY="$ID"
  fi
  LBL_SEEN["$REPLY"]="$ID"
  LBL="$REPLY"
}

# Locate a non-empty STAR BAM for a sample (used on re-runs to skip
# re-alignment); sets STAR_BAM.  Tries the FASTQ-derived label first
# (current naming), then the bare sample-dir name (predecessor runs).
find_star_bam() {
  local cand pref
  for pref in "${LBL:-}" "$ID"; do
    [[ -n "$pref" ]] || continue
    for cand in "$OUT/bam/${pref}"*.sortedByCoord*.bam "$OUT/bam/${pref}"*.out.bam "$OUT/bam/${pref}".bam; do
      [[ -f "$cand" && -s "$cand" ]] && { STAR_BAM="$cand"; return 0; }
    done
  done
  return 1
}

# --- stage 1: FastQC ----------------------------------------------------------
stage_fastqc() {
  local f base
  for f in "$R1" "$R2"; do
    base="$(basename "$f")"
    base="${base%.fastq.gz}"
    base="${base%.fq.gz}"
    base="${base%.fastq}"
    if [[ "$RESUME" -eq 1 ]] && ls "$OUT/fastqc/" 2>/dev/null | grep -qE "^${base}(_report|_fastqc)?\.zip$"; then
      info "$ID: FastQC present for $base (skipped)"
      continue
    fi
    if [[ "$DRY_RUN" == "1" ]]; then
      info "$ID [dry-run]: fastqc -f -o $OUT/fastqc '$f'"
      continue
    fi
    local rc=0
    fastqc -f fastq -o "$OUT/fastqc" "$f" >>"$RUN_LOG" 2>&1 || rc=$?
    if (( rc == 0 )); then info "$ID: FastQC done for $base"
    else warn "$ID: FastQC rc=$rc on $base (non-fatal)"; fi
  done
}

# --- stage 2: STAR ------------------------------------------------------------
stage_star() {
  if [[ -n "$STAR_BAM" && "$RESUME" -eq 1 ]]; then
    info "$ID: STAR BAM present (skipped): $(basename "$STAR_BAM")"
    return 0
  fi
  local STAR_CMD=("$STAR_BIN"
      "--runThreadN"         "$THREADS_PER"
      "--limitBAMsortRAM"    "20000000000"
      "--genomeLoad"         "LoadAndKeep"
      "--genomeDir"          "$STAR_INDEX"
      "--readFilesIn"        "$R1" "$R2"
      "--readFilesCommand"   "zcat"
      "--limitSjdbInsertNsj" "2000000"
      "--outFilterIntronMotifs"          "RemoveNoncanonicalUnannotated"
      "--outSAMtype"         BAM SortedByCoordinate
      "--outFileNamePrefix"  "$OUT/bam/${LBL}_")
  if [[ "$DRY_RUN" == "1" ]]; then
    info "$ID [dry-run]: ${STAR_CMD[*]}"
    return 0
  fi
  info "$ID: STAR starting (${THREADS_PER} threads)"
  "${STAR_CMD[@]}" >>"$RUN_LOG" 2>&1
  local rc=$?
  if (( rc != 0 )); then err "$ID: STAR failed rc=$rc — see $RUN_LOG"; return "$rc"; fi
  if ! find_star_bam; then
    err "$ID: STAR finished but no sorted BAM found in $OUT/bam"
    return 1
  fi
  info "$ID: STAR done -> $(basename "$STAR_BAM")"
}

# --- QC report: uniquely-mapped % from STAR Log.final.out ---
# Informational only (per Alithea 4.2.4): logged to _gates.csv and a [WARN]
# is emitted if <70%.  Samples are NOT filtered or dropped here.
log_mapping_qc() {
  local lf="$OUT/bam/${LBL}_Log.final.out"
  [[ -f "$lf" ]] || { warn "$ID: no Log.final.out — uniquely-mapped % not recorded"; return 0; }
  # STAR 2.7.11b Log.final.out format:
  #   Number of input reads |\t2287173
  #   Uniquely mapped reads number |\t2080181
  #   Uniquely mapped reads % |\t90.95%
  local pct um tot
  read -r pct um tot <<<"$(grep -E 'input reads|Uniquely mapped' "$lf" 2>/dev/null | awk -F'|' '{gsub(/%/, "", $2); gsub(/ /,"",$2)} /input reads/{t=$2} /number/{u=$2} /%/ {p=$2} END{print p, u, t}')"
  GATE_PCT="${pct:-NA}"
  if [[ "$GATE_PCT" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "$ID,$GATE_PCT,${um:-NA},${tot:-NA},$SPECIES" >>"$OUT/logs/_gates.csv"
    if awk -v p="$GATE_PCT" 'BEGIN{exit !(p+0 >= 70)}'; then
      info "$ID: ${GATE_PCT}% uniquely mapped (recorded)"
    else
      warn "$ID: ${GATE_PCT}% uniquely mapped (<70% reference; informational only)"
    fi
  else
    echo "$ID,NA,NA,NA,$SPECIES" >>"$OUT/logs/_gates.csv"
    warn "$ID: could not parse unique-mapped % from Log.final.out"
  fi
}

# --- stage 3: featureCounts ----------------------------------------------------
stage_count() {
  local out="$OUT/featurecounts/${LBL}.featureCounts.txt"
  # Resume by label; migrate a legacy table named after the sample dir so
  # the column name in the matrix follows the FASTQ-based label.
  if [[ -s "$out" && "$RESUME" -eq 1 ]]; then
    info "$ID: featureCounts present (skipped)"
    return 0
  fi
  if [[ "$RESUME" -eq 1 && "$LBL" != "$ID" && -s "$OUT/featurecounts/${ID}.featureCounts.txt" ]]; then
    info "$ID: renaming legacy count table ${ID}.featureCounts.txt -> ${LBL}.featureCounts.txt"
    mv "$OUT/featurecounts/${ID}.featureCounts.txt" "$out"
    return 0
  fi
  # Alithea counting convention (Part 4):
  #   featureCounts -T 1 -t exon -g gene_name --fracOverlap 0.25
  # -p is required: STAR emits a paired-end BAM (both mates carry mate-info);
  # -g gene_name doesn't work with Ensembl gtf. -g gene_id does the job
  # without -p featureCounts hard-errors on "Paired-end reads in a single-end library".
  local FC_CMD=(featureCounts -p -T 1 -t exon -g gene_id --fracOverlap 0.25
                -a "$GTF" -o "$out" "$STAR_BAM")
  if [[ "$DRY_RUN" == "1" ]]; then
    info "$ID [dry-run]: ${FC_CMD[*]}"
    return 0
  fi
  info "$ID: featureCounts starting"
  "${FC_CMD[@]}" >>"$RUN_LOG" 2>&1
  local rc=$?
  (( rc != 0 )) && err "$ID: featureCounts rc=$rc"
  info "$ID: featureCounts done"
  return $rc
}

# ---------------------------------------------------------------------------
# Per-sample driver: runs the full chain for one sample, aborts the sample
# at the first hard failure; other samples are unaffected.
# ---------------------------------------------------------------------------
run_sample() {
  local id="$1"
  ID="$id"
  FAILSFILE="$OUT/logs/_failed.samples"
  { echo; echo "==================== $ID  start: $(date -Is) ===================="; } >&3
  if ! find_pairs; then
    { err "$ID: no paired *_R1_001.fastq.gz / *_R2_001.fastq.gz in $FASTQ_DIR/$ID — skipped"; } >&3
    echo "$id" >>"$FAILSFILE"
    return 1
  fi
  derive_label
  info "$ID: matrix column / count-table name: $LBL"
  if [[ "$DO_FASTQC" -eq 1 && "$COUNT_ONLY" -eq 0 ]]; then stage_fastqc; fi

  STAR_BAM=""
  if [[ "$SKIP_STAR" -eq 0 && "$COUNT_ONLY" -eq 0 ]]; then
    find_star_bam || STAR_BAM=""
    stage_star || { err "$ID: STAR failed — see $RUN_LOG"; echo "$id" >>"$FAILSFILE"; return 1; }
    log_mapping_qc
  else
    if ! find_star_bam; then
      if [[ "$DRY_RUN" == "1" ]]; then
        info "$ID [dry-run]: (skip-star/count-only — would reuse an existing BAM in $OUT/bam)"
      else
        { err "$ID: --skip-star/--count-only but no BAM under $OUT/bam — skipped"; } >&3
        echo "$id" >>"$FAILSFILE"
        return 1
      fi
      STAR_BAM="$OUT/bam/${LBL}_Aligned.sortedByCoord.out.bam"
    fi
  fi

  if stage_count; then :; else err "$ID: featureCounts failed — see $RUN_LOG"; echo "$id" >>"$FAILSFILE"; return 1; fi
  { echo; echo "==================== $ID   done:  $(date -Is) ===================="; } >&3
  return 0
}

# ---------------------------------------------------------------------------
# Step 6: writes _merge.R and runs it to combine the per-sample
# featureCounts tables into one genes x samples matrix.
# ---------------------------------------------------------------------------
write_merge_r() {
  cat > "$OUT/logs/_merge.R" <<'R'
suppressMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
})
args     <- commandArgs(trailingOnly = TRUE)
countdir <- args[1]; outdir <- args[2]

files <- list.files(countdir, pattern = "\\.featureCounts\\.txt$", full.names = TRUE)
files <- files[grepl("\\.featureCounts\\.txt$", basename(files))]    # drop .S1_annotated_summary
files <- sort(files)
if (length(files) == 0) stop("no featureCounts files in ", countdir)
ids <- sub("\\.featureCounts\\.txt$", "", basename(files))   # strip suffix, keep prefix

read_one <- function(path, id) {
  df <- read_tsv(path, show_col_types = FALSE, comment = "#")
  if (ncol(df) < 2) stop("featureCounts file too narrow: ", path)
  df <- df[, c(1, ncol(df)), drop = FALSE]
  colnames(df) <- c("geneID", "count")
  tibble(geneID = as.character(df[["geneID"]]),
         count  = as.integer(df[["count"]]),
         sample = id)
}
long <- bind_rows(lapply(seq_along(files), function(i) read_one(files[i], ids[i])))

geneids <- unique(long$geneID)                       # union alignment by geneID
wide <- long %>%
  mutate(geneID = factor(geneID, levels = geneids)) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0L) %>%
  as_tibble()
wide[["geneID"]] <- as.character(wide[["geneID"]])
order_cols <- c("geneID", sort(setdiff(colnames(wide), "geneID")))
wide <- wide[, order_cols]

# --- Per-cell genes detected (two count thresholds: >=1 and >=5) -------------
# Informational only, same stance as _gates.csv: FLASH-seq is sparse, so low
# detection counts are expected and are never a hard gate. For each cell
# (column) count how many distinct genes meet each threshold.
cells     <- setdiff(colnames(wide), "geneID")
mat       <- as.matrix(wide[, cells, drop = FALSE])        # gene x cell, numeric
gpc       <- tibble(cell = colnames(mat),
                    genes_geq1 = as.integer(colSums(mat >= 1L)),
                    genes_geq5 = as.integer(colSums(mat >= 5L)))
gpc_out   <- file.path(outdir, "genes_per_cell.csv")
write.csv(gpc, gpc_out, row.names = FALSE)
cat(sprintf("GENES/CELL: %d cells -> %s\n", length(cells), gpc_out))
for (nm in c("genes_geq1", "genes_geq5")) {
  v <- gpc[[nm]]
  cat(sprintf("  %-10s min=%d  median=%d  max=%d\n", nm, min(v), as.integer(round(median(v))), max(v)))
}

out <- file.path(outdir, "counts_matrix.tsv")
write_tsv(wide, out)
cat(sprintf("MERGE: %d samples x %d genes -> %s\n", length(ids), nrow(wide) - 1, out))
R
}

merge_matrix() {
  write_merge_r
  Rscript "$OUT/logs/_merge.R" "$OUT/featurecounts" "$OUT/matrix"
}

# ---------------------------------------------------------------------------
# Step 6b (best-effort): derive a sample -> group metadata CSV straight from
# the merged matrix's header row (first column geneID, other columns = sample
# names; group = leading token, e.g. lung_001 -> lung). Writes
# matrix/counts_matrix_metadata.csv. Any failure (no python3, unreadable
# matrix, ...) is a WARNING only — the pipeline must never stop here.
# ---------------------------------------------------------------------------
write_metadata_py() {
  cat > "$OUT/logs/_make_metadata.py" <<'PY'
import csv
import os
import sys
from collections import Counter

def group_of(name):
    """Infer the sample group: the name minus a trailing '_<digits>' index."""
    parts = name.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return name

def load_matrix_labels(path):
    """Sample names = header row entries after the first column (geneID)."""
    with open(path, "r", newline="") as f:
        first = f.readline()
    if not first:
        raise ValueError("empty count matrix: %s" % path)
    first = first.rstrip("\n")
    parts = first.split("\t") if "\t" in first else first.split(",")
    labels = [p.strip() for p in parts[1:] if p.strip()]
    if not labels:
        raise ValueError("no sample columns found in header of: %s" % path)
    return labels

def main():
    count_matrix = sys.argv[1]
    stem, _ = os.path.splitext(count_matrix)
    if stem.endswith("_metadata"):
        stem = stem[: -len("_metadata")]
    out_csv = stem + "_metadata.csv"
    outdir = os.path.dirname(os.path.abspath(out_csv))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    labels = load_matrix_labels(count_matrix)
    pairs = {name: group_of(name) for name in labels}    # stable order
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "sample"])
        for n in labels:
            w.writerow([n, pairs[n]])
    counts = Counter(pairs.values())
    print("METADATA: %d samples -> %s" % (len(labels), out_csv))
    print("  group counts: %s" % dict(sorted(counts.items())))

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:      # report + non-zero; pipeline treats as warning
        sys.stderr.write("metadata CSV generation failed: %s\n" % exc)
        sys.exit(1)
PY
}

make_metadata() {
  local matrix="$OUT/matrix/counts_matrix.tsv"
  if [[ ! -s "$matrix" ]]; then
    warn "metadata step skipped: $matrix is missing or empty"
    return 0
  fi
  if ! command -v python3 >/dev/null 2>&1; then
    warn "metadata step skipped: python3 not on PATH"
    return 0
  fi
  write_metadata_py
  info "Writing sample -> group metadata CSV (best-effort) ..."
  local rc=0
  python3 "$OUT/logs/_make_metadata.py" "$matrix" 1>&3 2>&3 || rc=$?
  if (( rc != 0 )); then
    warn "metadata CSV was NOT generated — continuing (see $RUN_LOG)"
  else
    info "metadata : $OUT/matrix/counts_matrix_metadata.csv"
  fi
  return 0
}

# ===========================================================================
# MAIN
# ===========================================================================
info "=== FlashSeq pipeline (Alithea 2026 Part-4) ==="
info "species    : $SPECIES ($GENESET)"
info "fastq root : $FASTQ_DIR"
info "STAR index : $STAR_INDEX   (bin: $STAR_BIN)"
info "GTF        : $GTF"
info "out dir    : $OUT"
info "parallel   : ${MAX_JOBS} samples x ${THREADS_PER} threads"
info "stages     : fastqc=$DO_FASTQC  skip-star=$SKIP_STAR  count-only=$COUNT_ONLY  resume=$RESUME  dry-run=$DRY_RUN"
{
  printf '=== run %s ===\n' "$(date -Is)"
  printf 'species=%s(%s) avt=%s out=%s parallel=%sx%s\n' "$SPECIES" "$GENESET" "$AVT" "$OUT" "$MAX_JOBS" "$THREADS_PER"
  printf '%-21s: %s\n' 'flashseq_pipeline.sh' "$SCRIPT_VERSION"
  echo "--- tools ---"
  { command -v "$STAR_BIN" >/dev/null && "$STAR_BIN" --version 2>&1 | head -1; } || STAR --version 2>&1 | head -1
  featureCounts -v 2>&1 | head -1
  Rscript -e 'suppressMessages({library(dplyr);library(readr);library(stringr);library(tibble);library(tidyr)});
              cat("R", R.version.string, "| dplyr", as.character(packageVersion("dplyr")),
                  "| readr", as.character(packageVersion("readr")),
                  "| stringr", as.character(packageVersion("stringr")),
                  "| tibble", as.character(packageVersion("tibble")),
                  "| tidyr", as.character(packageVersion("tidyr")), "\n")' 2>&1 | tail -1
} 1>&3 2>&1

[[ -d "$FASTQ_DIR" ]] || die "fastq dir not found: $FASTQ_DIR"
[[ -f "$GTF" && -r "$GTF" ]] || die "GTF not readable: $GTF"

# Resolve the STAR index (use existing, or build from --ref-fasta).
# No-op when the STAR stage is bypassed (--skip-star / --count-only).
resolve_star_index

list_samples
[[ ${#SAMPLES_LIST[@]} -gt 0 ]] || die "no samples found (check --fastq-dir / --samples)"
N=${#SAMPLES_LIST[@]}
info "samples    : $N"

if [[ "$DRY_RUN" == "1" ]]; then
  check_genome_readable || warn "genome check failed (dry-run continuing regardless)"
  for n in $(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort); do
    run_sample "$n"
  done
  echo
  echo "=========================== DRY RUN complete — nothing was launched ==========================="
  exit 0
fi

# Only enforce the genome-dir check when the STAR alignment stage actually
# runs; --skip-star / --count-only bypass alignment entirely, so a missing
# or empty $STAR_INDEX must not abort the count/resume path.
if [[ "$SKIP_STAR" -eq 0 && "$COUNT_ONLY" -eq 0 ]]; then
  check_genome_readable || die "STAR genome not readable — nothing was launched (chmod hint printed above)."
fi
if [[ -S /dev/shm/shm  || -d /dev/shm ]]; then
  df -h /dev/shm | tail -1 | awk '{printf "shm: %s/%s\n", $2, $3}' >&3
  shms=$(df -k /dev/shm 2>/dev/null | awk 'NR==2{print $4}')
  if [[ -n "$shms" ]] && (( shms < 30000000 )); then
    warn "/dev/shm free = $((shms/1024)) MB < ~30 GB (STAR genome needs this to LoadAndKeep-share)"
  fi
fi

echo "sample,unique_mapped_pct,uni_reads,total_reads,species" >"$OUT/logs/_gates.csv"

START=$SECONDS
run_in_batch() {
  local id p
  for id in $(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort); do
    run_sample "$id" &
    PIDS+=("$!")
    if (( ${#PIDS[@]} >= MAX_JOBS )); then
      # Throttle: block until the OLDEST in-flight job finishes (FIFO).
      # Portable: works on bash < 4.3 (no `wait -n`) and < 5.1 (no `wait -p`).
      # `wait -n -p wpid` + set -u crashed here on older bash ("wpid: unbound
      # variable") and killed the whole run at the first throttle point.
      wait "${PIDS[0]}" 2>/dev/null || true
      PIDS=("${PIDS[@]:1}")
    fi
  done
  for p in "${PIDS[@]}"; do wait "$p" 2>/dev/null || true; done
  return 0
}
run_in_batch

ELAPSED=$(( SECONDS - START ))

# --- merge matrix -------------------------------------------------------------
if ls "$OUT/featurecounts"/*.featureCounts.txt >/dev/null 2>&1; then
  info "Merging gene x sample matrix ..."
  if merge_matrix 1>&3 2>&1; then
    M="$OUT/matrix/counts_matrix.tsv"
    info "Matrix : $M"
    info "shape  : $(awk -F'\t' 'NR==1{n=NF} END{print (NR-1)" genes x "(n-1)" samples"}' "$M")"
    # plausibility: distinct gene_id seen in featureCounts output vs matrix rows
    one_fc=$(ls -1 "$OUT/featurecounts" | grep -E '\.featureCounts\.txt$' | head -1)
    if [[ -n "$one_fc" && -f "$OUT/featurecounts/$one_fc" ]]; then
      fcgenes=$(awk -F'\t' 'NR>1{print $3}' "$OUT/featurecounts/$one_fc" | sort -u | wc -l)
      mgenes=$(( $(wc -l < "$M") - 1 ))
      if (( mgenes < fcgenes )); then
        info "note: matrix has $mgenes genes, < the $fcgenes distinct gene_id in a featureCounts output row set (featureCounts drops unassigned reads; confirm output not truncated)"
      else
        info "gene count plausible (matrix rows=$mgenes)"
      fi
    fi
    empty=$(awk -F'\t' 'NR==1{for(i=2;i<=NF;i++)h[i]=i; next} {for(i=2;i<=NF;i++) if($i+0>0)c[h[i]]=1} END{n=0; for(i in h) if(!(i in c)) n++; print n}' "$M")
    info "empty sample columns: $empty"
  else
    err "Matrix merge failed — see $OUT/logs/_merge.R and $RUN_LOG"
    exit 2
  fi
else
  warn "no featureCounts outputs found — matrix merge skipped"
fi

# Best-effort metadata from the (possibly pre-existing) matrix. Never fatal:
# make_metadata always returns 0 and downgrades its own failures to warnings.
make_metadata

# --- final summary --------------------------------------------------------------
below=$(awk -F',' 'NR>1 && $2 ~ /^[0-9.]+$/ && ($2+0) < 70 {print $1}' "$OUT/logs/_gates.csv" 2>/dev/null | tr '\n' ' ')
failed=""
[[ -f "$OUT/logs/_failed.samples" ]] && failed="$(sort -u "$OUT/logs/_failed.samples" 2>/dev/null | tr '\n' ' ')"
# Mirror the final summary to both the terminal and pipeline_run.log.
SUMMARY="$(
  {
    echo
    info "=== SUMMARY ==="
    info "elapsed : ${ELAPSED}s"
    info "species : $SPECIES"
    info "out dir : $OUT"
    info "samples : $N"
    info "failed/skipped : ${failed:-none}"
    info "uniquely mapped <70% (informational): ${below:-none}"
    info "matrix  : $OUT/matrix/counts_matrix.tsv"
    info "qc      : $OUT/matrix/genes_per_cell.csv  (genes/cell >=1 and >=5)"
    info "logs    : $OUT/logs/  (pipeline_run.log, _gates.csv, _merge.R)"
    info "DONE"
  } 2>&1
)"
if [[ -s "$RUN_LOG" ]]; then
  printf '%s\n' "$SUMMARY" | tee -a "$RUN_LOG"
else
  printf '%s\n' "$SUMMARY"
fi
