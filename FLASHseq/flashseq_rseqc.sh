#!/usr/bin/env bash
#=============================================================================
# flashseq_rseqc.sh -- RSeQC QC for FLASH-seq STAR BAMs
#                      (split out of flashseq_pipeline.sh, 2026-09-17)
#
# geneBody_coverage (--steps gbcov), TWO modes:
#   MERGED (default): PARALLEL per-sample geneBody_coverage.py --skip-plot
#                     calls (fanned out to --max-jobs); the per-sample tables
#                     are concatenated into a single combined coverage table,
#                     and the combined plots are then redrawn from that table
#                     with gbcov_redraw.py (no BAM re-read in the plot phase).
#                     Row names/order follow RSeQC (valid_name(bam stem),
#                     sorted members).
#                     Flat (all-zero) profiles do not abort the run -- the
#                     redraw step renders them as zero curves and flags the
#                     sample in _failed.samples.
#                     Artifacts: rseqc/geneBodyCombined.* + per-sample
#                     intermediates rseqc/geneBodyCombined.per_sample/ +
#                     logs/gbcov_parallel_<LBL>.rseqc.log.
#   PER-SAMPLE (--per-sample): one call per sample
#                              -- per-sample table + per-sample plots
#                              (heatMap + curves) under rseqc/<LBL>_gbcov.*,
#                              with per-sample <LBL>.rseqc.log entries.
# read_distribution (--steps rdist) is ALWAYS per-sample in either mode:
#                       per-sample <LBL>.readDistribution.txt
#
# This script is the ONLY thing that needs STAR BAMs (no alignment, no
# counting): point --outdir at an existing flashseq_pipeline.sh run and
# it picks up every bam/*.out.bam already there.
#
# Default REFERENCES (--species human / mouse) are PER-TOOL and DIFFERENT:
#   geneBody_coverage  -> the HOUSEKEEPING-GENE BED (12-column gene-model BED
#                         for RSeQC):
#       human  -> .../current_human_star/hg38.HouseKeepingGenes.bed
#       mouse  -> .../current_mouse_star/mm39.HouseKeepingGenes.bed
#   read_distribution  -> the GENOME-WIDE gene model, BAKED to BED12 with
#                         the UCSC genePred chain (gtfToGenePred |
#                         genePredToBed) from the species GTF (REFS table).
#                         NEVER the housekeeping set -- read_distribution
#                         exons + introns to split reads into CDS/5'UTR/
#                         3'UTR/intron/intergenic.
# Override with --refbed PATH (ready housekeeping BED) and/or --gtf PATH
# (genome-wide gene model; also the gbcov fallback when no housekeeping
# BED exists for the species).
#
# Artifacts -- standard output layout:
#   $OUT/rseqc/<gtf basename>.bed12
#   $OUT/rseqc/<LBL>_gbcov.*
#   $OUT/rseqc/<LBL>.readDistribution.txt
#   $OUT/logs/<LBL>.rseqc.log
#   $OUT/logs/_failed.samples
#
# Resume-safe: a sample is skipped when BOTH
#   rseqc/<LBL>_gbcov.geneBodyCoverage.txt and
#   rseqc/<LBL>.readDistribution.txt
# are non-empty.
#
# RSeQC 5.0.5 CLI (verified 2026-09-15, see
#   ~/.hermes/skills/research/flashseq-rnaseq/references/
#   rseqc-5.0.5-integration.md):
#   geneBody_coverage.py  -i BAM -r BED12 -o PREFIX [--skip-plot]
#     reference = the HOUSEKEEPING-GENE BED (per-species default, --refbed)
#   read_distribution.py  -i BAM -r BED12        (report printed to stdout)
#     reference = the GENOME-WIDE gene-model BED12 baked from the species GTF
#     (--gtf). The two tools use DIFFERENT references on purpose: gbcov is a
#     gene-body QC over the housekeeping set; rdist needs exons+introns to
#     bucket CDS/5'UTR/3'UTR/intron, so it always runs against the whole genome.
#   * both require an INDEXED BAM (.bam.bai present; auto-created with
#     samtools index if missing)
#   * reference MUST be a 12-column gene-model BED (BED12). Per-species gbcov
#     default = the housekeeping BED in its STAR-references dir (REFBEDS);
#     --refbed overrides it; --gtf bakes a BED12 from a GTF via the UCSC
#     genePred chain (gtfToGenePred -> genePredToBed) for rdist (.gz
#     references are decompressed first -- UCSC tools don't auto-inflate them).
#     A non-12-column reference is WARNed about because geneBody_coverage.py
#     skips such rows.
#=============================================================================
set -uo pipefail

# --- defaults ---------------------------------------------------------------
SCRIPT_VERSION="1.4.0"   # Bump on any change that alters outputs/behaviour.
SPECIES=""
GTF_ARG=""
OUT=""
SAMPLES=""
MAX_JOBS=4
DRY_RUN=0
STEPS="gbcov,rdist"   # default: both
MERGED=1              # gbcov default: COMBINED output (one combined table +
                      # one combined plot) covering all selected samples,
                      # built in parallel (per-sample --skip-plot table
                      # calls + concatenated table + gbcov_redraw.py).
                      # --per-sample (MERGED=0): one call per sample, with
                      # its own per-sample plots.
                      # read_distribution is ALWAYS per-sample.
GBCOV_FORMAT="pdf"    # MERGED combined-plot format (gbcov_redraw.py --format)
                      # = pdf|png|jpeg|html. html = self-contained interactive
                      # pages via plotly (no Rscript; the 300+-sample legend
                      # stays readable: hover + clickable legend). Per-sample
                      # (--per-sample) plots are always PDF from RSeQC itself
                      # and ignore this value.
# table->plot (no BAMs). Sits in the SAME FOLDER as this script, so
# resolve it relative to the script's own location (not the caller's CWD).
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" 2>/dev/null && pwd -P)"
GBCOV_REDRAW_PY="${SCRIPT_DIR:-.}/gbcov_redraw.py"

# Built-in reference BEDs -- HOUSEKEEPING-GENE set used by geneBody_coverage:
#   <species>|<ready reference BED path>
# The lab keeps a curated gene-region BED per species in its STAR-references
# dir; geneBody_coverage runs against it directly. If the file is
# missing or unreadable, gbcov falls back to the genome-wide BED12 baked
# from the built-in species GTF (a warning is emitted).
REFBEDS=(
  "human|/home/gecf/references/STAR_references/FLASHseq/current_human_star/hg38.HouseKeepingGenes.bed"
  "mouse|/home/gecf/references/STAR_references/FLASHseq/current_mouse_star/mm39.HouseKeepingGenes.bed"
)

# Built-in GTF references (kept in sync with flashseq_pipeline.sh) -- the
# GENOME-WIDE gene model, used to BAKE the read_distribution BED12 (UCSC
# genePred chain), or via --gtf as an override:
#   <species>|<gtf path>|<gene-set label>
REFS=(
  "human|/home/gecf/references/STAR_references/FLASHseq/current_human_star/Homo_sapiens.GRCh38.116.gtf|GRCh38.116"
  "mouse|/home/gecf/references/STAR_references/FLASHseq/current_mouse_star/Mus_musculus.GRCm39.116.gtf|GRCm39.116"
)

BEDREF_ARG=""   # --refbed: explicit ready housekeeping-gene BED for GBCOV

# Per-TOOL references (set during resolution below) -- they are DELIBERATELY
# separate because the two tools need different gene models:
#   geneBody_coverage  -> the HOUSEKEEPING-GENE set (GBCOV_REFBED)
#                         = REFBEDS per-species default, or --refbed override.
#                         Fallback if neither exists: genome-wide BED12 baked
#                         from GBCOV_FALLBACK_GTF (warned about).
#   read_distribution  -> the FULL-GENOME gene model (RDIST_BED12)
#                         = REFS per-species GTF (or --gtf) BAKED to BED12.
#                         NEVER the housekeeping set: read_distribution needs
#                         exons + introns to split CDS/5'UTR/3'UTR/intron.
#                          (RDIST_GTF = the GTF it was baked from.)
GBCOV_REFBED="";    GBCOV_FALLBACK_GTF=""
RDIST_BED12="";     RDIST_GTF=""
BAKED_BED=""        # set by bake_bed12 (last GTF->BED12 conversion)
# BAKED_BEDS: set of genome-wide BED12 intermediates produced (or reused) by
# bake_bed12 THIS run, so we can remove them after use (post-run cleanup).
# A ready housekeeping BED (GBCOV_REFBED) is NOT baked and is left untouched.
BAKED_BEDS=()
# GBCOV_BED12: the final reference actually passed to geneBody_coverage
#               (set at runtime; equals the ready housekeeping BED, or a
#               baked genome-wide BED12 if no ready housekeeping BED exists).
GBCOV_BED12=""
declare -A _BAKED_DRYRUN_ANNOUNCED=()   # path -> "printed" ; de-dups the dry-run log when both tools share the same GTF
RUN_TS=""
RUN_LOG=""
ID=""
LABEL=""
BAM=""
PIDS=()
declare -A SAMPLES_LIST

# ---------------------------------------------------------------------------
# Log helpers. Timestamped + tagged: [INFO]/[WARN] progress/advisory,
# [ERROR] + die() = aborted (mirrors flashseq_pipeline.sh).
# ---------------------------------------------------------------------------
# Append a line to RUN_LOG if it points at a real file (guard: RUN_LOG is
# unset during early validation; /dev/null during --dry-run).
_log_to_file() {
  if [[ -n "${RUN_LOG:-}" && "${RUN_LOG}" != "/dev/null" ]]; then
    printf '%s\n' "$1" >>"$RUN_LOG"
  fi
}
info() { local l; l="$(date +%H:%M:%S) [INFO]   $1"; echo "$l"; _log_to_file "$l"; }
warn() { local l; l="$(date +%H:%M:%S) [WARN]   $1"; echo "$l"; _log_to_file "$l"; }
err()  { local l; l="$(date +%H:%M:%S) [ERROR]  $1"; echo "$l"; _log_to_file "$l"; }
die()  { err "$1"; exit 1; }

# log_tool_versions -- capture the version of every external tool this script
# depends on, to the console (live + dry-run preview) and the run log (real
# runs). Runs ONCE at the start of the flow, before any sample work launches,
# so the provenance of every run is recorded. Each tool self-reports via its
# --version flag:
#   geneBody_coverage.py / read_distribution.py --version   (RSeQC)
#   gtfToGenePred / genePredToBed                    (UCSC; no --version flag --
#                                                     self-identified by no-arg usage)
#   samtools --version                    (optional -- only for .bam.bai)
#   Rscript (drives geneBody_coverage plotting)            (optional)
log_tool_versions() {
  local v lbl
  info "=== tool versions (provenance) ==="
  lbl=$(printf '%-21s: ' 'flashseq_rseqc.sh'); info "${lbl}${SCRIPT_VERSION}"
  # RSeQC python tools + UCSC genePred chain (required); samtools + R
  # (gbcov plotting) are optional here. Pad each label to 21 chars so the
  # colon + value align in one column no matter the label length.
  # gtfToGenePred/genePredToBed don't take --version (they print usage to
  # stderr), so capture the name line from a no-arg invocation instead --
  # proves they are present + runnable, which is what matters here.
  for t in geneBody_coverage.py read_distribution.py gtfToGenePred genePredToBed samtools; do
    if command -v "$t" >/dev/null 2>&1; then
      case "$t" in
        gtfToGenePred|genePredToBed) v=$("$t" 2>&1 | head -1) ;;
        *) v=$("$t" --version 2>&1 | head -1) ;;
      esac
      lbl=$(printf '%-21s: ' "$t"); info "${lbl}${v}"
    else
      lbl=$(printf '%-21s: ' "$t"); info "${lbl}not in PATH"
    fi
  done
  if command -v Rscript >/dev/null 2>&1; then
    v=$(Rscript -e 'cat(R.version.string)' 2>/dev/null | sed 's/^R version //')
    lbl=$(printf '%-21s: ' 'R (gbcov plotting)'); info "${lbl}${v:-unknown}"
  else
    lbl=$(printf '%-21s: ' 'R (gbcov plotting)'); info "${lbl}not in PATH"
  fi
  if command -v python3 >/dev/null 2>&1; then
    v=$(python3 -c 'import plotly;print(plotly.__version__)' 2>/dev/null)
    lbl=$(printf '%-21s: ' 'plotly (gbcov html)'); info "${lbl}${v:-not installed}"
  else
    lbl=$(printf '%-21s: ' 'plotly (gbcov html)'); info "${lbl}python3 not in PATH"
  fi
}

usage() {
cat <<'USAGE'
Per-sample RSeQC QC (geneBody_coverage + read_distribution) for
flashseq_pipeline.sh BAMs. Standalone: no FASTQ, no STAR, no counting --
just re-runs the QC stage on existing STAR BAMs.

REQUIRED
  --species   human | mouse | <name>: per-tool defaults
                (see "per-species defaults" block at the top of this file):
                  geneBody_coverage  -> HOUSEKEEPING-GENE BED
                                        (hg38 / mm39 .HouseKeepingGenes.bed)
                  read_distribution  -> GENOME-WIDE gene-model BED12 baked
                                        from the species GTF (NEVER the
                                        housekeeping set -- rdist needs
                                        exons + introns to bucket
                                        CDS/5'UTR/3'UTR/intron)
  --outdir    PATH   an existing flashseq_pipeline.sh run directory
                     (must contain bam/ and logs/).

OPTIONAL
  --steps L,L      which RSeQC tool(s) to run:
                       gbcov  = geneBody_coverage.py
                       rdist  = read_distribution.py
                       both   = both (default)
                   any comma-separated subset: both | gbcov,rdist |
                   gbcov | rdist
  --refbed PATH        ready-made reference BED for geneBody_coverage
                       (HOUSEKEEPING set), used as-is (overrides the
                       per-species default). RSeQC needs a 12-column gene-model
                       BED (BED12); a bed/bedGraph region list would be skipped.
                       Does NOT change the read_distribution reference.
  --gtf PATH           gene model GTF used by read_distribution (baked to
                       BED12 via UCSC genePred chain: gtfToGenePred |
                       genePredToBed); used for custom species or to
                       override the per-species default.
  --samples a,b        process only these sample labels (from bam/ names)
  --max-jobs N         samples in parallel (default 4)
  --per-sample         geneBody_coverage: PER-SAMPLE mode -- one
                       geneBody_coverage.py call per sample, each with its
                       own table + plots + log entry under rseqc/.
                       (Default = MERGED mode: one combined table for all
                       selected samples + two combined plots (heatMap +
                       curves), built in parallel.)
  --gbcov-format FMT   FMT (pdf|png|jpeg|html, default pdf) for the MERGED
                       combined plot, rendered by gbcov_redraw.py.
                       html = self-contained interactive browser pages via
                       plotly (no Rscript; the 300+-sample legend stays
                       readable: hover a curve, click a legend entry to
                       hide/show, zoom/pan). Affects the MERGED combined plot
                       only; --per-sample plots are always PDF (RSeQC's own).
  --dry-run            print planned commands, execute nothing
  --version            print the script version and exit

EXAMPLES
  # all samples of an existing run, both tools, MERGED gbcov (default)
  bash flashseq_rseqc.sh --species human --outdir /path/to/AVT0551_human_20260910_0930
  # per-sample mode (one geneBody_coverage.py call + plots per sample)
  bash flashseq_rseqc.sh --species human --outdir <run dir> --per-sample
  # one sample (merged pool of one = that one sample, same result as per-sample)
  bash flashseq_rseqc.sh --species human --outdir <run dir> --samples lung_052
  # geneBody_coverage only / read_distribution only
  bash flashseq_rseqc.sh --species human --outdir <run dir> --steps gbcov
  bash flashseq_rseqc.sh --species human --outdir <run dir> --steps rdist
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species)    SPECIES="$2"; shift 2;;
    --gtf)        GTF_ARG="$2"; shift 2;;
    --refbed)     BEDREF_ARG="$2"; shift 2;;
    --outdir)     OUT="$2"; shift 2;;
    --samples)    SAMPLES="$2"; shift 2;;
    --max-jobs)   MAX_JOBS="$2"; shift 2;;
    --dry-run)    DRY_RUN=1; shift;;
    --steps)      STEPS="$2"; shift 2;;
    --gbcov-format) GBCOV_FORMAT="$2"; shift 2;;
    --per-sample) MERGED=0; shift;;
    --version)    echo "flashseq_rseqc.sh $SCRIPT_VERSION"; exit 0;;
    -h|--help)    usage; exit 0;;
    *) die "Unknown option '$1' (see --help)";;
  esac
done

# --- validate --gbcov-format -------------------------------------------------
# Fail fast: these are the only values gbcov_redraw.py --format accepts, so a
# typo should not cost a full parallel batch that dies at the redraw step.
case "${GBCOV_FORMAT,,}" in
  pdf|png|jpeg|html)  GBCOV_FORMAT="${GBCOV_FORMAT,,}" ;;
  *) die "Invalid --gbcov-format '$GBCOV_FORMAT' (use: pdf | png | jpeg | html)" ;;
esac

# --- parse & validate --steps (which RSeQC tool(s) to run) ------------------
# Accepts "both", "gbcov", "rdist", or any comma combination of the two.
# NOTE: `IFS=',' read -r single_var` never tokenizes (bash hands the last var the
# whole remainder) -- so split into an ARRAY first, then validate per token.
RUN_GBCOV=0; RUN_RDIST=0
_norm=""
IFS=',' read -r -a _toks <<< "$STEPS"
for _tok in "${_toks[@]}"; do
  _tok="${_tok//[[:space:]]/}"            # trim
  case "${_tok,,}" in
    "" )  ;;                               # ignore empties from e.g. "gbcov, ,rdist"
    both) _norm="gbcov,rdist" ;;
    gbcov|genebody|genebodycoverage|1)  _norm="${_norm:+$_norm,}gbcov" ;;
    rdist|readdistribution|2)           _norm="${_norm:+$_norm,}rdist" ;;
    *) die "Invalid --steps value '${_tok}' (use: both | gbcov | rdist | gbcov,rdist)" ;;
  esac
done

# de-dupe in a stable order (gbcov first)
for _t in $(printf '%s' "$_norm" | tr ',' '\n' | awk 'NF' | awk '!s[$0]++'); do
  case "$_t" in gbcov) RUN_GBCOV=1;; rdist) RUN_RDIST=1;; esac
done
if [ "$RUN_GBCOV" -eq 0 ] && [ "$RUN_RDIST" -eq 0 ]; then
  die "--steps must produce at least one tool (got '$STEPS')"
fi
STEPS_LABEL="gbcov"
[ "$RUN_GBCOV" -eq 1 ] && [ "$RUN_RDIST" -eq 1 ] && STEPS_LABEL="gbcov+rdist"
[ "$RUN_RDIST" -eq 1 ] && [ "$RUN_GBCOV" -eq 0 ] && STEPS_LABEL="rdist"

[[ -n "$SPECIES" ]] || { err "--species is required (human|mouse|<custom>)"; echo; usage; exit 1; }
[[ -n "$OUT"     ]] || { err "--outdir is required (the flashseq_pipeline.sh run directory)"; echo; usage; exit 1; }

# --- tool resolution --------------------------------------------------------
# RSeQC 5.0.5 entry points must be on PATH (or the known envs get prepended).
# Rscript for plots is optional; samtools/UCSC-genePred/gzip are required-ish.
_prepend_if_present() {
  local d="$1"
  [[ -d "$d" ]] || return 0
  case ":$PATH:" in *":$d:"*) : ;; *) export PATH="$d:$PATH" ;; esac
}
if command -v geneBody_coverage.py >/dev/null 2>&1 \
 && command -v read_distribution.py  >/dev/null 2>&1 \
 && command -v gtfToGenePred         >/dev/null 2>&1 \
 && command -v genePredToBed         >/dev/null 2>&1; then
  :
else
  if   [[ -n "${CONDA_PREFIX:-}" ]]; then _prepend_if_present "$CONDA_PREFIX/bin"; fi
  _prepend_if_present "/home/hermes/.conda/envs/flashseq/bin"        # drafting box
  _prepend_if_present "/home/gecf/softwares/anaconda/envs/flashseq/bin"  # ptegsrv3
fi
# UCSC genePred chain (gtfToGenePred + genePredToBed): prepend a hand
# installation if present (harmless -- _prepend_if_present skips an already
# active dir). Covers envs built before the ucsc-* blocks were added. See
# flashseq-env.yaml.
for _ucsc_dir in "${UCSC_TOOLS_DIR:-}" "$HOME/.local/bin/ucsc"; do
  [[ -n "$_ucsc_dir" ]] && _prepend_if_present "$_ucsc_dir"
done
unset _ucsc_dir
# conda_env_hint -- shared guidance when an RSeQC-env tool is missing.
# The user should be in the 'flashseq' conda env; if it doesn't exist yet,
# build it from the env file that ships next to this script.
# Defined BEFORE its first caller below (die ... conda_env_hint) since bash
# resolves functions at call time.
conda_env_hint() {
  # Print guidance telling the user to activate (or create) the 'flashseq'
  # conda env, which provides geneBody_coverage.py / read_distribution.py /
  # gtfToGenePred / genePredToBed (UCSC chain). The caller prints the specific
  # missing tool on the [ERROR] line.
  cat <<HINT
You are probably not in the 'flashseq' conda env. Fix:

  1. Activate it:
       conda activate flashseq

  2. If that fails ('Channels not configured' / 'prefix not found'), create
     it first from the env file next to this script:
       conda env create -f /home/hermes/flashseq-env.yaml
       conda activate flashseq

Then re-run this script.
HINT
}
_missing=""
command -v geneBody_coverage.py >/dev/null 2>&1 || _missing="$_missing geneBody_coverage.py"
command -v read_distribution.py >/dev/null 2>&1 || _missing="$_missing read_distribution.py"
command -v gtfToGenePred        >/dev/null 2>&1 || _missing="$_missing gtfToGenePred"
command -v genePredToBed        >/dev/null 2>&1 || _missing="$_missing genePredToBed"
[[ -n "$_missing" ]] && {
  err "RSeQC tool(s) not found on PATH:$_missing"
  conda_env_hint
  exit 1
}
command -v Rscript   >/dev/null 2>&1 || warn "Rscript not found -- geneBody_coverage plots skipped (tabular output still produced)"
command -v samtools  >/dev/null 2>&1 || warn "samtools not found -- cannot auto-index missing .bam.bai (pipeline normally writes it)"
for _t in geneBody_coverage.py read_distribution.py gtfToGenePred genePredToBed; do
  info "${_t} -> $(command -v "$_t" 2>/dev/null || echo MISSING)"
done

# --- output dir validation ---------------------------------------------------
OUT="${OUT%/}"
[[ -d "$OUT" ]]     || die "outdir not found: $OUT (expected the flashseq_pipeline.sh run directory)"
[[ -d "$OUT/bam" ]] || die "no bam/ subdirectory in $OUT -- does this point at a flashseq_pipeline.sh run?"

# --- reference resolution ---------------------------------------------------
# References are per-TOOL: geneBody_coverage needs the HOUSEKEEPING-GENE BED;
# read_distribution needs a GENOME-WIDE BED12 baked from the GTF (it can't
# work off the housekeeping set -- rdist needs exons+introns to bucket reads).
#
#   geneBody_coverage (GBCOV_REFBED):
#     1. --refbed PATH           : ready BED, used as-is
#     2. per-species default     : REFBEDS table (housekeeping BED)
#     3. (fallback) GTF          : --gtf if given, else the REFS table GTF --
#                                  genome-wide, BAKED to BED12 (warned about
#                                  at run time)
#   read_distribution (RDIST_BED12):
#     1. --gtf PATH              : bake a BED12 from this GTF
#     2. per-species default     : REFS table GTF, BAKED to BED12
#     3. (hard-fail) neither     : err -- cannot proceed without a GTF
#
# Sets GBCOV_REFBED, GBCOV_FALLBACK_GTF, RDIST_GTF. The actual GTF->BED12
# conversion happens later (bake_bed12 called from prepare_gbcov_ref /
# prepare_rdist_ref), so dry-run and preflight checks see the intent.
_refbed_for() {   # $1=species -> housekeeping BED path or ""
  local b="" line sp p
  for line in "${REFBEDS[@]}"; do IFS='|' read -r sp p <<<"$line"; [[ "$sp" == "$1" ]] && b="$p"; done
  printf '%s' "$b"
}
_refgtf_for() {   # $1=species -> GTF path or ""
  local g="" line sp p label
  for line in "${REFS[@]}"; do IFS='|' read -r sp p label <<<"$line"; [[ "$sp" == "$1" ]] && g="$p"; done
  printf '%s' "$g"
}

# --- geneBody_coverage reference -------------------------------------------
# GBCOV_REFBED = ready housekeeping BED path; GBCOV_FALLBACK_GTF set (instead)
# only if no ready BED is available -- in that case we bake a genome-wide
# BED12 from that GTF and warn that the housekeeping default is missing.
if [[ -n "$BEDREF_ARG" ]]; then
  GBCOV_REFBED="$BEDREF_ARG"
else
  _dbed=$(_refbed_for "$SPECIES")
  if [[ -n "$_dbed" && -f "$_dbed" && -r "$_dbed" ]]; then
    GBCOV_REFBED="$_dbed"
  else
    if [[ -n "$_dbed" ]]; then
      warn "gbcov housekeeping BED missing/unreadable: $_dbed"
    fi
    # Fallback gene model: --gtf if the user gave one (custom species without
    # a housekeeping BED), else the per-species default GTF.
    _fgtf="${GTF_ARG:-$(_refgtf_for "$SPECIES")}"
    if [[ -n "$_fgtf" ]]; then
      warn "gbcov: falling back to GENOME-WIDE BED12 baked from GTF ($_fgtf) -- not the housekeeping set"
      GBCOV_FALLBACK_GTF="$_fgtf"
    else
      err "no reference for --species $SPECIES -- give --refbed PATH (ready housekeeping BED) or --gtf PATH (gene-model GTF, genome-wide BED12 baked from it)"
      echo; usage; exit 1
    fi
  fi
fi

# --- read_distribution reference -------------------------------------------
# RDIST_GTF = the GTF whose BED12 we will bake. Never the housekeeping BED.
# Priority: --gtf > per-species default.
if [[ -n "$GTF_ARG" ]]; then
  RDIST_GTF="$GTF_ARG"
else
  _sgtf=$(_refgtf_for "$SPECIES")
  if [[ -n "$_sgtf" ]]; then
    RDIST_GTF="$_sgtf"
  else
    err "read_distribution needs a genome-wide gene model GTF for --species $SPECIES"
    err "give --gtf PATH (a GTF with exon entries)"
    echo; usage; exit 1
  fi
fi

# --- output dir: done (validated above) --------------------------------------
# Create the working subdirs only for a REAL run. In --dry-run we leave the
# filesystem untouched (no logs/, no rseqc/, no files) so a preview has zero
# side effects.
if [[ "$DRY_RUN" != "1" ]]; then
  [[ -d "$OUT/logs"  ]] || { info "no logs/ in $OUT -- creating"; mkdir -p "$OUT/logs"; }
  [[ -d "$OUT/rseqc" ]] || { info "no rseqc/ in $OUT -- creating"; mkdir -p "$OUT/rseqc"; }
fi

# --- dry-run / real-run: run log --------------------------------------------
if [[ "$DRY_RUN" == "1" ]]; then
  RUN_LOG="/dev/null"
else
  RUN_LOG="$OUT/logs/rseqc_run.log"
  : > "$OUT/logs/_failed.samples"
fi
exec 3>>"$RUN_LOG"
log_tool_versions   # record versions of the tools we're about to call

# ============================================================================
# Helpers
# ============================================================================

# derive_label BAM -> LABEL
# Mirrors the pipeline's FASTQ-derived label so resume-checks line up with
# outputs flashseq_pipeline.sh --rseqc produced (which are named after the
# FASTQ-derived LBL, NOT the BAM file name).
# Pipeline BAM names: <LBL>_Aligned.sortedByCoord.out.bam (STAR prefix
# $OUT/bam/${LBL}_), plus the find_star_bam fallbacks *.out.bam / *.bam.
#   lung_001_Aligned.sortedByCoord.out.bam -> lung_001
#   SampleA_S3_R1_28_R2_90.out.bam         -> SampleA
derive_label() {
  local base="${BAM##*/}"
  local name="$base"
  name="${name%.sortedByCoord.out.bam}"
  # STAR appends a fixed "Aligned" stage before ".sortedByCoord...(out.)bam":
  # the pipeline writes --outFileNamePrefix "$OUT/bam/${LBL}_"  =>  file
  # <LBL>_Aligned.sortedByCoord.out.bam. Strip the literal _Aligned (never a
  # label, since LBL never ends in "_Aligned").
  name="${name%_Aligned}"
  name="${name%.out.bam}"
  name="${name%.bam}"
  if [[ "$name" =~ ^(.*)_S[0-9]+_ ]]; then LABEL="${BASH_REMATCH[1]}"
  elif [[ "$name" =~ ^(.*)_[Rr][12]_ ]]; then LABEL="${BASH_REMATCH[1]}"
  else LABEL="$name"; fi
}

# ensure_index BAM -- create .bam.bai if missing (samtools index).
# geneBody_coverage.py refuses an unindexed BAM. Quiet if already there.
ensure_index() {
  local bam="$1"
  [[ -f "${bam}.bai" ]] && return 0
  command -v samtools >/dev/null 2>&1 || {
    warn "no .bam.bai for $bam and samtools not in PATH -- geneBody_coverage will refuse it"
    return 1
  }
  info "$(basename "$bam"): creating .bam.bai (samtools index) ..."
  samtools index "$bam" >>"$RUN_LOG" 2>&1 \
    || warn "samtools index rc=$? on $(basename "$bam") (non-fatal; geneBody_coverage may still refuse)"
}

# rseqc_preflight -- hard-stop if any required tool is absent.
# Runs once in the parent before sample work (also during --dry-run).
# gtfToGenePred/genePredToBed/gzip are only needed if we will actually bake a
# BED12 (RDIST always bakes; GBCOV only if no ready housekeeping BED was
# found).
rseqc_preflight() {
  local miss=""
  command -v geneBody_coverage.py  >/dev/null 2>&1 || miss="$miss geneBody_coverage.py"
  command -v read_distribution.py  >/dev/null 2>&1 || miss="$miss read_distribution.py"
  [[ -n "$miss" ]] && { err "--rseqc stage missing: $miss"; conda_env_hint; exit 1; }
  # RDIST always needs to bake (we always have an GTF; GBCOV only if falling back).
  if [[ -n "$RDIST_GTF" ]]; then
    for _t in gtfToGenePred genePredToBed; do
      command -v "$_t" >/dev/null 2>&1 || { err "$_t not found -- needed to bake BED12 from GTF ($RDIST_GTF) for read_distribution"; conda_env_hint; exit 1; }
    done
    case "$RDIST_GTF" in
      *.gz)  command -v gzip  >/dev/null 2>&1 || die "gzip  not found but RDIST GTF is .gz"  ;;
      *.bz2) command -v bzcat >/dev/null 2>&1 || die "bzcat not found but RDIST GTF is .bz2" ;;
    esac
  fi
  if [[ -n "$GBCOV_FALLBACK_GTF" ]]; then
    for _t in gtfToGenePred genePredToBed; do
      command -v "$_t" >/dev/null 2>&1 || { err "$_t not found -- needed to bake BED12 from GTF ($GBCOV_FALLBACK_GTF) for geneBody_coverage (no ready housekeeping BED found)"; conda_env_hint; exit 1; }
    done
    case "$GBCOV_FALLBACK_GTF" in
      *.gz)  command -v gzip  >/dev/null 2>&1 || die "gzip  not found but GBCOV fallback GTF is .gz"  ;;
      *.bz2) command -v bzcat >/dev/null 2>&1 || die "bzcat not found but GBCOV fallback GTF is .bz2" ;;
    esac
  fi
}

# bake_bed12 GTF -- GTF -> BED12 (UCSC genePred chain). Returns the BED path
# in $BAKED_BED (also echoed). File-level guard is concurrent-safe.
# Shared by the two per-tool reference paths below; safe to call twice in the
# same run (e.g. once for RDIST, once for a GBCOV fallback) -- the second call
# will re-bake only if it's a DIFFERENT GTF (different output path), or reuse
# the existing BED12 if the target file is already there.
#
# Chain (minimal -- exactly what UCSC docs recommend, and what we verified
# deterministic on AVT0551):
#          gtfToGenePred GTF   gp   (transcript-level genePred, 8 cols)
#          genePredToBed gp    bed  (12-col BED12)
#
# NOT sorted: read_distribution reads the BED12 row-by-row and doesn't need
# bedtools ordering, so we don't add a third dependency. If you later feed
# this BED12 to a tool that does, prepend `bedSort` to the pipeline.
#
# Why UCSC chain, not gffread --bed: for transcripts whose GTF has no explicit
# CDS feature (many lncRNAs / pseudogenes), gffread fills thickStart/thickEnd
# over the WHOLE transcription-unit span -- RSeQC then classifies every exon
# as "CDS_exons". UCSC leaves thickStart=thickEnd (empty CDS) -- which is the
# biologically faithful reading of "no CDS in this transcript". Verified on
# AVT0551 lung_116: identical read totals, different CDS/UTR partition
# (gffread CDS 146M bases, UCSC CDS 38M -- the difference moved to 5'/3' UTR).
bake_bed12() {
  local gtf="$1"
  local bed="$OUT/rseqc/$(basename "${gtf%.gz}").bed12"
  BAKED_BED="$bed"
  # Track this genome-wide BED12 for post-run cleanup (dedup). It is a BAKED
  # intermediate in every case that reaches here (both rdist and the gbcov
  # fallback call this function) — never a ready housekeeping BED.
  if [[ " ${BAKED_BEDS[*]:-} " != *" $bed "* ]]; then BAKED_BEDS+=("$bed"); fi
  if [[ -s "$bed" ]]; then
    info "baking BED12: using existing (skip): $bed"
    return 0
  fi
  if [[ "$DRY_RUN" == "1" ]]; then
    # De-dup: MAIN calls prepare_rdist_ref first, then run_combined_gbcov
    # calls prepare_gbcov_ref; if both resolve to the same GTF (and thus the
    # same output BED) we only want this logged once.
    if [[ -z "${_BAKED_DRYRUN_ANNOUNCED[$bed]:-}" ]]; then
      if [[ "$gtf" == *.gz ]]; then
        info "  [dry-run] gzip -dc   $gtf | gtfToGenePred -  tmp.gp"
        info "             genePredToBed tmp.gp  $bed     (UCSC chain)"
      else
        info "  [dry-run] gtfToGenePred $gtf tmp.gp ; genePredToBed tmp.gp $bed   (UCSC chain)"
      fi
      _BAKED_DRYRUN_ANNOUNCED["$bed"]=1
    else
      info "  [dry-run] baking (reused from same GTF): $bed"
    fi
    return 0
  fi
  local tmpdir rc=1
  tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/rseqc_bed12.XXXXXX")
  local gp="$tmpdir/genes.gp" bedtmp="$tmpdir/genes.bed"
  # -- GTF -> genePred -------------------------------------------------
  if   [[ "$gtf" == *.gz  ]]; then gzip -dc  "$gtf"  | gtfToGenePred -  "$gp"; rc=$?
  elif [[ "$gtf" == *.bz2 ]]; then bzcat -d  "$gtf"  | gtfToGenePred -  "$gp"; rc=$?
  else                             gtfToGenePred     "$gtf"           "$gp"; rc=$?; fi
  (( rc == 0 )) || { rm -rf "$tmpdir"; die "gtfToGenePred failed rc=$rc on $gtf -- check GTF validity / ucsc-gtftogenepred version"; }
  [[ -s "$gp" ]]  || { rm -rf "$tmpdir"; die "gtfToGenePred produced EMPTY genePred for $gtf -- is it a valid GTF with exon entries?"; }
  # -- genePred -> BED12 (genePredToBed takes <in.gp> <out.bed>; it prints
  # usage to stderr and returns 0 in the usage-error case, so we MUST check
  # the output file size, not just the exit code) ----------------------
  if ! genePredToBed "$gp" "$bedtmp"; then
    rm -rf "$tmpdir"
    die "genePredToBed failed on $gtf (gp=$(du -h "$gp" 2>/dev/null | cut -f1)) -- is the genePred file consistent?";
  fi
  [[ -s "$bedtmp" ]] || { rm -rf "$tmpdir"; die "genePredToBed produced EMPTY/missing BED for $gtf (gp=$(du -h "$gp" 2>/dev/null | cut -f1) gp) -- genePred may be corrupt or the tool rejected it"; }
  mv -f "$bedtmp" "$bed"
  rm -rf "$tmpdir"
  info "$bed: $(wc -l <"$bed") BED12 rows (UCSC genePred chain: gtfToGenePred | genePredToBed)"
}

# prepare_rdist_ref -- resolve RDIST_BED12 (always: bake BED12 from RDIST_GTF).
# Idempotent (checks existing output).
prepare_rdist_ref() {
  [[ -n "$RDIST_BED12" && -s "$RDIST_BED12" ]] && return 0
  [[ -n "$RDIST_GTF" ]] || die "RDIST_GTF not set -- call prepare_rdist_ref after reference resolution"
  bake_bed12 "$RDIST_GTF"
  RDIST_BED12="$BAKED_BED"
}

# prepare_gbcov_ref -- resolve GBCOV_BED12.
#   If GBCOV_REFBED (ready housekeeping BED) is set and looks like a proper
#   BED12 -> use as-is.
#   If the housekeeping BED is missing/unreadable, fall back to baking a
#   genome-wide BED12 from GBCOV_FALLBACK_GTF (a warning is emitted).
# Idempotent.
prepare_gbcov_ref() {
  [[ -n "$GBCOV_BED12" && -s "$GBCOV_BED12" ]] && return 0
  if [[ -n "$GBCOV_REFBED" ]]; then
    GBCOV_BED12="$GBCOV_REFBED"
    local nc
    nc=$(awk -F'\t' 'NR==1{print NF; exit}' "$GBCOV_REFBED" 2>/dev/null || echo "?")
    info "gbcov reference: ready BED (as-is): $GBCOV_BED12"
    if   [[ "$nc" == "4" || "$nc" == "6" ]]; then
      warn "gbcov reference $GBCOV_REFBED looks like a region list ($nc cols) -- geneBody_coverage.py will likely skip all rows"
    elif [[ "$nc" != "12" && "$nc" != "?" ]]; then
      warn "gbcov reference $GBCOV_REFBED has $nc columns -- RSeQC expects 12 (BED12); non-12-column rows may be skipped"
    fi
    return 0
  fi
  [[ -n "$GBCOV_FALLBACK_GTF" ]] || die "GBCOV_REFBED and GBCOV_FALLBACK_GTF both empty -- call prepare_gbcov_ref after reference resolution"
  if [[ -z "${_BAKED_DRYRUN_ANNOUNCED[gbcov_fallback_$(basename "$GBCOV_FALLBACK_GTF")]:-}" ]]; then
    info "gbcov: no housekeeping BED found -- using GENOME-WIDE BED12 (baked from $GBCOV_FALLBACK_GTF)"
    _BAKED_DRYRUN_ANNOUNCED["gbcov_fallback_$(basename "$GBCOV_FALLBACK_GTF")"]=1
  fi
  bake_bed12 "$GBCOV_FALLBACK_GTF"
  GBCOV_BED12="$BAKED_BED"
}

# gbcov_ref_sanity BAM BED -- fail FAST (before any sample is processed) on the
# single root cause of "Cannot get coverage signal from <bam>! Skip" on EVERY
# sample: the reference gene-model BED shares NO contig name with the BAM, or
# has no parseable BED12 rows. Contig-name mismatches (M vs chrM vs NC_...,
# or an assembly STAR used that the BED's source assembly doesn't name
# identically) are invisible in the per-sample error -- RSeQC skips a contig
# it can't pileup -- so without this check the failure surfaces late and
# cryptically as "no valid coverage profiles were generated".
# GBCOV_SKIP_REF_CHECK=1 disables (e.g. contigs deliberately renamed by the caller).
gbcov_ref_sanity() {
  local bam="$1" bed="$2"
  [[ "${GBCOV_SKIP_REF_CHECK:-0}" == "1" ]] && return 0
  [[ -s "$bed" ]] || die "gbcov check: reference BED is empty: $bed"
  command -v samtools >/dev/null 2>&1 || { warn "gbcov check: samtools not found -- skipping contig-overlap sanity check"; return 0; }
  local bam_contigs bed_contigs n_common n_bed n_bam
  bam_contigs=$(samtools view -H "$bam" 2>/dev/null | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i; break}}' | sort -u)
  bed_contigs=$(awk '!/^#/{print $1}' "$bed" 2>/dev/null | sort -u)
  n_bam=$(printf '%s' "$bam_contigs" | grep -c . || true)
  n_bed=$(printf '%s' "$bed_contigs" | grep -c . || true)
  n_common=$(comm -12 <(printf '%s\n' "$bam_contigs") <(printf '%s\n' "$bed_contigs") 2>/dev/null | grep -c . || true)
  if (( n_bed == 0 )); then
    die "gbcov reference has no parseable rows (BED must have data rows; is it a BED12?) -- $bed"
  fi
  if (( n_common == 0 )); then
    die "gbcov reference $bed shares NO contig name with $bam.
  BAM contigs  : $(printf '%s ' $bam_contigs | cut -d' ' -f1-6) (of $n_bam)
  REF contigs  : $(printf '%s ' $bed_contigs | cut -d' ' -f1-6) (of $n_bed)
  This is the classic cause of 'Cannot get coverage signal ... Skip' on every
  sample. Fix: pass a --refbed whose contig names MATCH the contigs your BAMs
  were aligned to (build the BED12 from the exact GTF/assembly STAR used).
  Contig naming differs from assembly to assembly (M vs chrM vs NC_..., GRCh38
  vs T2T vs a draft). Override with GBCOV_SKIP_REF_CHECK=1 only if you know
  the reference is intentionally renamed."
  fi
  info "gbcov ref/BAM contig overlap: $n_common / $n_bed BED contigs present in BAM (OK)"
}

# stage_rseqc SKIP_GBCOV -- run the two RSeQC tools for one BAM (already resolved).
# SKIP_GBCOV=1 -> skip the geneBody_coverage block inside (used in MERGED mode,
# where the combined-BAM gbcov is handled separately by run_combined_gbcov).
# Uses LABEL + BAM set by the caller. Soft-fail per design: QC tool failures
# don't abort the sample; they're flagged and written to _failed.samples so
# the final summary can show what didn't produce output.
stage_rseqc() {
  local skip_gbcov="${1:-0}"
  local log="$OUT/logs/${LABEL}.rseqc.log"
  info "$LABEL: RSeQC starting ($STEPS_LABEL$([[ $skip_gbcov -eq 1 ]] && echo ", skip-gbcov"))"
  # Resolve per-tool references (idempotent -- safe to call per-sample).
  [[ "$RUN_GBCOV"  -eq 1 ]] && prepare_gbcov_ref
  [[ "$RUN_RDIST" -eq 1 ]] && prepare_rdist_ref
  # In --per-sample mode stage_rseqc (skip_gbcov=0) is the ONLY path into
  # geneBody_coverage.py, so this is where the pre-sample contig-overlap
  # sanity check lives (in MERGED mode run_combined_gbcov already does it
  # once for the whole batch, and stage_rseqc is called with skip_gbcov=1
  # so we skip double-checking here).
  if [[ "$skip_gbcov" -eq 0 && "$RUN_GBCOV" -eq 1 && "$DRY_RUN" != "1" ]]; then
    gbcov_ref_sanity "$BAM" "$GBCOV_BED12" || return 1
  fi
  ensure_index "$BAM"
  local plot_flag=()
  if [[ "$skip_gbcov" -eq 0 ]] && [ "$RUN_GBCOV" -eq 1 ]; then
    command -v Rscript >/dev/null 2>&1 || plot_flag=(--skip-plot)
  fi
  {
    date -Is
    echo "BAM                  : $BAM"
    echo "gbcov  reference     : $GBCOV_BED12"
    echo "rdist reference (BED12) : $RDIST_BED12"
    echo "steps                : $STEPS_LABEL"
    echo
    if [[ $skip_gbcov -eq 0 ]] && [ "$RUN_GBCOV" -eq 1 ]; then
      if [[ "${#plot_flag[@]}" -eq 0 ]]; then
        echo "+ geneBody_coverage.py -i $BAM -r $GBCOV_BED12 -o $OUT/rseqc/${LABEL}_gbcov"
      else
        echo "+ geneBody_coverage.py -i $BAM -r $GBCOV_BED12 -o $OUT/rseqc/${LABEL}_gbcov --skip-plot"
        echo "  (Rscript not in PATH; plots skipped; tabular output still written)"
      fi
      geneBody_coverage.py -i "$BAM" -r "$GBCOV_BED12" \
        -o "$OUT/rseqc/${LABEL}_gbcov" ${plot_flag[@]+"${plot_flag[@]}"}
      echo
    fi
    if [ "$RUN_RDIST" -eq 1 ]; then
      echo "+ read_distribution.py -i $BAM -r $RDIST_BED12"
      echo "  (stdout below -- 5.0.5 has no -o flag; we extract it to ${LABEL}.readDistribution.txt)"
      read_distribution.py -i "$BAM" -r "$RDIST_BED12"
    fi
  } >"$log" 2>&1
  local rc=$?
  # read_distribution.py 5.0.5 prints a 3-line summary + === separator + Group
  # table; awk stops at the SECOND === row. A range ending at the first ===
  # row silently drops the whole table (known bug, fixed in pipeline 2026-09-15).
  if [ "$RUN_RDIST" -eq 1 ]; then
    local rd="$OUT/rseqc/${LABEL}.readDistribution.txt"
    awk '/Total Reads/{f=1} f{print} f && /^=*$/ {n++; if(n==2)exit}' \
      "$log" >"$rd" 2>/dev/null || echo "(no report -- see $log)" >"$rd"
  fi
  # RSeQC 5.0.5 also writes a throwaway generated R script per sample
  # (<LBL>_gbcov.geneBodyCoverage.r). Its plots (.curves.pdf) are already
  # baked into the data before this point, and we don't want the .r kept, so
  # drop it right away (resume is gated on the .txt, not the .r). Only
  # relevant when geneBody_coverage actually ran for THIS sample.
  if [[ $skip_gbcov -eq 0 ]] && [ "$RUN_GBCOV" -eq 1 ]; then
    rm -f "$OUT/rseqc/${LABEL}_gbcov.geneBodyCoverage.r"
  fi
  if (( rc != 0 )); then
    warn "$LABEL: RSeQC rc=$rc (soft-fail per design -- sample skipped, see $log)"
    echo "$LABEL" >>"$OUT/logs/_failed.samples"
    return 1
  fi
  local done_msg="RSeQC done"
  if [[ $skip_gbcov -eq 0 ]] && [ "$RUN_GBCOV" -eq 1 ]; then
    done_msg="$done_msg -> $OUT/rseqc/${LABEL}_gbcov.geneBodyCoverage.txt"
  fi
  if [[ $skip_gbcov -eq 0 ]] && [ "$RUN_GBCOV" -eq 1 ] && [ "$RUN_RDIST" -eq 1 ]; then
    done_msg="$done_msg + "
  fi
  if [ "$RUN_RDIST" -eq 1 ]; then
    done_msg="$done_msg $OUT/rseqc/${LABEL}.readDistribution.txt"
  fi
  info "$LABEL: $done_msg"
  return 0
}

# is_resumable SAMPLE_LABEL -- true if this sample's outputs already exist,
# considering only the tool(s) selected via --steps. A sample is "done" when
# every selected tool already has its non-empty primary output.
# (In MERGED mode the per-sample gbcov output is produced by the COMBINED
# call, not per sample -- use is_resumable_rdist for the per-sample pass then.)
is_resumable() {
  local lbl="$1"
  if [ "$RUN_GBCOV" -eq 1 ]; then
    [[ -s "$OUT/rseqc/${lbl}_gbcov.geneBodyCoverage.txt" ]] || return 1
  fi
  if [ "$RUN_RDIST" -eq 1 ]; then
    [[ -s "$OUT/rseqc/${lbl}.readDistribution.txt" ]] || return 1
  fi
  return 0
}

# is_resumable_rdist LABEL -- true if the sample's read_distribution output
# already exists (merged mode: per-sample pass is rdist-only). Trivially true
# when rdist is not selected (nothing per-sample to do).
is_resumable_rdist() {
  local lbl="$1"
  if [ "$RUN_RDIST" -eq 1 ]; then
    [[ -s "$OUT/rseqc/${lbl}.readDistribution.txt" ]] || return 1
  fi
  return 0
}

# merged_gbcov_resumable -- true if the combined gbcov output exists AND was
# built from exactly the same set of sample BAMs as currently selected
# (members list match; a --samples subset change invalidates the pool).
merged_gbcov_resumable() {
  local gtxt="$OUT/rseqc/geneBodyCombined.geneBodyCoverage.txt"
  local members="$OUT/rseqc/geneBodyCombined.members.list"
  local want cur
  [[ -s "$gtxt" ]] || return 1
  want="$(for lbl in "${!SAMPLES_LIST[@]}"; do printf '%s\n' "${SAMPLES_LIST[$lbl]}"; done | sort)"
  cur="$( [[ -f "$members" ]] && sort "$members" || true)"
  [[ -n "$want" && "$want" == "$cur" ]]
}

# gbcov_row_ok FILE -- true if a per-sample geneBodyCoverage.txt actually
# holds >=1 sample data row (header-only files can be left behind by a
# failed geneBody_coverage.py call; the header is written before the loop).
gbcov_row_ok() {
  local f="$1"
  [[ -s "$f" ]] || return 1
  tail -n +2 "$f" | grep -q . || return 1
  return 0
}

# gbcov_parallel_one -- one sample of the PARALLEL merged gbcov phase.
# Runs geneBody_coverage.py --skip-plot against the SINGLE BAM and writes
# <per_sample dir>/<LBL>.geneBodyCoverage.txt (header + one row) -- the
# table only; no per-sample plotting (the combined plot is redrawn ONCE
# after concatenation, which is exactly what --per-sample uses -- except
# that it is rendered ONCE instead of once per sample).
# LABEL/BAM/GBCOV_BED12 are set by the caller; safe as a subshell worker.
gbcov_parallel_one() {
  local psdir="$OUT/rseqc/geneBodyCombined.per_sample"
  local prefix="$psdir/${LABEL}"
  local log="$OUT/logs/gbcov_parallel_${LABEL}.rseqc.log"
  {
    date -Is
    echo "+ geneBody_coverage.py -i $BAM -r $GBCOV_BED12 -o $prefix --skip-plot"
    geneBody_coverage.py -i "$BAM" -r "$GBCOV_BED12" -o "$prefix" --skip-plot
  } >"$log" 2>&1
  local rc=$?
  # --skip-plot still emits a throwaway per-sample .r we never use; drop it.
  rm -f "$prefix.geneBodyCoverage.r"
  (( rc != 0 )) && return 1
  gbcov_row_ok "$prefix.geneBodyCoverage.txt" || return 1
  return 0
}

# run_combined_gbcov -- MERGED mode (default): produce ONE combined gene
# body output (coverage table + heatMap + curves) over the selected samples.
#
#   v1.1.0 (only implementation):
#     Phase 1: geneBody_coverage.py --skip-plot PER sample, fanned out to
#              --max-jobs (the actual BAM reads happen here).
#     Phase 2: concatenate the per-sample tables (one row per sample, in
#              the sorted member order RSeQC would use), with one shared
#              header -- a flat (all-zero) sample simply has no row and is
#              flagged, instead of aborting the whole run.
#     Phase 3: redraw the combined heatMap+curves from the concatenated
#              table with gbcov_redraw.py (reuses RSeQC's own R template --
#              no BAMs are re-read).
#
# Artifacts (shared, under $OUT/):
#   rseqc/geneBodyCombined.geneBodyCoverage.txt     (coverage table, multi-row)
#   rseqc/geneBodyCombined.geneBodyCoverage.heatMap.pdf + .curves.pdf (plots)
#   rseqc/geneBodyCombined.geneBodyCoverage.r       (throwaway, removed)
#   rseqc/geneBodyCombined.members.list             (sorted input BAMs, resume)
#   rseqc/geneBodyCombined.per_sample/<LBL>.geneBodyCoverage.txt (Phase 1)
#   logs/gbcov_parallel_<LBL>.rseqc.log            (Phase 1 per-sample logs)
#
# Precondition: every selected BAM has a .bai in the same directory
# (RSeQC's getBamFiles.isbamfile hard-requires this). We index in
# parallel before the RSeQC phase so per-BAM index work isn't serial.
#
# Per-sample geneBody outputs are NOT produced in this mode (use --per-sample).
# The rdist pass still runs per-sample independently (see run_in_batch).
run_combined_gbcov() {
  local outprefix="$OUT/rseqc/geneBodyCombined"
  local gtxt="${outprefix}.geneBodyCoverage.txt"
  local members="$OUT/rseqc/geneBodyCombined.members.list"
  local want_members
  want_members="$(for lbl in "${!SAMPLES_LIST[@]}"; do printf '%s\n' "${SAMPLES_LIST[$lbl]}"; done | sort)"
  local n=0
  for lbl in "${!SAMPLES_LIST[@]}"; do n=$(( n + 1 )); done
  [[ "$n" -gt 0 ]] || die "run_combined_gbcov: no samples selected"

  # Resolve the GBCOV reference (idempotent). Must happen BEFORE the dry-run
  # branch since that branch prints the resolved path to the command line.
  prepare_gbcov_ref

  # One-time, pre-sample reference/BAM sanity check (see gbcov_ref_sanity): a
  # contig-name mismatch in the reference is the classic cause of
  # "Cannot get coverage signal ... Skip" on EVERY sample, and it is invisible
  # in the per-sample error. Fail here with an actionable message instead of
  # dying later with "no valid coverage profiles were generated".
  if [[ "$DRY_RUN" != "1" ]]; then
    local _sany_key=() _sany_bam=""
    for _k in "${!SAMPLES_LIST[@]}"; do _sany_key=("$_k"); break; done
    if [[ -n "${_sany_key[0]:-}" ]]; then
      _sany_bam="${SAMPLES_LIST[${_sany_key[0]}]}"
      gbcov_ref_sanity "$_sany_bam" "$GBCOV_BED12" || return 1
    fi
  fi

  if merged_gbcov_resumable; then
    info "combined gbcov already done for these $n samples -- reusing: $gtxt"
    return 0
  fi

  # Parallel index pass: RSeQC's isbamfile() skips any BAM without .bai.
  # stage_rseqc in --per-sample mode also indexes; here we must do it too so
  # combined mode works standalone. Parallelism = MAX_JOBS.
  local idx_targets=()
  for lbl in "${!SAMPLES_LIST[@]}"; do
    local b="${SAMPLES_LIST[$lbl]}"
    [[ -f "${b}.bai" ]] || idx_targets+=("$b")
  done
  if [[ ${#idx_targets[@]} -gt 0 ]]; then
    if [[ "$DRY_RUN" == "1" ]]; then
      info "  [dry-run] samtools index (parallel, max $MAX_JOBS) for ${#idx_targets[@]} sample BAM(s) without .bai"
    else
      command -v samtools >/dev/null 2>&1 \
        || die "samtools not found -- required to create .bai for the ${#idx_targets[@]} missing index(es)"
      info "indexing ${#idx_targets[@]} sample BAM(s) (parallel, max $MAX_JOBS) ..."
      local ipids=()
      for b in "${idx_targets[@]}"; do
        ( samtools index "$b" >>"$RUN_LOG" 2>&1 \
            || warn "samtools index rc=$? on $(basename "$b")" ) &
        ipids+=("$!")
        if (( ${#ipids[@]} >= MAX_JOBS )); then
          wait "${ipids[0]}" 2>/dev/null || true
          ipids=("${ipids[@]:1}")
        fi
      done
      local p; for p in "${ipids[@]}"; do wait "$p" 2>/dev/null || true; done
    fi
  fi

    # ------------------------------------------------------------------
    # PARALLEL (only implementation):
    #   1. geneBody_coverage.py --skip-plot per sample (fanned to MAX_JOBS)
    #   2. concatenate per-sample tables into the combined table
    #   3. redraw the combined heatMap+curves via gbcov_redraw.py (table
    #      in, R out -- no BAM re-read)
    # ------------------------------------------------------------------
    local psdir="$OUT/rseqc/geneBodyCombined.per_sample"
    local plog="$OUT/logs/gbcovParallel.rseqc.log"

    if [[ "$DRY_RUN" == "1" ]]; then
      info "  [dry-run] PARALLEL merged gbcov ($n sample(s), $MAX_JOBS workers):"
      info "            per sample : geneBody_coverage.py -i <BAM> -r $GBCOV_BED12 -l 100 -o $psdir/<LBL> --skip-plot \  (logs/gbcov_parallel_<LBL>.rseqc.log)"
      info "            concatenate: per-sample tables -> $gtxt"
      info "            redraw     : python3 $GBCOV_REDRAW_PY --table $gtxt --prefix $outprefix --format $GBCOV_FORMAT"
      return 0
    fi

    [[ -f "$GBCOV_REDRAW_PY" ]] \
      || die "gbcov_redraw.py not found at $GBCOV_REDRAW_PY (put it in the same folder as flashseq_rseqc.sh)"
    command -v python3 >/dev/null 2>&1 \
      || die "python3 not found -- required by gbcov_redraw.py"
    if [[ "$GBCOV_FORMAT" == "html" ]]; then
      # html is rendered in Python via plotly (no Rscript needed). Fail fast
      # if plotly is missing so we don't run the whole parallel batch only to
      # die at the redraw step.
      python3 -c "import plotly" 2>/dev/null \
        || die "plotly not importable by python3 -- required for --gbcov-format html (pip install plotly)"
    else
      command -v Rscript >/dev/null 2>&1 \
        || die "Rscript not found -- required to redraw the combined plot (gbcov_redraw.py, format=$GBCOV_FORMAT)"
    fi
    mkdir -p "$psdir"

    {
      date -Is
      echo "MODE    : PARALLEL merged (v1.1.0) -- $n sample(s), $MAX_JOBS workers"
      echo "PREFIX  : $outprefix"
      echo "MEMBERS :"
      printf '%s\n' "$want_members" | sed 's/^/  /'
      echo "gbcov ref : $GBCOV_BED12"
      echo "redraw  : $GBCOV_REDRAW_PY"
    } >"$plog" 2>&1

    # Phase 1: per-sample --skip-plot table calls, fanned out to MAX_JOBS.
    # Per-sample resume: a row file that already holds data skips the worker.
    local pids=() lbl_p
    for lbl_p in $(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort); do
      LABEL="$lbl_p"
      BAM="${SAMPLES_LIST[$lbl_p]}"
      if gbcov_row_ok "$psdir/$lbl_p.geneBodyCoverage.txt"; then
        info "  [skip] $lbl_p (gbcov table row already present)"
        continue
      fi
      gbcov_parallel_one &
      pids+=("$!")
      if (( ${#pids[@]} >= MAX_JOBS )); then
        wait "${pids[0]}" 2>/dev/null || true
        pids=("${pids[@]:1}")
      fi
    done
    local p
    for p in "${pids[@]}"; do wait "$p" 2>/dev/null || true; done

    # Collect the outcome per sample. A missing row (flat profile or tool
    # failure) is FLAGGED, not fatal -- the combined table is still valid
    # for the samples that do have rows.
    local have=() missing=()
    for lbl_p in "${!SAMPLES_LIST[@]}"; do
      if gbcov_row_ok "$psdir/$lbl_p.geneBodyCoverage.txt"; then have+=("$lbl_p")
      else missing+=("$lbl_p"); fi
    done
    if [[ ${#have[@]} -eq 0 ]]; then
      warn "parallel gbcov: NO sample produced a coverage row -- see logs/gbcov_parallel_<LBL>.rseqc.log"
      echo "combined_gbcov" >>"$OUT/logs/_failed.samples"
      return 1
    fi
    if [[ ${#missing[@]} -gt 0 ]]; then
      warn "parallel gbcov: ${#missing[@]} sample(s) have no row (flat profile or failure): $(printf '%s ' "${missing[@]}")"
      for lbl_p in "${missing[@]}"; do echo "gbcov_flat:$lbl_p" >>"$OUT/logs/_failed.samples"; done
    fi

    # Phase 2: concatenate the per-sample rows into the combined table.
    # Header once + rows in sorted member order == what a single
    # geneBody_coverage.py -i <all BAMs> call would write (same row
    # naming/order), but flat samples are simply omitted and flagged
    # rather than aborting the run.
    local tmp_gtxt concat_rc=0
    tmp_gtxt=$(mktemp "${TMPDIR:-/tmp}/gbcovCombined.XXXXXX")
    {
      head -n 1 "$psdir/${have[0]}.geneBodyCoverage.txt"
      for lbl_p in "${have[@]}"; do tail -n +2 "$psdir/$lbl_p.geneBodyCoverage.txt"; done
    } >"$tmp_gtxt" || concat_rc=$?
    if (( concat_rc == 0 )) && [[ -s "$tmp_gtxt" ]]; then
      mv -f "$tmp_gtxt" "$gtxt"
      info "parallel gbcov: concatenated ${#have[@]}/${n} sample tables -> $gtxt"
    else
      rm -f "$tmp_gtxt"
      warn "parallel gbcov: table concatenation failed (rc=$concat_rc)"
      echo "combined_gbcov" >>"$OUT/logs/_failed.samples"
      return 1
    fi

    # Phase 3: redraw the combined heatMap + curves from the concatenated
    # table (gbcov_redraw.py reuses RSeQC's own R plotting template -- the
    # BAMs are never re-read during plotting).
    info "parallel gbcov: redrawing combined plot (gbcov_redraw.py) ..."
    local drrc=0
    python3 "$GBCOV_REDRAW_PY" --table "$gtxt" --prefix "$outprefix" --format "$GBCOV_FORMAT" \
      >>"$plog" 2>&1 || drrc=$?
    if (( drrc != 0 )); then
      warn "gbcov_redraw.py rc=$drrc failed (see $plog) -- table still produced: $gtxt"
      echo "combined_gbcov_redraw" >>"$OUT/logs/_failed.samples"
      return 1
    fi
    rm -f "${outprefix}.geneBodyCoverage.r"
    # Members list only after success so a failure never poisons the resume
    # check (same invariant as the --per-sample path).
    printf '%s\n' "$want_members" | sort >"$members"
    info "parallel gbcov done -> $gtxt (+ $outprefix.heatMap.$GBCOV_FORMAT / $outprefix.curves.$GBCOV_FORMAT)"
    return 0
}

# discover_samples -- list BAMs under $OUT/bam/ (sortedByCoord preferred).
# Sets SAMPLES_LIST[label]=bam_path.
discover_samples() {
  local bam
  for bam in "$OUT"/bam/*.sortedByCoord.out.bam "$OUT"/bam/*.out.bam "$OUT"/bam/*.bam; do
    [[ -f "$bam" && -s "$bam" ]] || continue
    # Skip duplicates: prefer sortedByCoord; once we've seen it, don't add
    # the same logical sample from a different glob again.
    BAM="$bam" derive_label
    if [[ -z "${SAMPLES_LIST[$LABEL]:-}" ]]; then
      SAMPLES_LIST["${LABEL}"]="$bam"
    fi
  done
}

# ===========================================================================
# MAIN
# ===========================================================================
info "=== Flash-seq RSeQC QC (standalone) ==="
info "species    : $SPECIES"
[ "$RUN_GBCOV" -eq 1 ] && {
  if [[ -n "$GBCOV_REFBED" ]]; then
    info "gbcov  ref : $GBCOV_REFBED  (housekeeping-gene BED, used as-is)"
  else
    info "gbcov  ref : genome-wide BED12 baked from $GBCOV_FALLBACK_GTF  (housekeeping BED unavailable)"
  fi
}
[ "$RUN_RDIST" -eq 1 ] && info "rdist  ref : genome-wide BED12 baked from $RDIST_GTF"
info "out dir    : $OUT"
info "steps      : $STEPS_LABEL"
info "parallel   : ${MAX_JOBS} samples"
info "resume     : on (skip samples with existing rseqc outputs)  dry-run=$DRY_RUN"

discover_samples
[[ ${#SAMPLES_LIST[@]} -gt 0 ]] || die "no STAR BAMs found under $OUT/bam/ -- did the pipeline run?"
TOTAL=${#SAMPLES_LIST[@]}

# --samples filter -- warn (not die) on unmatched labels so a partial typo
# still processes the valid subset instead of halting the whole run.
if [[ -n "$SAMPLES" ]]; then
  # Capture the full discovered set NOW, before the filter unsets any keys.
  known_all="$(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort | tr '\n' ' ')"
  IFS=',' read -r -a FILTER <<<"$SAMPLES"
  for lbl in "${!SAMPLES_LIST[@]}"; do
    keep=0
    for f in "${FILTER[@]}"; do [[ "$lbl" == "$f" ]] && { keep=1; break; }; done
    (( keep )) || unset "SAMPLES_LIST[$lbl]"
  done
  unmatched=()
  for f in "${FILTER[@]}"; do
    [[ -n "${SAMPLES_LIST[$f]:-}" ]] || unmatched+=("$f")
  done
  if (( ${#unmatched[@]} > 0 )); then
    warn "--samples: label(s) NOT MATCHED under $OUT/bam/: ${unmatched[*]}"
    warn "--samples: available labels: $known_all"
  fi
  [[ ${#SAMPLES_LIST[@]} -gt 0 ]] || {
    die "--samples matched nothing under $OUT/bam/. Available labels: ${known_all:-<none>}"
  }
fi
info "samples    : ${#SAMPLES_LIST[@]} of $TOTAL discovered"
if [ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 0 ] && [ "$GBCOV_FORMAT" != "pdf" ]; then
  warn "--gbcov-format $GBCOV_FORMAT only affects the MERGED combined plot; --per-sample plots stay PDF (RSeQC's own output)"
fi

# ============================================================================
# Reference / sample pre-flight check (runs in BOTH real and dry-run modes).
# Consolidates the per-file existence/readability/non-empty checks into one
# clear [PASS]/[FAIL] list. Dies with a single actionable message at the end
# if anything is missing, so we don't walk forward into the UCSC chain /
# RSeQC only to discover the input was 0-byte or on a different machine.
# ============================================================================
check_ref_file() {
  # $1 = human label, $2 = path. Reports [ok]/[FAIL]; appends to _REF_FAILS.
  local label="$1" p="$2" reason
  if [[ ! -e "$p" ]]; then
    reason="MISSING"
  elif [[ ! -r "$p" ]]; then
    reason="UNREADABLE (permissions)"
  elif [[ ! -s "$p" ]]; then
    reason="EMPTY (0 bytes)"
  else
    local sz; sz=$(stat -c%s -- "$p" 2>/dev/null || stat -f%z -- "$p" 2>/dev/null || echo "?")
    info "[ok]   $label : $p  ($sz B)"
    return 0
  fi
  local l; l="[FAIL] $label : $p -- $reason"
  echo "$l"; _log_to_file "$l"
  _REF_FAILS+=("$label ($reason)")
}
_REF_FAILS=()
if [ "$RUN_GBCOV" -eq 1 ]; then
  if [[ -n "$GBCOV_REFBED" ]]; then
    check_ref_file "gbcov housekeeping BED" "$GBCOV_REFBED"
  elif [[ -n "$GBCOV_FALLBACK_GTF" ]]; then
    check_ref_file "gbcov genome-wide GTF (fallback)" "$GBCOV_FALLBACK_GTF"
  else
    _REF_FAILS+=("gbcov reference (neither GBCOV_REFBED nor GBCOV_FALLBACK_GTF set)")
  fi
fi
if [ "$RUN_RDIST" -eq 1 ]; then
  check_ref_file "rdist genome-wide GTF" "$RDIST_GTF"
else
  [[ -z "$RDIST_GTF" ]] || warn "rdist not in --steps, but RDIST_GTF=$RDIST_GTF (ignored)"
fi
if (( ${#_REF_FAILS[@]} > 0 )); then
  echo
  err "reference pre-flight FAILED -- missing/unusable input(s):"
  for r in "${_REF_FAILS[@]}"; do err "  - $r"; done
  err "fix these, then re-run. (Use --refbed / --gtf to supply alternative paths.)"
  exit 1
fi
info "[PASS] references: all required input files present and non-empty"

# (Legacy per-tool [[ -f ]] hard-fails removed: check_ref_file below/above
# covers missing/unreadable/empty in ONE pass and runs in BOTH dry-run and
# real modes; a dedicated post-flight re-check is therefore unnecessary.)

rseqc_preflight
[ "$RUN_GBCOV" -eq 1 ] && prepare_gbcov_ref
[ "$RUN_RDIST" -eq 1 ] && prepare_rdist_ref

# NOTE: discover_samples + `--samples` filter + reference pre-flight check run
# ONCE, earlier in the flow (before prepare_*_ref); not repeated here.

# --- DRY-RUN preview (merged mode first, then per-sample rdist if needed) ---
if [[ "$DRY_RUN" == "1" ]]; then
  if [ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 1 ]; then
    info "MODE     : MERGED -- one combined coverage table + two combined plots for ${#SAMPLES_LIST[@]} sample(s), built in parallel"
    # prepare_gbcov_ref resolves GBCOV_BED12 already (above).
    run_combined_gbcov   # prints the intended commands and returns 0
  fi
  if [ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 0 ]; then
    info "MODE     : PER-SAMPLE -- one geneBody_coverage.py call + per-sample plots, per sample"
  fi
  # Per-sample pass (rdist and/or per-sample gbcov in per-sample mode).
  local_skip_gbcov_for_presample=0
  [ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 1 ] && local_skip_gbcov_for_presample=1
  for lbl in $(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort); do
    BAM="${SAMPLES_LIST[$lbl]}"
    if [ "$local_skip_gbcov_for_presample" -eq 1 ]; then
      # merged mode: per-sample pass is rdist-only
      if is_resumable_rdist "$lbl"; then
        info "  [skip] $lbl (already has rdist output)"
      else
        _bn="$(basename "${SAMPLES_LIST[$lbl]}")"
        info "  [run]  $lbl -> $_bn  (rdist only)"
      fi
    else
      if is_resumable "$lbl"; then
        info "  [skip] $lbl (already has rseqc outputs)"
      else
        info "  [run]  $lbl -> $(basename "${SAMPLES_LIST[$lbl]}")"
      fi
    fi
  done
  echo
  echo "=========================== DRY RUN complete -- nothing was launched =========================== "
  exit 0
fi

START=$SECONDS
# Merged gbcov (default) -- one combined geneBody_coverage.py call, BEFORE the
# per-sample rdist pass. rdist is always per-sample (it reports read
# distribution of each individual BAM), so this split is clean:
#   MERGED=1:     [combined gbcov]  then [per-sample rdist]
#   MERGED=0:    [per-sample gbcov + per-sample rdist]
if [ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 1 ]; then
  run_combined_gbcov || warn "combined gbcov reported failure (soft-fail, continuing)"
fi

SKIP_GBCOV=0
[ "$RUN_GBCOV" -eq 1 ] && [ "$MERGED" -eq 1 ] && SKIP_GBCOV=1   # already done in the pool pass

run_in_batch() {
  local id
  for id in $(printf '%s\n' "${!SAMPLES_LIST[@]}" | sort); do
    # In merged mode the per-sample pass is rdist-only, so use rdist resume;
    # in per-sample mode the per-sample pass is the full pass (gbcov+rdist),
    # so use the full resume check.
    local resume_fn=is_resumable
    [ "$SKIP_GBCOV" -eq 1 ] && [ "$RUN_GBCOV" -eq 1 ] && resume_fn=is_resumable_rdist
    if "$resume_fn" "$id"; then
      info "  [skip] $id (already has rseqc outputs)"
      continue
    fi
    BAM="${SAMPLES_LIST[$id]}"
    LABEL="$id"
    stage_rseqc "$SKIP_GBCOV" &
    PIDS+=("$!")
    if (( ${#PIDS[@]} >= MAX_JOBS )); then
      wait "${PIDS[0]}" 2>/dev/null || true
      PIDS=("${PIDS[@]:1}")
    fi
  done
  for p in "${PIDS[@]}"; do wait "$p" 2>/dev/null || true; done
  return 0
}
run_in_batch


#--- post-run cleanup: remove baked genome-wide BED12 intermediates ----------
# read_distribution / geneBody_coverage have now consumed every reference we
# baked, so the (often hundreds-of-MB) genome-wide BED12 intermediates are
# removed to keep the out dir lean. Only files bake_bed12 produced THIS run
# are touched (BAKED_BEDS) — a ready housekeeping BED (GBCOV_REFBED) is a
# user asset and is never in this list. Resume-safe: the .bed12 is a pure
# function of the GTF, so re-baking on a future run is fast and identical.
if [[ ${#BAKED_BEDS[@]} -gt 0 ]]; then
  for _bed in "${BAKED_BEDS[@]}"; do
    if [[ -f "$_bed" ]]; then
      rm -f "$_bed"
      info "removed baked BED12 intermediate (cleaned up after use): $_bed"
    fi
  done
  unset _bed
fi

ELAPSED=$(( SECONDS - START ))
failed=""
[[ -f "$OUT/logs/_failed.samples" ]] && failed="$(sort -u "$OUT/logs/_failed.samples" 2>/dev/null | tr '\n' ' ')"
SUMMARY="$(
  {
    echo
    info "=== SUMMARY ==="
    info "elapsed : ${ELAPSED}s"
    info "species : $SPECIES"
    info "out dir : $OUT"
    info "samples : ${#SAMPLES_LIST[@]}"
    info "failed  : ${failed:-none}"
    info "rseqc   : $OUT/rseqc/  (geneBody_coverage, read_distribution; baked BED12 intermediates removed after use)"
    info "logs    : $OUT/logs/  (per-sample <label>.rseqc.log)"
    info "DONE"
  } 2>&1
)"
# NOTE: the info() lines above already appended to $RUN_LOG (via _log_to_file)
# during capture, so we only echo to stdout here — no second append.
printf '%s\n' "$SUMMARY"
exit 0
