#!/usr/bin/env bash
#
# build_star_index.sh — Build a STAR genome index (one-time step)
#
# One-time step for the MERCURIUS 3' RNA-seq (BRB-seq/DRUG-seq) pipeline:
# build the STAR genome index used by brbseq_pipeline.sh.
#
# You supply the reference FASTA and GTF; no genome is downloaded here.
# Use the primary assembly FASTA (no _sm/_rm masking) and the standard
# GTF (no chr/ab initio suffixes), with consistent chromosome naming
# between the two files.
#
# ERCC spike-ins (optional): pass --ercc-fasta and --ercc-gtf together and
# the script appends the ERCC sequences to the genome FASTA and the ERCC
# gene/exon entries to the genome GTF before building the index. The
# combined files keep your original names with a _with_ERCC suffix
# (e.g. genome.fa -> genome_with_ERCC.fa) so it is always obvious which
# reference the index was built from. The script errors out if only one
# of the two ERCC files is given.
#
# All reference files you passed in are moved into the index directory
# (original names kept, .gz decompressed in place there), so the index
# folder is self-contained and its provenance is readable at a glance:
#
#   <outdir>/
#     Genome  SA  SAindex  genomeParameters.txt  Log.out
#     genome.fa                          (decompressed your genome.fa)
#     annotation.gtf                     (decompressed your annotation.gtf)
#     (+ with the ERCC pair:)
#     genome_with_ERCC.fa                (genome.fa + ERCC92.fa)
#     annotation_with_ERCC.gtf           (annotation.gtf + ERCC92.gtf)
#
# Typical usage:
#   ./build_star_index.sh \
#       --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#       --gtf  Homo_sapiens.GRCh38.108.gtf \
#       --ercc-fasta ERCC92.fasta \
#       --ercc-gtf   ERCC92.gtf \
#       --outdir GRCh38_STAR_index \
#       --threads 8
#   # without spike-ins:
#   ./build_star_index.sh --fasta genome.fa --gtf annot.gtf --outdir idx --threads 8
#   # large genomes (e.g. some plant genomes):
#   ./build_star_index.sh --fasta genome.fa --gtf annot.gtf --ercc-fasta ERCC92.fa \
#       --ercc-gtf ERCC92.gtf --large-genome
#
set -uo pipefail

usage() {
  cat <<'EOF'
build_star_index.sh — build a one-time STAR genome index for the BRB-seq pipeline

PREREQUISITE
  STAR must be on PATH — activate the brb conda environment first, e.g.:
      conda activate brbseq

REQUIRED
  --fasta FILE        Reference genome FASTA (.fa / .fa.gz), primary assembly
  --gtf FILE          GTF annotation (.gtf / .gtf.gz)
  --outdir DIR        Output directory for the STAR index (genomeDir)

OPTIONAL
  --ercc-fasta FILE   ERCC spike-in FASTA to append to the reference
                      (must be used together with --ercc-gtf)
  --ercc-gtf FILE     ERCC spike-in GTF to append to the annotation
                      (must be used together with --ercc-fasta).
                      Combined files are named <original>_with_ERCC.<ext>
  --threads N         STAR --runThreadN (default: 8)
  --sjdbOverhang N    read 2 length minus 1 (default: 89, i.e. for 90 bp R2)
  --large-genome      add --limitSjdbInsertNsj 2000000 and
                      --limitGenomeGenerateRAM 50000000000 (for large genomes)
  -h, --help          show this help
EOF
  exit 0
}

FASTA="" GTF="" OUTDIR="" ERCC_FASTA="" ERCC_GTF="" THREADS=8 SJDB_OVERHANG=89
LARGE_GENOME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --ercc-fasta) ERCC_FASTA="$2"; shift 2 ;;
    --ercc-gtf) ERCC_GTF="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sjdbOverhang) SJDB_OVERHANG="$2"; shift 2 ;;
    --large-genome) LARGE_GENOME=1; shift ;;
    -h|--help) usage ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$FASTA" && -n "$GTF" && -n "$OUTDIR" ]] || { echo "ERROR: --fasta, --gtf and --outdir are required." >&2; usage; }
if [[ ( -n "$ERCC_FASTA" && -z "$ERCC_GTF" ) || ( -z "$ERCC_FASTA" && -n "$ERCC_GTF" ) ]]; then
  echo "ERROR: --ercc-fasta and --ercc-gtf must be passed together (got only one)." >&2
  usage
fi

# --- sanity checks -----------------------------------------------------------
for f in "$FASTA" "$GTF" ${ERCC_FASTA:+"$ERCC_FASTA"} ${ERCC_GTF:+"$ERCC_GTF"}; do
  [[ -f "$f" ]] || { echo "ERROR: input file not found: $f" >&2; exit 1; }
done
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || { echo "ERROR: --threads must be a positive integer." >&2; exit 1; }
[[ "$SJDB_OVERHANG" =~ ^[0-9]+$ && "$SJDB_OVERHANG" -ge 20 ]] || { echo "ERROR: --sjdbOverhang must be >= 20 (read 2 length minus 1)." >&2; exit 1; }
if [[ -e "$OUTDIR" && -n "$(ls -A "$OUTDIR" 2>/dev/null)" ]]; then
  echo "ERROR: outdir '$OUTDIR' already exists and is not empty. STAR refuses to build into an existing index directory — pick a new directory." >&2
  exit 1
fi

# --- locate STAR -------------------------------------------------------------
# The script expects STAR to be on PATH — activate the brb conda env first
# (e.g. conda activate brbseq), or make sure STAR is otherwise reachable
# on PATH.
STAR_BIN="$(command -v STAR || true)"
[[ -n "$STAR_BIN" ]] || { echo "ERROR: STAR not found on PATH. Activate the brb conda environment first (e.g. 'conda activate brbseq') and re-run." >&2; exit 1; }
[[ -x "$STAR_BIN" ]] || { echo "ERROR: STAR binary not executable: $STAR_BIN" >&2; exit 1; }
echo "Using STAR:  $STAR_BIN ($("$STAR_BIN" --version 2>&1 | head -n1))"

# Build into the final directory directly. There is no temporary workdir:
# the index folder is *created by* this run, so the provenance files (the
# FASTA/GTF STAR actually read) go right next to the index, under their
# real names.
mkdir -p "$OUTDIR" || { echo "ERROR: cannot create output directory: $OUTDIR" >&2; exit 1; }
rm -rf -- "${OUTDIR:?}"
mkdir -p "$OUTDIR"

# strip_suffix FILE → echoes the basename without a trailing .gz
strip_suffix() {
  local b
  b="$(basename "$1")"
  printf '%s' "${b%.gz}"
}

# materialize SRC DEST
#   If SRC is .gz, decompress into DEST; otherwise move SRC into DEST.
#   The moved/decompressed file keeps SRC's original filename, just placed
#   inside the index directory (which is exactly where it should live).
materialize() {
  local src="$1" dst
  dst="$OUTDIR/$(strip_suffix "$src")"
  if [[ -L "$src" ]]; then
    local real
    real="$(readlink -f "$src")"
    mv -f -- "$real" "$dst"
  elif [[ "$(readlink -f "$src")" == "$(readlink -f "$OUTDIR")/"* ]]; then
    mv -f -- "$src" "$dst"
  else
    case "$src" in
      *.gz) gzip -dc -- "$src" > "$dst" || { echo "ERROR: failed to decompress $src" >&2; exit 1; } ;;
      *)    mv -f -- "$src" "$dst" || { echo "ERROR: failed to move $src into $OUTDIR (check disk space / permissions)" >&2; exit 1; } ;;
    esac
  fi
}

echo
echo "Placing reference files in the index directory (original names kept):"
materialize "$FASTA"
materialize "$GTF"
REF_FASTA="$OUTDIR/$(strip_suffix "$FASTA")"
GTF_PLAIN="$OUTDIR/$(strip_suffix "$GTF")"

# --- merge ERCC spike-ins, if given ------------------------------------------
if [[ -n "$ERCC_FASTA" && -n "$ERCC_GTF" ]]; then
  materialize "$ERCC_FASTA"
  materialize "$ERCC_GTF"
  # Names keep the original basenames with a _with_ERCC suffix, so the
  # merged reference is always recognizable at a glance:
  #   genome.fa  + ERCC92.fa  ->  genome_with_ERCC.fa
  #   annot.gtf  + ERCC92.gtf ->  annot_with_ERCC.gtf
  REF_BASE="$(strip_suffix "$FASTA")"; REF_BASE="${REF_BASE%.fa}"
  GTF_BASE="$(strip_suffix "$GTF")";  GTF_BASE="${GTF_BASE%.gtf}"
  REF_FASTA="$OUTDIR/${REF_BASE}_with_ERCC.fa"
  GTF_PLAIN="$OUTDIR/${GTF_BASE}_with_ERCC.gtf"
  # Append ERCC spike-in sequences to the reference FASTA (plain concatenation;
  # STAR is fine with the two blocks back to back).
  cat "$OUTDIR/$(strip_suffix "$FASTA")" "$OUTDIR/$(strip_suffix "$ERCC_FASTA")" > "$REF_FASTA" \
    || { echo "ERROR: failed to build $REF_FASTA" >&2; exit 1; }
  # Append ERCC gene/transcript/exon records to the annotation GTF.
  cat "$OUTDIR/$(strip_suffix "$GTF")" "$OUTDIR/$(strip_suffix "$ERCC_GTF")" > "$GTF_PLAIN" \
    || { echo "ERROR: failed to build $GTF_PLAIN" >&2; exit 1; }
  echo "ERCC merged: $(basename "$FASTA") + $(basename "$ERCC_FASTA")   ->   $(basename "$REF_FASTA")"
  echo "ERCC merged: $(basename "$GTF") + $(basename "$ERCC_GTF")        ->   $(basename "$GTF_PLAIN")"
fi

# --- run STAR genomeGenerate --------------------------------------------------
STAR_ARGS=(--runMode genomeGenerate
           --runThreadN "$THREADS"
           --genomeDir "$OUTDIR"
           --genomeFastaFiles "$REF_FASTA"
           --sjdbGTFfile "$GTF_PLAIN"
           --sjdbOverhang "$SJDB_OVERHANG")
if [[ "$LARGE_GENOME" -eq 1 ]]; then
  STAR_ARGS+=(--limitSjdbInsertNsj 2000000 --limitGenomeGenerateRAM 50000000000)
fi

echo
echo "Running: $STAR_BIN ${STAR_ARGS[*]}"
"$STAR_BIN" "${STAR_ARGS[@]}"
star_rc=$?
if [[ "$star_rc" -ne 0 ]]; then
  echo "ERROR: STAR genomeGenerate failed (rc=$star_rc). Index directory with partial contents left at: $OUTDIR" >&2
  exit "$star_rc"
fi

# --- verification ------------------------------------------------------------
for f in Genome SA SAindex genomeParameters.txt Log.out; do
  [[ -e "$OUTDIR/$f" ]] || { echo "ERROR: expected index file missing: $OUTDIR/$f" >&2; exit 1; }
done
grep -q "finished successfully" "$OUTDIR/Log.out" || { echo "ERROR: STAR Log.out does not report a successful build." >&2; exit 1; }
echo
echo "OK — STAR genome index built in: $OUTDIR"
SJDB_USED="$(awk -F'\t' '$1=="sjdbOverhang"{print $2}' "$OUTDIR/genomeParameters.txt")"
echo "   reference FASTA:  $REF_FASTA"
echo "   annotation GTF:   $GTF_PLAIN"
echo "   sjdbOverhang in use: ${SJDB_USED:-<not found>}"
echo "   (provenance: $OUTDIR/genomeParameters.txt — exact command used to build it)"
echo "   (all input files live inside the index directory, under their original names)"
