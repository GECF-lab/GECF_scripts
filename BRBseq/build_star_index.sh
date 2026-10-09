#!/usr/bin/env bash
#
# build_star_index.sh — Build a STAR genome index (one-time step)
#
# Implements section "1.1. Building a genome index (one-time)" of the
# Alithea MERCURIUS 3' RNA-seq (BRB-seq/DRUG-seq) Data Analysis User Guide.
#
# You supply the reference FASTA and GTF; no genome is downloaded here.
# Use the primary assembly FASTA (no _sm/_rm masking) and the standard
# GTF (no chr/ab initio suffixes), with consistent chromosome naming
# between the two files.
#
# ERCC spike-ins (optional): pass --ercc-fasta and --ercc-gtf together and
# the script appends the ERCC sequences to the genome FASTA and the ERCC
# gene/exon entries to the genome GTF before building the index. The script
# errors out if only one of the two files is given.
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
#   ./build_star_index.sh --fasta genome.fa --gtf annot.gtf --ercc-fasta ERCC92.fasta \
#       --ercc-gtf ERCC92.gtf --large-genome
#
set -euo pipefail

usage() {
  cat <<'EOF'
build_star_index.sh — build a STAR genome index for BRB-seq (MERCURIUS guide 1.1)

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
                      (must be used together with --ercc-fasta)
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
[[ ( -n "$ERCC_FASTA" && -n "$ERCC_GTF" ) ]] || [[ ( -z "$ERCC_FASTA" && -z "$ERCC_GTF" ) ]] || { echo "ERROR: --ercc-fasta and --ercc-gtf must be passed together (got only one)." >&2; usage; }

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

# --- prepare reference files (decompress .gz if needed) ----------------------
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

decompress() {
  local src="$1" dst="$2"
  case "$src" in
    *.gz) gzip -dc "$src" > "$dst" ;;
    *) cp "$src" "$dst" ;;
  esac
  echo "$dst"
}

REF_FASTA="$(decompress "$FASTA" "$WORKDIR/reference.fa")"
GTF_PLAIN="$(decompress "$GTF" "$WORKDIR/annotation.gtf")"

if [[ -n "$ERCC_FASTA" && -n "$ERCC_GTF" ]]; then
  ERCC_FASTA_PLAIN="$(decompress "$ERCC_FASTA" "$WORKDIR/ercc.fa")"
  ERCC_GTF_PLAIN="$(decompress "$ERCC_GTF" "$WORKDIR/ercc.gtf")"
  # Append ERCC spike-in sequences to the reference FASTA (plain concatenation;
  # STAR is fine with the two blocks back to back).
  REF_FASTA="$WORKDIR/reference_with_ERCC.fa"
  cat "$WORKDIR/reference.fa" "$ERCC_FASTA_PLAIN" > "$REF_FASTA"
  # Append ERCC gene/transcript/exon records to the annotation GTF.
  GTF_PLAIN="$WORKDIR/annotation_with_ERCC.gtf"
  cat "$WORKDIR/annotation.gtf" "$ERCC_GTF_PLAIN" > "$GTF_PLAIN"
  echo "ERCC: appended $(awk -F'\t' '$3=="gene"{n++} END{print n+0}' "$ERCC_GTF_PLAIN") gene records and $(grep -c '^>' "$ERCC_FASTA_PLAIN") sequences to the reference."
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

echo "Running: $STAR_BIN ${STAR_ARGS[*]}"
"$STAR_BIN" "${STAR_ARGS[@]}"

# --- verification ------------------------------------------------------------
for f in Genome SA SAindex genomeParameters.txt Log.out; do
  [[ -e "$OUTDIR/$f" ]] || { echo "ERROR: expected index file missing: $OUTDIR/$f" >&2; exit 1; }
done
grep -q "finished successfully" "$OUTDIR/Log.out" || { echo "ERROR: STAR Log.out does not report a successful build." >&2; exit 1; }
echo
echo "OK — STAR genome index built in: $OUTDIR"
SJDB_USED="$(awk -F'\t' '$1=="sjdbOverhang"{print $2}' "$OUTDIR/genomeParameters.txt")"
echo "   sjdbOverhang in use: ${SJDB_USED:-<not found>}"
echo "   (provenance: $OUTDIR/genomeParameters.txt — exact command used to build it)"

# --- move the original input files into the index directory ------------------
# Move the actual --fasta/--gtf (and --ercc-fasta/--ercc-gtf, if given)
# files into the freshly built index folder, keeping their original names,
# so the index folder is self-contained. Files already inside the outdir
# are left in place.
for f in "$FASTA" "$GTF" ${ERCC_FASTA:+"$ERCC_FASTA"} ${ERCC_GTF:+"$ERCC_GTF"}; do
  base="$(basename "$f")"
  if [[ "$(dirname "$(readlink -f "$f")")" != "$(readlink -f "$OUTDIR")" ]]; then
    mv -f "$f" "$OUTDIR/$base"
    echo "Moved: $f -> $OUTDIR/$base"
  fi
done
echo "   (input FASTA/GTF are now inside the index directory)"
