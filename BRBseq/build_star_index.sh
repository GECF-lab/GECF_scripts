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
    *) echo "ERROR: I don't recognise the option '$1'. The full list of options is printed by:  build_star_index.sh --help" >&2; usage ;;
  esac
done

if [[ -z "$FASTA" || -z "$GTF" || -z "$OUTDIR" ]]; then
  echo "ERROR: This script needs three things to run, and one or more of them is missing: --fasta (the reference genome file), --gtf (the gene annotation file), and --outdir (where to save the built index). The full list of options, with a worked example, is printed by:  build_star_index.sh --help" >&2
  usage
fi
if [[ ( -n "$ERCC_FASTA" && -z "$ERCC_GTF" ) || ( -z "$ERCC_FASTA" && -n "$ERCC_GTF" ) ]]; then
  echo "ERROR: The two ERCC spike-in options come as a pair: --ercc-fasta and --ercc-gtf must be given together (one adds the spike-in DNA sequences to your genome file, the other adds the matching gene entries to your annotation — either one on its own would leave a mismatched reference). Give both, or omit both, and re-run." >&2
  usage
fi

# --- sanity checks -----------------------------------------------------------
for f in "$FASTA" "$GTF" ${ERCC_FASTA:+"$ERCC_FASTA"} ${ERCC_GTF:+"$ERCC_GTF"}; do
  [[ -f "$f" ]] || { echo "ERROR: One of the reference files you supplied could not be found: $f. Check that the file name and folder are spelled exactly right (file names are case-sensitive), then re-run." >&2; exit 1; }
done
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || { echo "ERROR: The value given for --threads must be a whole number greater than zero (it sets how many CPU cores STAR is allowed to use at once). For example: --threads 8. Re-run with a valid number." >&2; exit 1; }
[[ "$SJDB_OVERHANG" =~ ^[0-9]+$ && "$SJDB_OVERHANG" -ge 20 ]] || { echo "ERROR: The value given for --sjdbOverhang must be a whole number of at least 20. In plain terms: this number should be your R2 (genomic) read length minus 1 — for the 90 bp R2 reads used by BRB-seq, that is 89, which is also the default this script already uses, so you can usually just leave --sjdbOverhang out entirely. Re-run with a valid number (or omit the option) for a non-standard read length." >&2; exit 1; }
if [[ -e "$OUTDIR" && -n "$(ls -A "$OUTDIR" 2>/dev/null)" ]]; then
  echo "ERROR: The folder you asked to use as the new index location ($OUTDIR) already exists and already contains files. The index-building tool (STAR) refuses to build into a folder that is not brand-new or empty, to avoid mixing a fresh build with old files. Delete or rename that folder first, or pick a different folder name for --outdir, then re-run." >&2
  exit 1
fi

# --- locate STAR -------------------------------------------------------------
# The script expects STAR to be on PATH — activate the brb conda env first
# (e.g. conda activate brbseq), or make sure STAR is otherwise reachable
# on PATH.
STAR_BIN="$(command -v STAR || true)"
if [[ -z "$STAR_BIN" ]]; then
  echo "ERROR: The STAR program (the tool that both builds this index and, later, aligns your reads to it) is not installed or not currently active where this script can reach it. It ships inside the 'brbseq' Conda environment: activate it with 'conda activate brbseq' (create it first with 'conda env create -f brbseq-env.yaml' if you have never made it), then re-run this script. Everything else about the reference files you supplied is fine — this is purely an environment setup issue." >&2
  exit 1
fi
if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: A STAR program was found at $STAR_BIN, but it is not currently permitted to be run (its file permissions are set to 'not executable'). This is a file-permission problem on this computer, not a problem with your reference files. If you set up this install yourself, check the file's permissions; otherwise, re-install the 'brbseq' Conda environment and try again." >&2
  exit 1
fi
echo "Using STAR:  $STAR_BIN ($("$STAR_BIN" --version 2>&1 | head -n1))"

# Build into the final directory directly. There is no temporary workdir:
# the index folder is *created by* this run, so the provenance files (the
# FASTA/GTF STAR actually read) go right next to the index, under their
# real names.
mkdir -p "$OUTDIR" || { echo "ERROR: Could not create or open the output folder ($OUTDIR) that was to hold the new STAR index. This is a file-system permission problem on this computer (for example the parent folder is not writable by whoever is running this script), not a problem with your reference files. Check permissions on that path, or choose a different --outdir location, and re-run." >&2; exit 1; }
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
      *.gz) gzip -dc -- "$src" > "$dst" || { echo "ERROR: The zipped reference file '$src' could not be read/unpacked into the index folder. This is usually a file-system permission or disk-space problem on this computer, or the .gz file itself is damaged. Check those two things and re-run." >&2; exit 1; } ;;
      *)    mv -f -- "$src" "$dst" || { echo "ERROR: The reference file '$src' could not be moved into the index folder ($OUTDIR). Check that the computer has enough free disk space and that the file is not locked or protected, then re-run." >&2; exit 1; } ;;
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
    || { echo "ERROR: Could not combine your main reference genome and the ERCC spike-in sequences into the single merged file ($REF_FASTA) that STAR needs. This points to a file-system problem (permissions or disk space) on this computer, not to a problem with the two source files themselves (which are both readable, since they were placed into the index folder just now). Check permissions/disk space and re-run." >&2; exit 1; }
  # Append ERCC gene/transcript/exon records to the annotation GTF.
  cat "$OUTDIR/$(strip_suffix "$GTF")" "$OUTDIR/$(strip_suffix "$ERCC_GTF")" > "$GTF_PLAIN" \
    || { echo "ERROR: Could not combine your main gene annotation and the ERCC spike-in annotation into the single merged file ($GTF_PLAIN) that STAR needs. This points to a file-system problem (permissions or disk space) on this computer, not to a problem with the two source files themselves (which are both readable, since they were placed into the index folder just now). Check permissions/disk space and re-run." >&2; exit 1; }
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
  echo "ERROR: STAR finished the index build, but it reported that something went wrong part-way through (its own log files hold the specific reason). The folder $OUTDIR was left in place with everything STAR did manage to build before it stopped — you can open the file Log.out inside it (or Log.final.out, if present) to see STAR's own words on what failed, which usually pinpoints the real cause (for example an unrecognised chromosome line in your FASTA, or a gene/exon line in your GTF that STAR does not understand). Fix whatever that log points at, delete or rename the half-built $OUTDIR folder (STAR refuses to build into a folder that already has files in it), and re-run this script. Your original reference files are safe — they have only been copied into $OUTDIR, never modified in their original location." >&2
  exit "$star_rc"
fi

# --- verification ------------------------------------------------------------
for f in Genome SA SAindex genomeParameters.txt Log.out; do
  [[ -e "$OUTDIR/$f" ]] || { echo "ERROR: STAR said it finished, but one of the files a fully working STAR index must contain is missing from $OUTDIR (expected: $f). This means the build did not actually complete on STAR's side, even if no clear error was printed — the most reliable place to look for the real reason is STAR's own log file, Log.out, inside $OUTDIR. Read through it, fix whatever it points at (common culprits: an odd line in the FASTA or GTF, insufficient disk space), delete or rename the half-built folder, and re-run." >&2; exit 1; }
done
grep -q "finished successfully" "$OUTDIR/Log.out" || { echo "ERROR: The file inside $OUTDIR where STAR records how its build went (Log.out) does not say the build finished successfully — which means, even though the expected files are all there, STAR's own report does not confirm a clean finish. Open Log.out and read it for the specific reason, fix what it points at, delete or rename the half-built folder, and re-run." >&2; exit 1; }
echo
echo "OK — STAR genome index built in: $OUTDIR"
SJDB_USED="$(awk -F'\t' '$1=="sjdbOverhang"{print $2}' "$OUTDIR/genomeParameters.txt")"
echo "   reference FASTA:  $REF_FASTA"
echo "   annotation GTF:   $GTF_PLAIN"
echo "   sjdbOverhang in use: ${SJDB_USED:-<not found>}"
echo "   (provenance: $OUTDIR/genomeParameters.txt — exact command used to build it)"
echo "   (all input files live inside the index directory, under their original names)"
