#!/bin/bash
# compute_region_coverage.sh
#
# Per-region coverage for a set of samples (sorted pass BAMs) over a BED of
# target regions, using `samtools bedcov` (one row per region per sample).
#
# Produces 3 CSVs in <output_dir>:
#   <prefix>_region_coverage.csv         long: chrom start end len sample bases cov
#   <prefix>_region_coverage_matrix.csv  one row per region, one column per sample
#   <prefix>_region_coverage_summary.csv total mapped bases + mean coverage per sample
#
# Coverage = mapped base-pairs / region length (weighted by region size).
# Supplementary alignments (0x800) are included (bedcov default).

set -e -o pipefail
on_exit () { echo "This script has exited in error (on line $1)" >&2; }
trap 'on_exit $LINENO' ERR

usage() {
	echo "Usage: $0 --bed <regions.bed> --bam-dir <dir> [--prefix <name>] [--out <dir>] [--glob GLOB]"
	echo ""
	echo "Arguments:"
	echo "  --bed     <file>    BED file with target regions (chrom start end [name])"
	echo "  --bam-dir <dir>     Folder containing all sample BAMs (e.g. .../mapping_SOL0056)"
	echo "  --prefix <name>     Output file prefix (default: basename of --bam-dir)"
	echo "  --out <dir>         Output directory (default: current directory)"
	echo "  --glob GLOB         Pattern of sample BAMs inside --bam-dir (default: *.pass.sorted.bam)"
	echo "  -h                  Show this help"
	exit 0
}

########################################## CHECK DEPENDENCIES #################################
if [[ -z "${CONDA_PREFIX:-}" ]] || [[ "${CONDA_DEFAULT_ENV:-}" == "base" ]]; then
	echo "Please activate your Conda environment first."
	echo "-> conda activate basecaller"
	exit 1
fi
for tool in samtools awk python3; do
	if ! command -v "$tool" &> /dev/null; then
		echo "Error: '$tool' is required." >&2
		exit 1
	fi
done
##############################################################################

BED=""
BAM_DIR=""
PREFIX=""
OUT="."
GLOB="*.pass.sorted.bam"

while [[ $# -gt 0 ]]; do
	case "$1" in
	--bed)          BED="$2";     shift 2 ;;
	--bam-dir)      BAM_DIR="$2"; shift 2 ;;
	--prefix)       PREFIX="$2";  shift 2 ;;
	--out)          OUT="$2";     shift 2 ;;
	--glob)         GLOB="$2";    shift 2 ;;
	-h|--help)      usage ;;
	*)              echo "Unknown option: $1" >&2; usage ;;
	esac
done

[[ -z "$BED" ]]     && { echo "Error: --bed is required"   >&2; usage; }
[[ -z "$BAM_DIR" ]] && { echo "Error: --bam-dir is required" >&2; usage; }
[[ -f "$BED" ]]     || { echo "Error: BED file not found: $BED" >&2; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "Error: BAM dir not found: $BAM_DIR" >&2; exit 1; }
[[ -n "$PREFIX" ]] || PREFIX="$(basename "$BAM_DIR")"
mkdir -p "$OUT"

# Collect the sample BAMs (sorted for stable column order).
BAMS=( $(cd "$BAM_DIR" && ls $GLOB 2>/dev/null | sort) )
if [[ ${#BAMS[@]} -eq 0 ]]; then
	echo "Error: no BAM files matching '$GLOB' in $BAM_DIR" >&2
	exit 1
fi
echo "Found ${#BAMS[@]} sample BAM(s) in $BAM_DIR:"
printf '  %s\n' "${BAMS[@]}"
echo

LONG_TSV="$OUT/${PREFIX}_region_coverage_long.tsv"
: > "$LONG_TSV"

for b in "${BAMS[@]}"; do
	sample="${b%.bam}"
	sample="${sample%.pass.sorted}"
	bam="$BAM_DIR/$b"
	[[ -f "$bam" ]] || { echo "skip (not a file): $b"; continue; }

	# Index if missing (bedcov wants sorted, indexed input).
	if [[ ! -f "$bam.bai" ]]; then
		echo "indexing $b"
		samtools index -@8 "$bam"
	fi

	# Verify every BED chromosome exists in this BAM, else bedcov silently
	# emits zero rows and looks like the script "worked".
	bam_chroms="$(samtools view -H "$bam" | awk '/^@SQ/ {s=$2; sub(/^SN:/,"",s); sub(/^[^:]*:/,"",s); print s}')"
	bed_chroms="$(awk 'NR>1 {if(substr($1,1,1)=="#")next; print $1}' "$BED" | sort -u)"
	while IFS= read -r c; do
		[[ -z "$c" ]] && continue
		if ! printf '%s\n' "$bam_chroms" | grep -qx -- "$c"; then
			echo "ERROR: chromosome '$c' (from BED) not found in $b" >&2
			exit 1
		fi
	done <<< "$bed_chroms"

	echo "-- $sample: $(date +%H:%M:%S)"
	# bedcov output: chrom  start  end  total_mapped_bases
	# We add len(end-start) as col4 and keep sample as col5 so the downstream
	# python can compute cov = bases/len.
	samtools bedcov "$BED" "$bam" \
		| awk -v s="$sample" '{print $1"	"$2"	"$3"	"($3-$2)"	"s"	"$4}' >> "$LONG_TSV"
done

n_regions=$(wc -l < "$LONG_TSV")
echo "long table: $n_regions region x sample rows -> $LONG_TSV"

python3 - "$LONG_TSV" "$PREFIX" "$OUT" <<'PYEOF'
import csv, sys
long_tsv, prefix, out = sys.argv[1], sys.argv[2], sys.argv[3]

with open(long_tsv) as f:
    rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
# cols: 0 chrom, 1 start, 2 end, 3 len, 4 sample, 5 bases

regions = []
for r in rows:
    key = (r[0], int(r[1]), int(r[2]))
    if key not in regions:
        regions.append(key)
samples = []
for r in rows:
    if r[4] not in samples:
        samples.append(r[4])
samples.sort()

cov = {(r[0], int(r[1]), int(r[2]), r[4]) : int(r[5]) for r in rows}

def covval(ch, st, en, s):
    b = cov.get((ch, st, en, s), 0)
    L = en - st
    return (b / L) if L > 0 else 0.0

# long CSV
with open(f"{out}/{prefix}_region_coverage.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["chrom", "start", "end", "len", "sample", "bases", "cov"])
    for (ch, st, en) in regions:
        for s in samples:
            b = cov.get((ch, st, en, s), 0)
            w.writerow([ch, st, en, en - st, s, b, f"{covval(ch, st, en, s):.4f}"])

# matrix CSV
with open(f"{out}/{prefix}_region_coverage_matrix.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["chrom", "start", "end", "len"] + samples)
    for (ch, st, en) in regions:
        w.writerow([ch, st, en, en - st] + [f"{covval(ch, st, en, s):.4f}" for s in samples])

# summary CSV
with open(f"{out}/{prefix}_region_coverage_summary.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["sample", "target_length_bps", "target_length_MB", "total_mapped_bases", "mean_coverage"])
    tot_l = sum(en - st for (ch, st, en) in regions)
    for s in samples:
        tot_b = sum(cov.get((ch, st, en, s), 0) for (ch, st, en) in regions)
        w.writerow([s, tot_l, f"{tot_l/1e6:.4f}", tot_b, f"{(tot_b/tot_l if tot_l else 0):.2f}"])

print("wrote: %s_region_coverage.csv" % prefix)
print("wrote: %s_region_coverage_matrix.csv" % prefix)
print("wrote: %s_region_coverage_summary.csv" % prefix)
PYEOF

echo
echo "Long TSV (intermediate): $LONG_TSV"
