#!/usr/bin/env bash

set -e

on_exit () {
	echo "This script has exited in error (on line $1)"
}

trap 'on_exit $LINENO' ERR

# --- Helper function ---
usage() {
	echo "Usage: $0 [-o output] [-p pattern] <input> [<input> ...]"
	echo ""
	echo "Arguments:"
	echo "  <input>                Input to process: a .bam/.ubam file, OR a folder"
	echo "                         containing .bam/.ubam files. May be given multiple times."
	echo ""
	echo "Each input file is reported as a 2-column CSV (Metric,Value) with these rows in order:"
	echo "  Number of reads, Total bases (Gb), Mean read length, Median read length,"
	echo "  Read length N50, Mean read quality, Median read quality, Active channel,"
	echo "  and the 1st and 5th longest reads (combined '<len> bases (<q> Q)')."
	echo ""
	echo "Options:"
	echo "  -o output            Optional CSV output target."
	echo "                       -> a file path : write that one input's CSV there (single input only)"
	echo "                       -> a directory : write one <name>.csv per input file inside it"
	echo "                       (default when omitted: print CSV to stdout, one block per input file)"
	echo "  -p pattern           Optional: file pattern to match when expanding folders"
	echo "                       (find -name, e.g. '*.pass.bam'). Default: '*.bam' and '*.ubam'"
	echo "  -h                   Show this help message"
	echo ""
	echo "Note: options must come before the input file(s)/folder(s)."
	echo ""
	echo "Examples:"
	echo "  $0 run1/sample_name1.pass.bam                    # CSV to stdout (pipe with '>' to save)"
	echo "  $0 -o stats.csv run1/sample_name1.pass.bam       # single input -> stats.csv"
	echo "  $0 -o outdir /data1/ont_raw_data/SOL0057         # folder -> outdir/<name>.csv per file"
	echo "  $0 -p '*.pass.bam' -o outdir runA /data1/ont_raw_data/SOL0045"
	exit 0
}

########################################## CHECK DEPENDENCIES ###########################################

# Check if conda environment is activated
if [[ -z "$CONDA_PREFIX" ]] || [[ "$CONDA_DEFAULT_ENV" == "base" ]]; then
        echo "Please activate your Conda environment first."
        echo "-> conda activate basecaller"
        exit 1
fi

if ! command -v samtools &> /dev/null; then
	echo "Error: 'samtools' is required to parse uBAM files." >&2
	exit 1
fi

#########################################################################################################

OUT=""
PATTERN=""
while getopts "o:p:h" opt; do
	case $opt in
		o) OUT="$OPTARG" ;;
		p) PATTERN="$OPTARG" ;;
		h) usage ;;
		*) usage ;;
	esac
done
shift $((OPTIND - 1))

if [ $# -eq 0 ]; then
	echo "Error: no input file or folder provided." >&2
	usage
fi

# --- Expand each input into a concrete list of BAM/uBAM files ---
declare -A SEEN
FILES=()
add_file() {
	local f="$1"
	local key="${f#/}"
	if [ -z "${SEEN[$key]:-}" ]; then
		SEEN[$key]=1
		FILES+=("$f")
	fi
}

for input in "$@"; do
	if [ -f "$input" ]; then
		# An explicitly given file is always processed, regardless of -p.
		add_file "$input"
	elif [ -d "$input" ]; then
		if [ -n "$PATTERN" ]; then
			while IFS= read -r -d '' f; do
				add_file "$f"
			done < <(find "$input" -type f -name "$PATTERN" -print0 | sort -z)
		else
			while IFS= read -r -d '' f; do
				add_file "$f"
			done < <(find "$input" -type f \( -name '*.bam' -o -name '*.ubam' \) -print0 | sort -z)
		fi
	else
		echo "Error: '$input' is neither an existing file nor a folder." >&2
		exit 1
	fi
done

if [ ${#FILES[@]} -eq 0 ]; then
	echo "Error: no .bam/.ubam files found in the given input(s)." >&2
	exit 1
fi

NUMFILES=${#FILES[@]}

# Resolve the CSV output target for a given input file.
#   OF (set in the caller) = ""            -> write to stdout
#   OF (set in the caller) = "<path>.csv"  -> write to that file (overwrite)
resolve_out() {
	local infile="$1"
	OF=""
	# Name for the per-file CSV (bare basename minus .bam/.ubam, then .csv)
	local name
	name=$(basename "$infile")
	name="${name%.ubam}"
	name="${name%.bam}.csv"

	[ -z "$OUT" ] && return 0                          # no -o -> stdout (OF already "")

	case "$OUT" in
	*.csv)                                              # explicit .csv file target
		if [ "$NUMFILES" -gt 1 ]; then
			echo "Error: multiple input files but -o '$OUT' is a single .csv file." >&2
			echo "       Use a directory with -o for multiple inputs." >&2
			exit 1
		fi
		OF="$OUT" ;;
	*/)                                                # explicit directory (trailing slash)
		mkdir -p "$OUT"; OF="${OUT%/}/$name" ;;
	*)                                                  # bare path
		if [ -d "$OUT" ]; then
			OF="$OUT/$name"
		elif [ -f "$OUT" ]; then
			if [ "$NUMFILES" -gt 1 ]; then
				echo "Error: multiple input files but -o '$OUT' is an existing file." >&2
				exit 1
			fi
			OF="$OUT"
		else
			# Default: treat as an output directory and create it.
			mkdir -p "$OUT"; OF="${OUT%/}/$name"
		fi ;;
	esac
}

echo "Found ${NUMFILES} file(s) to process." >&2
echo "Extracting BAM tags (qs, ch) per file..." >&2

# --- Per-file statistics (2-column CSV: Metric,Value) ---
# $1 = path to the BAM file.
# Resolve $1 into OF first ("" for stdout) via resolve_out.
# Returns 0 on success, 1 if the file has no valid reads / cannot be read.
process_file() {
	local input_file="$1"

	# Extract read sequence length, mean Q-score tag (qs:f), and channel tag (ch:i).
	# Bypasses character-by-character calculation, speeding processing. Hardcoded thread count.
	local stream_cmd
	stream_cmd="samtools view -@ 10 -F 2304 \"$input_file\" | awk -F'	' '
	{
		len = length(\$10);
		# Extract Q-score float tag (qs:f:XX)
		match(\$0, /qs:f:[^	]+/);
		qs = (RSTART > 0) ? substr(\$0, RSTART+5, RLENGTH-5) : 0;
		# Extract Channel integer tag (ch:i:XX)
		match(\$0, /ch:i:[^	]+/);
		ch = (RSTART > 0) ? substr(\$0, RSTART+5, RLENGTH-5) : \"unknown\";
		print len \"	\" qs \"	\" ch
	}'"

	# Median Q. Also acts as an empty/unreadable-file probe (awk exits 1 on NR==0).
	local med_q
	if ! med_q=$(eval "$stream_cmd" | awk -F'\t' '{print $2}' | sort -n | awk -v fname="$input_file" '
	{
		scores[NR] = $1
	}
	END {
		if (NR == 0) {
			print "Error: No valid sequencing reads parsed from: " fname > "/dev/stderr"
			exit 1
		}
		if (NR % 2 == 1) printf "%.2f", scores[(NR+1)/2]
		else printf "%.2f", (scores[NR/2] + scores[(NR/2)+1]) / 2
	}'); then
		echo "WARNING: skipping '$input_file' (no valid sequencing reads)." >&2
		return 1
	fi

	resolve_out "$input_file"

	# For a multi-file stdout dump, prepend a comment line so the blocks are separable.
	if [ -z "$OF" ] && [ "$NUMFILES" -gt 1 ]; then
		echo "# $input_file"
	fi

	# Emit the 10-row, 2-column CSV to OF (file) or stdout.
	local fail_rc=0
	eval "$stream_cmd" | sort -k1,1n | awk -F'\t' -v of="$OF" -v medq="$med_q" '
	function emitrow(s) {  # write one CSV line to of (file) or stdout
		if (of != "") print s > of
		else          print s
	}
	function cell(a, b) { emitrow(a "," b) }
	{
		lens[NR] = $1
		qs[NR] = $2
		total_bases += $1
		q_sum += 10^(-$2/10)
		if ($1 > 40000) over40k++
		if ($1 > 100000) over100k++
		channels[$3]++
	}
	END {
		if (NR == 0) {
			print "Error: No valid sequencing reads parsed from BAM." > "/dev/stderr"
			emitrow("Metric,Value")
			exit 1
		}
		if (NR % 2 == 1) med_len = lens[(NR+1)/2]
		else             med_len = (lens[NR/2] + lens[(NR/2)+1]) / 2
		mean_len = total_bases / NR
		mean_q = -10 * log(q_sum / NR) / log(10)
		half_bases = total_bases / 2; running_sum = 0; n50 = 0
		for (i = NR; i >= 1; i--) {
			running_sum += lens[i]
			if (running_sum >= half_bases) { n50 = lens[i]; break; }
		}

		emitrow("Metric,Value")
		cell("Number of reads", NR)
		cell("Total bases (Gb)", sprintf("%.9f", total_bases / 1e9))
		cell("Mean read length", sprintf("%.1f", mean_len))
		cell("Median read length", sprintf("%.0f", med_len))
		cell("Read length N50", sprintf("%.0f", n50))
		cell("Mean read quality", sprintf("%.1f", mean_q))
		cell("Median read quality", sprintf("%.2f", medq))
		cell("Active channel", length(channels))
		cell("1st of longest reads (length and mean quality score)",
		     sprintf("%.0f bases (%.2f Q)", lens[NR], qs[NR]))
		cell("5th of longest reads (length and mean quality score)",
		     sprintf("%.0f bases (%.2f Q)", lens[NR-4], qs[NR-4]))
		emitrow("")
		emitrow("")
		emitrow("")
		cell("Reads > 40kb", over40k + 0)
		cell("Reads > 100kb", over100k + 0)
	}
' || fail_rc=1

	if [ "$fail_rc" -ne 0 ]; then
		echo "WARNING: no stats produced for '$input_file'." >&2
		return 1
	fi
	return 0
}

# --- Main: write a CSV for every discovered file ---
succeeded=0
skipped=0
for input_file in "${FILES[@]}"; do
	if process_file "$input_file"; then
		succeeded=$((succeeded + 1))
	else
		skipped=$((skipped + 1))
	fi
done

echo "Done: $succeeded CSV(s) written, $skipped file(s) skipped, of ${NUMFILES} discovered." >&2
if [ "$succeeded" -eq 0 ]; then
	echo "Error: no input file produced statistics." >&2
	exit 1
fi
exit 0
