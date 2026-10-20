#!/usr/bin/env bash
set -euo pipefail

script_name="$(basename "$0")"
input_file=""
whitelist_file=""
fastq_file=""
bam_file=""
has_header=false
declare -A seen_samples=()

usage() {
    cat <<EOF
Usage: ${script_name} -i INPUT_FILE -w WHITELIST_FILE -f FASTQ_FILE -b BAM_FILE [-H]

Validate a strict two-column sample/barcode manifest and derive downstream
barcode files.

Required options:
  -i INPUT_FILE     Two-column manifest with sample ID in column 1 and barcode in column 2.
  -w WHITELIST_FILE Output path for the STAR barcode whitelist.
  -f FASTQ_FILE     Output path for the FASTQ demux metadata TSV.
  -b BAM_FILE       Output path for the BAM demux barcode map.
  -H                Treat the first row as a header and skip it.
  -h                Show this help text and exit.
EOF
}

while getopts ":i:w:c:f:b:Hh" opt; do
    case "$opt" in
        i)
            input_file="$OPTARG"
            ;;
        w)
            whitelist_file="$OPTARG"
            ;;
        f)
            fastq_file="$OPTARG"
            ;;
        b)
            bam_file="$OPTARG"
            ;;
        H)
            has_header=true
            ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Option -${OPTARG} requires a value." >&2
            usage >&2
            exit 1
            ;;
        ?)
            echo "Unknown option: -${OPTARG}" >&2
            usage >&2
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))

if [[ $# -ne 0 ]]; then
    echo "Unexpected positional arguments: $*" >&2
    usage >&2
    exit 1
fi

required_args=(
    "$input_file"
    "$whitelist_file"
    "$fastq_file"
    "$bam_file"
)

for required_arg in "${required_args[@]}"; do
    if [[ -z "$required_arg" ]]; then
        echo "All output and input paths are required." >&2
        usage >&2
        exit 1
    fi
done

declare -A seen_filenames=()
declare -A seen_filename_flags=()

path_args=(
    "-i:$input_file"
    "-w:$whitelist_file"
    "-f:$fastq_file"
    "-b:$bam_file"
)

for path_arg in "${path_args[@]}"; do
    arg_flag="${path_arg%%:*}"
    arg_path="${path_arg#*:}"
    arg_name="$(basename "$arg_path")"
    if [[ -n "${seen_filenames[$arg_name]:-}" ]]; then
        echo "Conflicting file names provided to ${script_name}: ${arg_flag} and ${seen_filename_flags[$arg_name]} both use '${arg_name}'. Input and output file names must be unique to avoid overwriting files." >&2
        exit 1
    fi
    seen_filenames["$arg_name"]=1
    seen_filename_flags["$arg_name"]="$arg_flag"
done

if [[ ! -f "$input_file" ]]; then
    echo "Barcode manifest not found: $input_file" >&2
    exit 1
fi

: > "$whitelist_file"
printf 'sample_id\tbarcode\n' > "$fastq_file"
: > "$bam_file"

line_number=0
processed_samples=0

while IFS=$'\t' read -r sample_id barcode extra || [[ -n "$sample_id" || -n "$barcode" || -n "$extra" ]]; do
    line_number=$((line_number + 1))
    if [[ "$has_header" == true && $line_number -eq 1 ]]; then
        continue
    fi
    sample_id="${sample_id%$'\r'}"
    barcode="${barcode%$'\r'}"
    extra="${extra%$'\r'}"
    if [[ -z "$sample_id" && -z "$barcode" && -z "$extra" ]]; then
        echo "Invalid barcode manifest at line ${line_number}: empty lines are not allowed." >&2
        exit 1
    fi
    if [[ -n "$extra" ]]; then
        echo "Invalid barcode manifest at line ${line_number}: expected exactly 2 tab-delimited columns." >&2
        exit 1
    fi
    if [[ -z "$sample_id" || -z "$barcode" ]]; then
        echo "Invalid barcode manifest at line ${line_number}: sample ID and barcode are required." >&2
        exit 1
    fi
    if [[ "$sample_id" == */* ]]; then
        echo "Invalid barcode manifest at line ${line_number}: sample IDs must not contain '/'." >&2
        exit 1
    fi
    if [[ -n "${seen_samples[$sample_id]:-}" ]]; then
        echo "Invalid barcode manifest at line ${line_number}: duplicate sample ID '${sample_id}'." >&2
        exit 1
    fi
    normalized_barcode="$(printf '%s' "$barcode" | tr '[:lower:]' '[:upper:]')"
    if [[ ! "$normalized_barcode" =~ ^[ACGTN]+$ ]]; then
        echo "Invalid barcode manifest at line ${line_number}: barcode must contain only ACGTN characters." >&2
        exit 1
    fi
    seen_samples["$sample_id"]=1
    processed_samples=$((processed_samples + 1))
    printf '%s\n' "$normalized_barcode" >> "$whitelist_file"
    printf '%s\t%s\n' "$sample_id" "$normalized_barcode" >> "$fastq_file"
    printf '%s\t%s\n' "$sample_id" "$normalized_barcode" >> "$bam_file"
done < "$input_file"

if [[ $processed_samples -eq 0 ]]; then
    echo "Invalid barcode manifest: file is empty." >&2
    exit 1
fi