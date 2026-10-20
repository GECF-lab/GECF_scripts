#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <output_dir> <bam_file> <barcode_file> <tag_name>" >&2
    exit 1
fi

output_dir="$1"
bam_file="$2"
barcode_file="$3"
tag_name="$4"

mkdir -p "$output_dir"

while IFS=$'\t' read -r sample_id tag_value extra || [[ -n "$sample_id" || -n "$tag_value" || -n "$extra" ]]; do
    sample_id="${sample_id%$'\r'}"
    tag_value="${tag_value%$'\r'}"
    extra="${extra%$'\r'}"

    # Allow blank and comment lines in the barcode mapping file.
    if [[ -z "$sample_id" && -z "$tag_value" && -z "$extra" ]]; then
        continue
    fi

    if [[ "$sample_id" =~ ^[[:space:]]*# ]]; then
        continue
    fi

    picard FilterSamReads \
        I="${bam_file}" \
        O="${output_dir}/${sample_id}.bam" \
        TAG="${tag_name}" \
        TAG_VALUE="${tag_value}" \
        FILTER=includeTagValues
done < "${barcode_file}"