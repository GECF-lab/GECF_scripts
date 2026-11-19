#!/usr/bin/env bash
# Usage: ./count_reads.sh /path/to/folder
set -euo pipefail

find "$1" -type f -name '*.bam' -print0 |
while IFS= read -r -d '' bam; do
	bam_name=$(basename $bam)
	printf '%s: %s\n' "$bam_name" "$(samtools view -c "$bam")"
done
