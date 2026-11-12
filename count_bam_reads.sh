#!/usr/bin/env bash
# Usage: ./count_reads.sh /path/to/folder
set -euo pipefail

find "$1" -type f -name '*.bam' -print0 |
while IFS= read -r -d '' bam; do
    printf '%s: %s\n' "$bam" "$(samtools view -c "$bam")"
done
