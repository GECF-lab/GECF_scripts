#!/bin/bash
# TSO_percent.sh — GEX Read-2 library-quality metric (TSO / SMART / polyA / polyT)
#
# On a 5,000,000-read downsample of the Read-2 sequence lines, reports the
# fraction of reads that:
#   intact TSO   : contain  AAGCAGTGGTATCAACGCAGAGTACAT
#   intact SMART : contain  GGTATCAACGCAGAGTACTTTTTTTT
#   polyA tail   : START WITH >=5 As
#   polyT tail   : START WITH >=5 Ts
#
# Standalone usage:
#   bash TSO_percent.sh <read2.fastq.gz> [sample_name] [out_file]
#       <read2.fastq.gz>  the Read-2 fastq (zcat-able .fastq.gz)
#       [sample_name]     optional prefix on the report header
#       [out_file]        if given, the report is also written there
#
# Library usage: the GEX generator embeds the region above the
# "# %%EMBED-END%%" marker verbatim into the generated .sh and then calls:
#   tso_percent <read2.fastq.gz> [sample_name] [out_file]
# =============================================================================

tso_percent() {
    local R2="${1:-}"
    local SAMPLE="${2:-}"
    local OUT="${3:-}"

    if [ -z "${R2}" ] || [ ! -r "${R2}" ]; then
        echo "tso_percent: no readable R2 fastq given (got: '${R2}')" >&2
        return 2
    fi

    local DOWNSAMPLING=5000000
    local TSO="AAGCAGTGGTATCAACGCAGAGTACAT"
    local SMART="GGTATCAACGCAGAGTACTTTTTTTT"
    local POLYA="AAAAA"
    local POLYT="TTTTT"

    # Sequence line = the 2nd line of each 4-line record. Read the first
    # DOWNSAMPLING sequence lines ONCE into a temp file (avoids 4 passes over
    # the .gz and avoids stashing gigabytes of text in a shell variable).
    local TMPLINE
    TMPLINE="$(mktemp "${TMPDIR:-/tmp}/tso_percent.XXXXXXXX")"
    zcat "${R2}" 2>/dev/null | sed -n '2~4p' | head -n "${DOWNSAMPLING}" > "${TMPLINE}" || true

    local N_TSO N_SMART N_POLYA N_POLYT TOTAL D
    N_TSO="$(grep -c "${TSO}"    "${TMPLINE}" || true)"
    N_SMART="$(grep -c "${SMART}" "${TMPLINE}" || true)"
    N_POLYA="$(grep -c "^${POLYA}" "${TMPLINE}" || true)"
    N_POLYT="$(grep -c "^${POLYT}" "${TMPLINE}" || true)"
    TOTAL="$(grep -c '.'  "${TMPLINE}" || true)"
    D="${DOWNSAMPLING}"
    [ "${TOTAL}" -gt 0 ] && D="${TOTAL}"

    _pct() { awk -v a="$1" -v b="${D}" 'BEGIN{ if (b>0) printf "%.2f",(a*100)/b; else printf "NA" }'; }
    local P_TSO P_SMART P_POLYA P_POLYT
    P_TSO="$(_pct "${N_TSO}")";   P_SMART="$(_pct "${N_SMART}")"
    P_POLYA="$(_pct "${N_POLYA}")"; P_POLYT="$(_pct "${N_POLYT}")"

    local REPORT
    REPORT="${SAMPLE:-read2}  (read count sampled: ${D})
intact TSO   : ${P_TSO} %
intact SMART : ${P_SMART} %
polyA tail   : ${P_POLYA} %
polyT tail   : ${P_POLYT} %
TSO=${P_TSO} SMART=${P_SMART} POLYA=${P_POLYA} POLYT=${P_POLYT} n=${D}"

    printf '%s\n' "${REPORT}"
    [ -n "${OUT}" ] && printf '%s\n' "${REPORT}" > "${OUT}"

    rm -f "${TMPLINE}"
    return 0
}
# %%EMBED-END%%

# ---- standalone runner (NOT part of the embedded region above) ----
if [ "${BASH_SOURCE[0]:-}" = "${0}" ]; then
    if [ $# -lt 1 ] || [ ! -e "${1}" ]; then
        echo "Usage: bash TSO_percent.sh <read2.fastq.gz> [sample_name] [out_file]" >&2
        exit 1
    fi
    tso_percent "$@"
fi
