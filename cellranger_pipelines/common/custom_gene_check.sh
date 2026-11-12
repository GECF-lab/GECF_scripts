#!/usr/bin/env bash
# =============================================================================
# custom_gene_check.sh — report expression of appended "custom" genes in a
# cellranger count output
#
# Version: 1.0.0
#
# What it does
#   For a cellranger count output directory and a list of custom (transgene /
#   appended) genes, reports — per gene —
#       status          EXPRESSED if the transgene row has >=1 detected cell,
#                       else SILENT (a "detected" cell = UMI count >= threshold;
#                       default 3, change with --min-umi N)
#       tg_cells        cells with UMI >= threshold on the TRANSGENE row only
#       tg_UMI          total UMIs in those detected cells (transgene row)
#       tg_per_cell     tg_UMI / tg_cells (0 if none)
#       endo_cells      cells with UMI >= threshold on ENDOGENOUS rows (other
#                       feature rows that carry the same name, e.g. Ensembl-keyed)
#       endo_UMI        total UMIs in those detected cells (endogenous rows)
#   plus a summary line.
#
# How the split works (verified on AVT0571 / CRREF040)
#   cellranger counts are assigned by ALIGNED COORDINATE, never by name. In a
#   custom genome reference the appended transgene is its own contig; its
#   feature row is the one whose feature ID equals the gene name (e.g. ID
#   "ZNF721", name "ZNF721"), while the endogenous locus is keyed by Ensembl
#   ID (e.g. ID "ENSG00000182903", name "ZNF721"). So:
#       transgene row  = the feature whose ID == the listed gene name
#       endogenous row = every OTHER feature row with that name
#   Reads that map equally well to both loci are DISCARDED by cellranger
#   (multi-mappers to >1 gene are not counted, not split, not misassigned).
#   Consequence: both columns are independent, non-overlapping counts. A
#   SILENT transgene means "no uniquely-assigned reads to the contig" —
#   for transcripts that are sequence-identical copies of the endogenous
#   gene this can also mean the transgene is ON but indistinguishable; the
#   endogenous column is NOT transgene signal either.
#
# FEATURE LIST FORMAT
#   plain text, ONE custom gene name PER LINE. Each line:
#     - first whitespace-separated field is the gene name (rest ignored)
#     - lines starting with '#' and blank lines are ignored
#   The name is matched against the FEATURES table, NOT the cells:
#     - a "transgene row" is the feature whose ID == the listed name
#       (e.g. ID "KZFP-sfGFP", name "KZFP-sfGFP") -> tg_cells/tg_UMI columns
#     - other feature rows whose NAME == the listed name (e.g. Ensembl-ID'd
#       endogenous copies) -> endo_cells/endo_UMI columns
#   So the list file just needs the name you see in features.tsv.gz; if both
#   an Ensembl row and an appended row exist under that name, both columns
#   are populated from the right rows automatically.
#
# Input contract
#   <count_out>     the directory that CONTAINS filtered_feature_bc_matrix/
#                   (i.e. the cellranger 'outs' directory, e.g. NN_sample/outs)
#   <feature_list>  file with custom gene names, ONE PER LINE.
#                   Optional trailing columns (e.g. contig) are ignored.
#
# Detection
#   features:   filtered_feature_bc_matrix/{features.tsv.gz, features.tsv,
#                                   genes.tsv.gz,     genes.tsv}
#   barcodes:   filtered_feature_bc_matrix/{barcodes.tsv.gz, barcodes.tsv}
#   matrix:     filtered_feature_bc_matrix/{matrix.mtx.gz, matrix.mtx}
#
# Detection threshold
#   A cell counts as "detecting" a row when its UMI count for that row is
#   >= the threshold: default 3, change with --min-umi N (integer >= 1,
#   accepted at any CLI position). EXPRESSED = >=1 detected cell on the
#   transgene row.
#
# Usage:
#   bash custom_gene_check.sh --sample-dir  <sample_dir>    <genes> [report_path]
#   bash custom_gene_check.sh --project-dir <project_dir> <genes> [report_path]
#   bash custom_gene_check.sh --list        <sample_dir>    <gene_name>
#
#   <genes> = a gene name (matched exactly) OR a file with one gene name/line.
#
# The report is written to report_path IF given; ALWAYS echoed to stdout.
# Only bash + awk + (z)cat are required (mawk-compatible).
# =============================================================================

# ---- internal helpers ----

# print first existing file among args, or nothing on empty list
_exists() {
    local f
    for f in "$@"; do
        [ -n "$f" ] && [ -f "$f" ] && { printf '%s' "$f"; return 0; }
    done
    return 1
}

# gzip/cat a given path depending on extension
_zread() {
    local f="$1"
    case "$f" in
        *.gz) zcat "$f" 2>/dev/null ;;
        *)    cat  "$f" 2>/dev/null ;;
    esac
}

# default CSV report path for an input folder: <dir>/custom_genes_<dir>.csv
_default_report() {
    local D="${1%/}"
    printf '%s/custom_genes_%s.csv\n' "$D" "$(basename "$D")"
}

# convert a tab-separated report on stdin to RFC-4180 CSV on stdout.
# A field is quoted (embedded double-quotes doubled) only when it contains
# a comma or a double-quote; other fields are written verbatim.
_to_csv() {
    awk '
        BEGIN { FS = "\t" }
        {
            line = ""
            for (i = 1; i <= NF; i++) {
                f = $i
                if (f ~ /[",]/) { gsub(/"/, "\"\"", f); f = "\"" f "\"" }
                line = (i == 1 ? f : line "," f)
            }
            print line
        }
    '
}

# ---- main ----

check_custom_genes() {
    local OUT="${1:-}"
    local LIST="${2:-}"
    local REPORT="${3:-}"

    if [ -z "$OUT" ] || [ ! -d "$OUT" ]; then
        echo "check_custom_genes: <count_out> not a directory: '${OUT}'" >&2
        return 2
    fi
    if [ -z "$LIST" ] || [ ! -f "$LIST" ]; then
        echo "check_custom_genes: <feature_list> not a file: '${LIST}'" >&2
        return 2
    fi

    local FBM="${OUT}/filtered_feature_bc_matrix"
    if [ ! -d "$FBM" ]; then
        echo "check_custom_genes: filtered_feature_bc_matrix/ not found under ${OUT}" >&2
        return 2
    fi

    local FEATURES MATRIX BARCH
    FEATURES="$(_exists "${FBM}/features.tsv.gz" "${FBM}/features.tsv" \
                       "${FBM}/genes.tsv.gz"    "${FBM}/genes.tsv")"
    MATRIX="$(_exists     "${FBM}/matrix.mtx.gz" "${FBM}/matrix.mtx")"
    BARCH="$(_exists      "${FBM}/barcodes.tsv.gz" "${FBM}/barcodes.tsv")"

    if [ -z "$FEATURES" ] || [ -z "$MATRIX" ] || [ -z "$BARCH" ]; then
        echo "check_custom_genes: required files not found in ${FBM}" >&2
        [ -z "$FEATURES" ] && echo "  missing: features.tsv(.gz) or genes.tsv(.gz)" >&2
        [ -z "$MATRIX"   ] && echo "  missing: matrix.mtx(.gz)" >&2
        [ -z "$BARCH"    ] && echo "  missing: barcodes.tsv(.gz)" >&2
        return 2
    fi

    # Denominator for per-gene cell fractions is reported per column (cells
    # with that feature >= 1); TOTAL is kept only as the denominator for
    # pct columns.
    local TOTAL
    TOTAL="$(_zread "$BARCH" | awk 'NF > 0 { n++ } END { print n }')"
    if ! [ "$TOTAL" -gt 0 ] 2>/dev/null; then
        echo "check_custom_genes: could not read a positive total from ${BARCH}" >&2
        return 2
    fi

    # Sample label: sample_name.txt if present, else parent dir name, else basename.
    local SAMPLE
    if [ -f "${OUT}/sample_name.txt" ]; then
        SAMPLE="$(head -n1 "${OUT}/sample_name.txt" | tr -d '\r')"
    else
        SAMPLE="$(basename "$(dirname "$OUT")")"
    fi
    [ -z "$SAMPLE" ] && SAMPLE="$(basename "$OUT")"

    local TMPREPORT
    TMPREPORT="$(mktemp "${TMPDIR:-/tmp}/custom_gene_check.XXXXXXXX")"

    # One streaming pass. Input stream (explicit section separators make the
    # parse deterministic — robust under mawk):
    #   list lines | !!FEAT!! | feature lines | !!MTX!! | raw matrix lines
    # Matrix convention (MatrixMarket): '%' comment lines, then a 3-field
    # dimension line (Nfeat Ncells Nnz) which is the FIRST non-'%' non-blank
    # line, then 'i j v' entries (i indexes features, j indexes cells).
    #
    # Feature classification per listed gene name NM:
    #   TGSET(NM) = feature rows whose ID equals NM        (transgene row)
    #   ENSET(NM) = rows whose NAME equals NM, minus TGSET (endogenous rows)
    local AWK
    AWK='
    function trim(s) { gsub(/^[ \t\r]+/, "", s); gsub(/[ \t\r]+$/, "", s); return s }
    $0 == "!!FEAT!!" { ST = 2; next }
    $0 == "!!MTX!!"  { ST = 3; next }
    ST == 1 {
        # list: first field of each line
        split($0, a, /[ \t]+/)
        nm = trim(a[1])
        if (nm != "") wantorder[++nwant] = nm
        next
    }
    ST == 2 {
        # feature line: <id>\t<name>[\t<type>]; index = 1-based row number
        nfeat++
        nf = split($0, f, "\t")
        idv = trim(f[1]); name = (nf >= 2 ? trim(f[2]) : "")
        if (name != "" && name in namecnt) { namecnt[name] = namecnt[name] " " nfeat }
        else if (name != "") { namecnt[name] = nfeat }
        # transgene row = ID equals the name
        if (idv == name && name != "") {
            if (name in tgid) { tgid[name] = tgid[name] " " nfeat }
            else              { tgid[name] = nfeat }
        }
        next
    }
    ST == 3 {
        if ($0 ~ /^%/ || $0 ~ /^[ \t\r]*$/) next
        if (!dim_seen) { dim_seen = 1; next }   # dimension line: Nfeat Ncells Nnz
        if (NF >= 3 && $3 + 0 >= T) { fi = $1 + 0; cells[fi]++; umi[fi] += $3 }
        next
    }
    BEGIN { ST = 1 }
    END {
        printf "sample\tgene\tstatus\ttg_cells\ttg_UMI\ttg_per_cell\tendo_cells\tendo_UMI\tcells_total\n"
        n_exp = 0; n_silent = 0; n_duo = 0; n_notg = 0
        for (k = 1; k <= nwant; k++) {
            nm = wantorder[k]
            delete tseen; delete eseen
            # --- transgene row indices (id == name) ---
            ttn = 0
            if (nm in tgid) { tt = split(tgid[nm], TA, " "); for (r = 1; r <= tt; r++) { v = TA[r] + 0; if (v > 0 && !(v in tseen)) { tseen[v] = 1; TIX[++ttn] = v } } }
            # --- endogenous row indices (name match, minus transgene rows) ---
            een = 0
            if (nm in namecnt) {
                nn = split(namecnt[nm], NA, " ")
                for (r = 1; r <= nn; r++) {
                    v = NA[r] + 0
                    if (v > 0 && !(v in eseen) && !(v in tseen)) { eseen[v] = 1; EIX[++een] = v }
                }
            }
            tc = 0; tu = 0; ec = 0; eu = 0
            for (q = 1; q <= ttn; q++) { ii = TIX[q]; tc += cells[ii]; tu += umi[ii] }
            for (q = 1; q <= een; q++) { ii = EIX[q]; ec += cells[ii]; eu += umi[ii] }
            if (een > 0) n_duo++
            if (ttn == 0) n_notg++
            tgp = (tc > 0) ? (tu / tc) : 0
            if (tc > 0) { status = "EXPRESSED"; n_exp++ }
            else        { status = "SILENT";    n_silent++ }
            printf "%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\n", sample, nm, status, tc, tu, tgp, ec, eu, TOT
        }
        printf "%s\tCUSTOM SUMMARY: %d custom genes | EXPRESSED=%d | SILENT=%d | with-endogenous-row=%d | no-transgene-row=%d | min_UMI=%d\n", \
               sample, nwant, n_exp, n_silent, n_duo, n_notg, T
    }
    '

    local T="${MIN_UMI:-3}"
    local RC=0
    if [ -t 2 ]; then
        echo "check_custom_genes: scanning '${SAMPLE}' (min_UMI=${T}) -- reads the full filtered matrix;" >&2
        echo "                     large runs take a minute or two, hang tight..." >&2
    fi
    {
        _zread "$LIST" | awk 'NF > 0 && $1 !~ /^#/ { print $1 }'
        printf '!!FEAT!!\n'
        _zread "$FEATURES"
        printf '!!MTX!!\n'
        _zread "$MATRIX"
    } | awk -v sample="$SAMPLE" -v TOT="$TOTAL" -v T="$T" "$AWK" > "$TMPREPORT" 2>/dev/null
    RC=$?
    if [ $RC -ne 0 ]; then
        echo "check_custom_genes: awk failed (rc=$RC)" >&2
        rm -f "$TMPREPORT"
        return $((RC + 1))
    fi

    cat "$TMPREPORT"
    if [ -n "$REPORT" ]; then
        _to_csv < "$TMPREPORT" > "$REPORT"
    fi
    rm -f "$TMPREPORT"
    return 0
}
# ---- discovery/report helpers (self-contained: awk, grep, sed, mktemp only) ----

# _first_feature_file <sample_dir>
#   echo the first existing features/genes file under outs (shared by helpers)
_first_feature_file() {
    local OUTB="${1%/}"
    _exists "${OUTB}/filtered_feature_bc_matrix/features.tsv.gz" \
            "${OUTB}/filtered_feature_bc_matrix/features.tsv" \
            "${OUTB}/filtered_feature_bc_matrix/genes.tsv.gz" \
            "${OUTB}/filtered_feature_bc_matrix/genes.tsv"
}

# list_custom_genes <sample_dir> <gene_name>
#   Show what the run's reference contains for <gene_name>: every row of
#   features.tsv(.gz) mentioning it (case-insensitive), marking the
#   transgene row (feature ID == gene name), and printing the line to put
#   in a feature-list file. Does not run the expression check.
#   rc=0 transgene row exists | rc=3 no transgene row in this reference
list_custom_genes() {
    local OUT="${1%/}"
    local SYM="$2"
    local FEAT
    FEAT="$(_first_feature_file "$OUT")"
    if [ -z "$FEAT" ]; then
        echo "list_custom_genes: no features.tsv(.gz)/genes.tsv(.gz) under ${OUT}" >&2
        return 2
    fi
    echo "Reference feature table: $FEAT"
    echo "Rows mentioning '${SYM}' (case-insensitive):"
    echo "    feature_id                name                     type"
    _zread "$FEAT" | grep -i -- "$SYM" | sed 's/^/    /' || \
        echo "    (no rows mention it)"
    local TGROW
    TGROW="$(_zread "$FEAT" | awk -F'\t' '
        BEGIN { ts = tolower("'"$SYM"'") }
        tolower($1) == ts && tolower($2) == ts { print $1; exit }')"
    if [ -n "$TGROW" ]; then
        echo
        echo "  -> Transgene row (id==name): '${TGROW}'"
        echo "  -> Feature-list line:"
        echo "         $TGROW"
        echo "  (The other row(s) above are endogenous, Ensembl-keyed copies;"
        echo "   the check reports them in the endo_* columns.)"
        return 0
    fi
    echo
    echo "  -> No transgene row (id==name) for '${SYM}' in this run's reference."
    echo "     The transgene is not appended to the genome reference this run"
    echo "     was aligned to, so there is no separate count row for it."
    echo "     A check will report tg_cells=0 for it."
    return 3
}

# check_project_genes <project_dir> <gene_name> [report_path]
#   Run the check for <gene_name> on every outs/ directory found under
#   <project_dir> (recursively) and emit one combined table: header, then one
#   row per sample (sample name in the first column), then a project-level
#   summary line. Per-sample numbers come from check_custom_genes so they are
#   identical to a single-sample run.
#   rc 0 = at least one sample checked | 3 = no outs/ dirs found under project
check_project_genes() {
    local PROJ="${1%/}"
    local GENE="$2"
    local REPORT="${3:-}"

    if [ -z "$GENE" ]; then
        echo "check_project_genes: <gene_name> required" >&2
        return 2
    fi

    local FBM_DIRS=()
    mapfile -t FBM_DIRS < <(find "$PROJ" -type d -name filtered_feature_bc_matrix 2>/dev/null | sort)
    if [ "${#FBM_DIRS[@]}" -eq 0 ]; then
        echo "check_project_genes: no outs/filtered_feature_bc_matrix found under ${PROJ}" >&2
        return 3
    fi

    local LIST
    LIST="$(mktemp "${TMPDIR:-/tmp}/custom_gene_check.XXXXXX")" || return 2
    printf '%s\n' "$GENE" > "$LIST"

    local header="" rows="" RC=0
    local sname n_found=0 n_exp=0 n_silent=0 n_norows=0
    local d FBM OUT out grow st ntl
    local _NS="${#FBM_DIRS[@]}" _i=0

    if [ -t 2 ]; then
        echo "check_project_genes: ${_NS} sample(s) found under ${PROJ} -- scanning each matrix now (can take a while)" >&2
        echo "                     a line prints as each sample finishes; the table below is the final result..." >&2
    fi
    for FBM in "${FBM_DIRS[@]}"; do
        _i=$((_i + 1))
        OUT="$(dirname "$FBM")"
        if [ -t 2 ]; then
            echo "check_project_genes: [${_i}/${_NS}] scanning $(basename "$(dirname "$OUT")") ..." >&2
        fi
        out="$(check_custom_genes "$OUT" "$LIST" 2>/dev/null)" || RC=$?
        if [ -t 2 ]; then
            echo "check_project_genes: [${_i}/${_NS}] $(basename "$(dirname "$OUT")") done" >&2
        fi
        [ -z "$out" ] && continue
        if [ -z "$header" ]; then
            header="$(printf '%s\n' "$out" | sed -n '1p')"
        fi
        grow="$(printf '%s\n' "$out" | awk -F'\t' -v g="$GENE" '$2 == g { print; exit }')"
        [ -z "$grow" ] && continue
        sname="$(printf '%s\n' "$grow" | awk -F'\t' '{ print $1 }')"
        st="$(printf '%s\n' "$grow" | awk -F'\t' '{ print $3 }')"
        ntl="$(printf '%s\n' "$out" | grep -o 'no-transgene-row=[0-9]*' | grep -o '[0-9]*' || true)"
        rows+="${grow}"$'\n'
        n_found=$((n_found + 1))
        case "$st" in
            EXPRESSED) n_exp=$((n_exp + 1)) ;;
            SILENT)    n_silent=$((n_silent + 1)) ;;
        esac
        if [ "${ntl:-0}" -ge 1 ] 2>/dev/null; then n_norows=$((n_norows + 1)); fi
    done

    rm -f "$LIST"

    if [ ${n_found} -eq 0 ]; then
        echo "check_project_genes: no samples produced a result for '${GENE}' under ${PROJ}" >&2
        return $(( ${RC:-3} ))
    fi

    local RESULT
    RESULT="$(
        printf '%s\n' "${header}"
        printf '%s' "${rows}"
        printf '%s\tPROJECT SUMMARY: %d sample(s) | EXPRESSED=%d | SILENT=%d | without-transgene-row=%d | min_UMI=%s\n' \
            "PROJECT" "${n_found}" "${n_exp}" "${n_silent}" "${n_norows}" "${MIN_UMI:-3}" 
    )"

    printf '%s\n' "${RESULT}"
    if [ -n "${REPORT}" ]; then
        printf '%s\n' "${RESULT}" | _to_csv > "${REPORT}"
    fi
    return 0
}

# ---- sample / project / list runners ----

# run_sample_dir <sample_dir> <gene_name_or_file> [report_path]
#   Unified single-sample entry: <gene_name_or_file> is either a single gene
#   name (matched EXACTLY against the reference) or a file with one gene name
#   per line. Both route into check_custom_genes (exact match on the row name).
run_sample_dir() {
    local OUT="${1%/}"
    local TOKEN="${2:-}"
    local REPORT="${3:-}"
    if [ -z "$TOKEN" ]; then
        echo "run_sample_dir: <gene_name> or <gene_list_file> is required" >&2
        return 2
    fi
    if [ -f "$TOKEN" ]; then
        check_custom_genes "$OUT" "$TOKEN" "$REPORT"
        return $?
    fi
    if printf '%s' "$TOKEN" | grep -q '/'; then
        echo "run_sample_dir: '$TOKEN' contains '/' but is not an existing file" >&2
        echo "  - a single gene name must not contain '/'   - a list file must exist" >&2
        echo "  find the exact row name:  --list ${OUT} <gene_name>" >&2
        return 2
    fi
    local TMP
    TMP="$(mktemp "${TMPDIR:-/tmp}/custom_gene_check.XXXXXX")" || return 2
    printf '%s\n' "$TOKEN" > "$TMP"
    local RC=0
    check_custom_genes "$OUT" "$TMP" "$REPORT" || RC=$?
    rm -f "$TMP"
    return "$RC"
}

# run_project_genes <project_dir> <gene_name_or_file> [report_path]
#   Single gene  -> one combined table across every sample in the project
#                   (one row per sample) + PROJECT SUMMARY (reuses
#                   check_project_genes; values identical to single-sample).
#   List file    -> one table per gene (each with its own PROJECT SUMMARY).
#   Writes one combined report to [report_path] when given. rc = last failure.
run_project_genes() {
    local PROJ="${1%/}"
    local TOKEN="${2:-}"
    local REPORT="${3:-}"
    if [ -z "$TOKEN" ]; then
        echo "run_project_genes: <gene_name> or <gene_list_file> is required" >&2
        return 2
    fi
    local GLIST
    GLIST="$(mktemp "${TMPDIR:-/tmp}/custom_gene_check.XXXXXX")" || return 2
    if [ -f "$TOKEN" ]; then
        _zread "$TOKEN" | awk 'NF > 0 && $1 !~ /^#/ { print $1 }' > "$GLIST"
    else
        printf '%s\n' "$TOKEN" > "$GLIST"
    fi
    local NCOUNT
    NCOUNT="$(awk 'NF>0' "$GLIST" | wc -l | awk '{print $1}')"
    if [ "$NCOUNT" -eq 0 ]; then
        echo "run_project_genes: no gene names to check (empty file?)" >&2
        rm -f "$GLIST"
        return 2
    fi
    local COMBINED="" RC=0 IDX=0
    local g HDR BLK
    while IFS= read -r g; do
        [ -n "$g" ] || continue
        IDX=$((IDX+1))
        HDR=""
        if [ "$NCOUNT" -gt 1 ]; then
            HDR=$(printf '\n============================ gene %s (%s/%s) ============================' "$g" "$IDX" "$NCOUNT")
        fi
        if [ -t 2 ]; then
            echo "run_project_genes: checking gene '${g}' (${IDX}/${NCOUNT}) ..." >&2
        fi
        BLK="$(check_project_genes "$PROJ" "$g")" || RC=$?
        if [ -n "$HDR" ]; then
            printf '%s\n' "$HDR"
            COMBINED="${COMBINED}${HDR}"$'\n'
        fi
        if [ -n "$BLK" ]; then
            printf '%s\n' "$BLK"
            COMBINED="${COMBINED}${BLK}"$'\n\n'
        fi
    done < "$GLIST"
    rm -f "$GLIST"
    if [ -z "$COMBINED" ]; then
        echo "run_project_genes: no samples produced a result under ${PROJ}" >&2
        return "$RC"
    fi
    if [ -n "$REPORT" ]; then
        printf '%s\n' "$COMBINED" | _to_csv > "$REPORT"
        echo "CSV report written to: $REPORT" >&2
    fi
    return "$RC"
}

# ---- CLI runner ----

VERSION="1.0.0"
_usage() {
    cat <<'USAGE'
custom_gene_check.sh — report transgene expression in a cellranger count
output: per gene, transgene-row counts (tg_*) vs endogenous-row counts
(endo_*).
Version: 1.0.0

Usage:
  bash custom_gene_check.sh --sample-dir  <sample_dir>    <genes> [report_path]
  bash custom_gene_check.sh --project-dir <project_dir> <genes> [report_path]
  bash custom_gene_check.sh --list        <sample_dir>    <gene_name>
  bash custom_gene_check.sh --version     Print the script version and exit.
  bash custom_gene_check.sh --help        This text.

Global flag (works at any position with --sample-dir / --project-dir):
  --min-umi N   per-cell UMI detection threshold (default 3)

A CSV report is always saved next to the input folder:
  <input_dir>/custom_genes_<input_dir_name>.csv
Override its location with the optional [report_path]. The same report is
also printed to stdout.

<genes> is ONE argument, either a single gene name OR a list file:
  ZNF721          a gene name, matched EXACTLY against the reference
  genes.txt       a file with one gene name per line; each matched EXACTLY
  A reference that uses a prefixed name (e.g. GRCh38_ZNF721) must be passed
  exactly -- use --list to find that exact row name first.

Options:
  --sample-dir   <sample_dir> <genes> [report_path]
        Check <genes> on one sample (one 'outs' dir).
  --project-dir  <project_dir> <genes> [report_path]
        Check <genes> on every sample under <project_dir> (recursively finds
        outs/filtered_feature_bc_matrix/):
          - a single gene  -> one table (one row per sample) + PROJECT SUMMARY
          - a list file    -> one table per gene (each with its own summary)
        Per-sample values are identical to a single-sample run.
  --list         <sample_dir> <gene_name>
        Inspection only. Prints every row of the run's feature table
        (filtered_feature_bc_matrix/features.tsv.gz) mentioning <gene_name>
        and marks the transgene row (feature id == name) -- that exact name is
        the one to pass as <genes>. Does not run the expression check.
  --min-umi N    Global flag (any position) for --sample-dir / --project-dir.
        Minimum UMI count for a cell to count as detecting a row.
        Default: 3. A gene is EXPRESSED when >=1 cell on its transgene row
        reaches N. N must be an integer >= 1.
  --version      Print the script version and exit.
  --help         This text.

Arguments:
  <gene_name>         name of the transgene/gene as present in the reference
                      (e.g. ZNF721); matched exactly
  <genes>             one gene name, OR a file with one gene name per line
                      (first field of each line used; '#' comments and blank
                      lines ignored); names matched exactly
  <sample_dir>          cellranger 'outs' directory (contains
                      filtered_feature_bc_matrix/)
  <project_dir>       folder containing one or more samples (each with its
                      own outs/filtered_feature_bc_matrix/); searched
                      recursively for --project-dir
  [report_path]       path of the CSV report file. Defaults to
                      <input_dir>/custom_genes_<input_dir_name>.csv.
                      The report is always also printed to stdout.

Transgene vs endogenous (report columns):
  tg_*   transgene row: the feature whose id equals the gene name
         (e.g. 'ZNF721  ZNF721'); reads mapped to the appended contig
  endo_* endogenous row(s): same name, Ensembl id (e.g.
         'ENSG00000182903  ZNF721'); the genomic copy
  tg_*=0 with endo_>0: no uniquely-assigned reads to the appended contig
  tg_cells/endo_cells only count cells with UMI >= --min-umi (default 3)

Exit codes:
  0  check ran (per-gene status SILENT/EXPRESSED in the report)
  1  missing required argument
  2  missing/invalid input (files, dirs)
  3  <gene_name> has no transgene row in the run's reference, or the
     project contains no samples
USAGE
}
# ---- optional global flag: --min-umi N (per-cell UMI detection threshold) ----
MIN_UMI="3"
_args=()
while [ $# -gt 0 ]; do
    if [ "$1" = "--min-umi" ]; then
        if [ $# -lt 2 ]; then
            echo "ERROR: --min-umi needs a value (integer >= 1)" >&2
            exit 2
        fi
        MIN_UMI="$2"
        shift
    else
        _args+=("$1")
    fi
    shift
done
set -- "${_args[@]}"
if ! [[ "$MIN_UMI" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --min-umi must be an integer >= 1 (got '${MIN_UMI}')" >&2
    exit 2
fi

case "${1:-}" in
    --version|-v)  echo "${VERSION}" ;;
    --help|-h)     _usage ;;
    --list|-l)
        if [ $# -lt 3 ]; then
            echo "Usage: --list <sample_dir> <gene_name>" >&2
            exit 2
        fi
        list_custom_genes "${2%/}" "$3"
        ;;
    --sample-dir|-s)
        if [ $# -lt 3 ]; then
            echo "Usage: --sample-dir <sample_dir> <genes> [report_path]" >&2
            exit 2
        fi
        rep="${4:-}"
        [ -n "$rep" ] || rep="$(_default_report "${2%/}")"
        run_sample_dir "${2%/}" "$3" "$rep"
        ;;
    --project-dir|-p)
        if [ $# -lt 3 ]; then
            echo "Usage: --project-dir <project_dir> <genes> [report_path]" >&2
            exit 2
        fi
        rep="${4:-}"
        [ -n "$rep" ] || rep="$(_default_report "${2%/}")"
        run_project_genes "${2%/}" "$3" "$rep"
        ;;
    *)
        if [ $# -lt 2 ]; then
            _usage >&2
            echo "" >&2
            echo "ERROR: <mode> and <genes> are required." >&2
            echo "  sample mode :  --sample-dir  <sample_dir> <genes> [report]" >&2
            echo "  project mode:  --project-dir <project>  <genes> [report]" >&2
            exit 1
        fi
        echo "ERROR: unknown mode '${1}' (expected --sample-dir, --project-dir, or --list)." >&2
        echo "  Run 'custom_gene_check.sh --help' for details." >&2
        exit 2
        ;;
esac
