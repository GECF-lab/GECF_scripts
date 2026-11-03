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
#       status          EXPRESSED if the transgene row has >=1 expressed cell,
#                       else SILENT
#       tg_cells        cells with >=1 count on the TRANSGENE feature row only
#       tg_UMI          total UMIs on the TRANSGENE feature row only
#       tg_per_cell     tg_UMI / tg_cells (0 if none)
#       endo_cells      cells with >=1 count on ENDOGENOUS rows (other feature
#                       rows that carry the same name, e.g. Ensembl-keyed)
#       endo_UMI        total UMIs on those endogenous rows
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
# Library usage (embedded into a generated pipeline .sh)
#   The generator inlines everything ABOVE the `# %%EMBED-END%%` marker
#   verbatim, then calls:
#     check_custom_genes <count_out> <feature_list> [out_report_file]
#
# Standalone usage:
#   bash custom_gene_check.sh <count_out> <feature_list> [report_path]
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
        if (NF >= 3) { fi = $1 + 0; cells[fi]++; umi[fi] += $3 }
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
        printf "%s\tCUSTOM SUMMARY: %d custom genes | EXPRESSED=%d | SILENT=%d | with-endogenous-row=%d | no-transgene-row=%d\n", \
               sample, nwant, n_exp, n_silent, n_duo, n_notg
    }
    '

    local RC=0
    {
        _zread "$LIST" | awk 'NF > 0 { sub(/^#/, ""); print }' | awk 'NF > 0 { print $1 }'
        printf '!!FEAT!!\n'
        _zread "$FEATURES"
        printf '!!MTX!!\n'
        _zread "$MATRIX"
    } | awk -v sample="$SAMPLE" -v TOT="$TOTAL" "$AWK" > "$TMPREPORT" 2>/dev/null
    RC=$?
    if [ $RC -ne 0 ]; then
        echo "check_custom_genes: awk failed (rc=$RC)" >&2
        rm -f "$TMPREPORT"
        return $((RC + 1))
    fi

    cat "$TMPREPORT"
    if [ -n "$REPORT" ]; then
        cp "$TMPREPORT" "$REPORT"
    fi
    rm -f "$TMPREPORT"
    return 0
}
# ---- discovery/report helpers (embedded with the script; self-contained:
#      awk, grep, sed, mktemp only — safe to call from generated .sh files) ----

# _first_feature_file <outs_dir>
#   echo the first existing features/genes file under outs (shared by helpers)
_first_feature_file() {
    local OUTB="${1%/}"
    _exists "${OUTB}/filtered_feature_bc_matrix/features.tsv.gz" \
            "${OUTB}/filtered_feature_bc_matrix/features.tsv" \
            "${OUTB}/filtered_feature_bc_matrix/genes.tsv.gz" \
            "${OUTB}/filtered_feature_bc_matrix/genes.tsv"
}

# list_custom_genes <outs_dir> <gene_name>
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

# detect_custom_gene <outs_dir> <gene_name> [report_path]
#   One-shot check with no list file: find the transgene row for <gene_name>
#   in the reference (id==name), write it to a temp list file, run the
#   normal check, print the result table to stdout, and copy it to
#   [report_path] if given.
#   rc same as check_custom_genes (0 ok) | 3 no transgene row in reference
detect_custom_gene() {
    local OUT="${1%/}"
    local SYM="$2"
    local REPORT="${3:-}"
    local FEAT
    FEAT="$(_first_feature_file "$OUT")"
    if [ -z "$FEAT" ]; then
        echo "detect_custom_gene: no features.tsv(.gz)/genes.tsv(.gz) under ${OUT}" >&2
        return 2
    fi
    local TMP
    TMP="$(mktemp "${TMPDIR:-/tmp}/custom_gene_check.XXXXXX")" || return 2
    # transgene rows: id == name and row matches symbol (case-insensitive)
    _zread "$FEAT" | awk -F'\t' -v s="$SYM" '
        BEGIN { ts = tolower(s) }
        {
            id = tolower($1); nm = tolower($2)
            if (id == nm && (id == ts || index(nm, ts) > 0)) print $1
        }' > "$TMP"
    if [ ! -s "$TMP" ]; then
        echo "detect_custom_gene: no transgene feature row (id==name) matched '${SYM}' in $FEAT" >&2
        echo "  show all matching rows:  --list ${SYM} ${OUT}" >&2
        rm -f "$TMP"
        return 3
    fi
    echo "auto-selected feature row(s) from reference:"
    sed 's/^/  /' "$TMP"
    echo "----------------------------------------------"
    local RC=0
    check_custom_genes "$OUT" "$TMP" "$REPORT" || RC=$?
    rm -f "$TMP"
    return $RC
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

    for FBM in "${FBM_DIRS[@]}"; do
        OUT="$(dirname "$FBM")"
        out="$(check_custom_genes "$OUT" "$LIST" 2>/dev/null)" || RC=$?
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
        printf '%s\tPROJECT SUMMARY: %d sample(s) | EXPRESSED=%d | SILENT=%d | without-transgene-row=%d\n' \
            "PROJECT" "${n_found}" "${n_exp}" "${n_silent}" "${n_norows}"
    )"

    printf '%s\n' "${RESULT}"
    if [ -n "${REPORT}" ]; then
        printf '%s\n' "${RESULT}" > "${REPORT}"
    fi
    return 0
}

# make_list <out_file> <symbol> [symbol ...]
#   Convenience: write a ready-to-use feature-list file (one symbol per line)
#   so subsequent `--detect` or direct runs have consistent input.
make_list_file() {
    local OUTFILE="${1:-}"
    shift
    if [ -z "$OUTFILE" ] || [ $# -eq 0 ]; then
        echo "make_list_file: usage: --make-list <out_file> <symbol> [symbol ...]" >&2
        return 2
    fi
    : > "$OUTFILE"
    local s
    for s in "$@"; do
        printf '%s\n' "$s" >> "$OUTFILE"
    done
    echo "wrote feature list: $OUTFILE ($(wc -l < "$OUTFILE") lines)"
    return 0
}
# ---- the embedded region ends here (generator inlines everything above this) ----
# %%EMBED-END%%

# ---- standalone CLI runner (NOT embedded; only runs when executed directly) ----
if [ "${BASH_SOURCE[0]:-}" = "${0}" ]; then
    VERSION="1.0.0"
    _usage() {
        cat <<'USAGE'
custom_gene_check.sh — report transgene expression in a cellranger count
output: per gene, transgene-row counts (tg_*) vs endogenous-row counts
(endo_*).
Version: 1.0.0

Usage:
  bash custom_gene_check.sh --detect  <gene_name> <outs_dir> [report_path]
  bash custom_gene_check.sh --list    <gene_name> <outs_dir>
  bash custom_gene_check.sh --project <project_dir> <gene_name> [report_path]
  bash custom_gene_check.sh --make-list <out_file> <gene_name> [gene_name ...]
  bash custom_gene_check.sh <outs_dir> <feature_list_file> [report_path]
  bash custom_gene_check.sh --version
  bash custom_gene_check.sh --help

Options:
  --detect <gene_name> <outs_dir> [report_path]
      One-shot check on one sample (one 'outs' dir).
  --list <gene_name> <outs_dir>
      Inspection only. Prints every row of the run's feature table
      (filtered_feature_bc_matrix/features.tsv.gz) mentioning <gene_name>
      and prints the transgene row (feature id == name), which is the line
      to put in a feature-list file. Does not run the expression check.
  --project <project_dir> <gene_name> [report_path]
      Check one <gene_name> across all samples in a project: finds every
      outs/filtered_feature_bc_matrix under <project_dir> (recursively),
      runs the check on each sample, and prints one combined table
      (one row per sample) plus a PROJECT SUMMARY line (EXPRESSED / SILENT
      counts and how many samples lack the transgene row in their
      reference). Per-sample values are identical to a single-sample run.
  --make-list <out_file> <gene_name> ...
      Writes a feature-list file, one gene name per line.
  --version
      Print the script version and exit.
  --help
      This text.

  For --detect and --list, <gene_name> and <outs_dir> are accepted in
  either order.

Arguments:
  <gene_name>         name of the transgene/gene as present in the reference
                      (e.g. ZNF721)
  <outs_dir>          cellranger 'outs' directory (contains
                      filtered_feature_bc_matrix/)
  <project_dir>       folder containing one or more samples (each with its
                      own outs/filtered_feature_bc_matrix/); searched
                      recursively for --project
  <feature_list_file> plain text, one gene name per line; first field of
                      each line used; '#' comments and blank lines ignored
  [report_path]       file to also write the report to; table is always
                      printed to stdout

Transgene vs endogenous (report columns):
  tg_*   transgene row: the feature whose id equals the gene name
         (e.g. 'ZNF721  ZNF721'); reads mapped to the appended contig
  endo_* endogenous row(s): same name, Ensembl id (e.g.
         'ENSG00000182903  ZNF721'); the genomic copy
  tg_*=0 with endo_>0: no uniquely-assigned reads to the appended contig

Exit codes:
  0  check ran (per-gene status SILENT/EXPRESSED in the report)
  1  missing required argument
  2  missing/invalid input (files, dirs)
  3  <gene_name> has no transgene row in the run's reference, or the
     project contains no samples
USAGE
    }
    case "${1:-}" in
        --version|-v)  echo "${VERSION}" ;;
        --help|-h)  _usage ;;
        --list|-l)
            [ $# -ge 3 ] || { echo "Usage: --list <gene_name> <outs_dir>" >&2; exit 2; }
            if [ -d "${2%/}" ]; then list_custom_genes "$2" "$3"
            elif [ -d "${3%/}" ]; then list_custom_genes "$3" "$2"
            else list_custom_genes "$2" "$3"; fi
            ;;
        --detect|-d)
            [ $# -ge 3 ] || { echo "Usage: --detect <gene_name> <outs_dir> [report_path]" >&2; exit 2; }
            if [ -d "${2%/}" ]; then detect_custom_gene "$2" "$3" "${4:-}"
            elif [ -d "${3%/}" ]; then detect_custom_gene "$3" "$2" "${4:-}"
            else detect_custom_gene "$2" "$3" "${4:-}"; fi
            ;;
        --project|-p)
            [ $# -ge 3 ] || { echo "Usage: --project <project_dir> <gene_name> [report_path]" >&2; exit 2; }
            if [ -d "${2%/}" ]; then check_project_genes "$2" "$3" "${4:-}"
            elif [ -d "${3%/}" ]; then check_project_genes "$3" "$2" "${4:-}"
            else check_project_genes "$2" "$3" "${4:-}"; fi
            ;;
        --make-list|-m)
            [ $# -ge 3 ] || { echo "Usage: --make-list <out_file> <symbol> [symbol ...]" >&2; exit 2; }
            make_list_file "$2" "${@:3}"
            ;;
        *)
            if [ $# -lt 2 ]; then
                _usage
                echo "" >&2
                echo "ERROR: <feature_list_file> missing — one custom gene name per line." >&2
                echo "  build one with: --make-list <out_file> <gene_name> ..." >&2
                echo "  or run without a list file: --detect <gene_name> <outs_dir>" >&2
                exit 1
            fi
            check_custom_genes "$@"
            ;;
    esac
fi
