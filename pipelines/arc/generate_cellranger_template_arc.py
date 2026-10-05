#!/usr/bin/env python3
"""
generate_cellranger_template_arc.py

Usage:
    python3 generate_cellranger_template_arc.py file.xlsx [more.xlsx ...]
    python3 generate_cellranger_template_arc.py -n file.xlsx

Reads one or more AVT submission Excel files. From each workbook it
extracts the ``Sample name`` and ``GECF fastq ID (avidxxxx)`` columns,
a run name matching ``AVT0XXX`` (any cell or the file name), and the
requester (first word of the ``User:`` cell, if present).

Libraries are grouped by biologic sample: the base name is the sample
name with a trailing ``_GEX`` / ``_ATAC`` suffix removed, and each
library is classified Gene Expression or Chromatin Accessibility.

For every biologic sample the generated script writes one ``--libraries``
CSV (the same 3-column layout cellranger-arc consumes) and one
``cellranger-arc count`` invocation:

    fastqs,sample,library_type
    <fastq_leaf_dir>,<fastq_id>,Gene Expression
    <fastq_leaf_dir>,<fastq_id>,Chromatin Accessibility

    cellranger-arc count \\
        --id=20260310_mAP_CD28_4_ED --description=avid12500,avid12499 \\
        --reference=$REF --libraries=libraries_20260310_mAP_CD28_4_ED.csv \\
        --create-bam=true --localcores=44 --localmem=150 \\
        --cell-annotation-model=auto \\
        --tenx-cloud-token-path=/home/gecf/.config/txg/credentials

Fastq leaf directories are resolved at generation time under
``<run_dir>/fastq_bases2fastq/<batch>/Samples/<experiment>/<fastq_id>/``;
when they cannot be resolved (raw data not mounted here) a clearly marked
placeholder is emitted for the operator to fill in.
"""

import argparse
import os
import re
import sys
import warnings
from collections import OrderedDict

warnings.filterwarnings("ignore",
                        message="Data Validation extension is not supported")
warnings.filterwarnings("ignore",
                        message="Conditional Formatting extension is not supported")

try:
    import openpyxl
except ImportError:
    sys.exit("Error: the python module 'openpyxl' is required.\n"
             "Install it with: pip3 install openpyxl")

# Ensure the pipelines/ package root is on sys.path so ``common.intake``
# resolves.  This file lives at ``pipelines/arc/`` so its parent is the
# package root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PIPELINES_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir))
if _PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, _PIPELINES_ROOT)

# Shared intake atoms — see pipelines/common/intake.py.  All the low-level
# workbook helpers (norm, merged-cell lookup, cell read, the 'User:' label
# scan) are imported from there and used by both this generator and the
# ATAC generator on an identical code path.
from common.intake import (  # noqa: E402
    NGS_BASE_DIR_TEMPLATE,
    AVT_PATTERN,
    norm,
    build_merged_cell_lookup,
    get_cell_value,
    read_user_field,
    set_non_interactive,
    setup_readline_completer,
    ask_avt_run,
    ask_reference_genome,
)

NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/arc"

HUMAN_GENOME_PATH = ("/home/gecf/references/cellranger_references/"
                     "current_human_arc_ATAC/"
                     "refdata-cellranger-arc-GRCh38-2024-A")
MOUSE_GENOME_PATH = ("/home/gecf/references/cellranger_references/"
                     "current_mouse_arc_ATAC/"
                     "refdata-cellranger-arc-GRCm39-2024-A")

# (AVT_PATTERN is imported from common.intake — see the import block above.)

# Batch / control folders that are not real sample leaves.
EXCLUDED_BATCH_PATTERNS = (
    re.compile(r"phix", re.IGNORECASE),
    re.compile(r"\bunassigned\b", re.IGNORECASE),
    re.compile(r"for_check", re.IGNORECASE),
)

LOCAL_CORES = "44"
LOCAL_MEM = "150"
TENX_TOKEN_PATH = "/home/gecf/.config/txg/credentials"


# ─── input helpers ────────────────────────────────────────────────────


setup_readline_completer()  # tab-complete paths in the input() prompts


# (norm / get_cell_value / build_merged_cell_lookup are imported from
# common.intake — see the import block above.)


# A real GECF library id is ``avid`` followed by 4-6 digits (e.g. avid12499).
# Label/legend/placeholder rows in the same columns carry strings like
# "0", "GECF libr. ID (avidxxxx)" or "Loaded molarity (pM)" — this pattern
# distinguishes genuine sample rows from those.
FASTQ_ID_PATTERN = re.compile(r"^avid\d{4,6}$", re.IGNORECASE)


def _locate_header(ws, merged_lookup):
    """Return (fastq_col, sample_col, header_row) or (None, None, None)."""
    fastq_col = sample_col = header_row = None
    for row_idx in range(1, min(25, ws.max_row) + 1):
        for col_idx in range(1, ws.max_column + 1):
            v = norm(get_cell_value(ws, merged_lookup, row_idx, col_idx))
            if "fastq id" in v and "avid" in v and fastq_col is None:
                fastq_col = col_idx
            if "sample name" in v and sample_col is None:
                sample_col = col_idx
        if fastq_col is not None and sample_col is not None:
            header_row = row_idx
            break
    return fastq_col, sample_col, header_row


def _sheet_sample_rows(ws, merged_lookup):
    """Extract (sample, fastq) rows from one sheet, or [] if it is not a data sheet."""
    fastq_col, sample_col, header_row = _locate_header(ws, merged_lookup)
    if fastq_col is None or sample_col is None:
        return []

    data_start = header_row + 1
    for rng in ws.merged_cells.ranges:
        if (rng.min_row <= header_row <= rng.max_row
                and ((rng.min_col <= fastq_col <= rng.max_col)
                     or (rng.min_col <= sample_col <= rng.max_col))):
            data_start = max(data_start, rng.max_row + 1)

    rows = []
    blank_streak = 0
    for row_idx in range(data_start, ws.max_row + 1):
        sample_val = get_cell_value(ws, merged_lookup, row_idx, sample_col)
        fastq_val = get_cell_value(ws, merged_lookup, row_idx, fastq_col)
        sample_txt = str(sample_val).strip() if sample_val is not None else ""
        fastq_txt = str(fastq_val).strip() if fastq_val is not None else ""
        sample_txt = "" if sample_txt == "0" else sample_txt
        fastq_txt = "" if fastq_txt == "0" else fastq_txt

        if "run details" in norm(sample_txt) or "run details" in norm(fastq_txt):
            break

        if sample_txt and fastq_txt and FASTQ_ID_PATTERN.match(fastq_txt):
            rows.append((sample_txt, fastq_txt))
            blank_streak = 0
        else:
            blank_streak += 1
            if rows and blank_streak >= 2:
                break
    return rows


def best_data_sheet(wb):
    """
    Pick the sheet that carries the real sample table, and return
    (sheet_name, rows).

    Preference order:
      1. the ``User_x1x_copy`` sheet — the canonical single-user sample
         block of the GECF submission template;
      2. any other ``User_xNx_copy`` sheet, lowest N first;
      3. fallback: whichever of the remaining sheets (e.g. the
         multi-user ``whole_severaluser`` summary) yields the most valid
         rows. ``Spike_*`` sheets are listed for context only and ranked
         last, since spike controls belong to a different pipeline input.
    """
    parsed = []
    for name in wb.sheetnames:
        ws = wb[name]
        if ws.max_row is None or ws.max_row < 2:
            continue
        try:
            rows = _sheet_sample_rows(ws, build_merged_cell_lookup(ws))
        except Exception:
            rows = []
        if rows:
            parsed.append((name, rows))

    def score(name):
        low = name.lower()
        # Sheet names look like "User_x1x_copy" / "User_x2x_copy".
        m1 = re.match(r"^user_x(\d+)x?_copy$", low)
        if m1:
            tier = 0 if m1.group(1) == "1" else 1
            return (tier, int(m1.group(1)), low)
        if low.startswith("whole_severaluser"):
            return (2, 0, low)
        if low.startswith("spike"):
            return (4, 0, low)
        return (3, 0, low)

    parsed.sort(key=lambda c: (score(c[0]), -len(c[1])))
    if not parsed:
        return None, []
    return parsed[0]


def parse_xlsx(path):
    """
    Parse one AVT submission workbook.

    Returns (rows, detected_avt) where ``rows`` is a list of
    ``(sample_name, fastq_id)`` tuples (real library ids only) and
    ``detected_avt`` is the AVT run name found in a cell or the file
    name (or None).

    Handles both the single-sheet layout of older submissions
    (e.g. AVT0519) and the 15-sheet GECF submission workbooks where the
    sample table is NOT the first sheet (first sheet is GUIDELINES).
    """
    wb = openpyxl.load_workbook(path, data_only=True)

    sheet_name, rows = best_data_sheet(wb)
    if rows:
        # De-duplicate, preserving first-seen order.
        seen = set()
        deduped = []
        for sample, fastq in rows:
            key = (sample.lower(), fastq.lower())
            if key in seen:
                continue
            seen.add(key)
            deduped.append((sample, fastq))
        rows = deduped

    # Detect the AVT run name from any data row or, failing that, the file name.
    detected_avt = None
    for sample, fastq in rows:
        m = AVT_PATTERN.search(sample) or AVT_PATTERN.search(fastq)
        if m:
            detected_avt = m.group(0).upper()
            break
    if detected_avt is None:
        m = AVT_PATTERN.search(os.path.basename(path))
        if m:
            detected_avt = m.group(0).upper()

    return rows, detected_avt


def latest_avt(avt_names):
    """Return the detected AVT run with the largest numeric suffix."""
    return max(
        avt_names,
        key=lambda name: int(re.search(r"\d+", name).group()),
        default=None,
    )

# (read_user_field is imported from common.intake — see the import block above.)


# ─── grouping ─────────────────────────────────────────────────────────


def base_sample_name(name):
    """Strip a trailing ``_GEX`` / ``_ATAC`` (any case) from a sample name."""
    return re.sub(r"_(GEX|ATAC)$", "", name, flags=re.IGNORECASE)


def modality_of(full_name):
    if re.search(r"_GEX$", full_name, re.IGNORECASE):
        return "Gene Expression"
    if re.search(r"_ATAC$", full_name, re.IGNORECASE):
        return "Chromatin Accessibility"
    return "Gene Expression"


def group_libraries(rows):
    """
    Group libraries by biologic sample.

    Returns an OrderedDict: base_name -> list of
    ``(full_name, fastq_id, modality)`` in original order.
    """
    groups = OrderedDict()
    for sample, fid in rows:
        full = sample.strip()
        if not full:
            continue
        base = base_sample_name(full)
        groups.setdefault(base, []).append(
            (full, (fid or "").strip(), modality_of(full)))
    return groups


# ─── fastq leaf resolution ────────────────────────────────────────────


def is_excluded_batch(folder_name):
    return any(p.search(folder_name) for p in EXCLUDED_BATCH_PATTERNS)


def resolve_fastq_leaf(run_dir, fastq_id):
    """
    Resolve the fastq leaf directory for one fastq ID under the run
    folder's ``fastq_bases2fastq/<batch>/Samples`` structure.

    Returns the absolute leaf path, or None if not found.
    """
    if not fastq_id:
        return None
    root = os.path.join(run_dir, "fastq_bases2fastq")
    if not os.path.isdir(root):
        return None
    for batch in sorted(os.listdir(root)):
        if not batch or batch.startswith(".") or is_excluded_batch(batch):
            continue
        samples_dir = os.path.join(root, batch, "Samples")
        if not os.path.isdir(samples_dir):
            continue
        for dirpath, dirnames, _ in os.walk(samples_dir):
            if fastq_id in dirnames:
                leaf = os.path.join(dirpath, fastq_id)
                if os.path.isdir(leaf):
                    return leaf
    return None


# ─── prompts ──────────────────────────────────────────────────────────


# (ask_avt_run / ask_reference_genome come from common/intake — see the
# import block near the top of this file.)

# ─── display ──────────────────────────────────────────────────────────


def print_groups_table(groups):
    print("Biologic samples (one cellranger-arc count each):")
    print(f"{'#':<3} {'Sample (id)':<34} {'GEX':>4} {'ATAC':>5}  fastq IDs")
    for i, (base, libs) in enumerate(groups.items(), start=1):
        n_gex = sum(1 for _, _, m in libs if m == "Gene Expression")
        n_atac = sum(1 for _, _, m in libs if m == "Chromatin Accessibility")
        ids = ", ".join(fid for _, fid, _ in libs if fid) or "-"
        print(f"{i:<3} {base:<34} {n_gex:>4} {n_atac:>5}  {ids}")


# ─── shell script emission ────────────────────────────────────────────

# Library-type ordering: Gene Expression rows are emitted before Chromatin
# Accessibility rows so the generated libraries CSV and the --description
# value follow the arc convention.
_GEX, _ATAC = "Gene Expression", "Chromatin Accessibility"


def build_script(reference_genome, groups, avt_run, user, run_dir):
    L = []
    L.append("#!/bin/bash")
    L.append("set -euo pipefail")
    L.append("")
    L.append("# cellranger-arc — one `cellranger-arc count` per biologic sample.")
    L.append("# Each sample carries a Gene Expression and a Chromatin")
    L.append("# Accessibility library, resolved from the run's")
    L.append("# `fastq_bases2fastq/<batch>/Samples/<experiment>/<avid>/` leaves.")
    L.append(f"# Requester: {user or 'n/a'}")
    L.append("")

    # Fixed paths.
    L.append(f"avt_run='{avt_run}'")
    L.append(f'REFPATH="{reference_genome}"')
    L.append(f'NGSruns_folder="{NGSruns_ARCHIVE_DIR}"')
    # Local working folder where cellranger-arc writes its run dirs +
    # libraries CSVs (separate from the raw fastq source folder above).
    # Namespaced per run, same as CRarc_patched_full:
    # /mnt/data/CellRangerCountOutput/<AVT>_<user>_<timestamp>/
    L.append('scratch_root="/mnt/data/CellRangerCountOutput"')
    L.append("timestamp=$(date +%Y%m%d-%H%M%S)")
    L.append(f"user='{user or ''}'")
    L.append('run_subdir="${avt_run}${user:+_${user}}_${timestamp}"')
    L.append('Server_folder="${scratch_root}/${run_subdir}"')
    L.append('metrics_folder="${Server_folder}/raw_metrics_folder"')
    L.append('archive_dir="${NGSruns_folder}/${run_subdir}"')
    L.append('# set keepdata="yes" to skip local cleanup (CS tree + workdir) at the end')
    L.append('keepdata="no"')
    L.append("run_ids=(" + " ".join(f"'{b}'" for b in groups) + ")")
    L.append("")
    L.append('mkdir -p "${metrics_folder}" \\')
    L.append('    "${metrics_folder}/GEX_only" "${metrics_folder}/ATAC_only" \\')
    L.append('    "${archive_dir}"')
    L.append('cd "${Server_folder}"')
    L.append("")

    # One arc invocation per biologic sample.
    for idx, (base, libs) in enumerate(groups.items(), start=1):
        ordered = sorted(
            libs, key=lambda t: 0 if t[2] == _GEX else 1
        )
        fids = [fid for _, fid, _ in ordered if fid]
        desc = ",".join(fids) if fids else base
        n_gex = sum(1 for _, _, m in libs if m == _GEX)
        n_atac = sum(1 for _, _, m in libs if m == _ATAC)
        csv_name = f"libraries_{base}.csv"

        L.append(f"# ══ Sample {idx}: {base}  ({n_gex} GEX + {n_atac} ATAC) ══")
        L.append("")
        L.append(f'csv_path="${{Server_folder}}/{csv_name}"')
        L.append('echo "fastqs,sample,library_type" > "${csv_path}"')

        for full, fid, modality in ordered:
            leaf = resolve_fastq_leaf(run_dir, fid) if fid else None
            if leaf:
                L.append(
                    f'echo "{leaf},{full},{modality}" >> "${{csv_path}}"'
                )
            else:
                L.append(
                    f"# NOTE: fastq leaf for '{full}' "
                    f"(fid {fid or '-'}) was not resolved here; "
                    f"fill in the directory below."
                )
                L.append(
                    f'echo "/PATH/TO/Samples/.../{fid or full},{full},{modality}" '
                    f'>> "${{csv_path}}"'
                )

        L.append(f'cmd_log="${{Server_folder}}/{base}_script_log.txt"')
        # CRarc_patched_full (l.558) logs the exact command on ONE line. The
        # emitted `echo "…"` body must therefore carry the logical newlines of
        # the cellranger invocation as `\\` (escaped backslash) INSIDE the
        # double-quoted shell string — bash then prints a literal `\` at each
        # break, giving a single-line, readable log entry.
        log_line = (
            [
                "cellranger-arc count",
                f"--id={base}",
                f"--description={desc}",
                '--reference="$REFPATH"',
                f"--libraries={csv_name}",
                "--create-bam=true",
                f"--localcores={LOCAL_CORES}",
                f"--localmem={LOCAL_MEM}",
                "--cell-annotation-model=auto",
                f"--tenx-cloud-token-path={TENX_TOKEN_PATH}",
            ]
        )
        log_body = " \\\\\n    ".join(log_line)
        L.append(f'echo "{log_body}" > "${{cmd_log}}"')
        L.append("cellranger-arc count \\")
        L.append(f"    --id={base} \\")
        L.append(f"    --description={desc} \\")
        L.append('    --reference="$REFPATH" \\')
        L.append(f"    --libraries={csv_name} \\")
        L.append("    --create-bam=true \\")
        L.append(f"    --localcores={LOCAL_CORES} \\")
        L.append(f"    --localmem={LOCAL_MEM} \\")
        L.append("    --cell-annotation-model=auto \\")
        L.append(f"    --tenx-cloud-token-path={TENX_TOKEN_PATH}")
        # CRarc_patched_full (l.567-573): after the arc call, probe for
        # leftover `tmp` folders as an incomplete-count signal and persist
        # the result to `<sample>_error_log.txt` (empty if all good).
        L.append(
            f'error_log="${{Server_folder}}/{base}_error_log.txt"'
        )
        L.append(
            'find "${Server_folder}/' + base
            + '" -mindepth 2 -maxdepth 2 -type d -name "tmp" -print '
            '2>/dev/null > "${error_log}"'
        )
        L.append('if [ -s "${error_log}" ]; then')
        err_msg = (
            "    echo \"Error: cellranger-arc count for "
            + base
            + " didn't fully complete (leftover tmp folders detected). "
            + "See "
            + base
            + "_error_log.txt\" >&2"
        )
        L.append(err_msg)
        L.append('else')
        L.append('    rm -f "${error_log}"')
        L.append('fi')
        L.append("")

    # Post-processing: collect per-sample web summaries + summary CSV into
    # raw_metrics_folder (matches CRarc_patched_full extract_summary_metrics),
    # then archive everything to the NGSruns tree, then clean up locally.
    L.append("# ══ Post-processing ══")
    L.append("# Mirrors CRarc_patched_full::extract_summary_metrics (l.596-657).")
    L.append('for sd in "${run_ids[@]}"; do')
    L.append("")
    # 1) Joint web summary — always present.
    L.append('    if [ -f "${Server_folder}/${sd}/outs/web_summary.html" ]; then')
    L.append('        cp "${Server_folder}/${sd}/outs/web_summary.html" \\')
    L.append('            "${metrics_folder}/${sd}_web_summary.html"')
    L.append("    fi")
    # 2) Joint metrics CSV — name differs by pipeline layout (multiome/CR>=v8
    #    use summary.csv, older GEX-only uses metrics_summary.csv).
    L.append('    if [ -f "${Server_folder}/${sd}/outs/summary.csv" ]; then')
    L.append('        cp "${Server_folder}/${sd}/outs/summary.csv" \\')
    L.append('            "${metrics_folder}/${sd}_summary.csv"')
    L.append("    elif [ -f \"${Server_folder}/${sd}/outs/metrics_summary.csv\" ]; then")
    L.append('        cp "${Server_folder}/${sd}/outs/metrics_summary.csv" \\')
    L.append('            "${metrics_folder}/${sd}_metrics_summary.csv"')
    L.append('    else')
    L.append('        echo "  [warn] no summary.csv or metrics_summary.csv for ${sd}" >&2')
    L.append("    fi")
    # 3) Optional: cell-type web summary.
    L.append('    if [ -f "${Server_folder}/${sd}/outs/cell_types/web_summary_cell_types.html" ]; then')
    L.append('        cp "${Server_folder}/${sd}/outs/cell_types/web_summary_cell_types.html" \\')
    L.append('            "${metrics_folder}/${sd}_web_summary_cell_types.html"')
    L.append("    fi")
    # 4) Optional: GEX-only internal web summary — located by stage name
    #    under the arc CS tree, independent of fork/join/chnk depth.
    L.append(
        '    gex_web=$(find "${Server_folder}/${sd}/SC_ATAC_GEX_COUNTER_CS/'
        'SC_ATAC_GEX_COUNTER/GEX_SUMMARIZE_REPORTS" \\'
    )
    L.append('        -path "*/files/web_summary.html" -type f 2>/dev/null | head -1 || true)')
    L.append('    if [ -n "$gex_web" ]; then')
    L.append('        cp "$gex_web" "${metrics_folder}/GEX_only/${sd}_web_summary.html"')
    L.append("    else")
    L.append('        echo "  [warn] no GEX-only web_summary for ${sd} (GEX_SUMMARIZE_REPORTS absent)" >&2')
    L.append("    fi")
    # 5) Optional: ATAC-only internal web summary (same layout-robust lookup).
    L.append(
        '    atc_web=$(find "${Server_folder}/${sd}/SC_ATAC_GEX_COUNTER_CS/'
        'SC_ATAC_GEX_COUNTER/_SC_ATAC_REPORTER/CREATE_WEBSUMMARY" \\'
    )
    L.append('        -path "*/files/web_summary.html" -type f 2>/dev/null | head -1 || true)')
    L.append('    if [ -n "$atc_web" ]; then')
    L.append('        cp "$atc_web" "${metrics_folder}/ATAC_only/${sd}_web_summary.html"')
    L.append("    else")
    L.append('        echo "  [warn] no ATAC-only web_summary for ${sd} (CREATE_WEBSUMMARY absent)" >&2')
    L.append("    fi")
    # 6) md5 checksums over the sample run dir (CRarc_patched_full l.582-584).
    L.append('    find "${Server_folder}/${sd}" -type f 2>/dev/null | sort \\')
    L.append('        | xargs md5sum 2>/dev/null \\')
    L.append('        | sed "s|${Server_folder}||" \\')
    L.append('        > "${Server_folder}/${sd}_md5_checksums.txt"')
    L.append("done")
    L.append("")
    L.append('echo "=== archives ==="')
    L.append('if command -v rsync >/dev/null 2>&1; then')
    L.append('    for sd in "${run_ids[@]}"; do')
    L.append('        if [ -d "${Server_folder}/${sd}" ]; then')
    L.append('            rsync -avh "${Server_folder}/${sd}" "${archive_dir}/"')
    L.append("        fi")
    L.append('        if [ -f "${Server_folder}/libraries_${sd}.csv" ]; then')
    L.append('            rsync -avh "libraries_${sd}.csv" "${archive_dir}/"')
    L.append("        fi")
    L.append('        if [ -f "${Server_folder}/${sd}_script_log.txt" ]; then')
    L.append('            rsync -avh "${sd}_script_log.txt" "${archive_dir}/"')
    L.append("        fi")
    L.append('        if [ -f "${Server_folder}/${sd}_error_log.txt" ]; then')
    L.append('            rsync -avh "${sd}_error_log.txt" "${archive_dir}/"')
    L.append("        fi")
    L.append('        if [ -f "${Server_folder}/${sd}_md5_checksums.txt" ]; then')
    L.append('            rsync -avh "${sd}_md5_checksums.txt" "${archive_dir}/"')
    L.append("        fi")
    L.append("    done")
    L.append('    if [ -d "${metrics_folder}" ]; then')
    L.append('        rsync -avh "${metrics_folder}/" "${archive_dir}/raw_metrics_folder/"')
    L.append("    fi")
    L.append("fi")
    L.append("")
    # Local cleanup, same logic as CRarc_patched_full: drop the big
    # SC_ATAC_GEX_COUNTER_CS trees first (only after web summaries are copied),
    # then delete the local workdir unless keepdata=yes.
    L.append('if [ "${keepdata}" != "yes" ]; then')
    L.append('    for sd in "${run_ids[@]}"; do')
    L.append('        if [ -d "${Server_folder}/${sd}/SC_ATAC_GEX_COUNTER_CS" ]; then')
    L.append('            rm -r "${Server_folder}/${sd}/SC_ATAC_GEX_COUNTER_CS"')
    L.append("        fi")
    L.append("    done")
    L.append('    rm -r "${Server_folder}" 2>/dev/null')
    L.append("fi")
    L.append("")
    return "\n".join(L) + "\n"


# ─── main ─────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parses one or more AVT submission Excel files and generates a "
            "cellranger-arc_AVT0XXX.sh script with one `cellranger-arc count` "
            "per biologic sample (GEX + ATAC libraries grouped under the "
            "base sample name)."
        ),
        epilog=(
            "Example:\n"
            "  python3 %(prog)s samples.xlsx\n"
            "  python3 %(prog)s -n run1.xlsx run2.xlsx run3.xlsx"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "xlsx_files", nargs="+",
        help="One or more AVT submission Excel files.",
    )
    parser.add_argument(
        "--non-interactive", "-n", action="store_true", default=False,
        help="Skip prompts; use the detected AVT run name.",
    )
    parser.add_argument(
        "--genome", choices=["human", "mouse"], default="human",
        help="Reference genome for non-interactive runs (default: human). "
             "Interactive runs always prompt.",
    )
    args = parser.parse_args()

    set_non_interactive(args.non_interactive)
    default_genome = {
        "human": HUMAN_GENOME_PATH,
        "mouse": MOUSE_GENOME_PATH,
    }[args.genome]

    all_rows, detected_avts, any_user = [], [], ""
    for path in args.xlsx_files:
        if not os.path.isfile(path) or not path.lower().endswith((".xlsx", ".xlsm")):
            print(f"Warning: skipping '{path}' (not found or not .xlsx).")
            continue
        rows, avt = parse_xlsx(path)
        if rows:
            all_rows.extend(rows)
        if avt:
            detected_avts.append(avt)
        if not any_user:
            any_user = read_user_field(path)

    if not all_rows:
        sys.exit("Error: no sample rows could be parsed from the provided file(s).")

    groups = group_libraries(all_rows)
    print()
    print_groups_table(groups)
    print()

    any_avt = latest_avt(detected_avts)
    if not any_avt:
        any_avt = ask_avt_run(None)
    else:
        any_avt = ask_avt_run(any_avt)
    if not any_avt:
        sys.exit("Error: AVT run name cannot be empty.")

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=default_genome)
    user = any_user
    run_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=any_avt)

    script = build_script(reference_genome, groups, any_avt, user, run_dir)
    out_script = f"cellranger-arc_{any_avt}.sh"
    with open(out_script, "w") as f:
        f.write(script)
    os.chmod(out_script, 0o755)
    print(f"Generated: {out_script}")
    if "/PATH/TO/" in script:
        print("Note: some fastq leaves could not be resolved here — "
              "replace the /PATH/TO/... placeholders in the script before running.")


if __name__ == "__main__":
    main()
