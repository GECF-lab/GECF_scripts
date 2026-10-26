#!/usr/bin/env python3
"""
generate_cellranger_template_gex.py

Usage:
    python3 generate_cellranger_template_gex.py file.xlsx [more.xlsx ...]
    python3 generate_cellranger_template_gex.py -n file.xlsx

GEX-only pipeline (Elisa: `CRcount_GEXonly`). Reads one or more Excel
datasheets (the `CR_count_datasheet_XXX_template.xlsx` layout) and emits a
single self-contained bash script — `cellranger-count_AVT0XXX.sh` — with
one `cellranger count` block per unique sample.

Design follows the arc/atac generator pattern (common/intake atoms, shared
parse_xlsx, one bash script out); only the emitted command and its output
layout differ.

Key differences from the arc/atac generators:

  * Binary: plain `cellranger` (CellRanger 10.1.0), not `cellranger-arc` /
    `cellranger-atac`.
  * Reference flag: `--transcriptome` (arc/atac use `--reference`).
  * Reference paths: gene-expression references under
    `cellranger_references/current_human|current_mouse` (not the
    `current_human_arc_ATAC|current_mouse_arc_ATAC` dirs arc/atac use).
  * Per-sample output leaves: `NN_<sample>_GEX` (atac: `_atac`, arc:
    `GEX_only`/`ATAC_only` under `cellranger-arc`).
  * Post-processing: no web-summary HTML to collect (bare `cellranger count`
    does not emit one) and no SC_* tree to delete. Summary.csv (run-level QC)
    is still merged across samples.
  * Archive destination:
    `/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/count/` (arc: `arc/`,
    atac: `atac/`). The `count/` subfolder already exists on the server.
"""

import argparse
import os
import re
import sys
import unicodedata
import warnings
from collections import Counter, defaultdict

# Suppress harmless openpyxl warnings about unsupported Excel features
warnings.filterwarnings("ignore", message="Data Validation extension is not supported")
warnings.filterwarnings("ignore", message="Conditional Formatting extension is not supported")

try:
    import openpyxl
except ImportError:
    sys.exit(
        "Error: the python module 'openpyxl' is required.\n"
        "Install it with: pip3 install openpyxl"
    )

# Ensure the pipelines/ package root is on sys.path so ``common.intake``
# resolves. The generator lives at ``pipelines/<assay>/`` so its parent
# dir is the package root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PIPELINES_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir))
if _PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, _PIPELINES_ROOT)

from common.intake import (  # noqa: E402
    NGS_BASE_DIR_TEMPLATE,
    AVT_PATTERN,
    FASTQ_ID_PATTERN,
    LOCAL_CORES,        # cellranger --localcores (shared default "44")
    LOCAL_MEM,          # cellranger --localmem   (shared default "150")
    norm,
    build_merged_cell_lookup,
    get_cell_value,
    read_user_field,
    set_non_interactive,
    setup_readline_completer,
    ask_avt_run,
    ask_reference_genome,
)

HUMAN_GENOME_PATH = (
    "/home/gecf/references/cellranger_references/current_human/"
    "refdata-gex-GRCh38-2024-A"
)
MOUSE_GENOME_PATH = (
    "/home/gecf/references/cellranger_references/current_mouse/"
    "refdata-gex-GRCm39-2024-A"
)
NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/count"

# Root on the *analysis* server where cellranger count output + post-processing
# are written (matches arc's convention). On that host this is a real, mounted
# path; on the dev box it is not expected to exist. Never derive from $PWD —
# the script must land in the same place regardless of the caller's cwd.
SCRATCH_ROOT = "/mnt/data/CellRangerCountOutput"


def folder_safe(name: str) -> str:
    """Reduce *name* (a requester string from the datasheet) to a safe
    path-segment for the working-folder name.

    Unicode is NFKD-decomposed and non-ASCII leftovers dropped, so accented
    text like ``Cèsar`` becomes ``Cesar``; any other character that cannot
    live in a CIFS path / URL is replaced with ``_`` and a trailing dot/space
    is stripped (Windows).  Empty input yields ``""`` so the caller can omit
    the whole ``_<user>`` segment instead of emitting a stray underscore.
    """
    if not name:
        return ""
    s = unicodedata.normalize("NFKD", name)
    s = s.encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    s = re.sub(r"_{2,}", "_", s).strip("._ ")
    return s

setup_readline_completer()  # tab-complete paths in the input() prompts


def ranked_sheet_names(wb):
    """Sheet order for parsing: User_x1x_copy first, then workbook order.
    The first sheet that exposes the required header row wins."""
    ranked = []
    if "User_x1x_copy" in wb.sheetnames:
        ranked.insert(0, "User_x1x_copy")
    for name in wb.sheetnames:
        if name not in ranked:
            ranked.append(name)
    return ranked


def parse_xlsx(path):
    """
    Returns (sample_fastq_pairs, detected_avt_run) for the given Excel file.
    sample_fastq_pairs: list of (sample_name, fastq_id, avt_run) triples.

    Three-level run-name detection:
      1) Per-row: fastq_id column, then adjacent columns
      2) Sheet-wide: full scan of all cells
      3) Filename fallback: AVT pattern from the file name
    """
    wb = openpyxl.load_workbook(path, data_only=True)
    print(f"  Sheets: {', '.join(wb.sheetnames)}")

    def _find_header(ws, merged_lookup):
        max_header_scan = min(15, ws.max_row)
        for row_idx in range(1, max_header_scan + 1):
            row_vals = [
                get_cell_value(ws, merged_lookup, row_idx, c)
                for c in range(1, ws.max_column + 1)
            ]
            prev_row_vals = (
                [
                    get_cell_value(ws, merged_lookup, row_idx - 1, c)
                    for c in range(1, ws.max_column + 1)
                ]
                if row_idx > 1
                else [""] * ws.max_column
            )
            combined = [
                f"{norm(prev)} {norm(cur)}".strip()
                for prev, cur in zip(prev_row_vals, row_vals)
            ]
            f_idx = None
            s_idx = None
            for i, v in enumerate(combined):
                if "fastq id" in v and "avid" in v:
                    f_idx = i
                if "sample name" in v:
                    s_idx = i
            if f_idx is not None and s_idx is not None:
                return f_idx + 1, s_idx + 1, row_idx
        return None

    ws = None
    fastq_col = None
    sample_col = None
    header_row_idx = None
    for candidate in ranked_sheet_names(wb):
        cand_ws = wb[candidate]
        cand_lookup = build_merged_cell_lookup(cand_ws)
        found = _find_header(cand_ws, cand_lookup)
        if found is None:
            continue
        fastq_col, sample_col, header_row_idx = found
        ws = cand_ws
        merged_lookup = cand_lookup
        if candidate != "User_x1x_copy":
            print(f"  -> 'User_x1x_copy' not found with data; using '{candidate}' (fallback scan)")
        else:
            print("  -> Using 'User_x1x_copy' sheet")
        break

    if ws is None:
        sys.exit(
            "Error: could not find both 'GECF fastq ID (avidxxxx)' and "
            "'Sample name' columns in the first 15 rows of any sheet."
        )

    data_start_row = header_row_idx + 1
    for merged_range in ws.merged_cells.ranges:
        if (
            merged_range.min_row <= header_row_idx <= merged_range.max_row
            and (
                merged_range.min_col <= fastq_col <= merged_range.max_col
                or merged_range.min_col <= sample_col <= merged_range.max_col
            )
        ):
            data_start_row = max(data_start_row, merged_range.max_row + 1)

    print(
        f"Found header on row {header_row_idx}: "
        f"fastq column {fastq_col}, sample column {sample_col}"
    )

    pairs = []
    blank_streak = 0
    for row_idx in range(data_start_row, ws.max_row + 1):
        fastq_val = get_cell_value(ws, merged_lookup, row_idx, fastq_col)
        sample_val = get_cell_value(ws, merged_lookup, row_idx, sample_col)

        fastq_val = str(fastq_val).strip() if fastq_val is not None else ""
        sample_val = str(sample_val).strip() if sample_val is not None else ""

        fastq_val = "" if fastq_val == "0" else fastq_val
        sample_val = "" if sample_val == "0" else sample_val

        if "run details" in norm(fastq_val) or "run details" in norm(sample_val):
            break

        if fastq_val and sample_val and FASTQ_ID_PATTERN.match(fastq_val):
            row_avt = None
            fq_match = AVT_PATTERN.search(fastq_val)
            if fq_match:
                row_avt = fq_match.group(0).upper()
            else:
                for c in range(1, ws.max_column + 1):
                    if c == fastq_col or c == sample_col:
                        continue
                    val = get_cell_value(ws, merged_lookup, row_idx, c)
                    if val is not None:
                        cm = AVT_PATTERN.search(str(val))
                        if cm:
                            row_avt = cm.group(0).upper()
                            break

            pairs.append((sample_val, fastq_val, row_avt))
            blank_streak = 0
        elif not fastq_val and not sample_val:
            blank_streak += 1
            if pairs and blank_streak >= 2:
                break
        else:
            print(
                f"Warning: row {row_idx} has only one of "
                f"fastq ID / sample name filled in; skipping."
            )
            blank_streak = 0

    avt_found = None
    for row in ws.iter_rows():
        for cell in row:
            if cell.value is None:
                continue
            m = AVT_PATTERN.search(str(cell.value))
            if m:
                avt_found = m.group(0).upper()
                break
        if avt_found:
            break

    if not avt_found:
        fname_match = AVT_PATTERN.search(os.path.basename(path))
        if fname_match:
            avt_found = fname_match.group(0).upper()
            print(f"  -> Detected run '{avt_found}' from file name")

    return pairs, avt_found


def print_pairs_table(pairs):
    if not pairs:
        return
    fastq_hdr, sample_hdr, avt_hdr = "Fastq ID", "Sample name", "AVT run"
    fw = max(len(fastq_hdr), max(len(f) for _, f, _ in pairs))
    sw = max(len(sample_hdr), max(len(s) for s, _, _ in pairs))
    aw = max(len(avt_hdr), max(len(a or "(none)") for _, _, a in pairs))
    sep = f"+{'-' * (fw + 2)}+{'-' * (sw + 2)}+{'-' * (aw + 2)}+"
    print(sep)
    print(f"| {fastq_hdr:<{fw}} | {sample_hdr:<{sw}} | {avt_hdr:<{aw}} |")
    print(sep)
    for sample_name, fastq_id, avt in pairs:
        print(f"| {fastq_id:<{fw}} | {sample_name:<{sw}} | {(avt or '(none)'):<{aw}} |")
    print(sep)


def find_fastq_path(avt_run, fastq_id):
    """Search recursively under the run base dir for a directory named
    `fastq_id`. Returns full path if found, else None (so the generated
    script can fall back to a FIXME_LOCATE_ marker)."""
    base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=avt_run)
    if not os.path.isdir(base_dir):
        return None
    for root, dirnames, _ in os.walk(base_dir):
        for d in dirnames:
            if d == fastq_id:
                return os.path.join(root, d)
    return None


#: Per-sample knobs emitted at the top of the generated bash script.
#: Each entry is (var_suffix, default_value, short_comment).
#:  - default_value ""  -> user can fill, blank by default (flag omitted when blank)
#:  - default_value "true"/"false" -> both values accepted; user flips
#:  - default_value <number> -> user can change
#: Old CRcount parity: each of these was a per-sample column in the
#: ``CR_count_datasheet`` (expect/force-cells mutual exclusion preserved).
_PER_SAMPLE_PARAMS = (
    ("localcores",      LOCAL_CORES, "blank => Unrestricted (no --localcores flag)"),
    ("localmem",        LOCAL_MEM,   "blank => Unrestricted (no --localmem flag)"),
    ("include_introns", "true",      "blank => omit (auto) ;  true / false"),
    ("create_bam",      "true",      "blank => omit (auto) ;  true / false"),
    ("chemistry",       "",          "blank => auto-detect"),
    ("expect_cells",    "",          "blank => omit  (MUTUALLY EXCLUSIVE with force_cells)"),
    ("force_cells",     "",          "blank => omit  (MUTUALLY EXCLUSIVE with expect_cells)"),
    ("r1_length",       "",          "blank => auto-detect"),
    ("r2_length",       "",          "blank => auto-detect"),
    ("custom_options",  "",          "free cellranger flags (old: count_special_settings)"),
)


def build_script_full(reference_genome, pairs, user=""):
    """Build the complete `cellranger count` bash script (GEX).

    Layout (top -> bottom):
      1. Constants (timestamp, reference genome).
      2. Per-sample parameter variables (one block per sample),
         editable in place; blank a value to drop the flag for that sample.
      3. Validation: ``--expect-cells`` and ``--force-cells`` are mutually
         exclusive per CRcount v2.4; if both are non-blank for any
         sample the script exits 1 *before* any side effects.
      4. Side effects: working dir, NGS-runs archive dir, reference
         genome existence check, FASTQ path resolution.
      5. One `cellranger count` call per sample, each using its own
         *_<idx> variables.
      6. Post-processing: merge per-sample ``summary.csv`` into one
         ``all_summary.csv`` under ``<avt>_Summaries/``.
      7. rsync the working dir up to the shared NGS-runs archive.

    Old-CRcount parity (recovered):
      * ``--cell-annotation-model=auto``
        (the sibling arc/atac generators already emit this; GEX lost it
         when it was rewritten to consume the Excel intake)
      * ``--tenx-cloud-token-path=/home/gecf/.config/txg/credentials``
      * ``--sample=<name>`` alongside ``--id``
      * Per-sample ``--chemistry`` / ``--expect-cells`` /
        ``--force-cells`` / ``--include-introns`` / ``--create-bam`` /
        ``--localcores`` / ``--localmem`` / ``--r1-length`` /
        ``--r2-length`` and a free ``--flags`` string
        (old: ``count_special_settings``)
    """
    lines = []

    avt_counts = Counter(avt for _, _, avt in pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    samples_to_fastqs = defaultdict(list)
    for sn, fq, avt in pairs:
        samples_to_fastqs[sn].append((fq, avt))
    unique_samples = sorted(samples_to_fastqs.keys())
    n = len(unique_samples)

    # ---- 1) header / constants ----
    lines.append("#!/bin/bash")
    lines.append(f"# cellranger count — GEX-only ({main_avt})")
    lines.append("# Generated by generate_cellranger_template_gex.py "
                 f"({n} sample(s)).")
    lines.append(f"# Requester: {user or 'n/a'}")
    lines.append("# Restored from old CRcount: --cell-annotation-model=auto, "
                 "--tenx-cloud-token-path, --sample, per-sample knobs,")
    lines.append("# and the expect-cells / force-cells mutual-exclusion rule.")
    lines.append("")
    lines.append("TIMESTAMP=$(date +%Y%m%d-%H%M%S)")
    lines.append(f"reference_genome_used={reference_genome}")
    lines.append("")

    # ---- 2) per-sample parameter variables (editable) ----
    lines.append("# ============================================================")
    lines.append("# PER-SAMPLE PARAMETERS  (one block per sample)")
    lines.append("# ----------------------------------------------------------------")
    lines.append('#   * blank a value  ( value="" )  -> flag is dropped for that sample')
    lines.append("#   * number / text -> emitted as --flag=<value>")
    lines.append("#   * localcores / localmem : blank = Unrestricted  (cellranger default)")
    lines.append("#   * include_introns / create_bam : blank = omit  (cellranger auto);")
    lines.append("#                                     default is \"true\" (old behavior)")
    lines.append("#   * expect_cells and force_cells are MUTUALLY EXCLUSIVE,")
    lines.append("#     per the old CRcount rule. Both non-blank => this script exit 1s.")
    lines.append("#   * custom_options: free-form string of extra --flags, untouched")
    lines.append("#     (blank for none; otherwise passed as-is)")
    lines.append("# ============================================================")
    for idx, sn in enumerate(unique_samples, start=1):
        lines.append("")
        lines.append(f"# --- Sample {idx}: {sn} ---")
        for suffix, default, comment in _PER_SAMPLE_PARAMS:
            base = f'{suffix}_{idx}="{default}"'
            pad = max(1, 60 - len(base))
            lines.append(f"{base}   {' ' * pad}# {comment}")
    lines.append("")

    # ---- 3) validation (MUST run before any side effects) ----
    lines.append("# ===== VALIDATION: --expect-cells vs --force-cells are mutually exclusive (old CRcount rule) =====")
    lines.append("_vfail=0")
    for i in range(1, n + 1):
        lines.append(f'if [[ -n "${{expect_cells_{i}:-}}" && -n "${{force_cells_{i}:-}}" ]]; then')
        lines.append(f'    echo "ERROR: sample {i} has both expect_cells_{i} and force_cells_{i} set. Pick one; clear the other." >&2')
        lines.append("    _vfail=1")
        lines.append("fi")
    lines.append("if [[ ${_vfail} -ne 0 ]]; then")
    lines.append('    echo "Aborting: expect-cells / force-cells are mutually exclusive per sample (old CRcount rule). See errors above." >&2')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")

    # ---- 4) side effects (working dirs, archive, reference presence, FASTQ paths) ----
    # Requester tag, reduced to a safe path-segment (accented names like
    # "Cèsar" -> "Cesar").  Empty/sanitized-blank requester -> the whole
    # " _<user>" segment is omitted so we never emit a stray underscore.
    folder_user = folder_safe(user)
    user_tag = "_" + folder_user if folder_user else ""
    # Working dir: FIXED root on the analysis server (see SCRATCH_ROOT),
    # never derived from $PWD so the output lands identically whichever
    # directory the script is launched from.
    lines.append("# Working directory: cellranger count output + post-processing (fixed root, arc convention)")
    lines.append(f'scratch_root="{SCRATCH_ROOT}"')
    lines.append(f'server_folder_name="{main_avt}{user_tag}_$TIMESTAMP"')
    lines.append('Server_folder="${scratch_root}/${server_folder_name}"')
    lines.append('mkdir -vp "${Server_folder}"')
    lines.append("")
    lines.append("# Central archive destination (fixed target, shared with arc/atac)")
    lines.append("NGSruns_archive=" + NGSruns_ARCHIVE_DIR)
    lines.append('mkdir -p "${NGSruns_archive}/${server_folder_name}"')
    lines.append("NGSruns_folder=${NGSruns_archive}/${server_folder_name}")
    lines.append("")
    lines.append("# Local workdir cleanup (old CRcount keepdata semantics).")
    lines.append('#   After deploy to the archive, the local working folder is')
    lines.append('#   deleted UNLESS keepdata="yes".')
    lines.append('keepdata=""')
    lines.append("")
    lines.append("# Check that the reference genome directory exists")
    lines.append('if [[ ! -d "${reference_genome_used}" ]]; then')
    lines.append('    echo "Reference genome ${reference_genome_used} does not exist. Aborting." >&2')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")
    lines.append("# ---- FASTQ path resolution (auto-detected where possible) ----")
    lines.append("# Generated by the script at generation time:")
    lines.append("#   /mnt/PTEGraw/AA/NGSruns/AVT/AV240401/<run>/<avt_fastq_id>/")
    lines.append("# where <run> is the AVT run that carried the sample; if a")
    lines.append("# sample is in multiple AVT runs, the first one (by pair order)")
    lines.append("# is used. FIXME_LOCATE_ markers mean the generator could not")
    lines.append("# find the folder — set those path_*_N variables by hand.")

    any_missing = False
    for idx, sn in enumerate(unique_samples, start=1):
        fastqs = samples_to_fastqs[sn]
        lines.append("")
        lines.append(f"# --- Sample {idx}: {sn} ({len(fastqs)} FASTQ(s)) ---")
        for fq, avt in fastqs:
            run_part = avt or main_avt
            base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=run_part)
            resolved = find_fastq_path(run_part, fq)
            if resolved:
                lines.append(f"path_{idx}_{fq}={resolved}")
            else:
                lines.append(f"path_{idx}_{fq}={base_dir}/FIXME_LOCATE_{fq}")
                lines.append("# (generator could not resolve the folder — fix by hand)")
                if os.path.isdir(base_dir):
                    print(
                        f"  Warning: no folder named '{fq}' found under "
                        f"{base_dir} — path_{idx}_{fq} needs to be set manually.",
                        file=sys.stderr,
                    )
                else:
                    any_missing = True

    if any_missing:
        print("  Note: some base dirs are not on this machine —", file=sys.stderr)
        print("  all flagged path_* variables (FIXME_LOCATE_ markers) must be set before running.",
              file=sys.stderr)

    # ---- 5) one cellranger count per sample, each with its old-CRcount
    #          post-processing (SC tree cleanup, completion check, per-sample
    #          deploy, md5 checksums) ----
    _CONT = " \\"    # bash line-continuation used between flags
    for idx, sn in enumerate(unique_samples, start=1):
        fastqs = samples_to_fastqs[sn]
        numbered_dir = f"{idx:02d}_{sn}_GEX"
        fq_list_str = ",".join(fq for fq, _ in fastqs)
        fq_paths = ", ".join(f"${{path_{idx}_{fq}}}" for fq, _ in fastqs)
        SFR = "${Server_folder}"

        flag_lines = [
            'cellranger count --id="' + sn + '"',
            '    --sample="' + sn + '"',
            '    --description="' + fq_list_str + '"',
            '    --transcriptome="${reference_genome_used}"',
            '    --fastqs="' + fq_paths + '"',
            '    --cell-annotation-model=auto',
            '    --tenx-cloud-token-path="/home/gecf/.config/txg/credentials"',
            '    --output="' + SFR + '/' + numbered_dir + '"',
        ]
        # Per-sample conditional knobs. Blank value -> flag omitted for this
        # sample only. custom_options is deliberately UNQUOTED on the last
        # line so it word-splits into multiple --flags (old CRcount
        # count_special_settings behavior).
        for suffix, _default, _comment in _PER_SAMPLE_PARAMS:
            v = f"{suffix}_{idx}"
            if suffix == "custom_options":
                flag_lines.append(f"    ${{{v}}}")
            else:
                flag = "--" + suffix.replace("_", "-")
                flag_lines.append(f"    ${{{v}:+{flag}=${v}}}")

        lines.append("# ============================================================")
        lines.append(f"# EXECUTE  (one cellranger count per sample; each pulls from")
        lines.append("# its own *_<idx> parameters above; the flag is only emitted")
        lines.append("# when that variable is non-blank)")
        lines.append("# Per-sample post-processing replicates old CRcount:")
        lines.append("#   rm -r <sample>/SC_RNA_COUNTER_CS")
        lines.append("#   find <sample> -type d -name tmp  (completion check)")
        lines.append("#   rsync -avh <sample>/ + logs + md5 checksums -> archive")
        lines.append("# ============================================================")
        lines.append("")
        lines.append(f"# --- Run sample {idx}: {sn} ---")
        lines.append(
            'echo "cellranger count --id=' + sn + ' --sample=' + sn
            + ' ... (see full command below)" '
            + '> "' + SFR + '/' + numbered_dir + '_script_log.txt"'
        )
        # Join flag_lines with the line-continuation between them.
        for i, fl in enumerate(flag_lines):
            lines.append(fl + (_CONT if i < len(flag_lines) - 1 else ""))

        # Old CRcount per-sample post-processing.
        lines.append("")
        lines.append(f"# --- Post-processing sample {idx}: {sn} ---")
        lines.append(
            'rm -r "' + SFR + '/' + numbered_dir + '/SC_RNA_COUNTER_CS" 2>/dev/null || true'
        )
        lines.append(
            'find "' + SFR + '/' + numbered_dir
            + '" -type d -name "tmp" -print '
            '> "' + SFR + '/' + numbered_dir + '_error_log.txt"'
        )
        lines.append('if [ -s "' + SFR + '/' + numbered_dir + '_error_log.txt" ]; then')
        lines.append(
            '    echo "Error: CR count didn'"'"'t complete for sample '
            + sn + '"' + ' >&2'
        )
        lines.append("else")
        lines.append('    rm -f "' + SFR + '/' + numbered_dir + '_error_log.txt"')
        lines.append("fi")
        lines.append('rsync -avh "' + SFR + '/' + numbered_dir + '" "${NGSruns_folder}"')
        lines.append(
            'rsync -avh "' + SFR + '/' + numbered_dir + '_script_log.txt" "${NGSruns_folder}"'
        )
        lines.append(
            'rsync -avh "' + SFR + '/' + numbered_dir + '_error_log.txt" "${NGSruns_folder}" 2>/dev/null'
        )
        # md5 checksum (old CRcount: cd into sample dir, find -type f, xargs md5sum,
        # output file is SIBLING to the sample dir with relative ./ paths inside)
        lines.append("# md5 checksum (old CRcount)")
        lines.append(
            '( cd "' + SFR + '/' + numbered_dir
            + '" && find . -type f | sort'
            + ' | xargs md5sum > "../' + numbered_dir + '_md5_checksums.txt" )'
        )
        lines.append(
            'rsync -avh "' + SFR + '/' + numbered_dir + '_md5_checksums.txt" "${NGSruns_folder}"'
        )
        lines.append("")

    # ---- 6) extract_summary_metrics (old CRcount parity) ----
    #   copy per-sample web_summary.html / web_summary_cell_types.html /
    #   metrics_summary.csv into ${Server_folder}/raw_metrics_folder renamed
    #   as <sample>_... ; strip comma thousands-separators inside quoted
    #   fields (awk -F'"'); strip "N%" -> "N" (sed); build summary_metrics.csv
    #   with Samplename column; reorder to metric order 1,2,3,8,4,21,20,19,13,17.
    #   (Old CRcount::extract_summary_metrics, lines 651-717.)
    SFR = "${Server_folder}"
    # literal sample-name list (no dynamic expansion in the emitted bash)
    _sample_list = " ".join(
        str(i).zfill(2) + "_" + s + "_GEX" for i, s in enumerate(unique_samples, 1)
    )
    lines.append("# ------------------------------------ EXTRACT SUMMARY METRICS ------------------------------------ #")
    lines.append("# (replicates old CRcount::extract_summary_metrics)")
    lines.append("")
    lines.append('metrics_folder="' + SFR + '/raw_metrics_folder"')
    lines.append('mkdir -p "${metrics_folder}"')
    lines.append('for _sd in ' + _sample_list + '; do')
    lines.append("    # web_summary.html")
    lines.append('    cp "' + SFR + '/${_sd}/outs/web_summary.html" "${metrics_folder}/${_sd}_web_summary.html" 2>/dev/null')
    lines.append("    # web_summary_cell_types.html")
    lines.append('    cp "' + SFR + '/${_sd}/outs/cell_types/web_summary_cell_types.html" "${metrics_folder}/${_sd}_web_summary_cell_types.html" 2>/dev/null')
    lines.append("    # metrics_summary.csv")
    lines.append('    src="' + SFR + '/${_sd}/outs/metrics_summary.csv"')
    lines.append('    dst="${metrics_folder}/${_sd}_metrics_summary.csv"')
    lines.append('    if [ -f "${src}" ]; then')
    lines.append('        cp "${src}" "${dst}"')
    lines.append('    else')
    lines.append('        echo "WARNING: no metrics_summary.csv for ${_sd}" >&2')
    lines.append('        continue')
    lines.append('    fi')
    lines.append('    tmp="${dst}_tmp.csv"')
    lines.append("    # strip comma thousands-separators inside quoted fields (old CRcount)")
    lines.append(r'''    awk -F '"' 'BEGIN {OFS = "\""} { for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' "${dst}" > "${tmp}"''')
    lines.append("    # strip % symbols (old CRcount)")
    lines.append('    sed -E "s/([0-9]+)%/\\1/g" "${tmp}" > "${dst}"')
    lines.append('    rm -f "${tmp}"')
    lines.append("done")
    lines.append("")
    lines.append("# Build summary_metrics.csv with Samplename column (old CRcount)")
    lines.append("first_csv_processed=false")
    lines.append("for _sd in " + _sample_list + "; do")
    lines.append('    file="${metrics_folder}/${_sd}_metrics_summary.csv"')
    lines.append('    if [ ! -f "${file}" ]; then')
    lines.append('        echo "WARNING: ${file} not found, skipping" >&2')
    lines.append('        continue')
    lines.append('    fi')
    lines.append('    if [ "$first_csv_processed" = false ]; then')
    lines.append('        echo -n "Samplename," > summary_metrics_tmp.csv')
    lines.append('        sed -n "1p" "${file}" >> summary_metrics_tmp.csv')
    lines.append('        first_csv_processed=true')
    lines.append('    fi')
    lines.append('    echo -n "${_sd}," >> summary_metrics_tmp.csv')
    lines.append('    sed -n "2p" "${file}" >> summary_metrics_tmp.csv')
    lines.append("done")
    lines.append("# Reorder to classic CRcount metric order (old CRcount):")
    lines.append("#   1, 2, 3, 8, 4, 21, 20, 19, 13, 17")
    lines.append("awk -F',' '{ print $1, $2, $3, $8, $4, $21, $20, $19, $13, $17 }' OFS=',' \"summary_metrics_tmp.csv\" > \"summary_metrics.csv\"")
    lines.append("rm -f summary_metrics_tmp.csv")
    lines.append('rsync -avh "${metrics_folder}" "${NGSruns_folder}"')
    lines.append('rsync -avh "summary_metrics.csv" "${NGSruns_folder}"')
    lines.append("")

    # ---- 7) keepdata cleanup (old CRcount) ----
    # After all deploy, remove the local workdir unless keepdata=yes.
    lines.append("# ------------------------------------ LOCAL WORKDIR CLEANUP ------------------------------------ #")
    lines.append('# Old CRcount: `if test "$keepdata" != "yes"; then rm -r ${pathcount} 2> /dev/null; fi`')
    lines.append('if [[ "${keepdata}" != "yes" ]]; then')
    lines.append('    rm -r "${Server_folder}" 2>/dev/null || true')
    lines.append("fi")
    lines.append("")

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parses one or more GEX-only Excel datasheets and generates a "
            "cellranger-count_AVT0XXX.sh script with one `cellranger count` "
            "block per sample."
        ),
        epilog="Example:\n  python3 %(prog)s samples.xlsx\n  python3 %(prog)s run1.xlsx run2.xlsx",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "xlsx_files",
        nargs="+",
        help=(
            "One or more Excel files with the 'GECF fastq ID (avidxxxx)' and "
            "'Sample name' columns (CR_count_datasheet_*_template.xlsx layout). "
            "AVT run name (AVT0XXX) is also auto-detected."
        ),
    )
    parser.add_argument(
        "--non-interactive", "-n",
        action="store_true",
        default=False,
        help="Skip prompts; use detected values and human genome as default.",
    )
    args = parser.parse_args()

    set_non_interactive(args.non_interactive)
    xlsx_files = args.xlsx_files

    all_pairs = []
    for xlsx_file in xlsx_files:
        if not os.path.isfile(xlsx_file):
            sys.exit(f"Error: file '{xlsx_file}' does not exist.")
        if not xlsx_file.lower().endswith(".xlsx"):
            sys.exit(f"Error: '{xlsx_file}' does not look like an .xlsx file.")

        pairs, detected_avt = parse_xlsx(xlsx_file)
        if not pairs:
            print(f"Warning: no sample/fastq ID pairs found in '{xlsx_file}' — skipping.")
            continue
        all_pairs.extend(
            (s, f, a if a is not None else detected_avt)
            for s, f, a in pairs
        )

    if not all_pairs:
        sys.exit("Error: no sample/fastq ID pairs could be parsed from the provided file(s).")

    print()
    print(f"Total: {len(all_pairs)} sample(s) from {len(xlsx_files)} Excel file(s):")
    print_pairs_table(all_pairs)
    print()

    needs_avt = [i for i, (_, _, avt) in enumerate(all_pairs) if avt is None]
    if needs_avt:
        known_avts = list(dict.fromkeys(avt for _, _, avt in all_pairs if avt is not None))
        detected_hint = known_avts[0] if known_avts else None
        print(f"{len(needs_avt)} sample(s) have no AVT run assigned.")
        avt_run = ask_avt_run(detected_hint)
        if not avt_run:
            sys.exit("Error: AVT run name cannot be empty.")
        for i in needs_avt:
            all_pairs[i] = (all_pairs[i][0], all_pairs[i][1], avt_run)

    avt_counts = Counter(avt for _, _, avt in all_pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=HUMAN_GENOME_PATH
    )
    user = read_user_field(xlsx_files[0])

    script_content = build_script_full(reference_genome, all_pairs, user)
    out_script = f"cellranger-count_{main_avt}.sh"
    with open(out_script, "w") as f:
        f.write(script_content)
    os.chmod(out_script, 0o755)

    print(f"Generated: {out_script}")
    if "FIXME_LOCATE_" in script_content:
        print(
            "Note: some fastq paths could not be resolved — set the flagged "
            "path_* variables (FIXME_LOCATE_ markers) before running."
        )


if __name__ == "__main__":
    main()
