#!/usr/bin/env python3
"""
generate_cellranger_template_atac.py

Usage:
    python3 generate_cellranger_template_atac.py file.xlsx

Reads an Excel sheet containing (at least) these columns:
    - "GECF fastq ID (avidxxxx)"
    - "Sample name" (trailing spaces in the header are tolerated)
and a cell somewhere in the sheet containing the run name (AVT0XXX).

It then asks which reference genome to use, and generates a
cellranger-atac_AVT0XXX.sh script pre-filled with one
"cellranger-atac count" block per sample.
"""

import argparse
import os
import re
import sys
import warnings

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
# resolves.  The generator lives at ``pipelines/<assay>/`` so its parent
# directory is the package root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_PIPELINES_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir))
if _PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, _PIPELINES_ROOT)

# Shared intake atoms — see pipelines/common/intake.py.  The constants are
# re-exported here so any downstream code that imported them from this
# module keeps working.
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

HUMAN_GENOME_PATH = "/home/gecf/references/cellranger_references/current_human_arc_ATAC/refdata-cellranger-arc-GRCh38-2024-A"
MOUSE_GENOME_PATH = "/home/gecf/references/cellranger_references/current_mouse_arc_ATAC/refdata-cellranger-arc-GRCm39-2024-A"
# Central archive destination for cellranger-atac runs (fixed target).
NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/atac"

# A real GECF library id is `avid` followed by 4-6 digits (e.g. avid12499).
# Label/legend/placeholder rows in the same columns carry strings like
# "0", "GECF libr. ID (avidxxxx)" or "Loaded molarity (pM)" — this pattern
# distinguishes genuine sample rows from those.
FASTQ_ID_PATTERN = re.compile(r"^avid\d{4,6}$", re.IGNORECASE)


setup_readline_completer()  # tab-complete paths in the input() prompts


def ranked_sheet_names(wb):
    """Sheet order for parsing: User_x1x_copy first (the canonical
    single-user block), then all remaining sheets in workbook order.
    Used by parse_xlsx: the first sheet that exposes the required
    header row wins, so the canonical sheet is always preferred and
    the rest of the workbook is only fallback."""
    ranked = []
    if "User_x1x_copy" in wb.sheetnames:
        ranked.insert(0, "User_x1x_copy")
    for name in wb.sheetnames:
        if name not in ranked:
            ranked.append(name)
    return ranked


# ─── main parsing logic ───────────────────────────────────────────────

def parse_xlsx(path):
    """
    Returns (sample_fastq_pairs, detected_avt_run) for the given Excel file.
    
    sample_fastq_pairs: list of (sample_name, fastq_id, avt_run) triples
    
    Three-level run-name detection:
      1) Per-row: look in the fastq_id column and adjacent columns on each data row
      2) Sheet-wide: full scan of all cells 
      3) Filename fallback: extract AVT pattern from the file name itself
    """
    wb = openpyxl.load_workbook(path, data_only=True)
    print(f"  Sheets: {', '.join(wb.sheetnames)}")

    # Prefer the canonical single-user sheet (User_x1x_copy); otherwise scan
    # every worksheet in workbook order and take the first one that exposes
    # the required header row. This is resilient to multi-sheet GECF
    # workbooks where the data block is no longer the first sheet (e.g. the
    # 15-sheet layout with GUIDELINES first).
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
            print(f"  -> Using 'User_x1x_copy' sheet")
        break

    if ws is None:
        sys.exit(
            "Error: could not find both 'GECF fastq ID (avidxxxx)' and "
            "'Sample name' columns in the first 10 rows of any sheet."
        )

    # If the header cell(s) are merged across multiple rows, make sure we
    # start reading data *after* the full merged range, not right after
    # the row where we first detected the header text.
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

        # Excel returns 0 for empty numeric-formatted cells — treat those
        # the same as genuinely blank.
        fastq_val = "" if fastq_val == "0" else fastq_val
        sample_val = "" if sample_val == "0" else sample_val

        # Stop immediately if we hit a "Run details" marker row (a common
        # section break some sheets use after the sample table)
        if "run details" in norm(fastq_val) or "run details" in norm(sample_val):
            break

        if fastq_val and sample_val and FASTQ_ID_PATTERN.match(fastq_val):
            # Level 1 – per-row AVT detection: try fastq_id first, then other columns
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
            # Fully empty row: could be a stray blank row, or the real end
            # of the table. Tolerate isolated blank rows, but stop once we
            # see two blank rows in a row (after we've already collected
            # at least one pair).
            blank_streak += 1
            if pairs and blank_streak >= 2:
                break
        else:
            # One of the two cells is empty but not both: likely a
            # malformed/partial row. Skip it but keep going.
            print(
                f"Warning: row {row_idx} has only one of "
                f"fastq ID / sample name filled in; skipping."
            )
            blank_streak = 0

    # Level 2 – search whole sheet for an AVT run name
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

    # Level 3 – fallback: extract AVT from the file name itself
    if not avt_found:
        fname_match = AVT_PATTERN.search(os.path.basename(path))
        if fname_match:
            avt_found = fname_match.group(0).upper()
            print(f"  -> Detected run '{avt_found}' from file name")

    return pairs, avt_found

# (read_user_field is imported from common.intake — see the import block
# near the top of this file.)


def print_pairs_table(pairs):
    """
    Print a simple aligned table of (sample_name, fastq_id, avt_run) triples,
    with the fastq ID column first.
    """
    if not pairs:
        return

    fastq_header = "Fastq ID"
    sample_header = "Sample name"
    avt_header = "AVT run"

    fastq_width = max(len(fastq_header), max(len(f) for _, f, _ in pairs))
    sample_width = max(len(sample_header), max(len(s) for s, _, _ in pairs))
    avt_width = max(len(avt_header), max(len(a or "(none)") for _, _, a in pairs))

    sep = f"+{'-' * (fastq_width + 2)}+{'-' * (sample_width + 2)}+{'-' * (avt_width + 2)}+"

    print(sep)
    print(f"| {fastq_header:<{fastq_width}} | {sample_header:<{sample_width}} | {avt_header:<{avt_width}} |")
    print(sep)
    for sample_name, fastq_id, avt in pairs:
        print(f"| {fastq_id:<{fastq_width}} | {sample_name:<{sample_width}} | {(avt or '(none)'):<{avt_width}} |")
    print(sep)

# (ask_avt_run / ask_reference_genome come from common/intake — see the
# import block near the top of this file.)


def build_script(reference_genome, pairs):
    """Generate the cellranger-atac shell script content."""

    # Derive a main AVT for output filenames (most common one, or MIXED)
    from collections import Counter
    avt_counts = Counter(avt for _, _, avt in pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    # Build per-sample fastq path variables with timestamps and custom output dirs
    NGS_BASE_DIR_TEMPLATE = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/{avt_run}"

    lines.append(f"# Generated by pipelines/atac/generate_cellranger_template_atac.py")
    lines.append("# Review fastq paths and the reference genome before running.")
    lines.append("")

    # Build per-fastq variables. Each uses the AVT embedded in its path string.
    for idx, (sample_name, fastq_id, avt) in enumerate(pairs, start=1):
        run_part = avt or main_avt
        base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=run_part)
        lines.append(
            f"# --- {idx}. {sample_name} ({fastq_id}, AVT: {run_part}) ---"
        )
        # We will fill in the full path interactively or after review
        avt_path = base_dir  # Will need manual fine-tuning for exact subfolder
        lines.append(
            f"# fastq_path_{idx}={avt_path}/{fastq_id}"
        )

    return "\n".join(lines) + "\n"


# Global list for building output script
lines = []


def find_fastq_path(avt_run, fastq_id):
    """
    Search recursively under the run's base directory for a directory
    whose name is exactly `fastq_id` (it can live under any subfolder
    structure). Returns the full path if found, or None if not found
    (e.g. base dir doesn't exist on this machine, or no matching folder).
    """
    base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=avt_run)

    if not os.path.isdir(base_dir):
        return None

    for root, dirnames, _ in os.walk(base_dir):
        for d in dirnames:
            if d == fastq_id:
                return os.path.join(root, d)

    return None


def build_script_full(reference_genome, pairs, user=""):
    """Build complete cellranger-atac shell script.

    user is the requester read from the workbook (first word). Woven into
    the per-sample output folder; empty leaves the folder without a user
    segment (no placeholder).
    """
    global lines
    lines = []

    # Derive a main AVT for output filenames (most common one, or MIXED)
    from collections import Counter
    avt_counts = Counter(avt for _, _, avt in pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    lines.append("#!/bin/bash")
    lines.append("TIMESTAMP=$(date +%Y%m%d-%H%M%S)")
    lines.append("reference_genome_used=" + reference_genome)
    lines.append("")

    # Group fastq IDs per sample name so each unique sample gets ONE count call
    from collections import defaultdict, Counter as Ctr

    samples_to_fastqs = defaultdict(list)
    for sn, fq, avt in pairs:
        samples_to_fastqs[sn].append((fq, avt))

    unique_samples = sorted(samples_to_fastqs.keys())

    # Probe real filesystem for each fastq dir and build path variables
    any_missing = False
    for sn in unique_samples:
        fastqs = samples_to_fastqs[sn]
        lines.append("# --- " + sn + " (" + str(len(fastqs)) + " FASTQ(s)) ---")
        for fq, avt in fastqs:
            run_part = avt or main_avt
            base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=run_part)

            resolved_path = find_fastq_path(run_part, fq)
            if resolved_path:
                lines.append("path_" + sn + "_" + fq + "=" + resolved_path)
            else:
                lines.append("# path_" + sn + "_" + fq + "=" + base_dir + "/FIXME_LOCATE_" + fq)
                if os.path.isdir(base_dir):
                    print(
                        "  Warning: no folder named '" + fq + "' found under "
                        + base_dir + " - path_" + sn + "_" + fq + " needs to be set manually."
                    )
                else:
                    any_missing = True
        lines.append("")

    if any_missing:
        print("  Note: some base directories are not on this machine —")
        print("  all fastq_path variables under those runs need to be set manually.")
    lines.append("")
    lines.append("# Extra cellranger-atac options (leave empty if none)")
    lines.append("# Examples: --chemistry=ARC-v1  --no-bam  --min-atac-count=500")
    lines.append("CUSTOM_OPTIONS=")
    lines.append("")

    # Working folder name: <AVT>[_<user>]_<YYYYmmdd-HHMMSS>. Timestamp
    # expands at run time (bash, same slot as arc/multi do it); user is
    # baked in here because the workbook already told us who runs this.
    user_tag = "_" + user if user else ""
    # Server_folder_name must be a SHELL variable in the emitted script
    # (it is referenced below as ${server_folder_name} at run time).
    lines.append("")
    lines.append("# Working directory: cellranger-atac output and post-processing live here")
    lines.append(f'server_folder_name="{main_avt}{user_tag}_$TIMESTAMP"')
    lines.append('Server_folder="$PWD/${server_folder_name}"')
    lines.append('mkdir -vp "$Server_folder"')
    lines.append("")
    # Final archive destination (fixed central archive for all cellranger-atac runs)
    lines.append("# Central archive for cellranger-atac runs (fixed destination)")
    lines.append("NGSruns_archive=" + NGSruns_ARCHIVE_DIR)
    lines.append("mkdir -p \"${NGSruns_archive}/${server_folder_name}\"")
    lines.append("NGSruns_folder=${NGSruns_archive}/${server_folder_name}")
    lines.append("")
    lines.append("# Check that the reference genome directory exists")
    lines.append('if [[ ! -d "${reference_genome_used}" ]]; then')
    lines.append("    echo \"Reference genome '${reference_genome_used}' does not exist. Aborting.\" >&2")
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")
    lines.append("")

    for idx, sn in enumerate(unique_samples, start=1):
        fastqs = samples_to_fastqs[sn]
        numbered_dir = f"{idx:02d}_{sn}_atac"
        fq_list_str = ",".join(fq for fq, _ in fastqs)

        lines.append("# --- Sample " + str(idx) + ": " + sn + " (" + str(len(fastqs)) + " FASTQ(s)) ---")
        lines.append("cellranger-atac count --id=" + sn)
        lines.append("                       --description=\"" + fq_list_str + "\" \\")
        lines.append("                       --reference=${reference_genome_used} \\")

        # Build space-separated path list for all FASTQ dirs of this sample
        fq_paths = " ".join("${path_" + sn + "_" + fq + "}" for fq, _ in fastqs)
        lines.append("                       --fastqs=" + fq_paths + r" \\")
        lines.append("                       --output=${Server_folder}/" + numbered_dir + r" \\")
        lines.append("                       ${CUSTOM_OPTIONS}")
        lines.append("")

        lines.append("# Clean up SC_ATAC_COUNTER_CS (non-useful)")
        lines.append(
            "rm -r ${Server_folder}/" + numbered_dir + "/SC_ATAC_COUNTER_CS 2>/dev/null || true"
        )
        lines.append("")

    post_avt = main_avt or "MIXED_runs"
    # Post-processing in Server_folder; rsync everything to NGSruns_folder at the end
    lines.append("")
    lines.append("# ------------------------------------ POST-PROCESSING ------------------------------------ #")
    lines.append("")
    lines.append("# Collect per-sample web summary HTMLs into one folder for easy browsing")
    lines.append(f"mkdir -p ${{Server_folder}}/{post_avt}_Summaries/web_summaries")
    lines.append("")
    for idx, sn in enumerate(unique_samples, start=1):
        ndir = f'{idx:02d}_{sn}_atac/outs'
        lines.append(
            f"cp -n ${{Server_folder}}/{ndir}/web_summary.html "
            f"${{Server_folder}}/{post_avt}_Summaries/web_summaries/{sn}_WebSummary.html 2>/dev/null || true"
        )
    lines.append("")

    # Merge per-sample summary.csv (run-level QC metrics) across all samples into one file
    lines.append(
        "# Merge run-level QC metrics (summary.csv) from each sample into a single consolidated table"
    )
    lines.append(
        f'merged_csv="${{Server_folder}}/{post_avt}_Summaries/all_summary.csv"'
    )
    lines.append("header_csv_set=false")
    for idx, sn in enumerate(unique_samples, start=1):
        ndir = f'{idx:02d}_{sn}_atac/outs'
        src = f'${{Server_folder}}/{ndir}/summary.csv'
        lines.append(
            'if [ -f "' + src + '" ]; then '
            'if [ "$header_csv_set" = false ]; then '
            'cp "' + src + '" "$merged_csv"; '
            'header_csv_set=true; else '
            'tail -n +2 "' + src + '" >> "$merged_csv"; fi; '
            'else echo "WARNING: ' + src + ' not found" >&2; fi'
        )
    lines.append("")

    # Rsync everything from Server_folder to NGSruns_folder in one shot
    lines.append("")
    lines.append("# ------------------------------------ DEPLOY TO NGS RUNS FOLDER ------------------------------------ #")
    lines.append("")
    lines.append("# Rsync all results (sample dirs + Summaries folder) to the permanent location in one command")
    lines.append("rsync -a --progress ${Server_folder}/ ${NGSruns_folder}")
    lines.append("")
    return "\n".join(lines) + "\n"


# ─── main entry point ─────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parses one or more Excel sheets and generates a cellranger-atac_AVT0XXX.sh "
            "script with one \"cellranger-atac count\" block per sample."
        ),
        epilog="Example:\n  python3 %(prog)s samples.xlsx\n  python3 %(prog)s run1.xlsx run2.xlsx run3.xlsx",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "xlsx_files",
        nargs="+",
        help=(
            "One or more Excel files containing at least the columns "
            "'GECF fastq ID (avidxxxx)' and 'Sample name'. "
            "A cell containing the run name (e.g. AVT0226) is also "
            "searched for automatically."
        ),
    )
    parser.add_argument(
        "--non-interactive", "-n",
        action="store_true",
        default=False,
        help="Skip all prompts, use detected values and human genome as defaults.",
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

        # Backfill any pair whose per-row AVT is None with the sheet-level detection
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

    # Assign a run name to any pairs that have None
    needs_avt = [i for i, (_, _, avt) in enumerate(all_pairs) if avt is None]
    if needs_avt:
        known_avts = list(dict.fromkeys(avt for _, _, avt in all_pairs if avt is not None))
        detected_avt_hint = known_avts[0] if known_avts else None

        print(f"{len(needs_avt)} sample(s) have no AVT run assigned.")
        if detected_avt_hint:
            avt_run = ask_avt_run(detected_avt_hint)
        else:
            avt_run = ask_avt_run(None)

        if not avt_run:
            sys.exit("Error: AVT run name cannot be empty.")

        for i in needs_avt:
            all_pairs[i] = (all_pairs[i][0], all_pairs[i][1], avt_run)

    # Derive a representative AVT for the output filename: pick the most common one
    from collections import Counter
    avt_counts = Counter(avt for _, _, avt in all_pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=HUMAN_GENOME_PATH)

    # Read the requester from the workbook (first word, "" if absent)
    user = read_user_field(xlsx_files[0])

    script_content = build_script_full(reference_genome, all_pairs, user)

    out_script = f"cellranger-atac_{main_avt}.sh"
    with open(out_script, "w") as f:
        f.write(script_content)
    os.chmod(out_script, 0o755)

    print(f"Generated: {out_script}")
    if "FIXME_LOCATE_" in script_content:
        print("Note: some fastq paths could not be resolved — set the flagged "
              "path_* variables (FIXME_LOCATE_ markers) before running.")


if __name__ == "__main__":
    main()
