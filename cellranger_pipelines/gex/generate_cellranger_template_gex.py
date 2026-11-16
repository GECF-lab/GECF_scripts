#!/usr/bin/env python3
"""
generate_cellranger_template_gex.py

Usage:
    python3 generate_cellranger_template_gex.py file.xlsx [more.xlsx ...]
    python3 generate_cellranger_template_gex.py -n file.xlsx

GEX-only pipeline. Reads one or more Excel datasheets (the
`CR_count_datasheet_XXX_template.xlsx` layout) and emits a single
self-contained bash script (`cellranger-count_AVT0XXX.sh`) with one
`count_sample_<i>()` function per unique sample, plus the run-wide
dispatch / post-processing / archive block.
"""

import argparse
import os
import re
import sys
import time
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

# Ensure the cellranger_pipelines/ package root is on sys.path so ``common.intake``
# resolves. The generator lives at ``cellranger_pipelines/<assay>/`` so its parent
# dir is the package root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_CELLRANGER_PIPELINES_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir))
if _CELLRANGER_PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, _CELLRANGER_PIPELINES_ROOT)

from common.intake import (
    NGS_BASE_DIR_TEMPLATE,
    AVT_PATTERN,
    FASTQ_ID_PATTERN,
    LOCAL_CORES,        # cellranger --localcores (shared default "32")
    LOCAL_MEM,          # cellranger --localmem   (shared default "126")
    norm,
    build_merged_cell_lookup,
    get_cell_value,
    read_user_field,
    set_non_interactive,
    setup_readline_completer,
    ask_avt_run,
    ask_reference_genome,
    ask_sample_selection,
    choose_sample_dir,      # sample-level run-folder disambiguation (AVT0551)
    list_sample_dirs,       # ... list every candidate sample dir for an AVT
)

# Reference-genome root.  The genome folder inside each ``current_*``
# directory is never hard-coded here: each of them holds exactly one
# genome folder and that folder's name can change on the server
# (rename/replace), so the real path is resolved once at startup from
# the live filesystem (see the _resolve_reference_genome block below).
CELLRANGER_REFERENCES_ROOT = "/home/gecf/references/cellranger_references"
_HUMAN_CURRENT_DIR = os.path.join(CELLRANGER_REFERENCES_ROOT, "current_human")
_MOUSE_CURRENT_DIR = os.path.join(CELLRANGER_REFERENCES_ROOT, "current_mouse")
_HUMAN_MOUSE_CURRENT_DIR = os.path.join(CELLRANGER_REFERENCES_ROOT, "current_human_and_mouse")


def _is_version_check() -> bool:
    """True when the user is only asking for ``-v``/``--version`` so the
    import-time reference-genome scan can stay quiet (argparse hasn't run
    yet but ``sys.argv`` already carries the flag)."""
    return ("-v" in sys.argv or "--version" in sys.argv)


def _resolve_reference_genome(current_dir, label):
    """Return the genome path inside ``current_dir`` — its single subdirectory.

    Only *current_* is guaranteed to persist across the analysis servers;
    the genome folder it contains (refdata-gex-...) can be renamed or
    replaced without the generator being updated.  Non-recursive, so an
    inner ``old/`` tier never leaks in (only the current genome counts).
    """
    try:
        children = sorted(p for p in os.listdir(current_dir)
                          if os.path.isdir(os.path.join(current_dir, p)))
    except OSError:
        children = []
    if children:
        return os.path.join(current_dir, children[0])
    if not _is_version_check():
        print(f"WARNING: no genome folder found under {current_dir} "
              f"({label}); the reference prompt will fall back to it as-is")
    return current_dir


HUMAN_GENOME_PATH = _resolve_reference_genome(_HUMAN_CURRENT_DIR, "Human")
MOUSE_GENOME_PATH = _resolve_reference_genome(_MOUSE_CURRENT_DIR, "Mouse")
# Combined human+mouse reference (single reference usable for a mix of
# human and mouse samples) — offered as extra choice 3 in the genome
# prompt.  Same auto-detection as the others: current_human_and_mouse
# holds exactly one genome folder and its name is not hard-coded.
# The directory is a shared convention on the analysis servers; if it is
# absent on the caller's box, choice 3 is still listed and the generated
# script will fail its reference-genome existence check (same behaviour as
# choosing a custom path that does not exist yet).
HUMAN_MOUSE_GENOME_PATH = _resolve_reference_genome(
    _HUMAN_MOUSE_CURRENT_DIR, "Human+Mouse")
NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/count"

# "TO delete" folder where processed datasheets are parked for later
# cleanup. Two destinations, tried in order: Daniel's CR_pipelines tree
# first, the older CRcount_GEXonly layout as fallback.  The generated
# script uses the first folder that exists and creates neither (no
# forced mkdir).
DELETE_PRIMARY = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CR_pipelines_Daniel/TO_DELETE"
)
DELETE_FALLBACK = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CRcount_GEXonly/TO_DELETE"
)

# Root on the *analysis* server where cellranger count output + post-processing
# are written. On that host this is a real, mounted
# path; on the dev box it is not expected to exist. Never derive from $PWD —
# the script must land in the same place regardless of the caller's cwd.
SCRATCH_ROOT = "/mnt/data/CellRangerCountOutput"

# Provenance: version stamped into the generated .sh header
# (scriptversion=...) and into the script's own source as SCRIPT_VERSION.
# Generator version == generated-script version; the kept .sh in the
# archive tells you which one produced it.  Bump on any change to the
# generator.
SCRIPT_VERSION = "1.1.1"

# ---------------------------------------------------------------------------
# TSO% library-quality helper (cellranger_pipelines/common/TSO_percent.sh)
# The body of `tso_percent()` is inlined once per generated script so the
# analysis server needs no runtime dependency on this repo.
# ---------------------------------------------------------------------------
_TSO_HELPER_PATH = os.path.join(_CELLRANGER_PIPELINES_ROOT, "common", "TSO_percent.sh")


def _load_tso_helper_body() -> str:
    """Return the region of ``cellranger_pipelines/common/TSO_percent.sh`` up to the
    ``# %%EMBED-END%%`` marker (the ``tso_percent()`` function and its
    comments, WITHOUT the standalone-runner footer).

    Returns a warning string if the file is absent so the generated script
    still runs — the TSO% call will be skipped with a clear message.
    """
    if not os.path.isfile(_TSO_HELPER_PATH):
        return (
            "echo 'WARNING: cellranger_pipelines/common/TSO_percent.sh not found "
            f"(expected: {_TSO_HELPER_PATH}); tso_percent() will be a no-op.' >&2\n"
            "tso_percent() { :; }"
        )
    with open(_TSO_HELPER_PATH, encoding="utf-8") as fh:
        text = fh.read()
    marker = "# %%EMBED-END%%"
    pos = text.rfind(marker)  # last occurrence, in case a comment mentions it
    if pos == -1:
        return text.rstrip()
    return text[:pos].rstrip()


_TSO_HELPER_BODY = _load_tso_helper_body()


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


def find_fastq_path(avt_run, fastq_id, preferred_sample_dir=None):
    """Search under the run base dir for a directory named ``fastq_id``.

    ``preferred_sample_dir`` (a sample dir the user picked in the run-folder
    disambiguation menu, or ``None``) is honoured first: ``preferred/<id>``
    is checked before falling back to a full walk of the AVT tree.  This is
    what pins ``avidNNNNN`` into the correct run sub-folder when the same
    sample name exists under ``1_wrong_index_...`` AND
    ``2_corrected_indexes_...`` (the AVT0551 case).  Returns the full path
    if found, else None (so the generated script emits a FIXME_LOCATE_
    marker).
    """
    from common.intake import find_fastq_path_with_pref
    return find_fastq_path_with_pref(avt_run, fastq_id,
                                     preferred_sample_dir=preferred_sample_dir)


#: Per-sample knobs emitted at the top of the generated bash script.
#: Each entry is (var_suffix, default_value, short_comment).
#:  - default_value ""  -> user can fill, blank by default (flag omitted when blank)
#:  - default_value "true"/"false" -> both values accepted; user flips
#:  - default_value <number> -> user can change
#: Comments are user-facing: they describe what the user may type, so be
#: explicit about allowed values instead of referring to the flag.
_PER_SAMPLE_PARAMS = (
    ("include_introns", "true",
     "true = count intronic reads (default); false = exonic reads only"),
    ("create_bam",      "true",
     "true = emit aligned BAM to outs/ (default); false = faster + smaller output"),
    ("chemistry",       "",
     "blank = auto-detect (default); see header for explicit values"),
    ("expect_cells",    "",
     "expected number of cells; blank = cellranger estimates it (default)"),
    ("force_cells",     "",
     "keep exactly this many cells (overrides auto cell-calling); "
     "MUTUALLY EXCLUSIVE with expect_cells"),
    ("r1_length",       "",
     "trim Read 1 to this length in bases; blank = full read (default)"),
    ("r2_length",       "",
     "trim Read 2 to this length in bases; blank = full read (default)"),
    ("custom_options",  "",
     "extra flags passed through as-is; see header for examples"),
)


def build_script_full(reference_genome, pairs, user="", source_files=None,
                      chosen_sample_dirs=None, added_paths=None,
                      sample_order=False, number_prefix=True):
    """Build the complete `cellranger count` bash script (GEX).

    *sample_order* (default False): order of samples as they appear in the
    datasheet(s) — first-seen order preserved, NOT alphabetical sorting.
    Only pass True when the input is a SINGLE Excel file (with several
    files the cross-file order is ambiguous, so the alphabetical sort is
    kept as the unambiguous fallback).

    *number_prefix* (default True): zero-pad the sample index (01_, 02_,
    …) onto each per-sample output folder, keeping the folders in sample
    order even when the names are not alphabetical.  Pass False (CLI:
    --no-number-prefix) for bare sample names.

    Layout (top -> bottom):
      1. Header / constants (timestamp, reference genome, scriptversion,
         the source files list, the requester).
      2. RUN-WIDE PARAMETERS — one editable block at the top of the script:
         keepdata, max_parallel, localcores, localmem, scratch_root,
         server_folder_name, NGSruns_archive, and the per-sample FASTQ
         path_* variables (resolved at generation time).
      3. Per-sample parameter variables (one block per sample),
         editable in place; blank a value to drop the flag for that sample.
      4. Validation: ``--expect-cells`` and ``--force-cells`` are mutually
         exclusive per sample; if both are non-blank the script exits 1
         *before* any side effects.
      5. Side effects: working dir + archive dir creation, and a reference
         genome existence check (values come from the block in 2.).
      5a. The ``tso_percent()`` helper body inlined from
         ``cellranger_pipelines/common/TSO_percent.sh`` (self-contained .sh — no repo
         dependency on the analysis server).
      6. One `count_sample_<i>()` per sample: a `cellranger count`
         invocation using its `*_<idx>` variables, followed by that
         sample's own post-processing (SC_RNA_COUNTER_CS cleanup, tmp-dir
         completion check, TSO% metric to ``tso_report_parts/``, rsync of
         the sample dir + error log to the archive).
      6b/6c. DISPATCH (serial or a rolling pool of ``max_parallel``
         concurrent samples), the per-sample status report, and the single
         assembled ``tso_report.txt`` (concatenation of all parts, in
         sample order).
      7. ``extract_summary_metrics``: copy each sample's ``web_summary.html``,
         ``cell_types/web_summary_cell_types.html`` and
         ``metrics_summary.csv`` into ``raw_metrics_folder/`` (renamed
         ``<sample>_...``), and build a merged ``summary_metrics.csv`` with
         a Samplename column and every metric column (no filtering / reordering).
      7b. Copy this .sh into the working folder and the archive as the
         provenance record.
      7d. A single ``md5sums.txt`` over the entire run folder (replaces
         per-sample ``*_md5_checksums.txt``).
      8. Local workdir cleanup — unless ``keepdata=yes``.
      9. Park this run's source datasheets in the team "TO delete" folder
         (first existing of the two DELETE_* destinations; no mkdir -p).
    """
    lines = []

    avt_counts = Counter(avt for _, _, avt in pairs if avt is not None)
    main_avt = avt_counts.most_common(1)[0][0] if avt_counts else "MIXED_runs"

    samples_to_fastqs = defaultdict(list)
    for sn, fq, avt in pairs:
        samples_to_fastqs[sn].append((fq, avt))
    # Order: datasheet order (first-seen, insertion-preserving) or the
    # alphabetical sort fallback (multi-file input — see docstring).
    if sample_order:
        unique_samples = list(samples_to_fastqs.keys())
    else:
        unique_samples = sorted(samples_to_fastqs.keys())
    n = len(unique_samples)
    # Per-sample output folder name: 01_, 02_, … (ON by default) keeps the
    # folders in sample order even when the names are not alphabetical.
    # --no-number-prefix -> bare sample name (previous behaviour).
    dir_names = [f"{i:02d}_{s}" if number_prefix else s
                 for i, s in enumerate(unique_samples, 1)]

    # ---- 1) header / constants ----
    gen_stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    lines.append("#!/bin/bash")
    lines.append(f"# cellranger count — GEX-only ({main_avt})")
    lines.append("# Generated by generate_cellranger_template_gex.py "
                 f"({n} sample(s)).")
    lines.append(f"# Generated at: {gen_stamp}")
    _srcs = source_files or []
    if _srcs:
        lines.append("# Source file(s):")
        for _s in _srcs:
            lines.append(f"#   {os.path.basename(_s)}")
    else:
        lines.append("# Source file(s): (not provided)")
    lines.append(f"# Requester: {user or 'n/a'}")
    lines.append("")
    lines.append("TIMESTAMP=$(date +%Y%m%d_%H%M%S)")
    lines.append(f"reference_genome_used={reference_genome}")
    lines.append(
        # Provenance stamp in the kept .sh (copied into the run folder +
        # archive): which GENERATOR produced this script.  The generator
        # version == the generated-script version — edit SCRIPT_VERSION
        # above to bump it.
        "scriptversion=" + SCRIPT_VERSION
    )
    lines.append("")

    # ---- 2) run-wide parameters (defaults, editable at the top) ----
    # Requester tag, reduced to a safe path-segment (accented names like
    # "Cèsar" -> "Cesar").  Empty/sanitized-blank requester -> the whole
    # "_<user>" segment is omitted so we never emit a stray underscore.
    folder_user = folder_safe(user)
    user_tag = "_" + folder_user if folder_user else ""
    lines.append("# ============================================================")
    lines.append("# RUN-WIDE PARAMETERS  (defaults — edit here)")
    lines.append("# ----------------------------------------------------------------")
    lines.append('#   keepdata : "yes" = keep the local working folder after the run;')
    lines.append('#               "no"  = delete it (data already rsynced to the archive).')
    lines.append('keepdata="no"')
    lines.append('#   max_parallel : how many cellranger counts to run at the same time.')
    lines.append('#   1 = serial (one sample after another); N>1 = rolling pool of N jobs,')
    lines.append('#          the next sample starts as soon as a slot frees up.')
    lines.append('#   localcores : CPU cores per cellranger count (applies to ALL samples);')
    lines.append('#                blank = cellranger default.')
    lines.append('#   localmem   : RAM in GB per cellranger count (applies to ALL samples);')
    lines.append('#                blank = cellranger default.')
    lines.append('#   With max_parallel>1, total = max_parallel x (localcores cores, localmem GB),')
    lines.append("#   so all three must fit within the server's resources together.")
    lines.append('max_parallel=3')
    lines.append(f'localcores="{LOCAL_CORES}"')
    lines.append(f'localmem="{LOCAL_MEM}"')
    lines.append("#   scratch_root : root on the analysis server where cellranger count")
    lines.append("#                  output is written (fixed convention — change only if")
    lines.append("#                  the mount point differs on your server).")
    lines.append(f'scratch_root="{SCRATCH_ROOT}"')
    lines.append("#   server_folder_name : working + archive leaf folder name.")
    lines.append("#                  default is <AVT>[_user]_<timestamp>; change only if needed.")
    lines.append(f'server_folder_name="{main_avt}{user_tag}_$TIMESTAMP"')
    lines.append("#   NGSruns_archive  : central archive destination (fixed target).")
    lines.append(f'NGSruns_archive="{NGSruns_ARCHIVE_DIR}"')
    lines.append("#   FASTQ paths : one per (sample, FASTQ) pair, resolved by the generator")
    lines.append("#                  at generation time under the AVT run tree.")
    lines.append("#                  If a sample name exists in SEVERAL run sub-folders,")
    lines.append("#                  the generator asks which folder to use before")
    lines.append("#                  resolving — the FASTQ then lands inside the chosen")
    lines.append("#                  folder (e.g. corrected indexes vs other folders).")
    lines.append("#                  FIXME_LOCATE_ markers mean the generator could not find")
    lines.append("#                  the folder — set that path_* variable by hand.")

    any_missing = False
    per_sample_dir = {sn: chosen_sample_dirs.get(sn)
                      for sn in unique_samples} if chosen_sample_dirs else {}
    for idx, sn in enumerate(unique_samples, start=1):
        fastqs = samples_to_fastqs[sn]
        lines.append("")
        lines.append(f"# --- FASTQ paths — Sample {idx}: {sn} ({len(fastqs)} FASTQ(s)) ---")
        for fq, avt in fastqs:
            run_part = avt or main_avt
            base_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=run_part)
            # Per-sample chosen folder pins the run sub-folder (AVT0551's
            # corrected-vs-other case).  The helper below re-checks that the
            # chosen dir actually lives under this run's base before using
            # it, then falls back to a full walk otherwise — so we can pass
            # the sample's chosen dir unconditionally.
            chosen = per_sample_dir.get(sn)
            resolved = added_paths.get(fq) if added_paths else None
            if not resolved:
                resolved = find_fastq_path(run_part, fq, chosen)
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
    lines.append("")

    # ---- 3) per-sample parameter variables (editable) ----
    lines.append("# ============================================================")
    lines.append("# PER-SAMPLE PARAMETERS  (one block per sample)")
    lines.append("# ----------------------------------------------------------------")
    lines.append("#   * How to set:")
    lines.append("#       - blank a value ( value=\"\" )  -> flag DROPPED for that sample")
    lines.append("#       - number / text -> emitted as --flag=<value>")
    lines.append("#       - leave a default as-is to keep it")
    lines.append("#   * Allowed values (chemistry is the only one with a fixed list):")
    lines.append("#     include_introns: true | false")
    lines.append("#     create_bam     : true | false")
    lines.append("#     chemistry      : threeprime | fiveprime | SC3Pv1 | SC3Pv2 | SC3Pv3 |")
    lines.append("#                     SC3Pv4 | SC3Pv3HT | SC5P-PE | SC5P-PE-v3 | SC5P-R2 |")
    lines.append("#                     SC5P-R2-v3 | SC-FB | ARC-v1, or blank (auto-detect)")
    lines.append("#     expect_cells   : expected number of cells, or blank")
    lines.append("#     force_cells    : number of cells to keep, or blank")
    lines.append("#     r1_length      : length for Read 1 in bases, or blank")
    lines.append("#     r2_length      : length for Read 2 in bases, or blank")
    lines.append("#     custom_options : extra flags as-is, e.g. --nosecondary --localvmem=300")
    lines.append("#   * expect_cells and force_cells are MUTUALLY EXCLUSIVE.")
    lines.append("#     Both non-blank => this script exit 1s.")
    lines.append("#   * chemistry blank = cellranger auto-detects the assay (recommended).")
    lines.append("#     Set only if auto-detection fails. ARC (Multiome GEX) = ARC-v1.")
    lines.append("# ============================================================")
    for idx, sn in enumerate(unique_samples, start=1):
        lines.append("")
        lines.append(f"# --- Sample {idx}: {sn} ---")
        for suffix, default, comment in _PER_SAMPLE_PARAMS:
            base = f'{suffix}_{idx}="{default}"'
            pad = max(1, 60 - len(base))
            lines.append(f"{base}   {' ' * pad}# {comment}")
    lines.append("")

    # ---- 4) validation (MUST run before any side effects) ----
    lines.append("# ===== VALIDATION: --expect-cells vs --force-cells are mutually exclusive =====")
    lines.append("_vfail=0")
    for i in range(1, n + 1):
        lines.append(f'if [[ -n "${{expect_cells_{i}:-}}" && -n "${{force_cells_{i}:-}}" ]]; then')
        lines.append(f'    echo "ERROR: sample {i} has both expect_cells_{i} and force_cells_{i} set. Pick one; clear the other." >&2')
        lines.append("    _vfail=1")
        lines.append("fi")
    lines.append("if [[ ${_vfail} -ne 0 ]]; then")
    lines.append('    echo "Aborting: expect-cells / force-cells are mutually exclusive per sample. See errors above." >&2')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")

    # ---- 5) side effects (working dirs, archive, reference presence) ----
    # Values come from the RUN-WIDE PARAMETERS block at the top of the
    # script (edit there if a path differs on your server).
    lines.append("# Working directory: cellranger count output + post-processing (from RUN-WIDE PARAMETERS block)")
    lines.append('Server_folder="${scratch_root}/${server_folder_name}"')
    lines.append('mkdir -vp "${Server_folder}"')
    lines.append("# Parts dir: each sample (possibly running in parallel) writes its")
    lines.append("# TSO% report body here; the single centralized tso_report.txt is")
    lines.append("# assembled from them after DISPATCH (see TSO% REPORT section).")
    lines.append('mkdir -p "${Server_folder}/tso_report_parts"')
    lines.append("")
    lines.append('mkdir -p "${NGSruns_archive}/${server_folder_name}"')
    lines.append('NGSruns_folder=${NGSruns_archive}/${server_folder_name}')
    SFR = "${Server_folder}"
    lines.append("")
    lines.append("# Check that the reference genome directory exists")
    lines.append('if [[ ! -d "${reference_genome_used}" ]]; then')
    lines.append('    echo "Reference genome ${reference_genome_used} does not exist. Aborting." >&2')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append("")

    # ---- 5a) tso_percent() helper body — inlined from
    #          cellranger_pipelines/common/TSO_percent.sh so the generated .sh is
    #          self-contained on the analysis server (no repo dependency at
    #          run time).  Edit the .sh and regenerate to change behaviour ----
    lines.append("")
    lines.append("# ==============================================================================")
    lines.append("# tso_percent helper (from cellranger_pipelines/common/TSO_percent.sh, embedded verbatim)")
    lines.append("# Computes Read-2 percentage metrics: intact TSO, intact SMART, polyA/polyT tails.")
    lines.append("# ==============================================================================")
    lines.extend(_TSO_HELPER_BODY.splitlines())

    # ---- 6) one cellranger count per sample, each with its old-CRcount
    #          post-processing (SC tree cleanup, completion check, per-sample
    #          deploy, md5 checksums) ----
    _CONT = " \\"    # bash line-continuation used between flags
    lines.append("# ============================================================")
    lines.append("# PER-SAMPLE FUNCTIONS  (one count_sample_<i>() per sample)")
    lines.append("#   Each function is the old inline per-sample block:")
    lines.append("#     * cellranger count - flags from the *_<idx> variables above;")
    lines.append("#       a blank value drops that flag for that sample only")
    lines.append("#     * completion check: find <sample> -type d -name tmp")
    lines.append("#     * post-processing: rm -r <sample>/SC_RNA_COUNTER_CS")
    lines.append("#     * rsync -avh <sample>/ + logs -> the archive")
    lines.append("#     * TSO% library-quality metric (Read-2)")
    lines.append("#   (the single whole-folder md5sums.txt is written ONCE, after")
    lines.append("#   DISPATCH — see MD5 CHECKSUMS (WHOLE FOLDER), not per sample)")
    lines.append("#   The DISPATCH section below runs the functions either serially")
    lines.append("#   (max_parallel=1) or in a rolling pool of max_parallel jobs.")
    lines.append("# ============================================================")
    lines.append("")
    for idx, sn in enumerate(unique_samples, start=1):
        fastqs = samples_to_fastqs[sn]
        numbered_dir = dir_names[idx - 1]
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
            '    --output-dir="' + SFR + '/' + numbered_dir + '"',
            # CPU/memory caps are run-wide (localcores / localmem), not per-sample
            '    ${localcores:+--localcores=$localcores}',
            '    ${localmem:+--localmem=$localmem}',
        ]
        # Per-sample conditional knobs. Blank value -> flag omitted for this
        # sample only. custom_options is deliberately UNQUOTED on the last
        # line so it word-splits into multiple --flags.
        for suffix, _default, _comment in _PER_SAMPLE_PARAMS:
            v = f"{suffix}_{idx}"
            if suffix == "custom_options":
                flag_lines.append(f"    ${{{v}}}")
            else:
                flag = "--" + suffix.replace("_", "-")
                flag_lines.append(f"    ${{{v}:+{flag}=${v}}}")

        lines.append("")
        lines.append(f"# --- Run sample {idx}: {sn} ---")
        lines.append("function count_sample_" + str(idx) + "() {")
        # Emit the command for execution, line-continuation formatted,
        # from the SAME flag_lines list (so what runs is exactly what
        # appears in the .sh).  The .sh itself (copied into the run
        # folder + archive by KEEP A COPY OF THIS SCRIPT) is the
        # provenance record; cellranger's own command.log in
        # <sample>/outs/logs/ is the execution record — no separate
        # per-sample script_log needed.
        for i, fl in enumerate(flag_lines):
            lines.append(fl + (_CONT if i < len(flag_lines) - 1 else ""))

        # Per-sample post-processing.
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
            'rsync -avh "' + SFR + '/' + numbered_dir + '_error_log.txt" "${NGSruns_folder}" 2>/dev/null'
        )

        # ---- TSO% library-quality metric ------------------------------------
        # Compute intact-TSO / intact-SMART / polyA / polyT fractions on a
        # DOWNSAMPLING-sized sample of Read-2 sequences.  The R2 fastq is
        # auto-discovered from the resolved fastq path:
        #   * if it is already a file ending in *_R2_* — use it directly
        #   * otherwise take the first *_R2_* file in the dir (one avid
        #     folder holds exactly one R1 + one R2, so first == only)
        # The report body is echoed to the console AND written to a
        # PER-SAMPLE PART file under
        # ${Server_folder}/tso_report_parts/<NN_sample>.txt.  One centralized
        # tso_report.txt is assembled from all parts AFTER DISPATCH (samples
        # may run in parallel, so each part has exactly one writer);
        # see the "TSO% REPORT (SINGLE FILE)" section.
        if fastqs:
            _fq0_var = "path_%d_%s" % (idx, fastqs[0][0])
            # Part name stays 01_sn (matches the assembler below) regardless
            # of the per-sample folder naming.
            tags = f"{idx:02d}_{sn}"
            lines.append("")
            lines.append("# TSO% library-quality metric (Read-2 only) — intact TSO/SMART + poly tails")
            lines.append(f'_fq="${{{_fq0_var}}}"')
            lines.append('_TAGS="%s"' % tags)
            lines.append('if [ ! -f "${_fq}" ]; then')
            lines.append('    for _c in "${_fq}"/*_R2_*.fastq.gz "${_fq}"/*_R2.fastq.gz; do')
            lines.append('        [ -f "${_c}" ] && { _fq="${_c}"; break; }')
            lines.append("    done")
            lines.append("fi")
            lines.append('if [ -f "${_fq}" ]; then')
            lines.append('    _tso_part="' + SFR + '/tso_report_parts/' + tags + '.txt"')
            lines.append('    tso_percent "${_fq}" "${_TAGS}" | tee "${_tso_part}"')
            lines.append("else")
            lines.append(f'    echo "tso: no Read-2 fastq found near ${{{_fq0_var}}}; TSO% step skipped." >&2')
            lines.append("fi")
        lines.append("}")
        lines.append("")

    # ---- 6b) DISPATCH + per-sample status report ----
    _N = len(unique_samples)
    _DASH = "# --------------------------------------------------------"
    lines.append(_DASH)
    lines.append("# DISPATCH  (run count_sample_1 .. count_sample_" + str(_N) + ")")
    lines.append("#   max_parallel=1 : serial, one sample after another (previous behaviour)")
    lines.append("#   max_parallel=N : rolling pool of N concurrent cellranger counts; the")
    lines.append("#                    next sample starts as soon as one slot frees up")
    lines.append("#   Resource caps are run-wide (localcores / localmem), so")
    lines.append("#   max_parallel=3 @ defaults (32 cores / 126 GB) tops out at 96 cores /")
    lines.append("#   378 GB - within the 128-core / ~500 GB analysis box.")
    lines.append(_DASH)
    lines.append("NSAMPLES=" + str(_N))
    lines.append('case "${max_parallel}" in')
    lines.append("    ''|*[!0-9]*)")
    lines.append('        echo "ERROR: max_parallel must be a positive integer" >&2')
    lines.append("        exit 1")
    lines.append("    ;;")
    lines.append("esac")
    lines.append('if [ "${max_parallel}" -lt 1 ]; then')
    lines.append('    echo "ERROR: max_parallel must be >= 1" >&2')
    lines.append("    exit 1")
    lines.append("fi")
    lines.append('# Cap max_parallel at NSAMPLES: there is no point requesting more')
    lines.append('# concurrent cellranger counts than the number of samples, so a')
    lines.append('# small run (e.g. 1 sample) never hard-fails on a larger default.')
    lines.append('if [ "${max_parallel}" -gt "${NSAMPLES}" ]; then')
    lines.append('    echo "Note: ${NSAMPLES} sample(s) only — capping max_parallel ${max_parallel} -> ${NSAMPLES}." >&2')
    lines.append("    max_parallel=${NSAMPLES}")
    lines.append("fi")
    lines.append("")
    lines.append('if [ "${max_parallel}" -eq 1 ]; then')
    lines.append("    # --- SERIAL ---")
    lines.append("    for ((_s = 1; _s <= ${NSAMPLES}; _s++)); do")
    lines.append("        count_sample_${_s}")
    lines.append("    done")
    lines.append("else")
    lines.append("    # --- ROLLING POOL OF ${max_parallel} CONCURRENT JOBS ---")
    lines.append("    # Slot arrays: _pids[_i] holds the background PID, _idxs[_i] the sample.")
    lines.append("    declare -a _pids=() _idxs=()")
    lines.append("    for ((_i = 0; _i < ${max_parallel}; _i++)); do")
    lines.append('        _pids[$_i]=""')
    lines.append("        _idxs[$_i]=0")
    lines.append("    done")
    lines.append("    _next=1")
    lines.append("    while :; do")
    lines.append("        # 1) reap finished slots")
    lines.append("        for ((_i = 0; _i < ${max_parallel}; _i++)); do")
    lines.append('            if [ -n "${_pids[$_i]}" ] && ! kill -0 "${_pids[$_i]}" 2>/dev/null; then')
    lines.append('                wait "${_pids[$_i]}"')
    lines.append("                _rv=$?")
    lines.append('                echo "[REAP] sample #${_idxs[$_i]} finished (rc=${_rv})"')
    lines.append('                _pids[$_i]=""')
    lines.append("            fi")
    lines.append("        done")
    lines.append("        # 2) fill free slots with the next sample (if any left)")
    lines.append("        _launched=0")
    lines.append("        for ((_i = 0; _i < ${max_parallel}; _i++)); do")
    lines.append('            if [ -z "${_pids[$_i]}" ] && [ "${_next}" -le "${NSAMPLES}" ]; then')
    lines.append("                count_sample_${_next} &")
    lines.append("                _pids[$_i]=$!")
    lines.append("                _idxs[$_i]=${_next}")
    lines.append("                _next=$(( _next + 1 ))")
    lines.append("                _launched=1")
    lines.append("            fi")
    lines.append("        done")
    lines.append("        # 3) exit when all samples launched AND all jobs reaped")
    lines.append("        _busy=0")
    lines.append("        for ((_i = 0; _i < ${max_parallel}; _i++)); do")
    lines.append('            [ -n "${_pids[$_i]}" ] && _busy=1')
    lines.append("        done")
    lines.append('        if [ "${_next}" -gt "${NSAMPLES}" ] && [ "${_busy}" -eq 0 ]; then')
    lines.append("            break")
    lines.append("        fi")
    lines.append("        # all slots busy and more samples remain: give them a moment")
    lines.append('        if [ "${_launched}" -eq 0 ] && [ "${_next}" -le "${NSAMPLES}" ]; then')
    lines.append("            sleep 1")
    lines.append("        fi")
    lines.append("    done")
    lines.append("fi")
    lines.append("")
    # Per-sample status report - each sample's own completion markers:
    #   *_error_log.txt non-empty    -> cellranger left a "tmp" dir behind on failure
    #   no outs/metrics_summary.csv  -> cellranger count never completed
    lines.append("# Per-sample status report (from each sample's own completion markers)")
    lines.append("_overall_rc=0")
    _dirs_bash = " ".join('"' + d + '"' for d in dir_names)
    lines.append("for _d in " + _dirs_bash + "; do")
    lines.append('    if [ -s "${Server_folder}/${_d}_error_log.txt" ]; then')
    lines.append('        echo "[FAIL] ${_d}  (see ${Server_folder}/${_d}_error_log.txt)" >&2')
    lines.append("        _overall_rc=1")
    lines.append('    elif [ ! -f "${Server_folder}/${_d}/outs/metrics_summary.csv" ]; then')
    lines.append('        echo "[FAIL] ${_d}  (no outs/metrics_summary.csv; cellranger count did not complete)" >&2')
    lines.append("        _overall_rc=1")
    lines.append("    else")
    lines.append('        echo "[ OK ] ${_d}"')
    lines.append("    fi")
    lines.append("done")
    lines.append("")

    # ---- 6c) TSO% — assemble the ONE centralized report -------------------
    # Samples (serially or in a parallel pool) each wrote their part under
    # tso_report_parts/; now concatenate them, in sample order, into a single
    # tso_report.txt and rsync that one file to the archive.  The working
    # copy is removed (the archive copy is the deliverable; with
    # keepdata=no the whole working dir is deleted below anyway).
    lines.append("# ------------------------------------ TSO% REPORT (SINGLE FILE) ------------------------------------ #")
    lines.append("#   Centralized TSO% report: one tso_report.txt gathering every sample's")
    lines.append("#   per-sample section (from tso_report_parts/), in sample order.")
    lines.append('_tso_out="' + SFR + '/tso_report.txt"')
    lines.append("{")
    lines.append('    echo "TSO% library-quality report — one section per sample"')
    lines.append('    echo "Generated: $(date)"')
    lines.append('    echo "Server folder: ' + SFR + '"')
    lines.append('    echo "Samples: ' + str(n) + '"')
    lines.append('    echo ""')
    for idx_t, sn_t in enumerate(unique_samples, start=1):
        part = f"{idx_t:02d}_{sn_t}"
        lines.append(f'    if [ -f "{SFR}/tso_report_parts/{part}.txt" ]; then')
        lines.append(f'        cat "{SFR}/tso_report_parts/{part}.txt"')
        lines.append("    else")
        lines.append(f'        echo "[{part}] (no TSO% report — no Read-2 fastq found or tso_percent skipped)" >&2')
        lines.append(f'        echo "[{part}] (no report: Read-2 fastq not found or step skipped)"')
        lines.append("    fi")
        lines.append('    echo ""')
    lines.append('} > "${_tso_out}"')
    lines.append('rm -rf "' + SFR + '/tso_report_parts"')
    lines.append('echo "TSO% report (single file): ${_tso_out}"')
    lines.append('rsync -avh "${_tso_out}" "${NGSruns_folder}" 2>/dev/null || true')
    lines.append("")

    # ---- 7) extract_summary_metrics ----
    #   copy per-sample web_summary.html / web_summary_cell_types.html /
    #   metrics_summary.csv into ${Server_folder}/raw_metrics_folder renamed
    #   as <sample>_... ; strip comma thousands-separators inside quoted
    #   fields (awk -F'"'); strip "N%" -> "N" (sed); build summary_metrics.csv
    #   with Samplename column and all metric columns (no filtering / reordering).
    SFR = "${Server_folder}"
    # literal sample-name list (no dynamic expansion in the emitted bash)
    # — must match the per-sample dir names used in the cellranger count
    # section above (01_, 02_, … or bare names with --no-number-prefix).
    _sample_list = " ".join(dir_names)
    lines.append("# ------------------------------------ EXTRACT SUMMARY METRICS ------------------------------------ #")
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
    lines.append("    # strip comma thousands-separators inside quoted fields")
    lines.append(r'''    awk -F '"' 'BEGIN {OFS = "\""} { for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' "${dst}" > "${tmp}"''')
    lines.append("    # strip % symbols")
    lines.append('    sed -E "s/([0-9]+)%/\\1/g" "${tmp}" > "${dst}"')
    lines.append('    rm -f "${tmp}"')
    lines.append("done")
    lines.append("")
    lines.append("# Build summary_metrics.csv with Samplename column")
    lines.append("first_csv_processed=false")
    lines.append('SM_TMP="' + SFR + '/summary_metrics_tmp.csv"')
    lines.append('SM_OUT="' + SFR + '/summary_metrics.csv"')
    lines.append("for _sd in " + _sample_list + "; do")
    lines.append('    file="${metrics_folder}/${_sd}_metrics_summary.csv"')
    lines.append('    if [ ! -f "${file}" ]; then')
    lines.append('        echo "WARNING: ${file} not found, skipping" >&2')
    lines.append('        continue')
    lines.append('    fi')
    lines.append('    if [ "$first_csv_processed" = false ]; then')
    lines.append('        echo -n "Samplename," > "${SM_TMP}"')
    lines.append('        sed -n "1p" "${file}" >> "${SM_TMP}"')
    lines.append('        first_csv_processed=true')
    lines.append('    fi')
    lines.append('    echo -n "${_sd}," >> "${SM_TMP}"')
    lines.append('    sed -n "2p" "${file}" >> "${SM_TMP}"')
    lines.append("done")
    # Promote the full merged table. SM_TMP already holds Samplename plus
    # every metric column from the per-sample summaries (10X natural order),
    # so no column is filtered or reordered.
    lines.append("# Promote the full merged table (Samplename + all metric columns, none filtered).")
    lines.append('if [ -f "${SM_TMP}" ]; then')
    lines.append('    mv "${SM_TMP}" "${SM_OUT}"')
    lines.append('else')
    lines.append('    echo "WARNING: no per-sample metrics_summary.csv found; summary_metrics.csv not created" >&2')
    lines.append('fi')
    lines.append('rsync -avh "${metrics_folder}" "${NGSruns_folder}"')
    lines.append('rsync -avh "${SM_OUT}" "${NGSruns_folder}"')
    lines.append("")

    # ---- 7b) copy the generating .sh into the project folders ----
    # The user runs this script from their working directory (e.g.
    # `./cellranger-count_AVT0571.sh`), so $0 points at it.  Save a copy in
    # the local project folder AND the central archive, which both live
    # under paths that outlive the caller's cwd.
    lines.append("# ------------------------------------ KEEP A COPY OF THIS SCRIPT ------------------------------------ #")
    lines.append("#   Copy the .sh used to generate this project into the project folders")
    lines.append("#   (local working dir + central archive) so the exact pipeline is")
    lines.append("#   preserved with the results.")
    lines.append('if [ -n "${BASH_SOURCE:-}" ]; then _self="${BASH_SOURCE[0]}"')
    lines.append('elif [ -n "${0:-}" ]; then _self="${0}"; fi')
    lines.append('if [ -f "${_self:-/nonexistent}" ]; then')
    lines.append('    cp "${_self}" "${Server_folder}/$(basename "${_self}")" 2>/dev/null || true')
    lines.append('    cp "${_self}" "${NGSruns_folder}/$(basename "${_self}")" 2>/dev/null || true')
    lines.append("fi")
    lines.append("")
    lines.append("")
    # ---- 7d) single md5 checksums file for the ENTIRE run folder -------
    # One md5sums.txt covering every file under ${Server_folder} (sample dirs,
    # logs, TSO report, summary metrics, the .sh copy — the whole tree), with
    # paths relative to the folder.  Computed after every artifact is in
    # place and rsync'd to the archive.  Placed after KEEP A COPY OF THIS
    # SCRIPT so the .sh is included too; before LOCAL WORKDIR CLEANUP so
    # the local files still exist for the checksum to be computed from.
    lines.append("# ------------------------------------ MD5 CHECKSUMS (WHOLE FOLDER) ------------------------------------ #")
    lines.append("#   One md5sums.txt covering the ENTIRE run folder (every file,")
    lines.append("#   relative ./ paths).  Replace the old per-sample *_md5_checksums.txt files.")
    lines.append('_md5_out="' + SFR + '/md5sums.txt"')
    lines.append('( cd "' + SFR + '" && find . -type f ! -name "md5sums.txt" | sort | xargs md5sum > "md5sums.txt" )')
    lines.append('echo "MD5 checksums (whole folder): ${_md5_out}"')
    lines.append('rsync -avh "${_md5_out}" "${NGSruns_folder}" 2>/dev/null || true')
    lines.append("")
    lines.append("")
    # ---- 8) keepdata cleanup ----
    lines.append("# ------------------------------------ LOCAL WORKDIR CLEANUP ------------------------------------ #")
    lines.append('#   After deploy to the archive, the local working folder is deleted')
    lines.append('#   UNLESS keepdata="yes".')
    lines.append('if [[ "${keepdata}" != "yes" ]]; then')
    lines.append('    rm -r "${Server_folder}" 2>/dev/null || true')
    lines.append("fi")
    lines.append("")

    # ---- 9) park processed datasheets in the team "TO delete" folder ----
    # Shared destination instead of a per-script ${mypwd}.  No mkdir -p:
    # the target is created only as a side effect of the mv, and we fall
    # back to the CRcount_GEXonly destination if Daniel's CR_pipelines
    # tree is absent.
    lines.append("# ------------------------------------ DATASHEET CLEANUP ------------------------------------ #")
    lines.append("#   Park the Excel datasheets USED TO GENERATE THIS SCRIPT here for later")
    lines.append("#   deletion — only these exact files; any other .xlsx in the working dir")
    lines.append("#   is left untouched.  Absolute paths (captured at generation time), so")
    lines.append("#   this works no matter where the script is launched from.")
    lines.append(f"#   Primary: {DELETE_PRIMARY}")
    lines.append(f"#   Fallback: {DELETE_FALLBACK}")
    lines.append(f'if [ -d "{DELETE_PRIMARY}" ]; then _delfolder="{DELETE_PRIMARY}"')
    lines.append(f'elif [ -d "{DELETE_FALLBACK}" ]; then _delfolder="{DELETE_FALLBACK}"')
    lines.append("fi")
    lines.append('if [ -n "${_delfolder:-}" ]; then')
    _srcs = []
    for _s in (source_files or []):
        _s_abs = os.path.abspath(_s)
        if _s_abs not in _srcs:
            _srcs.append(_s_abs)
    for _s in _srcs:
        # `${_delfolder}` must survive the f-string as literal bash:
        # doubled braces `${{_delfolder}}` -> `${_delfolder}`.
        lines.append(f'    mv "{_s}" "${{_delfolder}}" 2>/dev/null || true')
    lines.append("fi")
    lines.append("")

    lines.append('exit "${_overall_rc}"')

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parses one or more GEX-only Excel datasheets and generates a "
            "cellranger-count_AVT0XXX.sh script with one "
            "`count_sample_<i>()` function per sample, plus the run-wide "
            "dispatch / post-processing / archive block."
        ),
        epilog="Example:\n  python3 %(prog)s samples.xlsx\n  python3 %(prog)s run1.xlsx run2.xlsx",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"generate_cellranger_template_gex.py {SCRIPT_VERSION}",
        help="Print the generator version and exit.",
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
    parser.add_argument(
        "--no-number-prefix",
        action="store_true",
        default=False,
        help=(
            "Do not number the per-sample output folders (ON by default): "
            "the folder becomes '<sample>' instead of '01_<sample>'. "
            "Numbering keeps the folders in sample order even when the "
            "names are not alphabetical."
        ),
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

        print(f"Reading '{os.path.basename(xlsx_file)}' …", flush=True)
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

    main_avt = (Counter(avt for _, _, avt in all_pairs if avt is not None)
                .most_common(1)[0][0] if any(a for _, _, a in all_pairs) else "MIXED_runs")

    # Pre-build sample -> fastqs (same shape `build_script_full` uses
    # internally) so the selection prompt can display fastq counts for each
    # sample; used below to re-filter the pair list to only the kept samples.
    _prefilter_s2f = defaultdict(list)
    for sn, fq, avt in all_pairs:
        _prefilter_s2f[sn].append((fq, avt))

    # ── Sample selection (keep/drop + add sheet-less samples) ──────────────
    # Mirrors the reference-genome prompt pattern: interactive by default,
    # skipped under `-n` / non-tty, in which case every detected sample is
    # kept and none are added (identical to the old behaviour).  Dropped
    # samples are removed from the generated script entirely — no per-sample
    # block, no cellranger invocation.  Added samples are each specified by
    # the user as a FULL path to an avid folder (whose basename matches
    # ``avid\d{4,6}``); these paths are used verbatim in the generated
    # script and their AVT run is left ``None`` (unused, the path is
    # self-referential).
    kept_names, added_raw = ask_sample_selection(
        _prefilter_s2f,
        default_root="/mnt/PTEGraw/AA/NGSruns/AVT/AV240401",
    )
    if dropped_names := [sn for sn in _prefilter_s2f if sn not in kept_names]:
        print(f"  Dropped from script: {', '.join(dropped_names)}")
    new_pairs = [(sn, fq, avt) for sn in kept_names
                 for fq, avt in _prefilter_s2f[sn]]
    # Samples ADDED off-prompt carry a list of FULL avid-folder paths
    # (self-referential — no AVT-tree walk needed).  The pair's
    # fastq_id is the folder's basename (e.g. "avid13144"); the
    # generator consumes the full path directly when it is known.
    added_paths = {}   # {fastq_id: full_path}
    for sn, info in added_raw.items():
        for p in info:
            added_paths[os.path.basename(p)] = p
        new_pairs.extend((sn, os.path.basename(p), None) for p in info)
    if not new_pairs:
        sys.exit("Error: no samples to add to the script.")
    all_pairs = new_pairs

    # ── Run-folder disambiguation (sample-level pick) ───────────────────────
    # A sample name can exist in SEVERAL run sub-folders under the AVT tree
    # (AVT0551: same 384 avid folders under both 1_wrong_index_… and
    # 2_corrected_indexes_…).  A bare avid id cannot pin down a location, so
    # before we resolve any FASTQ we pin the SAMPLE's folder here.  We only
    # ask when there really is more than one candidate for that AVT; a
    # single candidate resolves silently (the 99% case).  `choose_sample_dir`
    # also handles non-interactive mode (first candidate + warning).
    # Keyed by (sample, avt) because a sample could in principle appear
    # under two AVTs; for the usual one-AVT run every entry shares the
    # same avt.
    #
    # Samples the user ADDED on the command prompt (not from the Excel
    # sheet) carry a user-supplied FULL path to each avid folder — that
    # is the location and no disambiguation prompt is needed.
    selected_dirs = {}
    for sn in dict.fromkeys(sn for sn, _, _ in all_pairs):
        if sn in added_raw:
            continue  # user typed the full avid path(s); no anchor to pick
        # Each (sample, fq, avt) row's avt may be None; resolve to main_avt.
        row_avts = {avt for sn2, _, avt in all_pairs if sn2 == sn and avt}
        avt_use = next(iter(row_avts)) if row_avts else main_avt
        cands = list_sample_dirs(avt_use, sn)
        if len(cands) > 1:
            chosen = choose_sample_dir(sn, cands)
            if chosen:
                selected_dirs[(sn, avt_use)] = chosen

    # Build the sample_name -> sample_dir map the builder consumes.
    # Added samples are NOT included here — their paths are passed
    # separately via `added_paths` and used verbatim.
    chosen_for_build = None
    if selected_dirs:
        flat = {}
        for (sn, _), dirp in selected_dirs.items():
            flat.setdefault(sn, dirp)
        chosen_for_build = flat
        n_asked = len(flat)
        if n_asked:
            print(f"  Run-folder disambiguation picked for {n_asked} sample(s).")

    # ── Final sample summary (shown before the genome prompt) ──────────────
    # Recap of everything that will appear in the generated script: the
    # kept samples (with their FASTQs) plus any sheet-less samples the user
    # added by full path.  Always prints, in interactive and non-interactive
    # mode alike.
    _final_s2f = defaultdict(list)
    for sn, fq, _ in all_pairs:
        _final_s2f[sn].append(fq)
    _n_fq = len(all_pairs)
    _kept = [sn for sn in _final_s2f if sn not in added_raw]
    _added = [sn for sn in dict.fromkeys(added_raw.keys())]
    print()
    print("  ── Final samples (genome selection next) ──")
    for _sn in _kept:
        print(f"    keep   {_sn}  [{', '.join(_final_s2f[_sn])}]")
    for _sn in _added:
        print(f"    added  {_sn}  [{', '.join(_final_s2f[_sn])}]  (added by path)")
    print(f"    Total: {len(_kept) + len(_added)} sample(s), "
          f"{_n_fq} FASTQ ID(s)")
    print()

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=HUMAN_GENOME_PATH,
        human_label="Human", mouse_label="Mouse",
        extra_choices=[("Human + mouse (combined)", HUMAN_MOUSE_GENOME_PATH)],
    )
    print("Reading the requester (User) field and generating the script …", flush=True)
    user = read_user_field(xlsx_files[0])

    script_content = build_script_full(reference_genome, all_pairs, user,
                                       source_files=xlsx_files,
                                       chosen_sample_dirs=chosen_for_build,
                                       added_paths=added_paths,
                                       # Datasheet order for the 01_/02_/…
                                       # numbering: only when a SINGLE Excel
                                       # file drives the generation (with
                                       # several files the cross-file order
                                       # is ambiguous -> alphabetical).
                                       sample_order=(len(xlsx_files) == 1),
                                       # Per-sample folder numbering (01_,
                                       # 02_, …) is ON by default;
                                       # --no-number-prefix disables it.
                                       number_prefix=not args.no_number_prefix)
    out_script = f"cellranger-count_{main_avt}.sh"
    with open(out_script, "w") as f:
        f.write(script_content)
    os.chmod(out_script, 0o755)

    print(f"Generated: {out_script}")
    # Only flag when a path was genuinely left unresolved.  Real unresolved
    # markers are always emitted in *path position* — a trailing
    # "/FIXME_LOCATE_<id>" segment (e.g. "…/<run>/FIXME_LOCATE_<id>").  The
    # header help text uses the same token as a noun ("FIXME_LOCATE_ markers
    # mean …"), so scanning for the bare token fires even when every path
    # resolved (which is why this note used to show up on clean runs).
    n_unresolved = script_content.count("/FIXME_LOCATE_")
    if n_unresolved:
        print(
            f"Note: {n_unresolved} FASTQ path(s) could not be resolved — "
            "set the flagged path_* variable(s) (FIXME_LOCATE_ markers) "
            "before running."
        )


if __name__ == "__main__":
    main()
