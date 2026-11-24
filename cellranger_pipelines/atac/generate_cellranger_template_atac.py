#!/usr/bin/env python3
"""
generate_cellranger_template_atac.py

Usage:
    python3 generate_cellranger_template_atac.py file.xlsx [more.xlsx ...]
    python3 generate_cellranger_template_atac.py -n file.xlsx

Generate a cellranger-atac count script from an AVT submission Excel file.
The generator reads the sheet named after ``run_details`` (default:
"run details"), detects the AVT run, asks which sample(s) to keep (and lets
you add sample(s) that are NOT in the sheet), resolves each FASTQ leaf
directory under the AVT run tree, and writes a self-contained
``cellranger-atac_<AVT>.sh`` next to the generator.  When the inputs span
several AVT runs, every run is listed in the name
(``cellranger-atac_AVT0001_AVT0002.sh``; the working/archive folder
``AVT0001_AVT0002_user_$TIMESTAMP`` carries the same tag) — a single-run
input keeps the plain name.

Excel layout (one row per line):
    Column A = sample name (required)
    Column B = fastq ID (optional when a FASTQ path is given; e.g. "avid12543").
               A row that gives the full FASTQ path but no fastq ID also works.
    Column C = full FASTQ path (optional; the row must give at least a
               fastq ID or a FASTQ path — a row with neither is skipped)
    Column D = run name, header 'AVT run' (optional; an AVT/MSH/MSQ/NVQ
               code; when blank the run detected from the sheet/file name
               is used)

Sheets:
    The first sheet that carries the data rows (``run_details``-style header)
    wins; a ``User_x1x_copy`` sheet - if present - is always preferred.
    Sheets whose name contains 'spike', 'several user', 'guideline', or
'index table' (GECF legend/reference sheets) are never picked.

Run name detection (AVT/MSH/MSQ/NVQ codes):
    * Excel file name:  e.g. "012534_AVT0543_ChangXu.xlsx" -> AVT0543
    * run_details sheet: value in a cell (column A-D) of the first 10 rows

FASTQ leaf resolution (AVT0551 case, shared with GEX/ARC generators):
    * ``list_sample_dirs()`` lists every sample dir under the AVT tree
      (one per demux attempt in a ``fastq_bases2fastq/`` style layout).
    * A single candidate resolves silently (the 99% case); more than one
      prompts via ``choose_sample_dir()`` (auto-pick + warning under -n).
    * ``find_fastq_path_with_pref()`` then resolves the avid id inside the
      chosen dir first and falls back to the classic full-tree walk
      (``find_fastq_path()``).
    * A row that names no AVT inherits the run named by its OWN file, so a
      worksheet mixing two runs resolves each FASTQ under its own run tree.
      Rows whose file names no run at all are explicitly prompted for one
      (GEX 'needs_avt' parity, shared ``resolve_unassigned_avt()``) instead
      of silently inheriting the representative run.

Reference genome:
    Like the GEX/ARC generators, the genome folder name is NOT hard-coded.
    Each ``current_<species>_arc_ATAC`` directory under
    /home/gecf/references/cellranger_references is read at startup and its
    single genome folder (refdata-cellranger-arc-...) is used — rename or
    replace the folder on the server and the generator follows.  A custom
    path can still be typed at the reference prompt.

Generated script (cellranger-atac_<AVT>.sh; multi-run inputs:
cellranger-atac_AVT0001_AVT0002.sh), top to bottom:
    1.  Header / constants: requester, timestamp, reference, AVT run,
        scriptversion stamp, editable RUN-WIDE PARAMETERS (keepdata,
        max_parallel, localcores, localmem, scratch_root,
        server_folder_name, NGSruns_archive).
    2.  FASTQ PATHS: one variable per sample FASTQ, resolved at generation
        time (empty + FIXME_LOCATE marker when unresolved).
    3.  PER-SAMPLE PARAMETERS: custom_options, one editable block per
        sample (blank a value to drop the flag). CPU/memory caps are
        run-wide (localcores / localmem in RUN-WIDE PARAMETERS), matching
        the GEX/ARC generators — they are NOT per-sample.
    4.  One ``atac_sample_<i>()`` function per sample: runs
        ``cellranger-atac count`` (flags from its parameter block), checks
        for leftover tmp folders, deletes the non-essential SC_ATAC_COUNTER_CS
        scratch tree, harvests the web summary into ``<AVT>_Summaries/``.
    5.  DISPATCH: serial (max_parallel=1) or a rolling pool of max_parallel
        concurrent counts (default 3).
    6.  Per-sample STATUS SWEEP: [OK]/[FAIL] per sample, ``_overall_rc``.
    7.  MERGED SUMMARY: ``all_summary.csv`` (Samplename column first).
    8.  rsync of the working folder to the central archive.
    9.  WHOLE-FOLDER ``md5sums.txt`` + rsync of the single checksums file.
    10. A copy of the .sh used to generate the project is kept in the work
        folder and the archive.
    11. LOCAL WORKDIR CLEANUP (keepdata=no deletes it after the rsync).
    12. DATASHEET CLEANUP: the source .xlsx file(s) are moved to the team
        "TO delete" folder for later cleanup (only the exact source file(s);
        any other .xlsx in the working dir is left untouched).
    13. ``exit ${_overall_rc}`` — any failed sample makes the script end non-zero.

Parity note: this generator mirrors the GEX/ARC parity-pass (sample
selection, run-folder disambiguation, per-sample parameters, parallel
dispatch, status sweep, merged summary CSV, single md5sums file, .sh
provenance copy, TO_DELETE datasheet parking, scriptversion stamp), adapted
to single-ATAC runs — there is no Multiome grouping, no libraries CSV, and
no TSO% report (ATAC has no Read-2 intact-TSO to measure).
"""

import argparse
import os
import re
import sys
import time
import unicodedata
import warnings
from collections import OrderedDict

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
    run_base_dir,        # run name -> base folder (all sequencers: AVT/MSH/MSQ/NVQ…)
    FASTQ_ID_PATTERN,
    AVT_PATTERN,
    avt_run_tag,         # multi-run tag (AVT0001_AVT0002) for script/folder names
    inherit_file_avt,    # rows inherit the run named in their OWN file
    resolve_unassigned_avt,  # GEX 'needs_avt' prompt (shared) for run-less rows
    norm,
    build_merged_cell_lookup,
    get_cell_value,
    read_user_field,
    set_non_interactive,
    setup_readline_completer,
    ask_avt_run,
    ask_reference_genome,
    LOCAL_CORES,        # cellranger-atac --localcores (shared default "32")
    LOCAL_MEM,          # cellranger-atac --localmem (shared default "126")
    # The same shared helpers the GEX/ARC generators use for their
    # sample-selection prompt, run-folder disambiguation, and preferred-dir
    # fastq resolution (single source of truth in common.intake).
    ask_sample_selection,
    NGS_RUNS_ROOT,          # prefill root for the add-sample prompt
    list_sample_dirs,            # candidate dirs for one (AVT, sample)
    choose_sample_dir,           # pick ONE when >1 candidate
    find_fastq_path_with_pref,   # resolve an avid id honouring the pick
)

# The central archive where every ATAC run's folder is kept.
# (Fixed target; the generated script rsyncs the whole working folder there.)
NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/ZZ_10XG/atac"

# Root on the *analysis* server where cellranger-atac output + post-processing
# are written. On that host this is a real, mounted path; on the dev box it is
# not expected to exist. Never derive from $PWD — the script must land in the
# same place regardless of the caller's cwd.
SCRATCH_ROOT = "/mnt/data/CellRangerCountOutput"

# "TO delete" folder where processed datasheets are parked for later cleanup.
# Two destinations, tried in order (same folders as the GEX/ARC generators);
# the generated script uses the first folder that exists and creates neither.
DELETE_PRIMARY = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CR_pipelines_Daniel/TO_DELETE"
)
DELETE_FALLBACK = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CRcount_GEXonly/TO_DELETE"
)

# Provenance: version stamped into the generated .sh header
# (scriptversion=...) and into the script's own source as SCRIPT_VERSION.
# Generator version == generated-script version; the kept .sh in the archive
# tells you which one produced it.  Bump on any change to the generator.
# 1.3.5 = archive root moved off the sequencer folder: results now deploy
#         to the shared, sequencer-agnostic /mnt/PTEGraw/AA/NGSruns/ZZ_10XG/.
SCRIPT_VERSION = "1.3.5"

# ---------------------------------------------------------------------------
# Reference genomes (auto-detected from the live filesystem, GEX/ARC pattern)
# ---------------------------------------------------------------------------
# Reference-genome root.  The genome folder inside each ``current_*_arc_ATAC``
# directory is never hard-coded here: each of them holds exactly one genome
# folder and that folder's name can change on the server (rename/replace), so
# the real path is resolved once at startup from the live filesystem.
CELLRANGER_REFERENCES_ROOT = "/home/gecf/references/cellranger_references"
_HUMAN_CURRENT_DIR = os.path.join(
    CELLRANGER_REFERENCES_ROOT, "current_human_arc_ATAC")
_MOUSE_CURRENT_DIR = os.path.join(
    CELLRANGER_REFERENCES_ROOT, "current_mouse_arc_ATAC")


def _is_version_check() -> bool:
    """True when the user is only asking for ``-v``/``--version`` so the
    import-time reference-genome scan can stay quiet (argparse hasn't run
    yet but ``sys.argv`` already carries the flag).

    (ported from the completed GEX generator)
    """
    return ("-v" in sys.argv or "--version" in sys.argv)


def _resolve_reference_genome(current_dir, label):
    """Return the genome path inside ``current_dir`` — its single subdirectory.

    Only *current_* is guaranteed to persist across the analysis servers;
    the genome folder it contains (refdata-cellranger-arc-...) can be
    renamed or replaced without the generator being updated.  Non-recursive,
    so an inner ``old/`` tier never leaks in (only the current genome counts).

    (ported from the completed GEX generator)
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


def folder_safe(name: str) -> str:
    """Reduce *name* (a requester string from the datasheet) to a safe
    path-segment for the working-folder name.

    Unicode is NFKD-decomposed and non-ASCII leftovers dropped, so accented
    text like ``Cèsar`` becomes ``Cesar``; any other character that cannot
    live in a CIFS path / URL is replaced with ``_`` and a trailing dot/space
    is stripped (Windows).  Empty input yields ``""`` so the caller can omit
    the whole ``_<user>`` segment instead of emitting a stray underscore.

    (ported from the completed GEX/ARC generators — single convention)
    """
    if not name:
        return ""
    s = unicodedata.normalize("NFKD", name)
    s = s.encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    s = re.sub(r"\.{2,}", ".", s)
    s = re.sub(r"_{2,}", "_", s)
    return s.strip("._-")


def latest_avt(avt_names):
    """Return the detected AVT run with the largest numeric suffix.

    (ported from the completed ARC generator)
    """
    return max(
        avt_names,
        key=lambda name: int(re.search(r"\d+", name).group()),
        default=None,
    )


# PER-SAMPLE PARAMETERS for cellranger-atac count (parity: same "one editable
# block per sample, blank = drop the flag" pattern as GEX/ARC — the knobs are
# the ATAC-relevant ones; cellranger-atac has no --create-bam /
# --cell-annotation-model (GEX/ARC-only). ATAC-specific tuning goes in
# custom_options — only pass flags this cellranger-atac build accepts.
# CPU/memory caps are NOT per-sample (GEX/ARC convention): localcores /
# localmem are run-wide, emitted in the RUN-WIDE PARAMETERS block and baked
# into the generated script from the shared imports below.
_ATAC_PER_SAMPLE_PARAMS = (
    ("custom_options", "",
     "extra flags as-is, e.g. --chemistry=ARC-v1  --peaks=/refs/peaks.bed"),
)


# ---------------------------------------------------------------------------
# Excel parsing (sheet preference / blanks / AVT detection — unchanged)
# ---------------------------------------------------------------------------

def ranked_sheet_names(wb) -> list:
    """Order workbook sheets by data-sheet precedence.

    Precedence rules:
      1. A ``User_x1x_copy`` sheet (if present) is always first.
      2. Sheets whose names contain 'spike', 'several user', 'guideline',
         or 'index table' (case-insensitive) are never picked — they are
         legend/distractor sheets (GECF GUIDELINES legend, or the
         'Index table MM.YY' control-library reference).
      3. The rest, in workbook order (the first one that exposes the data
         header row wins).
    """
    names = list(wb.sheetnames)
    preferred, rest = [], []
    for n in names:
        lo = n.lower()
        if ("spike" in lo or "several user" in lo or "guideline" in lo
                or "index table" in lo):
            continue
        (preferred if _run_details_name(n) else rest).append(n)
    return preferred + rest


def _run_details_name(name: str) -> bool:
    return name.lower().replace("-", "_").strip() == "user_x1x_copy"


def _locate_header(ws, merged):
    """Return ``(fastq_col, sample_col, header_row)`` or ``(None, None, None)``.

    The GECF submission layout marks the sample table with a header row that
    names both the fastq-ID column (a ``GECF fastq ID`` header cell) and the
    ``Sample name`` column.  Column positions therefore come from the header
    cells (case/whitespace-insensitive via ``norm``), never from fixed
    A/B/C/D assumptions — ARC-parity.
    """
    max_row = ws.max_row or 0
    for row_idx in range(1, min(25, max_row) + 1):
        fastq_col = sample_col = None
        for col_idx in range(1, ws.max_column + 1):
            v = norm(get_cell_value(ws, merged, row_idx, col_idx))
            if fastq_col is None and "fastq id" in v:
                fastq_col = col_idx
            if sample_col is None and "sample name" in v:
                sample_col = col_idx
        if fastq_col is not None and sample_col is not None:
            return fastq_col, sample_col, row_idx
    return None, None, None


def _row_avt(ws, merged, row_idx, fastq_col, sample_col):
    """Per-row run name (parity item 11, shared with the ARC generator).

    The run name may sit in any cell of the table row except the fastq-ID /
    sample-name columns (which carry the avid id / the sample label, not
    the run).  Returns the upper-cased run name (``AVT0XXX``) or ``None``
    when this row does not name a run.
    """
    for col_idx in range(1, ws.max_column + 1):
        if col_idx in (fastq_col, sample_col):
            continue
        val = get_cell_value(ws, merged, row_idx, col_idx)
        if val is None:
            continue
        m = AVT_PATTERN.search(str(val))
        if m:
            return m.group(0).upper()
    return None


def _sheet_rows_header_driven(ws, merged):
    """Extract (sample, fastq, avt) triples from one GECF sample table.

    Returns ``[]`` when the sheet does not carry the header row (fastq-ID
    column + Sample name).  Rows are collected below the header until the
    table ends: a ``Run details`` row, or two consecutive non-data rows.
    """
    fastq_col, sample_col, header_row = _locate_header(ws, merged)
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
        sample_val = get_cell_value(ws, merged, row_idx, sample_col)
        fastq_val = get_cell_value(ws, merged, row_idx, fastq_col)
        sample_txt = str(sample_val).strip() if sample_val is not None else ""
        fastq_txt = str(fastq_val).strip() if fastq_val is not None else ""
        sample_txt = "" if sample_txt == "0" else sample_txt
        fastq_txt = "" if fastq_txt == "0" else fastq_txt

        if "run details" in norm(sample_txt) or "run details" in norm(fastq_txt):
            break

        if sample_txt and fastq_txt and FASTQ_ID_PATTERN.match(fastq_txt):
            rows.append((sample_txt, fastq_txt,
                         _row_avt(ws, merged, row_idx, fastq_col, sample_col)))
            blank_streak = 0
        else:
            if sample_txt and not fastq_txt:
                print(f"  Warning: sample '{sample_txt}' has no valid fastq "
                      f"ID — row skipped.")
            blank_streak += 1
            if rows and blank_streak >= 2:
                break
    return rows


def _parse_positional(wb):
    """Legacy positional fallback: col A sample, col B fastq, col C FASTQ
    path, col D AVT — the pre-GECF-template layout.  Used only when no
    sheet exposes the fastq-ID / Sample-name header row (keeps the old
    behaviour byte-for-byte for such workbooks).
    """

    def is_blank(v) -> bool:
        return v is None or (isinstance(v, str) and v.strip() == "")

    def _is_header_text(text: str) -> bool:
        n = norm(text)
        return ("sample" in n) or ("run details" in n)

    all_pairs = []
    seen_pairs = set()  # (sample, fastq) to avoid duplicate runs
    best_sheet = None
    best_score = -1
    avt_from_sheet = None

    for sheet_name in ranked_sheet_names(wb):
        ws = wb[sheet_name]
        rows = list(ws.iter_rows(values_only=True))
        score = 0
        for r in rows[:10]:
            vals = [c for c in r if c is not None]
            if vals:
                s = " ".join(str(v) for v in vals)
                if _is_header_text(s) and any(
                    AVT_PATTERN.search(str(v)) for v in vals
                ):
                    m = AVT_PATTERN.search(s)
                    if m and not avt_from_sheet:
                        avt_from_sheet = m.group()
        # Count candidate data rows (sample cell non-blank, not the word
        # "sample", and at least one of fastq/path is present)
        for r in rows:
            c_sample = r[0] if len(r) > 0 else None
            c_fastq = r[1] if len(r) > 1 else None
            c_path = r[2] if len(r) > 2 else None
            if is_blank(c_sample):
                continue
            if str(c_sample).strip().lower() in {"sample", "name", ""}:
                continue
            if is_blank(c_fastq) and is_blank(c_path):
                continue  # no fastq at all -> not a data row
            score += 1
        if score > best_score:
            best_sheet = ws
            best_score = score

    if best_sheet is None:
        return [], None

    rows = list(best_sheet.iter_rows(values_only=True))
    # find header row for column positions
    hidx = 0
    for i, r in enumerate(rows):
        v0 = r[0] if r else None
        if v0 is not None and _is_header_text(str(v0)):
            hidx = i
            break
    for r in rows[hidx + 1:]:
        vals = list(r)
        c_sample = vals[0] if len(vals) > 0 else None
        c_fastq = vals[1] if len(vals) > 1 else None
        c_path = vals[2] if len(vals) > 2 else None
        c_avt = vals[3] if len(vals) > 3 else None
        # stop at a blank sample cell after we've started collecting data
        if is_blank(c_sample):
            if all_pairs:
                break
            continue
        sample = str(c_sample).strip()
        if sample.lower() in {"name", "sample"}:
            continue
        fastq = ""
        if not is_blank(c_fastq):
            fastq = str(c_fastq).strip()
        elif not is_blank(c_path):
            # Use the basename as the fastq id; the full path is kept for
            # resolution (it may already be the leaf directory).
            fastq = os.path.basename(str(c_path).strip())
        else:
            print(f"  Warning: sample '{sample}' has no fastq ID and no "
                  f"FASTQ path — row skipped.")
            continue
        # Per-row AVT (optional; blank -> inherit the run named in this file)
        row_avt = None
        if not is_blank(c_avt):
            m = AVT_PATTERN.search(str(c_avt))
            if m:
                row_avt = m.group()
        key = (sample, fastq)
        if key in seen_pairs:
            continue
        seen_pairs.add(key)
        all_pairs.append((sample, fastq, row_avt))

    return all_pairs, avt_from_sheet


def _sheet_avt(ws):
    """Run name from a sheet's first rows (e.g. the
    ``AVITI24 run AVT0565`` title row), or ``None``."""
    for r in list(ws.iter_rows(values_only=True))[:12]:
        for c in r:
            if c is None:
                continue
            m = AVT_PATTERN.search(str(c))
            if m:
                return m.group(0).upper()
    return None


def parse_xlsx(path: str):
    """Return (list_of_pairs, detected_avt).

    ``list_of_pairs`` items are ``(sample, fastq, per_row_avt)`` —
    ``per_row_avt`` is the AVT named on that row (or None when the row
    names no run).

    Parsing rules:
      * GECF submission layout (preferred, ARC-parity): the sheet whose
        header row names both a fastq-ID column (a ``GECF fastq ID`` header
        cell) and a ``Sample name`` column.  Column positions come
        from the header cells, data rows are collected below the header
        until the ``Run details`` block (or two consecutive non-data rows).
        When several sheets qualify, the one with the most valid rows wins;
        ranked_sheet_names precedence breaks ties (User_x1x_copy first).
      * Legacy layout fallback (no such header anywhere): positional col A
        sample, col B fastq, col C FASTQ path, col D AVT — the old
        behaviour, unchanged.
      * detected_avt: the first per-row AVT, else the sheet title / early
        rows, else the file name.
    """
    wb = openpyxl.load_workbook(path, data_only=True)

    best_rows, best_sheet = [], None
    for sheet_name in ranked_sheet_names(wb):
        ws = wb[sheet_name]
        try:
            rows = _sheet_rows_header_driven(ws, build_merged_cell_lookup(ws))
        except Exception:
            rows = []
        if len(rows) > len(best_rows):
            best_rows, best_sheet = rows, ws

    if best_rows:
        # De-duplicate, preserving first-seen order.
        seen = set()
        deduped = []
        for sample, fastq, avt in best_rows:
            key = (sample.lower(), fastq.lower())
            if key in seen:
                continue
            seen.add(key)
            deduped.append((sample, fastq, avt or None))
        sheet_avt = next((a for (_s, _f, a) in deduped if a), None)
        if sheet_avt is None and best_sheet is not None:
            sheet_avt = _sheet_avt(best_sheet)
        if sheet_avt is None:
            m = AVT_PATTERN.search(os.path.basename(path))
            if m:
                sheet_avt = m.group(0).upper()
        if sheet_avt:
            return ([
                (s, f, (a or sheet_avt)) for (s, f, a) in deduped
            ], sheet_avt)
        return deduped, None

    return _parse_positional(wb)


def print_pairs_table(pairs):
    """Human-readable preview: one line per (sample, fastq) pair."""
    if not pairs:
        return
    print("")
    print("Detected samples (one `cellranger-atac count` per sample):")
    shown = OrderedDict()
    for sample, fastq, _avt in pairs:
        shown.setdefault(sample, []).append(fastq)
    for i, (sample, fqs) in enumerate(shown.items(), start=1):
        print(f"  {i}. {sample}:  " + ", ".join(f for f in fqs if f))


# ---------------------------------------------------------------------------
# FASTQ leaf resolution (classic full-tree walk — unchanged behaviour)
# ---------------------------------------------------------------------------

def find_fastq_path(avt_run: str, fastq_id: str):
    """Find ``fastq_id`` directory anywhere under the AVT run tree, or None.

    (Kept as the classic fallback walk; ``find_fastq_path_with_pref()`` from
    common.intake is preferred when a sample dir was disambiguated.)
    """
    base = run_base_dir(avt_run)
    if not os.path.isdir(base):
        return None
    for r, dirs, _ in os.walk(base):
        for d in dirs:
            if d == fastq_id:
                return os.path.join(r, d)
    return None


# ---------------------------------------------------------------------------
# Shell script emission
# ---------------------------------------------------------------------------

def merge_sample_order(orders):
    """Merge the sample orders of several datasheets into one order.

    *orders* is a list of sample-name lists, one per input Excel file, in
    the order the user listed the files on the command line (each list may
    be raw row order — duplicates are collapsed to first-seen).

    Result:
      * if all files AGREE on the relative order of the samples they share
        (no pair of samples is ordered oppositely in two files), the single
        common order across all files — samples that only some files list
        still land at their agreed position;
      * otherwise the FIRST file's order wins, with samples that appear
        only in later files appended after it (first-seen order in the
        remaining files).

    Every sample present in any input file appears exactly once in the
    result — files may cover different sample subsets, nothing is dropped.
    """
    orders = [list(dict.fromkeys(o)) for o in orders]
    orders = [o for o in orders if o]
    if len(orders) <= 1:
        return orders[0] if orders else []

    union_list = list(dict.fromkeys(s for o in orders for s in o))

    # Constraint graph: each file is a chain — consecutive samples a, b in
    # the same file must keep that relative order in the merged result.
    preds = {s: set() for s in union_list}
    succs = {s: set() for s in union_list}
    for o in orders:
        for a, b in zip(o, o[1:]):
            succs[a].add(b)
            preds[b].add(a)

    # Deterministic tie-break: prefer the sample that sits earliest in the
    # file containing it (file #, then row) — keeps the first file's order
    # at the front of the merged result.
    rank = {}
    for s in union_list:
        for i, o in enumerate(orders):
            if s in o:
                rank[s] = (i, o.index(s))
                break

    # Topological sort.  A deadlock (ready set empty while samples remain)
    # means two files order some pair of samples oppositely -> they
    # DISAGREE -> fall through to the first-file fallback below.
    live = set(union_list)
    emitted = []
    while live:
        ready = [s for s in live if not preds[s]]
        if not ready:
            break
        s = min(ready, key=lambda s: rank[s])
        emitted.append(s)
        live.discard(s)
        for t in succs[s]:
            preds[t].discard(s)

    if len(emitted) == len(union_list):
        return emitted  # all files agree -> common order

    # Disagreement -> first file's order + first-seen order of the rest.
    merged = list(orders[0])
    seen = set(merged)
    for o in orders[1:]:
        for s in o:
            if s not in seen:
                merged.append(s)
                seen.add(s)
    return merged


def build_script_full(reference_genome, pairs, user="", main_avt=None,
                      avt_tag=None, source_files=None, resolve=None,
                      number_prefix=True, ordered_samples=None):
    """Emit the ``cellranger-atac_<AVT>.sh`` text (multi-run inputs:
    ``cellranger-atac_AVT0001_AVT0002.sh``).

    ``pairs`` is a list of ``(sample, fastq_id, per_row_avt)``; ``resolve``
    maps ``(sample, fastq_id)`` to a leaf directory (or None) and is called
    once per pair at generation time, baking literal paths into the .sh.
    When ``resolve`` is None the classic walk (``find_fastq_path``) under the
    main AVT run tree is used.
    *number_prefix* (default True): zero-pad the sample index (01_, 02_,
    …) onto each per-sample run folder (cellranger-atac count --id and
    --output-dir) and its companion files, keeping the folders in sample
    order even when the names are not alphabetical.  Pass False (CLI:
    --no-number-prefix) for bare sample names (the previous behaviour).

    *ordered_samples* (optional): the merged cross-file sample order
    (``merge_sample_order``).  When given, it renumbers the samples;
    any sample it does not list (e.g. added by full path in the sample
    prompt) keeps first-seen order and is appended at the end.


    Layout (top -> bottom) — the same section order as the GEX/ARC generated
    scripts, adapted to single-ATAC runs:
      1.  Header / constants (timestamp, reference, source files, requester,
          scriptversion stamp, avt_run, user).
      2.  RUN-WIDE PARAMETERS (editable at the top): keepdata, max_parallel,
          localcores, localmem, scratch_root, server_folder_name,
          NGSruns_archive, then all FASTQ path variables (resolved at
          generation time), then derived vars (Server_folder, summaries
          folder, archive_dir), mkdir, ref check.
      3.  PER-SAMPLE PARAMETERS: one block per sample — custom_options;
          blank a value to drop that flag for that sample. CPU/memory caps
          are run-wide (localcores / localmem in RUN-WIDE PARAMETERS), per
          the GEX/ARC convention.
      4.  One ``atac_sample_<i>()`` function per sample: runs
          ``cellranger-atac count`` (flags from the ``*_<i>`` variables),
          checks for leftover tmp folders, deletes the SC_ATAC_COUNTER_CS
          scratch tree, harvests the web summary into ``<AVT>_Summaries/``.
      5.  DISPATCH: serial (max_parallel=1) or rolling pool of max_parallel
          concurrent atac samples (default 3).
      6.  Per-sample status sweep: [OK]/[FAIL] per sample,
          ``_overall_rc`` accumulator.
      7.  Merged summary metrics: one ``all_summary.csv`` with a Samplename
          column first (``<AVT>_Summaries/all_summary.csv``).
      8.  rsync of the whole working folder to the central archive.
      9.  Whole-folder ``md5sums.txt`` + rsync of the single checksums file.
      10. ``.sh`` copy to work folder + archive (provenance).
      11. keepdata cleanup.
      12. TO_DELETE parking of the source .xlsx files.
      13. ``exit ${_overall_rc}``.
    """
    # Group by sample (first-seen order), keeping each fastq only once.
    groups = OrderedDict()
    for sample, fastq, _row_avt in pairs:
        fqs = groups.setdefault(sample, [])
        if not any(f == fastq for f in fqs):
            fqs.append(fastq)
    samples = list(groups)
    if ordered_samples:
        # Cross-file consensus order (merge_sample_order) renumbers the
        # samples; anything it does not list is appended after, first-seen.
        ordered_set = set(ordered_samples)
        samples = [si for si in ordered_samples if si in groups] + \
                  [si for si in samples if si not in ordered_set]
    N = len(samples)
    # Per-sample run folder name: 01_, 02_, … (ON by default) keeps the
    # folders in sample order even when the names are not alphabetical.
    # --no-number-prefix -> bare sample name (previous behaviour).
    run_names = [f"{i:02d}_{s}" if number_prefix else s
                 for i, s in enumerate(samples, 1)]
    _resolve = resolve if resolve is not None else (
        lambda sn, fid: find_fastq_path(main_avt or "", fid))
    # Identity of the script name + working/archive folder: the combined
    # multi-run tag (AVT0001_AVT0002, common.intake.avt_run_tag) when the
    # inputs span several runs, else the representative AVT.  Per-FASTQ
    # resolution still uses each row's own AVT (main_avt only as fallback).
    tag = avt_tag or main_avt
    L = []

    # ---- 1) header / constants -----------------------------------------------------------------
    L.append("#!/bin/bash")
    L.append("# cellranger-atac — one `cellranger-atac count` per sample (ATAC only).")
    L.append(f"# Generated by generate_cellranger_template_atac.py ({N} sample(s)).")
    L.append(f"# Generated at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    _srcs = list(dict.fromkeys(source_files or []))
    if _srcs:
        L.append("# Source file(s):")
        for _s in _srcs:
            L.append(f"#   {os.path.basename(_s)}")
    else:
        L.append("# Source file(s): (not provided)")
    L.append(f"# Requester: {user or 'n/a'}")
    # Provenance stamp in the kept .sh: which GENERATOR produced this script.
    L.append("scriptversion=" + SCRIPT_VERSION)
    L.append("")
    L.append("TIMESTAMP=$(date +%Y%m%d_%H%M%S)")
    L.append('reference_genome="' + reference_genome + '"')
    L.append(f"avt_run='{tag}'")
    # The requester string (datasheet ``User:`` cell) is not validated; clean
    # it to a CIFS-/URL-safe path segment HERE at generation time so an
    # accented / punctuation-laden name cannot later break the working-folder
    # name.  If the result is empty, the ``${user:+_${user}}`` conditional
    # below simply omits the ``_user`` segment.
    L.append("user='" + folder_safe(user or "") + "'")
    L.append("")

    # ---- 2) RUN-WIDE PARAMETERS (editable block) -----------------------------------------------
    L.append("# ============================================================")
    L.append("# RUN-WIDE PARAMETERS  (defaults — edit here)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   keepdata : "yes" = keep the local working folder after the run;')
    L.append('#               "no"  = delete it (data already rsynced to the archive).')
    L.append('keepdata="no"')
    L.append('#   max_parallel : how many cellranger-atac counts to run at the same time.')
    L.append('#   1 = serial (one sample after another); N>1 = rolling pool of N jobs,')
    L.append('#          the next sample starts as soon as a slot frees up.')
    L.append('#   localcores : CPU cores per cellranger-atac count (applies to ALL samples);')
    L.append('#                blank = cellranger default.')
    L.append('#   localmem   : RAM in GB per cellranger-atac count (applies to ALL samples);')
    L.append('#                blank = cellranger default.')
    L.append('#   With max_parallel>1, total = max_parallel x (localcores cores, localmem GB),')
    L.append("#   so all three must fit within the server's resources together.")
    L.append('max_parallel=3')
    L.append(f'localcores="{LOCAL_CORES}"')
    L.append(f'localmem="{LOCAL_MEM}"')
    L.append('#   scratch_root : root on the analysis server where cellranger-atac output is written')
    L.append('#                  (fixed convention — change only if the mount differs on your server).')
    L.append('scratch_root="' + SCRATCH_ROOT + '"')
    L.append('#   server_folder_name : working + archive leaf folder name.')
    L.append('#                  default is <AVT>[_user]_<TIMESTAMP>; change only if needed.')
    L.append('server_folder_name="${avt_run}${user:+_${user}}_$TIMESTAMP"')
    L.append('#   NGSruns_archive: central archive destination (fixed target).')
    L.append('NGSruns_archive="' + NGSruns_ARCHIVE_DIR + '"')
    L.append("")

    # ---- FASTQ PATHS (one variable per sample FASTQ, resolved at generation time) ------
    L.append("# ============================================================")
    L.append("# FASTQ PATHS  (one per sample FASTQ — resolved at generation time)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   If a variable is set to the empty string ("") or has a FIXME_LOCATE_')
    L.append('#   marker, the generator could not resolve the leaf directory — set the')
    L.append('#   real path by hand before running.')
    any_missing = False
    for idx, sample in enumerate(samples, start=1):
        fqs = groups[sample]
        L.append("")
        L.append(f"# --- Sample {idx}: {sample}  ({len(fqs)} FASTQ) ---")
        for j, fid in enumerate(fqs, start=1):
            leaf = _resolve(sample, fid)
            if leaf:
                L.append(f'fastq_{idx}_{j}="{leaf}"')
            else:
                L.append(f'fastq_{idx}_{j}=""   '
                         f'# FIXME_LOCATE_{fid or j}: leaf not resolved — set by hand before running')
                any_missing = True
    if any_missing:
        print("  Warning: some fastq leaves could not be resolved at generation time.",
              file=sys.stderr)
        print("  The flagged fastq_<idx>_<j> variables (FIXME_LOCATE_ markers) "
              "must be set before running.", file=sys.stderr)
    L.append("")

    # ---- derived working-dir variables ---------------------------------------------------------
    L.append('Server_folder="${scratch_root}/atac/${server_folder_name}"')
    L.append('Summaries_folder="${Server_folder}/' + (tag or "AVT")
             + '_Summaries"')
    L.append('Web_summaries="${Summaries_folder}/web_summaries"')
    L.append('archive_dir="${NGSruns_archive}/${server_folder_name}"')
    L.append('mkdir -p "${Server_folder}" "${Web_summaries}" "${archive_dir}"')
    L.append("")

    # ---- reference-genome existence check ------------------------------------------------------
    L.append("# Check that the reference genome directory exists")
    L.append('if [[ ! -d "${reference_genome}" ]]; then')
    L.append('    echo "Reference genome ${reference_genome} does not exist. Aborting." >&2')
    L.append("    exit 1")
    L.append("fi")
    L.append("")

    # ---- 3) PER-SAMPLE PARAMETERS (one block per sample) ---------------------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE PARAMETERS  (one block per sample — edit values in place)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   * blank a value ( value="" )  -> that flag is DROPPED for that sample')
    L.append("#   * custom_options : extra flags as-is, e.g. --peaks=/refs/peaks.bed")
    L.append("#   * localcores / localmem are RUN-WIDE (see RUN-WIDE PARAMETERS);")
    L.append("#     they are not per-sample — same convention as the GEX/ARC scripts.")
    L.append("# ============================================================")
    for idx, sample in enumerate(samples, start=1):
        L.append(f"# Sample {idx}: {sample}")
        for suffix, default, comment in _ATAC_PER_SAMPLE_PARAMS:
            var = f"{suffix}_{idx}"
            base_str = f'{var}="{default}"'
            pad = max(1, 56 - len(base_str))
            L.append(f"    {base_str}{' ' * pad}# {comment}")
    L.append("")

    # Boundary marker: emitted below all editable blocks (RUN-WIDE / FASTQ /
    # PER-SAMPLE) so the user/code boundary in the generated script is
    # unambiguous.
    L.append("# " + "=" * 73)
    L.append("#   ***  END OF USER-EDITABLE SECTION  ***")
    L.append("#")
    L.append("#   Everything ABOVE is yours to edit: RUN-WIDE PARAMETERS, FASTQ PATHS,")
    L.append("#   PER-SAMPLE PARAMETERS.")
    L.append("#")
    L.append("#   Everything BELOW is pipeline machinery — pure code, DO NOT EDIT.")
    L.append("# " + "=" * 73)
    L.append("")

    # ---- 4) atac_sample_<i>() functions ----------------------------------------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE FUNCTIONS  (one atac_sample_<i> per sample)")
    L.append("#   Each function:")
    L.append("#     * runs cellranger-atac count (per-sample flags from the *_<i> variables above;")
    L.append("#       a blank value drops that flag for that sample only — CPU/memory caps")
    L.append("#       localcores / localmem are run-wide)")
    L.append("#     * checks for leftover tmp folders -> <sample>_error_log.txt")
    L.append("#     * deletes the non-essential SC_ATAC_COUNTER_CS scratch tree")
    L.append("#     * harvests the web summary into <AVT>_Summaries/web_summaries/")
    L.append("# ============================================================")

    for idx, sample in enumerate(samples, start=1):
        fqs = groups[sample]
        n_fq = len(fqs)
        desc = ",".join(f for f in fqs if f) or sample
        run_name = run_names[idx - 1]
        L.append(f"# --- atac_sample_{idx}() — Sample {idx}: {sample}  ({n_fq} FASTQ) ---")
        L.append(f"atac_sample_{idx}() {{")

        # Fastq value(s) for --fastqs (comma-joined when the sample has more than one)
        fq_vars = [f"fastq_{idx}_{j}" for j in range(1, n_fq + 1)]
        fq_arg = ",".join("${" + v + "}" for v in fq_vars)

        flag_lines = [
            'cellranger-atac count',
            f'    --id={run_name}',
            f'    --description={desc}',
            '    --reference="${reference_genome}"',
            f'    --fastqs={fq_arg}',
            f'    --output-dir=${{Server_folder}}/{run_name}',
            # CPU/memory caps are run-wide (localcores / localmem), not per-sample
            '    ${localcores:+--localcores=$localcores}',
            '    ${localmem:+--localmem=$localmem}',
        ]
        # Per-sample conditional knobs. Blank value -> flag omitted for this
        # sample only. custom_options is deliberately UNQUOTED on the last
        # line so it word-splits into multiple --flags.
        for suffix, _default, _comment in _ATAC_PER_SAMPLE_PARAMS:
            v = f"{suffix}_{idx}"
            if suffix == "custom_options":
                flag_lines.append(f"    ${{{v}}}")
            else:
                flag = "--" + suffix.replace("_", "-")
                flag_lines.append(f"    ${{{v}:+{flag}=${v}}}")

        # Write the exact command to <sample>_script_log.txt (one readable line)
        L.append(f'    _cmd_log="${{Server_folder}}/{run_name}_script_log.txt"')
        _log_body = " ".join(fl.strip() for fl in flag_lines)
        L.append(f"    echo {_log_body!r} > \"${{_cmd_log}}\"")
        L.append("")
        # Execute the command (same flag_lines list, line-continuation formatted)
        for i, fl in enumerate(flag_lines):
            cont = " \\" if i < len(flag_lines) - 1 else ""
            L.append(f"    {fl}{cont}")
        L.append("")

        # Completion check: leftover tmp folders indicate cellranger-atac count
        # did not finish (the same completion marker the GEX/ARC scripts use).
        L.append(f'    _err_log="${{Server_folder}}/{run_name}_error_log.txt"')
        L.append(f'    find "${{Server_folder}}/{run_name}" -mindepth 2 -maxdepth 2 '
                 f'-type d -name "tmp" -print 2>/dev/null > "${{_err_log}}"')
        L.append('    if [ -s "${_err_log}" ]; then')
        L.append(f'        echo "Error: cellranger-atac count for {sample} did not complete '
                 f'(leftover tmp folders detected).\\nSee ${{Server_folder}}/{run_name}_error_log.txt" >&2')
        L.append('    else')
        L.append('        rm -f "${_err_log}"')
        L.append('    fi')
        L.append("")

        # Delete the non-essential SC_ATAC_COUNTER_CS intermediate tree
        L.append(f'    rm -r "${{Server_folder}}/{run_name}/SC_ATAC_COUNTER_CS" 2>/dev/null || true')
        L.append("")

        # Harvest the web summary into <AVT>_Summaries
        L.append(f'    if [ -f "${{Server_folder}}/{run_name}/outs/web_summary.html" ]; then')
        L.append(f'        cp "${{Server_folder}}/{run_name}/outs/web_summary.html" '
                 f'"${{Web_summaries}}/{run_name}_WebSummary.html"')
        L.append('    else')
        L.append(f'        echo "  [warn] no web summary for {sample}" >&2')
        L.append('    fi')

        L.append("}")
        L.append("")

    # ---- 5) DISPATCH (serial or rolling pool of max_parallel) -----------------------------------
    L.append("# ============================================================")
    L.append("# DISPATCH  (run atac_sample_1 .. atac_sample_" + str(N) + ")")
    L.append("#   max_parallel=1 : serial, one sample after another")
    L.append("#   max_parallel=N : rolling pool of N concurrent cellranger-atac counts —")
    L.append("#                    the next sample starts as soon as a slot frees up")
    L.append("#   (default 3 in the generated script — see RUN-WIDE PARAMETERS above)")
    L.append("# ============================================================")
    L.append("NSAMPLES=" + str(N))
    L.append('case "${max_parallel}" in')
    L.append("    ''|*[!0-9]*)")
    L.append('        echo "ERROR: max_parallel must be a positive integer" >&2')
    L.append("        exit 1")
    L.append("    ;;")
    L.append("esac")
    L.append('if [ "${max_parallel}" -lt 1 ]; then')
    L.append('    echo "ERROR: max_parallel must be >= 1" >&2')
    L.append("    exit 1")
    L.append("fi")
    L.append("# Cap max_parallel at NSAMPLES (no point requesting more slots than samples)")
    L.append('if [ "${max_parallel}" -gt "${NSAMPLES}" ]; then')
    L.append('    echo "Note: only ${NSAMPLES} sample(s) — capping max_parallel ${max_parallel} -> ${NSAMPLES}." >&2')
    L.append("    max_parallel=${NSAMPLES}")
    L.append("fi")
    L.append("")
    L.append('if [ "${max_parallel}" -eq 1 ]; then')
    L.append("    # --- SERIAL ---")
    L.append("    for (( _s = 1; _s <= ${NSAMPLES}; _s++ )); do")
    L.append("        atac_sample_${_s}")
    L.append("    done")
    L.append("else")
    L.append("    # --- ROLLING POOL OF ${max_parallel} CONCURRENT JOBS ---")
    L.append("    # Slot arrays: _pids[_i] = PID, _idxs[_i] = sample number")
    L.append("    declare -a _pids=() _idxs=()")
    L.append("    for (( _i = 0; _i < ${max_parallel}; _i++ )); do")
    L.append('        _pids[$_i]=""')
    L.append("        _idxs[$_i]=0")
    L.append("    done")
    L.append("    _next=1")
    L.append("    while :; do")
    L.append("        # 1) reap finished slots")
    L.append("        for (( _i = 0; _i < ${max_parallel}; _i++ )); do")
    L.append('            if [ -n "${_pids[$_i]}" ] && ! kill -0 "${_pids[$_i]}" 2>/dev/null; then')
    L.append('                wait "${_pids[$_i]}"')
    L.append("                _rv=$?")
    L.append('                echo "[REAP] sample #${_idxs[$_i]} finished (rc=${_rv})"')
    L.append('                _pids[$_i]=""')
    L.append("            fi")
    L.append("        done")
    L.append("        # 2) fill free slots with the next sample (if any remain)")
    L.append("        _launched=0")
    L.append("        for (( _i = 0; _i < ${max_parallel}; _i++ )); do")
    L.append('            if [ -z "${_pids[$_i]}" ] && [ "${_next}" -le "${NSAMPLES}" ]; then')
    L.append("                atac_sample_${_next} &")
    L.append("                _pids[$_i]=$!")
    L.append("                _idxs[$_i]=${_next}")
    L.append("                _next=$(( _next + 1 ))")
    L.append("                _launched=1")
    L.append("            fi")
    L.append("        done")
    L.append("        # 3) exit when all samples launched AND all jobs reaped")
    L.append("        _busy=0")
    L.append("        for (( _i = 0; _i < ${max_parallel}; _i++ )); do")
    L.append('            [ -n "${_pids[$_i]}" ] && _busy=1')
    L.append("        done")
    L.append('        if [ "${_next}" -gt "${NSAMPLES}" ] && [ "${_busy}" -eq 0 ]; then')
    L.append("            break")
    L.append("        fi")
    L.append('        if [ "${_launched}" -eq 0 ] && [ "${_next}" -le "${NSAMPLES}" ]; then')
    L.append("            sleep 1")
    L.append("        fi")
    L.append("    done")
    L.append("fi")
    L.append("")

    # ---- 6) per-sample status sweep -------------------------------------------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE STATUS REPORT")
    L.append("#   Each sample's own completion markers determine the outcome:")
    L.append("#     *_error_log.txt non-empty  -> cellranger-atac count left tmp dirs on failure")
    L.append("#     no outs/summary.csv        -> cellranger-atac count never completed")
    L.append("# ============================================================")
    L.append("_overall_rc=0")
    _dirs_bash = " ".join('"' + d + '"' for d in run_names)
    L.append("for _d in " + _dirs_bash + "; do")
    L.append('    if [ -s "${Server_folder}/${_d}_error_log.txt" ]; then')
    L.append('        echo "[FAIL] ${_d}  (see ${Server_folder}/${_d}_error_log.txt)" >&2')
    L.append("        _overall_rc=1")
    L.append('    elif [ ! -f "${Server_folder}/${_d}/outs/summary.csv" ] && [ ! -f "${Server_folder}/${_d}/outs/metrics_summary.csv" ]; then')
    L.append('        echo "[FAIL] ${_d}  (no summary output — cellranger-atac count did not complete)" >&2')
    L.append("        _overall_rc=1")
    L.append("    else")
    L.append('        echo "[ OK ] ${_d}"')
    L.append("    fi")
    L.append("done")
    L.append("")

    # ---- 7) merged summary metrics (one all_summary.csv, Samplename column first) --------------
    L.append("# ============================================================")
    L.append("# MERGED SUMMARY METRICS  (one all_summary.csv, Samplename column first)")
    L.append("#   Gathers each sample's outs/summary.csv into a single CSV so all")
    L.append("#   samples' key metrics can be read in one table.")
    L.append("#   (created BEFORE the archive rsync so it lands in the central copy too)")
    L.append("# ============================================================")
    L.append('SM_TMP="${Summaries_folder}/all_summary_tmp.csv"')
    L.append('SM_OUT="${Summaries_folder}/all_summary.csv"')
    L.append(': > "${SM_TMP}"')
    L.append('_sm_hdr_written=0')
    for _i, _s in enumerate(samples, 1):
        _dn = run_names[_i - 1]
        L.append('    _sm_src="${Server_folder}/' + _dn + '/outs/summary.csv"')
        L.append('    [ -f "${_sm_src}" ] || _sm_src="${Server_folder}/' + _dn + '/outs/metrics_summary.csv"')
        L.append('    if [ -f "${_sm_src}" ]; then')
        L.append('        _sm_hdr=$(head -n1 "${_sm_src}" | tr -d "\\t")')
        L.append('        _sm_dat=$(tail -n +2 "${_sm_src}" | head -n1 | tr -d "\\t")')
        L.append('        if [ "${_sm_hdr_written}" -eq 0 ] && [ -n "${_sm_hdr}" ]; then')
        L.append('            echo "Samplename,${_sm_hdr}" > "${SM_TMP}"')
        L.append('            _sm_hdr_written=1')
        L.append('        fi')
        L.append('        if [ -n "${_sm_dat}" ]; then echo "' + _dn + ',${_sm_dat}" >> "${SM_TMP}"; fi')
        L.append("    fi")
    L.append('if [ -s "${SM_TMP}" ]; then')
    L.append('    cp "${SM_TMP}" "${SM_OUT}"')
    L.append('    rm -f "${SM_TMP}"')
    L.append('    echo "Merged summary metrics: ${SM_OUT}"')
    L.append("else")
    L.append('    rm -f "${SM_TMP}"')
    L.append('    echo "  [warn] no per-sample summary.csv found; all_summary.csv not created" >&2')
    L.append("fi")
    L.append("")

    # ---- 8) rsync working folder to the central archive -----------------------------------------
    L.append("# ============================================================")
    L.append("# DEPLOY  (rsync the whole working folder to the central archive)")
    L.append("# ============================================================")
    L.append('if command -v rsync >/dev/null 2>&1; then')
    L.append('    rsync -a --progress "${Server_folder}/" "${archive_dir}"')
    L.append('else')
    L.append('    cp -r "${Server_folder}"/* "${archive_dir}"')
    L.append("fi")
    L.append("")

    # ---- 9) whole-folder md5sums.txt (single file) -----------------------------------------------
    L.append("# ============================================================")
    L.append("# MD5 CHECKSUMS  (WHOLE FOLDER)")
    L.append("#   One md5sums.txt covering every file under ${Server_folder},")
    L.append("#   relative ./ paths.")
    L.append("# ============================================================")
    L.append('_md5_out="${Server_folder}/md5sums.txt"')
    L.append('( cd "${Server_folder}" && find . -type f ! -name "md5sums.txt" | sort | xargs md5sum > "md5sums.txt" )')
    L.append('echo "MD5 checksums (whole folder): ${_md5_out}"')
    L.append('if command -v rsync >/dev/null 2>&1; then rsync -avh "${_md5_out}" "${archive_dir}/" 2>/dev/null || true; fi')
    L.append("")

    # ---- 10) copy the .sh used to generate this project ------------------------------------------
    L.append("# ============================================================")
    L.append("# KEEP A COPY OF THIS SCRIPT")
    L.append("#   Copy the .sh used to generate this project into the project folders")
    L.append("#   (local working dir + central archive) so the exact pipeline is")
    L.append("#   preserved with the results.")
    L.append("# ============================================================")
    # Self-resolve THIS script to an absolute path right where it is used —
    # this far down the script nothing has `cd`'d in the main shell (the
    # per-sample functions call cellranger-atac with --output-dir and the
    # md5sum step runs in a subshell), so a relative BASH_SOURCE[0] still
    # means "launched from here".  Kept here (not at the top) because it is
    # plumbing, not a user-editable parameter.  NOTE: if a bare `cd` is ever
    # added to the main flow ABOVE this point, the relative name breaks —
    # then move this capture above that `cd` instead.
    L.append('# Path of this script as an absolute path, resolved against the launch directory')
    L.append('_SELF_PATH="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)/$(basename "${BASH_SOURCE[0]:-$0}")"')
    L.append('if [ ! -f "${_SELF_PATH}" ]; then _SELF_PATH="${BASH_SOURCE[0]:-$0}"; fi')
    L.append('if [ -n "${_SELF_PATH:-}" ] && [ -f "${_SELF_PATH}" ]; then')
    L.append('    cp "${_SELF_PATH}" "${Server_folder}/$(basename "${_SELF_PATH}")" 2>/dev/null || true')
    L.append('    cp "${_SELF_PATH}" "${archive_dir}/$(basename "${_SELF_PATH}")" 2>/dev/null || true')
    L.append("fi")
    L.append("")

    # ---- 11) keepdata cleanup --------------------------------------------------------------------
    L.append("# ============================================================")
    L.append("# LOCAL WORKDIR CLEANUP")
    L.append('#   After deploy to the archive, the local working folder is deleted')
    L.append('#   UNLESS keepdata="yes".  (MD5 checksums and the .sh copy are done first.)')
    L.append("# ============================================================")
    L.append('if [[ "${keepdata}" != "yes" ]]; then')
    L.append('    rm -r "${Server_folder}" 2>/dev/null || true')
    L.append("fi")
    L.append("")

    # ---- 12) park processed datasheets in the team "TO delete" folder ----------------------------
    L.append("# ============================================================")
    L.append("# DATASHEET CLEANUP")
    L.append("#   Park the Excel datasheets USED TO GENERATE THIS SCRIPT for later deletion —")
    L.append("#   only these exact files; any other .xlsx in the working dir is left untouched.")
    L.append("#   Primary:  " + DELETE_PRIMARY)
    L.append("#   Fallback: " + DELETE_FALLBACK)
    L.append("# ============================================================")
    L.append('if [ -d "' + DELETE_PRIMARY + '" ]; then _delfolder="' + DELETE_PRIMARY + '"')
    L.append('elif [ -d "' + DELETE_FALLBACK + '" ]; then _delfolder="' + DELETE_FALLBACK + '"')
    L.append("fi")
    L.append('if [ -n "${_delfolder:-}" ]; then')
    _srcs_abs = list(dict.fromkeys(os.path.abspath(s) for s in (source_files or [])))
    for _s in _srcs_abs:
        L.append(f'    mv "{_s}" "${{_delfolder}}" 2>/dev/null || true')
    L.append("fi")
    L.append("")

    # ---- 13) exit --------------------------------------------------------------------------------
    L.append('exit "${_overall_rc}"')

    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parses one or more AVT submission Excel files and generates a "
            "cellranger-atac_<AVT>.sh script (multi-run inputs: "
            "cellranger-atac_AVT0001_AVT0002.sh) with one `cellranger-atac "
            "count` per sample (ATAC only — no Multiome grouping, no TSO%)."
        ),
        epilog=(
            "Example:\n"
            "  python3 %(prog)s samples.xlsx\n"
            "  python3 %(prog)s -n run1.xlsx run2.xlsx run3.xlsx\n"
            "#   -> one merged script, named after every run it contains:\n"
            "#   cellranger-atac_AVT0001_AVT0002.sh"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"generate_cellranger_template_atac.py {SCRIPT_VERSION}",
        help="Print the generator version and exit.",
    )
    parser.add_argument(
        "xlsx_files", nargs="+",
        help="One or more AVT submission Excel files.",
    )
    parser.add_argument(
        "--non-interactive", "-n", action="store_true", default=False,
        help="Skip prompts; use the detected run name (AVT/MSH/MSQ/NVQ code).",
    )
    parser.add_argument(
        "--genome", choices=["human", "mouse"], default="human",
        help="Reference genome for non-interactive runs (default: human). "
             "Interactive runs always prompt.",
    )
    parser.add_argument(
        "--no-number-prefix",
        action="store_true",
        default=False,
        help=(
            "Do not number the per-sample run folders (ON by default): "
            "cellranger-atac count --id/--output-dir become '<sample>' "
            "instead of '01_<sample>', and the folder keeps the bare "
            "sample name.  Numbering keeps the folders in sample order "
            "even when the names are not alphabetical."
        ),
    )
    args = parser.parse_args()

    set_non_interactive(args.non_interactive)
    default_genome = {
        "human": HUMAN_GENOME_PATH,
        "mouse": MOUSE_GENOME_PATH,
    }[args.genome]
    setup_readline_completer()

    all_pairs, detected_avts = [], []
    per_file_orders = []  # per-file row order — feeds merge_sample_order
    any_user = ""
    source_files = []  # accepted Excel file(s) — parked in TO_DELETE by the generated .sh
    for path in args.xlsx_files:
        if not os.path.isfile(path) or not path.lower().endswith((".xlsx", ".xlsm")):
            print(f"Warning: skipping '{path}' (not found or not .xlsx).")
            continue
        source_files.append(path)
        rows, avt = parse_xlsx(path)
        if rows:
            # Per-file run inheritance (GEX parity): rows that name no run
            # take the run THIS FILE names (row-level names win) — instead
            # of one representative run stamped onto every row.
            all_pairs.extend(inherit_file_avt(rows, avt))
            # Row order of the samples in this datasheet — feeds the
            # cross-file order merge in merge_sample_order() below.
            per_file_orders.append([si for si, _, _ in rows])
        if avt:
            detected_avts.append(avt)
        if not any_user:
            any_user = read_user_field(path)

    if not all_pairs:
        sys.exit("Error: no sample rows could be parsed from the provided file(s).")

    # Cross-file sample order: the common order of all files when they
    # agree, else the FIRST file's order (later files' new samples are
    # appended after it).  See merge_sample_order.
    merged_sample_order = merge_sample_order(per_file_orders)

    print()
    print_pairs_table(all_pairs)
    print()

    # Run name (largest numeric suffix wins; prompt to confirm/override)
    main_avt = latest_avt(detected_avts)
    if not main_avt:
        main_avt = ask_avt_run(None)
    else:
        main_avt = ask_avt_run(main_avt)
    if not main_avt:
        sys.exit("Error: run name cannot be empty.")

    # GEX parity (shared helper): rows whose file named no run get ONE
    # explicit prompt (hint = first detected run) instead of silently
    # inheriting the representative run.  Rows that already carry a run —
    # stamped per-file by inherit_file_avt() or named on the row — are
    # left untouched.  Skipped when NO run is detected anywhere: main_avt
    # above already holds the answer and the 'or main_avt' fallbacks cover
    # these rows (single prompt, like GEX in that case).
    resolve_unassigned_avt(all_pairs)

    # ── Parity item 3: sample selection (keep/drop + add sheet-less) ────────
    # Same shared helper + prompt as the GEX/ARC generators.
    # ``ask_sample_selection`` wants ``{sample: [(fastq_id, avt_run)]}``;
    # atac keys it by the (un-suffixed) sample name from the sheet.  Returns
    # (kept_names, added) where added maps a sample name to the FULL
    # avid-folder path(s) the operator typed — those are self-referential and
    # used verbatim (no AVT-tree walk).  Under -n the helper returns every
    # sample kept and adds none (identical to pre-parity behaviour).
    sample_fastqs = OrderedDict()
    for _sample, _fid, _ravt in all_pairs:
        sample_fastqs.setdefault(_sample, []).append(
            (_fid or "", (_ravt or main_avt)))
    kept_names, added_raw = ask_sample_selection(sample_fastqs,
                                                 default_root=NGS_RUNS_ROOT)
    if dropped := [sn for sn in sample_fastqs if sn not in kept_names]:
        print(f"  Dropped from script: {', '.join(dropped)}")
    all_pairs = [p for p in all_pairs if p[0] in kept_names]
    # Added (sheet-less) samples: carry a full avid-folder path.  The row's
    # fastq id is the folder's basename (e.g. avid13144); the generator uses
    # the full path verbatim.
    added_paths = {}   # {fastq id (avid basename): full avid-folder path}
    for name, folders in added_raw.items():
        for p in folders:
            added_paths[os.path.basename(p)] = p
            # Per-row avt = None: a folder added by full path keeps the run
            # implied by the folder the operator typed.
            all_pairs.append((name, os.path.basename(p), None))
    if not all_pairs:
        sys.exit("Error: no samples selected — nothing to generate.")

    # ── Per-row AVT — resolve each FASTQ under its own run tree ─────────────
    # Each sample's FASTQ is resolved under the AVT run named on its own sheet
    # row, so a worksheet mixing two runs resolves correctly.  A row that
    # names no run inherits the run named by its own file; rows whose file
    # names no run were explicitly prompted for one first
    # (resolve_unassigned_avt, GEX parity).
    fid_avt = {}
    for _sample, _fid, _ravt in all_pairs:
        _fid = (_fid or "").strip()
        if _fid:
            fid_avt.setdefault(_fid, _ravt or main_avt)
    _distinct_avts = sorted(set(fid_avt.values()))
    if len(_distinct_avts) > 1:
        print(f"  Note: this workspace spans AVT runs {', '.join(_distinct_avts)}; "
              f"each FASTQ resolves under the run named on its own row.")

    print()
    print_pairs_table(all_pairs)
    print()

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=default_genome,
        human_label="Human ATAC", mouse_label="Mouse ATAC")
    user = any_user or ""
    run_dir = run_base_dir(main_avt)  # noqa: F841 (kept for reference)

    # ── Parity item 4: run-folder disambiguation (AVT0551 case) ─────────────
    # Before: find_fastq_path() blind-walked the run tree and took the FIRST
    # hit — the AVT0551 case, where the same sample recurs under
    # 1_wrong_index_… AND 2_corrected_indexes_… and silently resolved to the
    # wrong library.  Now we pin the SAMPLE's directory via list_sample_dirs/
    # choose_sample_dir before resolving; choose_sample_dir prompts only when
    # >1 candidate exists (the 99% single-candidate case resolves silently)
    # and also no-ops under -n (first candidate + a warning naming the
    # alternatives).  Samples added by full path have no candidate set —
    # they are used verbatim.
    sample_names = list(dict.fromkeys(p[0] for p in all_pairs))
    chosen_dirs = {}
    for sn in sample_names:
        _sn_fid = next((p[1] for p in all_pairs
                        if p[0] == sn and (p[1] or "")), "")
        _sn_avt = fid_avt.get(_sn_fid, main_avt) or main_avt
        cands = list_sample_dirs(_sn_avt, sn)
        if len(cands) > 1:
            chosen = choose_sample_dir(sn, cands)
            if chosen:
                chosen_dirs[sn] = chosen
    if chosen_dirs:
        print(f"  Run-folder disambiguation pinned {len(chosen_dirs)} sample(s).")

    # Emission-time resolver: build_script_full calls resolve(sample, fid).
    # Added-by-path samples come back verbatim; the rest honour the chosen
    # sample dir first (find_fastq_path_with_pref) and fall back to the
    # classic full-tree walk under the FASTQ's own run tree.
    def resolve(sample, fid):
        _fid = (fid or "").strip()
        if not _fid:
            return None
        if _fid in added_paths:            # added-by-path sample → verbatim
            return added_paths[_fid]
        _f_avt = fid_avt.get(_fid, main_avt)
        pref = chosen_dirs.get(sample)
        if pref:                          # preferred sample dir first
            hit = find_fastq_path_with_pref(_f_avt, _fid,
                                            preferred_sample_dir=pref)
            if hit:
                return hit
        return find_fastq_path(_f_avt, _fid)   # classic full-tree walk

    # Multi-run inputs: the script name and the working/archive folder
    # enumerate every AVT run the kept FASTQs resolve under (sorted
    # numerically), e.g. cellranger-atac_AVT0001_AVT0002.sh and
    # AVT0001_AVT0002_user_$TIMESTAMP; a single-run input keeps the plain
    # name.
    avt_tag = avt_run_tag(set(fid_avt.values())) or main_avt
    script = build_script_full(
        reference_genome, all_pairs, user=user, main_avt=main_avt,
        avt_tag=avt_tag, source_files=source_files, resolve=resolve,
        number_prefix=not args.no_number_prefix,
        ordered_samples=merged_sample_order)
    out_script = f"cellranger-atac_{avt_tag}.sh"
    with open(out_script, "w") as f:
        f.write(script)
    os.chmod(out_script, 0o755)
    print(f"Generated: {out_script}")
    # NB: match only actual flagged variables (fastq_<i>_<j>="" # FIXME_LOCATE_),
    # NOT the static FIXME_LOCATE_ mention in the FASTQ PATHS help comment.
    n_unresolved = len(re.findall(
        r'fastq_\d+_\d+=""\s*#\s*FIXME_LOCATE_', script))
    if n_unresolved:
        print(f"Note: {n_unresolved} fastq leaf(s) could not be resolved here — "
              "replace the empty fastq_<idx>_<j> variables flagged "
              "FIXME_LOCATE_ in the script before running.")


if __name__ == "__main__":
    main()
