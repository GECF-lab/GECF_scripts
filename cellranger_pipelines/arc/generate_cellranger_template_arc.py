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
        --create-bam=true --localcores=32 --localmem=126 \\
        --cell-annotation-model=auto \\
        --tenx-cloud-token-path=/home/gecf/.config/txg/credentials

Fastq leaf directories are resolved at generation time under
``<run_dir>/fastq_bases2fastq/<batch>/Samples/<experiment>/<fastq_id>/``
(or — when a run has no ``fastq_bases2fastq`` level — directly under
``<run_dir>/<batch>/Samples/<experiment>/<fastq_id>/``); when they
cannot be resolved (raw data not mounted here) a clearly marked
placeholder is emitted for the operator to fill in.
"""

import argparse
import os
import re
import subprocess
import sys
import time
import unicodedata
import warnings
from collections import OrderedDict, defaultdict

warnings.filterwarnings("ignore",
                        message="Data Validation extension is not supported")
warnings.filterwarnings("ignore",
                        message="Conditional Formatting extension is not supported")

try:
    import openpyxl
except ImportError:
    sys.exit("Error: the python module 'openpyxl' is required.\n"
             "Install it with: pip3 install openpyxl")

# Ensure the cellranger_pipelines/ package root is on sys.path so ``common.intake``
# resolves.  This file lives at ``cellranger_pipelines/arc/`` so its parent is the
# package root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_CELLRANGER_PIPELINES_ROOT = os.path.normpath(os.path.join(_HERE, os.pardir))
if _CELLRANGER_PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, _CELLRANGER_PIPELINES_ROOT)

# Shared intake atoms — see cellranger_pipelines/common/intake.py.  All the low-level
# workbook helpers (norm, merged-cell lookup, cell read, the 'User:' label
# scan) are imported from there and used by both this generator and the
# ATAC generator on an identical code path.
from common.intake import (
    NGS_BASE_DIR_TEMPLATE,
    AVT_PATTERN,
    FASTQ_ID_PATTERN,   # row gate: ^avid\d{4,6}$ — single source of truth
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
    # Same shared helpers the GEX generator uses for its sample-selection
    # prompt and run-folder disambiguation (single source of truth in
    # common.intake).
    ask_sample_selection,        # keep/drop + add sheet-less samples
    list_sample_dirs,            # candidate dirs for one (AVT, sample)
    choose_sample_dir,           # pick ONE when >1 candidate
    find_fastq_path_with_pref,   # resolve an avid id honouring the pick
)

# Provenance: version stamped into the generated .sh header
# (``scriptversion=...``) and reportable on the generator itself via
# ``-v``/``--version``.  Generator version == generated-script version;
# the .sh kept in the archive tells you which one produced it.
# Bump on any change to the generator.
#
# 1.3.1 = one TSO% section per GEX-leg run: every GEX fastq (Read-2) of a
#         sample is measured, not just the first leg; report sections are
#         labeled <NN_sample>__<avid>.
SCRIPT_VERSION = "1.3.1"

NGSruns_ARCHIVE_DIR = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/arc"

# Reference-genome root.  The genome folder inside each ``current_*``
# directory is never hard-coded here: each of them holds exactly one
# genome folder and that folder's name can change on the server
# (rename/replace), so the real path is resolved once at startup from
# the live filesystem (see the _resolve_reference_genome block below).
CELLRANGER_REFERENCES_ROOT = "/home/gecf/references/cellranger_references"
_HUMAN_CURRENT_DIR = os.path.join(CELLRANGER_REFERENCES_ROOT, "current_human_arc_ATAC")
_MOUSE_CURRENT_DIR = os.path.join(CELLRANGER_REFERENCES_ROOT, "current_mouse_arc_ATAC")


# ---------------------------------------------------------------------------
# Requester -> safe path-segment (parity item 12).  Ported from the finished
# GEX generator.  The requester string comes from the ``User:`` cell of the
# datasheet and is not validated; it is sanitised to a CIFS-/URL-safe segment
# before it is baked into the working-folder name, so an accented or
# punctuation-laden requester cannot break the path.
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# TSO% library-quality helper (cellranger_pipelines/common/TSO_percent.sh).
# The body of `tso_percent()` is inlined once per generated script so the
# analysis server needs no runtime dependency on this repo.  Ported from the
# finished GEX generator (parity item 9).
# ---------------------------------------------------------------------------
_TSO_HELPER_PATH = os.path.join(_CELLRANGER_PIPELINES_ROOT, "common", "TSO_percent.sh")


def _load_tso_helper_body() -> str:
    """Return the region of ``cellranger_pipelines/common/TSO_percent.sh`` up to
    the ``# %%EMBED-END%%`` marker (the ``tso_percent()`` function and its
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


def _is_version_check() -> bool:
    """True when the user is only asking for ``-v``/``--version`` so the
    import-time reference-genome scan can stay quiet (argparse hasn't run
    yet but ``sys.argv`` already carries the flag)."""
    return ("-v" in sys.argv or "--version" in sys.argv)


def _resolve_reference_genome(current_dir, label):
    """Return the genome path inside ``current_dir`` — its single subdirectory.

    Only *current_* is guaranteed to persist across the analysis servers;
    the genome folder it contains (refdata-cellranger-arc-...) can be
    renamed or replaced without the generator being updated.
    Non-recursive, so an inner ``old/`` tier never leaks in (only the
    current genome counts).
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

# (AVT_PATTERN is imported from common.intake — see the import block above.)

# Batch / control folders that are not real sample leaves.
EXCLUDED_BATCH_PATTERNS = (
    re.compile(r"phix", re.IGNORECASE),
    re.compile(r"\bunassigned\b", re.IGNORECASE),
    re.compile(r"for_check", re.IGNORECASE),
)
TENX_TOKEN_PATH = "/home/gecf/.config/txg/credentials"


# Analysis-server working root (scratch) where cellranger-arc output is
# written; emitted into the generated .sh's run-wide parameters block.
SCRATCH_ROOT = "/mnt/data/CellRangerCountOutput"

# GEX parity item 8c: park the Excel datasheet(s) used to generate the
# project .sh in a team "TO DELETE" folder for later cleanup.  Primary
# first, fallback if it isn't on this server; whichever directory exists
# at run time wins.
DELETE_PRIMARY = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CR_pipelines_Daniel/TO_DELETE"
)
DELETE_FALLBACK = (
    "/home/gecf/NGS_pipelines/CR_pipeline/CRarc_Multiome/TO_DELETE"
)

# GEX parity item 5: per-sample editable parameters emitted as one block
# per sample at the top of the generated .sh (mirrors gex._PER_SAMPLE_PARAMS,
# restricted to cellranger-arc's knobs: no include_introns/chemistry/
# expect_cells/force_cells/r1_length/r2_length).
#   * create_bam            : true | false — emit aligned fastq BAM
#   * cell_annotation_model : auto | <model name>, or blank (=auto)
#   * custom_options        : extra flags as-is, e.g. --nosecondary --localvmem=300
# A blank value ("") drops that flag for that sample.
# CPU/memory limits are NOT per-sample: localcores / localmem are run-wide
# parameters (see the RUN-WIDE PARAMETERS block in build_script).
_ARC_PER_SAMPLE_PARAMS = (
    ("create_bam",            "true",
     "true = emit aligned BAM (default); false = faster + smaller output"),
    ("cell_annotation_model", "auto",
     "auto | <cellranger-arc annotation model name>; blank = auto (default)"),
    ("custom_options",        "",
     "extra flags passed through as-is; e.g. --nosecondary --localvmem=300"),
    ("gex_exclude_introns",   "",
     "any non-blank value adds --gex-exclude-introns for this sample (skip GEX "
     "intronic reads); blank (default) = include introns.  Feeds the 'Include "
     "introns?' column of the stats sheet."),
)



# ─── input helpers ────────────────────────────────────────────────────


setup_readline_completer()  # tab-complete paths in the input() prompts


# (norm / get_cell_value / build_merged_cell_lookup are imported from
# common.intake — see the import block above.)



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


def _row_avt(ws, merged_lookup, row_idx, fastq_col, sample_col, fastq_txt):
    """Per-row AVT run-name detection (parity item 11).

    Mirrors the finished GEX generator: the AVT run name may sit in the
    fastq-ID cell itself or in any *adjacent* cell of the same table row.
    The fastq-ID column and the sample-name column are skipped (they carry
    the avid id / the sample label, not the run).  Returns the upper-cased
    run name (``AVT0XXX``) or ``None`` when this row does not name a run.
    """
    m = AVT_PATTERN.search(fastq_txt or "")
    if m:
        return m.group(0).upper()
    for c in range(1, ws.max_column + 1):
        if c == fastq_col or c == sample_col:
            continue
        val = get_cell_value(ws, merged_lookup, row_idx, c)
        if val is None:
            continue
        cm = AVT_PATTERN.search(str(val))
        if cm:
            return cm.group(0).upper()
    return None


def _sheet_sample_rows(ws, merged_lookup):
    """
    Extract (sample, fastq, avt) triples from one sheet, or [] if it is not a
    data sheet.  Each row's ``avt`` is the per-run AVT name (``AVT0XXX``)
    detected from the fastq-ID cell or an adjacent cell of that row, or
    ``None`` when the row does not name a run (parity item 11).
    """
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
            rows.append((sample_txt, fastq_txt,
                         _row_avt(ws, merged_lookup, row_idx,
                                  fastq_col, sample_col, fastq_txt)))
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
    ``(sample_name, fastq_id, avt_run)`` triples (real library ids only) and
    ``detected_avt`` is the AVT run name found in a data row, any sheet cell,
    or the file name (or None).

    Handles both the single-sheet layout of older submissions
    (e.g. AVT0519) and the 15-sheet GECF submission workbooks where the
    sample table is NOT the first sheet (first sheet is GUIDELINES).
    """
    wb = openpyxl.load_workbook(path, data_only=True)

    sheet_name, rows = best_data_sheet(wb)
    if rows:
        # De-duplicate, preserving first-seen order.  Rows are
        # (sample, fastq, avt) triples; the per-row avt rides along and is not
        # part of the duplicate key (parity item 11).
        seen = set()
        deduped = []
        for sample, fastq, avt in rows:
            key = (sample.lower(), fastq.lower())
            if key in seen:
                continue
            seen.add(key)
            deduped.append((sample, fastq, avt))
        rows = deduped

    # Sheet-wide fallback: rows that do not name a run inherit the AVT run
    # found in any data cell of the table (a run name written once in the
    # sheet), so per-sample resolution still has a value to hang on.
    sheet_avt = None
    for sample, fastq, avt in rows:
        m = (AVT_PATTERN.search(avt) if avt else None) \
            or AVT_PATTERN.search(sample) or AVT_PATTERN.search(fastq)
        if m:
            sheet_avt = m.group(0).upper()
            break
    if sheet_avt:
        rows = [(s, f, (a or sheet_avt)) for (s, f, a) in rows]

    # Final fallback: the file name.
    if sheet_avt is None:
        m = AVT_PATTERN.search(os.path.basename(path))
        if m:
            sheet_avt = m.group(0).upper()
            rows = [(s, f, (a or sheet_avt)) for (s, f, a) in rows]

    return rows, sheet_avt


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
    for row in rows:
        # Rows are (sample, fid, [avt]); the per-row avt rides separately via
        # the fid->avt map built in main(), so it is not stored here.
        sample, fid = row[0], row[1]
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
    folder.

    Two known demux layouts are accepted (mirroring the depths that
    ``common.intake.list_sample_dirs`` honours):

      1. ``<run_dir>/fastq_bases2fastq/<batch>/Samples/<experiment>/<fastq_id>``
         — the canonical layout (e.g. AVT0660);
      2. ``<run_dir>/<batch>/Samples/<experiment>/<fastq_id>``
         — demux output directly under the run folder, without the
         ``fastq_bases2fastq`` level (e.g. AVT0654); the flat
         ``<run_dir>/Samples/<experiment>/<fastq_id>`` form is preferred.

    Layout 2 is only searched when ``fastq_bases2fastq`` is absent, so
    runs with the canonical layout keep the exact previous behaviour.
    Batch folders matching EXCLUDED_BATCH_PATTERNS (phix / unassigned /
    for_check) are skipped.  Returns the absolute leaf path, or None if
    the id cannot be resolved.
    """
    if not fastq_id:
        return None

    def _in(samples_dir):
        if not os.path.isdir(samples_dir):
            return None
        for dirpath, dirnames, _ in os.walk(samples_dir):
            if fastq_id in dirnames:
                leaf = os.path.join(dirpath, fastq_id)
                if os.path.isdir(leaf):
                    return leaf
        return None

    strict_root = os.path.join(run_dir, "fastq_bases2fastq")
    if os.path.isdir(strict_root):
        # Canonical layout — unchanged from the previous behaviour.
        for batch in sorted(os.listdir(strict_root)):
            if not batch or batch.startswith(".") or is_excluded_batch(batch):
                continue
            hit = _in(os.path.join(strict_root, batch, "Samples"))
            if hit:
                return hit
        return None

    # No fastq_bases2fastq level: the demux output sits directly under the
    # run folder.  Try the flat form first, then each top-level batch folder
    # (same exclusions as the canonical layout).
    try:
        tops = sorted(os.listdir(run_dir))
    except OSError:
        return None
    hit = _in(os.path.join(run_dir, "Samples"))
    if hit:
        return hit
    for top in tops:
        if not top or top.startswith(".") or is_excluded_batch(top):
            continue
        hit = _in(os.path.join(run_dir, top, "Samples"))
        if hit:
            return hit
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


def _probe_arc_version():
    """Best-effort generation-time probe of ``cellranger-arc --version`` so the
    generated .sh carries a sane fallback for the 'Cell Ranger arc version'
    stats column.  At run time the .sh re-queries the binary and only falls
    back to this value.  Returns "" on any failure — never required."""
    try:
        proc = subprocess.run(["cellranger-arc", "--version"],
                              capture_output=True, text=True, timeout=30)
        raw = (proc.stdout or "") + " " + (proc.stderr or "")
        m = re.search(r"(\d+\.\d+(?:\.\d+)?)", raw)
        return ("v" + m.group(1)) if m else ""
    except Exception:
        return ""


# Embedded Python that assembles the two-sheet stats workbook (or a two-CSV
# fallback when openpyxl is unavailable).  Emitted verbatim into the .sh
# inside a quoted heredoc; kept free of backslash escapes so it survives
# here-doc processing intact.
_PY_STATS = [
    "import sys, os",
    "gex_tsv, atac_tsv, genome, arc_ver, flags, outxlsx = sys.argv[1:7]",
    "gex_hdr = ['Official sample name', 'GECF library label (L_xxx)', 'GECF library exp code', 'seq. run', 'fastq code GEX (nxid, nvid)', 'R2 length', 'Cell Ranger arc version', 'genome assembly*', 'Include introns?']",
    "atac_hdr = ['Official sample name', 'GECF library label (L_xxx)', 'GECF library exp code', 'seq. Run ATAC', 'fastq code ATAC (nxid, nvid)', 'R2 length ATAC', 'Cell Ranger arc version', 'genome assembly*', 'Include introns?']",
    "sep = chr(9)",
    "def _rows(path):",
    "    if not os.path.isfile(path):",
    "        return []",
    "    out = []",
    "    with open(path) as fh:",
    "        for ln in fh:",
    "            ln = ln.rstrip(chr(10))",
    "            if not ln.strip():",
    "                continue",
    "            p = ln.split(sep)",
    "            if len(p) < 6:",
    "                continue",
    "            while len(p) < 7:",
    "                p.append('')",
    "            i = int(p[0])",
    "            inc = 'no' if (len(flags) >= i and flags[i - 1] == 'n') else 'yes'",
    "            out.append([p[1], p[2], p[3], p[4], p[5], p[6], arc_ver, genome, inc])",
    "    return out",
    "gexrows = _rows(gex_tsv)",
    "atarows = _rows(atac_tsv)",
    "def _write_csv(path, hdr, rows):",
    "    import csv",
    "    with open(path, 'w', newline='') as fh:",
    "        w = csv.writer(fh)",
    "        w.writerow(hdr)",
    "        w.writerows(rows)",
    "try:",
    "    import openpyxl",
    "    wb = openpyxl.Workbook()",
    "    ws = wb.active",
    "    ws.title = 'GEX'",
    "    ws.append(gex_hdr)",
    "    for r in gexrows:",
    "        ws.append(r)",
    "    ws2 = wb.create_sheet('ATAC')",
    "    ws2.append(atac_hdr)",
    "    for r in atarows:",
    "        ws2.append(r)",
    "    wb.save(outxlsx)",
    "    print('STATS: wrote ' + outxlsx + '  (GEX %d row(s), ATAC %d row(s))' % (len(gexrows), len(atarows)))",
    "except Exception as e:",
    "    gex_csv = os.path.splitext(outxlsx)[0] + '_GEX.csv'",
    "    atac_csv = os.path.splitext(outxlsx)[0] + '_ATAC.csv'",
    "    _write_csv(gex_csv, gex_hdr, gexrows)",
    "    _write_csv(atac_csv, atac_hdr, atarows)",
    "    print('STATS: WARNING openpyxl unavailable (%s) - wrote a CSV fallback instead of the xlsx:' % (e,))",
    "    print('       ' + gex_csv)",
    "    print('       ' + atac_csv)",
]


def build_stats_block(groups, fid_avt, avt_run, reference_genome):
    """Emit the 'STATS' section of the generated .sh.

    Produces a two-sheet ``<avt>_stats.xlsx`` (GEX + ATAC) in ``Server_folder``;
    one row per sample per modality with columns:

      official sample name | GECF library label | GECF library exp code
      seq. run             | fastq code         | R2 length
      cellranger-arc version*  | genome assembly* | Include introns?

    (*resolved at run time on the analysis host; the first six columns are
    baked in at generation time from the submission sheet).  The GECF library
    exp code and R2 length are left blank — they are not present in the run
    data (LPX codes live only in the deleted submission xlsx, and R2 length is
    derived from the manifest mask, not the CR output).

    Non-fatal: without python3 the step warns and skips; with python3 but
    without openpyxl it writes a two-CSV fallback.  Either way the run
    proceeds, and the workbook (or CSVs) is md5'd + rsync'd to the archive.
    """
    fid_avt = fid_avt or {}
    gen_arc = _probe_arc_version()
    gen_genome = os.path.basename((reference_genome or "").rstrip("/"))
    n = len(groups)

    # ---- static rows (known at generation time) ---------------------------------------------
    # Row = [sample_idx, official, label, exp_code, seq_run, fastq, r2_length]
    gex_rows, atac_rows = [], []
    gcnt = acnt = 0
    for idx, (_base, libs) in enumerate(groups.items(), start=1):
        for full, fid, modality in libs:
            if modality == _GEX:
                gcnt += 1
                gex_rows.append([idx, full, "GEX_L_%d" % gcnt, "",
                                 fid_avt.get(fid, avt_run) or "", fid, ""])
                break
        for full, fid, modality in libs:
            if modality == _ATAC:
                acnt += 1
                atac_rows.append([idx, full, "AT_L_%d" % acnt, "",
                                  fid_avt.get(fid, avt_run) or "", fid, ""])
                break

    def _tsv(rows):
        return "\n".join("\t".join(str(c) for c in r) for r in rows)

    out_name = "%s_stats.xlsx" % avt_run
    L = []
    L.append("# ============================================================")
    L.append("# STATS -- two-sheet %s  (one row per sample per modality)" % out_name)
    L.append("#   official name | GECF library label | exp code | seq. run | fastq | R2 length")
    L.append("#   + (run-time)  cellranger-arc version | genome assembly | Include introns?")
    L.append("#   Baked-in columns from the submission sheet; the last three resolved at run")
    L.append("#   time.  Written into ${Server_folder} so md5 + rsync below ship it.")
    L.append("#   Non-fatal: no python3 -> warn + skip; no openpyxl -> two CSVs.  Run proceeds.")
    L.append("# ============================================================")
    L.append('STATS_OUT="${Server_folder}/%s"' % out_name)

    if gex_rows:
        L.append("cat > \"${Server_folder}/.stats_gex.tsv\" <<'STATS_GEX_EOF'")
        L.extend(_tsv(gex_rows).split("\n"))
        L.append("STATS_GEX_EOF")
    else:
        L.append(': > "${Server_folder}/.stats_gex.tsv"')
    if atac_rows:
        L.append("cat > \"${Server_folder}/.stats_atac.tsv\" <<'STATS_ATAC_EOF'")
        L.extend(_tsv(atac_rows).split("\n"))
        L.append("STATS_ATAC_EOF")
    else:
        L.append(': > "${Server_folder}/.stats_atac.tsv"')

    # ---- runtime-resolved values -----------------------------------------------------------
    L.append("_ARC_VERSION=\"\"")
    L.append("if command -v cellranger-arc >/dev/null 2>&1; then")
    L.append("  _av=\"$(cellranger-arc --version 2>/dev/null)\"")
    L.append("  _v=\"\"")
    L.append("  for _t in ${_av}; do")
    L.append("    case \"${_t}\" in")
    L.append("      *.* ) _d=\"${_t//[!0-9]/}\"; if [ -n \"${_d}\" ]; then _v=\"${_t}\"; break; fi ;;")
    L.append("    esac")
    L.append("  done")
    L.append("  [ -n \"${_v}\" ] && _ARC_VERSION=\"v${_v}\"")
    L.append("fi")
    L.append("[ -z \"${_ARC_VERSION}\" ] && _ARC_VERSION=\"%s\"" % (gen_arc or ""))
    L.append("_GENOME_ASSEMBLY=\"$(basename \"${REFPATH}\" 2>/dev/null)\"")
    L.append("[ -z \"${_GENOME_ASSEMBLY}\" ] && _GENOME_ASSEMBLY=\"%s\"" % (gen_genome or ""))
    L.append("_INTRON_FLAGS=\"\"")
    _intr_line = '  if [ -n "${gex_exclude_introns_%d:-}" ]; then _INTRON_FLAGS="${_INTRON_FLAGS}n"; else _INTRON_FLAGS="${_INTRON_FLAGS}y"; fi'
    for _i in range(1, n + 1):
        L.append(_intr_line % _i)

    # ---- write the workbook (xlsx preferred; CSV fallback) --------------------------------
    L.append("if command -v python3 >/dev/null 2>&1; then")
    L.append('  echo "STATS: building ${STATS_OUT} ..."')
    L.append("  python3 - \"${Server_folder}/.stats_gex.tsv\" \"${Server_folder}/.stats_atac.tsv\" \"${_GENOME_ASSEMBLY}\" \"${_ARC_VERSION}\" \"${_INTRON_FLAGS}\" \"${STATS_OUT}\" <<'STATS_PY_EOF' || echo \"  [warn] STATS: could not write the stats workbook (see above); run continues.\" >&2")
    L.extend(_PY_STATS)
    L.append("STATS_PY_EOF")
    L.append("else")
    L.append("  echo '  [warn] STATS: python3 not found on this host - the stats workbook' >&2")
    L.append("  echo '         was not generated.  Install python3 + openpyxl on the analysis server.' >&2")
    L.append("fi")

    # ---- ship whichever artifact was produced to the archive --------------------------------
    L.append("if command -v rsync >/dev/null 2>&1; then")
    L.append("  for _f in \"${Server_folder}/%s_stats.xlsx\" \"${Server_folder}/%s_stats_GEX.csv\" \"${Server_folder}/%s_stats_ATAC.csv\"; do" % (avt_run, avt_run, avt_run))
    L.append("    [ -f \"${_f}\" ] && rsync -avh \"${_f}\" \"${archive_dir}/\" 2>/dev/null || true")
    L.append("  done")
    L.append("fi")
    L.append('rm -f "${Server_folder}/.stats_gex.tsv" "${Server_folder}/.stats_atac.tsv"')
    L.append("")
    return L


def build_script(reference_genome, groups, avt_run, user, run_dir,
                 source_files=None, resolve=None, fid_avt=None,
                 number_prefix=True):
    """Emit the ``cellranger-arc_AVT0XXX.sh`` text.

    ``resolve`` maps ``(sample_base, fastq_id)`` to a leaf directory (or
    None); resolved at generation time and baked into the .sh as literal
    path assignments.  When ``resolve`` is None the classic full-tree walk
    (``resolve_fastq_leaf``) is used directly.
    *number_prefix* (default True): zero-pad the sample index (01_, 02_,
    …) onto each per-sample run folder (--id / --output-dir) and its
    companion files, keeping the folders in sample order even when the
    names are not alphabetical.  Pass False (CLI: --no-number-prefix)
    for bare sample names (the previous behaviour).


    Layout (top -> bottom) — items 5–8 of the arc/GEX parity list:
      1.  Header / constants (timestamp, reference, source files, requester,
          scriptversion stamp, REFPATH, avt_run, user).
      2.  RUN-WIDE PARAMETERS (editable at the top): keepdata, max_parallel,
          localcores, localmem, scratch_root, server_folder_name,
          NGSruns_archive, then all FASTQ path variables (resolved at
          generation time).
      3.  PER-SAMPLE PARAMETERS (item 5 — one block per sample): create_bam,
          cell_annotation_model, custom_options, gex_exclude_introns;
          blank a value to drop that flag for that sample.
          All user-editable blocks (2 + 3) stay contiguous at the top of the
          script; the non-editable housekeeping sits just below them.
      3b. PIPELINE PLUMBING: reference-genome existence check (fails fast,
          before any mkdir), then derived vars (Server_folder,
          metrics_folder, archive_dir) and the output-folder mkdirs.
      4.  One ``arc_sample_<i>()`` function per sample: builds the libraries
          CSV, runs ``cellranger-arc count`` (flags from the ``*_<i>``
          variables in block 3), checks for leftover tmp folders, harvests
          web summaries (joint + GEX_only + ATAC_only), deletes the
          SC_ATAC_GEX_COUNTER_CS scratch tree, rsync's the sample dir and
          auxiliary files to the archive.
      5.  DISPATCH (item 6): serial (max_parallel=1) or rolling pool of
          max_parallel concurrent arc samples.
      6.  Per-sample status sweep (item 7): [OK]/[FAIL] per sample,
          ``_overall_rc`` accumulator.
      7.  rsync ``raw_metrics_folder`` to the archive (global — arc-specific).
      8.  Whole-folder ``md5sums.txt`` (item 8a; replaces per-sample
          ``*_md5_checksums.txt``).
      9.  ``.sh`` copy to work folder + archive (item 8b).
      10. keepdata cleanup.
      11. TO_DELETE parking of the source .xlsx files (item 8c).
      12. ``exit ${_overall_rc}``.
    """
    _resolve = (resolve if resolve is not None
                else (lambda base, fid: resolve_fastq_leaf(run_dir, fid)))
    # Per-sample run folder name: 01_, 02_, … (ON by default) keeps the
    # folders in sample order even when the names are not alphabetical.
    # --no-number-prefix -> bare sample name (previous behaviour).
    run_names = [f"{i:02d}_{b}" if number_prefix else b
                 for i, b in enumerate(groups, 1)]
    L = []

    # ---- 1) header / constants -----------------------------------------------------------------
    L.append("#!/bin/bash")
    L.append("# cellranger-arc — one `cellranger-arc count` per biologic sample (GEX + ATAC Multiome).")
    L.append(f"# Generated by generate_cellranger_template_arc.py ({len(groups)} sample(s)).")
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
    L.append('REFPATH="' + reference_genome + '"')
    L.append(f"avt_run='{avt_run}'")
    # Item 12: the requester string (datasheet ``User:`` cell) is not
    # validated; clean it to a CIFS-/URL-safe path segment HERE at generation
    # time so an accented / punctuation-laden name cannot later break the
    # working-folder name.  If the result is empty, the ``${user:+_${user}}``
    # conditional below simply omits the ``_user`` segment.
    L.append("user='" + folder_safe(user or "") + "'")
    L.append("")

    # ---- 2) RUN-WIDE PARAMETERS (editable block) -----------------------------------------------
    L.append("# ============================================================")
    L.append("# RUN-WIDE PARAMETERS  (defaults — edit here)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   keepdata : "yes" = keep the local working folder after the run;')
    L.append('#               "no"  = delete it (data already rsynced to the archive).')
    L.append('keepdata="no"')
    L.append('#   max_parallel : how many cellranger-arc counts to run at the same time.')
    L.append('#   1 = serial (one sample after another); N>1 = rolling pool of N jobs,')
    L.append('#         the next sample starts as soon as a slot frees up.  Default is 3,')
    L.append('#         i.e. 3 concurrent runs.')
    L.append('max_parallel=3')
    L.append('#   localcores : CPU cores per cellranger-arc count (applies to ALL samples);')
    L.append('#                blank = cellranger default.')
    L.append('#   localmem   : RAM in GB per cellranger-arc count (applies to ALL samples);')
    L.append('#                blank = cellranger default.')
    L.append('#   With max_parallel>1, total = max_parallel x (localcores cores, localmem GB),')
    L.append('#   so all three must fit within the server\'s resources together.')
    L.append(f'localcores="{LOCAL_CORES}"')
    L.append(f'localmem="{LOCAL_MEM}"')
    L.append('#   scratch_root : root on the analysis server where cellranger-arc output is written')
    L.append('#                  (fixed convention — change only if the mount differs on your server).')
    L.append('scratch_root="' + SCRATCH_ROOT + '"')
    L.append('#   server_folder_name : working + archive leaf folder name.')
    L.append('#                  default is <AVT>[_user]_<TIMESTAMP>; change only if needed.')
    L.append('server_folder_name="${avt_run}${user:+_${user}}_$TIMESTAMP"')
    L.append('#   NGSruns_archive: central archive destination (fixed target).')
    L.append('NGSruns_archive="' + NGSruns_ARCHIVE_DIR + '"')
    L.append("")

    # Fastq paths — one variable per library, resolved at generation time
    L.append("# ============================================================")
    L.append("# FASTQ PATHS  (one per sample library — resolved at generation time)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   If a variable is set to the empty string (""), the generator could not')
    L.append('#   resolve the leaf directory — set the real path by hand before running.')
    any_missing = False
    for idx, (base, libs) in enumerate(groups.items(), start=1):
        ordered = sorted(libs, key=lambda t: 0 if t[2] == _GEX else 1)
        L.append("")
        L.append(f"# --- Sample {idx}: {base}  ({len(ordered)} library(ies)) ---")
        for full, fid, modality in ordered:
            if not fid:
                var = f"fastq_{idx}_{modality.replace(' ', '_')}"
                L.append(f'{var}=""   # (no fastq ID parsed — set the leaf path by hand)')
                any_missing = True
                continue
            var = f"fastq_{idx}_{fid}"
            leaf = _resolve(base, fid)
            if leaf:
                L.append(f'{var}="{leaf}"')
            else:
                L.append(f'{var}=""   # FIXME: leaf for \'{full}\' (fid {fid}) — set by hand before running')
                any_missing = True
    if any_missing:
        L.append("")
        L.append("# NOTE: at least one fastq variable is empty — set it before running.")
        print("  Warning: some fastq leaves could not be resolved at generation time.",
              file=sys.stderr)
        print("  The flagged fastq_<idx>_<fid> variables must be set before running.",
              file=sys.stderr)
    L.append("")

    # ---- 3) PER-SAMPLE PARAMETERS (item 5) ------------------------------------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE PARAMETERS  (one block per sample — edit values in place)")
    L.append("# ----------------------------------------------------------------")
    L.append('#   * blank a value ( value="" )  -> that flag is DROPPED for that sample')
    L.append("#   * create_bam             : true | false")
    L.append("#   * cell_annotation_model  : auto | <model name>, or blank (=auto)")
    L.append("#   * custom_options         : extra flags as-is, e.g. --nosecondary --localvmem=300")
    L.append("#   * gex_exclude_introns    : any value adds --gex-exclude-introns (skip GEX introns);")
    L.append("#                             blank (default) = include introns.  Feeds 'Include introns?'.")
    L.append("# ============================================================")
    for idx, (base, libs) in enumerate(groups.items(), start=1):
        L.append(f"# Sample {idx}: {base}")
        for suffix, default, comment in _ARC_PER_SAMPLE_PARAMS:
            base_str = f'{suffix}_{idx}="{default}"'
            pad = max(1, 56 - len(base_str))
            L.append(f"    {base_str}{' ' * pad}# {comment}")
    L.append("")

    # ---- 3b) PIPELINE PLUMBING (derived from the parameters above) ----------------------------
    # Emitted BELOW all editable blocks (RUN-WIDE / FASTQ / PER-SAMPLE) so the
    # user-facing "edit here" region stays contiguous at the top of the
    # script.  The reference-genome check fails fast BEFORE any mkdir, so a
    # bad REFPATH leaves no stray directories behind.  A banner marks the
    # user/code boundary unambiguously.
    L.append("")
    L.append("# " + "=" * 73)
    L.append("#   ***  END OF USER-EDITABLE SECTION  ***")
    L.append("#")
    L.append("#   Everything ABOVE is yours to edit: RUN-WIDE PARAMETERS, FASTQ PATHS,")
    L.append("#   PER-SAMPLE PARAMETERS.")
    L.append("#")
    L.append("#   Everything BELOW is pipeline machinery — pure code, DO NOT EDIT.")
    L.append("# " + "=" * 73)
    L.append("")
    L.append("# PIPELINE PLUMBING  (derived from the parameters above)")
    # Reference genome existence check (MUST run before any side effects)
    L.append("# Check that the reference genome directory exists")
    L.append('if [[ ! -d "${REFPATH}" ]]; then')
    L.append('    echo "Reference genome ${REFPATH} does not exist. Aborting." >&2')
    L.append("    exit 1")
    L.append("fi")
    L.append("")
    # Derived working-dir variables (used throughout the sample functions)
    L.append('Server_folder="${scratch_root}/${server_folder_name}"')
    L.append('metrics_folder="${Server_folder}/raw_metrics_folder"')
    L.append('archive_dir="${NGSruns_archive}/${server_folder_name}"')
    L.append('mkdir -p "${metrics_folder}/GEX_only" "${metrics_folder}/ATAC_only"')
    L.append('mkdir -p "${archive_dir}"')
    L.append('mkdir -p "${Server_folder}/tso_report_parts"')
    L.append("")

    # ---- tso_percent() library-quality helper (inlined from common/TSO_percent.sh) ----
    # (parity item 9) The body of `tso_percent()` is inlined once per generated
    # script, so the run needs no runtime dependency on the cellranger_pipelines
    # repo.  Called per sample (from the arc_sample_<i>() functions) on the
    # sample's GEX-leg Read-2 fastq.
    #
    # Emitted AFTER both parameter blocks (run-wide + per-sample) so every
    # editable parameter is listed first and the tso_percent helper follows  —
    # the same top-to-bottom order the GEX generated script uses.
    L.append("# ==============================================================================")
    L.append("# tso_percent helper (from cellranger_pipelines/common/TSO_percent.sh, embedded verbatim)")
    L.append("# Computes Read-2 percentage metrics: intact TSO, intact SMART, polyA/polyT tails.")
    L.append("# ==============================================================================")
    L.append(_TSO_HELPER_BODY)
    L.append("")

    # ---- 4) arc_sample_<i>() functions (item 5 — wrapped in functions) ---------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE FUNCTIONS  (one arc_sample_<i> per biologic sample)")
    L.append("#   Each function:")
    L.append("#     * builds the libraries CSV for this sample")
    L.append("#     * runs cellranger-arc count (flags from the *_<i> variables above)")
    L.append("#     * checks for leftover tmp folders -> <sample>_error_log.txt")
    L.append("#     * harvests web summaries (joint + GEX_only + ATAC_only)")
    L.append("#     * deletes the SC_ATAC_GEX_COUNTER_CS scratch tree")
    L.append("#     * rsync's the sample dir + libraries CSV + logs to the archive")
    L.append("# ============================================================")

    for idx, (base, libs) in enumerate(groups.items(), start=1):
        ordered = sorted(libs, key=lambda t: 0 if t[2] == _GEX else 1)
        fids = [fid for _, fid, _ in ordered if fid]
        desc = ",".join(fids) if fids else base
        n_gex = sum(1 for _, _, m in libs if m == _GEX)
        n_atac = sum(1 for _, _, m in libs if m == _ATAC)
        run_name = run_names[idx - 1]
        csv_name = f"libraries_{run_name}.csv"

        L.append(f"# --- arc_sample_{idx}() — Sample {idx}: {base}  ({n_gex} GEX + {n_atac} ATAC) ---")
        L.append(f"arc_sample_{idx}() {{")

        # Build the libraries CSV (one row per library)
        L.append(f'    csv_path="${{Server_folder}}/{csv_name}"')
        L.append('    echo "fastqs,sample,library_type" > "${csv_path}"')
        for full, fid, modality in ordered:
            var = (f"fastq_{idx}_{fid}" if fid
                   else f"fastq_{idx}_{modality.replace(' ', '_')}")
            L.append(f'    echo "${{{var}}},{full},{modality}" >> "${{csv_path}}"')
        L.append("")

        # cellranger-arc count — flag_lines list = what runs AND appears in the .sh
        flag_lines = [
            'cellranger-arc count',
            f'    --id={run_name}',
            f'    --description={desc}',
            '    --reference="${REFPATH}"',
            # absolute --libraries path (the CSV is written into Server_folder
            # just above this command) — no cd needed
            f'    --libraries="${{Server_folder}}/{csv_name}"',
            # --output-dir pins the run folder to the exact layout a cd-based
            # invocation used to produce (run root = Server_folder/<id>),
            # without ever moving the main shell's CWD
            f'    --output-dir="${{Server_folder}}/{run_name}"',
            # CPU/memory caps are run-wide (localcores / localmem), not per-sample
            '    ${localcores:+--localcores=$localcores}',
            '    ${localmem:+--localmem=$localmem}',
        ]
        for suffix, _default, _comment in _ARC_PER_SAMPLE_PARAMS:
            v = f"{suffix}_{idx}"
            if suffix == "custom_options":
                flag_lines.append(f"    ${{{v}}}")
            elif suffix == "gex_exclude_introns":
                flag_lines.append(f"    ${{{v}:+--gex-exclude-introns}}")
            else:
                flag = "--" + suffix.replace("_", "-")
                flag_lines.append(f"    ${{{v}:+{flag}=${v}}}")
        flag_lines.append(f'    --tenx-cloud-token-path={TENX_TOKEN_PATH}')

        # Log the exact command, WITH RESOLVED VALUES, to <base>_script_log.txt
        # (one readable line).  A single bash array (_cmd) is built from the
        # SAME flag_lines list the execution line uses, so the log can never
        # drift from what actually runs: `echo "${_cmd[@]}"` expands
        # ${REFPATH}, ${localcores:+--localcores=$localcores},
        # ${localmem:+--localmem=$localmem} and the per-sample flag values
        # exactly as the execution sees them.
        # (The old version printf'd a Python-repr()'d string — a
        # single-quoted echo — which wrote the unexpanded ${...} template
        # into the log: parameter NAMES without VALUES.)
        L.append('    _cmd=( ' + " ".join(fl.strip() for fl in flag_lines) + ' )')
        L.append(f'    _cmd_log="${{Server_folder}}/{run_name}_script_log.txt"')
        L.append('    echo "${_cmd[@]}" > "${_cmd_log}"')
        L.append("")
        # Execute the command (same flag_lines list, line-continuation formatted)
        for i, fl in enumerate(flag_lines):
            cont = " \\" if i < len(flag_lines) - 1 else ""
            L.append(f"    {fl}{cont}")
        L.append("")

        # Completion check: leftover tmp folders indicate cellranger-arc count did not finish
        L.append(f'    _err_log="${{Server_folder}}/{run_name}_error_log.txt"')
        L.append(f'    find "${{Server_folder}}/{run_name}" -mindepth 2 -maxdepth 2 '
                 f'-type d -name "tmp" -print 2>/dev/null > "${{_err_log}}"')
        L.append('    if [ -s "${_err_log}" ]; then')
        L.append(f'        echo "Error: cellranger-arc count for {base} did not complete '
                 f'(leftover tmp folders detected).\\nSee ${{Server_folder}}/{run_name}_error_log.txt" >&2')
        L.append('    else')
        L.append('        rm -f "${_err_log}"')
        L.append('    fi')
        L.append("")

        # Harvest web summaries (joint + GEX_only + ATAC_only)
        L.append("    # Harvest web summaries (joint + GEX_only + ATAC_only)")
        L.append(f'    if [ -f "${{Server_folder}}/{run_name}/outs/web_summary.html" ]; then')
        L.append(f'        cp "${{Server_folder}}/{run_name}/outs/web_summary.html" '
                 f'"${{metrics_folder}}/{run_name}_web_summary.html"')
        L.append("    fi")
        L.append(f'    if [ -f "${{Server_folder}}/{run_name}/outs/summary.csv" ]; then')
        L.append(f'        cp "${{Server_folder}}/{run_name}/outs/summary.csv" '
                 f'"${{metrics_folder}}/{run_name}_summary.csv"')
        L.append(f'    elif [ -f "${{Server_folder}}/{run_name}/outs/metrics_summary.csv" ]; then')
        L.append(f'        cp "${{Server_folder}}/{run_name}/outs/metrics_summary.csv" '
                 f'"${{metrics_folder}}/{run_name}_summary.csv"')
        L.append('    else')
        L.append(f'        echo "  [warn] no summary.csv/metrics_summary.csv for {base}" >&2')
        L.append('    fi')
        L.append(f'    if [ -f "${{Server_folder}}/{run_name}/outs/cell_types/web_summary_cell_types.html" ]; then')
        L.append(f'        cp "${{Server_folder}}/{run_name}/outs/cell_types/web_summary_cell_types.html" '
                 f'"${{metrics_folder}}/{run_name}_web_summary_cell_types.html"')
        L.append("    fi")
        # GEX-only web summary (from the CS tree — harvest before CS-tree deletion)
        L.append(f'    _gex_web=$(find "${{Server_folder}}/{run_name}/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/GEX_SUMMARIZE_REPORTS" \\')
        L.append('        -path "*/files/web_summary.html" -type f 2>/dev/null | head -1 || true)')
        L.append('    if [ -n "${_gex_web}" ]; then')
        L.append(f'        cp "${{_gex_web}}" "${{metrics_folder}}/GEX_only/{run_name}_web_summary.html"')
        L.append('    else')
        L.append(f'        echo "  [warn] no GEX-only web summary for {base}" >&2')
        L.append('    fi')
        # ATAC-only web summary (from the CS tree)
        L.append(f'    _atac_web=$(find "${{Server_folder}}/{run_name}/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/_SC_ATAC_REPORTER/CREATE_WEBSUMMARY" \\')
        L.append('        -path "*/files/web_summary.html" -type f 2>/dev/null | head -1 || true)')
        L.append('    if [ -n "${_atac_web}" ]; then')
        L.append(f'        cp "${{_atac_web}}" "${{metrics_folder}}/ATAC_only/{run_name}_web_summary.html"')
        L.append('    else')
        L.append(f'        echo "  [warn] no ATAC-only web summary for {base}" >&2')
        L.append('    fi')
        L.append("")

        # Delete the SC_ATAC_GEX_COUNTER_CS scratch tree (summaries already harvested)
        L.append("    rm -r \"${Server_folder}/" + run_name + "/SC_ATAC_GEX_COUNTER_CS\" 2>/dev/null || true")
        L.append("")

        # rsync sample dir + auxiliary files to the archive
        L.append("    # Archive this sample + auxiliary files")
        L.append('    if command -v rsync >/dev/null 2>&1; then')
        L.append(f'        if [ -d "${{Server_folder}}/{run_name}" ]; then')
        L.append(f'            rsync -avh "${{Server_folder}}/{run_name}" "${{archive_dir}}/"')
        L.append('        fi')
        L.append(f'        if [ -f "${{Server_folder}}/{csv_name}" ]; then')
        L.append(f'            rsync -avh "${{Server_folder}}/{csv_name}" "${{archive_dir}}/"')
        L.append('        fi')
        L.append(f'        if [ -f "${{Server_folder}}/{run_name}_script_log.txt" ]; then')
        L.append(f'            rsync -avh "${{Server_folder}}/{run_name}_script_log.txt" "${{archive_dir}}/"')
        L.append('        fi')
        L.append(f'        if [ -f "${{Server_folder}}/{run_name}_error_log.txt" ]; then')
        L.append(f'            rsync -avh "${{Server_folder}}/{run_name}_error_log.txt" "${{archive_dir}}/" 2>/dev/null || true')
        L.append('        fi')
        L.append('    fi')

        # ---- TSO% library-quality metric (parity item 9) — GEX leg only ----
        # intact-TSO / intact-SMART / polyA+polyT fractions, measured on a
        # downsampled Read-2 subset of EACH GEX-leg fastq — one section per
        # avid run, not just the first leg.  The ATAC leg is chromatin and
        # carries no meaningful TSO% reading, so it is never measured.
        # Each section is tee'd to a part under tso_report_parts/; the
        # single centralized tso_report.txt is assembled AFTER DISPATCH
        # (samples may run concurrently, so each part has exactly one
        # writer).  See the "TSO% REPORT (SINGLE FILE)" section below.
        gex_fids = [fid for _, fid, m in libs if m == _GEX and fid]
        if gex_fids:
            for fid in gex_fids:
                _tso_var = f"fastq_{idx}_{fid}"
                _tso_tags = f"{idx:02d}_{base}__{fid}"
                L.append("")
                L.append(f"    # TSO% ({fid}) — intact TSO/SMART + poly tails on Read-2")
                L.append(f'    _tso_fq="${{{_tso_var}}}"')
                L.append(f'    _tso_tags="{_tso_tags}"')
                L.append('    if [ ! -f "${_tso_fq}" ] && [ -d "${_tso_fq}" ]; then')
                L.append('        for _c in "${_tso_fq}"/*_R2_*.fastq.gz "${_tso_fq}"/*_R2.fastq.gz; do')
                L.append('            [ -f "${_c}" ] && { _tso_fq="${_c}"; break; }')
                L.append("        done")
                L.append("    fi")
                L.append('    if [ -f "${_tso_fq}" ]; then')
                L.append(f'        _tso_part="${{Server_folder}}/tso_report_parts/{_tso_tags}.txt"')
                L.append('        tso_percent "${_tso_fq}" "${_tso_tags}" | tee "${_tso_part}"')
                L.append("    else")
                L.append(f'        echo "tso: no Read-2 fastq near {_tso_var}; TSO% step skipped ({fid})." >&2')
                L.append("    fi")
        else:
            L.append("")
            L.append(f'    echo "[{idx:02d}_{base}] no GEX leg in this sample; TSO% step skipped (ATAC-only)" >&2')

        L.append("}")
        L.append("")

    # ---- 5) DISPATCH (item 6) -------------------------------------------------------------------
    N = len(groups)
    L.append("# ============================================================")
    L.append("# DISPATCH  (run arc_sample_1 .. arc_sample_" + str(N) + ")")
    L.append("#   max_parallel=1 : serial, one sample after another")
    L.append("#   max_parallel=N : rolling pool of N concurrent cellranger-arc counts —")
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
    L.append("        arc_sample_${_s}")
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
    L.append("                arc_sample_${_next} &")
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
    L.append("        if [ \"${_launched}\" -eq 0 ] && [ \"${_next}\" -le \"${NSAMPLES}\" ]; then")
    L.append("            sleep 1")
    L.append("        fi")
    L.append("    done")
    L.append("fi")
    L.append("")

    # ---- 6) per-sample status sweep (item 7) -----------------------------------------------------
    L.append("# ============================================================")
    L.append("# PER-SAMPLE STATUS REPORT")
    L.append("#   Each sample's own completion markers determine the outcome:")
    L.append("#     *_error_log.txt non-empty  -> cellranger-arc count left tmp dirs on failure")
    L.append("#     no outs/summary.csv        -> cellranger-arc count never completed")
    L.append("# ============================================================")
    L.append("_overall_rc=0")
    _dirs_bash = " ".join('"' + d + '"' for d in run_names)
    L.append("for _d in " + _dirs_bash + "; do")
    L.append('    if [ -s "${Server_folder}/${_d}_error_log.txt" ]; then')
    L.append('        echo "[FAIL] ${_d}  (see ${Server_folder}/${_d}_error_log.txt)" >&2')
    L.append("        _overall_rc=1")
    L.append('    elif [ ! -f "${Server_folder}/${_d}/outs/summary.csv" ] && [ ! -f "${Server_folder}/${_d}/outs/metrics_summary.csv" ]; then')
    L.append('        echo "[FAIL] ${_d}  (no summary output — cellranger-arc count did not complete)" >&2')
    L.append("        _overall_rc=1")
    L.append("    else")
    L.append('        echo "[ OK ] ${_d}"')
    L.append("    fi")
    L.append("done")
    L.append("")

    # ---- 6b) TSO% REPORT (single file) — parity item 9 -------------------------------------------
    L.append("# ============================================================")
    L.append("# TSO% REPORT  (SINGLE FILE)")
    L.append("#   Assembles one tso_report.txt from the (sample, avid) parts under")
    L.append("#   tso_report_parts/ (per GEX leg: intact-TSO / intact-SMART /")
    L.append("#   polyA+polyT block), in sample order then avid order.  ATAC-only")
    L.append("#   samples write no part (no Read-2 fastq to measure) and are")
    L.append("#   simply absent.")
    L.append("# ============================================================")
    L.append('_tso_out="${Server_folder}/tso_report.txt"')
    L.append("{")
    L.append('    echo "# ============================================================"')
    L.append('    echo "# TSO% REPORT — one section per (sample, avid) — ' + str(len(groups)) + ' sample(s) — AVT run: ${avt_run}"')
    L.append('    echo "# Generated: $(date -u)  —  generator ' + SCRIPT_VERSION + '"')
    L.append('    echo "# (per GEX leg: intact TSO % / intact SMART % / polyA tail % / polyT tail %)"')
    L.append('    echo "# ============================================================"')
    L.append('    echo ""')
    for idx, (base, _libs) in enumerate(groups.items(), start=1):
        _gex_fids = [fid for _, fid, m in _libs if m == _GEX and fid]
        if not _gex_fids:
            continue      # ATAC-only sample — no part to list
        for _fid in _gex_fids:
            _tag = "%02d_%s__%s" % (idx, base, _fid)
            L.append('    if [ -f "${Server_folder}/tso_report_parts/%s.txt" ]; then' % _tag)
            L.append('        echo "### Sample %d: %s — %s (TSO%%) ###"' % (idx, base, _fid))
            L.append('        cat "${Server_folder}/tso_report_parts/%s.txt"' % _tag)
            L.append('        echo ""')
            L.append("    else")
            L.append('        echo "### Sample %d: %s — %s (TSO%%) ###"' % (idx, base, _fid))
            L.append('        echo "[no TSO%% report — no Read-2 fastq found or tso_percent skipped]"')
            L.append('        echo ""')
            L.append("    fi")
    L.append('} > "${_tso_out}"')
    L.append('rm -rf "${Server_folder}/tso_report_parts"')
    L.append('echo "TSO%% report (single file): ${_tso_out}"')
    L.append('if command -v rsync >/dev/null 2>&1; then rsync -avh "${_tso_out}" "${archive_dir}/" 2>/dev/null || true; fi')
    L.append("")

    # ---- 6c) MERGED SUMMARY METRICS (one summary_metrics.csv) — parity item 10 ---------------------
    L.append("# ============================================================")
    L.append("# MERGED SUMMARY METRICS  (one summary_metrics.csv, Samplename column first)")
    L.append("#   Gathers each sample's outs/summary.csv into a single CSV so all")
    L.append("#   samples' key metrics can be read in one table.")
    L.append("# ============================================================")
    L.append('SM_TMP="${Server_folder}/summary_metrics_tmp.csv"')
    L.append('SM_OUT="${Server_folder}/summary_metrics.csv"')
    L.append(': > "${SM_TMP}"')
    L.append('_sm_hdr_written=0')
    for _i, _base in enumerate(groups, 1):
        _dn = run_names[_i - 1]
        L.append('    _sm_src="${Server_folder}/%s/outs/summary.csv"' % _dn)
        L.append('    [ -f "${_sm_src}" ] || _sm_src="${Server_folder}/%s/outs/metrics_summary.csv"' % _dn)
        L.append('    if [ -f "${_sm_src}" ]; then')
        L.append('        _sm_hdr=$(head -n1 "${_sm_src}" | tr -d "\\t")')
        L.append('        _sm_dat=$(tail -n +2 "${_sm_src}" | head -n1 | tr -d "\\t")')
        L.append('        if [ "${_sm_hdr_written}" -eq 0 ] && [ -n "${_sm_hdr}" ]; then')
        L.append('            echo "Samplename,${_sm_hdr}" > "${SM_TMP}"')
        L.append('            _sm_hdr_written=1')
        L.append('        fi')
        L.append('        if [ -n "${_sm_dat}" ]; then echo "%s,${_sm_dat}" >> "${SM_TMP}"; fi' % _dn)
        L.append("    fi")
    L.append('if [ -s "${SM_TMP}" ]; then')
    L.append('    cp "${SM_TMP}" "${SM_OUT}"')
    L.append('    rm -f "${SM_TMP}"')
    L.append('    echo "Merged summary metrics: ${SM_OUT}"')
    L.append('    if command -v rsync >/dev/null 2>&1; then rsync -avh "${SM_OUT}" "${archive_dir}/" 2>/dev/null || true; fi')
    L.append("else")
    L.append('    rm -f "${SM_TMP}"')
    L.append('    echo "  [warn] no per-sample summary.csv found; summary_metrics.csv not created" >&2')
    L.append("fi")
    L.append("")

    # ---- 6b) per-run stats workbook (GEX + ATAC) -> Server_folder ----
    #   Two-sheet xlsx (9 cols, one row per sample per modality), built at run
    #   time from the run data.  Written into ${Server_folder} so the
    #   whole-folder md5 below (item 8a) covers it; its own rsync ships it to
    #   the archive.  Non-fatal: no python3 -> warn+skip, no openpyxl -> CSVs.
    L.extend(build_stats_block(groups, fid_avt, avt_run, reference_genome))

    # ---- 7) rsync raw_metrics_folder -------------------------------------------------------------
    L.append("# ============================================================")
    L.append("# POST-PROCESSING  (rsync raw_metrics_folder to the archive)")
    L.append("# ============================================================")
    L.append('if [ -d "${metrics_folder}" ] && command -v rsync >/dev/null 2>&1; then')
    L.append('    rsync -avh "${metrics_folder}/" "${archive_dir}/raw_metrics_folder/"')
    L.append("fi")
    L.append("")

    # ---- 8a) whole-folder md5sums.txt (item 8) ---------------------------------------------------
    L.append("# ============================================================")
    L.append("# MD5 CHECKSUMS  (WHOLE FOLDER)")
    L.append("#   One md5sums.txt covering every file under ${Server_folder},")
    L.append("#   relative ./ paths.")
    L.append("# ============================================================")
    L.append('_md5_out="${Server_folder}/md5sums.txt"')
    L.append('( cd "${Server_folder}" && find . -type f ! -name "md5sums.txt" | sort | xargs md5sum > "md5sums.txt" )')
    L.append('echo "MD5 checksums (whole folder): ${_md5_out}"')
    L.append('rsync -avh "${_md5_out}" "${archive_dir}/" 2>/dev/null || true')
    L.append("")

    # ---- 8b) copy the .sh used to generate this project ------------------------------------------
    L.append("# ============================================================")
    L.append("# KEEP A COPY OF THIS SCRIPT")
    L.append("#   Copy the .sh used to generate this project into the project folders")
    L.append("#   (local working dir + central archive) so the exact pipeline is")
    L.append("#   preserved with the results.")
    L.append("# ============================================================")
    # Self-resolve THIS script to an absolute path right where it is used —
    # this far down the script nothing has `cd`'d in the main shell (the
    # arc_sample_<i> functions call cellranger-arc count with --output-dir
    # and the md5sum step runs in a subshell), so a relative BASH_SOURCE[0]
    # still means "launched from here".  Kept here (not at the top) because
    # it is plumbing, not a user-editable parameter.  NOTE: if a bare `cd`
    # is ever added to the main flow ABOVE this point, the relative name
    # breaks — then move this capture above that `cd` instead.
    L.append('# Path of this script as an absolute path, resolved against the launch directory')
    L.append('_SELF_PATH="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)/$(basename "${BASH_SOURCE[0]:-$0}")"')
    L.append('if [ ! -f "${_SELF_PATH}" ]; then _SELF_PATH="${BASH_SOURCE[0]:-$0}"; fi')
    L.append('if [ -n "${_SELF_PATH:-}" ] && [ -f "${_SELF_PATH}" ]; then')
    L.append('    cp "${_SELF_PATH}" "${Server_folder}/$(basename "${_SELF_PATH}")" 2>/dev/null || true')
    L.append('    cp "${_SELF_PATH}" "${archive_dir}/$(basename "${_SELF_PATH}")" 2>/dev/null || true')
    L.append("fi")
    L.append("")

    # ---- 10) keepdata cleanup ---------------------------------------------------------------------
    L.append("# ============================================================")
    L.append("# LOCAL WORKDIR CLEANUP")
    L.append('#   After deploy to the archive, the local working folder is deleted')
    L.append('#   UNLESS keepdata="yes".  (MD5 checksums and the .sh copy are done first.)')
    L.append("# ============================================================")
    L.append('if [[ "${keepdata}" != "yes" ]]; then')
    L.append('    rm -r "${Server_folder}" 2>/dev/null || true')
    L.append("fi")
    L.append("")

    # ---- 11) park processed datasheets in the team "TO delete" folder -----------------------------
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
    if not _srcs_abs:
        # bash rejects an empty `then ... fi` body — keep one (no-op) command
        L.append('    :   # (no source datasheets to park)')
    L.append("fi")
    L.append("")

    # ---- 12) exit ---------------------------------------------------------------------------------
    L.append('exit "${_overall_rc}"')

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
        "-v", "--version",
        action="version",
        version=f"generate_cellranger_template_arc.py {SCRIPT_VERSION}",
        help="Print the generator version and exit.",
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
    parser.add_argument(
        "--no-number-prefix",
        action="store_true",
        default=False,
        help=(
            "Do not number the per-sample run folders (ON by default): "
            "cellranger-arc count --id becomes '<sample>' instead of "
            "'01_<sample>', and the folder keeps the bare sample name. "
            "Numbering keeps the folders in sample order even when the "
            "names are not alphabetical."
        ),
    )
    args = parser.parse_args()

    set_non_interactive(args.non_interactive)
    default_genome = {
        "human": HUMAN_GENOME_PATH,
        "mouse": MOUSE_GENOME_PATH,
    }[args.genome]

    all_rows, detected_avts, any_user = [], [], ""
    source_files = []  # accepted Excel file(s) — parked in TO_DELETE by the generated .sh
    for path in args.xlsx_files:
        if not os.path.isfile(path) or not path.lower().endswith((".xlsx", ".xlsm")):
            print(f"Warning: skipping '{path}' (not found or not .xlsx).")
            continue
        source_files.append(path)
        print(f"Reading '{os.path.basename(path)}' …", flush=True)
        rows, avt = parse_xlsx(path)
        if rows:
            all_rows.extend(rows)
        if avt:
            detected_avts.append(avt)
        if not any_user:
            print("Reading the requester (User) field …", flush=True)
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

    # ── Parity item 3: sample selection (keep/drop + add sheet-less) ───
    # Same shared helper + prompt as the GEX generator.  ``ask_sample_selection``
    # wants ``{sample: [(fastq_id, avt_run)]}`` (the defaultdict shape the
    # generators build); arc keys it by the BIologic base name (GEX/_ATAC legs
    # collapse to one entry).  Returns (kept_names, added) where added maps a
    # sample name to the FULL avid-folder path(s) the operator typed — those
    # are self-referential and used verbatim (no AVT-tree walk).  Under -n the
    # helper returns every sample kept and adds none (identical to old behaviour).
    sample_fastqs = defaultdict(list)
    for _row in all_rows:
        # (sample, fid, [per-row avt]); the per-row avt (parity item 11) is
        # what the selection helper sees per FASTQ — any_avt is only the
        # fallback for rows that name no run.
        _sample, _fid = _row[0], _row[1]
        _ravt = _row[2] if len(_row) > 2 else None
        sample_fastqs[base_sample_name(_sample)].append((
            (_fid or ""), (_ravt or any_avt)))
    kept_names, added_raw = ask_sample_selection(
        sample_fastqs, default_root=NGSruns_ARCHIVE_DIR)
    if dropped := [sn for sn in sample_fastqs if sn not in kept_names]:
        print(f"  Dropped from script: {', '.join(dropped)}")
    all_rows = [r for r in all_rows
                if base_sample_name(r[0]) in kept_names]
    # Added (sheet-less) samples: carry a full avid-folder path per leg.  The
    # row's fastq id is the folder's basename (e.g. avid13144); the generator
    # uses the full path verbatim.  A sample added via two folders (GEX + ATAC)
    # contributes one row per folder.
    added_paths = {}   # {fastq id (avid basename): full avid-folder path}
    for name, folders in added_raw.items():
        for p in folders:
            added_paths[os.path.basename(p)] = p
            # Per-row avt = None: a folder added by full path keeps the run
            # implied by the folder the operator typed (parity item 11).
            all_rows.append((name, os.path.basename(p), None))
    if not all_rows:
        sys.exit("Error: no samples selected — nothing to generate.")

    # ── Parity item 11: per-row AVT — resolve each FASTQ under its own run tree ──
    # Each sample's FASTQ is resolved under the AVT run named on its own sheet
    # row (the run may sit in the fastq-ID cell or an adjacent cell), so a
    # worksheet mixing two runs resolves correctly.  A row that names no run
    # inherits the sheet-wide / file-name-detected run; a sample added by full
    # path keeps the run implied by the folder that was typed.  When every
    # FASTQ names one run, the outcome is identical to the pre-item-11
    # generator (single-run workspace).
    fid_avt = {}
    base_avt = {}
    for _row in all_rows:
        _sample, _fid, _ravt = _row[0], _row[1], (_row[2] if len(_row) > 2 else None)
        _fid = (_fid or "").strip()
        if _fid:
            fid_avt.setdefault(_fid, _ravt or any_avt)
        if _ravt:
            base_avt.setdefault(base_sample_name(_sample), _ravt)
    _distinct_avts = sorted(set(fid_avt.values()) | set(base_avt.values()))
    if len(_distinct_avts) > 1:
        print(f"  Note: this workspace spans AVT runs {', '.join(_distinct_avts)}; "
              f"each FASTQ resolves under the run named on its own row.")

    groups = group_libraries(all_rows)
    print()
    print_groups_table(groups)

    reference_genome = ask_reference_genome(
        HUMAN_GENOME_PATH, MOUSE_GENOME_PATH, default_path=default_genome,
        human_label="Human", mouse_label="Mouse")
    user = any_user
    run_dir = NGS_BASE_DIR_TEMPLATE.format(avt_run=any_avt)

    # ── Parity item 4: run-folder disambiguation (AVT0551 case) ─────────
    # Before: resolve_fastq_leaf() blind-walked the run tree and took the
    # FIRST hit — the AVT0551 case, where the same sample recurs under
    # 1_wrong_index_… AND 2_corrected_indexes_… and silently resolved to the
    # wrong library.  Now we pin the SAMPLE's directory (keyed by base) via
    # list_sample_dirs/choose_sample_dir before resolving; choose_sample_dir
    # prompts only when >1 candidate exists (the 99% single-candidate case
    # resolves silently) and also no-ops under -n (first candidate + a warning
    # naming the alternatives).  Samples added by full path have no candidate
    # set — they are used verbatim.
    chosen_dirs = {}
    for base in kept_names:
        cands = list_sample_dirs(base_avt.get(base, any_avt), base)
        if len(cands) > 1:
            chosen = choose_sample_dir(base, cands)
            if chosen:
                chosen_dirs[base] = chosen
    if chosen_dirs:
        print(f"  Run-folder disambiguation pinned {len(chosen_dirs)} sample(s).")

    # Emission-time resolver: build_script calls resolve(base, fid).  Each
    # FASTQ walks the run tree of the AVT run named on its own row (parity
    # item 11); `any_avt` is only the fallback for rows that name no run.
    def resolve(base, fid):
        _fid = (fid or "").strip()
        _f_avt = fid_avt.get(_fid, any_avt) if _fid else any_avt
        if _fid and _fid in added_paths:          # added-by-path sample → verbatim
            return added_paths[_fid]
        pref = chosen_dirs.get(base)
        if pref:                                  # preferred sample dir first
            hit = find_fastq_path_with_pref(_f_avt, _fid,
                                            preferred_sample_dir=pref)
            if hit:
                return hit
        _run = NGS_BASE_DIR_TEMPLATE.format(avt_run=_f_avt)
        return resolve_fastq_leaf(_run, _fid)     # full-tree walk (parity item 11)

    script = build_script(reference_genome, groups, any_avt, user, run_dir,
                          resolve=resolve, source_files=source_files,
                          fid_avt=fid_avt,
                          number_prefix=not args.no_number_prefix)
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
