"""Shared GECF Excel intake atoms.

Every 10XG cellranger template generator (ATAC, ARC, and — eventually — Multi
and Visium) shares the same low-level workbook plumbing:

* merged-cell resolution
* cell value reading with graceful fallback
* header normalisation
* the ``User:`` requester field
* the shared interactive prompts (confirm AVT run name, choose
  reference genome) and the readline path-completion setup

The per-assay generators (``pipelines/{atac,arc}/generate_cellranger_template_*.py``)
import *these* helpers from here and layer their own sheet-selection, row-gate,
AVT-extraction, and data-model logic on top.

This module is the single source of truth for the atoms.  If any of the helpers
have a bug or need a format change, the fix lives in one place and propagates
to all generators on their next import — no cross-file syncing required.

Import convention (from a generator file at ``pipelines/<assay>/...py``):

    import os, sys
    _HERE      = os.path.dirname(os.path.abspath(__file__))
    _PIPELINES = os.path.normpath(os.path.join(_HERE, ".."))
    if _PIPELINES not in sys.path:
        sys.path.insert(0, _PIPELINES)

    from common.intake import (
        NGS_BASE_DIR_TEMPLATE,     # template literal w/ {avt_run} placeholder
        AVT_PATTERN,               # compiled re for 'AVT0519'-style names
        norm,                      # lowercased/stripped cell value
        build_merged_cell_lookup,  # (row,col)->value dict for merged ranges
        get_cell_value,            # merged-aware cell read
        read_user_field,           # 'User:' label -> first word
        set_non_interactive,       # flag prompt helpers
        is_non_interactive,        # read prompt-flag state
        ask_avt_run,               # confirm / ask the AVT run name
        ask_reference_genome,      # choose a reference genome path
        setup_readline_completer,  # tab-completion for input() prompts
    )

Non-interactive mode (CI / headless / ``--non-interactive`` CLI flag):
call ``set_non_interactive(True)`` once in the generator's ``main()`` before
any prompt helper is invoked.
"""

from __future__ import annotations

import glob
import os
import re

# ── shared path constants ─────────────────────────────────────────────────

#: Base directory under which per-AVT FASTQ trees live on the GECF server.
#: The generator interpolates the AVT run name in at use-time.
NGS_BASE_DIR_TEMPLATE = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/{avt_run}"

#: Anchored match for a bare AVT run name (e.g. ``AVT0519``).
AVT_PATTERN = re.compile(r"AVT0?\d{3,5}")

#: A real GECF library id is ``avid`` followed by 4-6 digits (e.g. ``avid12499``).
#: Label / legend / placeholder rows in the same columns carry strings like
#: ``"0"``, ``"GECF libr. ID (avidxxxx)"`` or ``"Loaded molarity (pM)"`` — this
#: pattern distinguishes genuine sample rows from those.  Both the ARC and ATAC
#: generators use this as their row gate, so the single source of truth lives
#: here rather than in each generator.
FASTQ_ID_PATTERN = re.compile(r"^avid\d{4,6}$", re.IGNORECASE)

#: Default ``--localcores`` / ``--localmem`` for cellranger invocations emitted
#: by the generators.  Shared so the numbers stay in lockstep across assays;
#: per-assay overrides remain possible via the generator's own constant.
LOCAL_CORES = "44"
LOCAL_MEM = "150"


# ── non-interactive mode ──────────────────────────────────────────────────

_NON_INTERACTIVE = False


def set_non_interactive(flag: bool) -> None:
    """Set the module-level non-interactive flag (affects prompt helpers)."""
    global _NON_INTERACTIVE
    _NON_INTERACTIVE = bool(flag)


def is_non_interactive() -> bool:
    """Return the current non-interactive flag."""
    return _NON_INTERACTIVE


# ── cell-level helpers ────────────────────────────────────────────────────


def norm(value) -> str:
    """Return a lowercased, whitespace-normalised string form of *value*.

    ``None`` → ``""``.  Every run of whitespace and any ``_ . -`` separator
    is collapsed to a single space, so header cells such as
    ``Sample  Name`` / ``Sample_Name`` / ``sample-name`` all match the same
    way.  Used for case-insensitive header / alias matching in both
    generators — this is the canonical form (the ARC variant; the ATAC
    variant, which kept ``_`` intact, is a strict subset for the substring
    matches those generators perform).
    """
    if value is None:
        return ""
    return re.sub(r"[\s._-]+", " ", str(value)).strip().lower()


def build_merged_cell_lookup(ws) -> dict:
    """Map ``(row, col) → top_left_value`` for every cell inside any merged range.

    Header rows that merge a label across two rows (or a column) still resolve
    correctly through this lookup, which is how ATAC and ARC workbooks both
    survive openpyxl's otherwise unreliable merged-cell detection.
    """
    lookup: dict = {}
    for merged_range in ws.merged_cells.ranges:
        min_row, min_col = merged_range.min_row, merged_range.min_col
        max_row, max_col = merged_range.max_row, merged_range.max_col
        top_left_value = ws.cell(row=min_row, column=min_col).value
        for r in range(min_row, max_row + 1):
            for c in range(min_col, max_col + 1):
                lookup[(r, c)] = top_left_value
    return lookup


def get_cell_value(ws, merged_lookup, row_idx: int, col_idx: int):
    """Return the value at ``(row_idx, col_idx)`` honouring merged ranges.

    Falls back to a plain ``ws.cell(...)`` read when the position is not part
    of a merge, and swallows ``Exception`` (malformed cells, out-of-bounds)
    returning ``None`` rather than raising mid-parse.
    """
    try:
        if (row_idx, col_idx) in merged_lookup:
            val = merged_lookup[(row_idx, col_idx)]
            if val is not None:
                return val
        return ws.cell(row=row_idx, column=col_idx).value
    except Exception:
        return None


# ── 'User:' requester field ───────────────────────────────────────────────


def read_user_field(path: str) -> str:
    """Locate a ``User:`` label (any row, any sheet) in the workbook at
    *path* and return the first whitespace-delimited word of the first
    non-empty cell to its right.

    Returns ``""`` when the label is absent, the adjacent cell is blank,
    or the workbook cannot be opened.

    GECF submission sheets carry the requester at approximately row 2
    (label col 1, value col 2) but the exact position has moved between
    templates, so a full sheet scan is used rather than a hard-coded
    row/column.  Cell content is normalised by stripping every
    character outside ``[a-z]`` so header variants such as ``User:``  /
    ``User``  /  ``USER:`` are all matched.
    """
    try:
        import openpyxl
        wb = openpyxl.load_workbook(path, data_only=True)
    except Exception:
        return ""
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for cell in row:
                if cell.value is None:
                    continue
                label = re.sub(r"[^a-z]", "", str(cell.value).lower())
                if label != "user":
                    continue
                for c in range(cell.column + 1, ws.max_column + 1):
                    v = ws.cell(row=cell.row, column=c).value
                    if v is not None and str(v).strip():
                        return str(v).strip().split()[0]
                return ""
    return ""


# ── shared interactive prompts ─────────────────────────────────────────


def setup_readline_completer() -> None:
    """Enable tab-completion of file paths in interactive ``input()`` prompts."""
    import readline

    def _completer(text, state):
        matches = glob.glob(os.path.expanduser(text) + "*")
        matches = [m + "/" if os.path.isdir(m) else m for m in matches]
        return matches[state] if state < len(matches) else None

    readline.set_completer_delims(" \t\n;")
    readline.set_completer(_completer)
    readline.parse_and_bind("tab: complete")


def ask_avt_run(detected_avt) -> str:
    """Confirm the detected AVT run name, or ask for one.

    In non-interactive mode returns ``detected_avt`` (or ``"AVT0000"``)
    upper-cased, without prompting.
    """
    if is_non_interactive():
        return (detected_avt or "AVT0000").upper()
    if detected_avt:
        while True:
            print(f"Detected run name: {detected_avt}")
            answer = input("Use this run name? [Y/n]: ").strip()
            if not answer or answer[0].lower() == "y":
                return detected_avt.upper()
            if answer[0].lower() == "n":
                break
            print("Answer can only be Y or n. Please try again.")
    while True:
        value = input("Enter AVT run name (e.g. AVT0519): ").strip()
        if value:
            return value.upper()
        print("Run name cannot be empty. Please try again.\n")


def ask_reference_genome(human_path, mouse_path, default_path=None):
    """Prompt for a reference genome path.

    ``default_path`` is returned immediately in non-interactive mode
    (falls back to the human path); interactive mode offers the
    human/mouse choices or accepts a custom path.
    """
    if is_non_interactive():
        return default_path or human_path
    choices = [("Human ARC ATAC", human_path),
               ("Mouse ARC ATAC", mouse_path)]
    while True:
        print("\nAvailable reference genomes:")
        for i, (label, p) in enumerate(choices, start=1):
            print(f"  {i}. {label} — {p}")
        print("  Custom: type a full path")
        answer = input("Choice [1]: ").strip() or "1"
        if answer == "1":
            return human_path
        if answer == "2":
            return mouse_path
        custom = os.path.expanduser(answer)
        if os.path.isdir(custom):
            return custom
        print(f"Path '{custom}' does not exist. Try again.\n")
