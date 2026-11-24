"""Shared GECF Excel intake atoms.

Every 10XG cellranger template generator (ATAC, ARC, and — eventually — Multi
and Visium) shares the same low-level workbook plumbing:

* merged-cell resolution
* cell value reading with graceful fallback
* header normalisation
* the ``User:`` requester field
* the shared interactive prompts (confirm the run name — AVT/MSH/MSQ/NVQ…,
  choose reference genome, select / add samples) and the readline
  path-completion setup

The per-assay generators (``cellranger_pipelines/{atac,arc}/generate_cellranger_template_*.py``)
import *these* helpers from here and layer their own sheet-selection, row-gate,
run-name-extraction (AVT/MSH/MSQ/NVQ… — any code in KNOWN_RUN_CODES), and
data-model logic on top.

This module is the single source of truth for the atoms.  If any of the helpers
have a bug or need a format change, the fix lives in one place and propagates
to all generators on their next import — no cross-file syncing required.

Import convention (from a generator file at ``cellranger_pipelines/<assay>/...py``):

    import os, sys
    _HERE      = os.path.dirname(os.path.abspath(__file__))
    _CELLRANGER_PIPELINES = os.path.normpath(os.path.join(_HERE, ".."))
    if _CELLRANGER_PIPELINES not in sys.path:
        sys.path.insert(0, _CELLRANGER_PIPELINES)

    from common.intake import (
        run_base_dir,            # run name → base folder (all sequencers)
        AVT_PATTERN,             # compiled re for run names (AVT/MSH/MSQ/NVQ + digits)
        NGS_BASE_DIR_TEMPLATE,   # legacy AVT-only template literal
        avt_run_tag,             # multi-run tag: AVT0001_AVT0002 for names
        norm,                    # lowercased/stripped cell value
        build_merged_cell_lookup,  # (row,col)->value dict for merged ranges
        get_cell_value,          # merged-aware cell read
        read_user_field,         # 'User:' label -> first word
        set_non_interactive,     # flag prompt helpers
        is_non_interactive,      # read prompt-flag state
        ask_avt_run,             # confirm / ask the run name (ANY sequencer)
        inherit_file_avt,        # stamp rows with their own file's run
        resolve_unassigned_avt,  # GEX-style prompt for run-less rows
        ask_reference_genome,    # choose a reference genome path
        ask_sample_selection,    # keep a subset of samples / add more
        list_sample_dirs,        # every sample dir under a run (disambig)
        choose_sample_dir,       # pick ONE of those when >1 exists
        find_fastq_path_with_pref, # resolve a fastq id honouring the pick
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
import sys

# ── shared path constants ─────────────────────────────────────────────────

#: Root of the GECF FASTQ trees on the server.  Every sequencer keeps its
#: run folders under its own code folder here (see RUN_BASE_SUBPATH).
NGS_RUNS_ROOT = "/mnt/PTEGraw/AA/NGSruns"

#: Known GECF sequencer codes — the prefixes a run name carries
#: (``AVT0519``, ``MSH0024``, ``MSQ016``, ``NVQ0312`` …).
KNOWN_RUN_CODES = ("AVT", "MSH", "MSQ", "NVQ")

#: Per-code subpath below NGS_RUNS_ROOT where the run folders live.  The
#: AVT tree keeps one extra machine level (AV240401); MSH/MSQ (MiSeq) and
#: NVQ (NovaSeq) hang their runs directly under the code folder.  Onboard a
#: new sequencer by adding its code here AND to KNOWN_RUN_CODES — the
#: pattern below picks the run names up automatically.
RUN_BASE_SUBPATH = {
    "AVT": "AVT/AV240401",
    "MSH": "MSH",
    "MSQ": "MSQ",
    "NVQ": "NVQ",
}

#: Legacy AVT-only template literal — kept so existing scripts and notes that
#: interpolate the run name themselves keep working.  New code should call
#: :func:`run_base_dir`, which honours RUN_BASE_SUBPATH per sequencer.
NGS_BASE_DIR_TEMPLATE = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/{avt_run}"

#: A bare run name on any known GECF sequencer: the code + 2-6 digits
#: (optional leading zero), case-insensitive — ``AVT0519``, ``MSH0024``,
#: ``MSQ016`` (MSQ runs are 3 digits), ``NVQ0312``.
RUN_PATTERN = re.compile(
    "(" + "|".join(KNOWN_RUN_CODES) + r")0?\d{2,6}", re.IGNORECASE
)

#: Legacy name for RUN_PATTERN — imported by the generators and the smoke
#: test.  Same object; run names are now recognised on every sequencer,
#: not only AVT.
AVT_PATTERN = RUN_PATTERN

#: A real GECF library id is a leaf named ``<prefix>id<4-6 digits>`` where
#: ``<prefix>`` is the letters of the instrument's namespace only (``avid``
#: AVT, ``mhid``/``msid`` MSQ/MSH, ``nvid`` NVQ, and whatever a future sequencer uses)
#: — never digits.  Requiring only the ``id<digits>`` tail keeps the row gate
#: and the add-sample folder-basename validation portable across platforms,
#: and keeps the id shell-safe (ARC embeds it in ``fastq_<idx>_<fid>``
#: variable names in the generated script).
#: Label / legend / placeholder rows in the same columns carry strings like
#: ``"0"``, ``"GECF libr. ID (avidxxxx)"`` or ``"Loaded molarity (pM)"`` — this
#: pattern distinguishes genuine sample rows from those (they contain spaces,
#: punctuation, or no ``id<digits>`` tail).  All generators use this as their
#: row gate, so the single source of truth lives here rather than in each
#: generator.
FASTQ_ID_PATTERN = re.compile(r"^[a-z]*id\d{4,6}$", re.IGNORECASE)


def run_base_dir(run_name, runs_root=NGS_RUNS_ROOT):
    """Resolve a GECF run name to its on-disk base folder (any sequencer).

    ``run_base_dir('AVT0519')``  -> ``<root>/AVT/AV240401/AVT0519``
    ``run_base_dir('MSH0024')``  -> ``<root>/MSH/MSH0024``
    ``run_base_dir('MSQ016')``   -> ``<root>/MSQ/MSQ016`` (earlier MSQ runs
                                    sit in ``…/MSQ016metr`` on disk and are
                                    found by the prefix match below)
    ``run_base_dir('NVQ0312')``  -> ``<root>/NVQ/NVQ0312``
    (see RUN_BASE_SUBPATH)

    Resolution order:
      1. the literal path, when it exists on this machine;
      2. a prefix match among the run folders of that code — so ``MSQ016``
         finds the on-disk ``MSQ016metr``;
      3. the literal path anyway, so the caller can report "not found" at
         the expected location (e.g. the ``AVT0000`` placeholder draft).

    A value that is not a known run code is treated as a path: absolute
    paths pass through untouched, anything else is taken relative to
    ``runs_root``.  With several prefix matches the first (alphabetical)
    wins with a warning — never a crash, never a silent wrong tree.
    """
    name = str(run_name or "").strip()
    if not name:
        name = "UNKNOWN_RUN"
    m = RUN_PATTERN.match(name)
    if not m:
        if os.path.isabs(name):
            return name
        return os.path.join(runs_root, name)
    code = name[:3].upper()
    parent = os.path.join(runs_root, *RUN_BASE_SUBPATH.get(code, code).split("/"))
    literal = os.path.join(parent, name)
    if os.path.isdir(literal):
        return literal
    try:
        cands = sorted(
            e for e in os.listdir(parent)
            if e.upper().startswith(name.upper())
            and os.path.isdir(os.path.join(parent, e))
        )
    except OSError:
        cands = []
    if cands:
        if len(cands) > 1:
            shown = ", ".join(cands[:5]) + ("…" if len(cands) > 5 else "")
            print(
                f"  Warning: run '{name}' matches {len(cands)} folders under "
                f"{parent} ({shown}) — using first: '{cands[0]}'.",
                file=sys.stderr,
            )
        return os.path.join(parent, cands[0])
    return literal

#: Default ``--localcores`` / ``--localmem`` for cellranger invocations emitted
#: by the generators.  Shared so the numbers stay in lockstep across assays;
#: per-assay overrides remain possible via the generator's own constant.
#: 32 cores per the 10x recommendation ("negligible return beyond 32 cores
#: or 256G").
LOCAL_CORES = "32"
#: 126 GB = observed high-water mark for a 20k-cell GEX run (10x EC2 trials).
LOCAL_MEM = "126"


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
    """Confirm a detected run name (any GECF sequencer), or ask for one.

    In non-interactive mode a detected name is returned uppercased
    unchanged.  When *no* run name could be detected, a loud WARNING is
    printed and the placeholder ``"AVT0000"`` is returned — matching the
    existing FIXME_LOCATE_ pattern (generate a draft, edit the unknowns
    before running) — instead of failing hard, which would leave the
    ``-n`` workflow with no way to draft a script for a not-yet-named run.
    """
    if is_non_interactive():
        if detected_avt:
            return detected_avt.upper()
        print(
            "\nWARNING: no run name detected in the workbook or\n"
            "         filename (no AVT/MSH/MSQ/NVQ code found) —\n"
            "         generated script uses placeholder run 'AVT0000'.\n"
            "         Edit 'avt_run=' (and any dependent paths) in the\n"
            "         script before running.\n"
            "         To avoid: name the workbook file with the run\n"
            "         (e.g. 'AVT0123_...xlsx' or 'MSH0024_...xlsx'),\n"
            "         or run interactively.\n"
        )
        return "AVT0000"
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
        value = input(
            "Enter run name (e.g. AVT0519, MSH0024, MSQ016, NVQ0312): "
        ).strip()
        if value:
            return value.upper()
        print("Run name cannot be empty. Please try again.\n")


def inherit_file_avt(rows, file_avt):
    """Normalise parse rows to 3-tuples and stamp in the run name the FILE
    itself carries (per-file parity for the ATAC/ARC generators).

    GEX resolves each pair under its own file's run while parsing; ATAC/ARC
    return the per-file ``avt`` separately and historically stamped ONE
    representative run onto every row.  This gives each row the run named in
    its own sheet (row-level names win), so a workbook mixing e.g. an AVT
    run and an MSH run is not silently mislabelled.

    ``rows``:      list of (sample, fastq_id[, per-row avt or None]) tuples
    ``file_avt``:  run detected for the file (sheet/filename chain), or None
    Returns a new list of 3-element lists.
    """
    out = []
    for r in rows:
        r = list(r)
        while len(r) < 3:
            r.append(None)
        if not r[2] and file_avt:
            r[2] = file_avt
        out.append(r)
    return out


def resolve_unassigned_avt(rows, require_known=True):
    """Give every row whose run name could not be detected an explicit run —
    the shared GEX 'needs_avt' prompt (parity across all three generators).

    * rows that name a run are left untouched;
    * rows without a run get ONE user decision for all of them:
      interactive mode asks once (the first detected run is the default
      hint — 'Y' confirms, typing replaces);  non-interactive (``-n``)
      keeps the ``ask_avt_run()`` contract (returned 'AVT0000' placeholder
      / loud WARNING when no run name is detected anywhere).

    ``require_known=True`` (ATAC/ARC): when NO run is detected in any file,
    the caller's own run-name answer (main_avt / any_avt) already covers
    these rows via its existing fallback — no second prompt (matches the
    GEX single-prompt result in that case).
    ``require_known=False`` (GEX): always prompt (GEX has no separate
    run-name confirmation, and the prompt also fills main_avt).

    ``rows``: list of (sample, fastq_id, avt_or_None); mutated in place.
    Returns the number of rows backfilled.
    """
    unassigned = [i for i, r in enumerate(rows)
                  if not (len(r) > 2 and r[2])]
    if not unassigned:
        return 0
    known = [r[2] for r in rows if len(r) > 2 and r[2]]
    if require_known and not known:
        return 0
    hint = known[0] if known else None
    print(f"{len(unassigned)} sample(s) have no run name assigned.")
    run = ask_avt_run(hint)
    if not run:
        raise SystemExit("Error: run name cannot be empty.")
    for i in unassigned:
        r = list(rows[i])
        while len(r) < 3:
            r.append(None)
        r[2] = run
        rows[i] = r
    return len(unassigned)


def avt_run_tag(avt_names) -> str:
    """Join a set of AVT run names into one filename-safe tag.

    ``avt_run_tag(['AVT0002', 'AVT0001'])`` -> ``'AVT0001_AVT0002'``:
    deduplicated, normalised to uppercase, sorted by numeric suffix, so the
    tag is deterministic regardless of input or CLI argument order.  A
    single run returns that run name alone, so single-AVT inputs keep
    exactly the same script/folder names as before this helper existed.

    All three generators (GEX/ATAC/ARC) use it as the name suffix of the
    generated script AND of the working/archive folder whenever the inputs
    span several AVT runs (e.g. two .xlsx files, or a sheet whose rows
    name different runs) — so the provenance is readable from the name
    alone: ``cellranger-count_AVT0001_AVT0002.sh`` /
    ``AVT0001_AVT0002_user_$TIMESTAMP``.  Returns ``''`` for empty/None
    input so the caller keeps its own fallback.
    """
    def _num(name: str) -> int:
        m = re.search(r"(\d+)$", name)
        return int(m.group(1)) if m else -1
    if not avt_names:
        return ""

    uniq = sorted(
        {str(a).upper() for a in avt_names if a},
        key=lambda name: (_num(name), name),
    )
    return "_".join(uniq)


def ask_reference_genome(human_path, mouse_path, default_path=None,
                         human_label="Human ARC ATAC",
                         mouse_label="Mouse ARC ATAC",
                         extra_choices=None):
    """Prompt for a reference genome path.

    ``default_path`` is returned immediately in non-interactive mode
    (falls back to the human path); interactive mode offers the
    human/mouse choices or accepts a custom path. ``human_label`` /
    ``mouse_label`` are the display names shown for the two built-in
    choices (assay-specific, since the same function serves all the
    10xG pipelines).

    ``extra_choices`` is an optional list of ``(label, path)`` tuples
    appended after the mouse entry — e.g. the gex pipeline offers the
    combined human+mouse reference as choice 3.  Non-interactive mode
    still honors ``default_path`` (set it to an extra path to default to
    it); extras only affect the interactive menu.
    """
    if is_non_interactive():
        return default_path or human_path
    choices = [(human_label, human_path),
               (mouse_label, mouse_path)]
    for label, p in (extra_choices or []):
        choices.append((label, p))
    while True:
        print("\nAvailable reference genomes:")
        for i, (label, p) in enumerate(choices, start=1):
            print(f"  {i}. {label} — {p}")
        print("  Custom: type a full path")
        answer = input("Choice [1]: ").strip() or "1"
        if answer.isdigit() and 1 <= int(answer) <= len(choices):
            return choices[int(answer) - 1][1]
        custom = os.path.expanduser(answer)
        if os.path.isdir(custom):
            return custom
        print(f"Path '{custom}' does not exist. Try again.\n")


def _parse_keep_selection(tokens: list, n: int) -> tuple:
    """Parse a keep-selection line into (indices, names, invalid_count).

    *tokens* is the whitespace-split user input.  Purely numeric tokens are
    treated as 1-based list indices; anything else is treated as an exact
    sample name (matched by the caller).  Returns the index list (validated
    1..n, order preserved, duplicates removed), the name list (verbatim,
    deduplicated), and the number of rejected index tokens — so the caller
    can distinguish "all rejected" (show the full list again) from
    "some bad names" (show partial feedback).
    """
    idx, names, invalid = [], [], 0
    for t in tokens:
        if t.isdigit():
            v = int(t)
            if 1 <= v <= n and v not in idx:
                idx.append(v)
            else:
                invalid += 1
        else:
            normed = t.strip()
            if normed and normed not in names:
                names.append(normed)
    return idx, names, invalid


def _prompt_prefilled(prompt: str, default: str = "") -> str:
    """input() with the line pre-filled (the bash ``read -e -i`` equivalent).

    Uses the standard readline startup-hook idiom: the hook inserts *default*
    into the edit buffer just before the prompt is drawn, so the field is
    already filled (``prompt: DEFAULT`` with the cursor at the end) and the
    user just edits it and presses Enter.  The hook is cleared again in
    ``finally`` so it never leaks into the next prompt.

    Relies on a functioning readline.  On builds where readline is a stub /
    does not honour insert_text (or on non-tty stdin), we fall back to a
    bracketed-default prompt ``prompt [DEFAULT]:`` — a bare Enter ends the
    per-folder loop and a relative tail can still be joined onto the base.
    Zero external dependencies (unlike prompt_toolkit).

    ``insert_text`` + ``get_line_buffer`` are both present in standard CPython
    (Linux / macOS, libedit-backed where applicable) and in ``pyreadline3``;
    they are not, in sandbox/stub readlines, which is what the fallback
    covers.
    """
    if default and sys.stdin.isatty():
        try:
            import readline
            readline.set_startup_hook(lambda: readline.insert_text(default))
            try:
                return input(prompt + ": ")
            finally:
                readline.set_startup_hook()  # clear — no sticky prefill
        except Exception:
            pass
    return input(f"{prompt} [{default}]: ")


def ask_sample_selection(sample_fastqs: dict,
                         default_root: str = "",
                         ) -> tuple:
    """Interactive keep/drop + add-more sample selection.

    *sample_fastqs* maps sample name -> list of (fastq_id, avt_run)
    (exactly the ``defaultdict(list)`` shape the generators build), in the
    order the samples should be presented (normally sorted).  Returns
    ``(kept_names, added)`` where:

      * ``kept_names`` is a list of sample names from *sample_fastqs*, in the
        order the user confirmed them;
      * ``added`` is ``{name: [<full fastq folder path>, ...]}`` for sample
        names NOT found in the workbook.  Each entry is a FULL path to one
        sample FASTQ folder (the directory holding R1/R2 FASTQs, whose
        basename matches ``FASTQ_ID_PATTERN`` — ``avid13144``,
        ``mhid00491``, ``nvid06797``, …); a bare id is rejected at the
        prompt.  These paths are self-referential — the generator uses
        them as-is, no run-tree walk involved.

    *default_root* is prefilled into the "full path to the sample FASTQ
    folder" prompt (bash ``read -e -i`` style).  Generators pass
    ``NGS_RUNS_ROOT`` (``/mnt/PTEGraw/AA/NGSruns``): an added sample can
    come from ANY sequencer tree, so the prefill stops at the level all of
    them share (``…/AVT/…``, ``…/MSH/…``, ``…/NVQ/…`` live below it) —
    never the archive/rsync destination.  Blank input simply stops the
    per-sample prompt (no sample is added for that name).

    Non-interactive mode (``set_non_interactive(True)`` or a non-tty stdin)
    returns everything and skips all prompts, so ``-n`` runs stay CI-safe
    (the ``added`` dict is empty in that path).
    """
    if is_non_interactive():
        return list(sample_fastqs.keys()), {}

    unique = list(dict.fromkeys(sample_fastqs.keys()))  # input order, dedup
    n = len(unique)

    nmap = {i: name for i, name in enumerate(unique, 1)}
    iname = {name: i for i, name in nmap.items()}

    # ---- keep/drop ----
    print("\nDetected samples (type 's' to see all with their FASTQs):")
    for i, name in enumerate(unique, 1):
        print(f"  {i}. {name}  ({len(sample_fastqs[name])} FASTQ)")

    kept = []
    while True:
        answer = input(
            "Keep which samples? (numbers, e.g. '1 3'; exact names; blank = all): "
        ).strip()
        if answer == "s":
            for i, name in enumerate(unique, 1):
                print(f"  {i}. {name}: " + ", ".join(fq for fq, _ in sample_fastqs[name]))
            continue
        tokens = answer.split()
        if not tokens:
            kept = list(unique)
            print(f"All {n} sample(s) will be kept.")
            break
        idxs, names, invalid = _parse_keep_selection(tokens, n)
        sel = [nmap[i] for i in idxs]
        unmatched = [x for x in names if x not in iname]
        sel += [x for x in names if x in iname]
        seen = set()
        sel = [x for x in sel if not (x in seen or seen.add(x))]
        if unmatched or invalid:
            bad = list(unmatched) + [f"index {i}" for i in (set(int(t) for t in tokens if t.isdigit()) - set(range(1, n + 1)))]
            print("Could not match: " + ", ".join(bad) + "\n", file=sys.stderr)
            for x in unmatched:
                print(f"  unknown sample name: {x}", file=sys.stderr)
            continue
        if not sel:
            print("Nothing selected — enter at least one number or name.\n")
            continue
        kept = sel
        for name in kept:
            print(f"  keep: {name}")
        n_dropped = n - len(kept)
        if n_dropped:
            print(f"  drop: {n_dropped} sample(s) — they will not appear in the generated script.")
        else:
            print("  drop: (none)")
        print()
        while True:
            conf = input("Confirm this sample set? [y/N]: ").strip()
            if conf and conf[0].lower() == "y":
                break
            if conf and conf[0].lower() == "n":
                break
            print("Answer can only be y or N. Please try again.")
        if conf and conf[0].lower() == "y":
            break

    # ---- add more (optional, repeatable) ----
    added = {}
    while True:
        while True:
            cont = input("Add another sample that is NOT in the Excel sheet? [y/N]: ").strip()
            c = cont[0].lower() if cont else ""
            if c == "y" or c == "n":
                break
            print("Answer can only be y or N. Please try again.")
        if c == "y":
            while True:
                sname = input("Sample name: ").strip()
                if not sname:
                    print("Sample name cannot be empty.\n")
                    continue
                if sname in sample_fastqs or sname in added:
                    print(f"Sample '{sname}' is already in the selection.\n")
                    continue
                break
            # ── Full path to each sample's FASTQ folder ──────────────────
            # Samples from ANY sequencer work here: the folder just has to
            # be the one holding R1/R2 FASTQs (basename ``avidNNNNN`` or
            # ``mhidNNNNN`` / ``nvidNNNNN`` / …, i.e. FASTQ_ID_PATTERN).
            # A bare id is NEVER accepted — it is not a location, and
            # tree-walking to find it fails for samples whose demux does
            # not follow the standard layout — the user points at the
            # folder itself (e.g. ``…/AVT0548/…/avid13144`` or
            # ``…/MSH0024/fastq_BCL-Convert/<attempt>/mhid00491``).
            # ``default_root`` is the prefill (generators pass NGS_RUNS_ROOT
            # — the added sample may come from ANY sequencer tree, so we
            # never guess deeper); relative tails are completed against it,
            # so a full path from any other location is accepted verbatim.
            base_hint = (default_root or "").rstrip(os.path.sep)
            paths = []
            while True:
                # Field comes pre-filled (bash read -e -i style) so the user
                # types only the tail: ``…/AVT0548/Samples/name/avidNNNNN``.
                # A bare Enter (or Enter on the untouched prefill) ends this
                # sample's folder list — same as the old blank stop.
                if base_hint:
                    p = _prompt_prefilled(
                        "Full path to the sample's FASTQ folder",
                        base_hint).strip()
                    if p == base_hint:
                        p = ""
                else:
                    p = input(
                        "Full path to the sample's FASTQ folder "
                        "(a dir like avid13144 / mhid00491 that holds the "
                        "R1/R2 fastqs; blank ends this sample): "
                    ).strip()
                if not p:
                    break
                p = os.path.expanduser(p)
                if (base_hint and os.path.sep in p
                        and not p.startswith(os.path.sep)):
                    # Relative tail — complete it against the run base.
                    p = os.path.normpath(os.path.join(base_hint, p))
                else:
                    p = os.path.abspath(p)
                if not os.path.isdir(p):
                    print(f"  '{p}' is not an existing directory — try again.\n")
                    continue
                if FASTQ_ID_PATTERN.fullmatch(os.path.basename(p)):
                    paths.append(p)
                    print(f"  added: {p}")
                else:
                    print(f"  The folder name '{os.path.basename(p)}' does "
                          "not look like a fastq ID (letters + 'id' + "
                          "4-6 digits, e.g. avid13144, mhid00491) — try "
                          "again.\n")
            if not paths:
                print("No FASTQ folder given — sample not added.\n")
                continue
            added[sname] = paths
            print(f"  {sname}: {len(paths)} FASTQ folder(s)")
            print()
        elif c == "n":
            break
        else:
            print("Answer can only be y or N. Please try again.")

    return kept, added


# ── FASTQ path resolution (sample-level disambiguation) ─────────────────────────
#
# Run folders on every GECF sequencer frequently land the same sample name
# in *multiple* sibling sub-folders (the AVT tree under fastq_bases2fastq/
# holds one sub-folder per demux attempt — corrected indexes, wrong indexes,
# …; MiSeq/NovaSeq trees use their own per-run layouts).  A bare
# fastq id alone does not pin down a physical location.  The helpers below
# let the caller:
#
#   * ``list_sample_dirs``         — return every sample dir that exists under
#                                    the run for a given sample name (ALL
#                                    of them, never "first walk hit");
#   * ``choose_sample_dir``        — collapse that candidate list to ONE dir
#                                    (no-op when exactly one candidate;
#                                    interactive lettered menu when more than
#                                    one; first candidate + warning in
#                                    non-interactive mode);
#   * ``find_fastq_path_with_pref``— resolve one fastq id to a full
#                                    path under the run base (see
#                                    ``run_base_dir``), honouring the chosen
#                                    sample dir first, falling back to a walk.


def list_sample_dirs(avt_run, sample_name):
    """Return every sample dir that exists under the run tree (any sequencer).

    *``avt_run``* is any GECF run name (AVT/MSH/MSQ/NVQ…); its base folder
    is resolved via :func:`run_base_dir`, which knows each sequencer's
    tree layout (including the ``MSQ016`` → ``MSQ016metr`` on-disk form).

    The layout can be either of two depths:

    * ``base/Samples/<name>``            (simple / older runs)
    * ``base/<folder>/<run>/Samples/<name>``  (the AVT0551 style, where
      ``base/fastq_bases2fastq/`` holds one ``<run>`` sub-folder per
      demux attempt — corrected indexes, wrong indexes, …)

    Returns a sorted list of absolute paths (deterministic order across
    runs / machines).  Empty list if the run base does not exist on this
    machine or the sample name is not present.  A sample whose name recurs
    under several ``<run>`` sub-folders — the wrong-index /
    corrected-index case — gets one path per occurrence.
    """
    base = run_base_dir(avt_run)
    if not os.path.isdir(base):
        return []

    def _try(path):
        cand = os.path.join(path, sample_name)
        if os.path.isdir(cand):
            out.append(cand)

    out: list = []
    # Depth A: base/Samples/<name>
    _try(os.path.join(base, "Samples"))
    # Depth B: base/<mid>/<run>/Samples/<name>  and
    #          base/<run>/Samples/<name>   (flat two-level layout)
    try:
        mid_entries = sorted(os.listdir(base))
    except OSError:
        return sorted(out)
    seen_mid: set = set()
    for mid in mid_entries:
        mid_path = os.path.join(base, mid)
        if not os.path.isdir(mid_path):
            continue
        # flat two-level: base/<run>/Samples/<name>
        _try(os.path.join(mid_path, "Samples"))
        # three-level: base/<mid>/<run>/Samples/<name>
        try:
            run_entries = sorted(os.listdir(mid_path))
        except OSError:
            continue
        for run in run_entries:
            run_path = os.path.join(mid_path, run)
            if os.path.isdir(run_path):
                _try(os.path.join(run_path, "Samples"))
    # de-dup (same absolute path could be reached twice via symlinks)
    dedup = []
    for p in out:
        rp = os.path.realpath(p)
        if rp not in seen_mid:
            seen_mid.add(rp)
            dedup.append(p)
    return sorted(dedup)


def choose_sample_dir(sample_name, candidates):
    """Pick one dir from *candidates* for *sample_name*.

    * Empty list → ``None`` (caller reports unresolved).
    * Exactly one candidate → returned immediately, no prompt, no output —
      the 99% case where the sample name pins the location uniquely.
    * More than one:
        - non-interactive: first candidate returned, with a warning on
          stderr naming all alternatives (so the user can grep the
          generated script and flip to a different candidate);
        - interactive: lettered menu (a/b/c…) with full paths plus the
          parent run folder so the choice is unambiguous.
    """
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    if is_non_interactive():
        chosen = candidates[0]
        alts = " ; ".join(candidates[1:])
        print(
            f"Warning: sample '{sample_name}' exists in {len(candidates)} "
            f"locations; auto-picked '{chosen}'.\n"
            f"If that is the wrong run, edit the path variable(s) to one of:\n"
            f"    {alts}",
            file=sys.stderr,
        )
        return chosen
    print(
        f"\nSample '{sample_name}' exists in {len(candidates)} locations — "
        f"pick the run sub-folder to use:"
    )
    for i, p in enumerate(candidates):
        letter = chr(ord("a") + i)
        parent = os.path.basename(os.path.dirname(os.path.dirname(p)))
        print(f"  {letter}. {p}")
        print(f"      (parent run folder: {parent})")
    while True:
        answer = input("Pick (a/b/c… or 1/2, or a full path): ").strip()
        if not answer:
            continue
        a = answer.lower()
        if len(a) == 1 and "a" <= a <= "z":
            i = ord(a) - ord("a")
            if 0 <= i < len(candidates):
                return candidates[i]
        if answer.isdigit():
            i = int(answer) - 1
            if 0 <= i < len(candidates):
                return candidates[i]
        raw = os.path.expanduser(answer)
        raw_real = os.path.realpath(raw)
        for p in candidates:
            if os.path.realpath(p) == raw_real:
                return p
        print(f"  No candidate matches '{answer}'. Try again.\n")


#: Directory names that never hold demux'd fastq leaves — skipped during
#: full-tree walks so we do not descend into raw BCL / raw basecall / QC
#: trees (huge, network-mounted, no demux'd samples inside).  The dated
#: raw-flowcell folders (``20260630_AV240401_AVT0551``-style, AVT/NextSeq
#: layout) are pruned by pattern instead.
_WALK_PRUNE_NAMES = frozenset({
    "bcl", "basecalls", "filter", "stats", "alignment", "location",
    "samplemetrics_and_jsonfiles", "script_logfiles", "multiqc_data",
    "info",
})
_RAW_FLOWCELL_DIR = re.compile(r"^\d{8}")


def find_fastq_path_with_pref(avt_run, fastq_id, preferred_sample_dir=None):
    """Resolve ``fastq_id`` to a real directory path.

    *``avt_run``* is any GECF run name (AVT/MSH/MSQ/NVQ…) resolved to its
    base folder via :func:`run_base_dir`.

    Search order:
      1. ``<preferred>/<fastq_id>`` (direct hit — the 99% case).
      2. Walk under *preferred* (handles ``<root>/<sub>/<fastq_id>`` and
         any layout the user pointed at by hand, including roots that live
         OUTSIDE the run tree).
      3. Walk under the run base (``run_base_dir(avt_run)``) — covers all
         known layouts:
         * AVT/NextSeq:  ``<run>/fastq_bases2fastq/<attempt>/Samples/<sample>``
         * MiSeq/NovaSeq: ``<run>/fastq_BCL-Convert/<attempt>/<project>/<mid>``
         * anything else that keeps the fastq-id leaf under the run folder.

    If the fastq id occurs in MORE THAN ONE place (the AVT0551
    ``1_wrong_index…``/``2_corrected_indexes…`` pair, or a MiSeq run with
    ``1_original_wrong_index``/``2_corrected_index_given``), the first two
    candidates are disambiguated with :func:`choose_sample_dir` (lettered
    menu when interactive; auto-pick of the first + warning otherwise)
    instead of silently taking the first walk hit.

    Walks early-exit as soon as two hits are found (a third copy of the
    same id is a pathological case), so single-hit runs cost only the walk
    up to the first match, and skip raw BCL / basecall / QC trees
    (see ``_WALK_PRUNE_NAMES``) to stay fast on the network-mounted run
    folders.

    Returns the resolved path or ``None`` (caller emits a FIXME_LOCATE_
    marker in the generated script).
    """
    preferred = (preferred_sample_dir or "").strip()

    def _direct(d):
        t = os.path.join(d, fastq_id)
        return t if os.path.isdir(t) else None

    def _find_hits(root, limit=2):
        """Collect up to *limit* dirs named ``fastq_id`` under *root*."""
        found = []
        for r, dirnames, _ in os.walk(root):
            kept = []
            for d in dirnames:
                if d in _WALK_PRUNE_NAMES:
                    continue
                if _RAW_FLOWCELL_DIR.match(d):
                    continue
                if d == fastq_id:
                    found.append(os.path.join(r, d))
                    if len(found) >= limit:
                        return found
                    continue  # leaf: do not descend into it
                kept.append(d)
            dirnames[:] = kept
        return found

    if preferred and os.path.isdir(preferred):
        hit = _direct(preferred)
        if hit:
            return hit
        hits = _find_hits(preferred)
        if hits:
            if len(hits) == 1:
                return hits[0]
            return choose_sample_dir(fastq_id, hits)

    base = run_base_dir(avt_run)
    if not os.path.isdir(base):
        return None
    if preferred and os.path.realpath(preferred) == os.path.realpath(base):
        # already searched; skip
        return None
    hits = _find_hits(base)
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    return choose_sample_dir(fastq_id, hits)
