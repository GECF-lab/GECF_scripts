"""Shared GECF Excel intake atoms.

Every 10XG cellranger template generator (ATAC, ARC, and — eventually — Multi
and Visium) shares the same low-level workbook plumbing:

* merged-cell resolution
* cell value reading with graceful fallback
* header normalisation
* the ``User:`` requester field
* the shared interactive prompts (confirm AVT run name, choose
  reference genome, select / add samples) and the readline
  path-completion setup

The per-assay generators (``cellranger_pipelines/{atac,arc}/generate_cellranger_template_*.py``)
import *these* helpers from here and layer their own sheet-selection, row-gate,
AVT-extraction, and data-model logic on top.

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
        ask_sample_selection,      # keep a subset of samples / add more
        list_sample_dirs,          # every sample dir under an AVT (disambig)
        choose_sample_dir,         # pick ONE of those when >1 exists
        find_fastq_path_with_pref, # resolve an avid id honouring the pick
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
      * ``added`` is ``{name: [<full avid folder path>, ...]}`` for sample
        names NOT found in the workbook.  Each entry is the FULL path to
        one avid folder (the directory that holds R1/R2 FASTQs, whose
        basename matches ``avid + 4-6 digits``); a bare ``avidNNNNN`` is
        rejected at the prompt.  These paths are self-referential — the
        generator uses them as-is, no AVT-tree walk involved.

    *default_root* is shown as a hint above the "full path to the avid
    folder" prompt, so the user can see the common AVT base for the main
    run without having to type it from memory.  Blank input simply stops
    the per-sample prompt (no sample is added for that name).

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
            # ── Full path to each avid folder (hint: common AVT base) ──
            # A bare ``avidNNNNN`` is NEVER accepted here — it is not a
            # location, and tree-walking to find it fails for samples whose
            # demux does not sit under the standard
            # ``<AVT base>/<…>/Samples/<name>/avidNNNNN`` layout.  The user
            # must point at the folder itself (e.g.
            # ``/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/AVT0548/…/avid13144``).
            base_hint = default_root or "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401"
            paths = []
            while True:
                # Field comes pre-filled (bash read -e -i style) so the user
                # types only the tail: ``…/AVT0548/Samples/name/avidNNNNN``.
                # A bare Enter (or Enter on the untouched base) ends this
                # sample's folder list — same as the old blank stop.
                p = _prompt_prefilled("Full path to the avid folder",
                                     base_hint).strip()
                if not p or p == base_hint:
                    break
                p = os.path.expanduser(p)
                if os.path.sep in p and not p.startswith(os.path.sep):
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
                    print(f"  The folder name '{os.path.basename(p)}' does not "
                          "look like an avid ID (avid + 4-6 digits) — try "
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
# AVT runs frequently land the same sample name in *multiple* sibling run
# sub-folders under ``fastq_bases2fastq/`` — e.g. AVT0551 carries the same
# 384 ``avidNNNNN`` folder names under BOTH ``1_wrong_index_fromAlithea_...``
# AND ``2_corrected_indexes_...``.  A bare ``avidNNNNN`` alone does not pin
# down a physical location.  The helpers below let the caller:
#
#   * ``list_sample_dirs``         — return every sample dir that exists under
#                                    the AVT run for a given sample name (ALL
#                                    of them, never "first walk hit");
#   * ``choose_sample_dir``        — collapse that candidate list to ONE dir
#                                    (no-op when exactly one candidate;
#                                    interactive lettered menu when more than
#                                    one; first candidate + warning in
#                                    non-interactive mode);
#   * ``find_fastq_path_with_pref``— resolve one ``avidNNNNN`` to a full
#                                    path, honouring the chosen sample dir
#                                    first, falling back to a full walk.


def list_sample_dirs(avt_run, sample_name):
    """Return every sample dir that exists under the AVT tree.

    The layout can be either of two depths:

    * ``base/Samples/<name>``            (simple / older runs)
    * ``base/<folder>/<run>/Samples/<name>``  (the AVT0551 style, where
      ``base/fastq_bases2fastq/`` holds one ``<run>`` sub-folder per
      demux attempt — corrected indexes, wrong indexes, …)

    Returns a sorted list of absolute paths (deterministic order across
    runs / machines).  Empty list if the AVT base does not exist on this
    machine or the sample name is not present.  A sample whose name recurs
    under several ``<run>`` sub-folders — the wrong-index /
    corrected-index case — gets one path per occurrence.
    """
    base = NGS_BASE_DIR_TEMPLATE.format(avt_run=avt_run)
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


def find_fastq_path_with_pref(avt_run, fastq_id, preferred_sample_dir=None):
    """Resolve ``fastq_id`` to a real directory path.

    Search order:
      1. ``<preferred>/<fastq_id>`` (direct hit — the 99% case).
      2. ``os.walk`` under *preferred* (handles ``<root>/<sub>/<fastq_id>``
         and any layout the user pointed at by hand, including roots that
         live OUTSIDE the AVT run tree).
      3. ``os.walk`` under the AVT run base
         (``NGS_BASE_DIR_TEMPLATE.format(avt_run=avt_run)``) — the historical
         fallback that finds samples which sit in non-Samples layouts.

    Returns the resolved path or ``None`` (caller emits a FIXME_LOCATE_
    marker in the generated script).
    """
    preferred = (preferred_sample_dir or "").strip()

    def _direct(d):
        t = os.path.join(d, fastq_id)
        return t if os.path.isdir(t) else None

    def _walk_under(root):
        for r, dirnames, _ in os.walk(root):
            for d in dirnames:
                if d == fastq_id:
                    return os.path.join(r, d)
        return None

    if preferred and os.path.isdir(preferred):
        hit = _direct(preferred) or _walk_under(preferred)
        if hit:
            return hit

    base = NGS_BASE_DIR_TEMPLATE.format(avt_run=avt_run)
    if not os.path.isdir(base):
        return None
    if preferred and os.path.realpath(preferred) == os.path.realpath(base):
        # already searched; skip
        return None
    return _walk_under(base)
