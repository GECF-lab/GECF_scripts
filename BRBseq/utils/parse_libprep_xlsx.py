#!/usr/bin/env python3
# /// script
# # requires-python = ">=3.9"
# # dependencies = ["openpyxl>=3.1"]
# ///
"""parse_libprep_xlsx.py — build a BRB-seq sample/barcode manifest from a
Mercurius HS BRB-seq lib-prep Excel workbook.

The lib-prep sheet records, per sample, the well of the Mercurius
Oligo-dT (barcoded adapter) plate where that sample was barcoded (e.g. 'A11').
This script joins that well to its 14-nt CB sequence from a well-matched
barcode plate file and emits the two-column manifest that
prepare_barcodes.sh expects:

    Sample_ID<TAB>Barcode        (header row included)

Which column becomes Sample_ID:
  * The SAMPLE NAME column is preferred — a header containing 'sample name'
    (e.g. 'GECF RNA sample name' / 'RNA sample name').
  * If no name column exists in the chosen table, the sample CODE column is
    used instead (header containing 'sample code' / 'sample id', or any
    header containing 'sample') — with a loud WARNING, so the fallback is
    never silent.
  * Names are SANITISED into safe IDs before use as column headers of the
    count matrix:
      - trim; runs of whitespace -> single '_'
      - common separators '/', '\\', '-', '.' -> '_'
      - everything outside [A-Za-z0-9_] removed
      - collapsed '_++' -> '_', leading/trailing '_' stripped, cap 64 chars
    e.g. 'Patient 1 / cohort A' -> 'Patient_1_cohort_A',
         'PBMC (fresh), batch 2' -> 'PBMC_fresh_batch_2'
  * Two different raw names that sanitise to the SAME id are FATAL (it is
    then impossible to tell wells apart in the count matrix). The original
    raw name is preserved in the provenance file either way.

Plate file format (two columns, tab or comma separated, with or without a
header):  well position + barcode sequence.  Well cells are matched
permissively: 'A01' == 'A1' == 'a01'; barcodes are upper-cased and must be
[ACGT]{14}.  Any layout with one position + one 14-nt barcode per line is
accepted (e.g. the 'barcodes_96_*_brb_STAR_solo.txt' well-matched files).

The workbook's 'Kit ref' row records the plate FAMILY: the same well
encodes a different 14-nt CB on a V5A plate vs a V5D plate.  If the file
named with --plate disagrees with the Kit ref, that is a LOUD WARNING (not
an abort): an explicit --plate is authoritative — whoever names the file
is expected to know which plate was used.  Real-world context (2026-10-16):
the PHX537/Eleni workbook shipped with its kit-ref cell mislabeled 'V5A,
Part #10513' even though the V5D plate was used — the cell was corrected
the same day (original backup in ~/.hermes/cache/scratch).  The provenance
TSV records the file actually used, for post-run audit.

Sheet selection:  the first sheet (top to bottom) containing a header row
with BOTH a sample column (name or code) and a well-position column.  This
finds the 'Libraries analysis' table even when other sample tables appear
earlier. Empty sample rows after the data are ignored; a sample whose
NAME is missing is fatal (no way to name its column); a sample whose WELL
CELL is empty or not a valid position is FATAL — the lib-prep workbook must
contain a verified well position for every sample.

Kit ref cross-check (prevents using a mis-labelled plate file):

    The template's "Kit ref:" row (e.g. "Barcoded Oligo-dT Adapters Module,
    V5A, Part #10513") records which barcode family was actually used.
    The same well position encodes a different CB in V5A vs V5D, so the
    wrong plate file silently tags every sample to the wrong barcode:

  * --plate omitted (default): the script reads the Kit ref, determines the
    barcode family (V5A / V5D), and scans the plate directory for the
    matching barcodes_96_<family>_*.txt file.  The directory is the
    hardcoded PLATE_DIR (prod: /home/gecf/NGS_pipelines/BRB-seq/BRB_plates);
    pass --plate explicitly to use a file in some other location.

  * --plate given: an explicit selection is AUTHORITATIVE — the file is
    used as-is.  If its family disagrees with the workbook's Kit ref this
    is a LOUD WARNING on stderr (visible in the run log), not an abort:
    whoever names the plate explicitly is expected to know which plate
    was actually used.
    WARNING: the same well position encodes a DIFFERENT CB on the two
    families, so a wrongly-chosen file mis-tags every sample.

Exit codes: 0 ok, 1 bad usage/input, 2 parse failure, 3 bad plate file.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

WELL_RE = re.compile(r"^\s*([A-Ha-h])\s*(\d{1,2})\s*$")
BARCODE_RE = re.compile(r"^[ACGT]{14}$")
MAX_ID_LEN = 64
# Matches a barcode family token (V5A10, V5A, V5D, ...) and captures the
# A/D discriminator.  Used for both the workbook Kit ref and the plate file
# name, which use the same vocabulary.
KIT_FAMILY_RE = re.compile(r"(?<![A-Za-z])V5([AD])\d*", re.IGNORECASE)
# Plate directory (prod server): the well-matched barcode plates
# (barcodes_96_V5A10_brb.txt / barcodes_96_V5D_brb_STAR_solo.txt) live here,
# not next to the lib-prep workbooks.  Used when --plate is omitted and the
# plate is inferred from the workbook's 'Kit ref' row.
PLATE_DIR = os.environ.get(
    "BRB_PLATE_DIR", "/home/gecf/NGS_pipelines/BRB-seq/BRB_plates"
)


# ---------------------------------------------------------------- helpers

def _norm_well(text: str) -> str:
    m = WELL_RE.match(text)
    if not m:
        return ""
    return f"{m.group(1).upper()}{int(m.group(2)):02d}"


def sanitize_sample_id(raw: str) -> str:
    """Turn a free-text sample name into a safe, unique-ish column id.

    Rules (see module docstring): whitespace runs and common separators
    ('/', '\\', '-', '.') become '_', anything outside [A-Za-z0-9_] is
    dropped, '_' runs collapse, ends are stripped, length capped at 64.
    """
    s = re.sub(r"\s+", "_", str(raw).strip())
    s = re.sub(r"[/\\\-\.]+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_]", "", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:MAX_ID_LEN].rstrip("_")


def _header_score(key: str) -> int:
    """Rank a sample-header candidate: name > code/id > plain 'sample'."""
    if "name" in key:
        return 0
    if "code" in key or "id" in key:
        return 1
    return 2


def _find_header_row(ws) -> tuple[int, dict]:
    """Return (header_row, info) for the first header row that has BOTH a
    sample column AND a well-position column, else (-1, {}).

    info = {sample_col, sample_header, well_col, well_header, kind}, where
    kind is 'name' or 'code' (drives the fallback warning).
    Within one row, the header containing 'name' wins over code/id/plain.
    Matching rules (case-insensitive substring):
      sample col : 'GECF RNA sample name' > 'GECF RNA sample code' > 'sample'
      well col   : cell text contains 'well' (e.g. 'Mercurius Oligo-dT primer
                   well'); kit-ref rows without 'well' never match.
    """
    for row in ws.iter_rows():
        best_sample: tuple[int, int] | None = None  # (score, column)
        best_sample_text = ""
        well_col = 0
        well_text = ""
        for cell in row:
            if cell.value is None:
                continue
            key = str(cell.value).strip().lower()
            if not key:
                continue
            if "sample" in key:
                score = _header_score(key)
                if best_sample is None or score < best_sample[0]:
                    best_sample = (score, cell.column)
                    best_sample_text = key
            if "well" in key and not well_col:
                well_col = cell.column
                well_text = key
        if best_sample is not None and well_col:
            return row[0].row, {
                "sample_col": best_sample[1],
                "sample_header": best_sample_text,
                "well_col": well_col,
                "well_header": well_text,
                "kind": "name" if best_sample[0] == 0 else "code",
            }
    return -1, {}


def extract_researcher_first_name(ws) -> str:
    """Scan a workbook sheet for a 'Researcher:' label in the
    'General description' block (or anywhere on the sheet — first hit wins)
    and return the value cell's first word (first name).  Empty string when
    no such row is present.
    """
    def _normkey(v):
        return re.sub(r"[^a-z]", "", str(v).lower()) if v is not None else ""

    for r in range(1, ws.max_row + 1):
        for c in range(1, min(ws.max_column, 20) + 1):
            if _normkey(ws.cell(r, c).value) != "researcher":
                continue
            # value: first non-empty cell to the right in this row
            for c2 in range(c + 1, min(ws.max_column, 20) + 1):
                v = ws.cell(r, c2).value
                if v is not None and str(v).strip():
                    return re.split(r"\s+", str(v).strip())[0]
            return ""
    return ""


# ------------------------------------------------------- kit ref extraction

def _family_from_text(text: str) -> str:
    """Return 'V5A'/'V5D' if text contains a barcode-family token, else ''."""
    m = KIT_FAMILY_RE.search(text)
    return f"V5{m.group(1).upper()}" if m else ""


def extract_kit_family(wb) -> tuple[str, str]:
    """Read the 'Kit ref:' value from any sheet and return
    (family, raw_text).

    family is 'V5A' or 'V5D', or '' if no Kit ref row was found or the
    value carries no recognisable barcode-family token.
    raw_text is the first non-empty value cell next to the label (for
    diagnostics; '' if none).
    """
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for cell in row:
                if cell.value is None:
                    continue
                key = re.sub(r"[^a-z]", "", str(cell.value).lower())
                if key != "kitref" and not (key.startswith("kit") and "ref" in key):
                    continue
                # value: first non-empty cell to the right (within 4 cols)
                for c2 in range(cell.column + 1,
                                min(cell.column + 5, ws.max_column + 1)):
                    v = ws.cell(cell.row, c2).value
                    if v is not None and str(v).strip():
                        text = str(v).strip()
                        return _family_from_text(text), text
                # look a few rows below (label/value may be wrapped)
                for r2 in range(cell.row + 1, min(cell.row + 5, ws.max_row + 1)):
                    for c2 in range(cell.column + 1,
                                    min(cell.column + 5, ws.max_column + 1)):
                        v = ws.cell(r2, c2).value
                        if v is not None and str(v).strip():
                            text = str(v).strip()
                            fam = _family_from_text(text)
                            if fam:
                                return fam, text
                            break
    return "", ""


def _infer_plate_from_kit(xlsx_dir: str, family: str) -> str:
    """Find a plate file in xlsx_dir matching the kit family, or '' if none.

    family is 'V5A' or 'V5D'.  Prefers a plate file whose name contains the
    exact family token ('barcodes_96_V5A10_brb.txt' -> 'V5A' matches).
    """
    try:
        entries = sorted(os.listdir(xlsx_dir))
    except OSError:
        return ""
    for e in entries:
        if not e.endswith(".txt"):
            continue
        if not (e.lower().startswith("barcodes_96_") or e.lower().startswith("barcode_")):
            continue
        if family in e:  # 'V5A' / 'V5D' as it appears in the filename
            return os.path.join(xlsx_dir, e)
    return ""


# --------------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--xlsx", required=True,
                   help="Mercurius lib-prep .xlsx workbook")
    p.add_argument("--plate", default=None,
                   help="well-matched barcode file (well + 14-nt CB per line). "
                        "If omitted (default), the plate is inferred from the "
                        "workbook's 'Kit ref' row by scanning the plate "
                        "directory (BRB_PLATE_DIR; prod: "
                        "/home/gecf/NGS_pipelines/BRB-seq/BRB_plates) for a "
                        "matching barcodes_96_<family>_*.txt file.")
    p.add_argument("--out", required=True,
                   help="output manifest TSV (Sample_ID<TAB>Barcode, with header)")
    p.add_argument("--provenance", default="",
                   help=("optional provenance TSV (Sample_ID, Sample_raw, Well, Barcode); "
                         "default: same dir as --out, name '<out>_libprep_provenance.tsv'"))
    p.add_argument("--run-name-out", default="",
                   help=("optional path: the researcher first name (from the workbook's "
                         "'Researcher' row) is written to this file (one line, no header) "
                         "for the pipeline to use as the run name. No line is written "
                         "when the row is missing."))
    args = p.parse_args(argv)

    if not os.path.isfile(args.xlsx):
        print(f"ERROR: The library-prep workbook this step needs could not be found at: {args.xlsx}. Check that the file name and folder are spelled exactly right (file names are case-sensitive) and that it is still there, then re-run. Nothing else needs to change — the workbook is the only missing piece.", file=sys.stderr)
        return 1

    # ---- open workbook (needed before plate selection for kit ref) -------
    try:
        import openpyxl
    except ImportError:
        print(
            "ERROR: This step reads the .xlsx workbook with a program called openpyxl, and it is not available in the Python that is currently active. This is an environment-setup issue, not a problem with your workbook. Activate (or create) the pipeline's Conda environment — for example 'conda activate brbseq' — and re-run; running the whole pipeline the usual way already puts you in the right environment.",
            file=sys.stderr,
        )
        return 2

    wb = openpyxl.load_workbook(args.xlsx, data_only=True)

    # ---- kit ref → barcode family (drives plate selection) ---------------
    kit_family, kit_raw = extract_kit_family(wb)
    if kit_raw:
        fam_note = f"  [family: {kit_family}]" if kit_family else "  [no barcode family recognised]"
        print(f"parse_libprep_xlsx: kit ref   = {kit_raw}{fam_note}")
    else:
        if args.plate:
            # Pipeline behaviour (v1.4.0+): --plate is mandatory and
            # authoritative, so a missing 'Kit ref' row is NOT an error —
            # just skip the family cross-check and a one-line record,
            # nothing else. (No stderr, no NOTE.)
            print("parse_libprep_xlsx: kit ref   = (row not found — no family cross-check possible; using the given --plate file as-is)")
        else:
            print("NOTE: no 'Kit ref' row was found in this workbook, so the pipeline cannot cross-check that the barcode plate being used matches the plate family that was actually used in the lab (for example V5A vs V5D). If you are telling the pipeline which plate file to use yourself, that is fine and you can keep going. If you normally rely on the pipeline picking the plate automatically from the workbook, add or check the workbook's 'Kit ref' row so it records the plate family used.",
                  file=sys.stderr)

    # ---- determine plate file -------------------------------------------
    plate_path = args.plate
    if not plate_path:
        # infer from kit family
        if not kit_family:
            print(
                "ERROR: I could not work out which barcode plate to use: no plate file was named explicitly, and the workbook has no usable 'Kit ref' row that says which plate family (e.g. V5A or V5D) was used. Either tell the pipeline the plate file directly with --plate <file>, or add/fix the workbook's 'Kit ref' row so the plate family is recorded, then re-run.",
                file=sys.stderr,
            )
            return 1
        # the plate directory is hardcoded (PLATE_DIR), not the xlsx's own dir
        plate_path = _infer_plate_from_kit(PLATE_DIR, kit_family)
        if not plate_path:
            try:
                candidates = sorted(
                    e for e in os.listdir(PLATE_DIR)
                    if e.endswith(".txt") and (e.lower().startswith("barcodes_96_")
                                               or e.lower().startswith("barcode_")))
            except OSError:
                candidates = []
            print(
                f"ERROR: The workbook's 'Kit ref' says this run used the {kit_family} plate family, but a matching barcode plate file (a file named like barcodes_96_{kit_family}_*.txt) was not found in the folder the pipeline looks in: {PLATE_DIR}.\n"
                f"  Files that do sit there: {candidates or '(none — the folder appears empty or unreadable)'}\n"
                f"  Either put the correct {kit_family} plate file into that folder, or point the pipeline at it directly with --plate <file>, then re-run.",
                file=sys.stderr,
            )
            return 1
        print(f"parse_libprep_xlsx: plate     = {plate_path}  (inferred from kit ref {kit_family})")
    else:
        if not os.path.isfile(plate_path):
            print(f"ERROR: The barcode plate file you asked to use could not be found at: {plate_path}. Check that the file name and folder are spelled exactly right (file names are case-sensitive) and that it is still there, then re-run.", file=sys.stderr)
            return 1
        # cross-check the explicit plate file against the kit ref.
        # An explicit --plate is AUTHORITATIVE: whoever names the file is
        # expected to know which plate was used.  A family disagreement is
        # a loud warning (not an abort) so a truly mislabeled Kit ref does
        # not block a correct run; the provenance TSV records the file
        # actually used, for post-run audit.
        if kit_family:
            plate_fam = _family_from_text(os.path.basename(plate_path))
            if plate_fam and plate_fam != kit_family:
                print(
                    f"!! WARNING — KIT / PLATE FAMILY MISMATCH (run continues, "
                    f"because an explicit --plate is authoritative) !!\n"
                    f"  Your plate file is from the '{plate_fam}' family "
                    f"({os.path.basename(plate_path)}), but the workbook's "
                    f"'Kit ref' row says this preparation used the "
                    f"'{kit_family}' family:\n"
                    f"    Workbook Kit ref: {kit_raw}\n"
                    f"    Plate file used : {plate_path}\n"
                    f"  RISK: the same well position is printed with a "
                    f"DIFFERENT 14-base barcode on a V5A plate than on a "
                    f"V5D plate, so if the file you named was picked by "
                    f"mistake EVERY sample will silently be tagged with "
                    f"the wrong barcode.\n"
                    f"  ACTION REQUIRED: confirm which plate family really "
                    f"was used in the lab and that {os.path.basename(plate_path)} "
                    f"is the matching file. If it genuinely is a deliberate "
                    f"re-prep where the two differ, no flag is needed — "
                    f"this warning is only a heads-up; the chosen file is "
                    f"used as-is and the provenance file records it.",
                    file=sys.stderr,
                )
        print(f"parse_libprep_xlsx: plate     = {plate_path}")

    # ---- load plate: well -> barcode ------------------------------------
    plate: dict[str, str] = {}
    with open(plate_path, encoding="utf-8", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip().rstrip("\r")
            if not line:
                continue
            parts = [c.strip() for c in re.split(r"[\t,; ]+", line) if c.strip()]
            if len(parts) < 2:
                continue  # header or malformed: skip silently (validated below)
            well, barcode = _norm_well(parts[0]), parts[-1].upper()
            if not well or not BARCODE_RE.match(barcode):
                continue
            if well in plate and plate[well] != barcode:
                print(
                    f"ERROR: The barcode plate file is internally inconsistent: well {well} is "
                    f"listed against two different barcodes ({plate[well]} on an earlier "
                    f"line, {barcode} on line {ln}). A plate file should give exactly one "
                    f"barcode per well, so this one appears to be damaged or was written "
                    f"by a different tool than expected. Please obtain the correct, complete "
                    f"plate file for this plate family and try again.",
                    file=sys.stderr,
                )
                return 3
            plate[well] = barcode
    if not plate:
        print(
            f"ERROR: The barcode plate file ({plate_path}) did not contain any lines that look like "
            "'a well position and a 14-base barcode' — it needs one such line per well, for "
            "example:  A01\tACGT...  (a single 'Well\\tBarcode' header row at the top is "
            "fine). It looks like this may not be a well-matched barcode list at all, or it "
            "is empty/corrupt. Please use the correct well-matched barcode file for this "
            "plate family, or point the pipeline at the right file with --plate <file>.",
            file=sys.stderr,
        )
        return 3

    # ---- find the sample/well table --------------------------------------
    header_row, info, chosen_sheet = -1, {}, None
    for ws in wb.worksheets:
        hr, cm = _find_header_row(ws)
        if hr > 0:
            header_row, info, chosen_sheet = hr, cm, ws.title
            break
    if header_row < 0:
        print(
            "ERROR: I could not find the sample/well table inside this workbook. It looks for "
            "any table whose header row names a 'sample' column (sample name or code) AND "
            "a 'well' column together, and found none. If this is genuinely the right "
            "lib-prep workbook, the table's column headings may have been renamed or "
            "rearranged — check that a column still says something like 'sample name' and "
            "another still says 'well', then re-run (or open the workbook and put the "
            "standard headers back).",
            file=sys.stderr,
        )
        return 2
    ws = wb[chosen_sheet]
    col_sample = info["sample_col"]
    col_well = info["well_col"]
    if info["kind"] == "code":
        print(
            f"NOTE: in the table I found (sheet '{chosen_sheet}'), I could not find a column "
            f"for the sample NAME (I found the column headed '{info['sample_header']}' instead), "
            f"so the pipeline will use each sample's CODE as its name in the results. Nothing is "
            f"lost and the run proceeds — it simply means the columns of the count table will "
            f"be labelled with sample codes rather than the friendlier, human sample names. If "
            f"you'd like to see the actual sample names instead, add an 'RNA sample name' column "
            f"to the workbook and re-run.",
            file=sys.stderr,
        )

    # ---- researcher name (only when the pipeline asked for it) -----------
    if args.run_name_out:
        researcher_first_name = extract_researcher_first_name(ws)
        parent = os.path.dirname(os.path.abspath(args.run_name_out))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(args.run_name_out, "w", encoding="utf-8") as fh:
            if researcher_first_name:
                fh.write(researcher_first_name + "\n")
                print(f"parse_libprep_xlsx: run name   = {researcher_first_name} (from 'Researcher:' row)")
        if not researcher_first_name:
            print("NOTE: this workbook has no 'Researcher' row, so there is no name to use when labelling this run automatically. The run still proceeds; if you want a specific run name, supply one with --run-name <your name>.",
                  file=sys.stderr)

    # ---- read the sample rows --------------------------------------------
    rows: list[dict] = []
    for r in range(header_row + 1, ws.max_row + 1):
        name_raw = ws.cell(r, col_sample).value
        well_raw = ws.cell(r, col_well).value
        name_s = str(name_raw).strip() if name_raw is not None else ""
        well_s = _norm_well(str(well_raw)) if well_raw is not None else ""
        if not name_s and not well_s:
            if rows:
                break  # end-of-table marker: first fully-empty row after data
            continue
        rows.append({
            "line": r,
            "name_raw": name_s,
            "well": well_s,
            "well_raw": str(well_raw).strip() if well_raw is not None else "",
        })
    if not rows:
        print(f"ERROR: In the workbook's sample table (sheet '{chosen_sheet}'), the column headers were found (row {header_row}) but there are no sample rows below them — the table appears to be empty. Fill in the sample/well rows for this run (or confirm you are using the right workbook) and re-run.",
              file=sys.stderr)
        return 2

    # ---- sanitise + join well -> barcode ---------------------------------
    used_ids: dict[str, int] = {}          # lower(id) -> line of first occurrence
    seen_id_raw: dict[str, str] = {}       # lower(id) -> raw name of first occurrence
    used_barcodes: set[str] = set()
    barcode_owner: dict[str, str] = {}     # barcode -> sample id (1 well = 1 sample)
    manifest_lines: list[str] = []
    prov_lines: list[str] = ["Sample_ID\tSample_raw\tWell\tBarcode"]
    missing_wells: list[str] = []
    for row in rows:
        raw_name = row["name_raw"]
        if not raw_name:
            print(
                f"ERROR: Row {row['line']} of the sample table (sheet '{chosen_sheet}') has a well "
                f"position filled in but no sample name (well '{row['well_raw'] or 'not filled in'}'). "
                f"Every well in the workbook must belong to a named sample — it cannot be left "
                f"bluntly empty, because that well's barcode would otherwise be unaccounted for "
                f"in the results. Please either give this sample its name, or delete this row "
                f"if it is not part of this run, then re-run.",
                file=sys.stderr,
            )
            return 2
        sid = sanitize_sample_id(raw_name)
        if not sid:
            print(
                f"ERROR: Row {row['line']} of the sample table (sheet '{chosen_sheet}') has the name "
                f"'{raw_name}', but after the pipeline strips out characters that are not allowed "
                f"in a file/column name (spaces, slashes, and similar), nothing is left. Please "
                f"rename this sample using at least one letter or number (for example a patient "
                f"or cohort identifier) — symbols like '#', '!', or '... ' on their own cannot "
                f"become a usable name — then re-run.",
                file=sys.stderr,
            )
            return 2
        if sid.lower().startswith("empty_well"):
            print(
                f"ERROR: Row {row['line']} of the sample table (sheet '{chosen_sheet}') has the name "
                f"'{raw_name}', which after cleaning becomes '{sid}'. That starts with "
                f"'empty_well', a label the pipeline reserves for the blank plate wells it "
                f"generates automatically, so using it for a real sample would make the "
                f"two indistinguishable. Please rename this sample (drop the 'empty_well' part).",
                file=sys.stderr,
            )
            return 2
        key = sid.lower()
        if key in used_ids:
            print(
                f"ERROR: Two samples in the sheet '{chosen_sheet}' end up with the SAME cleaned name: "
                f"'{seen_id_raw[key]}' (row {used_ids[key]}) and '{raw_name}' (row {row['line']}) "
                f"both become '{sid}'. The pipeline needs each sample to have its own distinct "
                f"name, because once the run has finished you can no longer tell these two apart "
                f"in the results just from their well. Please give the two samples clearly "
                f"different names (for example append a batch number or cell source to each) "
                f"and re-run.",
                file=sys.stderr,
            )
            return 2
        used_ids[key] = row["line"]
        seen_id_raw[key] = raw_name
        if not row["well"]:
            detail = "its well cell is empty" if not row["well_raw"] else \
                f"its well cell ('{row['well_raw']}') is not recognised as a valid well position — expected a letter A-H followed by a number, for example 'A01'"
            print(
                f"ERROR: Row {row['line']} of the sample table (sheet '{chosen_sheet}'), for sample "
                f"'{sid}', {detail}. The workbook must record the verified well position for "
                f"every sample (where it was barcoded on the plate). Please correct this cell "
                f"in the workbook, or delete this row if the sample is not part of this run, "
                f"then re-run.",
                file=sys.stderr,
            )
            return 2
        barcode = plate.get(row["well"])
        if barcode is None:
            missing_wells.append(f"{sid} (well {row['well_raw'] or row['well']})")
            continue
        # one well is one sample: two different samples mapped to the same
        # well are FATAL here, before alignment (rename_matrix.py would only
        # catch it after STAR, and the failure would be late and confusing)
        previous = barcode_owner.get(barcode)
        if previous is not None and previous != sid:
            print(
                f"ERROR: Well {row['well']} (barcode {barcode}) is listed for two different samples: "
                f"'{previous}' and '{sid}'. In this assay one well holds exactly one sample, "
                f"so two samples sharing a well is not a valid setup — this would make it "
                f"impossible to know which sample that well's molecules actually come from. "
                f"Please give one of these two samples its own, separate well on the plate "
                f"(and record that correctly in the workbook), then re-run.",
                file=sys.stderr,
            )
            return 2
        barcode_owner[barcode] = sid
        used_barcodes.add(barcode)
        manifest_lines.append(f"{sid}\t{barcode}")
        prov_lines.append(f"{sid}\t{raw_name}\t{row['well']}\t{barcode}")

    if not manifest_lines:
        print(f"ERROR: None of the samples in the workbook could be matched to a well on the barcode plate, so no usable sample list could be produced. "
              f"Unmatched wells: {', '.join(missing_wells) or 'none were read at all'}.\n"
              f"  This almost always means the well letters/numbers written in the workbook do not line up with the wells that exist on the plate file actually being used "
              f"(for example a V5A workbook paired with a V5D plate). The plate's wells (first 8 shown): {sorted(plate)[:8]}{'...' if len(plate)>8 else ''}.\n"
              f"  Please check that the plate file and the workbook describe the same physical plate, fix whichever side is wrong, and re-run.",
              file=sys.stderr)
        return 2
    if missing_wells:
        print(f"NOTE: the following sample(s) named in the workbook had a well position that does not exist on the plate file, so they were NOT included in the final sample list:",
              file=sys.stderr)
        for m in missing_wells:
            print(f"    {m}", file=sys.stderr)
        print(f"  The run still proceeds with the {len(manifest_lines)} sample(s) that could be matched, and nothing else is lost. "
              f"The plate's wells (first 12 shown): {sorted(plate)[:12]}. "
              f"This almost always means the well letters/numbers written in the workbook do not line up with the wells on the plate file actually being used "
              f"(for example a V5A workbook paired with a V5D plate). Check that those two describe the same physical plate, fix whichever side is wrong, and re-run to include these samples.",
              file=sys.stderr)

    # ---- empty wells (unconditional: always emitted) ---------------------
    # Every UNOCCUPIED plate well is emitted as empty_well_<WELL> (e.g.
    # empty_well_A01) so the whole plate appears in the STARsolo matrix:
    # empty reads are never dropped by the whitelist and the run gets both a
    # samples-only and a full-plate (debug) count matrix.  Samples occupying
    # a well (duplicates included) never become empty wells; the plate's
    # well ordering is preserved for the generated rows.
    n_samples = len(manifest_lines)
    n_empty = 0
    occupied = set(used_barcodes)
    for well, bc in plate.items():  # dict preserves the plate file's order
        if bc in occupied:
            continue
        name = f"empty_well_{well}"
        # belt & braces: the 'empty_well' prefix is already reserved for
        # generated IDs above, but an alias could still collide in theory
        if name in {m.split("\t")[0] for m in manifest_lines}:
            print(
                f"ERROR: There is a name clash: a real sample in the workbook has already taken the label '{name}', which the pipeline uses for one of the automatically generated blank-well columns. It cannot create a blank-well column under a name that already belongs to a sample. Please rename that sample so it no longer starts with 'empty_well'.",
                file=sys.stderr,
            )
            return 2
        manifest_lines.append(f"{name}\t{bc}")
        n_empty += 1
    print(f"parse_libprep_xlsx: empty wells -> {n_empty} "
          f"(plate {len(plate)} wells, {len(occupied)} occupied by samples; "
          "the 'empty_well' name prefix is owned by these generated IDs)")

    out_path = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Sample_ID\tBarcode\n")
        fh.write("\n".join(manifest_lines) + "\n")

    prov_path = os.path.abspath(args.provenance) if args.provenance else \
        os.path.splitext(out_path)[0] + "_libprep_provenance.tsv"
    if prov_path == out_path:
        prov_path = out_path + "_libprep_provenance.tsv"
    with open(prov_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(prov_lines) + "\n")

    print(f"parse_libprep_xlsx: {n_samples} sample(s) from sheet '{chosen_sheet}' "
          f"(header row {header_row}, sample col '{info['sample_header']}' [{info['kind']}], "
          f"well col '{info['well_header']}') + {n_empty} empty-well row(s), "
          f"plate={plate_path}")
    print(f"parse_libprep_xlsx: manifest     -> {out_path}")
    print(f"parse_libprep_xlsx: provenance   -> {prov_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
