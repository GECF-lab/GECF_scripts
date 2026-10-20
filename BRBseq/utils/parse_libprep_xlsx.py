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
[ACGT]{12,}.  Any layout with one position + one barcode per line is
accepted (e.g. the 'barcodes_96_*_brb_STAR_solo.txt' well-matched files).

The workbook's 'Kit ref' row cross-checks the plate FAMILY: the same well
encodes a different 14-nt CB on a V5A plate vs a V5D plate, so a workbook
label and plate file that disagree on the <family> token (V5A, V5D,
V5F ...) would silently mis-assign every barcode.  parse_libprep refuses
to emit a manifest in that case (exit 4) unless --force-plate is given
(deliberate re-prep / kit swap).  Real-world case this guards (2026-10-16):
the PHX537/Eleni workbook shipped with its kit-ref cell mislabeled 'V5A,
Part #10513' even though the V5D plate was used — the cell was corrected
the same day (original backup in ~/.hermes/cache/scratch).

Sheet selection:  the first sheet (top to bottom) containing a header row
with BOTH a sample column (name or code) and a well-position column.  This
finds the 'Libraries analysis' table even when other sample tables appear
earlier.  Empty sample rows after the data are ignored; a sample whose
NAME is missing is fatal (no way to name its column); a sample with an
empty well cell is skipped with a warning.

Kit ref cross-check (prevents using a mis-labelled plate file):

    The template's "Kit ref:" row (e.g. "Barcoded Oligo-dT Adapters Module,
    V5A, Part #10513") records which barcode family was actually used.
    The same well position encodes a different CB in V5A vs V5D, so the
    wrong plate file silently tags every sample to the wrong barcode:

  * --plate omitted (default): the script reads the Kit ref, determines the
    barcode family (V5A / V5D), and scans --xlsx's directory for the
    matching barcodes_96_<family>_*.txt file.  This is the recommended
    path — you can no longer pass the wrong plate file by mistake.

  * --plate given: the plate file name's family (V5A / V5D) must be
    compatible with the workbook's Kit ref, else the script exits with
    code 4 (family is a hard constraint).  --force-plate suppresses the
    check for a deliberate re-prep where the plate family intentionally
    differs from the Kit ref.

Exit codes: 0 ok, 1 bad usage/input, 2 parse failure, 3 bad plate file,
             4 kit/plate family mismatch.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

WELL_RE = re.compile(r"^\s*([A-Ha-h])\s*(\d{1,2})\s*$")
BARCODE_RE = re.compile(r"^[ACGT]{12,}$")
MAX_ID_LEN = 64
# Matches a barcode family token (V5A10, V5A, V5D, ...) and captures the
# A/D discriminator.  Used for both the workbook Kit ref and the plate file
# name, which use the same vocabulary.
KIT_FAMILY_RE = re.compile(r"(?<![A-Za-z])V5([AD])\d*", re.IGNORECASE)


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
                        "workbook's 'Kit ref' row by scanning --xlsx's directory "
                        "for a matching barcodes_96_<family>_*.txt file.")
    p.add_argument("--force-plate", action="store_true",
                   help="suppress the kit/plate family cross-check (use this to "
                        "override a deliberate re-prep where the plate family "
                        "intentionally differs from the workbook's Kit ref).")
    p.add_argument("--out", required=True,
                   help="output manifest TSV (Sample_ID<TAB>Barcode, with header)")
    p.add_argument("--empty-wells", action="store_true",
                   help="also emit every unoccupied plate well as an "
                        "empty_well_<WELL> column (old-pipeline convention, "
                        "e.g. empty_well_A01) so all plate wells appear in "
                        "the STARsolo count matrix.  In this mode the "
                        "'empty_well' name prefix is reserved for those "
                        "generated IDs.")
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
        print(f"ERROR: workbook not found: {args.xlsx}", file=sys.stderr)
        return 1

    # ---- open workbook (needed before plate selection for kit ref) -------
    try:
        import openpyxl
    except ImportError:
        print(
            "ERROR: openpyxl is not installed in the active python.  "
            "(on this host: /home/hermes/.conda/envs/brbseq/bin/python3 has it, "
            "or pass 'python' via the pipeline.)",
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
        print("WARNING: no 'Kit ref' row found in the workbook — "
              "kit/plate family cross-check and auto-plate-inference are disabled.",
              file=sys.stderr)

    # ---- determine plate file -------------------------------------------
    plate_path = args.plate
    if not plate_path:
        # infer from kit family
        if not kit_family:
            print(
                "ERROR: --plate not supplied and no usable 'Kit ref' row found in the "
                "workbook, so the barcode family cannot be determined.\n"
                "  Pass --plate explicitly, or add a 'Kit ref' row to the workbook.",
                file=sys.stderr,
            )
            return 1
        xlsx_dir = os.path.dirname(os.path.abspath(args.xlsx))
        plate_path = _infer_plate_from_kit(xlsx_dir, kit_family)
        if not plate_path:
            try:
                candidates = sorted(
                    e for e in os.listdir(xlsx_dir)
                    if e.endswith(".txt") and (e.lower().startswith("barcodes_96_")
                                               or e.lower().startswith("barcode_")))
            except OSError:
                candidates = []
            print(
                f"ERROR: no barcodes_96_{kit_family}_*.txt plate file found in "
                f"{xlsx_dir}.\n  Candidates there: {candidates or '(none)'}",
                file=sys.stderr,
            )
            return 1
        print(f"parse_libprep_xlsx: plate     = {plate_path}  (inferred from kit ref {kit_family})")
    else:
        if not os.path.isfile(plate_path):
            print(f"ERROR: barcode plate file not found: {plate_path}", file=sys.stderr)
            return 1
        # cross-check the explicit plate file against the kit ref
        if kit_family and not args.force_plate:
            plate_fam = _family_from_text(os.path.basename(plate_path))
            if plate_fam and plate_fam != kit_family:
                print(
                    f"ERROR: plate file family '{plate_fam}' "
                    f"({os.path.basename(plate_path)}) disagrees with the workbook's "
                    f"kit ref family '{kit_family}'.\n"
                    f"  Kit ref: {kit_raw}\n"
                    f"  Plate file: {plate_path}\n"
                    f"  The same well encodes a different CB in V5A vs V5D — the wrong "
                    f"plate file silently mis-assigns every barcode.\n"
                    f"  If this mismatch is intentional (e.g. a re-prep), re-run with --force-plate.",
                    file=sys.stderr,
                )
                return 4
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
                    f"ERROR: plate file: well {well} maps to both {plate[well]} and {barcode} (line {ln})",
                    file=sys.stderr,
                )
                return 3
            plate[well] = barcode
    if not plate:
        print(
            f"ERROR: plate file {plate_path}: no usable 'well + [ACGT]{{12,}}' lines found\n"
            "Expected one line per well, e.g.  A01<TAB>ACGT...  (or with header 'Position\\tBarcode').",
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
            "ERROR: no sheet found that contains BOTH a sample-name/code column "
            "and a well-position column (case-insensitive header match).",
            file=sys.stderr,
        )
        return 2
    ws = wb[chosen_sheet]
    col_sample = info["sample_col"]
    col_well = info["well_col"]
    if info["kind"] == "code":
        print(
            f"WARNING: no 'sample name' column found in sheet '{chosen_sheet}' "
            f"(sample header: '{info['sample_header']}', well header: '{info['well_header']}'). "
            "Falling back to the sample CODE as Sample_ID. Add an 'RNA sample name' "
            "column to use human names instead.",
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
            print("WARNING: no 'Researcher' row found in the workbook — "
                  "run name not inferred (pass --run-name to set it explicitly).",
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
        print(f"ERROR: sheet '{chosen_sheet}' header found (row {header_row}) but no sample rows below it.",
              file=sys.stderr)
        return 2

    # ---- sanitise + join well -> barcode ---------------------------------
    used_ids: dict[str, int] = {}          # lower(id) -> line of first occurrence
    seen_id_raw: dict[str, str] = {}       # lower(id) -> raw name of first occurrence
    used_barcodes: set[str] = set()
    manifest_lines: list[str] = []
    prov_lines: list[str] = ["Sample_ID\tSample_raw\tWell\tBarcode"]
    missing_wells: list[str] = []
    for row in rows:
        raw_name = row["name_raw"]
        if not raw_name:
            print(
                f"ERROR: sheet '{chosen_sheet}' row {row['line']}: sample name is empty "
                f"(well '{row['well_raw'] or 'n/a'}').  Provide the RNA sample name or "
                "remove the row — samples cannot be skipped silently, their well would be orphaned.",
                file=sys.stderr,
            )
            return 2
        sid = sanitize_sample_id(raw_name)
        if not sid:
            print(
                f"ERROR: sheet '{chosen_sheet}' row {row['line']}: sample name {raw_name!r} "
                f"sanitises to an empty ID.  Use names with at least one letter or digit.",
                file=sys.stderr,
            )
            return 2
        if args.empty_wells and sid.lower().startswith("empty_well"):
            print(
                f"ERROR: sheet '{chosen_sheet}' row {row['line']}: sample name {raw_name!r} "
                f"sanitises to '{sid}', which starts with the reserved 'empty_well' prefix "
                "used for the generated empty-well columns in --empty-wells mode.  "
                "Rename the sample or drop --empty-wells.",
                file=sys.stderr,
            )
            return 2
        key = sid.lower()
        if key in used_ids:
            print(
                f"ERROR: ambiguous sample names on sheet '{chosen_sheet}': "
                f"rows {used_ids[key]} and {row['line']}: {seen_id_raw[key]!r} (row {used_ids[key]}) "
                f"and {raw_name!r} (row {row['line']}) both sanitise to '{sid}'. "
                "Give them distinct names (the well alone is not recoverable downstream).",
                file=sys.stderr,
            )
            return 2
        used_ids[key] = row["line"]
        seen_id_raw[key] = raw_name
        if not row["well"]:
            detail = "empty well cell" if not row["well_raw"] else \
                f"well cell {row['well_raw']!r} is not a valid well position (expected e.g. 'A01')"
            print(f"WARNING: {sid} (row {row['line']}): {detail} — skipped (no barcode).",
                  file=sys.stderr)
            continue
        barcode = plate.get(row["well"])
        if barcode is None:
            missing_wells.append(f"{sid} (well {row['well_raw'] or row['well']})")
            continue
        # two samples sharing one well (duplicate barcoded in two wells of the
        # same sample) is allowed: STARsolo will emit two CB columns and both
        # carry this sample id
        used_barcodes.add(barcode)
        manifest_lines.append(f"{sid}\t{barcode}")
        prov_lines.append(f"{sid}\t{raw_name}\t{row['well']}\t{barcode}")

    if not manifest_lines:
        print(f"ERROR: no sample rows could be joined to the plate.  "
              f"Unmatched wells: {', '.join(missing_wells) or 'none read at all'}.",
              file=sys.stderr)
        print("  Check that the well column uses positions that exist in the "
              f"plate file ({sorted(plate)[:8]}{'...' if len(plate)>8 else ''}).",
              file=sys.stderr)
        return 2
    if missing_wells:
        print(f"WARNING: {len(missing_wells)} sample(s) had wells not in the plate file — "
              f"excluded from the manifest:", file=sys.stderr)
        for m in missing_wells:
            print(f"    {m}", file=sys.stderr)
        print(f"  Continuing with {len(manifest_lines)} sample(s).  "
              f"Plate wells (first 12): {sorted(plate)[:12]}", file=sys.stderr)

    # ---- empty wells (optional: --empty-wells) ---------------------------
    # Emit every UNOCCUPIED plate well as empty_well_<WELL> (old-pipeline
    # convention, e.g. empty_well_A01) so the whole plate appears in the
    # STARsolo matrix as zero-count columns.  Samples occupying a well
    # (duplicates included) never become empty wells; the plate's well
    # ordering is preserved for the generated rows.
    empty_lines: list[str] = []
    if args.empty_wells:
        occupied = set(used_barcodes)
        for well, bc in plate.items():  # dict preserves the plate file's order
            if bc in occupied:
                continue
            name = f"empty_well_{well}"
            # belt & braces: the 'empty_well' prefix is already reserved for
            # generated IDs above, but an alias could still collide in theory
            if name in {m.split("\t")[0] for m in manifest_lines}:
                print(
                    f"ERROR: generated empty-well ID {name!r} collides with a "
                    "real sample ID — cannot emit it.",
                    file=sys.stderr,
                )
                return 2
            empty_lines.append(f"{name}\t{bc}")
        print(f"parse_libprep_xlsx: empty wells -> {len(empty_lines)} "
              f"(plate {len(plate)} wells, {len(occupied)} occupied by samples; "
              "reserved 'empty_well' prefix is owned by these generated IDs)")
    manifest_lines.extend(empty_lines)

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

    print(f"parse_libprep_xlsx: {len(manifest_lines)} sample(s) from sheet '{chosen_sheet}' "
          f"(header row {header_row}, sample col '{info['sample_header']}' [{info['kind']}], "
          f"well col '{info['well_header']}'), plate={plate_path}")
    print(f"parse_libprep_xlsx: manifest     -> {out_path}")
    print(f"parse_libprep_xlsx: provenance   -> {prov_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
