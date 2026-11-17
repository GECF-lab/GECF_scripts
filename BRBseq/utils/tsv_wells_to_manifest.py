#!/usr/bin/env python3
"""tsv_wells_to_manifest.py — BRB-seq pipeline helper (utils/).

Turns the new-style ``--barcode-tsv`` input — a 2-column TSV with
``well <TAB> sample_id`` (an optional header row is auto-detected and
skipped) — into the full-plate sample→barcode manifest the rest of the
BRB-seq pipeline consumes (the SAME contract parse_libprep_xlsx.py
produces in workbook mode).  Sample names are canonicalised to a
column-safe id with the SAME rules as workbook mode (see
``sanitize_sample_id`` below, which MUST stay in sync with
``parse_libprep_xlsx.py:sanitize_sample_id``) — so both manifest sources
give the SAME Sample_ID for the SAME name.  The original name is kept in
the provenance file's ``sample_raw`` column.

    manifest:    sample_id <TAB> barcode   (header row; plus
                    empty_well_<WELL> rows for every plate well that no
                    sample occupies, so unoccupied wells stay trackable)
    provenance:  sample_id <TAB> sample_raw <TAB> well <TAB> barcode
                    (used by the report grid to fill the plate_well column)

Barcodes come from the barcode plate file given by the pipeline's
``--plate`` (well + 14-nt CB per line, optional header row).

The script is intentionally dependency-free (stdlib only) and fail-hard:
any structural problem in the TSV or the plate file stops the run with a
message that names the offending line(s)/well(s).

Usage:
    python tsv_wells_to_manifest.py \
        --tsv        /path/sample_wells.tsv \
        --plate      /path/barcodes_96_V5D_brb_STAR_solo.txt \
        --manifest   /path/barcode_manifest.tsv \
        --provenance /path/barcode_manifest_provenance.tsv
"""

import argparse
import re
import sys

WELL_RE = re.compile(r"^[A-H][0-9]{1,2}$")
BARCODE_RE = re.compile(r"^[ACGTN]{14}$")
MAX_ID_LEN = 64


def sanitize_sample_id(raw):
    """Turn a free-text sample name into a safe, unique-ish column id.

    Rules (MUST stay in sync with utils/parse_libprep_xlsx.py:
    sanitize_sample_id, so TSV and workbook mode name the same sample the
    same way): whitespace runs and common separators ('/', '\\', '-', '.')
    become '_', anything outside [A-Za-z0-9_] is dropped, '_' runs
    collapse, ends are stripped, length capped at 64.
    """
    s = re.sub(r"\s+", "_", str(raw).strip())
    s = re.sub(r"[/\\\-\.]+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_]", "", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:MAX_ID_LEN].rstrip("_")


def norm_well(well):
    """Accept A1 as well as A01; return the canonical 2-digit form (A01),
    or None if the token is not a plate well."""
    if not WELL_RE.match(well):
        return None
    row = int(well[1:])
    if not 1 <= row <= 96:
        return None
    return f"{well[0]}{row:02d}"


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def read_plate(path):
    """Return an ordered list of (well, barcode) pairs from a plate file
    (well + 14-nt barcode per line).  A leading header row (any row whose
    first column is not a well ID, e.g. 'Name') is skipped."""
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        die(f"could not read the barcode plate file ({path}): {exc}")

    plate = []
    seen = {}
    for lineno, line in enumerate(lines, start=1):
        line = line.strip()
        if not line:
            continue
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) < 2:
            if len(fields) == 1:
                continue  # header row, e.g. "Name"
            continue
        well, barcode = norm_well(fields[0].upper()), fields[1].upper()
        if well is None:
            if lineno == 1:
                continue  # header row, e.g. "Name<TAB>B1"
            die(f"plate file line {lineno} of {path}: first column '{fields[0]}' "
                f"is not a well ID (expected one of A01..H12). The plate file must "
                f"have one line per well: well <TAB> 14-nt barcode.")
        if not BARCODE_RE.match(barcode):
            die(f"plate file line {lineno} of {path}: barcode '{fields[1]}' for well "
                f"{well} is not a valid 14-nt cell barcode (expected only A, C, G, T, "
                f"N — e.g. 'TACGTTATTCCGAA').")
        if well in seen:
            die(f"plate file {path} lists well {well} more than once "
                f"(first at line {seen[well]}, again at line {lineno}); "
                f"each well must appear exactly once.")
        seen[well] = lineno
        plate.append((well, barcode))

    if not plate:
        die(f"the barcode plate file ({path}) contains no usable well/barcode rows "
            f"(expected one line per well: well <TAB> 14-nt barcode).")
    return plate


def parse_tsv_rows(path):
    """Read the well/sample TSV.  Returns (rows, header_row_str_or_None) with
    rows = list of (lineno, well, sample_id).  Fails hard on:
    * wrong column count, empty cells,
    * a header row in any other position than line 1,
    * the deprecated sample_id<TAB>barcode layout,
    * anything else that is not a well ID in column 1 (lists the lines)."""
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        die(f"could not read the well/sample list ({path}): {exc}")

    # physical lines that are not blank
    phys = [(i, ln) for i, ln in enumerate(lines, start=1) if ln.strip()]
    if not phys:
        die(f"the well/sample list ({path}) is empty — expected 2 columns per line: "
            f"well (A01..H96) in column 1, sample_id in column 2"
            f" (an optional header row on line 1 is fine).")

    def split_fields(i, ln):
        fields = [f.strip() for f in ln.split("\t")]
        if len(fields) != 2:
            die(f"well/sample list line {i}: expected exactly 2 tab-separated columns "
                f"(well <TAB> sample_id), found {len(fields)}"
                f"{' (an empty TSV column? check for double tabs)' if len(fields) > 2 else ''}: "
                f"'{ln.strip()}'")
        return fields[0], fields[1]

    def is_old_format(i, c1, c2):
        # column 1 is not a well and column 2 reads as a 14-nt barcode
        # => the deprecated 'sample_id <TAB> barcode' layout.
        return (not WELL_RE.match(c1.upper()) and BARCODE_RE.match(c2.upper()))

    first_i = phys[0][0]
    header = None
    rows = []
    bad = []
    for i, ln in phys:
        c1, c2 = split_fields(i, ln)
        well = norm_well(c1.upper())
        if well is not None:
            if c2 == "":
                die(f"well/sample list line {i}: sample_id (column 2) is empty "
                    f"(well '{well}').")
            rows.append((i, well, c2))
        elif i == first_i and not rows and not bad:
            header = ln.strip()  # line-1 header row: allowed, skipped
        else:
            bad.append((i, c1, c2))

    if not rows:
        # Nothing mapped to a well: either an old-format file or all rows bad.
        first_i, first_ln = phys[0]
        fc1, fc2 = [f.strip() for f in first_ln.split("\t")]
        if BARCODE_RE.match(fc2.upper()):
            die("this TSV looks like the OLD 2-column layout "
                "'sample_id <TAB> barcode' (column 2 of the first row reads as a "
                f"14-nt barcode: '{fc2}'). That format is no longer supported: the "
                "barcode is now derived from the --plate file, joined by well. "
                "Convert the file to the new layout — column 1 = well (A01..H96), "
                "column 2 = sample_id — and pass the plate file with --plate. "
                "(A well→sample_id pair per line; wells left out of this list "
                "become empty_well_<WELL> rows in the output.)")
        if bad:
            first = bad[0]
            die(f"well/sample list line {first[0]}: column 1 '{first[1]}' is not a "
                f"well ID (expected one of A01..H12), and the file contains no "
                f"usable rows. Each data line must be  well <TAB> sample_id.")
        die(f"the well/sample list ({path}) contains no data rows (only a header?). "
            f"Each data line must be  well <TAB> sample_id.")

    if bad:
        listing = "; ".join(f"line {i}: column 1 '{c1}'" for i, c1, c2 in bad[:8])
        more = f" (and {len(bad) - 8} more line(s))" if len(bad) > 8 else ""
        die(f"well/sample list column 1 must be a well ID (A01..H12) — offending "
            f"line(s): {listing}{more}. If these rows are 'sample_id <TAB> "
            f"barcode', the file is in the OLD format that is no longer supported "
            f"(column 1 = well, column 2 = sample_id instead).")

    return rows, header


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True, help="well/sample list (well <TAB> sample_id)")
    ap.add_argument("--plate", required=True, help="barcode plate file (well + 14-nt CB)")
    ap.add_argument("--manifest", required=True, help="output manifest (sample_id <TAB> barcode)")
    ap.add_argument("--provenance", required=True, help="output provenance (sample_id <TAB> sample_raw <TAB> well <TAB> barcode)")
    args = ap.parse_args()

    plate = read_plate(args.plate)
    bar_of = dict(plate)
    rows, header = parse_tsv_rows(args.tsv)

    # --- canonicalise sample names (same rules as workbook mode) -----------
    clean_rows = []  # (lineno, well, raw_name, cleaned_sid)
    for i, well, raw in rows:
        sid = sanitize_sample_id(raw)
        if not sid:
            die(f"well/sample list line {i}: sample_id '{raw}' leaves no letters "
                f"or numbers after cleaning (spaces and symbols like '/', '#', "
                f"'!' become '_' or are dropped), so it cannot become a usable "
                f"name. Rename that sample with at least one letter or number "
                f"(e.g. a patient or cohort identifier) and re-run.")
        clean_rows.append((i, well, raw, sid))

    # --- reserved sample names ---------------------------------------------
    for i, well, raw, sid in clean_rows:
        if sid.lower().startswith("empty_well"):
            die(f"well/sample list line {i}: sample_id '{raw}' (cleaned: '{sid}') "
                f"starts with the reserved prefix 'empty_well' — unoccupied plate "
                f"wells are named empty_well_<WELL> automatically, so real samples "
                f"may not use it. Rename that sample and re-run.")
        if sid == "Unknown":
            die(f"well/sample list line {i}: sample_id '{raw}' (cleaned: '{sid}') "
                f"is reserved for the unassigned-reads bucket of the "
                f"read-distribution report — a real sample may not use that exact "
                f"name. Rename it and re-run.")

    # --- duplicates (well and sample), and unmappable wells -----------------
    seen_wells = {}
    unmapped = []
    for i, well, raw, sid in clean_rows:
        if well not in bar_of:
            unmapped.append((i, well))
            continue
        if well in seen_wells:
            j, other_raw = seen_wells[well]
            die(f"well {well} appears more than once: line {j} (sample "
                f"'{other_raw}') and line {i} (sample '{raw}'). One sample per "
                f"well — merge the two rows into one or move one sample to a "
                f"different, unused well, then re-run.")
        seen_wells[well] = (i, raw)

    seen_ids = {}
    for i, well, raw, sid in clean_rows:
        key = sid.lower()
        if key in seen_ids:
            prev_line, prev_well, prev_raw = seen_ids[key]
            die(f"two samples clean to the SAME name '{sid}': line {prev_line} "
                f"('{prev_raw}', well {prev_well}) and line {i} ('{raw}', well "
                f"{well}). Once cleaned (spaces → '_', other symbols dropped) the "
                f"pipeline can no longer tell these two apart. Give the two samples "
                f"clearly different names and re-run.")
        seen_ids[key] = (i, well, raw)

    if unmapped:
        listing = "; ".join(f"line {i}: well {w}" for i, w in unmapped[:12])
        more = f" (and {len(unmapped) - 12} more)" if len(unmapped) > 12 else ""
        die(f"the following well(s) from the TSV are not present in the plate "
            f"file ({args.plate}), so no barcode can be assigned: {listing}{more}. "
            f"Check the well positions in the TSV against the plate (wells are "
            f"named A01..H96), then re-run.")

    # --- emit manifest + provenance ----------------------------------------
    occupied = {w for _, w, _, _ in clean_rows}
    manifest_lines = ["sample_id\tbarcode"]
    prov_lines = ["sample_id\tsample_raw\twell\tbarcode"]
    n_empty = 0
    for i, well, raw, sid in clean_rows:  # TSV order first (samples)
        barcode = bar_of[well]
        manifest_lines.append(f"{sid}\t{barcode}")
        prov_lines.append(f"{sid}\t{raw}\t{well}\t{barcode}")
    for well, barcode in plate:  # then empty wells, in plate-file order
        if well not in occupied:
            n_empty += 1
            name = f"empty_well_{well}"
            manifest_lines.append(f"{name}\t{barcode}")
            prov_lines.append(f"{name}\t{name}\t{well}\t{barcode}")

    try:
        with open(args.manifest, "w", encoding="utf-8") as fh:
            fh.write("\n".join(manifest_lines) + "\n")
        with open(args.provenance, "w", encoding="utf-8") as fh:
            fh.write("\n".join(prov_lines) + "\n")
    except OSError as exc:
        die(f"could not write the manifest/provenance files: {exc}")

    if header is not None:
        print(f"note: header row auto-detected in the TSV and skipped: '{header}'")
    print(f"ok: {len(rows)} sample row(s) mapped to plate barcodes (plate: "
          f"{len(plate)} wells); {n_empty} unoccupied well(s) emitted as "
          f"empty_well_<WELL> rows.")
    print(f"  manifest   = {args.manifest}")
    print(f"  provenance = {args.provenance}")


if __name__ == "__main__":
    main()
