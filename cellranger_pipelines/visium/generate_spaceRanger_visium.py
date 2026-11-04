#!/usr/bin/env python3
"""generate_spaceRanger_visium.py (v1 — spaceranger 4.1.0 compatible)

Reads a Visium HD samples workbook and emits one `spaceranger count` bash block
per sample, targeting spaceranger 4.1.0.

TWO MAIN APPLICATIONS:
  1. Visium HD WITH probes    -> --transcriptome + --probe-set
  2. Visium HD 3' (WITHOUT)   -> --transcriptome only

Usage
-----
  python3 generate_spaceRanger_visium.py samples.xlsx --species human \
      --cytaimage path_to_image.tif
  python3 generate_spaceRanger_visium.py samples.xlsx --species human \
      --slide V10J25-015 --area A1
  python3 generate_spaceRanger_visium.py samples.xlsx --species human \
      --unknown-slide visium-hd
  # Preview the commands without writing anything (non-coder review step):
  python3 generate_spaceRanger_visium.py samples.xlsx --species human \
      --cytaimage img.tif --print-commands

For non-coders, run_visium.sh wraps this: validate, preview, confirm,
launch, tail, and report PASS/FAIL.
"""

import argparse
import os
import re
import shutil
import sys
import warnings

warnings.filterwarnings("ignore", message="Data Validation extension is not supported")
warnings.filterwarnings("ignore", message="Conditional Formatting extension is not supported")

try:
    import openpyxl
except ImportError:
    sys.exit("Error: 'openpyxl' is required. pip3 install openpyxl")


# On ptegsrv3 the transcriptomes and the probes live under ONE shared root,
# /home/gecf/references/spaceranger_references :
#   human/  refdata-gex-GRCh38-2024-A  (current) + old/refdata-gex-GRCh38-2020-A
#   mouse/  refdata-gex-GRCm39-2024-A  (current) + old/refdata-gex-mm10-2020-A
#   Probes/human (current CSV)  + Probes/human/OLD
#   Probes/mouse/HD (current)  + Probes/mouse/NON_HD (v1.0) + Probes/mouse/OLD (v1.0+v2.0)
#
# Canonical transcriptome root + per-species subdir (ptegsrv3 layout).
# Override with $REF_ROOT when generating for a server with a different layout
# (e.g. SV-76L-002: /home/gecf/Documents/references).
REF_ROOT = os.environ.get("REF_ROOT", "/home/gecf/references/spaceranger_references")

# --transcriptome path = REF_ROOT/<species>/refdata-gex-*-2024-A. Same
# 2024-A refdata serves BOTH applications — the probe-set vs no-probe
# distinction only changes the --probe-set argument, not the reference.
#   mouse HD_PROBES + mouse HD_3P -> <root>/mouse/refdata-gex-GRCm39-2024-A
#   human HD_PROBES + human HD_3P -> <root>/human/refdata-gex-GRCh38-2024-A
REF = {
    ("human", "HD_PROBES"): REF_ROOT + "/human/refdata-gex-GRCh38-2024-A",   # 1: Visium HD with probes
    ("mouse", "HD_PROBES"): REF_ROOT + "/mouse/refdata-gex-GRCm39-2024-A",   # 1: Visium HD with probes
    ("human", "HD_3P"):     REF_ROOT + "/human/refdata-gex-GRCh38-2024-A",   # 2: Visium HD 3' (no probes)
    ("mouse", "HD_3P"):     REF_ROOT + "/mouse/refdata-gex-GRCm39-2024-A",   # 2: Visium HD 3' (no probes)
}

# --probe-set CSV relpath (application 1 only). Application 2 (Visium HD 3')
# rows resolve to None -> the --probe-set line is NOT emitted (per 10x).
# PROBE values are RELATIVE PATHS under PROBE_ROOT (not bare basenames),
# mirroring the ptegsrv3 filesystem layout:
#   Probes/human/<csv>          — human v2.1.0 (current) is top-level human/
#   Probes/human/OLD/<csv>      — human v2.0 (stale)
#   Probes/mouse/HD/<csv>       — mouse v2.1.0 (current — mouse needs the HD/ class folder)
#   Probes/mouse/NON_HD/<csv>   — mouse v1.0 (retired chemistry, kept on ptegsrv3)
#   Probes/mouse/OLD/<csv>      — mouse v1.0 + v2.0 (stale)
# _probe_path() joins these onto the runtime root, so the canonical emitted
# path is PROBE_ROOT/human/<csv> (human) or PROBE_ROOT/mouse/HD/<csv> (mouse).
# Override PROBE_ROOT via --probe-root or $PROBE_ROOT when the probes live
# elsewhere (e.g. SV-76L-002's flat bundle ~/Documents/spaceranger-4.1.0/
# probe_sets — the _probe_path flat-join fallback handles that case).
PROBE_ROOT = os.environ.get("PROBE_ROOT", "/home/gecf/references/spaceranger_references/Probes")
PROBE = {
    ("human", "HD_PROBES"): "human/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv",
    ("mouse", "HD_PROBES"): "mouse/HD/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv",
    ("human", "HD_3P"):     None,   # Visium HD 3' -> no --probe-set
    ("mouse", "HD_3P"):     None,   # Visium HD 3' -> no --probe-set
}

CHEM_ALIASES = {
    # Application 1: Visium HD WITH probes -> HD_PROBES
    # (nonHD alias removed 2026-09-03: nonHD no longer supported, see
    #  resolve_sample error)
    "hd_probes": "HD_PROBES", "hd_probeset": "HD_PROBES",
    "visium_hd_probes": "HD_PROBES", "visium_hd_with_probes": "HD_PROBES",
    # Application 2: Visium HD 3' WITHOUT probes -> HD_3P
    # (bare "visium hd", "visium", "gex", "3p", "3'", standard, etc. all mean
    #  the no-probe application)
    "hd": "HD_3P", "gex_hd": "HD_3P", "visium_hd": "HD_3P", "visiumhd": "HD_3P",
    "gex": "HD_3P", "gex_": "HD_3P", "visium": "HD_3P", "visium_": "HD_3P",
    "3prime": "HD_3P", "3p": "HD_3P", "gex3": "HD_3P", "gex3prime": "HD_3P",
    "gex_3prime": "HD_3P", "gex_3": "HD_3P", "std": "HD_3P", "std_gex": "HD_3P",
    "std_3p": "HD_3P", "standard": "HD_3P", "standard_gex": "HD_3P",
    "visium_3prime": "HD_3P", "visium3p": "HD_3P", "visium3prime": "HD_3P",
    "visium_3p": "HD_3P", "visium_3": "HD_3P",
    "10xg_visium_hd_3": "HD_3P", "10xvisiumhd3": "HD_3P",
    "visium_hd_3": "HD_3P", "visiumhd3": "HD_3P",
    "10xg_visium_hd": "HD_3P", "10xg_visium": "HD_3P", "10xvisium": "HD_3P",
    "spatial_gex_3prime": "HD_3P", "gex_spatial_3prime": "HD_3P",
}
CANONICAL_CHEM = ("HD_PROBES", "HD_3P")

NGS_BASE_DIR_TEMPLATE = "/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/{avt_run}"
AVT_PATTERN = re.compile(r"AVT(\d{3,5})", re.IGNORECASE)

SPACERANGER_BIN = os.environ.get("SPACERANGER_BIN", "spaceranger")

unknown_slide_choices = ["visium-1", "visium-2", "visium-2-large",
                         "visium-hd", "visium-hd-large"]


# ────────────────────────────────────────────────────────────────────────
# Column aliases
# ────────────────────────────────────────────────────────────────────────

COLUMN_ALIASES = {
    "samplename": {
        "samplename", "sample", "sampleid", "samplelabel",
    },
    "tubelabel": {
        "librarytubelabel", "tubelabel", "librarylabel", "libraryname",
    },
    "fastqs": {
        "fastqs", "fastqsdir", "fastq", "fastqid",
        "gcffastqid", "gecffastqid", "gcffastqidavidxxxx",
        "gecffastqidavidxxxx", "fastqdirectories",
    },
    "chemistry": {
        "chemistry", "librarymethod", "method", "chem", "assay",
    },
    "species":     {"species", "organism", "genus_species"},
    "cytaimage":   {"cytaimage", "celltreatmentimage"},
    "image":       {"image", "omeimage", "darkimage", "colorizedimage"},
    "slide":       {"slide", "slideid"},
    "area":        {"area", "region", "roi", "foci"},
    "loupe":       {"loupe", "fiducials", "fiducialsjson", "loupealign"},
    "description": {"description", "comments", "comment"},
    "output_dir":  {"outputdir", "out", "outdir"},
    "seq_dir":     {"seqdir", "seqpath", "seq"},
    "avt_run":     {"avtrun", "avtid", "runname", "gcfindivid"},
    # S1-S5 nucleus / binning / annotation / image quality (spaceranger 4.1.0)
    "nucleus_segmentation": {"nucleussegmentation", "nucleusseg", "nucleus_seg"},
    "custom_bin_size":      {"custombinsize", "binsize", "bin_size", "bin_um", "binum"},
    "max_nucleus_diameter_px": {"maxnucleusdiameterpx", "maxnucleusdiameter",
                                "nucleusdiameterpx", "nucleusdiameter", "maxnucleuspx"},
    "nucleus_expansion_distance_micron": {"nucleusexpansiondistancemicron",
                                          "expansiondistancemicron", "nucleusexpansion",
                                          "expansiondist", "expansionmicron"},
    "custom_segmentation_file": {"customsegmentationfile", "segmentationfile",
                                 "segmentation", "segfile"},
    "cell_annotation_model":    {"cellannotationmodel", "cellannotation",
                                 "annotationmodel", "annotation"},
    "disable_cell_annotation":  {"disablecellannotation", "nocellannotation",
                                 "disableannotation"},
    "dapi_index":               {"dapiindex", "dapi", "dapichannel", "dapiidx"},
    "reorient_images":          {"reorientimages", "reorient", "reorientimg"},
    "image_scale":              {"imagescale", "micronperpixel",
                                 "micronsperpixel", "scale"},
    "override_id":              {"overrideid", "overridemetadata", "override"},
    "filter_probes":            {"filterprobes"},
    "exp_code":    {"expcode", "expc"},
    "user":        {"user", "username"},
}

_CANON_BY_NORM = {}
for _canon, _set in COLUMN_ALIASES.items():
    for _h in _set:
        _CANON_BY_NORM.setdefault(_h, _canon)


# ────────────────────────────────────────────────────────────────────────
# Normalisation helpers
# ────────────────────────────────────────────────────────────────────────

def norm_header(h):
    if h is None:
        return ""
    return re.sub(r"[^a-z0-9]", "", str(h).lower())


def norm_value(v):
    if v is None:
        return ""
    s = str(v)
    if s.strip() in ("", "0"):
        return ""
    s = s.replace("\r", " ").replace("\n", " ")
    return re.sub(r"\s+", " ", s).strip()


def canon_for_header(h):
    n = norm_header(h)
    if not n:
        return None
    return _CANON_BY_NORM.get(n)


# ────────────────────────────────────────────────────────────────────────
# Merged-cell support (ported from generate_cellranger_template_atac.py)
# ────────────────────────────────────────────────────────────────────────

def build_merged_cell_lookup(ws):
    lookup = {}
    try:
        ranges = list(ws.merged_cells.ranges)
    except Exception:
        return lookup
    for mr in ranges:
        try:
            tr, tc, br, bc = mr.min_row, mr.min_col, mr.max_row, mr.max_col
        except Exception:
            continue
        try:
            v = ws.cell(row=tr, column=tc).value
        except Exception:
            v = None
        for r in range(tr, br + 1):
            for c in range(tc, bc + 1):
                lookup[(r, c)] = v
    return lookup


def get_cell(ws, merged, r, c):
    if (r, c) in merged:
        v = merged[(r, c)]
        if v is not None:
            return v
    try:
        return ws.cell(row=r, column=c).value
    except Exception:
        return None


def find_sheet_sheet_header(path, verbose=True):
    wb = openpyxl.load_workbook(path, data_only=True)
    if verbose:
        print("  sheets:", ", ".join(wb.sheetnames), file=sys.stderr)
    if "User_x1x_copy" in wb.sheetnames:
        chosen = ["User_x1x_copy"] + wb.sheetnames
    else:
        chosen = wb.sheetnames

    best = None
    for sheetname in chosen:
        ws = wb[sheetname]
        merged = build_merged_cell_lookup(ws)
        max_row = ws.max_row or 0
        max_col = ws.max_column or 0
        if max_row == 0 or max_col == 0:
            continue
        scan_to = min(15, max_row)

        for r in range(1, scan_to + 1):
            mapping = {}
            hits = 0
            for c in range(1, max_col + 1):
                v_top = get_cell(ws, merged, r, c)
                v_prev = get_cell(ws, merged, r - 1, c) if r > 1 else None
                combined = "{} {}".format(v_prev or "", v_top or "").strip()
                canon = canon_for_header(combined)
                if canon and canon not in mapping:
                    mapping[canon] = c
                    hits += 1
            if hits >= 2:
                if best is None or hits > best[3]:
                    best = (sheetname, r, mapping, hits)
    wb.close()
    if not best:
        sys.exit("FATAL: could not detect a header row in %s" % path)
    sheetname, row_idx, mapping, hits = best
    if verbose:
        print("  using sheet %r, header row %d, %d cols resolved" % (
            sheetname, row_idx, len(mapping)), file=sys.stderr)
        for k in sorted(mapping):
            print("   %s -> col %d" % (k, mapping[k]), file=sys.stderr)
    return sheetname, row_idx, mapping


def read_data_rows(path, sheetname, header_row_idx, mapping):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[sheetname]
    merged = build_merged_cell_lookup(ws)
    max_row = ws.max_row or 0

    data_start = header_row_idx + 1
    for mr in ws.merged_cells.ranges:
        if (mr.min_row <= header_row_idx <= mr.max_row and
                any(mr.min_col <= mapping.get(k, 0) <= mr.max_col
                    for k in mapping)):
            data_start = max(data_start, mr.max_row + 1)

    out = []
    blank_streak = 0
    for r in range(data_start, max_row + 1):
        row = {}
        for canon, c in mapping.items():
            v = get_cell(ws, merged, r, c)
            row[canon] = norm_value(v)
        for k in ("samplename", "fastqs"):
            if "run details" in row.get(k, "").lower():
                blank_streak = 99
                break
        if blank_streak == 99:
            break

        if not any(row.values()):
            blank_streak += 1
            if out and blank_streak >= 2:
                break
            continue
        blank_streak = 0

        if not row.get("samplename") and not row.get("fastqs"):
            continue
        out.append(row)
    wb.close()
    return out


def parse_xlsx(path, verbose=True):
    sheetname, header_row, mapping = find_sheet_sheet_header(path, verbose=verbose)
    rows = read_data_rows(path, sheetname, header_row, mapping)
    avt_m = None
    m = AVT_PATTERN.search(path)
    if m:
        avt_m = m.group(1)
    for r in rows:
        r["_detected_avt"] = avt_m
    return rows


# ────────────────────────────────────────────────────────────────────────
# Per-row validation
# ────────────────────────────────────────────────────────────────────────

def _chem_canon(c):
    c = (c or "").strip().lower().replace(" ", "_").replace("-", "_")
    c = re.sub(r"[^a-z0-9_]", "", c)
    return CHEM_ALIASES.get(c, c)


def _probe_path(relpath, args):
    """Resolve a probe CSV relpath (species[/class/]<csv>) against the
    runtime probe root.

    Priority: --probe-root arg > $PROBE_ROOT env > PROBE_ROOT default.
    Returns None for Visium HD 3' rows (no probe set).

    PROBE values are the canonical ptegsrv3 relpaths (human/<csv>,
    mouse/HD/<csv>). The root is usually a box-local or remote dir we can
    stat, and the layout there is either:
      (a) nested (ptegsrv3): <root>/<relpath>  — the CSV lives under the
          species[/class/] subfolders;
      (b) flat (SV-76L-002 bundle override): <root>/<basename>  — the CSV
          lives directly in the root, no subpath.
    Resolution order:
      1. If the nested path exists locally -> emit it (ptegsrv3 layout).
      2. Else if the flat-basename path exists locally -> emit it (flat
         bundle override the user pointed --probe-root at).
      3. Else (neither is stat-able, e.g. generator running on
         SV-76L-002 with PROBE_ROOT targeting ptegsrv3) -> emit the nested
         canonical path, trusting the PROBE table as the layout source.
    """
    if not relpath:
        return None
    root = (getattr(args, "probe_root", None) or PROBE_ROOT).strip() or PROBE_ROOT
    root = os.path.expanduser(root)
    nested = os.path.join(root, relpath)
    if os.path.exists(nested):
        return nested
    basename = relpath.rsplit("/", 1)[-1]
    flat = os.path.join(root, basename)
    if os.path.exists(flat):
        return flat
    # Neither present locally — target is canonical layout on a remote box
    # we can't stat from here (default ptegsrv3 target).
    return nested


def _split_punct(s):
    return [t for t in re.split(r"[,; ]+", (s or "").strip()) if t]


def _resolve_fastqs(raw_fastq, detected_avt):
    parts = _split_punct(raw_fastq)
    out = []
    if not parts:
        return out
    for p in parts:
        if os.path.isabs(p):
            out.append(p)
            continue
        if not os.path.sep in p and re.match(r"^[a-z0-9_-]{3,40}$", p):
            if detected_avt:
                base = NGS_BASE_DIR_TEMPLATE.format(avt_run=detected_avt)
                if os.path.isdir(base):
                    found = False
                    for root, dirs, _ in os.walk(base):
                        for d in dirs:
                            if d == p:
                                out.append(os.path.join(root, d))
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        out.append("%s/%s" % (base, p))
                else:
                    out.append(p)
            else:
                out.append(os.path.normpath(p))
            continue
        out.append(os.path.normpath(p))
    return out


def slide_serial_ok(v):
    return bool(re.match(r"^[A-Za-z0-9]{3,10}(?:-[A-Za-z0-9]{3,12})?$", v or ""))


def area_ok(v):
    return bool(re.match(r"^[A-Za-z]\d{0,2}$", v or ""))


def require_path(p, what):
    if p and not os.path.exists(p):
        raise ValueError("%s path does not exist: %s" % (what, p))


def _parse_bool(raw):
    """Return True / False / None (None = column absent / empty)."""
    s = (raw or "").strip().lower()
    if s in ("", "na", "n/a", "none", "null", "nan"):
        return None
    if s in ("0", "false", "no", "f", "off"):
        return False
    return True


def _parse_num(raw):
    """Return a normalised numeric string, or None if absent / not numeric."""
    s = (raw or "").strip()
    if not s:
        return None
    try:
        f = float(s)
    except ValueError:
        return None
    if f == int(f):
        return str(int(f))
    return ("%f" % f).rstrip("0").rstrip(".")


def resolve_sample(raw, args):
    # "Sample name" column is the sample; the tube/library label is only a
    # fallback for workbooks that do not carry a dedicated Sample Name column.
    name = (raw.get("samplename") or raw.get("tubelabel") or "").strip()
    if not name:
        raise ValueError("missing sample name (row %r)" % (dict(raw),))

    species = (raw.get("species") or "").strip().lower()
    if not species:
        species = (args.species or "").lower()
    if species not in ("human", "mouse"):
        raise ValueError("unknown species %r (want mouse|human, or --species)" % raw.get("species"))

    chem_in = (raw.get("chemistry") or raw.get("librarymethod") or "").strip()
    if not chem_in:
        chem_in = args.chemistry or "3p"
    chem = _chem_canon(chem_in)
    if chem in ("nonhd", "nonhd_probeset", "non_probeset", "nonprobeset"):
        raise ValueError(
            "nonHD is no longer supported (removed 2026-09-03: refs never on "
            "this box). This pipeline covers exactly two main applications — "
            "('1') Visium HD WITH probes (hd_probes) and "
            "('2') Visium HD 3' WITHOUT probes (3p). Use one of those.")
    if chem not in CANONICAL_CHEM:
        raise ValueError(
            "unknown chemistry %r — this pipeline covers exactly two main "
            "applications: 'Visium HD with probes' (hd_probes) or "
            "'Visium HD 3'' (3p, no probe set)." % raw.get("chemistry"))

    spec_chem = (species, chem)

    ref = REF.get(spec_chem)
    probe = _probe_path(PROBE.get(spec_chem), args)
    if not ref:
        raise ValueError("no --transcriptome resolved for %s %s — map it in REF" % (species, chem))
    if probe and not os.path.exists(probe):
        # We emit the path as-is regardless. Two legitimate cases:
        #   (a) Generator runs on the execution box (ptegsrv3) and the file
        #       is right there — _probe_path's os.path.exists check should
        #       have caught it already, so this case doesn't reach us.
        #   (b) Generator runs on SV-76L-002 and the target is ptegsrv3's
        #       nested tree — the file exists on ptegsrv3, not locally. This
        #       is the expected cross-box setup (user choice 2026-09-03).
        # Either way, surfacing the path is the correct output — but we flag
        # it on stderr so a human eyeballing the log sees "this one isn't
        # resolvable from where I'm standing."
        sys.stderr.write(
            "note: probe-set path not locally resolvable "
            "(correct if executing on the target box): %s\n" % (probe,))

    fastqs = _resolve_fastqs(raw.get("fastqs") or raw.get("fastqid"), raw.get("_detected_avt"))
    if not fastqs and (args.fastqs or "").strip():
        fastqs = _split_punct(args.fastqs)
        for p in list(fastqs):
            if not os.path.isabs(p):
                fastqs[fastqs.index(p)] = os.path.abspath(p)
    if not fastqs:
        raise ValueError("no fastq path resolved for sample %r (row %r)" % (name, dict(raw)))

    slide = (raw.get("slide") or args.slide or "").strip()
    area = (raw.get("area") or args.area or "").strip()
    unknown_slide = args.unknown_slide

    cytaimage = (raw.get("cytaimage") or args.cytaimage or "").strip()
    image = (raw.get("image") or args.image or "").strip()
    if not (cytaimage or image):
        raise ValueError(
            "spaceranger 4.1.0 requires at least one of --cytaimage / --image / --darkimage "
            "/ --colorizedimage for sample %s (row %r has none)." % (name, dict(raw)))
    require_path(cytaimage, "cytaimage")
    require_path(image, "image")

    # slide / area validation. Per the 4.1.0 binary (`spaceranger count --help`),
    # --slide/--area carry the capture area only when the CytAssist instrument
    # metadata does not; a CytAssist image (--cytaimage) has the slide serial
    # and area embedded in its metadata, in which case --slide/--area are
    # optional (spaceranger takes them from the cytaimage -- and
    # --override-id is the knob that re-applies user values on top of that
    # metadata rather than the other way round). A plain --image (microscopy,
    # no instrument metadata) has NO equivalent, so --slide && --area are
    # required whenever only an image group member without cytaimage is used.
    if unknown_slide:
        # --override-id is not compatible with --unknown-slide per 4.1.0
        # (--unknown-slide is incompatible with --slide, --area, or
        # --slidefile), and spaceranger only accepts --override-id with a
        # known slide serial, so this path is the one that is unambiguous:
        if slide:
            raise ValueError("slide %r set but unknown_slide=%s — spaceranger 4.1.0 "
                             "does not accept --slide with --unknown-slide" % (slide, unknown_slide))
        if area:
            raise ValueError("area %r set but unknown_slide=%s — spaceranger 4.1.0 "
                             "does not accept --area with --unknown-slide" % (area, unknown_slide))
    elif cytaimage:
        # CytAssist path: slide/area optional; if present, must still be well-formed.
        if slide and not slide_serial_ok(slide):
            raise ValueError("slide serial %r not a Visium serial (e.g. 'V10J25-015')" % slide)
        if area and not area_ok(area):
            raise ValueError("area id %r not a Visium area (e.g. 'A1' or 'A')" % area)
    else:
        # microscopy-only path: no instrument metadata → slide && area mandatory.
        if not slide:
            raise ValueError(
                "spaceranger 4.1.0 requires --slide for sample %s (row slide=%r, "
                "area=%r, no cytaimage to carry it). Pass --slide/--area, or pass "
                "--unknown-slide with one of: %s " % (name, slide, area, " | ".join(unknown_slide_choices)))
        if not area:
            raise ValueError(
                "spaceranger 4.1.0 requires --area for sample %s (row slide=%r, "
                "area=%r, no cytaimage to carry it). Pass --slide/--area, or pass "
                "--unknown-slide with one of: %s " % (name, slide, area, " | ".join(unknown_slide_choices)))
        if not slide_serial_ok(slide):
            raise ValueError("slide serial %r not a Visium serial (e.g. 'V10J25-015')" % slide)
        if not area_ok(area):
            raise ValueError("area id %r not a Visium area (e.g. 'A1' or 'A')" % area)

    loupe = (raw.get("loupe") or args.loupe or "").strip()
    require_path(loupe, "loupe-alignment file")

    # S1-S5 nucleus / binning / annotation / image quality (spaceranger 4.1.0).
    # These are optional; emit only when the row provides a value. Absent
    # column (None) -> spaceranger default applies.
    nucleus_seg = _parse_bool(raw.get("nucleus_segmentation"))
    custom_bin = _parse_num(raw.get("custom_bin_size"))
    max_nuc_px = _parse_num(raw.get("max_nucleus_diameter_px"))
    nuc_expum = _parse_num(raw.get("nucleus_expansion_distance_micron"))
    seg_file = (raw.get("custom_segmentation_file") or "").strip()
    ann_model = (raw.get("cell_annotation_model") or "").strip()
    disable_ann = _parse_bool(raw.get("disable_cell_annotation"))
    dapi_idx = _parse_num(raw.get("dapi_index"))
    reorient = _parse_bool(raw.get("reorient_images"))
    img_scale = _parse_num(raw.get("image_scale"))
    override_id = _parse_bool(raw.get("override_id"))
    filter_probes = _parse_bool(raw.get("filter_probes"))

    # 10x constraint: --custom-segmentation-file requires
    # --nucleus-expansion-distance-micron (both must be provided together).
    if seg_file and nuc_expum is None:
        raise ValueError(
            "--custom-segmentation-file is set but "
            "--nucleus-expansion-distance-micron is missing for sample %s "
            "(4.1.0 requires both together)" % name)
    # 10x constraint: --max-nucleus-diameter-px and --custom-segmentation-file
    # are mutually exclusive.
    if seg_file and max_nuc_px is not None:
        raise ValueError(
            "--max-nucleus-diameter-px and --custom-segmentation-file are "
            "mutually exclusive for sample %s" % name)

    # --probe-set rows: --filter-probes defaults to true; honour an explicit
    # false/true from the row.
    if filter_probes is None:
        filter_probes = bool(probe)

    return {
        "name": name,
        "species": species,
        "chem": chem,
        "ref": ref,
        "probe": probe,
        "fastqs": fastqs,
        "slide": slide,
        "area": area,
        "unknown_slide": unknown_slide or "",
        "cytaimage": cytaimage,
        "image": image,
        "loupe": loupe,
        "description": (raw.get("description") or "").strip(),
        "spaceranger_bin": args.spaceranger_bin,
        # S1-S5 optional flags (None / "" -> not emitted, spaceranger default)
        "nucleus_segmentation": nucleus_seg,
        "custom_bin_size": custom_bin,
        "max_nucleus_diameter_px": max_nuc_px,
        "nucleus_expansion_distance_micron": nuc_expum,
        "custom_segmentation_file": seg_file,
        "cell_annotation_model": ann_model,
        "disable_cell_annotation": disable_ann,
        "dapi_index": dapi_idx,
        "reorient_images": reorient,
        "image_scale": img_scale,
        "override_id": override_id,
        "filter_probes": filter_probes if probe else None,
    }


# ────────────────────────────────────────────────────────────────────────
# bash emitter (spaceranger 4.1.0 surface)
# ────────────────────────────────────────────────────────────────────────
# Use a single backslash-continuation string throughout to avoid the
# shell escaping traps (a single '\\ ' inside a Python string literal is
# the same char on either side, and it stays a literal backslash when
# printed to bash).

_BQ = '"'
_BS = "\\"   # bash line continuation


def command_lines(s, create_bam=False):
    """Return the bash lines of the `spaceranger count` invocation for one
    resolved sample, in emission order.

    Shared by the script emitter (build_bash_script) and the
    --print-commands preview so the two can never drift apart: whatever the
    reviewer reads is exactly what the generated .sh executes.
    """
    C = []
    C.append(_CONT('    %s count' % s["spaceranger_bin"]))
    # --id <ID>: spaceranger 4.1.0 help text — "A unique run id and
    # output folder name [a-zA-Z0-9_-]+"  (the sample name is used)
    C.append(_CONT('        --id "%s"' % s["name"]))
    if s["description"]:
        C.append(_CONT('        --description=%s' % LQ(s["description"])))
    C.append(_CONT('        --transcriptome=%s' % LQ(s["ref"])))
    if s["probe"]:
        C.append(_CONT('        --probe-set=%s' % LQ(s["probe"])))
    # --filter-probes: only meaningful on --probe-set rows (defaults true in
    # 4.1.0); emit false only when the row explicitly asked for it.
    if s["probe"] and s["filter_probes"] is False:
        C.append(_CONT('        --filter-probes=false'))
    # image-type / CytAssist flags (4.1.0 requires at least one of the image
    # group; cytaimage/image already validated by resolve_sample)
    if s["cytaimage"]:
        C.append(_CONT('        --cytaimage=%s' % LQ(s["cytaimage"])))
    if s["image"]:
        C.append(_CONT('        --image=%s' % LQ(s["image"])))
    # image-quality controls (only relevant with --darkimage/--fluoro;
    # harmless to emit)
    if s["dapi_index"] is not None:
        C.append(_CONT('        --dapi-index=%s' % s["dapi_index"]))
    if s["reorient_images"] is not None:
        C.append(_CONT('        --reorient-images=%s' %
                       ("true" if s["reorient_images"] else "false")))
    if s["image_scale"] is not None:
        C.append(_CONT('        --image-scale=%s' % s["image_scale"]))
    # slide / area, or unknown-slide fallback. 4.1.0: --slide/--area are
    # optional when a CytAssist image is present (its metadata carries the
    # slide serial and capture area); they stay mandatory for a plain
    # --image with no instrument metadata.
    if s["unknown_slide"]:
        C.append(_CONT('        --unknown-slide %s' % s["unknown_slide"]))
    if s["slide"]:
        C.append(_CONT('        --slide=%s' % LQ(s["slide"])))
        if s["area"]:
            C.append(_CONT('        --area=%s' % LQ(s["area"])))
    # --override-id: 4.1.0 says use it with --slide/--area or with
    # --unknown-slide (it overrides the value read from the CytAssist image
    # metadata). It is only emitted when a slide/unknown-slide context exists.
    if s["override_id"] is True and (s["unknown_slide"] or s["slide"]):
        C.append(_CONT('        --override-id'))
    # nucleus / segmentation family. --nucleus-segmentation is a VALUED flag
    # in 4.1.0 (true|false), not a bare flag.
    if s["nucleus_segmentation"] is not None:
        C.append(_CONT('        --nucleus-segmentation=%s' %
                       ("true" if s["nucleus_segmentation"] else "false")))
    if s["custom_bin_size"] is not None:
        C.append(_CONT('        --custom-bin-size=%s' % s["custom_bin_size"]))
    if s["max_nucleus_diameter_px"] is not None:
        C.append(_CONT('        --max-nucleus-diameter-px=%s' %
                       s["max_nucleus_diameter_px"]))
    if s["nucleus_expansion_distance_micron"] is not None:
        C.append(_CONT('        --nucleus-expansion-distance-micron=%s' %
                       s["nucleus_expansion_distance_micron"]))
    if s["custom_segmentation_file"]:
        C.append(_CONT('        --custom-segmentation-file=%s' %
                       LQ(s["custom_segmentation_file"])))
    # cell annotation family
    if s["cell_annotation_model"]:
        C.append(_CONT('        --cell-annotation-model=%s' %
                       LQ(s["cell_annotation_model"])))
    if s["disable_cell_annotation"] is True:
        C.append(_CONT('        --disable-cell-annotation'))
    # loupe alignment
    if s["loupe"]:
        C.append(_CONT('        --loupe-alignment=%s' % LQ(s["loupe"])))
    # --fastqs is REQUIRED by spaceranger 4.1.0. Emit one --fastqs per
    # resolved path (repeatable).
    for fq in s["fastqs"]:
        C.append(_CONT('        --fastqs=%s' % LQ(fq)))
    C.append(_CONT('        --create-bam=%s' % ("true" if create_bam else "false")))
    C.append('        --output-dir="$OUTPUT_DIR/${sn}"')
    return C


def print_commands(resolved, max_jobs, output_dir, create_bam=False):
    """Human-readable review of every `spaceranger count` invocation, built
    from command_lines() so it matches the generated script exactly. Prints to
    stdout and returns 0. No file is written and nothing is launched."""
    W = 76
    bar = "=" * W
    print("")
    print(bar)
    print("SpaceRanger  %-9s — %d sample(s) — up to %d concurrent" % (
        "count", len(resolved), max_jobs))
    print("Reference and transcriptome paths resolved on this machine.")
    print(bar)
    for i, s in enumerate(resolved, 1):
        print("")
        print(" [%02d] %s    (%s %s, %d fastq)" % (
            i, s["name"], s["species"], s["chem"], len(s["fastqs"])))
        for line in command_lines(s, create_bam):
            t = line.rstrip()
            if "--output-dir=" in t:
                t = t[:t.index("--output-dir=")] + "--output-dir=" + \
                    os.path.join(output_dir, s["name"])
            print("     " + t.strip())
    print("")
    print(bar)
    print("Preview only — nothing was written and nothing was launched.")
    print("Generate the script:  python3 %s <workbook> ..." %
        os.path.basename(__file__))
    print("Run it end to end:    bash run_visium.sh <workbook> <output_dir> ...")


def build_bash_script(samples, max_jobs, output_dir, create_bam=False):
    L = []
    L.append("#!/usr/bin/env bash")
    L.append("# auto-generated by generate_spaceRanger_visium.py (spaceranger 4.1.0)")
    L.append("set -uo pipefail")
    L.append("IFS=$'\\n\\t'")
    L.append("")
    L.append("OUTPUT_DIR=%s" % LQ(output_dir))
    L.append("LOG_DIR=\"$OUTPUT_DIR/logs\"")
    L.append('mkdir -p "$LOG_DIR" "$OUTPUT_DIR"')
    L.append("MAX_JOBS=%d" % max_jobs)
    L.append("FAILURES=0")
    L.append("FAIL_LOG=\"$OUTPUT_DIR/.failures.txt\"")
    L.append(': > "$FAIL_LOG"')
    L.append("")
    L.append("# Launch a sample in the background, capped at MAX_JOBS concurrent.")
    L.append("launch() {")
    L.append("    local fn=$1")
    L.append('    local running; running=$(jobs -rp | wc -l)')
    L.append('    if (( running >= MAX_JOBS )); then wait -n; fi')
    L.append('    "$fn" &')
    L.append("}")

    for i, s in enumerate(samples, 1):
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", s["name"]).replace('"', '_')
        L.append("")
        L.append("# %02d  %s  (%s %s, %d fastq)" % (
            i, s["name"], s["species"], s["chem"], len(s["fastqs"])))
        L.append("run_%02d() {" % i)
        L.append('    local sn=%s' % LQ(s["name"]))
        L.append('    local log="$LOG_DIR/${sn}.log"')
        L.append('    : > "$log"')
        L.append('    echo ">> $sn" | tee -a "$log"')
        L.append('    cd "$OUTPUT_DIR"')
        # spaceranger count invocation (single source of truth)
        L.extend(command_lines(s, create_bam))
        L.append("    local rc=$?")
        L.append('    echo "run rc=$rc" | tee -a "$log"')
        L.append("    if (( rc != 0 )); then")
        L.append('        echo "FAIL $sn rc=$rc" >> "$FAIL_LOG"')
        L.append("        FAILURES=$((FAILURES + 1))")
        L.append("        return 1")
        L.append("    fi")
        L.append("    # Remove SPATIAL_RNA_COUNTER_CS (martian runtime CS stages, not a deliverable)")
        L.append('    rm -r "$OUTPUT_DIR/${sn}/SPATIAL_RNA_COUNTER_CS" 2>/dev/null || true')
        L.append("    return 0")
        L.append("}")
        L.append("launch run_%02d" % i)

    L.append("")
    L.append("wait")
    L.append("")
    L.append("# ------------------------------------ POST-PROCESSING ------------------------------------ #")
    L.append("")
    L.append('SUMMARY_DIR="$OUTPUT_DIR/Summaries"')
    L.append('mkdir -p "$SUMMARY_DIR/web_summaries"')
    L.append("")
    L.append("# Collect per-sample web summary HTMLs into one folder for easy browsing")
    for s in samples:
        L.append(
            'cp -n "$OUTPUT_DIR/%s/outs/web_summary.html" '
            '"$SUMMARY_DIR/web_summaries/%s_WebSummary.html" 2>/dev/null || true'
            % (s["name"], s["name"]))
    L.append("")
    L.append("# Pull each sample outs/metrics_summary.csv into Summaries/ and sanitize it:")
    L.append("#   - drop thousands-separator commas inside quoted numeric fields")
    L.append("#   - strip a percent sign that follows a number (e.g. 75% -> 75)")
    for s in samples:
        nm = s["name"]
        L.append('if [ -f "$OUTPUT_DIR/%s/outs/metrics_summary.csv" ]; then' % nm)
        L.append('    mf="$SUMMARY_DIR/%s_metrics_summary.csv"' % nm)
        L.append('    cp "$OUTPUT_DIR/%s/outs/metrics_summary.csv" "$mf"' % nm)
        # (awk/sed sanitize; sample-independent, operates on the $mf set above)
        L.append(
            r'''    awk -F'"' 'BEGIN {OFS = "\""} { for (i=2; i<=NF; i+=2) gsub(",", "", $i) } 1' "$mf" > "$mf.tmp" && mv "$mf.tmp" "$mf"''')
        L.append(
            r'''    sed -E 's/([0-9]+)%/\1/g' "$mf" > "$mf.tmp" && mv "$mf.tmp" "$mf"''')
        L.append('else')
        L.append('    echo "WARNING: %s metrics_summary.csv not found" >&2' % nm)
        L.append('fi')
    L.append("")
    L.append("# Merge sanitized per-sample metrics into one table with a Samplename column")
    L.append("# (10X metrics_summary.csv is exactly a header + one data row)")
    L.append('merged_csv="$SUMMARY_DIR/all_metrics_summary.csv"')
    L.append(': > "$merged_csv"')
    L.append('header_set=false')
    for s in samples:
        nm = s["name"]
        L.append('mf="$SUMMARY_DIR/%s_metrics_summary.csv"' % nm)
        L.append('if [ -f "$mf" ]; then')
        L.append('    if [ "$header_set" = false ]; then')
        L.append('        printf "Samplename," > "$merged_csv"')
        L.append('        sed -n 1p "$mf" >> "$merged_csv"')
        L.append('        header_set=true')
        L.append('    fi')
        L.append("    printf '%s,' >> \"$merged_csv\"" % nm)
        L.append('    sed -n 2p "$mf" >> "$merged_csv"')
        L.append('fi')
    L.append("")
    L.append("# Deploy the finished run to the central sequence archive (fixed target)")
    L.append('SEQ_DIR="/mnt/PTEGraw/AA/NGSruns/AVT/AV240401/ZZ_10XG/spaceranger"')
    L.append('deploy_ok=true')
    L.append('if [ -d "$SEQ_DIR" ]; then')
    L.append('    echo ">> deploying $OUTPUT_DIR -> $SEQ_DIR"')
    L.append('    rsync -a "$OUTPUT_DIR" "$SEQ_DIR/" || deploy_ok=false')
    L.append('else')
    L.append('    echo "ERROR: SEQ_DIR $SEQ_DIR is missing; cannot deploy" >&2')
    L.append('    deploy_ok=false')
    L.append('fi')
    L.append("")
    L.append('count=$(ls -d "$OUTPUT_DIR"/*/outs 2>/dev/null | wc -l)')
    L.append('echo "done: $count sample output dirs under $OUTPUT_DIR"')
    L.append("")
    L.append("if (( FAILURES > 0 )); then")
    L.append('    echo "some runs failed — see $FAIL_LOG" >&2')
    L.append("    exit 1")
    L.append("fi")
    L.append('if [ "$deploy_ok" != true ]; then')
    L.append('    echo "ERROR: deploy to $SEQ_DIR failed" >&2')
    L.append("    exit 1")
    L.append("fi")
    L.append("exit 0")
    return "\n".join(L) + "\n"


def LQ(s):
    """Bash single-quote a string (safe for any text)."""
    s = str(s)
    return "'" + s.replace("'", "'\\''") + "'"


def _CONT(s):
    """Append a bash line-continuation backslash to a command fragment."""
    return s + " " + _BS

# ────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Generate a spaceranger 4.1.0 `count` runner for Visium / "
                    "SpaceRanger, from a VL or GECF sample workbook.")
    ap.add_argument("sheets", nargs="+",
                    help="one or more .xlsx workbooks to parse")
    ap.add_argument("--species", default=None,
                    help="default species (human|mouse) for rows without one. "
                         "Required for the GECF submission form.")
    ap.add_argument("--chemistry", default="3p",
                    help="default application for rows without one. Two main "
                         "applications: 'Visium HD with probes' (hd_probes) "
                         "or 'Visium HD 3'' (HD_3P, no probe set). "
                         "Default: 3p (Visium HD 3').")
    ap.add_argument("--probe-root", default=os.environ.get("PROBE_ROOT", PROBE_ROOT),
                    help="directory containing the Visium HD probe CSV (the "
                         "'Visium HD WITH probes' application only; the 3' "
                         "application emits no probe set). Default: %s, or "
                         "$PROBE_ROOT. PROBE carries the species[/class/] "
                         "subpath (ptegsrv3 canonical: human/<csv>, "
                         "mouse/HD/<csv>) which is joined onto this dir; the "
                         "flat SV-76L-002 bundle layout also resolves via "
                         "the _probe_path flat-join fallback." % PROBE_ROOT)
    ap.add_argument("--slide", default="",
                    help="default slide serial (e.g. 'V10J25-015') for rows "
                         "without one. Optional when the row also carries a "
                         "CytAssist image (--cytaimage, its metadata carries "
                         "the slide serial); required for a plain --image row "
                         "unless --unknown-slide is passed. Not compatible "
                         "with --unknown-slide.")
    ap.add_argument("--area", default="",
                    help="default area id (e.g. 'A1') for rows without one. "
                         "Optional with --cytaimage; required for a plain "
                         "--image row unless --unknown-slide is passed. Not "
                         "compatible with --unknown-slide.")
    ap.add_argument("--cytaimage", default="",
                    help="default --cytaimage path for rows without one (the "
                         "CytAssist/cyte image). spaceranger 4.1.0 requires at "
                         "least one of --cytaimage/--image/--darkimage/"
                         "--colorizedimage; the GECF submission workbook has "
                         "none of these columns, so pass it here when they are "
                         "shared across the run.")
    ap.add_argument("--image", default="",
                    help="default --image path for rows without one (OME-TIFF "
                         "spatial image). Use if the run has a plain spatial "
                         "image instead of a CytAssist image.")
    ap.add_argument("--unknown-slide", default=None,
                    choices=unknown_slide_choices,
                    help="use --unknown-slide <one of these> instead of "
                         "--slide/--area when the serial/area were unknown.")
    ap.add_argument("--loupe", default="",
                    help="global default --loupe-alignment path; per-row value wins.")
    ap.add_argument("--fastqs", default="",
                    help="global default FASTQ dir(s) (space/comma separated) for rows "
                         "without a resolvable fastq path. Overrides the per-row fastq "
                         "ID when that cannot be resolved on this machine.")
    ap.add_argument("--max-jobs", type=int, default=4,
                    help="max concurrent spaceranger runs (default 4)")
    ap.add_argument("--create-bam", default="false",
                    choices=["true", "false"],
                    help="--create-bam value (default false)")
    ap.add_argument("--spaceranger-bin", default=os.environ.get(
        "SPACERANGER_BIN", "spaceranger"),
        help="spaceranger binary path (default: spaceranger on PATH, "
             "or $SPACERANGER_BIN)")
    ap.add_argument("--output-dir", default=".",
                    help="output dir for generated artifacts (default: cwd)")
    ap.add_argument("--output-name", default=None,
                    help="basename of the generated shell script (without .sh)")
    ap.add_argument("--print", action="store_true",
                    help="print the generated script to stdout instead of "
                         "writing a file")
    ap.add_argument("--print-commands", action="store_true",
                    help="parse the workbook, validate the samples, and print "
                         "the exact spaceranger count command for each sample "
                         "for quick review; do not write or launch anything")
    args = ap.parse_args(argv)

    bin_path = shutil.which(args.spaceranger_bin)
    if not bin_path:
        sys.stderr.write(
            "WARN: %r not on PATH; the generated script will fail until it is.\n"
            % args.spaceranger_bin)

    resolved = []
    for sheet in args.sheets:
        rows = parse_xlsx(sheet)
        if not rows:
            print("FATAL: no data rows found in %s" % sheet, file=sys.stderr)
            sys.exit(2)
        for raw in rows:
            try:
                r = resolve_sample(raw, args)
                resolved.append(r)
            except ValueError as e:
                print("SKIPPED %s : %s" % (sheet, e), file=sys.stderr)
    if not resolved:
        print("FATAL: no usable samples after validation — at least one row "
              "must provide an image group (cytaimage/image), with --slide AND "
              "--area (or --unknown-slide) unless a CytAssist image carries the "
              "slide/area in its metadata.", file=sys.stderr)
        sys.exit(2)

    output_dir = args.output_dir or "."
    create_bam = (args.create_bam == "true")

    if args.print_commands:
        return print_commands(resolved, args.max_jobs, output_dir, create_bam)

    script = build_bash_script(resolved, args.max_jobs, output_dir, create_bam)

    if args.print:
        sys.stdout.write(script)
        return 0

    base = args.output_name
    if not base:
        base = "spaceranger_visium_" + os.path.splitext(
            os.path.basename(args.sheets[0]))[0]
    path = os.path.join(output_dir, base + ".sh")
    os.makedirs(output_dir, exist_ok=True)
    with open(path, "w") as f:
        f.write(script)
    os.chmod(path, 0o755)

    print("Wrote %d sample blocks -> %s" % (len(resolved), path))
    print("")
    for r in resolved:
        print("  %-32s  %-7s  %-3s  %d fastq" % (
            r["name"], r["species"], r["chem"], len(r["fastqs"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
