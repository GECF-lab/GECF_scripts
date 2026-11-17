#!/usr/bin/env python3
"""custom_gene_report.py -- HTML report for custom-gene (transgene) expression.

Reads cellranger count `outs/` dir(s) (filtered_feature_bc_matrix) for a list
of custom genes and writes a self-contained HTML report with:

  * Stats table (one row per gene):
      - detected cells (UMI per cell >= --min-umi) and % of all cells
      - cells at >=1 UMI and cells at >=3 UMI (fixed sensitivity columns;
        the one equal to --min-umi is omitted)
      - total UMIs on the transgene row
      - read count of the gene (primary alignments from the BAM, called
        cells only)
      - mean and median UMI among positive (>=1 UMI) cells
      - endogenous cell/UMI context
  * % positive per sample table (project mode).
  * Plots (Plotly figures embedded inline, no CDN; offline in the lab):
      - UMI histogram of positive cells (>= 1 UMI), per sample
      - % positive cells per sample (bars) + mean-UMI dot plot (project mode)
      - positive-vs-negative total-UMIs-per-cell depth boxplot, per sample
      - UMAP coloured by the gene + binary positive/negative panel
      - coverage track along the transgene contig
        (needs BAM + samtools + coordinate-bearing feature rows)

Missing optional inputs (BAM, samtools, UMAP coordinates)
degrade to "n/a" panels. Nothing is required beyond the filtered matrix; Plotly
is the only hard Python dependency (BAM/samtools are best-effort only).
numpy/scipy are optional accelerators: when importable the matrix scan runs
at C speed, otherwise the script silently falls back to a pure-stdlib line
scan (same result, slower).

The transgene row is the feature whose id equals the gene name (the appended
contig); endogenous rows are other rows carrying the same name (e.g. Ensembl).
Status: EXPRESSED = >=1 detected cell (UMI >= --min-umi, default 3);
SILENT = transgene row present, 0 detected cells;
NO-ROW = no transgene row in the feature table at all.

Usage:
  python3 custom_gene_report.py --sample-dir <sample_dir> --genes genes.txt [-o OUT_FOLDER]
  python3 custom_gene_report.py --project-dir <project_dir> --genes genes.txt [-o OUT_FOLDER]
  python3 custom_gene_report.py --help

Python environment: no conda activation is needed -- if the interpreter you
launched the script with does not provide numpy/scipy/plotly, the script
re-runs itself (transparently, same arguments) under the lab's conda env,
brbseq first, then basecaller (prod paths, then draft-box paths).
CGR_NO_ENV_SWITCH=1 disables the switch.
"""
import argparse
import errno
import gzip
import html
import math
import os
import shutil
import statistics
import subprocess
import sys
import time


# --------------------------------------------------------------------------- #
# conda-env auto-resolution (no activation needed)
#
# The lab runs this on the prod server with plain `python` in a shell where no
# conda env is active. numpy/scipy (fast matrix scan) and plotly all come from
# the lab's brbseq / basecaller envs, so if the current interpreter lacks them
# we re-exec the script under that env's python (brbseq preferred). The new
# process is identical (same script, same args, same cwd) -- nothing to
# activate, nothing to source.
#
#   CGR_NO_ENV_SWITCH=1        disable the auto switch entirely
# --------------------------------------------------------------------------- #
_CANDIDATE_ENVS = [
    # prod server -- NOTE: the env dir sits in a literal bracket-named
    # subdirectory on this host (yes, the brackets are part of the path);
    # python strings are not shell-globbed, so quote nothing, escape nothing.
    "/home/gecf/softwares/anaconda/[/home/gecf/softwares/anaconda]/envs/brbseq",
    "/home/gecf/softwares/anaconda/[/home/gecf/softwares/anaconda]/envs/basecaller",
    # prod server -- plain layout, kept in case a host does not have the
    # bracket dir (missing dirs are skipped in one isfile check)
    "/home/gecf/softwares/anaconda/envs/brbseq",
    "/home/gecf/softwares/anaconda/envs/basecaller",
    # draft box
    "/home/hermes/.conda/envs/brbseq",
    "/home/hermes/miniconda3/envs/basecaller",
]

_SAMTOOLS_RESOLVED = None


def samtools_cmd():
    """Resolved `samtools` executable: PATH first, then the known lab envs.

    Re-exec into a conda env's python does NOT put that env's bin/ on PATH,
    so a plain `shutil.which("samtools")` (or `subprocess.run(["samtools"])`)
    silently fails under env pythons even when the env ships samtools.
    Without this, every BAM-derived column degrades to 0 / n-a while the
    matrix columns still populate -- confusing, and it is exactly what made
    the 'reads (called cells)' column read zero.
    """
    global _SAMTOOLS_RESOLVED
    if _SAMTOOLS_RESOLVED is not None:
        return _SAMTOOLS_RESOLVED
    s = shutil.which("samtools")
    if s is None:
        for env in _CANDIDATE_ENVS:
            p = os.path.join(env, "bin", "samtools")
            try:
                if os.path.isfile(p) and os.access(p, os.X_OK):
                    s = p
                    break
            except OSError:
                continue
    _SAMTOOLS_RESOLVED = s
    return s


def _interpreter_ok(py):
    r = subprocess.run(
        [py, "-c", "import plotly, numpy, scipy"],
        timeout=60, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def _maybe_use_conda_env():
    if os.environ.get("CGR_NO_ENV_SWITCH"):
        return
    if len(sys.argv) < 2 or sys.argv[0] in ("", "-c"):
        return  # not running as this script file
    cur = os.path.abspath(sys.executable)
    # Already running under one of the target envs (user activated it), or
    # the current interpreter already has what we need: keep it.
    if any(cur.startswith(os.path.abspath(e) + os.sep)
           for e in _CANDIDATE_ENVS):
        return
    if cur:
        try:
            if _interpreter_ok(cur):
                return
        except Exception:
            pass
    for py in [os.path.join(e, "bin", "python") for e in _CANDIDATE_ENVS]:
        if not os.path.isfile(py) or not os.access(py, os.X_OK):
            continue
        try:
            if not _interpreter_ok(py):
                continue
        except Exception:
            continue
        env_name = os.path.basename(os.path.dirname(os.path.dirname(py)))
        print("[custom_gene_report] %s lacks numpy/scipy/plotly "
              "-- re-running under %s (%s)" % (cur or "current python",
                                              env_name, py),
              file=sys.stderr)
        os.environ["CGR_NO_ENV_SWITCH"] = "1"
        os.execv(py, [py, os.path.abspath(__file__)] + sys.argv[1:])
        return
    # No usable lab env found: keep the current interpreter (the pure-stdlib
    # fallbacks keep the script functional either way).


_maybe_use_conda_env()

VERSION = "1.1.0"


# --------------------------------------------------------------------------- #
# file helpers
# --------------------------------------------------------------------------- #
def opener(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") \
        if path.endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def read_lines(path):
    with opener(path) as f:
        for line in f:
            if line.strip():
                yield line.rstrip("\n")


def find_file(base_dir, names):
    for n in names:
        p = os.path.join(base_dir, n)
        if os.path.isfile(p):
            return p
    return None


FEATURE_NAMES = ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]
MATRIX_NAMES = ["matrix.mtx.gz", "matrix.mtx"]
BARCODE_NAMES = ["barcodes.tsv.gz", "barcodes.tsv"]


def discover_sample(out_dir):
    """Locate filtered_feature_bc_matrix core files.

    Accepts either the cellranger `outs` dir, the sample dir that contains an
    `outs/` subfolder, or the matrix dir itself."""
    candidates = [
        out_dir,
        os.path.join(out_dir, "filtered_feature_bc_matrix"),
        os.path.join(out_dir, "outs", "filtered_feature_bc_matrix"),
    ]
    for fbm in candidates:
        features = find_file(fbm, FEATURE_NAMES)
        barcodes = find_file(fbm, BARCODE_NAMES)
        matrix = find_file(fbm, MATRIX_NAMES)
        if features and barcodes and matrix:
            return {"fbm": fbm, "features": features, "barcodes": barcodes,
                    "matrix": matrix}
    raise SystemExit(f"ERROR: no filtered_feature_bc_matrix found under '{out_dir}'")


def read_features(path):
    """-> list[(idx1based, feature_id, name)]."""
    out = []
    for i, line in enumerate(read_lines(path), start=1):
        parts = line.split("\t")
        fid = parts[0]
        name = parts[1] if len(parts) > 1 and parts[1].strip() else fid
        out.append((i, fid, name))
    return out


def iter_matrix_entries(path):
    """Yield (row1, col1, val) from a MatrixMarket coordinate matrix.

    First non-comment line (dimension line) is skipped.
    """
    with opener(path) as f:
        seen_dims = False
        for line in f:
            if not line or line.startswith("%"):
                continue
            if not seen_dims:
                seen_dims = True
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    yield int(parts[0]), int(parts[1]), int(parts[2])
                except ValueError:
                    continue


def read_genes(genes_arg):
    if os.path.isfile(genes_arg):
        genes = []
        for line in read_lines(genes_arg):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.append(line.split()[0])
        return genes
    return [g for g in genes_arg.replace(",", " ").split() if g]


# --------------------------------------------------------------------------- #
# optional data sources (best-effort, one scan each)
# --------------------------------------------------------------------------- #
def bam_contigs(bam):
    """contig name -> length from the BAM @SQ header; None if not usable.

    Uses `samtools view -H` (any samtools in PATH, e.g. from a conda env)."""
    if not bam or not os.path.isfile(bam):
        return None
    st = samtools_cmd()
    if st is None:
        return None
    try:
        out = subprocess.check_output(
            [st, "view", "-H", bam],
            stderr=subprocess.DEVNULL, timeout=120, text=True)
    except Exception:
        return None
    contigs = {}
    for line in out.splitlines():
        if not line.startswith("@SQ"):
            continue
        name = length = None
        for tok in line.split("\t"):
            if tok.startswith("SN:"):
                name = tok[3:]
            elif tok.startswith("LN:"):
                try:
                    length = int(tok[3:])
                except ValueError:
                    pass
        if name and length is not None:
            contigs[name] = length
    return contigs or None


def bam_gene_depth(bam, contig):
    """Base coverage of one contig: {pos1: depth}, exactly as reported by
    `samtools depth -r <contig>` (no CIGAR re-interpretation on our side).

    None when the BAM, samtools, or the scan is unavailable."""
    if not bam or not os.path.isfile(bam):
        return None
    st = samtools_cmd()
    if st is None:
        return None
    try:
        out = subprocess.check_output(
            [st, "depth", "-r", contig, bam],
            stderr=subprocess.DEVNULL, timeout=900, text=True)
    except Exception:
        return None
    depth = {}
    for line in out.splitlines():
        p = line.split()
        if len(p) < 3:
            continue
        try:
            depth[int(p[1])] = int(p[2])
        except ValueError:
            continue
    return depth or None


def bam_gene_reads(bam, contig, genes, called):
    """Primary-alignment read count per gene: aligned on `contig` (RNAME),
    not secondary/supplementary, GN:<gene> tagged, from a called cell
    (CB tag matching a matrix barcode, with or without a -N duplicate
    suffix).

    One `samtools view -r <contig>` pass. Returns {gene: n}; {} when the
    scan could not run."""
    if not bam or not os.path.isfile(bam):
        return {}
    st = samtools_cmd()
    if st is None:
        return {}
    try:
        out = subprocess.check_output(
            [st, "view", bam, contig],   # region is a positional arg for `view`
            stderr=subprocess.DEVNULL, timeout=900, text=True)
    except Exception:
        return {}
    geneset = set(genes)
    reads = {g: 0 for g in genes}
    for line in out.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        if len(f) < 11 or f[2] != contig:   # RNAME guard: some samtools builds
            continue                        # can bleed other contigs into -r
        try:
            flags = int(f[1])
        except ValueError:
            continue
        if flags & 0x900:              # secondary / supplementary
            continue
        gn = None
        cb = None
        # tag area (after QUAL): some samtools builds emit tags space-
        # separated, so split on any whitespace
        for tok in " ".join(f[11:]).split():
            parts = tok.split(":", 2)
            if len(parts) != 3:
                continue
            if parts[0] == "GN":
                gn = parts[2]
            elif parts[0] == "CB":
                cb = parts[2]
        if gn in geneset and cb is not None \
                and (cb in called or cb.split("-", 1)[0] in called):
            reads[gn] += 1
    return reads


def sample_bam_summary(data, genes):
    """Per-sample BAM summary: transgene contig, base coverage, and reads
    per gene from called cells.

    Two lightweight passes, each opening the BAM once: `samtools depth -r`
    for the coverage track (numbers exactly as `depth` reports them) and
    `samtools view -r` for the primary-alignment read count.

    Returns a dict {gene, contig, length, depth, reads, contigs_n, no_contig}
    or None when no usable BAM/header is available for this sample.
    """
    bam = data.get("bam")
    if not bam:
        return None
    contigs = bam_contigs(bam)
    if contigs is None:
        return None
    called = set()
    for bc in data.get("barcodes") or []:      # called cells
        called.add(bc)
        called.add(bc.split("-", 1)[0])        # CB tag without a -N dedup suffix
    for g in genes:
        m = match_contig([g], contigs)
        if m is None:
            continue
        gene, cname = m
        length = contigs.get(cname) or 0
        depth = bam_gene_depth(bam, cname)
        reads = bam_gene_reads(bam, cname, genes, called)
        if depth is None and not reads:
            continue
        return {"gene": gene, "contig": cname, "length": length,
                "depth": depth, "reads": reads,
                "contigs_n": len(contigs), "no_contig": False}
    return {"gene": genes[0], "contig": None, "length": 0, "depth": None,
            "reads": None, "contigs_n": len(contigs),
            "no_contig": True}


# --------------------------------------------------------------------------- #
# per-gene stats for one sample
# --------------------------------------------------------------------------- #
def compute_sample(sample_dir, bam, genes, threshold):
    fbm = discover_sample(sample_dir)
    barcodes = list(read_lines(fbm["barcodes"]))
    cells_total = len(barcodes)
    wanted = set(genes)

    tg_rows = {}    # gene -> [feature idx, ...]   (transgene rows: id == name)
    endo_rows = {}  # gene -> [feature idx, ...]   (id != name, context rows)
    for idx, fid, name in read_features(fbm["features"]):
        if name not in wanted:
            continue
        if fid == name:
            tg_rows.setdefault(name, []).append(idx)
        else:
            endo_rows.setdefault(name, []).append(idx)

    stats = {g: {"n_pos": 0, "ge1": 0, "ge3": 0, "det": 0,
                 "umi": 0, "reads": None, "pos_vals": [],
                 "endo_cells": 0, "endo_umi": 0,
                 "row_present": g in tg_rows,
                 "mean": 0.0, "median": 0.0,
                 "status": "SILENT", "cells_total": cells_total}
             for g in genes}

    # per-barcode tallies (for plots): cell total UMI + transgene UMI per gene
    bc_total = {}          # barcode -> total UMI (all genes)
    bc_gene = {}           # barcode -> {gene -> UMI}
    try:
        # fast path: scipy C parser + numpy sums (C speed even on big runs)
        _scan_matrix_numpy(fbm["matrix"], barcodes, tg_rows, endo_rows,
                           threshold, stats, bc_total, bc_gene)
    except Exception as exc:
        # no scipy, or a matrix mmread can't handle -> stdlib line scan
        try:
            sys.stderr.write(
                "note: fast matrix scan unavailable (%s); "
                "falling back to the line-by-line scan\n"
                % exc.__class__.__name__)
        except Exception:
            pass
        _scan_matrix_python(fbm["matrix"], barcodes, tg_rows, endo_rows,
                            threshold, stats, bc_total, bc_gene)

    for g in genes:
        s = stats[g]
        if s["row_present"]:
            s["status"] = "EXPRESSED" if s["det"] > 0 else "SILENT"
        else:
            s["status"] = "NO-ROW"
        if s["n_pos"]:
            s["mean"] = s["umi"] / s["n_pos"]
            s["median"] = statistics.median(s["pos_vals"])

    return {"sample_dir": sample_dir, "bam": bam, "cells_total": cells_total,
            "barcodes": barcodes, "bc_total": bc_total,
            "bc_gene": bc_gene, "stats": stats}


def _scan_matrix_numpy(matrix_path, barcodes, tg_rows, endo_rows,
                       threshold, stats, bc_total, bc_gene):
    """Fast matrix scan: scipy C parser (+ numpy) instead of a Python loop.

    The whole sparse matrix is read once by scipy's C code, then:
      * per-cell total UMI      -> one C column sum over all nonzeros
      * each transgene row      -> a single row slice (a handful of rows)
    so the Python loop touches only (a) barcodes of cells with any UMI and
    (b) the positive hits on the wanted transgene rows -- not every matrix
    entry. Raises (caller falls back to _scan_matrix_python) if scipy is
    missing or the matrix dimensions disagree with barcodes.tsv.
    """
    import numpy as np
    from scipy.io import mmread

    X = mmread(str(matrix_path), spmatrix=True)
    ncells = len(barcodes)
    if X.shape[1] != ncells:
        raise ValueError("matrix has %d cells, barcodes.tsv has %d"
                         % (X.shape[1], ncells))
    csr = X if X.format == "csr" else X.tocsr()

    # total UMI per cell: one C pass over every nonzero in the matrix
    colsum = np.asarray(csr.sum(axis=0)).ravel()
    for i in np.flatnonzero(colsum).tolist():
        bc_total[barcodes[i]] = int(colsum[i])

    # transgene rows: the only rows whose per-cell values we need
    for g, rows in tg_rows.items():
        s = stats[g]
        vals = []
        for r1 in rows:
            row = csr.getrow(r1 - 1)
            data = row.data
            vals.extend(data.tolist())
            for ci, v in zip(row.indices.tolist(), data.tolist()):
                bc = barcodes[int(ci)]
                cd = bc_gene.get(bc)
                if cd is None:
                    bc_gene[bc] = {g: int(v)}
                else:
                    cd[g] = cd.get(g, 0) + int(v)
        if vals:
            a = np.asarray(vals, dtype=np.int64)
            s["n_pos"] = int(a.size)
            s["ge1"] = int((a >= 1).sum())
            s["ge3"] = int((a >= 3).sum())
            s["det"] = int((a >= threshold).sum())
            s["umi"] = int(a.sum())
            s["pos_vals"] = a.tolist()

    # endogenous context rows: total cells/UMIs touching the name
    for g, rows in endo_rows.items():
        s = stats[g]
        n = 0
        tot = 0
        for r1 in rows:
            d = csr.getrow(r1 - 1).data
            n += int(d.size)
            tot += int(d.sum())
        s["endo_cells"] = n
        s["endo_umi"] = tot


def _scan_matrix_python(matrix_path, barcodes, tg_rows, endo_rows,
                        threshold, stats, bc_total, bc_gene):
    """Line-by-line pure-stdlib scan (original algorithm).

    Kept as the fallback for environments without scipy/numpy, or for a
    matrix scipy.io.mmread refuses. Slow on large runs, always works.
    """
    tg_idx = {}
    for g, rows in tg_rows.items():
        for r1 in rows:
            tg_idx[r1] = g
    endo_idx = {}
    for g, rows in endo_rows.items():
        for r1 in rows:
            endo_idx[r1] = g
    ncells = len(barcodes)

    for r, c, v in iter_matrix_entries(matrix_path):
        bc = barcodes[c - 1] if 1 <= c <= ncells else None
        if bc is not None:
            bc_total[bc] = bc_total.get(bc, 0) + v
            if r in tg_idx:
                g = tg_idx[r]
                bcd = bc_gene.setdefault(bc, {})
                bcd[g] = bcd.get(g, 0) + v
        g = tg_idx.get(r)
        if g is None:
            g = endo_idx.get(r)
            if g is None:
                continue
            s = stats[g]
            s["endo_cells"] += 1
            s["endo_umi"] += v
            continue
        s = stats[g]
        s["n_pos"] += 1
        if v >= 1:
            s["ge1"] += 1
        if v >= 3:
            s["ge3"] += 1
        if v >= threshold:
            s["det"] += 1
        s["umi"] += v
        s["pos_vals"].append(v)


# --------------------------------------------------------------------------- #
# discovery for project mode
# --------------------------------------------------------------------------- #
BAM_NAMES = ["possorted_genome_bam.bam", "possorted_bam.bam",
             "gex_possorted_bam.bam"]


def detect_bam(sample_dir):
    """Find the aligned BAM in the sample dir (or its outs/); None if absent."""
    for sub in (".", "outs"):
        for n in BAM_NAMES:
            p = os.path.join(sample_dir, sub, n)
            if os.path.isfile(p):
                return p
    return None


def discover_samples(project_dir):
    samples = []
    for root, _dirs, _files in os.walk(project_dir):
        fbm = os.path.join(root, "filtered_feature_bc_matrix")
        if not os.path.isdir(fbm):
            continue
        if find_file(fbm, FEATURE_NAMES) and find_file(fbm, MATRIX_NAMES):
            sample_dir = root
            # report keyed by the sample folder, not its 'outs' subdir
            if os.path.basename(os.path.normpath(root)) == "outs":
                parent = os.path.dirname(root)
                if os.path.isdir(os.path.join(parent, "outs")):
                    sample_dir = parent
            samples.append({"sample_dir": sample_dir,
                            "bam": detect_bam(sample_dir)})
    uniq = []
    seen = set()
    for s in sorted(samples, key=lambda x: x["sample_dir"]):
        key = os.path.normpath(s["sample_dir"])
        if key not in seen:
            seen.add(key)
            uniq.append(s)
    return uniq


def find_umap_csv(sample_dir):
    """Locate a UMAP projection CSV under the analysis dir; None if absent.

    PCA projections share the 'projection.csv' filename, so a candidate is
    only accepted if its header truly is a UMAP (UMAP-1 / UMAP-2 columns)."""
    cands = []
    for sub in (os.path.join("outs", "analysis"), "analysis"):
        base = os.path.join(sample_dir, sub)
        if not os.path.isdir(base):
            continue
        for root, _dirs, files in os.walk(base):
            for fn in files:
                if fn.lower().endswith(".csv"):
                    p = os.path.join(root, fn)
                    if "umap" in p.lower() or "projection" in fn.lower():
                        cands.append(p)
    cands.sort(key=lambda p: (len(p), p))
    for p in cands:
        try:
            with open(p, "r", errors="replace") as f:
                head = f.readline().strip()
        except OSError:
            continue
        cols = [c.strip().lower()
                for c in (head.replace(",", "\t").split("\t"))]
        if len(cols) >= 3 and cols[0] == "barcode" \
                and "umap" in cols[1] and "umap" in cols[2]:
            return p
    return None


def read_umap_coords(path):
    """barcode -> (x, y) from a UMAP CSV/TSV (header row skipped); None if bad."""
    if not path or not os.path.isfile(path):
        return None
    try:
        text = open(path, "r", errors="replace").read()
    except OSError:
        return None
    coords = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t") if "\t" in line else line.split(",")
        if len(parts) < 3:
            continue
        try:
            x, y = float(parts[1]), float(parts[2])
        except ValueError:
            continue  # header or malformed line
        coords[parts[0].strip()] = (x, y)
    return coords or None

import plotly.graph_objects as go
# --------------------------------------------------------------------------- #
# Design system (dark-first token CSS) + lightweight JS for theme toggle
# --------------------------------------------------------------------------- #
CSS = """
:root{
  --bg:#0a0f1c; --bg-2:#0c1322; --panel:#121b30; --panel-2:#0f1830; --inset:#0b1424;
  --line:#22314f; --line-2:#33476f; --text:#e8eefb; --text-2:#aab8d6; --mut:#7e8fb2;
  --acc:#6366f1; --acc-2:#8b5cf6; --teal:#2dd4bf; --green:#22c55e; --red:#fb7185;
  --amber:#f5b13d; --rose:#f472b6;
  --shadow:0 18px 42px -16px rgba(2,6,20,.75);
  --rad:16px; --rad-sm:11px;
  --font:-apple-system,'Inter','Segoe UI',Roboto,system-ui,sans-serif;
  --mono:ui-monospace,'SF Mono','JetBrains Mono',Menlo,Consolas,monospace;
}
[data-theme="light"]{
  --bg:#eef2fa; --bg-2:#e7edf8; --panel:#ffffff; --panel-2:#f5f8ff; --inset:#eef3fc;
  --line:#dbe4f3; --line-2:#c3d1ea; --text:#131b30; --text-2:#3f4f74; --mut:#6c7ea8;
  --shadow:0 16px 38px -18px rgba(30,50,100,.28);
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;color:var(--text);font-family:var(--font);line-height:1.5;
  -webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility;
  background:
    radial-gradient(1100px 520px at 88% -8%, rgba(99,102,241,.18), transparent 60%),
    radial-gradient(900px 480px at -6% 2%, rgba(45,212,191,.11), transparent 55%),
    var(--bg);
  background-attachment:fixed;background-repeat:no-repeat;}
.wrap{max-width:1220px;margin:0 auto;padding:30px 22px 56px}
h2,h3{margin:0}
.mut{color:var(--mut)}
code,mono{font-family:var(--mono)}

/* ---- hero ---- */
.hero{position:relative;overflow:hidden;border:1px solid var(--line);
  border-radius:20px;box-shadow:var(--shadow);
  background:linear-gradient(135deg, rgba(99,102,241,.20), rgba(18,27,48,.15) 42%, rgba(45,212,191,.14));
  padding:26px 28px;margin-bottom:20px}
.hero-grid{display:flex;justify-content:space-between;gap:22px;align-items:flex-start;flex-wrap:wrap}
.eyebrow{font-size:12px;letter-spacing:.16em;text-transform:uppercase;color:var(--teal);
  font-weight:700;margin:0 0 9px;display:flex;align-items:center;gap:10px}
.eyebrow .tag{font-size:10.5px;letter-spacing:.12em;padding:3px 8px;border-radius:6px;
  border:1px solid var(--line-2);color:var(--text-2)}
h1{font-size:clamp(27px,4.6vw,42px);line-height:1.04;margin:0;font-weight:800;letter-spacing:-.025em}
h1 .sub{font-weight:500;color:var(--text-2);font-size:.58em;letter-spacing:0}
.chips{display:flex;flex-wrap:wrap;gap:9px;margin-top:18px}
.chip{display:inline-flex;align-items:center;gap:7px;padding:6px 12px;border:1px solid var(--line);
  border-radius:999px;background:rgba(255,255,255,.03);font-size:12.5px;color:var(--text-2)}
.chip b{color:var(--text);font-weight:700}
.hero-summ{display:flex;gap:12px;flex-wrap:wrap}
.bigstat{min-width:104px;border:1px solid var(--line);border-radius:var(--rad-sm);
  background:rgba(255,255,255,.03);padding:14px}
.bigstat .n{font-size:30px;font-weight:800;letter-spacing:-.02em;line-height:1}
.bigstat .n small{font-size:15px;color:var(--text-2);font-weight:600}
.bigstat .n .ok{color:var(--teal)}
.bigstat .n .no{color:var(--red)}
.bigstat .l{font-size:11.5px;color:var(--mut);margin-top:7px;letter-spacing:.01em}
.ghost{appearance:none;border:1px solid var(--line-2);background:rgba(255,255,255,.04);
  color:var(--text-2);border-radius:10px;padding:9px 15px;font-size:12.5px;font-weight:600;
  cursor:pointer;font-family:var(--font);transition:.15s;letter-spacing:.01em;display:inline-flex;
  align-items:center;gap:8px}
.ghost:hover{background:rgba(255,255,255,.09);color:var(--text)}
.ghost svg{width:16px;height:16px}
[data-theme="light"] .ghost{background:#fff;color:#3f4f74}

/* ---- sections + cards ---- */
.sec-h{display:flex;align-items:baseline;justify-content:space-between;gap:12px;margin:32px 2px 14px;flex-wrap:wrap}
.sec-h h2{font-size:19px;font-weight:700;letter-spacing:-.01em;display:flex;align-items:center;gap:11px}
.sec-h h2 .ic{width:24px;height:24px;padding:5px;border:1px solid var(--line-2);border-radius:8px;
  display:inline-grid;place-items:center;color:var(--teal);background:rgba(45,212,191,.10)}
.card{border:1px solid var(--line);border-radius:var(--rad);
  background:linear-gradient(180deg,var(--panel-2),var(--panel));box-shadow:var(--shadow);overflow:hidden}
.sample-card{scroll-margin-top:20px}
.card-head{display:flex;align-items:center;justify-content:space-between;gap:14px;
  padding:18px 20px;border-bottom:1px solid var(--line);flex-wrap:wrap}
.card-head .title{font-size:15.5px;font-weight:800;letter-spacing:-.01em;display:flex;align-items:center;gap:12px;flex-wrap:wrap}
.card-head .idx{font-family:var(--mono);font-size:11.5px;color:var(--mut);
  padding:3px 8px;border:1px solid var(--line);border-radius:6px;background:rgba(255,255,255,.02)}
.card-sub{padding:12px 20px 13px;font-size:12.5px;line-height:1.6;color:var(--mut);
  border-bottom:1px solid var(--line)}
.reading-note{margin:26px 0 18px;padding:14px 18px 15px;font-size:12.5px;line-height:1.65;
  color:var(--mut);border:1px solid var(--line);border-radius:var(--rad);
  background:linear-gradient(180deg,var(--panel-2),var(--panel))}
.reading-note .rn-t{display:block;margin-bottom:6px;font-size:11.5px;font-weight:700;
  letter-spacing:.05em;text-transform:uppercase;color:var(--teal)}

/* ---- pill status ---- */
.pill{display:inline-flex;align-items:center;gap:7px;padding:5px 12px;border-radius:999px;
  font-size:12px;font-weight:700;letter-spacing:.03em;white-space:nowrap}
.pill svg{width:14px;height:14px;flex:none}
.pill.exp{color:#b8f5e8;background:rgba(34,197,94,.16);border:1px solid rgba(45,212,191,.55)}
.pill.silent{color:#ffe0ad;background:rgba(245,177,61,.14);border:1px solid rgba(245,177,61,.5)}
.pill.norow{color:#ffc9d3;background:rgba(251,113,133,.14);border:1px solid rgba(251,113,133,.48)}
[data-theme="light"] .pill.exp{color:#0b5c3c;background:#e2faf0;border-color:#82e2c1}
[data-theme="light"] .pill.silent{color:#6e4405;background:#fff2dc;border-color:#f0c98a}
[data-theme="light"] .pill.norow{color:#8a1a2f;background:#ffe6ec;border-color:#f5a9ba}

/* ---- KPI strip ---- */
.statstrip{display:grid;grid-template-columns:repeat(6,1fr);gap:1px;
  background:var(--line);border-top:1px solid var(--line);border-bottom:1px solid var(--line)}
.statstrip .kpi{background:var(--panel);padding:15px 16px}
@media(max-width:940px){.statstrip{grid-template-columns:repeat(3,1fr)}}
@media(max-width:560px){.statstrip{grid-template-columns:repeat(2,1fr)}}
.kpi .v{font-size:23px;font-weight:800;letter-spacing:-.02em;line-height:1.05}
.kpi .v small{font-size:13px;color:var(--text-2);font-weight:600}
.kpi .v .ok{color:var(--teal)} .kpi .v .warn{color:var(--amber)}
.kpi .k{font-size:11px;color:var(--mut);margin-top:6px;letter-spacing:.02em;line-height:1.35}

/* ---- table ---- */
.tablewrap{overflow:auto;max-height:560px}
table{width:100%;border-collapse:separate;border-spacing:0;font-size:13.5px;min-width:820px}
thead th{position:sticky;top:0;z-index:2;background:var(--inset);color:var(--text-2);
  font-weight:600;text-align:left;padding:12px 14px;letter-spacing:.02em;font-size:11px;
  text-transform:uppercase;border-bottom:1px solid var(--line-2);white-space:nowrap}
thead th.num{text-align:right}
td{padding:12px 14px;border-bottom:1px solid var(--line);color:var(--text-2);white-space:nowrap}
td.num{text-align:right;font-family:var(--mono);font-variant-numeric:tabular-nums}
td.gene{font-weight:700;color:var(--text)}
tbody tr:not(.tot):hover td{background:rgba(99,102,241,.06)}
tr.tot td{background:rgba(99,102,241,.10);color:var(--text);font-weight:700;border-bottom:none}
tr.tot td.pill{background:transparent}
.flag{color:var(--amber);font-weight:800}
.na{color:var(--mut);opacity:.72;font-style:italic}

/* ---- chart grid + panels ---- */
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}
@media(max-width:900px){.grid2{grid-template-columns:1fr}}
@media(max-width:760px){.grid3{grid-template-columns:1fr}}
.panel{border:1px solid var(--line);border-radius:var(--rad);background:var(--panel-2);
  box-shadow:var(--shadow);overflow:hidden}
.panel-h{display:flex;align-items:center;justify-content:space-between;gap:10px;padding:14px 16px 4px}
.panel-h .t{font-size:14px;font-weight:700;display:flex;align-items:center;gap:9px;letter-spacing:-.01em}
.panel-h .t .sw{width:11px;height:11px;border-radius:4px;background:var(--teal);flex:none;
  box-shadow:0 0 10px rgba(45,212,191,.5)}
.panel-h .m{font-size:11.5px;color:var(--mut);font-family:var(--mono)}
.panel .js-plotly-plot{width:100%!important}
.panel .plotly .modebar{opacity:.25;transition:.15s}.panel .plotly .modebar:hover{opacity:1}
.note{padding:2px 17px 13px;font-size:12px;color:var(--mut)}
.empty{padding:34px 20px;color:var(--mut);font-size:13px;text-align:center}
.legend-inline{display:flex;gap:14px;align-items:center;padding:2px 16px 10px;font-size:12px;color:var(--text-2);flex-wrap:wrap}
.legend-inline .lg{display:inline-flex;align-items:center;gap:7px}
.legend-inline .dot{width:11px;height:11px;border-radius:3px}

/* ---- footer ---- */
.foot{margin-top:42px;padding-top:18px;border-top:1px solid var(--line);color:var(--mut);
  font-size:12px;display:flex;justify-content:space-between;gap:14px;flex-wrap:wrap}
.foot code{font-family:var(--mono);color:var(--text-2)}
"""

JS = """
function gcrToggleTheme(){
  var cur = document.documentElement.getAttribute('data-theme') || 'dark';
  var next = cur === 'dark' ? 'light' : 'dark';
  try { localStorage.setItem('gcr-theme', next); } catch (e) {}
  document.documentElement.setAttribute('data-theme', next);
  var font  = next === 'dark' ? '#e8eefb' : '#182338';
  var ticks = next === 'dark' ? '#aab8d6' : '#4a5c82';
  var grid  = next === 'dark' ? 'rgba(148,163,184,0.14)' : 'rgba(40,70,130,0.16)';
  var line  = next === 'dark' ? '#2e4067' : '#c3d1ea';
  if (window.Plotly) {
    document.querySelectorAll('.js-plotly-plot').forEach(function (el) {
      if (!el.data) return;
      try {
        Plotly.relayout(el, {
          'font.color': font,
          'xaxis.gridcolor': grid, 'yaxis.gridcolor': grid,
          'xaxis.linecolor': line, 'yaxis.linecolor': line,
          'xaxis.tickfont.color': ticks, 'yaxis.tickfont.color': ticks,
          'xaxis.title.font.color': ticks, 'yaxis.title.font.color': ticks,
          'legend.font.color': font, 'legend.title.font.color': font
        });
      } catch (e) {}
    });
  }
}
(function gcrInitTheme() {
  var m = 'dark';
  try { m = localStorage.getItem('gcr-theme') || 'dark'; } catch (e) {}
  document.documentElement.setAttribute('data-theme', m);
})();
"""

# Palette for the charts
_EXPR_SCALE = ["#14162e", "#3b2a76", "#7c3aed", "#c026d3", "#f07bb0", "#facc15"]
_POS = "#2dd4bf"
_NEG = "#3a4a66"
_BAR = "#2dd4bf"
_DOT = "#f472b6"


def _p_icon(status):
    if status == "EXPRESSED":
        return ('<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.3" '
                'stroke-linecap="round" stroke-linejoin="round"><path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/>'
                '<polyline points="22 4 12 14.01 9 11.01"/></svg>')
    if status == "SILENT":
        return ('<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.3" '
                'stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86 1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/>'
                '<line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/></svg>')
    return ('<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.3" '
            'stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="10"/>'
            '<line x1="8" y1="8" x2="16" y2="16"/><line x1="16" y1="8" x2="8" y2="16"/></svg>')


# --------------------------------------------------------------------------- #
# Plotly figure builders (dark-default layout, themeable via gcrToggleTheme)
# --------------------------------------------------------------------------- #
def _plotly_theme(fig, title, x_title="", y_title="", height=360, legend=False,
                  legend_title=""):
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#e8eefb", size=11, family="-apple-system,Inter,system-ui,sans-serif"),
        margin=dict(l=46, r=20, t=18, b=34), height=height, autosize=True,
        showlegend=legend,
        legend=dict(orientation="h", yanchor="bottom", y=-0.24, xanchor="left", x=0,
                    font=dict(size=11), bgcolor="rgba(0,0,0,0)"),
        hoverlabel=dict(bgcolor="#0c1526", bordercolor="#2e4067",
                        font=dict(color="#e8eefb", size=12,
                                 family="-apple-system,Inter,sans-serif")),
    )
    fig.update_xaxes(gridcolor="rgba(148,163,184,0.09)", linecolor="#2e4067",
                     zeroline=False, ticks="outside", ticklen=5, tickwidth=1,
                     title_text=x_title, title_font=dict(color="#8a99bd", size=11),
                     tickfont=dict(color="#aab8d6", size=10))
    fig.update_yaxes(gridcolor="rgba(148,163,184,0.09)", linecolor="#2e4067",
                     zeroline=False, ticks="outside", ticklen=5, tickwidth=1,
                     title_text=y_title, title_font=dict(color="#8a99bd", size=11),
                     tickfont=dict(color="#aab8d6", size=10))
    return fig


def _f2html(fig):
    cfg = dict(displaylogo=False,
               modeBarButtonsToRemove=["toImage", "hoverCompareCartesian"],
               scrollZoom=False, showTips=False)
    return fig.to_html(full_html=False, include_plotlyjs=False, config=cfg)


def _panel(title, fig, meta="", cap=""):
    if fig is None:
        return ('<div class="panel"><div class="panel-h"><span class="t">'
                '<span class="sw"></span>' + html.escape(title) + '</span>'
                + ('<span class="m">' + html.escape(meta) + '</span>' if meta else '')
                + '</div><div class="empty">not available</div></div>')
    note = ('<div class="note">' + html.escape(cap) + '</div>' if cap else '')
    meta_h = ('<span class="m">' + html.escape(meta) + '</span>' if meta else '')
    return ('<div class="panel"><div class="panel-h"><span class="t">'
            '<span class="sw"></span>' + html.escape(title) + '</span>' + meta_h
            + '</div>' + _f2html(fig) + note + '</div>')


def _fig_umap_expr(xs, ys, umi, label):
    if len(umi) < 2:
        return None
    t = [math.log10(v + 1) for v in umi]
    cmax = max(t) if max(t) > 0 else 1.0
    fig = go.Figure(go.Scattergl(
        x=list(xs), y=list(ys), mode="markers",
        marker=dict(size=5, opacity=0.82, color=t, colorscale=_EXPR_SCALE,
                    cmin=0, cmax=cmax, showscale=True,
                    colorbar=dict(
                        title=dict(text=label + "  log\u00b9\u2070 UMI", side="right",
                                  font=dict(size=10, color="#c9d3e6")),
                        len=0.92, thickness=11, outlinewidth=1,
                        outlinecolor="#22314f",
                        tickfont=dict(color="#c9d3e6", size=10), tickformat=".1f")),
        customdata=list(umi),
        hovertemplate="x %{x:.2f}<br>y %{y:.2f}<br><b>UMI %{customdata}</b>"
                      "<extra>single cell</extra>",
    ))
    return _plotly_theme(fig, label + "  \u00b7  expression on UMAP",
                         x_title="UMAP-1", y_title="UMAP-2", height=440)


def _fig_umap_bin(dx, nx, dy, ny, label, threshold):
    fig = go.Figure()
    if nx:
        fig.add_trace(go.Scattergl(x=list(nx), y=list(ny), mode="markers",
                                   name="not detected",
                                   marker=dict(size=3.5, color=_NEG, opacity=0.5),
                                   hovertemplate="not detected \u00b7 x %{x:.2f}<br>y %{y:.2f}"
                                                 "<extra></extra>"))
    if dx:
        fig.add_trace(go.Scattergl(x=list(dx), y=list(dy), mode="markers",
                                   name=label + "  detected (\u2265 %g UMI)" % threshold,
                                   marker=dict(size=4.6, color=_POS, opacity=0.92),
                                   hovertemplate="detected \u00b7 x %{x:.2f}<br>y %{y:.2f}"
                                                 "<extra></extra>"))
    return _plotly_theme(fig, label + "  \u00b7  detected vs not detected on UMAP",
                         x_title="UMAP-1", y_title="UMAP-2", height=440,
                         legend=True)


def _fig_histogram(values, label):
    if not values:
        return None
    t = [math.log10(v + 1) for v in values]
    tmin, tmax = min(t), max(t)
    if (tmax - tmin) < 1e-9:
        tmax = tmin + 1e-9
    nb = 28
    counts = [0] * nb
    for v in values:
        tt = math.log10(v + 1)
        counts[min(nb - 1, max(0, int((tt - tmin) / (tmax - tmin) * nb)))] += 1
    labels = [str(max(0, int(round(10 ** (tmin + (i + 0.5) / nb * (tmax - tmin)) - 1))))
              for i in range(nb)]
    fig = go.Figure(go.Bar(
        x=labels, y=counts,
        marker=dict(color="#38bdf8", line=dict(width=0)),
        hovertemplate="\u2248 %{x} UMI<br>%{y} cells<extra>expressing</extra>",
    ))
    fig.update_xaxes(type="category", tickvals=labels[::max(1, nb // 6)])
    fig.update_layout(bargap=0.05)
    return _plotly_theme(fig, "UMI distribution of expressing cells",
                         x_title="UMI per cell", y_title="cells", height=300)


def _fig_boxplot(pos_vals, neg_vals, label):
    fig = go.Figure()
    if neg_vals:
        fig.add_trace(go.Box(y=list(neg_vals), name="not detected",
                             marker=dict(color=_NEG, size=3, opacity=0.5,
                                         ),
                             line=dict(width=1.1),
                             fillcolor="rgba(58,74,102,0.18)"))
    if pos_vals:
        fig.add_trace(go.Box(y=list(pos_vals), name=label + "  detected",
                             marker=dict(color=_POS, size=3, opacity=0.6,
                                         ),
                             line=dict(width=1.2),
                             fillcolor="rgba(45,212,191,0.16)"))
    fig.update_yaxes(type="log")
    return _plotly_theme(fig, "Total UMI per cell  \u00b7  detected vs not detected",
                         y_title="total UMI per cell (log)", height=320,
                         legend=True)


def _fig_coverage(depth, label):
    if not depth:
        return None
    pos = sorted(depth)
    vals = [depth[p] for p in pos]
    hover = "position %%{x} bp \u00b7 coverage %%{y}\u00d7<extra>%s</extra>" % html.escape(label)
    fig = go.Figure(go.Scatter(
        x=pos, y=vals, mode="lines", fill="tozeroy",
        fillcolor="rgba(99,102,241,0.34)", line=dict(color="#818cf8", width=1.5),
        hovertemplate=hover,
    ))
    return _plotly_theme(fig, "Read-start coverage on the transgene contig",
                         x_title="position on contig (bp)", y_title="coverage (\u00d7)",
                         height=320)


def _fig_bars(labels, values, label):
    fig = go.Figure(go.Bar(
        x=labels, y=values,
        marker=dict(color=_BAR, line=dict(width=0)),
        text=["%.1f%%" % v for v in values], textposition="outside",
        textfont=dict(color="#c9d3e6", size=11),
        hovertemplate="%{x}<br>%{y:.1f}%% cells detected<extra></extra>",
    ))
    fig.update_layout(bargap=0.34)
    return _plotly_theme(fig, "% cells detected" + (("\u00b7 " + label) if label else ""),
                         x_title="sample", y_title="% cells (\u2265 threshold)", height=320)


def _fig_dots(labels, means, label):
    order = list(range(len(labels)))
    fig = go.Figure(go.Scatter(
        x=list(means), y=[labels[i] for i in order],
        mode="markers + text",
        marker=dict(size=11, color=_DOT, line=dict(color="#0c1526", width=1.5)),
        text=["%g" % v for v in means], textfont=dict(color="#c9d3e6", size=11),
        hovertemplate="%{y}<br>mean %{x:g} UMI<extra>expressing cells</extra>",
    ))
    fig.update_yaxes(categoryorder="array", categoryarray=labels, autorange="reversed")
    return _plotly_theme(fig, "Mean UMI per expressing cell" + (("\u00b7 " + label) if label else ""),
                         x_title="mean UMI per expressing cell", y_title="", height=320)


def match_contig(genes, contigs):
    if not genes or not contigs:
        return None
    exact = {c.lower(): c for c in contigs}
    for g in genes:
        gl = str(g).lower()
        if gl in exact:
            return (g, exact[gl])
        cands = [c for c in contigs if gl in c.lower()]
        if cands:
            cands = sorted(cands, key=lambda c: len(c))
            return (g, cands[0])
    return None


def _chip(key, val=""):
    v = html.escape(val) if val else ""
    body = ("<b>" + html.escape(key) + "</b>" + ((" " + v) if v else ""))
    return "<span class=\"chip\">" + body + "</span>"


def _num_cell(v, fn):
    if v is None:
        return "<td class=\"num na\">n/a</td>"
    try:
        return "<td class=\"num\">" + fn(v) + "</td>"
    except Exception:
        return "<td class=\"num na\">n/a</td>"


def _pct_cell(v):
    if v is None:
        return "<td class=\"num na\">n/a</td>"
    if isinstance(v, (int, float)):
        return "<td class=\"num\">%.1f%%</td>" % v
    return "<td class=\"num na\">n/a</td>"


# --------------------------------------------------------------------------- #
# public render entry points (keep the original signatures)
# --------------------------------------------------------------------------- #
# Shown once at the top of the report, above the per-sample tables (not repeated
# for every sample).
READING_NOTE = ("A gene is reported as &ldquo;detected&rdquo; in a cell when the cell &rsquo;s "
                "UMI count for that gene is equal to or greater than the detection threshold. "
                "This rule is applied everywhere the term "
                "&ldquo;detected&rdquo; appears in this report. "
                "The extra &ldquo;cells &ge; N UMI&rdquo; columns count cells at fixed low "
                "cutoffs, as a sensitivity check independent of the detection threshold. "
                "Mean and median UMI are computed over expressing cells only "
                "(cells with at least 1 UMI). "
                "If a gene&rsquo;s row has no numbers (status NO-ROW), the "
                "transgene row was not found in this sample&rsquo;s feature "
                "table &mdash; the gene was not part of the genome used for "
                "this run. "
                "The &ldquo;reads (called cells)&rdquo; column counts "
                "transcriptome reads aligned on the gene&rsquo;s contig "
                "(primary alignments only) that come from called cells, taken "
                "from the aligned read file (BAM); it is n/a when no BAM "
                "is available.")


def render_stats_table(sid, genes, stats, threshold):
    s = html.escape(str(sid))
    n_exp = sum(1 for g in genes if (stats.get(g) or {}).get("status") == "EXPRESSED")
    n_sil = sum(1 for g in genes if (stats.get(g) or {}).get("status") == "SILENT")
    n_norow = sum(1 for g in genes if (stats.get(g) or {}).get("status") == "NO-ROW")
    ct = next((st.get("cells_total") for st in stats.values() if st.get("cells_total")), 0) or 0

    def num(v, spec="{:,}"):
        return "<td class=num>" + (spec.format(v) if v is not None else "n/a") + "</td>"
    def pcell(p):
        return "<td class=num>" + (f"{p:.1f}%" if p is not None else "n/a") + "</td>"
    def pct(n, denom):
        return (100.0 * n / denom) if (denom and n is not None) else None

    # Sensitivity columns at fixed 1 and 3 UMI cutoffs; a column equal to the
    # active detection threshold is skipped so a count is never shown twice.
    thrs = []
    if threshold != 1:
        thrs.append(("ge1", "cells &ge; 1 UMI"))
    if threshold != 3:
        thrs.append(("ge3", "cells &ge; 3 UMI"))

    head = ("<thead><tr><th>gene</th><th>status</th>"
            f"<th class=num>cells &ge; {int(threshold)} UMI</th><th class=num>% cells</th>")
    for _k, lbl in thrs:
        head += f"<th class=num>{lbl}</th><th class=num>% cells</th>"
    head += ("<th class=num>total UMIs</th><th class=num>reads (called cells)</th>"
             "<th class=num>mean UMI</th>"
             "<th class=num>median UMI</th>"
             "<th class=num>endogenous cells</th><th class=num>endogenous UMIs</th>"
             "</tr></thead>")
    colcount = 2 + 2 + 2 * len(thrs) + 6

    rows = []
    for g in genes:
        d = stats.get(g)
        if d is None or d.get("status") == "NO-ROW":
            rows.append(f"<tr><td class=gene>{html.escape(g)}</td>"
                        f"<td><span class='pill norow'>{_p_icon('NO-ROW')} NO-ROW</span></td>"
                        f"<td colspan={colcount} class=na>no transgene feature row in this output "
                        f"(check the reference genome)</td></tr>")
            continue
        st = d["status"]
        cls = {"EXPRESSED": "exp", "SILENT": "silent"}.get(st, "norow")
        det = d.get("det") or 0
        cells = ["<tr>", f"<td class=gene>{html.escape(g)}</td>",
                 f"<td><span class='pill {cls}'>{_p_icon(st)}{st}</span></td>",
                 num(det), pcell(pct(det, ct))]
        for k, _lbl in thrs:
            v = d.get(k)
            cells += [num(v), pcell(pct(v, ct))]
        cells.append(num(d.get("umi")))
        cells.append(num(d.get("reads")))
        cells.append(num(d.get("mean"), "{:.2f}"))
        cells.append(num(d.get("median"), "{:.2f}"))
        cells.append(num(d.get("endo_cells")))
        cells.append(num(d.get("endo_umi")))
        cells.append("</tr>")
        rows.append("".join(cells))

    return (f"<section class=card sample-card id=s-{s}>"
            "<div class=card-head><div class=title><span class=idx>" + _p_icon("EXPRESSED") +
            f"</span>" + s + " &nbsp;"
            f"<span class=chip><b>{n_exp}</b> expressed / {len(genes)} gene(s)</span>"
            + (f"<span class=chip>{n_sil} silent</span>" if n_sil else "")
            + (f"<span class=chip>{n_norow} no-row</span>" if n_norow else "")
            + f"</div></div>"
            + (f"<div class=card-sub>{ct:,} cells total &middot; "
               f"detection threshold &ge; {int(threshold)} UMI/cell</div>" if ct else "")
            + "<div class=tablewrap><table>" + head + "<tbody>"
            + "".join(rows) + "</tbody></table></div>"
            "</section>")



def render_per_sample_table(results, genes):
    if not results or not genes or len(results) < 2:
        return ""
    labels = [sid for sid, _st, _ct in results]
    first = genes[0]
    dets, means = [], []
    expr_cnt = 0
    for _sid, stats, _ct in results:
        d = stats.get(first) or {}
        cells = d.get("cells_total") or 0
        det = d.get("det") or 0
        dets.append((100.0 * det / cells) if cells else 0.0)
        means.append(d.get("mean") or 0.0)
        if d.get("status") == "EXPRESSED":
            expr_cnt += 1
    avg = sum(dets) / len(dets) if dets else 0.0
    mx = max(means) if means else 0.0
    cards = ("<div class=grid3 style=margin-bottom:16px>"
             f"<div class=bigstat><div class=n>{expr_cnt}/{len(results)}</div>"
             f"<div class=l>samples where {html.escape(first)} is expressed</div></div>"
             f"<div class=bigstat><div class=n>{avg:.1f}%</div>"
             "<div class=l>mean share of cells that are detected (transgene UMI &ge; threshold)</div></div>"
             f"<div class=bigstat><div class=n>{mx:.2f}</div>"
             "<div class=l>max mean UMI per expressing cell</div></div></div>")
    bar = _panel("% cells detected", _fig_bars(labels, dets, first), "per sample, at the detection threshold")
    dot = _panel("Mean UMI per expressing cell", _fig_dots(labels, means, first), "expressing cells only (UMI \u2265 1)")
    return ('<div class=sec-h><h2><span class=ic>' + _p_icon("EXPRESSED") + "</span>Across samples</h2></div>"
            + cards + "<div class=grid2>" + bar + dot + "</div>")


def render_plots(plot_data, genes, threshold, coverage=True):
    if not plot_data:
        return ""
    out = []
    for p in plot_data:
        sid = p.get("id") or ""
        cells = p.get("cells_total") or 0
        coords = p.get("umap_coords") or None
        bc_total = p.get("bc_total") or {}
        bc_gene = p.get("bc_gene") or {}
        stats = p.get("stats") or {}
        ncoord = "available" if coords else "not available"
        out.append('<div class="sec-h"><h2><span class="ic">' + _p_icon("EXPRESSED") + "</span>"
                   + html.escape(sid) + "</h2>"
                   + _chip("cells", "{:,}".format(cells))
                   + _chip("UMAP", ncoord)
                   + _chip("detection threshold", "\u2265 %g UMI per cell" % threshold) + "</div>")

        # shared UMAP frame for this sample (barcode order is stable)
        if coords:
            order = list(coords.keys())
            cx = [coords[bc][0] for bc in order]
            cy = [coords[bc][1] for bc in order]
        else:
            order = []
            cx = cy = []

        for g in genes:
            st = stats.get(g) or {}
            pos_vals = st.get("pos_vals") or []
            if coords:
                cmi = [max(0.0, bc_gene.get(bc, {}).get(g, 0)) for bc in order]
                d_x, n_x, d_y, n_y = [], [], [], []
                for bc in order:
                    u = bc_gene.get(bc, {}).get(g, 0)
                    if u >= threshold:
                        d_x.append(coords[bc][0]); d_y.append(coords[bc][1])
                    else:
                        n_x.append(coords[bc][0]); n_y.append(coords[bc][1])
            else:
                cmi, d_x, n_x, d_y, n_y = [], [], [], [], []
            pv = [bc_total[bc] for bc in bc_total if bc_gene.get(bc, {}).get(g, 0) >= threshold]
            nv = [bc_total[bc] for bc in bc_total if bc_gene.get(bc, {}).get(g, 0) < threshold]

            if st.get("row_present") is False:
                norow_msg = ("this gene is not present in this sample \u2019s reference "
                             "(NO-ROW), so there are no cells to show")
                out.append('<div class="grid2" style="margin-bottom:16px">')
                out.append(_panel(g + " \u00b7 expression on UMAP", None, norow_msg))
                out.append(_panel(g + " \u00b7 detected vs not detected", None, norow_msg))
                out.append("</div>")
                continue

            out.append('<div class="grid2" style="margin-bottom:16px">')
            out.append(_panel(
                g + " \u00b7 expression on UMAP", _fig_umap_expr(cx, cy, cmi, g),
                meta="{}, cells".format(len(cx)),
                cap="Colour = transgene UMI per cell, on a log\u00b9\u2070 scale; dark = low, "
                    "bright = high. Drag to zoom, shift-drag to pan, double-click to reset, "
                    "hover a point for its UMI."))
            out.append(_panel(
                g + " \u00b7 detected vs not detected", _fig_umap_bin(d_x, n_x, d_y, n_y, g, threshold),
                meta="{} detected / {} not detected".format(len(d_x), len(n_x)),
                cap="A cell is shown as detected when its transgene UMI count for this gene "
                    "is equal to or greater than %g \u00b7 all other cells are shown as not "
                    "detected." % threshold))
            out.append("</div>")
            out.append('<div class="grid2" style="margin-bottom:16px">')
            out.append(_panel(
                "UMI distribution of expressing cells", _fig_histogram(pos_vals, g),
                meta="{}, expressing cells".format(len(pos_vals)),
                cap="Each bar is the number of expressing cells at that transgene UMI "
                    "level (bins are logarithmic)."))
            out.append(_panel(
                "Total UMI per cell \u00b7 detected vs not detected", _fig_boxplot(pv, nv, g),
                meta="{} detected vs {} not detected cells".format(len(pv), len(nv)),
                cap="Total UMI of the whole transcriptome per cell. Left box: cells where "
                    "the gene is not detected; right box: cells where it is detected "
                    "(log scale)."))
            out.append("</div>")

    if coverage:
        out.extend(_cov_blocks(plot_data))
    return "\n".join(out)


def _cov_blocks(plot_data):
    """Transgene-locus coverage panels (a short contig per sample; kept light
    on purpose so it can also live in the index file without adding weight)."""
    out = []
    cov_done = False
    for p in plot_data:
        cov = p.get("cov")
        if cov is None:
            continue
        g = cov.get("gene") or ""
        if cov.get("no_contig"):
            out.append(_panel("Coverage \u00b7 " + g, None,
                              "no contig matching any of the listed genes "
                              "was found among the {} contigs in the aligned "
                              "read file (BAM)".format(cov.get("contigs_n") or 0)))
            continue
        if cov.get("depth"):
            out.append(_panel(
                "Coverage \u00b7 " + g, _fig_coverage(cov["depth"], g),
                meta="contig {} \u00b7 {:,} bp \u00b7 sample {}".format(
                    html.escape(cov.get("contig") or ""), cov.get("length") or 0,
                    html.escape(p.get("id") or "")),
                cap="Number of reads covering each position along the transgene; "
                    "a flat, even plateau across the whole length indicates the "
                    "transgene is transcribed as a full-length construct."))
            cov_done = True
        else:
            out.append(_panel("Coverage \u00b7 " + g, None,
                              "no aligned reads were found on contig " + (cov.get("contig") or "")))
    if not cov_done and not any(p.get("bam") for p in plot_data):
        out.append('<p class="note">Coverage plot skipped: no aligned read file (BAM) '
                   'was found next to the sample outputs.</p>')
    return out


def render_coverage(plot_data):
    return "\n".join(_cov_blocks(plot_data))


def render_detail_nav(parts, total_genes):
    """Landing-page block: one link per plot detail file (its transgene range)."""
    items = []
    for p in parts:
        g = p["genes"]
        rng = g[0] if len(g) == 1 else "{} \u2013 {}".format(g[0], g[-1])
        items.append(
            '<li><a href="{}"><b>{}</b> &middot; {}\u00a0</a> '
            '<span class="rn-sub">{}/{} expressed</span></li>'.format(
                html.escape(p["fname"]), p["label"], html.escape(rng),
                p["n_exp"], len(g)))
    return ('<div class="card sample-card" id="detail-nav">'
            '<div class=card-head><div class=title><span class=ic>' + _p_icon("EXPRESSED") + '</span>'
            'Interactive plots \u00b7 detail files</div></div>'
            '<p>All {} transgenes are summarised in the tables on this page. The '
            'interactive plots are split into {} files (batches) so no single file is '
            'too large to open. Pick a file below to view its plots.</p>'.format(
                total_genes, len(parts))
            + '<ul class=detail-list>' + "".join(items) + '</ul></div>')


def render_back_note(index_fname, label, lo, hi, n_parts, k):
    """First line of a detail file: link back to the index + its own range."""
    return ('<p class=note><a href="{}">&larr; back to index</a> '
            '(tables for all transgenes) &middot; transgenes {}&ndash;{} &middot; '
            'plot file {}/{} &middot; {} transgene{pl}</p>'.format(
                html.escape(index_fname), lo, hi,
                int(label.replace("part", "").strip()) if label else 0,
                n_parts, k, pl="" if k == 1 else "s"))


def render_html(title, meta, sections, plots):
    m = {}
    for line in meta:
        if ": " in line:
            k, v = line.split(": ", 1)
            m[k.strip().lower()] = v.strip()
    gene = m.get("genes", "")
    source = m.get("mode", "")
    gen_count = m.get("samples", "")
    threshold = m.get("detection threshold", m.get("thresholddetected", m.get("threshold", "")))
    generated = m.get("generated", "")
    gene_main = gene.split(",")[0].strip() if gene else ""
    # counts (reliable from pills in the stat tables)
    body = "\n".join(sections)
    expr = body.count("pill exp")

    def _meta_chip(key, label=None):
        if key not in m:
            return ""
        lab = label or key
        return _chip(lab, m[key])

    chips = "".join([
        _chip("genes", gene),
        _chip("samples", gen_count),
        _chip("detection threshold", threshold),
        _chip("generated", generated),
    ])

    # plotly JS once — only when THIS file actually has interactive figures (a
    # pure-table index stays small); fallback to CDN if import fails so the file
    # still renders
    if "Plotly.newPlot" not in plots:
        js_block = ""
    else:
        try:
            import plotly.offline as _po
            _plotly_js = _po.get_plotlyjs()
            has_js = True
        except Exception:
            _plotly_js = ""
            has_js = False
        if not has_js:
            js_block = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>\n'
        else:
            js_block = "<script>" + _plotly_js + "</script>\n"

    hero_cards = (
        "<div class=\"hero-summ\">"
        "<div class=\"bigstat\"><div class=\"n\">"
        + (f"{expr}<small class=\"ok\">/{gen_count}</small>" if gen_count else f"{expr}")
        + "</div><div class=\"l\">expressed samples \u00b7 " +
          (html.escape(gene_main) if gene_main else "") + "</div></div>"
        "<div class=\"bigstat\"><div class=\"n\">" + (gen_count or "\u2013") +
        "</div><div class=\"l\">samples analysed</div></div>"
        "<div class=\"bigstat\"><div class=\"n\">" +
        (html.escape(threshold.split()[0]) if threshold else "\u2013") +
        "</div><div class=\"l\">detection threshold (UMI/cell)</div></div>"
        "</div>"
    )

    return (
        "<!doctype html>\n<html lang=\"en\" data-theme=\"dark\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<meta name=\"color-scheme\" content=\"dark light\">\n"
        "<title>" + html.escape(title) + " \u00b7 " + (html.escape(gene_main) if gene_main else "") + "</title>\n"
        "<script>(function(){var m='dark';try{m=localStorage.getItem('gcr-theme')||'dark';catch(e){}document.documentElement.setAttribute('data-theme',m);})();</script>\n"
        "<style>" + CSS + "</style>\n"
        + js_block + "\n"
        "</head>\n<body>\n<div class=\"wrap\">\n"
        # hero
        "<header class=\"hero\">\n<div class=\"hero-grid\">\n"
        "<div>\n"
        "<p class=\"eyebrow\">Custom / transgene expression</p>"
        "<h1>" + (html.escape(gene_main) if gene_main else html.escape(title)) +
        " <span class=\"sub\">&#183; expression report</span></h1>\n"
        "<div class=\"chips\">" + chips + "</div>\n"
        "</div>\n"
        + hero_cards +
        "</div>\n"
        "<div style=\"display:flex;justify-content:flex-end;margin-top:4px\">\n"
        "<button class=\"ghost\" onclick=\"gcrToggleTheme()\">"
        "<svg viewBox=\"0 0 24 24\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.8\" "
        "stroke-linecap=\"round\" stroke-linejoin=\"round\"><path d=\"M12 3a6 6 0 0 0 9 9 9 9 0 1 1-9-9Z\"/></svg>"
        "Toggle theme</button>\n"
        "</div>\n</header>\n"
        # sections
        + "\n".join(sections) + "\n"
        # plots
        + plots + "\n"
        # footer
        + "<footer class=\"foot\">"
        + "<div><code>custom_gene_report.py</code> &middot; Plotly &middot; self-contained &middot; offline"
        + (("&middot; source: <code>" + html.escape(source) + "</code>") if source else "") + "</div>"
        + ("<div><code>python custom_gene_report.py --project-dir ... --genes " +
          (html.escape(gene) if gene else "") + "</code></div>" if gene else "")
        + "</footer>\n"
        "</div>\n\n<script>" + JS + "</script>\n</body>\n</html>\n"
    )

# --------------------------------------------------------------------------- #
# safe write (overwrites any existing report)
# --------------------------------------------------------------------------- #
def safe_write(path, content):
    if os.path.isdir(path):
        raise OSError(f"{path} is a directory -- it should be a file "
                      f"(--out takes a folder; the report name is derived "
                      f"from the input folder name)")
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return path


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="HTML report of custom-gene (transgene) expression from "
                    "cellranger count output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--sample-dir", help="one sample folder (parent of the "
                                           "cellranger count `outs` directory)")
    grp.add_argument("--project-dir", help="project root; all sample `outs` under it")
    ap.add_argument("--genes", required=True,
                    help="gene name(s): file (one per line) or inline "
                         "comma/space-separated exact names")
    ap.add_argument("--min-umi", type=int, default=3,
                    help="per-cell UMI threshold for 'detected'/status (>= 1)")
    ap.add_argument("--bam", default=None,
                    help="BAM for the coverage track (default: auto-detect "
                         "possorted_genome_bam.bam / possorted_bam.bam / "
                         "gex_possorted_bam.bam in the sample dir)")
    ap.add_argument("--out", "-o", default=None,
                    help="FOLDER to save the report into (default: the input folder); "
                         "the main file is always named custom_genes_<input_dir_name>.html, "
                         "like custom_gene_check.sh's CSV")
    ap.add_argument("--per-file", type=int, default=50,
                    help="max transgenes per plot detail file. When the gene list is "
                         "longer than this, the report is split into an index (tables for "
                         "all genes + coverage tracks) plus N detail files (each with the "
                         "tables for its batch plus the plots, this many per file). "
                         "0 = one combined file (never split)")
    ap.add_argument("--version", action="version", version=VERSION)
    ap.add_argument("-h", "--help", action="help",
                    help="show this help message and exit")
    args = ap.parse_args(argv)

    if args.min_umi < 1:
        raise SystemExit("ERROR: --min-umi must be an integer >= 1")
    genes = list(dict.fromkeys(read_genes(args.genes)))
    if not genes:
        raise SystemExit("ERROR: no gene names found in --genes")

    if args.sample_dir:
        bam = args.bam or detect_bam(args.sample_dir)
        samples = [{"sample_dir": args.sample_dir, "bam": bam}]
    else:
        samples = discover_samples(args.project_dir)
        if args.bam:
            for s in samples:
                s["bam"] = args.bam
    if samples:
        # --out names the FOLDER the report is saved into; the report itself is
        # always custom_genes_<input_dir_name>.html (same name as
        # custom_gene_check.sh's CSV in the input folder).
        base_dir = os.path.normpath(args.sample_dir or args.project_dir)
        fname = "custom_genes_%s.html" % os.path.basename(base_dir)
        if args.out is None:
            out_dir = base_dir                     # default: inside the input folder
        else:
            if os.path.exists(args.out):
                if not os.path.isdir(args.out):
                    raise SystemExit(
                        f"ERROR: --out {args.out} is a file, not a folder. "
                        f"Pass the folder to save into, e.g. "
                        f"--out {os.path.dirname(os.path.normpath(args.out))}")
                out_dir = os.path.normpath(args.out)
            else:
                try:
                    os.makedirs(args.out)
                except OSError as e:
                    raise SystemExit(
                        f"ERROR: could not create --out folder {args.out} ({e})")
                out_dir = os.path.normpath(args.out)
        args.out = os.path.join(out_dir, fname)
    else:
        raise SystemExit("ERROR: no samples with a filtered_feature_bc_matrix found")

    results = []     # (sample_id, stats, cells_total) for the tables
    plot_data = []   # full per-sample dicts for render_plots
    for _n, s in enumerate(samples, 1):
        sid = os.path.basename(os.path.normpath(s["sample_dir"]))
        if sys.stderr.isatty():
            n_g = len(genes)
            print(f"[{_n}/{len(samples)}] {sid}: scanning matrix for {n_g} gene(s) "
                  f"(min_umi={args.min_umi})...",
                  file=sys.stderr, flush=True)
        data = compute_sample(s["sample_dir"], s.get("bam"), genes, args.min_umi)
        data["id"] = sid
        # per sample: `samtools depth` for the coverage track, plus one
        # `samtools view -r` pass for the read count in called cells
        cov = sample_bam_summary(data, genes)
        data["cov"] = cov
        if cov and cov.get("reads"):
            for g, n in cov["reads"].items():
                data["stats"][g]["reads"] = n
        found = find_umap_csv(s["sample_dir"])
        data["umap_coords"] = read_umap_coords(found) if found else None
        plot_data.append(data)
        results.append((sid, data["stats"], data["cells_total"]))

    # the "detected"/column definitions are shown once here, above the
    # sample tables -- not repeated inside every sample card
    sections = [
        ('<div class=reading-note><span class=rn-t>How to read the tables</span>'
         + READING_NOTE + "</div>")
    ]
    meta = [f"mode: {'project' if args.project_dir else 'sample'}",
            f"genes: {', '.join(genes)}",
            f"samples: {len(results)}",
            f"detection threshold: {args.min_umi} UMI/cell",
            f"generated: {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    for sid, stats, _ct in results:
        sections.append(render_stats_table(sid, genes, stats, args.min_umi))
    sections.append(render_per_sample_table(results, genes))

    per_file = args.per_file
    if not isinstance(per_file, int) or per_file < 0:
        raise SystemExit("ERROR: --per-file must be an integer >= 0")
    per_file = int(per_file)
    index_path = args.out                     # always the index / landing file
    out_dir = os.path.dirname(index_path) or "."
    stem = os.path.basename(index_path)
    if stem.endswith(".html"):
        stem = stem[:-len(".html")]
    n_g = len(genes)
    do_split = (per_file is not None) and (per_file > 0) and (n_g > per_file)
    num_parts = (n_g + per_file - 1) // per_file if (per_file is not None and per_file > 0) else 1

    if sys.stderr.isatty():
        label = "split report" if do_split else "report"
        print(f"rendering {label} ({len(results)} sample(s))...",
              file=sys.stderr, flush=True)

    title = "Custom gene expression report"
    texts = []
    if not do_split:
        # single combined file (current behaviour): tables + plots + coverage
        plots = render_plots(plot_data, genes, args.min_umi)
        texts.append((index_path, render_html(title, meta, sections, plots)))
    else:
        width = max(2, len(str(num_parts)))
        parts = []
        for b in range(num_parts):
            genes_b = genes[b * per_file:(b + 1) * per_file]
            n_exp = sum(
                sum(1 for g in genes_b
                    if (stats.get(g) or {}).get("status") == "EXPRESSED")
                for _sid, stats, _ct in results)
            parts.append({
                "fname": f"{stem}_part{b + 1:0{width}d}.html",
                "genes": genes_b,
                "n_exp": n_exp,
                "k": len(genes_b),
                "label": f"part {b + 1:0{width}d}",
                "lo": b * per_file + 1,
                "hi": b * per_file + len(genes_b),
            })
        # index: full tables (all genes) + light coverage + links to the parts
        sections_idx = list(sections)
        sections_idx.insert(1, render_detail_nav(parts, n_g))
        texts.append((index_path,
                      render_html(title, meta, sections_idx,
                                  render_coverage(plot_data))))
        # detail files: the interactive plots, batched, each self-contained
        for p in parts:
            genes_b = p["genes"]
            sections_b = ([render_back_note(os.path.basename(index_path), p["label"],
                                            p["lo"], p["hi"], num_parts, p["k"])]
                          + [render_stats_table(sid, genes_b, stats, args.min_umi)
                             for sid, stats, _ct in results])
            meta_b = [ln if not ln.startswith("genes:")
                      else "genes: " + ", ".join(genes_b) for ln in meta]
            title_b = (f"Custom gene expression \u00b7 {p['label']} "
                       f"({p['lo']}\u2013{p['hi']} of {n_g})")
            plots_b = render_plots(plot_data, genes_b, args.min_umi, coverage=False)
            texts.append((os.path.join(out_dir, p["fname"]),
                          render_html(title_b, meta_b, sections_b, plots_b)))

    # write the index first, then the detail files (a permission error hits the
    # index path, and the existing read-only/ path hints apply unchanged)
    written = []
    for path, text in texts:
        try:
            written.append(safe_write(path, text))
        except OSError as e:
            if e.errno in (errno.EACCES, errno.EROFS):
                hint = ("  The output folder is probably read-only for this user.\n"
                        "  Re-run with an explicit writable location, e.g.:\n"
                        "    ... --out /writable/path")
            else:
                hint = ("  Check the --out folder: it must be a folder the user "
                        "can write to (the report name is derived automatically).")
            raise SystemExit(f"ERROR: could not save report to {path} ({e}).\n{hint}")

    per = ";  ".join(
        f"{sid}: " + " ".join(f"{g}={stats[g]['status']}" for g in genes)
        for sid, stats, _ct in results)
    if len(written) == 1:
        print(f"HTML report written to: {written[0]}")
        print(f"  {per}")
    else:
        print(f"HTML report written to folder: {out_dir}")
        for p in written:
            print(f"  - {os.path.basename(p):34s} ({os.path.getsize(p) // 1024} KB)")
        print(f"  open the index first ({os.path.basename(index_path)}): all {n_g} "
              f"transgenes are summarised there, and each part file holds the "
              f"interactive plots for its batch")
        print(f"  {per}")


if __name__ == "__main__":
    main()
