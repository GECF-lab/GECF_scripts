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
      - one interactive UMAP per sample: one canvas re-coloured client-side
        as you pick a gene (dropdown; expression view, log10 UMI colour scale)
      - coverage track on the transgene contig, one per gene: it follows the
        Gene selector of the interactive plot set (needs BAM + samtools)

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
import json
import os
import random
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
    """-> list[(idx1based, feature_id, name, feature_type_or_None)].

    feature_type is None for lines without a 3rd column. Multi-modal
    3' runs add a feature_type column ('Gene Expression', 'Antisense',
    ...); the caller may then restrict per-cell total UMI to the
    'Gene Expression' rows. If every line lacks the column, feature_type
    is None everywhere and the caller behaves exactly as before.
    """
    out = []
    for i, line in enumerate(read_lines(path), start=1):
        parts = line.split("\t")
        fid = parts[0]
        name = parts[1] if len(parts) > 1 and parts[1].strip() else fid
        ftype = parts[2].strip() if len(parts) > 2 and parts[2].strip() else None
        out.append((i, fid, name, ftype))
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
        # -a reports zero-depth positions (gaps are not drawn as ramps) and
        # -d 0 removes the default 8000x cap, so the true maximum depth shows.
        out = subprocess.check_output(
            [st, "depth", "-a", "-d", "0", "-r", contig, bam],
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

    One `samtools view -r <contig>` pass. Returns {gene: n}; None when the
    scan could not run at all (the caller must then show n/a, never a
    silent 0)."""
    if not bam or not os.path.isfile(bam):
        return None
    st = samtools_cmd()
    if st is None:
        return None
    try:
        out = subprocess.check_output(
            [st, "view", bam, contig],   # region is a positional arg for `view`
            stderr=subprocess.DEVNULL, timeout=900, text=True)
    except Exception:
        return None
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
    """Per-sample BAM summary: reads per gene (called cells) and a base
    coverage track on EACH gene's own transgene contig.

    Every gene that matches a contig (BAM @SQ name) is scanned on its own
    contig: one `samtools view -r <contig>` pass each, counting primary
    alignments tagged GN:<gene> from called cells, and one
    `samtools depth -a` pass for its coverage track. So ALL genes in the
    list get their own count AND their own track -- earlier versions
    stopped after the first matched gene, leaving every other gene's
    'reads (called cells)' at 0 and only its contig ever plotted.

    Returns a dict
      {reads: {gene: int|None},
       contigs_n: int,
       no_contig: bool,                # no listed gene matches any contig
       genes: {gene: {contig: str, length: int, pos0: int,
                      depth: [int, ...]}}}
                      # dense depth from pos0 (empty list = contig present
                      # but with no aligned reads at all)
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
    reads = {g: None for g in genes}
    tracks = {}
    for g in genes:
        m = match_contig([g], contigs)
        if m is None:
            continue
        gene, cname = m
        # A failed BAM scan returns None: keep the gene at None ("n/a" in
        # the table) instead of a silent 0. Genes with no exactly-named
        # contig also stay None.
        gr = bam_gene_reads(bam, cname, [gene], called)
        if gr is not None:
            reads[gene] = gr.get(gene, 0)
        d = bam_gene_depth(bam, cname)
        if d:
            lo, hi = min(d), max(d)
            # -a makes positions contiguous first->last; encode the dense
            # stretch as (pos0, [depth...]). Transgene contigs are short;
            # the cap only guards against a same-named reference contig.
            if hi - lo < 2_000_000:
                tracks[gene] = {"contig": cname,
                                "length": contigs.get(cname) or 0,
                                "pos0": lo,
                                "depth": [d.get(i, 0) for i in range(lo, hi + 1)]}
        else:
            # contig matched but zero coverage: keep an (empty) track so the
            # panel can say "present in the BAM, no aligned reads"
            tracks[gene] = {"contig": cname,
                            "length": contigs.get(cname) or 0,
                            "pos0": 1, "depth": []}
    return {"reads": reads,
            "contigs_n": len(contigs),
            "no_contig": not tracks,
            "genes": tracks}


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
    gex_rows = None  # multi-modal runs (features.tsv 3rd column): 1-based rows
                     # that are 'Gene Expression'; None when the column is
                     # absent, in which case per-cell total UMI counts all rows
    for idx, fid, name, ftype in read_features(fbm["features"]):
        if ftype is not None:
            if gex_rows is None:
                gex_rows = set()
            if ftype.lower() == "gene expression":
                gex_rows.add(idx)
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
                           threshold, stats, bc_total, bc_gene, gex_rows)
    except Exception as exc:
        # no scipy, or a matrix mmread can't handle -> stdlib line scan
        try:
            sys.stderr.write(
                "note: fast matrix scan unavailable (%s); "
                "falling back to the line-by-line scan\n"
                % exc.__class__.__name__)
        except Exception:
            pass
        # The numpy scan may have partially filled these before it failed;
        # the stdlib scan recomputes from scratch, so reset first or every
        # count would be double-counted.
        bc_total.clear()
        bc_gene.clear()
        for g in genes:
            s = stats[g]
            s["n_pos"] = 0
            s["ge1"] = 0
            s["ge3"] = 0
            s["det"] = 0
            s["umi"] = 0
            s["pos_vals"] = []
            s["endo_cells"] = 0
            s["endo_umi"] = 0
        _scan_matrix_python(fbm["matrix"], barcodes, tg_rows, endo_rows,
                            threshold, stats, bc_total, bc_gene, gex_rows)

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
                       threshold, stats, bc_total, bc_gene, gex_rows=None):
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

    try:
        X = mmread(str(matrix_path), spmatrix=True)
    except TypeError:
        # SciPy < 1.12: mmread has no `spmatrix` keyword. Read without it
        # and convert whatever (sparse) matrix we get to CSR ourselves.
        X = mmread(str(matrix_path))
        if hasattr(X, "tocsr"):
            X = X.tocsr()
    ncells = len(barcodes)
    if X.shape[1] != ncells:
        raise ValueError("matrix has %d cells, barcodes.tsv has %d"
                         % (X.shape[1], ncells))
    csr = X if X.format == "csr" else X.tocsr()

    # total UMI per cell: one C pass over every nonzero in the matrix.
    # Multi-modal runs (gex_rows given): count only the 'Gene Expression'
    # rows; otherwise all rows, as before.
    if gex_rows is not None:
        rowkeep = np.zeros(csr.shape[0], dtype=bool)
        for r1 in gex_rows:
            if 1 <= r1 <= csr.shape[0]:
                rowkeep[r1 - 1] = True
        base = csr[rowkeep]
        if hasattr(base, "tocsr"):
            base = base.tocsr()
    else:
        base = csr
    colsum = np.asarray(base.sum(axis=0)).ravel()
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
                        threshold, stats, bc_total, bc_gene, gex_rows=None):
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
            # multi-modal runs: per-cell total UMI counts only 'Gene
            # Expression' rows; gex_rows None -> all rows, as before
            if gex_rows is None or r in gex_rows:
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
# Design system (light-theme token CSS) -- fixed day theme, no toggle
# --------------------------------------------------------------------------- #
CSS = """
:root{
  --bg:#0a0f1c; --bg-2:#0c1322; --panel:#121b30; --panel-2:#0f1830; --inset:#0b1424;
  --line:#22314f; --line-2:#33476f; --text:#e8eefb; --text-2:#aab8d6; --mut:#7e8fb2;
  --acc:#6366f1; --acc-2:#8b5cf6; --teal:#2dd4bf; --green:#22c55e; --red:#fb7185;
  --field:#0b1424;
  --amber:#f5b13d; --rose:#f472b6;
  --shadow:0 18px 42px -16px rgba(2,6,20,.75);
  --rad:16px; --rad-sm:11px;
  --font:-apple-system,'Inter','Segoe UI',Roboto,system-ui,sans-serif;
  --mono:ui-monospace,'SF Mono','JetBrains Mono',Menlo,Consolas,monospace;
}
[data-theme="light"]{
  --bg:#eef2fa; --bg-2:#e7edf8; --panel:#ffffff; --panel-2:#f5f8ff; --inset:#eef3fc;
  --line:#dbe4f3; --line-2:#c3d1ea; --text:#131b30; --text-2:#3f4f74; --mut:#6c7ea8;
  --field:#ffffff;
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
.hero{margin:0 0 20px}
h1{font-size:clamp(30px,5vw,46px);line-height:1.03;margin:0;font-weight:800;letter-spacing:-.025em}
.hero-sub{margin:4px 0 0;font-size:14.5px;color:var(--mut);font-weight:500}
.hero-meta{margin:15px 0 0;padding-top:14px;border-top:1px solid var(--line);
  font-size:12.5px;color:var(--mut)}
/* ---- detail-nav card (index links to plot files) — note: #detail-nav is an ID ---- */
#detail-nav>p{margin:0;padding:13px 20px;border-bottom:1px solid var(--line);
  font-size:12.5px;line-height:1.6;color:var(--mut)}
.detail-list{list-style:none;margin:0;padding:6px 0 10px}
.detail-list li{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap;
  padding:10px 20px;border-bottom:1px solid var(--line);font-size:13px}
.detail-list li:last-child{border-bottom:none}
.detail-list a{text-decoration:none;color:var(--text)}
.detail-list a:hover{text-decoration:underline}
.rn-sub{font-size:11.5px;color:var(--mut)}
.chip{display:inline-flex;align-items:center;gap:7px;padding:6px 12px;border:1px solid var(--line);
  border-radius:999px;background:rgba(255,255,255,.03);font-size:12.5px;color:var(--text-2)}
.chip b{color:var(--text);font-weight:700}
.bigstat{min-width:104px;border:1px solid var(--line);border-radius:var(--rad-sm);
  background:rgba(255,255,255,.03);padding:14px}
.bigstat .n{font-size:30px;font-weight:800;letter-spacing:-.02em;line-height:1}
.bigstat .n small{font-size:15px;color:var(--text-2);font-weight:600}
.bigstat .n .ok{color:var(--teal)}
.bigstat .n .no{color:var(--red)}
.bigstat .l{font-size:11.5px;color:var(--mut);margin-top:7px;letter-spacing:.01em}
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

/* ---- interactive UMAP (gene dropdown + view) ---- */
.umap-bar{display:flex;align-items:center;gap:18px;flex-wrap:wrap;padding:6px 16px 4px}
.umap-bar label{display:inline-flex;align-items:center;gap:8px;font-size:12px;font-weight:600;
  color:var(--text-2);text-transform:none}
.umap-sel{appearance:none;-webkit-appearance:none;background:var(--field);color:var(--text);
  border:1px solid var(--line-2);border-radius:9px;padding:7px 30px 7px 11px;font-size:12.5px;font-family:inherit;
  cursor:pointer;background-image:linear-gradient(45deg,transparent 50%,var(--mut) 50%),linear-gradient(135deg,var(--mut) 50%,transparent 50%);
  background-position:calc(100% - 14px) 15px,calc(100% - 9px) 15px;background-size:5px 5px,5px 5px;background-repeat:no-repeat}
.umap-sel:focus{outline:2px solid color-mix(in srgb,var(--teal) 55%,transparent);border-color:var(--teal)}
.umap-status{margin-left:auto;font-size:11.5px;color:var(--mut);font-family:var(--mono)}
.umap-host{width:100%}
.umap-host .js-plotly-plot{margin:2px 4px 0}

/* ---- footer ---- */
.foot{margin-top:42px;padding-top:18px;border-top:1px solid var(--line);color:var(--mut);
  font-size:12px;display:flex;justify-content:space-between;gap:14px;flex-wrap:wrap}
.foot code{font-family:var(--mono);color:var(--text-2)}
"""

JS = """
/* ---- interactive UMAP: ONE canvas per sample, re-coloured client-side ---- */
var GCR_UMAP_CFG = { displaylogo: false, scrollZoom: false, showTips: false,
  modeBarButtonsToRemove: ['toImage', 'hoverCompareCartesian'] };
/* box figures: also drop the "Autoscale" button. In this plotly build the
   Autoscale action recomputes box-trace auto-range (broken: first point +/-
   one order of magnitude), and its two-step relayout cannot be safely
   intercepted. Double-click (redirected to the data-derived range) and the
   "Reset axes" button cover the reset use case instead. */
var GCR_BOXCFG = { displaylogo: false, scrollZoom: false, showTips: false,
  modeBarButtonsToRemove: ['toImage', 'hoverCompareCartesian', 'autoScale2d'] };
function gcrUmapText() { return '#182338'; }
function gcrUmapAxis() { return '#4a5c82'; }
function gcrUmapGrid() { return 'rgba(100,116,139,0.14)'; }
function gcrUmapBaseLayout(legend) {
  return {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: gcrUmapText(), size: 11, family: "-apple-system,Inter,system-ui,sans-serif" },
    margin: { l: 46, r: 20, t: 38, b: 34 }, height: 452, autosize: true,
    showlegend: legend,
    legend: (legend ? { orientation: 'h', yanchor: 'bottom', y: -0.32, xanchor: 'left', x: 0,
      font: { size: 11 }, bgcolor: 'rgba(0,0,0,0)' } : {}),
    hoverlabel: { bgcolor: '#0c1526', bordercolor: '#2e4067',
      font: { color: '#e8eefb', size: 12, family: "-apple-system,Inter,system-ui,sans-serif" } },
    xaxis: { gridcolor: gcrUmapGrid(), linecolor: '#2e4067', zeroline: false, showline: true,
      ticks: 'outside', ticklen: 5, tickwidth: 1,
      title: { text: 'UMAP-1', font: { color: gcrUmapAxis(), size: 11 } },
      tickfont: { color: gcrUmapAxis(), size: 10 } },
    yaxis: { gridcolor: gcrUmapGrid(), linecolor: '#2e4067', zeroline: false, showline: true,
      ticks: 'outside', ticklen: 5, tickwidth: 1,
      title: { text: 'UMAP-2', font: { color: gcrUmapAxis(), size: 11 } },
      tickfont: { color: gcrUmapAxis(), size: 10 } }
  };
}
function gcrUmapFigure(d) {
  var gs = document.getElementById('umapgene-' + d.uid);
  var gene = gs ? gs.value : (d.genes[0] || '');
  var xs = d.xs, ys = d.ys, n = xs.length,
      umi = (d.umap_umis && d.umap_umis[gene]) ? d.umap_umis[gene] : (d.umis[gene] || []);
  var cols = new Float64Array(n), cmax = 0, j, v;
  for (j = 0; j < n; j++) { v = Math.log10((umi[j] || 0) + 1); cols[j] = v; if (v > cmax) cmax = v; }
  if (cmax <= 0) cmax = 1.0;
  var exprData = [{ type: 'scatter', mode: 'markers', x: xs, y: ys, name: gene,
    marker: { size: 5, opacity: 0.82, color: cols, colorscale: d.cscale, cmin: 0, cmax: cmax, showscale: true,
      colorbar: { title: { text: gene + '  log10 UMI', side: 'right', font: { size: 10, color: gcrUmapAxis() } },
        len: 0.9, thickness: 11, outlinewidth: 1, outlinecolor: '#22314f',
        tickfont: { color: gcrUmapAxis(), size: 10 }, tickformat: '.1f' } },
    customdata: umi,
    hovertemplate: 'x %{x:.2f}<br>y %{y:.2f}<br><b>UMI %{customdata}</b><extra>single cell</extra>' }];
  var exprLayout = gcrUmapBaseLayout(false);
  exprLayout.title = { text: gene + ' \u00b7 expression on UMAP',
    font: { color: gcrUmapText(), size: 14, family: "-apple-system,Inter,system-ui,sans-serif" } };
  return { data: exprData, layout: exprLayout, meta: n + ' cells' };
}
/* ---- shared light-theme base layout for the histogram / box / bars / dots ---- */
function gcrAxisBase() {
  return { gridcolor: gcrUmapGrid(), linecolor: '#c3d1ea', showline: true, zeroline: false,
    ticks: 'outside', ticklen: 5, tickwidth: 1,
    tickfont: { color: gcrUmapAxis(), size: 10 },
    title: { font: { color: gcrUmapAxis(), size: 11 } } };
}
function gcrCartBase(legend) {
  var layout = {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: gcrUmapText(), size: 11, family: '-apple-system,Inter,system-ui,sans-serif' },
    margin: { l: 46, r: 20, t: 30, b: 34 }, autosize: true, showlegend: !!legend,
    legend: { orientation: 'h', yanchor: 'bottom', y: -0.28, xanchor: 'left', x: 0,
      font: { size: 11 }, bgcolor: 'rgba(0,0,0,0)' } };
  /* xaxis and yaxis MUST be separate objects. A shared reference makes a
     "layout.yaxis.type = ..." leak into xaxis and makes each axis-title
     overwrite the other (all axes ended up with the last-set title and
     the box plot got a log-scale x axis of 1e-30..1e30). */
  layout.xaxis = gcrAxisBase(); layout.yaxis = gcrAxisBase();
  return layout;
}
function gcrTitle(text) {
  return { text: text, font: { color: gcrUmapText(), size: 14, family: '-apple-system,Inter,system-ui,sans-serif' } };
}
function gcrSetFig(host, fig) {
  if (!host) return;
  if (fig.none) {
    host.innerHTML = '<div class="empty">' + (fig.msg || 'not available') + '</div>';
    return;
  }
  /* remember the data-derived range so the "Autoscale" modebar button can
     never fall back to plotly's auto-range (box traces compute a broken
     auto-range in this plotly build: first data point +/- one order of
     magnitude, which is what made the axes jump and stretch on reset) */
  var rl = null;
  if (fig.layout) {
    if ((fig.layout.yaxis && fig.layout.yaxis.range) || (fig.layout.xaxis && fig.layout.xaxis.range)) {
      rl = {};
      if (fig.layout.yaxis && fig.layout.yaxis.range) {
        rl['yaxis.range'] = fig.layout.yaxis.range.slice(); rl['yaxis.autorange'] = false;
      }
      if (fig.layout.xaxis && fig.layout.xaxis.range) {
        rl['xaxis.range'] = fig.layout.xaxis.range.slice(); rl['xaxis.autorange'] = false;
      }
    }
  }
  host._gcrRestore = rl;
  if (rl && !host._gcrBound) {
    host._gcrBound = true;
    try {
      /* this plotly build has no plotly_autoscale event: the Autoscale
         button (and dblclick) arrive as a relayout with autorange:true,
         which for box traces recomputes the broken auto-range (first
         point +/- one order of magnitude).  Redirect any such request
         back to the data-derived range, so every reset lands on the
         same correct range instead of stretching. */
      host.on('plotly_relayout', function (data) {
        if (!host._gcrRestore || !data) return;
        var wantAutorange = (data['yaxis'] && data['yaxis'].autorange === true) ||
                            (data['xaxis'] && data['xaxis'].autorange === true);
        if (wantAutorange) window.Plotly.relayout(host, host._gcrRestore);
      });
    } catch (e2) {}
  }
  try { window.Plotly.react(host, fig.data, fig.layout, fig.cfg || GCR_UMAP_CFG); } catch (e) {}
}
/* UMI histogram of expressing cells, mirroring the old static figure:
   log10(UMI+1) binned to 28 bars, bar labels = the UMI at each bin centre. */
function gcrHistFig(d, gene) {
  if (d.norow && d.norow[gene]) {
    return { none: true, msg: 'this gene is not present in this sample\u2019s reference (NO-ROW), so there are no cells to show', meta: '' };
  }
  var umi = d.umis[gene] || [], t = [], tmin = Infinity, tmax = -Infinity, i, j;
  for (i = 0; i < umi.length; i++) {
    var v = umi[i] | 0;
    if (v >= 1) {
      var tv = Math.log10(v + 1);
      t.push(tv);
      if (tv < tmin) tmin = tv;
      if (tv > tmax) tmax = tv;
    }
  }
  var layout = gcrCartBase(false);
  layout.title = gcrTitle(gene + ' \u00b7 UMI distribution of expressing cells');
  layout.xaxis.title.text = 'UMI per cell';
  layout.yaxis.title.text = 'cells';
  layout.bargap = 0.05;
  if (!t.length) return { none: true, msg: 'no expressing cells for ' + gene, meta: '0 expressing cells' };
  if (tmax - tmin < 1e-9) tmax = tmin + 1e-9;
  var nb = 28, span = tmax - tmin, counts = [];
  for (j = 0; j < nb; j++) counts.push(0);
  for (j = 0; j < t.length; j++) {
    counts[Math.min(nb - 1, Math.max(0, Math.floor((t[j] - tmin) / span * nb)))] += 1;
  }
  var labels = [];
  for (j = 0; j < nb; j++) {
    labels.push(String(Math.max(0, Math.round(10 ** (tmin + (j + 0.5) / nb * span) - 1))));
  }
  layout.xaxis.type = 'category';
  layout.xaxis.tickvals = labels.filter(function (_v, i2) { return i2 % Math.max(1, nb / 6 | 0) === 0; });
  var cmax = 0;
  for (j = 0; j < counts.length; j++) if (counts[j] > cmax) cmax = counts[j];
  /* explicit range so a gene switch never inherits the previous gene's
     auto-range (plotly caches it across Plotly.react) */
  if (cmax > 0) { layout.yaxis.range = [0, cmax * 1.12]; layout.yaxis.autorange = false; }
  return { data: [{ type: 'bar', x: labels, y: counts,
    marker: { color: '#38bdf8', line: { width: 0 } },
    hovertemplate: '\u2248 %{x} UMI<br>%{y} cells<extra>expressing</extra>' }],
    layout: layout, meta: t.length + ' expressing cells' };
}
/* total UMI per cell, box split by the selected gene's detection status */
function gcrBoxFig(d, gene) {
  if (d.norow && d.norow[gene]) {
    return { none: true, msg: 'this gene is not present in this sample\u2019s reference (NO-ROW), so there are no cells to show', meta: '' };
  }
  var umi = d.umis[gene] || [], th = d.threshold, tot = d.totals || [];
  var pos = [], neg = [], npos0 = 0, nneg0 = 0, i, tv, vmin = Infinity, vmax = -Infinity;
  for (i = 0; i < tot.length; i++) {
    tv = tot[i];
    if ((umi[i] | 0) >= th) {
      if (tv > 0) { pos.push(tv); if (tv < vmin) vmin = tv; if (tv > vmax) vmax = tv; } else { npos0++; }
    } else {
      if (tv > 0) { neg.push(tv); if (tv < vmin) vmin = tv; if (tv > vmax) vmax = tv; } else { nneg0++; }
    }
  }
  var negName = 'not detected', posName = gene + '  detected';
  var layout = gcrCartBase(true);
  layout.title = gcrTitle('Total UMI per cell  \u00b7  detected vs not detected');
  layout.xaxis.type = 'category';
  layout.xaxis.categoryorder = 'array';
  layout.xaxis.categoryarray = [negName, posName];
  layout.xaxis.title.text = 'sample group';
  layout.yaxis.type = 'log';
  layout.yaxis.title.text = 'total UMI per cell (log)';
  /* Set an explicit range on every render. Plotly keeps the previous
     auto-range across Plotly.react, so after switching genes the axis kept
     the old gene's range and a double-click "reset" jumped to yet another
     range (the "reset only gets longer" bug). With an explicit range and
     autorange left true, zoom + double-click/modebar-reset returns to
     exactly this range. */
  if (vmin !== Infinity) {
    var lmin = Math.log10(vmin), lmax = Math.log10(vmax), pad = Math.max(0.12, (lmax - lmin) * 0.06);
    layout.yaxis.range = [lmin - pad, lmax + pad];
    layout.yaxis.autorange = false;
  }
  /* This plotly build breaks a box trace's statistics when x is shorter than
     y (a single category string makes plotly pair x[0] against y and the box
     collapses to a flat line). Repeat the category once per y value instead. */
  function fillX(nm, ary) { var x = new Array(ary.length), k; for (k = 0; k < ary.length; k++) x[k] = nm; return x; }
  var data = [];
  if (neg.length) data.push({ type: 'box', x: fillX(negName, neg), y: neg, name: negName,
    marker: { color: '#3a4a66', size: 3, opacity: 0.5 }, line: { width: 1.1 },
    fillcolor: 'rgba(58,74,102,0.18)' });
  if (pos.length) data.push({ type: 'box', x: fillX(posName, pos), y: pos, name: posName,
    marker: { color: '#2dd4bf', size: 3, opacity: 0.6 }, line: { width: 1.2 },
    fillcolor: 'rgba(45,212,191,0.16)' });
  /* Invisible 2-point anchor spanning the true data min/max. In this plotly
     build a box trace's auto-range is computed from its first point only
     (+/- one order of magnitude), so every dblclick / "Autoscale" reset
     landed on a nonsense axis (the "reset only gets longer" bug). A proper
     scatter anchor forces every auto-range pass to span the real data. */
  if (vmin !== Infinity) data.push({ type: 'scatter', mode: 'markers',
    x: [negName, negName], y: [vmin, vmax], marker: { size: 1, opacity: 0 },
    line: { width: 0 }, hoverinfo: 'skip', showlegend: false });
  /* keep a labelled slot for an empty group instead of silently dropping it,
     so a 100%-detected (or 100%-undetected) gene still explains where its
     "not detected" (or "detected") box went */
  var ann = [];
  if (!neg.length) ann.push({ x: negName, xref: 'x', yref: 'paper', y: 0.45,
    text: '0 cells', showarrow: false, font: { color: gcrUmapAxis(), size: 12 } });
  if (!pos.length) ann.push({ x: posName, xref: 'x', yref: 'paper', y: 0.45,
    text: '0 cells', showarrow: false, font: { color: gcrUmapAxis(), size: 12 } });
  if (ann.length) layout.annotations = ann;
  if (!data.length) return { none: true, msg: 'no cells to plot', meta: '' };
  var meta = pos.length + ' detected vs ' + neg.length + ' not detected cells';
  var hidden = npos0 + nneg0;
  if (hidden) meta += ' \u00b7 ' + hidden + ' cell(s) with 0 total UMI are not visible on a log axis';
  return { data: data, layout: layout, meta: meta, cfg: GCR_BOXCFG };
}
/* coverage track on the selected gene's own transgene contig: a dense depth
   array per gene, served in the payload (one `samtools depth -a` pass on its
   contig at report time). Follows the same Gene selector as the UMAP, the
   histogram and the box plot. */
function gcrCovFig(d, gene) {
  var cov = d.cov;
  if (!cov) {
    return { none: true,
      msg: 'no aligned read file (BAM) found next to this sample\u2019s outputs',
      meta: '' };
  }
  if (cov.no_contig) {
    return { none: true,
      msg: 'no contig matching any of the listed genes was found among the ' +
           (cov.contigs_n || 0) + ' contigs in the BAM',
      meta: '' };
  }
  var t = cov.genes[gene];
  if (!t) {
    return { none: true,
      msg: 'no contig named ' + gene + ' in this sample\u2019s BAM (of ' +
           (cov.contigs_n || 0) + ' contigs)',
      meta: '' };
  }
  if (!t.depth || !t.depth.length) {
    return { none: true,
      msg: 'contig ' + t.contig + ' is present in the BAM but has no aligned ' +
           'reads (the whole contig is at zero coverage)',
      meta: 'no reads on ' + t.contig };
  }
  var x = new Array(t.depth.length), i, dmax = 0;
  for (i = 0; i < t.depth.length; i++) {
    x[i] = t.pos0 + i;
    if (t.depth[i] > dmax) dmax = t.depth[i];
  }
  var layout = gcrCartBase(false);
  layout.title = gcrTitle(gene + ' \u00b7 coverage on contig ' + t.contig);
  layout.xaxis.title.text = 'position on contig (bp)';
  layout.yaxis.title.text = 'coverage (\u00d7)';
  /* explicit ranges on every render: plotly keeps the previous gene's
     cached auto-range across Plotly.react (same gene-switch bug the box
     and histogram work around) */
  if (t.depth.length) {
    layout.xaxis.range = [t.pos0, t.pos0 + t.depth.length - 1];
    layout.xaxis.autorange = false;
  }
  if (dmax > 0) { layout.yaxis.range = [0, dmax * 1.12]; layout.yaxis.autorange = false; }
  var r = cov.reads ? cov.reads[gene] : null;
  var meta = t.contig + ' \u00b7 ' + t.length + ' bp';
  if (r != null) meta += ' \u00b7 ' + r + ' reads (called cells)';
  return { data: [{ type: 'scatter', mode: 'lines', x: x, y: t.depth, name: gene,
    line: { color: '#818cf8', width: 1.2 },
    fill: 'tozeroy', fillcolor: 'rgba(99,102,241,0.34)',
    hovertemplate: 'position %{x} bp \u00b7 coverage %{y}\u00d7<extra>' + gene + '</extra>' }],
    layout: layout, meta: meta };
}
/* across-samples: % cells detected bars + mean UMI dots for the picked gene */
function gcrBarsFig(a, gene) {
  var vals = (a.det[gene] || []).map(Number);
  var layout = gcrCartBase(false);
  layout.title = gcrTitle('% cells detected \u00b7 ' + gene);
  layout.xaxis.type = 'category';
  layout.xaxis.title.text = 'sample';
  layout.yaxis.title.text = '% cells (\u2265 detection threshold)';
  layout.bargap = 0.34;
  var vmax = 0;
  for (var i = 0; i < vals.length; i++) if (vals[i] > vmax) vmax = vals[i];
  /* explicit range (25% headroom for the outside bar labels) so switching
     genes never inherits the previous gene's cached auto-range */
  if (vmax > 0) { layout.yaxis.range = [0, vmax * 1.25]; layout.yaxis.autorange = false; }
  return { data: [{ type: 'bar', x: a.labels, y: vals,
    marker: { color: '#2dd4bf', line: { width: 0 } },
    text: vals.map(function (v) { return v.toFixed(1) + '%'; }),
    textposition: 'outside', textfont: { color: '#4a5c82', size: 11 },
    hovertemplate: '%{x}<br>%{y:.1f}% cells detected<extra></extra>' }],
    layout: layout, meta: vals.length + ' samples' };
}
function gcrDotsFig(a, gene) {
  var raw = a.mean[gene] || [];
  /* null = a sample with no expressing cells for this gene: shows as a blank
     label (no dot, no "0.00"); the axis range below ignores nulls */
  var vals = raw.map(function (v) { return (v == null) ? null : Number(v); });
  var layout = gcrCartBase(false);
  layout.title = gcrTitle('Mean UMI per expressing cell \u00b7 ' + gene);
  layout.xaxis.title.text = 'mean UMI per expressing cell';
  layout.yaxis.type = 'category';
  layout.yaxis.categoryorder = 'array';
  layout.yaxis.categoryarray = a.labels.slice();
  layout.yaxis.autorange = 'reversed';
  var nv = vals.filter(function (v) { return v !== null; });
  /* explicit range (headroom for the value labels) so a gene switch never
     inherits the previous gene's cached auto-range */
  if (nv.length) {
    var xmin = Math.min.apply(null, nv), xmax = Math.max.apply(null, nv);
    var pad = Math.max(1, (xmax - xmin) * 0.15);
    layout.xaxis.range = [Math.max(0, xmin - pad), xmax + pad];
    layout.xaxis.autorange = false;
  }
  return { data: [{ type: 'scatter', x: vals, y: a.labels.slice(), mode: 'markers+text',
    marker: { size: 11, color: '#f472b6', line: { color: '#0c1526', width: 1.5 } },
    text: vals.map(function (v) { return v ? Number(v).toFixed(2) : ''; }),
    textposition: 'middle right',
    textfont: { color: '#4a5c82', size: 11 },
    hovertemplate: '%{y}<br>mean %{x:.2f} UMI<extra>expressing cells</extra>' }],
    layout: layout, meta: 'expressing cells only (UMI \u2265 1)' };
}
function acrossRender() {
  var a = window.__gcrAcross;
  if (!a || !window.Plotly) return;
  var sel = document.getElementById('gcracrosssel');
  var gene = (sel && sel.value) || (a.genes[0] || '');
  var bf = gcrBarsFig(a, gene);
  var bm = document.getElementById('gcracrossdm');
  if (bm) bm.textContent = bf.meta;
  gcrSetFig(document.getElementById('gcracrossbars'), bf);
  var df = gcrDotsFig(a, gene);
  var dm = document.getElementById('gcracrossmm');
  if (dm) dm.textContent = df.meta;
  gcrSetFig(document.getElementById('gcracrossdots'), df);
  /* stat cards follow the same gene selector (per-gene values in the payload) */
  var ev = a.expr ? (a.expr[gene] != null ? a.expr[gene] : 0) : 0;
  var e1 = document.getElementById('gcracrossexpr');
  if (e1) e1.textContent = ev + '/' + a.labels.length;
  var e1l = document.getElementById('gcracrossexpl');
  if (e1l) e1l.textContent = 'samples where ' + gene + ' is expressed';
  var dv = (a.det[gene] || []).map(Number);
  var dvsum = 0;
  for (var di = 0; di < dv.length; di++) dvsum += dv[di];
  var e2 = document.getElementById('gcracrossavg');
  if (e2) e2.textContent = ((dv.length ? dvsum / dv.length : 0).toFixed(1)) + '%';
  var mvv = [];
  var msrc = a.mean[gene] || [];
  for (var mi = 0; mi < msrc.length; mi++) {
    var mv = msrc[mi];
    if (mv != null && !isNaN(Number(mv))) mvv.push(Number(mv));
  }
  var e3 = document.getElementById('gcracrossmx');
  if (e3) e3.textContent = (mvv.length ? Math.max.apply(null, mvv).toFixed(2) : '0.00');
}
/* one UMAP + one histogram + one box plot per sample, all re-rendered in place
   by the gene <select> (umapSwitch). One canvas per plot type, always. */
function umapRender(uid) {
  var d = window.__umap && window.__umap[uid];
  if (!d || !window.Plotly) return;
  d.uid = uid;
  var sel = document.getElementById('umapgene-' + uid);
  var gene = (sel && sel.value) || (d.genes[0] || '');
  if (d.umap) {
    var f = gcrUmapFigure(d);
    var st = document.getElementById('umapstatus-' + uid);
    if (st) st.textContent = f.meta;
    gcrSetFig(document.getElementById('umaphost-' + uid), f);
  }
  var h = gcrHistFig(d, gene);
  var hm = document.getElementById('umhistmeta-' + uid);
  if (hm) hm.textContent = h.meta;
  gcrSetFig(document.getElementById('umhist-' + uid), h);
  var b = gcrBoxFig(d, gene);
  var bm = document.getElementById('umboxmeta-' + uid);
  if (bm) bm.textContent = b.meta;
  gcrSetFig(document.getElementById('umbox-' + uid), b);
  var c = gcrCovFig(d, gene);
  var cm = document.getElementById('umcovmeta-' + uid);
  if (cm) cm.textContent = c.meta;
  gcrSetFig(document.getElementById('umcov-' + uid), c);
}
function umapSwitch(uid) { umapRender(uid); }
/* the per-sample inline <script> tags in the body (parsed before this) put each
   payload into window.__umap and its uid in window.__umapInit; now that
   umapRender / umapSwitch / acrossRender are defined, render the registered ones. */
if (window.__umapInit) {
  for (var _i = 0; _i < window.__umapInit.length; _i++) { umapRender(window.__umapInit[_i]); }
}
if (window.__gcrAcross) { acrossRender(); }
"""

# Palette for the charts
_EXPR_SCALE = ["#14162e", "#3b2a76", "#7c3aed", "#c026d3", "#f07bb0", "#facc15"]


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
# Plotly figure builders (fixed light-theme layout)
# --------------------------------------------------------------------------- #
def _plotly_theme(fig, title, x_title="", y_title="", height=360, legend=False,
                  legend_title=""):
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#182338", size=11, family="-apple-system,Inter,system-ui,sans-serif"),
        margin=dict(l=46, r=20, t=18, b=34), height=height, autosize=True,
        showlegend=legend,
        legend=dict(orientation="h", yanchor="bottom", y=-0.24, xanchor="left", x=0,
                    font=dict(size=11), bgcolor="rgba(0,0,0,0)"),
        hoverlabel=dict(bgcolor="#0c1526", bordercolor="#2e4067",
                        font=dict(color="#e8eefb", size=12,
                                 family="-apple-system,Inter,sans-serif")),
    )
    fig.update_xaxes(gridcolor="rgba(40,70,130,0.16)", linecolor="#c3d1ea",
                     zeroline=False, ticks="outside", ticklen=5, tickwidth=1,
                     title_text=x_title, title_font=dict(color="#4a5c82", size=11),
                     tickfont=dict(color="#4a5c82", size=10))
    fig.update_yaxes(gridcolor="rgba(40,70,130,0.16)", linecolor="#c3d1ea",
                     zeroline=False, ticks="outside", ticklen=5, tickwidth=1,
                     title_text=y_title, title_font=dict(color="#4a5c82", size=11),
                     tickfont=dict(color="#4a5c82", size=10))
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


# unique, ever-increasing id so each sample's UMAP controls get a distinct DOM id
_UMAP_SEQ = [0]


def _cov_seed_meta(games, cov):
    """Server-rendered placeholder for the coverage meta line, using the first
    gene (the initial dropdown value). The JS replaces it on load and on every
    Gene-selector change (gcrCovFig keeps the same wording)."""
    if cov is None:
        return "no aligned read file (BAM)"
    tracks = cov.get("genes") or {}
    if cov.get("no_contig") or not tracks:
        return "no matching contig ({} contigs in BAM)".format(
            cov.get("contigs_n") or 0)
    for g in games:
        t = tracks.get(g)
        if t is None:
            continue
        if not t.get("depth"):
            return "no reads on " + t["contig"]
        meta = g + " \u00b7 contig " + t["contig"] + " \u00b7 {} bp".format(t["length"])
        n = (cov.get("reads") or {}).get(g)
        if n is not None:
            meta += " \u00b7 {:,} reads (called cells)".format(n)
        return meta
    return "no matching contig ({} contigs in BAM)".format(
        cov.get("contigs_n") or 0)


def render_sample_interactive(sample_name, coords, barcodes, bc_total, bc_gene,
                              genes, threshold, norow=None, cov=None):
    """One interactive plot set per sample (UMAP + UMI histogram + total-UMI box
    + coverage track on the gene's transgene contig).

    In place of one static figure per gene per plot type: one shared payload per
    sample (the UMAP coordinate frame, one UMI vector per gene, per-cell total
    UMI) is embedded ONCE as JSON; the gene <select> drives the client-side JS
    (umapRender / gcrHistFig / gcrBoxFig in the JS constant) which re-renders one
    canvas per plot type in place. So the DOM holds exactly one figure of each
    type no matter how many genes are inspected — what exhausts browser memory on
    big multi-gene runs — and the report stays fully offline.

    The histogram, box plot and the "· N detected" dropdown labels are computed
    from ALL barcodes of the sample (the same population the stats table uses).
    When the UMAP only covers a subset of cells, that coordinate subset is used
    for the UMAP trace alone, and the panel says so (plus a stderr note).

    cov (the sample_bam_summary output, or None without a BAM): the per-gene
    coverage tracks for this gene batch are embedded in the payload and the
    coverage panel follows the same Gene selector (gcrCovFig on the JS side).
    """
    have_umap = len(coords or {}) >= 2
    order = [bc for bc in (barcodes or []) if bc]
    if have_umap:
        uorder = [bc for bc in order if bc in coords]
        cx = [coords[bc][0] for bc in uorder]
        cy = [coords[bc][1] for bc in uorder]
        if len(uorder) != len(order):
            sys.stderr.write(
                "note: %s: UMAP covers %d of %d cells; histogram and box plot "
                "use all %d cells\n" % (sample_name, len(uorder), len(order),
                                       len(order)))
    else:
        uorder = []
        cx, cy = [], []
    # one aligned UMI vector per gene (ints), keyed by the (string) gene name,
    # plus the per-cell total UMI in the same (full) barcode order — the box and
    # the histogram are split from these two vectors on the client side
    umis = {}
    for g in genes:
        key = str(g)
        umis[key] = [int(bc_gene.get(bc, {}).get(g, 0)) for bc in order]
    totals = [int(bc_total.get(bc, 0)) for bc in order]
    # UMAP-trace UMI vectors, aligned to the (possibly smaller) UMAP subset
    umap_umis = {}
    if have_umap:
        for g in genes:
            umap_umis[str(g)] = [int(bc_gene.get(bc, {}).get(g, 0))
                                 for bc in uorder]

    _UMAP_SEQ[0] += 1
    uid = _UMAP_SEQ[0]

    umap_subset_note = ""
    if have_umap and len(uorder) != len(order):
        umap_subset_note = ('<div class="note">UMAP covers %d of %d cells; '
                            'histogram and box plot use all %d cells.</div>'
                            % (len(uorder), len(order), len(order)))

    opts = []
    for g in genes:
        key = str(g)
        det = sum(1 for v in umis[key] if v >= threshold)
        opts.append('<option value="%s">%s · %d detected</option>'
                    % (html.escape(key, quote=True), html.escape(key), det))

    ncells = len(order)
    # ensure_ascii keeps the embedded JSON plain-ASCII (gene names are
    # identifiers), so it can never contain the '</script>' terminator.
    pl = {"umap": have_umap, "genes": [str(g) for g in genes], "umis": umis,
          "totals": totals, "threshold": int(threshold),
          "norow": {str(g): 1 for g in (norow or {})}}
    if have_umap:
        pl["xs"], pl["ys"], pl["cscale"] = cx, cy, list(_EXPR_SCALE)
        pl["n_umap"] = len(cx)
        pl["umap_umis"] = umap_umis
    if cov is not None:
        tracks = cov.get("genes") or {}
        pl["cov"] = {
            "no_contig": bool(cov.get("no_contig")),
            "contigs_n": cov.get("contigs_n") or 0,
            "reads": {str(g): (cov.get("reads") or {}).get(str(g)) for g in genes},
            "genes": {str(g): tracks[str(g)] for g in genes if str(g) in tracks},
        }
    # seed the coverage meta line server-side (same wording gcrCovFig uses);
    # the JS overwrites it on load and on every gene switch
    cov_seed = _cov_seed_meta(genes, cov)
    payload = json.dumps(pl, separators=(",", ":"), ensure_ascii=True)

    if have_umap:
        umap_host = ('<div id="umaphost-' + str(uid) + '" class="umap-host"></div>'
                     + umap_subset_note)
    else:
        umap_host = ('<div class="empty">No interactive UMAP for '
                     + html.escape(sample_name)
                     + ' (needs a UMAP_1 / UMAP_2 coordinate file next to the '
                     'matrix — the UMI plots below still work).</div>')
    return (
        '<div class="panel umap-panel" style="margin-bottom:16px">'
        '<div class="panel-h"><span class="t"><span class="sw"></span>'
        'Per-gene plots · pick a gene</span>'
        '<span class="m">' + format(ncells, ",") + ' cells</span></div>'
        '<div class="umap-bar">'
        '<label>Gene&nbsp;<select class="umap-sel" id="umapgene-' + str(uid) + '" '
        'onchange="umapSwitch(' + str(uid) + ')">' + "".join(opts) + '</select></label>'
        '<span class="umap-status" id="umapstatus-' + str(uid) + '"></span>'
        '</div>'
        + umap_host +
        '<div class="note">The gene selector above re-renders all four plots in place (one of '
        'each type, never one per gene): the UMAP (dark to bright = low to high UMI, '
        'log10 scale), the histogram of expressing cells, the total-UMI box split by '
        'detection status, and the coverage track on the gene\u2019s transgene contig. '
        'Drag to zoom, shift-drag to pan, double-click to reset; hover a point for '
        'its UMI.</div>'
        '<script type="application/json" id="umapdata-' + str(uid) + '">' + payload + '</script>'
        '<script>(function(){'
        'var el=document.getElementById("umapdata-' + str(uid) + '");if(!el)return;'
        'var d=JSON.parse(el.textContent);'
        '(window.__umap=window.__umap||{})[' + str(uid) + ']=d;'
        '(window.__umapInit=window.__umapInit||[]).push(' + str(uid) + ');'
        '})();</script>'
        '</div>'
        '<div class="grid2" style="margin-bottom:16px;align-items:start">'
        '<div class="panel">'
        '<div class="panel-h"><span class="t"><span class="sw"></span>'
        'UMI distribution of expressing cells</span>'
        '<span class="m" id="umhistmeta-' + str(uid) + '"></span></div>'
        '<div id="umhist-' + str(uid) + '"></div>'
        '<div class="note">Expressing cells = cells with at least 1 UMI of the selected '
        'gene (bins are logarithmic). Follows the Gene selector above.</div>'
        '</div>'
        '<div class="panel">'
        '<div class="panel-h"><span class="t"><span class="sw"></span>'
        'Total UMI per cell · detected vs not detected</span>'
        '<span class="m" id="umboxmeta-' + str(uid) + '"></span></div>'
        '<div id="umbox-' + str(uid) + '"></div>'
        '<div class="note">Total UMI of the whole transcriptome per cell. Left box: cells '
        'where the selected gene is not detected; right box: cells where it is detected '
        '(log scale). Follows the Gene selector above.</div>'
        '</div>'
        '</div>'
        '<div class="panel" style="margin-bottom:16px">'
        '<div class="panel-h"><span class="t"><span class="sw"></span>'
        'Coverage on the transgene contig</span>'
        '<span class="m" id="umcovmeta-' + str(uid) + '">' + html.escape(cov_seed) + '</span></div>'
        '<div id="umcov-' + str(uid) + '"></div>'
        '<div class="note">Follows the Gene selector: aligned-read base coverage on the '
        "selected gene's own transgene contig (from the BAM, not the matrix). "
        "10x 3' data is expected to pile up toward the 3' end; coverage counts all "
        'barcodes and excludes duplicate-flagged reads, so it will not reconcile '
        'with the "reads (called cells)" column.</div>'
        '</div>'
    )


def match_contig(genes, contigs):
    if not genes or not contigs:
        return None
    # Exact (case-insensitive) match only — deliberately NO substring fallback:
    # a gene like "Cre" must never be labelled with the "CreERT2" contig. A
    # gene with no exactly-named contig is treated as "no contig".
    exact = {c.lower(): c for c in contigs}
    for g in genes:
        gl = str(g).lower()
        if gl in exact:
            return (g, exact[gl])
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
    # per-gene across-sample values (bars = % cells detected, dots = mean UMI) —
    # small (samples x genes numbers), so every gene's panels can live in the DOM
    # and the Gene selector below re-renders them client-side
    det_all, mean_all = {}, {}
    for g in genes:
        det_all[str(g)] = []
        mean_all[str(g)] = []
    expr_all = {str(g): 0 for g in genes}
    for _sid, stats, _ct in results:
        for g in genes:
            d = stats.get(g) or {}
            cells = d.get("cells_total") or 0
            det_all[str(g)].append((100.0 * (d.get("det") or 0) / cells) if cells else 0.0)
            # samples with no expressing cells for this gene send None; the JS
            # shows a blank label, not a misleading "0.00"
            mean_all[str(g)].append(float(d["mean"]) if d.get("n_pos") else None)
            if d.get("status") == "EXPRESSED":
                expr_all[str(g)] += 1
    dets = det_all[str(first)]
    means = [m for m in mean_all[str(first)] if m is not None]
    avg = sum(dets) / len(dets) if dets else 0.0
    mx = max(means) if means else 0.0
    opts = []
    for g in genes:
        n_exp = sum(1 for _sid, st, _ct in results
                    if (st.get(g) or {}).get("status") == "EXPRESSED")
        opts.append('<option value="%s">%s · %d/%d expressed</option>'
                    % (html.escape(str(g), quote=True), html.escape(str(g)),
                       n_exp, len(labels)))
    # ensure_ascii keeps the embedded JSON plain-ASCII (gene names are
    # identifiers), so it can never contain the '</script>' terminator.
    payload = json.dumps(
        {"labels": labels, "genes": [str(g) for g in genes],
         "expr": expr_all, "det": det_all, "mean": mean_all},
        separators=(",", ":"), ensure_ascii=True)
    cards = ("<div class=grid3 style=margin-bottom:16px>"
             f"<div class=bigstat><div class=n id=gcracrossexpr>{expr_all[str(first)]}/{len(results)}</div>"
             f"<div class=l id=gcracrossexpl>samples where {html.escape(str(first))} is expressed</div></div>"
             f"<div class=bigstat><div class=n id=gcracrossavg>{avg:.1f}%</div>"
             "<div class=l>mean share of cells that are detected (transgene UMI &ge; threshold)</div></div>"
             f"<div class=bigstat><div class=n id=gcracrossmx>{mx:.2f}</div>"
             "<div class=l>max mean UMI per expressing cell</div></div></div>")
    return ('<div class=sec-h><h2><span class=ic>' + _p_icon("EXPRESSED") + "</span>Across samples</h2></div>"
            + cards
            + '<div class=umap-bar style="margin:0 0 14px 6px">'
            '<label>Gene&nbsp;<select class="umap-sel" id="gcracrosssel" onchange="acrossRender()">'
            + "".join(opts)
            + '</select></label>'
            '<span class=umap-status>both panels below follow this selector</span>'
            '</div>'
            '<div class=grid2>'
            '<div class=panel><div class=panel-h><span class=t><span class=sw></span>% cells detected</span>'
            '<span class=m id=gcracrossdm></span></div><div id=gcracrossbars></div>'
            '<div class=note>per sample, at the detection threshold</div></div>'
            '<div class=panel><div class=panel-h><span class=t><span class=sw></span>Mean UMI per expressing cell</span>'
            '<span class=m id=gcracrossmm></span></div><div id=gcracrossdots></div>'
            '<div class=note>expressing cells only (UMI \u2265 1)</div></div>'
            '</div>'
            + '<script type="application/json" id="gcracrossdata">' + payload + '</script>'
            + '<script>(function(){var el=document.getElementById("gcracrossdata");if(!el)return;'
              'window.__gcrAcross=JSON.parse(el.textContent);})();</script>')


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

        # one gene-driven plot set per sample (one UMAP, one hist of expressing
        # cells, one total-UMI box, one coverage track on the transgene
        # contig): the client-side Gene selector re-renders all four in place,
        # so the DOM keeps exactly one figure per plot type no matter how many
        # genes are inspected (never one figure per gene again).
        norow = {str(g) for g in genes
                 if (stats.get(g) or {}).get("row_present") is False}
        out.append(render_sample_interactive(sid, coords, p.get("barcodes") or [],
                                             bc_total, bc_gene,
                                             genes, threshold, norow,
                                             cov=p.get("cov") if coverage else None))
    return "\n".join(out)


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
    gene_list = [g.strip() for g in gene.split(",") if g.strip()]
    multi = len(gene_list) > 1

    # hero title block — one gene: the gene itself; many genes: a neutral
    # heading with the count (and the part range when this is a batch file)
    if multi:
        head_text = "Gene Expression"
        sub_parts = []
        if "\u00b7" in title:
            sub_parts.append(title.rsplit("\u00b7", 1)[1].strip())
        sub_parts.append("{} genes".format(len(gene_list)))
        hero_sub = " \u00b7 ".join(sub_parts)
    else:
        head_text = gene_main or title
        hero_sub = "expression report"

    meta_bits = []
    if gen_count:
        meta_bits.append(html.escape(gen_count) +
                         (" sample" if str(gen_count) == "1" else " samples"))
    if threshold:
        meta_bits.append("detection threshold " + html.escape(threshold))
    if generated:
        meta_bits.append("generated " + html.escape(generated))
    hero_meta = ("<div class=\"hero-meta\">" + " &middot; ".join(meta_bits) + "</div>"
                 if meta_bits else "")

    title_tag = html.escape(title)
    if gene_main and not multi:
        title_tag += " \u00b7 " + html.escape(gene_main)
    elif multi:
        title_tag += " \u00b7 " + "{} genes".format(len(gene_list))

    # plotly JS once — only when THIS file actually has interactive figures (a
    # pure-table index stays small); fallback to CDN if import fails so the file
    # still renders. The token check must scan BOTH sections and plots: in
    # split mode the "Across samples" charts live in sections while `plots`
    # can be empty for the index file.
    blob = "\n".join(sections) + plots
    if not any(k in blob for k in ("Plotly.newPlot", "umhist-",
                                   "umaphost-", "gcracrossbars")):
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

    return (
        "<!doctype html>\n<html lang=\"en\" data-theme=\"light\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<meta name=\"color-scheme\" content=\"light\">\n"
        "<title>" + title_tag + "</title>\n"
        "<style>" + CSS + "</style>\n"
        + js_block + "\n"
        "</head>\n<body>\n<div class=\"wrap\">\n"
        # hero
        "<header class=\"hero\">\n"
        "<h1>" + html.escape(head_text) + "</h1>\n"
        "<p class=\"hero-sub\">" + html.escape(hero_sub) + "</p>\n"
        + hero_meta + "\n"
        "</header>\n"
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
                         "all genes + links) plus N detail files (each with the tables for "
                         "its batch plus the interactive plots, incl. the per-gene "
                         "coverage tracks, this many per file). 0 = one combined file "
                         "(never split)")
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
        # per matched gene (on its own contig): one `samtools view -r` pass for
        # the read count in called cells, and one `samtools depth -a` pass for
        # the coverage track shown under the Gene selector
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
        # index: full tables (all genes) + links to the parts. The per-gene
        # coverage tracks live in the detail files (inside the interactive plot
        # set, driven by the Gene selector), so the index stays tables-only.
        sections_idx = list(sections)
        sections_idx.insert(1, render_detail_nav(parts, n_g))
        texts.append((index_path, render_html(title, meta, sections_idx, "")))
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
            plots_b = render_plots(plot_data, genes_b, args.min_umi, coverage=True)
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
