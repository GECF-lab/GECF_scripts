#!/usr/bin/env python3
"""prepare_for_pca.py — emit one PCA-ready CSV from the BRB-seq count matrix.

The STARsolo + rename_matrix output is an MTX triplet (rows = genes,
columns = samples), which most downstream tools (PCA in R / Excel / CLC)
want as a single CSV.  The legacy BRB pipeline produced this as
``<pool>_counts.detailed_forPCA.csv``; this script is the new pipeline's
equivalent.

Columns are the Sample_IDs from ``barcodes.tsv`` (already remapped by
rename_matrix.py — well labels, project codes, ... anything the manifest
used).  Row order and contents follow ``features.tsv`` exactly, so the
table lines up gene-for-gene with the matrix.

Usage:
    prepare_for_pca.py --samples-dir <OUT>/star/samples \
                       --out-csv   <pool>_counts_for_pca.csv \
                       [--dedup 1MM_Directional]

Exit codes: 0 ok; 1 missing/unreadable input; 2 structural mismatch.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

TRIPLE_RE = re.compile(
    r"^(\d+)\s+(\d+)\s+([+-]?(?:\d+\.?\d*|\.?\d+)(?:[eE][+-]?\d+)?)\s*$"
)


def _next_meaningful(lines, i):
    """Index of the next non-blank, non-comment line (or len(lines))."""
    while i < len(lines):
        s = lines[i].strip()
        if s and not s.startswith("%"):
            return i
        i += 1
    return i


def read_mtx(path: str):
    """Parse a MatrixMarket matrix coordinate file. Return (rows, cols, triples).

    ``triples`` is a list of (row0, col0, value_string), 0-based, in file order.
    Header comments before the %%MatrixMarket line are tolerated.  Raises
    ValueError on anything that is not a valid coordinate file.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    n = len(lines)
    i = 0
    while i < n:
        s = lines[i].strip().lower()
        if s.startswith("%%matrixmarket matrix"):
            break
        if s and not s.startswith("%"):
            raise ValueError(f"{path}: data before the %%MatrixMarket header")
        i += 1
    else:
        raise ValueError(f"{path}: no '%%MatrixMarket matrix' header")

    if not s.startswith("%%matrixmarket matrix coordinate"):
        raise ValueError(f"{path}: expected a 'matrix coordinate' file, got {s!r}")

    i = _next_meaningful(lines, i + 1)
    if i >= n:
        raise ValueError(f"{path}: no dimension line after header")
    fmt = lines[i].split()
    if len(fmt) < 3:
        raise ValueError(f"{path}: malformed dimension line {lines[i]!r}")
    rows, cols, nnz = int(fmt[0]), int(fmt[1]), int(fmt[2])
    if rows < 1 or cols < 1 or nnz < 0:
        raise ValueError(f"{path}: impossible dimensions {rows}x{cols}, nnz={nnz}")

    triples = []
    for ln in lines[i + 1:]:
        s = ln.strip()
        if not s or s.startswith("%"):
            continue
        m = TRIPLE_RE.match(s)
        if not m:
            raise ValueError(f"{path}: bad coordinate line {ln!r}")
        r, c = int(m.group(1)) - 1, int(m.group(2)) - 1
        if not (0 <= r < rows and 0 <= c < cols):
            raise ValueError(
                f"{path}: coordinate out of range {ln!r} (rows={rows} cols={cols})"
            )
        triples.append((r, c, m.group(3)))
    if len(triples) != nnz:
        raise ValueError(
            f"{path}: header declares nnz={nnz} but the file has {len(triples)} entries"
        )
    return rows, cols, triples


def read_barcodes(path: str) -> list[str]:
    with open(path, encoding="utf-8", errors="replace") as fh:
        names = [ln.strip() for ln in fh]
    names = [x for x in names if x]
    if not names:
        raise ValueError(f"{path}: no sample names")
    if len(names) != len(set(names)):
        raise ValueError(f"{path}: duplicate sample names: {sorted(set(n for n in names if names.count(n) > 1))[:5]} …")
    return names


def read_features(path: str):
    """features.tsv -> (n_rows, list of (gene_id, gene_symbol, feature_type))."""
    feats: list[tuple[str, str, str]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            s = ln.rstrip("\n").strip()
            if not s:
                continue
            parts = s.split("\t") if "\t" in s else s.split()
            gid = parts[0].strip()
            sym = parts[1].strip() if len(parts) > 1 else ""
            typ = parts[2].strip() if len(parts) > 2 else ""
            feats.append((gid, sym, typ))
    return len(feats), feats


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--samples-dir", required=True,
                   help="rename_matrix.py output dir (barcodes.tsv, features.tsv, umiDedup-*.mtx)")
    p.add_argument("--out-csv", required=True, help="output CSV path (gene rows x sample cols)")
    p.add_argument("--dedup", default="1MM_Directional",
                   help="which umiDedup-<method>.mtx to source (default: 1MM_Directional, "
                        "the pipeline's dedup count mode)")
    p.add_argument("--all-features", action="store_true",
                   help="include non-'Gene Expression' feature rows "
                        "(default: only 'Gene Expression' rows, matching a "
                        "--soloFeatures Gene STARsolo run)")
    args = p.parse_args(argv)

    sdir = args.samples_dir
    bar_path = os.path.join(sdir, "barcodes.tsv")
    feat_path = os.path.join(sdir, "features.tsv")
    mtx_path = os.path.join(sdir, f"umiDedup-{args.dedup}.mtx")
    for req in (bar_path, feat_path, mtx_path):
        if not os.path.isfile(req) or os.path.getsize(req) == 0:
            print(f"ERROR: missing or empty required input: {req}", file=sys.stderr)
            return 1

    try:
        samples = read_barcodes(bar_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    try:
        n_rows, features = read_features(feat_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    try:
        rows, cols, triples = read_mtx(mtx_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    if rows != n_rows:
        print(f"ERROR: matrix has {rows} rows but features.tsv has {n_rows} lines "
              f"({mtx_path} vs {feat_path})", file=sys.stderr)
        return 2
    if cols != len(samples):
        print(f"ERROR: matrix has {cols} columns but barcodes.tsv has {len(samples)} "
              f"sample names ({mtx_path} vs {bar_path})", file=sys.stderr)
        return 2

    # Per-row sparse counts: counts[row_idx][sample_idx] = count.
    row_counts: list[dict[int, float]] = [dict() for _ in range(rows)]
    for r, c, v in triples:
        row_counts[r][c] = row_counts[r].get(c, 0.0) + float(v)

    # Fold duplicate gene symbols into one row.  The matrix row list may
    # contain a symbol more than once (mouse GTFs: e.g. two ENSMUSG entries
    # both named "Snhg4"); STARsolo emits a matrix row per GTF entry, and
    # HTSeq's default 'union' mode is the canonical way to collapse them.
    # We mirror that: one output row per (symbol, feature_type), summing the
    # counts over the rows that share the symbol, keeping the id of the first
    # entry. Empty symbols fall back to the gene id as their label.
    merged: list[tuple[str, str, str, str, dict[int, float]]] = []
    index: dict[tuple[str, str], tuple[str, str, str, str, dict[int, float]]] = {}
    n_dup_rows = 0
    for idx, (gid, sym, typ) in enumerate(features):
        label = sym or gid
        key = (label, typ)
        ent = index.get(key)
        if ent is None:
            ent = (label, gid, sym, typ, dict(row_counts[idx]))
            index[key] = ent
            merged.append(ent)
        else:
            n_dup_rows += 1
            for c, v in row_counts[idx].items():
                ent[4][c] = ent[4].get(c, 0.0) + v

    def fmt(v: float) -> str:
        f = float(v)
        return str(int(f)) if f.is_integer() else str(f)

    # First output column is always the gene symbol (empty symbols already
    # fell back to the id when the rows were folded).
    out_path = args.out_csv
    parent = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(parent, exist_ok=True)

    written = 0
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(",".join(["gene_symbol"] + samples) + "\n")
        for label, gid, sym, typ, vals in merged:
            if not args.all_features and typ and typ != "Gene Expression":
                continue
            if not label:
                continue  # unlabelled feature row: nothing safe to call it
            out.write(label + "," + ",".join(
                fmt(vals[c]) if c in vals else "0" for c in range(cols)
            ) + "\n")
            written += 1

    dropped = len(merged) - written
    if written == 0:
        print("ERROR: no feature rows passed the filter — check --all-features "
              "and the feature_type column in features.tsv", file=sys.stderr)
        return 2

    dup_note = f", {n_dup_rows} folded into duplicate-label rows" if n_dup_rows else ""
    print(
        f"OK  {out_path}\n"
        f"    rows (genes)  = {written} of {rows} matrix-rows{dup_note}"
        + (f", {dropped} dropped by feature filter" if dropped else "")
        + f"\n    cols (samples)= {cols}\n"
        f"    matrix source = {os.path.basename(mtx_path)} (nnz={len(triples)})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
