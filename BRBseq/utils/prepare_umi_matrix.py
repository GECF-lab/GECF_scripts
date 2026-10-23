#!/usr/bin/env python3
"""prepare_umi_matrix.py — emit a gene x sample matrix CSV from the BRB-seq MTX.

The STARsolo + rename_matrix output is an MTX triplet (rows = genes,
columns = samples), which most downstream tools (PCA in R / Excel / CLC)
want as a single CSV.  The legacy BRB pipeline produced this as
``<pool>_counts.detailed_forPCA.csv``; this script is the new pipeline's
equivalent.  brbseq_pipeline.sh calls it twice: with `--dedup
1MM_Directional` it writes the UMI matrix (``<pool>_umi_matrix.csv``,
the canonical table for PCA / DE), and with `--dedup NoDedup` the reads
matrix (``<pool>_reads_matrix.csv``, all PCR copies — kept for
duplication context, matches the old delivered read count).

Columns are the Sample_IDs from ``barcodes.tsv`` (already remapped by
rename_matrix.py — well labels, project codes, ... anything the manifest
used).  Row order and contents follow ``features.tsv`` exactly; the row key
is the unique gene id (column ``Gene_id``), with the symbol as an annotation
column (``gene_symbol``) — duplicates are impossible because every row has
its own id.

Usage:
    prepare_umi_matrix.py --samples-dir <OUT>/star/samples \
                          --out-csv   <pool>_umi_matrix.csv \
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
            raise ValueError("this is not recognisable as a well-formed BRB-seq/STAR 'MatrixMarket' count table — there is unexpected content before the part of the file that identifies it as such. This file was supposed to be produced by the STAR step; if it looked intact when you checked it, the problem is elsewhere in the pipeline and worth re-running that step")
        i += 1
    else:
        raise ValueError("this is not recognisable as a well-formed BRB-seq/STAR 'MatrixMarket' count table — the identifying header line is missing entirely. This file was supposed to be produced by the STAR step; the most likely explanation is that STAR's write-out of this file was interrupted, so re-running the earlier alignment step is the right move")

    if not s.startswith("%%matrixmarket matrix coordinate"):
        raise ValueError("this file identifies itself as a MatrixMarket file, but not in the 'coordinate' sub-format this step expects to read — STAR writes out the coordinate sub-format for BRB-seq, so a different sub-format here points to file corruption or a copy/move that went wrong")

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
            print(f"ERROR: The step that turns the BRB-seq count table into a plain gene-by-sample CSV is missing one of its required inputs: {req} is either not present or completely empty. This points to a problem in one of the earlier steps (the STAR counting step or the sample-renaming step) rather than in this step — see the messages printed above for where the earlier problem started. Fix that, then re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
            return 1

    try:
        samples = read_barcodes(bar_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: Could not build the finished CSV of gene counts: the file listing which sample name belongs to which column ({os.path.basename(bar_path)}) could not be read or contains a problem: {e}. This is generated by the earlier sample-renaming step, so the root cause is most likely there — check that step's messages (above this one) and its result folder, fix the problem, and re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
        return 2

    try:
        n_rows, features = read_features(feat_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: Could not build the finished CSV of gene counts: the list of gene names/IDs this step would use as row labels ({os.path.basename(feat_path)}) could not be read or contains a problem: {e}. This is generated by the earlier STAR / sample-renaming steps, so the root cause is most likely there — check those steps' messages (above this one) and their result folders, fix the problem, and re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
        return 2

    try:
        rows, cols, triples = read_mtx(mtx_path)
    except (OSError, ValueError) as e:
        print(f"ERROR: Could not build the finished CSV of gene counts: the count table itself ({os.path.basename(mtx_path)}) could not be read or its internal layout is inconsistent: {e}. This is generated by the earlier STAR / sample-renaming steps, so the root cause is most likely there — check those steps' messages (above this one) and their result folders, fix the problem, and re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
        return 2

    if rows != n_rows:
        print(f"ERROR: Refusing to build the finished CSV of gene counts: the count table has {rows} gene rows, but the list of gene names/IDs says there should be {n_rows}. These two should always agree, so a disagreement means one of the two files was produced by a different (partial) run than the other. This is not a problem with your sample names — it points to the earlier STAR / sample-renaming steps. Check those steps' messages (above this one) and their result folders, fix the problem, and re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
        return 2
    if cols != len(samples):
        print(f"ERROR: Refusing to build the finished CSV of gene counts: the count table has {cols} sample columns, but the list of sample names says there should be {len(samples)}. These two should always agree, so a disagreement means one of the two files was produced by a different (partial) run than the other. This is not a problem with your sample names — it points to the earlier STAR / sample-renaming steps. Check those steps' messages (above this one) and their result folders, fix the problem, and re-run the pipeline (using --resume-dir pointing at this same folder, so steps that already finished are not redone).", file=sys.stderr)
        return 2

    # Per-row sparse counts: counts[row_idx][sample_idx] = count.
    row_counts: list[dict[int, float]] = [dict() for _ in range(rows)]
    for r, c, v in triples:
        row_counts[r][c] = row_counts[r].get(c, 0.0) + float(v)

    # ONE output row per matrix row, in features.tsv order.  Cell values are
    # written exactly as they appear in the matrix — no summing, no folding.
    # The row key is the gene ID (unique: one ENSMUSG id per matrix row),
    # followed by the gene symbol as a human-readable annotation column —
    # the same two-key layout as the old FastReadCounter deliveries
    # (Gene_id, Gene_name).  Duplicate symbols (mouse GTFs carry ~316: e.g.
    # two entries both named "Snhg4") are therefore disambiguated by id and
    # need NO folding; downstream tools keep the id as the unique key.

    def fmt(v: float) -> str:
        f = float(v)
        return str(int(f)) if f.is_integer() else str(f)

    # Two leading key columns: Gene_id (unique) + gene symbol (annotation).
    # Rows carry a gene id for every well-behaved feature list, so the skip
    # guard below only drops a row when BOTH id and symbol are empty.
    out_path = args.out_csv
    parent = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(parent, exist_ok=True)

    written = 0
    dup_symbols = 0
    sym_seen: dict[str, int] = dict()
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(",".join(["Gene_id", "gene_symbol"] + samples) + "\n")
        for idx, (gid, sym, typ) in enumerate(features):
            if not args.all_features and typ and typ != "Gene Expression":
                continue
            if not gid and not sym:
                continue  # unlabelled feature row: nothing safe to call it
            vals = row_counts[idx]
            out.write(gid + "," + sym + "," + ",".join(
                fmt(vals[c]) if c in vals else "0" for c in range(cols)
            ) + "\n")
            if sym:
                sym_seen[sym] = sym_seen.get(sym, 0) + 1
            written += 1

    dropped = len(features) - written
    if written == 0:
        print("ERROR: Could not build the finished CSV of gene counts: after applying the filter for which gene rows to include, no rows were left at all — the gene list this step reads (features.tsv, produced by the earlier STAR step) either has no entries matching the 'Gene Expression' type this step expects, or every single entry is missing both an ID and a name. This is not a problem with your sample names. If you really do need non-standard gene rows included in the final table too, re-run this step's call with the --all-features switch (this normally means re-running the pipeline with different settings, not re-running this step by hand).", file=sys.stderr)
        return 2

    dup_symbols = sum(1 for c in sym_seen.values() if c > 1)
    dup_note = (f", {dup_symbols} symbols shared by >1 gene id "
                f"(disambiguated by Gene_id, not folded)" if dup_symbols else "")
    print(
        f"OK  {out_path}\n"
        f"    rows (genes)  = {written} of {rows} matrix-rows{dup_note}"
        + (f", {dropped} dropped by feature filter" if dropped else "")
        + f"\n    cols (samples)= {cols}\n"
        f"    key             = Gene_id (+ gene_symbol annotation)\n"
        f"    matrix source = {os.path.basename(mtx_path)} (nnz={len(triples)})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
