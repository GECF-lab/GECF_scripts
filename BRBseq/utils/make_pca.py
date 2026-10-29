#!/usr/bin/env python
# =============================================================================
#  make_pca.py  (v1.0.0)
# ----------------------------------------------------------------------------
#  PCA visualization of a gene x well UMI count matrix.
#
#  INPUT:   gene x well UMI count CSV (genes as rows, one column per well):
#           Gene_id,gene_symbol,<well1>,<well2>,...
#           also accepts: Gene_name,<sample1>,<sample2>,...
#
#  OUTPUT:  (with --outdir DIR, else <input_dir>/pca_<matrixname>/):
#           pca_<pool>.png  pca_<pool>_depth.png  pca_<pool>_elbow.png
#           well_qc.csv     (well, total_umi, genes_detected - one row per well)
#
#  BEHAVIOR: all wells from the input matrix are kept and plotted; every well
#            is annotated by its full name on both PC1/PC2 plots (no legend).
#
#  DEPS:    pip install scanpy anndata matplotlib adjustText   (Python >= 3.9, no R needed)
#
#  USAGE:   python make_pca.py <umi_matrix.csv> --outdir DIR
#          [--n-hvg 1500] [--n-pc 15] [--min-total-counts 5]
# =============================================================================
import argparse, os, sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

VERSION = "1.0.0"

# Gene-filter rule: keep a gene only if its TOTAL count across ALL samples
# (wells) is greater than this threshold.
MIN_TOTAL_COUNTS = 5

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_csv", help="genes x wells UMI count CSV")
    ap.add_argument("--version", action="version", version="make_pca.py " + VERSION)
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: <input_dir>/pca_<matrixname>/)")
    ap.add_argument("--min-total-counts", type=int, default=MIN_TOTAL_COUNTS,
                    help="keep genes whose TOTAL count across all wells > this (default 5)")
    ap.add_argument("--n-hvg",        type=int,   default=1500)
    ap.add_argument("--n-pc",         type=int,   default=15)
    a = ap.parse_args()

    if not os.path.exists(a.input_csv):
        sys.exit("Input not found: " + a.input_csv)
    outdir = a.outdir or os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(a.input_csv)),
        "pca_viz_" + os.path.basename(a.input_csv).replace(".csv", "")))
    os.makedirs(outdir, exist_ok=True)
    # pool name for output filenames: the input matrix file's basename, with
    # the extension and a trailing *_umi_matrix / *_reads_matrix / _matrix
    # removed (e.g. Lib_01_umi_matrix.csv -> Lib_01)
    pool = os.path.basename(a.input_csv)
    pool = pool[:-4] if pool.endswith(".csv") else pool
    for _suffix in ("_umi_matrix", "_reads_matrix", "_matrix"):
        if pool.endswith(_suffix):
            pool = pool[:-len(_suffix)]
            break
    print("INPUT ", a.input_csv); print("OUTPUT", outdir)
    print("make_pca.py v" + VERSION)

    # ---- 1. load + transpose (wells = rows, genes = cols) ------------------
    # Two accepted layouts:
    #   (a) new GECF BRB-seq pipeline: Gene_id, gene_symbol, <well1>, <well2>, ...
    #   (b) GECF BRB-seq pipeline / CLC:  Gene_name, <sample1>, <sample2>, ...
    raw = pd.read_csv(a.input_csv, dtype=str)
    if {"Gene_id", "gene_symbol"}.issubset(raw.columns):
        key, symb = "Gene_id", "gene_symbol"
    elif "Gene_name" in raw.columns:
        key, symb = "Gene_name", "Gene_name"   # gene name used as both id and label
    else:
        sys.exit("Unknown layout. First cols: " + ",".join(list(raw.columns)[:4]))
    wells = [c for c in raw.columns if c not in (key, symb)]
    gid   = list(raw[key])
    symb  = list(raw[symb])
    if len(set(gid)) != len(gid):
        # duplicate gene labels — fold by summing (HTSeq-union-ish). Common in
        # older GECF BRB-seq outputs (gene-name-only files).
        dup = raw[raw[key].duplicated(keep=False)][key].value_counts()
        print(f"WARNING: {len(dup)} gene labels appear >1x. Folding (summing) "
              f"their columns: top-5 {list(dup.head(5).items())}")
        raw[wells] = raw[wells].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        raw = raw.groupby(key, as_index=False, dropna=False)[wells].sum(numeric_only=True)
        gid  = list(raw[key])
        symb = list(raw[key])
    else:
        raw[wells] = raw[wells].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X = raw[wells].to_numpy(dtype=np.float32)
    ad_ = sc.AnnData(X=X.T,
                     obs=pd.DataFrame(index=wells),
                     var=pd.DataFrame({"gene_symbol": symb}, index=gid, columns=["gene_symbol"]))
    ad_.var.index.name = "gene_id"

    # ---- 2. well QC (per-well totals for the well_qc.csv log) --------------
    ad_.obs.index.name = "well"
    ad_.obs["total_umi"]      = ad_.X.sum(axis=1).astype(int)
    ad_.obs["genes_detected"] = (ad_.X > 0).sum(axis=1).astype(int)
    if ad_.n_obs < 2:
        sys.exit("Fewer than 2 wells in the matrix - nothing to plot.")
    ad_.obs[["total_umi", "genes_detected"]].to_csv(
        os.path.join(outdir, "well_qc.csv"), index=True)

    # ---- 3. gene filter -----------------------------------------------------
    # Rule: keep a gene only if its TOTAL count across ALL wells > --min-total-counts
    # (default 5), regardless of how many wells it appears in.
    totals = ad_.X.sum(axis=0)   # total count of each gene across ALL wells
    keep_g = totals > a.min_total_counts
    print(f"genes: {ad_.n_vars} -> {int(keep_g.sum())} "
          f"(dropped {int((~keep_g).sum())} with total count <= {a.min_total_counts})")
    ad_ = ad_[:, keep_g].copy()

    # ---- 4. normalise + HVG -------------------------------------------------
    # Normalise first, then HVG — this is the standard order: HVGs are chosen on
    # the log-scale which is where the biological signal lives, and it avoids
    # integer-count overflow (zero-variance genes -> inf dispersion ->
    # pandas.cut failure) that happens on raw counts for sparse matrices.
    sc.pp.normalize_total(ad_, target_sum=1e4)
    sc.pp.log1p(ad_)
    # Drop any gene still with non-finite dispersion (shouldn't happen after
    # normalise, but be defensive).
    _finite = np.isfinite(ad_.X).all(axis=0)
    ad_ = ad_[:, _finite].copy()
    sc.pp.highly_variable_genes(ad_, n_top_genes=a.n_hvg, flavor="seurat",
                                min_mean=0.01, subset=False)
    hvg = ad_.var["highly_variable"]
    print("HVGs:", int(hvg.sum()), "/", ad_.n_vars)

    # ---- 5. centre-scale on HVGs + PCA (data is already log-norm'd) --------
    # sc.pp.scale: zero-centers each gene, scales each to unit variance, then
    # clips to ±10 so a single outlier well cannot dominate the SVD.
    # (Replaces the old hand-rolled (Xh - mu)/sd.)
    adata_hvg = ad_[:, ad_.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, zero_center=True, max_value=10)
    M = np.asarray(adata_hvg.X, dtype=np.float64)
    # NB: ad_.X stays log-norm (full gene set) so rank_genes_groups / t-test see real means.
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    npca = max(1, min(a.n_pc, U.shape[1], M.shape[1] - 1))
    ad_.obsm["X_pca"] = U[:, :npca] * S[:npca]
    ad_.uns["pca"] = {
        "variance": (S**2) / (M.shape[0] - 1),
        "variance_ratio_": (S**2) / (S**2).sum(),
    }

    # ---- 6. plots (matplotlib only; no scanpy-plot dep needed) -------------
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    Xp = ad_.obsm["X_pca"]
    # variance explained by each PC, for the axis labels
    _s2 = S**2
    def _pct(i):  # i = 0 -> PC1, 1 -> PC2
        return 100.0 * (_s2[i] / _s2.sum()) if i < len(_s2) else 0.0
    lbl1 = f"PC 1 ({_pct(0):.1f}%)"
    lbl2 = f"PC 2 ({_pct(1):.1f}%)"

    # (a) PC1 vs PC2, one colour per well
    cmap = plt.get_cmap("tab20")

    def _annotate(ax, xs, ys):
        """Write the full sample name next to each dot. Uses adjustText for
        de-overlapping when it is installed; otherwise a simple fixed offset."""
        names = list(ad_.obs.index)
        try:
            from adjustText import adjust_text
            texts = [ax.text(x, y, s, fontsize=7, alpha=0.9)
                     for x, y, s in zip(xs, ys, names)]
            adjust_text(texts, x_coords=xs, y_coords=ys,
                        expand_points=(1.25, 1.25), only_move={"texts": "xy"},
                        arrowprops=dict(arrowstyle="-", color="0.55", lw=0.5,
                                        shrinkA=2, shrinkB=2))
        except Exception:
            for x, y, s in zip(xs, ys, names):
                ax.annotate(s, (x, y), textcoords="offset points", xytext=(5, 4),
                            fontsize=7, alpha=0.85)

    def _save(name):
        fig = plt.gcf()
        fig.savefig(os.path.join(outdir, name + ".png"), dpi=150)
        plt.close(fig)

    plt.figure(figsize=(8, 7))
    for i, _w in enumerate(ad_.obs.index.tolist()):
        plt.scatter(Xp[i, 0], Xp[i, 1], c=[cmap(i / max(1, ad_.n_obs))], s=12, linewidth=0)
    plt.xlabel(lbl1); plt.ylabel(lbl2); plt.title("PC1 vs PC2")
    _annotate(plt.gca(), Xp[:, 0], Xp[:, 1])
    plt.tight_layout()
    _save("pca_" + pool)

    # (c) by depth (total UMI per well)
    plt.figure(figsize=(8, 7))
    sc_ = plt.scatter(Xp[:, 0], Xp[:, 1], c=ad_.obs["total_umi"].values, cmap="inferno", s=12, linewidth=0)
    plt.colorbar(sc_, label="total UMI per well")
    plt.xlabel(lbl1); plt.ylabel(lbl2); plt.title("PC1 vs PC2 [well depth]")
    _annotate(plt.gca(), Xp[:, 0], Xp[:, 1])
    plt.tight_layout(); _save("pca_" + pool + "_depth")

    # (d) elbow
    vr = np.array([(s**2)/(S**2).sum() for s in S[:npca]])
    plt.figure(figsize=(9, 4))
    plt.bar(range(1, len(vr)+1), vr*100)
    plt.xlabel("PC"); plt.ylabel("variance explained (%)")
    plt.title("PCA variance"); plt.tight_layout()
    _save("pca_" + pool + "_elbow")
    print("DONE ->", outdir)

if __name__ == "__main__":
    main()
