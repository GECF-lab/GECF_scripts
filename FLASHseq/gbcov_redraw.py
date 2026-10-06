#!/usr/bin/env python
"""gbcov_redraw.py -- regenerate RSeQC geneBody_coverage plots from an EXISTING
table (no BAMs, no numpy, no RSeQC install needed).

Use when a run (per-sample or a big combined `-i bam1,bam2,...` call) finished
writing the coverage table but then crashed before producing the .r/PDFs --
e.g. the RSeQC 5.0.5 ZeroDivisionError on a flat (all-zero) profile -- or when
the plots were skipped (--skip-plot / no Rscript at the time).

  INPUT  : a *.geneBodyCoverage.txt table (header row "Percentile 1 2 ... 100",
           one row per sample, 100 values each)
  OUTPUT : <prefix>.heatMap.<fmt>   (>= 3 samples; heatmap ranked by skewness)
           <prefix>.curves.<fmt>    (line plot of all profiles)
           <prefix>.r               (generated R script; pdf/png/jpeg only,
                                    not written for html)

  NOTE   : --format html renders BOTH plots directly in Python (plotly) as
           self-contained interactive .html files -- no R needed. This is the
           readable 300+-sample case: hover any curve for its name, click a
           legend entry to hide/show, and zoom/pan.

Math + R output are identical to RSeQC 5.0.5, with two guarded exceptions
only for CONSTANT rows (RSeQC 5.0.5 dies on these):
  * skewness of a zero-variance row   -> 0.0  (was: float division by zero)
  * min-max normalization of a flat row -> 0.0 (was: float division by zero)
  All other rows are computed exactly as before.

Usage:
  python3 gbcov_redraw.py TABLE.geneBodyCoverage.txt [options]

Options:
  --format {pdf,png,jpeg,html}
                           plot format                (default: pdf; html =
                           interactive browser output via plotly, no R)
  --prefix PATH              output prefix              (default: derived from
                                               the table path)
  --rscript PATH             Rscript executable         (default: PATH, then
                                               common conda locations)
  --dry-run                  parse table + write .r only; do not run R
  --limit N                  plot first N samples only  (quick sanity check)

Exit codes: 0 ok | 2 bad/empty table | 3 --format html: plotly missing
(pip install plotly) | 5 Rscript missing for pdf/png/jpeg (.r still written;
run it manually) | R's own exit code on R failure.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
from operator import itemgetter

RSCRIPT_CANDIDATES = (
    None,  # "Rscript" on PATH
    "/home/gecf/softwares/anaconda3/envs/flashseq/bin/Rscript",
    "/home/gecf/softwares/anaconda/envs/flashseq/bin/Rscript",
    os.path.expanduser("~/.conda/envs/flashseq/bin/Rscript"),
)


def skewness(values):
    """Pearson-moment skewness (RSeQC algorithm); 0.0 for a constant row."""
    middle = values[len(values) // 2]
    n = len(values)
    mean = sum(values) / n
    if n > 1:
        sigma = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1))
    else:
        sigma = 0.0
    if sigma == 0.0:
        return 0.0
    cubes = [((v - middle) / sigma) ** 3 for v in values]
    return sum(cubes) / n


def load_dataset(path):
    """Parse the table; return rows sorted by skewness (most 3'-biased first)."""
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("Percentile"):
                continue
            fields = stripped.split()
            name = fields[0]
            values = [float(v) for v in fields[1:]]
            if not values:
                continue
            minimum, maximum = min(values), max(values)
            if maximum == minimum:  # flat profile: no range -> zero curve
                norm = [0.0] * len(values)
            else:                  # RSeQC min-max normalization
                norm = [(v - minimum) / (maximum - minimum) for v in values]
            rows.append((name, norm, skewness(values)))
    rows.sort(key=itemgetter(2), reverse=True)
    return rows


def write_r_code(dataset, file_prefix, output_format, column_count=100):
    """Generate the R script, mirroring RSeQC 5.0.5 write_r_code()."""
    r_script = f"{file_prefix}.r"
    names = [name for name, _d, _s in dataset]
    fmt = output_format.lower()
    ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    base_plot = (
        f"plot(x,{names[0]},type='l',"
        "xlab=\"Gene body percentile (5'->3')\","
        "ylab=\"Coverage\",lwd=0.8,col=icolor[1])"
    )
    legend = (
        f"legend(0,1,fill=icolor[1:{len(names)}],legend=c("
        + ",".join(f"'{n}'" for n in names)
        + "))"
    )
    lines_block = [
        f"lines(x,{name},type='l',col=icolor[{i}])"
        for i, name in enumerate(names[1:], start=2)
    ]

    with open(r_script, "w", encoding="utf-8") as output:
        for name, data, _skewness in dataset:
            output.write(
                f"{name} <- c({','.join(str(v) for v in data)})\n"
            )

        if len(names) >= 3:
            output.write(
                f"data_matrix <- matrix(c({','.join(names)}), byrow=T, "
                f"ncol={column_count})\n"
            )
            output.write(
                "rowLabel <- c("
                + ",".join(f'"{n}"' for n in names)
                + ")\n\n"
            )
            output.write(f'{fmt}("{file_prefix}.heatMap.{fmt}")\n')
            output.write("rc <- cm.colors(ncol(data_matrix))\n")
            output.write(
                "heatmap(data_matrix, scale=c(\"none\"),keep.dendro=F, "
                "labRow=rowLabel,Colv=NA,Rowv=NA,labCol=NA,"
                "col=cm.colors(256),margins=c(6,8),ColSideColors=rc,"
                "cexRow=1,cexCol=1,"
                "xlab=\"Gene body percentile (5'->3')\","
                "add.expr=x_axis_expr <- axis("
                f"side=1,at=c({','.join(str(i) for i in ticks)}),labels=c("
                + ",".join(f'"{i}"' for i in ticks)
                + ")))\n"
            )
            output.write("dev.off()\n")
        output.write("\n")

        output.write(f'{fmt}("{file_prefix}.curves.{fmt}")\n')
        output.write(f"x=1:{column_count}\n")
        output.write(
            "icolor = colorRampPalette(c("
            '"#7fc97f","#beaed4","#fdc086","#ffff99",'
            '"#386cb0","#f0027f"))'
            f"({len(names)})\n"
        )

        if len(names) == 1:
            output.write(base_plot + "\n")
        elif len(names) <= 6:
            output.write(base_plot + "\n")
            output.writelines(l + "\n" for l in lines_block)
            output.write(legend + "\n")
        else:
            output.write(
                "layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2),"
                "4,4,byrow=TRUE))\n"
            )
            output.write(base_plot + "\n")
            output.writelines(l + "\n" for l in lines_block)
            output.write("par(mar=c(1,0,2,1))\n")
            output.write("plot.new()\n")
            output.write(legend + "\n")

        output.write("dev.off()\n")
    return r_script


def _ramp(n, ramp_colors):
    """n evenly-spaced colors linearly interpolated over the ramp endpoint
    colors; mimics R's colorRampPalette(colors)(n). Pure RGB math, returns a
    list of '#rrggbb' hex strings (no plotly color-API dependency)."""
    def _rgb(hexcolor):
        h = hexcolor.lstrip("#")
        return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))

    stops = [_rgb(c) for c in ramp_colors]
    seg = len(stops) - 1
    out = []
    for i in range(n):
        if n <= 1:
            t_idx, local = 0, 0.0
        else:
            pos = (i / (n - 1)) * seg
            t_idx = min(int(pos), seg - 1)
            local = pos - t_idx
        a, b = stops[t_idx], stops[t_idx + 1]
        rgb = tuple(round(a[k] + local * (b[k] - a[k])) for k in range(3))
        out.append("#%02x%02x%02x" % rgb)
    return out


def write_html(dataset, file_prefix):
    """Render the two plots as self-contained interactive HTML via plotly,
    and return the list of written file paths (heatMap first when present).

    No Rscript needed. The legend problem of the static RSeQC output is
    solved: each curve is a hoverable / toggable trace in an interactive
    legend and heatmap rows carry their sample name on hover.
    """
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise SystemExit(
            "error: --format html requires the 'plotly' dependency "
            "(not installed):  pip install plotly"
        ) from exc

    names = [name for name, _d, _s in dataset]
    values = [vals for _n, vals, _s in dataset]
    n = len(names)
    column_count = len(values[0]) if n else 100
    x = list(range(1, column_count + 1))
    ramp = _ramp(n, ['#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f'])

    paths = []

    # --- heatMap (only when RSeQC itself would emit one: >= 3 samples) ---
    if n >= 3:
        heat_path = f"{file_prefix}.heatMap.html"
        fig = go.Figure(
            data=go.Heatmap(
                z=values,
                x=x,
                y=names,
                colorscale="Jet",
                zmin=0,
                zmax=1,
                colorbar=dict(title="Coverage (min-max normalized)"),
                hovertemplate="(sample, percentile)\n"
                             "Sample: %{y}<br>Pctile: %{x}<br>"
                             "Coverage: %{z:.3f}"
                             "<extra></extra>",
            )
        )
        fig.update_layout(
            title="gbcov heatMap (ranked by skewness, 3'-biased on top)",
            xaxis=dict(title="Gene body percentile (5'->3')"),
            yaxis_autorange="reversed",  # first sample (most skewed) on top, like R
            width=1000,
            height=min(400 + 6 * n, 2600),
            margin=dict(t=70, r=30, b=80, l=240),
        )
        fig.write_html(heat_path, include_plotlyjs="inline")
        paths.append(heat_path)

    # --- curves ---
    curves_path = f"{file_prefix}.curves.html"
    fig = go.Figure()
    for i, (name, vals) in enumerate(zip(names, values)):
        fig.add_scatter(
            x=x,
            y=vals,
            mode="lines",
            name=name,
            line=dict(width=1.2, color=ramp[i]),
            hovertemplate="Sample: %{fullData.name}<br>"
                          "Pctile: %{x}<br>Cov: %{y:.3f}"
                          "<extra></extra>",
        )
    fig.update_layout(
        title="gbcov curves (hover to identify, click legend to hide/show, "
              "scroll to zoom)",
        xaxis=dict(
            title="Gene body percentile (5'->3')",
            tickmode="array",
            tickvals=[1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
        ),
        yaxis=dict(title="Coverage (min-max normalized)", range=[0, 1]),
        legend=dict(title=dict(text="Sample"), bgcolor="rgba(255,255,255,0.85)"),
        height=720,
        width=1100,
        hovermode="closest",
    )
    fig.write_html(curves_path, include_plotlyjs="inline")
    paths.append(curves_path)
    return paths


def find_rscript(preferred=None):
    candidates = (preferred,) if preferred else RSCRIPT_CANDIDATES
    for cand in candidates:
        if cand is None:
            found = shutil.which("Rscript")
            if found:
                return found
        elif os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Regenerate RSeQC geneBody plots from an existing "
        "coverage table (no BAM re-read).",
    )
    parser.add_argument(
        "table",
        nargs="?",
        help="*.geneBodyCoverage.txt table file (positional; or use --table)",
    )
    parser.add_argument(
        "--table",
        dest="table_opt",
        default=None,
        help="alias for the positional table argument "
        "(flashseq_rseqc.sh >= v1.1.0 calls gbcov_redraw.py this way)",
    )
    parser.add_argument(
        "--format",
        choices=("pdf", "png", "jpeg", "html"),
        default="pdf",
    )
    parser.add_argument(
        "--prefix", help="output prefix (default: table path minus .txt)"
    )
    parser.add_argument("--rscript", help="Rscript executable path")
    parser.add_argument("--dry-run", action="store_true", help="write the .r only")
    parser.add_argument("--limit", type=int, help="plot first N samples only")
    args = parser.parse_args()

    # accept both "TABLE" (positional) and "--table TABLE"; prefer --table when given
    table = args.table_opt or args.table
    if not table:
        parser.error("table not specified (positional argument or --table)")
    args.table = table

    if not os.path.isfile(args.table):
        parser.error(f"table not found: {args.table}")

    prefix = args.prefix
    if prefix is None:
        prefix = args.table
        for suffix in (".geneBodyCoverage.txt", ".txt"):
            if prefix.endswith(suffix):
                prefix = prefix[: -len(suffix)]
                break

    if not os.path.isdir(os.path.dirname(os.path.abspath(prefix))):
        parser.error(f"output directory does not exist for prefix: {prefix}")

    dataset = load_dataset(args.table)
    if args.limit:
        dataset = dataset[: args.limit]
    if not dataset:
        print(f"error: no valid coverage rows in {args.table}", file=sys.stderr)
        return 2

    flat = [n for n, vals, _s in dataset if max(vals) == 0.0]
    print(f"samples in table  : {len(dataset)}", file=sys.stderr)
    if flat:
        print(
            f"flat (all-zero) profiles ({len(flat)}) -- treated as skewness 0 / "
            f"zero curve (RSeQC 5.0.5 would crash on these):",
            file=sys.stderr,
        )
        for n in flat:
            print(f"  {n}", file=sys.stderr)

    fmt = args.format.lower()
    if fmt == "html":
        try:
            html_files = write_html(dataset, prefix)
        except SystemExit as exc:
            print(str(exc), file=sys.stderr)
            return 3
        for p in html_files:
            print(f"HTML written      : {p}", file=sys.stderr)
        print("                  (open in a browser -- interactive legend)", file=sys.stderr)
        return 0

    r_script = write_r_code(dataset, prefix, fmt)
    print(f"R script written  : {r_script}", file=sys.stderr)

    if not args.dry_run:
        rscript = find_rscript(args.rscript)
        if rscript is None:
            print(
                "Rscript not found. The .r is written -- run it manually:\n"
                f"  Rscript {r_script}",
                file=sys.stderr,
            )
            return 5
        print(f"running R       : {rscript} {r_script}", file=sys.stderr)
        rc = subprocess.run([rscript, r_script]).returncode
        if rc != 0:
            print(f"error: Rscript exit {rc}", file=sys.stderr)
            return rc

    print(
        f"done            : {prefix}.heatMap.{fmt}, {prefix}.curves.{fmt}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
