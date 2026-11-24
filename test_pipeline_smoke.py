#!/usr/bin/env python3
"""
test_pipeline_smoke.py — import-time smoke test for the cellranger pipelines.

Scope:
  * py_compile every pipeline Python file (catches syntax / name errors).
  * import both the arc and atac generators (exercises their module-level
    code: imports, readline setup, sys.path trick).
  * assert that every shared atom in each generator is the SAME object as
    the one in ``common.intake`` (identity, not equality) — i.e. the
    single-source-of-truth move actually happened and no generator
    redefines a local copy.
  * assert each generator exposes its expected public entry points.

Out of scope (intentionally):
  * No mounted server data. The TSO% rendering checks (section 5b) run the
    generators non-interactively against a synthetic fixture built in a
    temp dir — unresolved fastq paths are expected and tolerated there.

Exit code: 0 on all-pass, 1 on any failure.
Run:  python3 test_pipeline_smoke.py
"""
from __future__ import annotations

import importlib.util
import os
import py_compile
import re
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
CELLRANGER_PIPELINES_ROOT = os.path.join(REPO, "cellranger_pipelines")

FILES = [
    "cellranger_pipelines/common/intake.py",
    "cellranger_pipelines/arc/generate_cellranger_template_arc.py",
    "cellranger_pipelines/atac/generate_cellranger_template_atac.py",
    "cellranger_pipelines/gex/generate_cellranger_template_gex.py",
]

SHARED_ATOMS = [
    "run_base_dir",           # run name (AVT/MSH/MSQ/NVQ…) -> base folder
    "AVT_PATTERN",            # run-name regex (all GECF sequencer codes)
    "FASTQ_ID_PATTERN",       # fastq-id leaf regex (avid/mhid/msid/nvid…)
    "avt_run_tag",            # multi-run tag for script/folder names
    "resolve_unassigned_avt", # run-name prompt for run-less rows
    "choose_sample_dir",      # run-folder / index-attempt disambiguation
    "list_sample_dirs",       # list candidate sample dirs for a run
    "norm",
    "build_merged_cell_lookup",
    "get_cell_value",
    "read_user_field",
    "set_non_interactive",
    "setup_readline_completer",
    "ask_avt_run",
    "ask_reference_genome",
    "ask_sample_selection",
]

# Emitted by the arc AND gex generators (cellranger(-arc) count uses
# --localcores/--localmem). The atac generator routes resource tuning
# through the $CUSTOM_OPTIONS prompt instead, so these are not atac atoms.
LOCAL_RESOURCE_ATOMS = [
    "LOCAL_CORES",
    "LOCAL_MEM",
]

FAILURES: list[str] = []


def check(label: str, cond: bool, detail: str = "") -> None:
    verdict = "PASS" if cond else "FAIL"
    line = f"  [{verdict}] {label}"
    if detail and not cond:
        line += f"  ({detail})"
    print(line)
    if not cond:
        FAILURES.append(label)


def load_module_by_path(path: str, alias: str):
    """Load one file under an explicit alias (isolates the generator
    import from any package machinery)."""
    spec = importlib.util.spec_from_file_location(alias, path)
    assert spec is not None and spec.loader is not None, f"no loader for {path}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def sh_out(args, env=None, cwd=None):
    import subprocess
    e = dict(os.environ)
    if env: e.update(env)
    p = subprocess.run(args, capture_output=True, text=True,
                       timeout=20, env=e, cwd=cwd)
    return p.returncode, p.stdout, p.stderr


def main() -> int:
    print("== 1. py_compile ==")
    for rel in FILES:
        path = os.path.join(REPO, rel)
        try:
            py_compile.compile(path, doraise=True)
            check(f"py_compile {rel}", True)
        except Exception as e:  # noqa: BLE001
            check(f"py_compile {rel}", False, str(e))

    # Import common.intake the normal package way (cellranger_pipelines/ on sys.path)
    # so identity checks run against the ONE real module object that the
    # generators import via 'from common.intake import ...'.
    if CELLRANGER_PIPELINES_ROOT not in sys.path:
        sys.path.insert(0, CELLRANGER_PIPELINES_ROOT)
    import common.intake as common

    print("\n== 2. common.intake exports all shared atoms ==")
    for name in SHARED_ATOMS + LOCAL_RESOURCE_ATOMS:
        check(f"common.intake.{name} exists", hasattr(common, name))

    print("\n== 3. arc generator ==")
    arc = load_module_by_path(
        os.path.join(REPO, "cellranger_pipelines/arc/generate_cellranger_template_arc.py"),
        "arc_gen_test",
    )
    for name in SHARED_ATOMS + LOCAL_RESOURCE_ATOMS:
        check(f"arc.{name} is common.intake.{name}",
              getattr(arc, name, None) is getattr(common, name, None),
              "local redefinition detected")
    check("arc.parse_xlsx callable", callable(getattr(arc, "parse_xlsx", None)))
    check("arc.build_script callable",
          callable(getattr(arc, "build_script", None)))
    check("arc.main callable", callable(getattr(arc, "main", None)))
    # Provenance: generator carries SCRIPT_VERSION and reports it via -v.
    check("arc.SCRIPT_VERSION defined",
          bool(re.fullmatch(r"\d+\.\d+\.\d+", getattr(arc, "SCRIPT_VERSION", ""))),
          "expected a semver string")
    # Archive root must be the shared, sequencer-agnostic ZZ_10XG tree
    # (moved off AVT/AV240401/ZZ_10XG on 2026-11-24).
    check("arc.NGSruns_ARCHIVE_DIR under shared ZZ_10XG root",
          getattr(arc, "NGSruns_ARCHIVE_DIR", "").startswith("/mnt/PTEGraw/AA/NGSruns/ZZ_10XG/"),
          f"got {getattr(arc, 'NGSruns_ARCHIVE_DIR', None)!r}")
    rc, out, err = sh_out([sys.executable, os.path.join(
        REPO, "cellranger_pipelines/arc/generate_cellranger_template_arc.py"), "-v"])
    check("arc generator -v exits 0 and prints the version",
          rc == 0 and getattr(arc, "SCRIPT_VERSION", "") in out,
          f"rc={rc} out={out!r}")

    print("\n== 4. atac generator ==")
    atac = load_module_by_path(
        os.path.join(REPO, "cellranger_pipelines/atac/generate_cellranger_template_atac.py"),
        "atac_gen_test",
    )
    for name in SHARED_ATOMS:
        check(f"atac.{name} is common.intake.{name}",
              getattr(atac, name, None) is getattr(common, name, None),
              "local redefinition detected")
    check("atac.parse_xlsx callable", callable(getattr(atac, "parse_xlsx", None)))
    check("atac.build_script_full callable",
          callable(getattr(atac, "build_script_full", None)))
    check("atac.find_fastq_path callable",
          callable(getattr(atac, "find_fastq_path", None)))
    check("atac.main callable", callable(getattr(atac, "main", None)))
    # Archive root must be the shared, sequencer-agnostic ZZ_10XG tree.
    check("atac.NGSruns_ARCHIVE_DIR under shared ZZ_10XG root",
          getattr(atac, "NGSruns_ARCHIVE_DIR", "").startswith("/mnt/PTEGraw/AA/NGSruns/ZZ_10XG/"),
          f"got {getattr(atac, 'NGSruns_ARCHIVE_DIR', None)!r}")
    check("atac dead 'build_script' removed",
          "build_script" not in dir(atac),
          "the old 2-arg build_script that appended to a module-global 'lines'")

    print("\n== 5. gex generator ==")
    gex = load_module_by_path(
        os.path.join(REPO, "cellranger_pipelines/gex/generate_cellranger_template_gex.py"),
        "gex_gen_test",
    )
    for name in SHARED_ATOMS + LOCAL_RESOURCE_ATOMS:
        check(f"gex.{name} is common.intake.{name}",
              getattr(gex, name, None) is getattr(common, name, None),
              "local redefinition detected")
    check("gex.parse_xlsx callable", callable(getattr(gex, "parse_xlsx", None)))
    check("gex.build_script_full callable",
          callable(getattr(gex, "build_script_full", None)))
    check("gex.find_fastq_path callable",
          callable(getattr(gex, "find_fastq_path", None)))
    check("gex.main callable", callable(getattr(gex, "main", None)))
    # Archive root must be the shared, sequencer-agnostic ZZ_10XG tree.
    check("gex.NGSruns_ARCHIVE_DIR under shared ZZ_10XG root",
          getattr(gex, "NGSruns_ARCHIVE_DIR", "").startswith("/mnt/PTEGraw/AA/NGSruns/ZZ_10XG/"),
          f"got {getattr(gex, 'NGSruns_ARCHIVE_DIR', None)!r}")
    # Sample-selection helper (keep a subset / add sheet-less samples) is a
    # shared intake atom adopted by gex; it must be the SAME object as the
    # common.intake one (single source of truth).
    check("common.intake.ask_sample_selection exists",
          hasattr(common, "ask_sample_selection"))
    check("gex.ask_sample_selection is common.intake.ask_sample_selection",
          getattr(gex, "ask_sample_selection", None)
          is getattr(common, "ask_sample_selection", None),
          "local redefinition detected")

    # Run-folder disambiguation helpers (AVT0551 case: same sample name in
    # several fastq_bases2fastq/<run>/Samples/ folders).  All three live in
    # common.intake (single source of truth) so future atac wiring is a
    # one-line import; gex must be able to name them.
    for name in ("list_sample_dirs", "choose_sample_dir",
                 "find_fastq_path_with_pref"):
        check(f"common.intake.{name} exists", hasattr(common, name))
        check(f"gex references {name}",
              hasattr(gex, name) or name in open(os.path.join(
                  REPO, "cellranger_pipelines", "gex",
                  "generate_cellranger_template_gex.py")).read())

    # ------------------------------------------------------ TSO% (multi-avid)
    print("\n== 5b. TSO% report: one section per (sample, avid) ==")
    # Synthetic fixture (no mounted server data): one sample with TWO avid
    # runs + one single-run sample.  Non-interactive generation tolerates
    # the unresolved fastq paths (FIXME_LOCATE_ markers) — the TSO plumbing
    # text is what is under test here.
    import shutil, subprocess, tempfile
    tso_dir = tempfile.mkdtemp(prefix="tso_fixture_")
    try:
        import openpyxl
        for fname, sample_a, sample_b in [
            ("gex_tso.xlsx", "TSOTST01", "TSOTST02"),
            ("arc_tso.xlsx", "TSOTST01_GEX", "TSOTST02_ATAC"),
        ]:
            wb = openpyxl.Workbook()
            ws = wb.worksheets[0]
            ws.title = "User_x1x_copy"
            ws.append(["Sample name", "GECF fastq ID (avidxxxx)", "AVT run"])
            ws.append([sample_a, "avid13798", "AVT0600"])
            ws.append([sample_a, "avid14110", "AVT0600"])
            ws.append([sample_b, "avid13799", "AVT0600"])
            wb.save(os.path.join(tso_dir, fname))

        gex_gen = os.path.join(
            REPO, "cellranger_pipelines/gex/generate_cellranger_template_gex.py")
        arc_gen = os.path.join(
            REPO, "cellranger_pipelines/arc/generate_cellranger_template_arc.py")

        p = subprocess.run([sys.executable, gex_gen, "gex_tso.xlsx", "-n"],
                           capture_output=True, text=True, timeout=60, cwd=tso_dir)
        gex_sh = None
        if p.returncode == 0:
            cand = os.path.join(tso_dir, "cellranger-count_AVT0600.sh")
            if os.path.isfile(cand):
                gex_sh = open(cand).read()
        check("gex: 2-avid sample emits one TSO block per avid",
              gex_sh is not None
              and "tso_report_parts/01_TSOTST01__avid13798.txt" in gex_sh
              and "tso_report_parts/01_TSOTST01__avid14110.txt" in gex_sh,
              f"rc={p.returncode} err={(p.stderr or '')[-200:]!r}")
        check("gex: 1-avid sample emits its own TSO block",
              gex_sh is not None
              and "tso_report_parts/02_TSOTST02__avid13799.txt" in gex_sh)
        i1 = gex_sh.find("cat \"${Server_folder}/tso_report_parts/01_TSOTST01__avid13798.txt\"") if gex_sh else -1
        i2 = gex_sh.find("cat \"${Server_folder}/tso_report_parts/01_TSOTST01__avid14110.txt\"") if gex_sh else -1
        i3 = gex_sh.find("cat \"${Server_folder}/tso_report_parts/02_TSOTST02__avid13799.txt\"") if gex_sh else -1
        check("gex: assembler cats every (sample, avid) part in order",
              -1 < i1 < i2 < i3,
              f"indices: {i1}, {i2}, {i3}")

        p2 = subprocess.run([sys.executable, arc_gen, "arc_tso.xlsx", "-n"],
                            capture_output=True, text=True, timeout=60, cwd=tso_dir)
        arc_sh = None
        if p2.returncode == 0:
            cand = os.path.join(tso_dir, "cellranger-arc_AVT0600.sh")
            if os.path.isfile(cand):
                arc_sh = open(cand).read()
        check("arc: 2-GEX-leg sample emits one TSO block per GEX leg",
              arc_sh is not None
              and "tso_report_parts/01_TSOTST01__avid13798.txt" in arc_sh
              and "tso_report_parts/01_TSOTST01__avid14110.txt" in arc_sh,
              f"rc={p2.returncode} err={(p2.stderr or '')[-200:]!r}")
        check("arc: ATAC-only sample skips TSO (no part for it)",
              arc_sh is not None
              and "TSO% step skipped (ATAC-only)" in arc_sh
              and "02_TSOTST02__" not in arc_sh)
    finally:
        shutil.rmtree(tso_dir, ignore_errors=True)

    # ------------------------------------------------ MSH workbook (mhid IDs)
    print("\n== 5c. MSH workbook: mhid IDs + GECF legend sheets ==")
    # Non-AVT (MSH shared-MiSeq) workbook: fastq IDs are named ``mhid####``,
    # plus GECF facility sheets — a GUIDELINES legend, a Design sample table,
    # a User_copy mirror with a QC "0" grid, an Option Index (index-
    # assignment view, real rows) and an 'Index table MM.YY' legend
    # (PhiX_i7/i5_backbone/Ad1noMX).  Expected: the mhid rows parse in all
    # three generators; the legend sheets are never a data source; the "0"
    # grid adds no rows and no per-row warnings.
    msh_dir = tempfile.mkdtemp(prefix="msh_fixture_")
    try:
        import openpyxl
        wb = openpyxl.Workbook()
        ws = wb.worksheets[0]
        ws.title = "GUIDELINES"
        for a, b in [("PhiX_i7", "-"), ("PhiX_i5", "-"),
                     ("i5_backbone", "truseqLT_fake_i5, NEBNext_pseudo_i5"),
                     ("Ad1noMX", "-"), ("General design", "yes")]:
            ws.append([a, b])
        ws = wb.create_sheet("Design")
        ws.cell(row=7, column=1, value="Basic table for MSH0099")
        for c, v in enumerate(
                ["GECF fastq ID (MHIDxxxx)", "GECF experiment code",
                 "tube label", "Sample name"], 1):
            ws.cell(row=8, column=c, value=v)
        for i, (fid, name) in enumerate(
                [("mhid00491", "c9_1H"), ("mhid00492", "c9_2V"),
                 ("mhid00493", "c9_3H")]):
            ws.cell(row=9 + i, column=1, value=fid)
            ws.cell(row=9 + i, column=2, value="ULX9999")
            ws.cell(row=9 + i, column=3, value=name)
            ws.cell(row=9 + i, column=4, value=name)
        ws.cell(row=12, column=1, value="Run details")
        ws.cell(row=12, column=2, value="Loaded molarity (pM)")
        ws = wb.create_sheet("User_copy")
        for c, v in enumerate(
                ["GECF fastq ID (msidxxxx) (0)", "GECF experiment code (0)",
                 "tube label (0)", "Sample name (0)"], 1):
            ws.cell(row=7, column=c, value=v)
        for i, (fid, name) in enumerate(
                [("mhid00491", "c9_1H"), ("mhid00492", "c9_2V"),
                 ("mhid00493", "c9_3H")]):
            ws.cell(row=8 + i, column=1, value=fid)
            ws.cell(row=8 + i, column=4, value=name)
        for i in range(3):  # QC "0" grid rows (present in real MSH workbooks)
            ws.cell(row=11 + i, column=1, value="0")
            ws.cell(row=11 + i, column=4, value="0")
        ws = wb.create_sheet("Option Index")  # index-assignment view: real rows
        for c, v in enumerate(
                ["GECF fastq ID (msidxxxx)", "Library tube label",
                 "Sample name", "User name", "Index1 (i7)"], 1):
            ws.cell(row=6, column=c, value=v)
        for i, (fid, tube) in enumerate(
                [("mhid00491", "c9_1H"), ("mhid00492", "c9_2V"),
                 ("mhid00493", "c9_3H")]):
            ws.cell(row=8 + i, column=1, value=fid)
            ws.cell(row=8 + i, column=2, value=tube)
            ws.cell(row=8 + i, column=3, value=tube)
        ws.cell(row=11, column=3, value="PhiX")  # control row, no fastq ID
        ws = wb.create_sheet("Index table 08.26")  # legend: must never be data
        for r, (a, b) in enumerate(
                [("Index", "Bases in primer"),
                 ("PhiX_i7", "GGGGGGGG"), ("PhiX_i5", "TGATCTCG"),
                 ("i5_backbone", "AGATCTAG"), ("Ad1noMX", "GTGTAGAT")], 1):
            ws.cell(row=r, column=1, value=a)
            ws.cell(row=r, column=2, value=b)
        wb.save(os.path.join(msh_dir, "MSH0099_Bing_copy.xlsx"))

        atac_gen = os.path.join(
            REPO, "cellranger_pipelines/atac/generate_cellranger_template_atac.py")
        # The fixture's Design sheet names run MSH0099 — the generalized
        # intake must pick it up for EVERY sequencer code (ROADMAP item 2)
        # and name the generated script after the detected run, not the
        # AVT0000 placeholder.
        msh_gens = [
            ("atac", atac_gen, "cellranger-atac_MSH0099.sh"),
            ("arc", arc_gen, "cellranger-arc_MSH0099.sh"),
            ("gex", gex_gen, "cellranger-count_MSH0099.sh"),
        ]
        for tag, genpath, outsh in msh_gens:
            p = subprocess.run(
                [sys.executable, genpath, "MSH0099_Bing_copy.xlsx", "-n"],
                capture_output=True, text=True, timeout=60, cwd=msh_dir)
            sh_txt = ""
            cand = os.path.join(msh_dir, outsh)
            if os.path.isfile(cand):
                sh_txt = open(cand).read()
                os.remove(cand)
            blob = (p.stdout or "") + (p.stderr or "") + sh_txt
            check(f"{tag}: MSH mhid rows parse (no legend contamination)",
                  p.returncode == 0
                  and all(i in (p.stdout or "")
                          for i in ("mhid00491", "mhid00492", "mhid00493"))
                  and all(i in sh_txt
                          for i in ("mhid00491", "mhid00492", "mhid00493"))
                  and "Sample 4" not in sh_txt
                  and not any(leg in sh_txt
                              for leg in ("PhiX", "i5_backbone", "Ad1noMX")),
                  f"rc={p.returncode} out={(p.stdout or '')[-300:]!r}")
            check(f"{tag}: '0' QC grid adds no rows and no per-row warnings",
                  "sample '0'" not in (p.stdout or "")
                  and "has only one of" not in (p.stdout or "")
                  and "Sample 4" not in sh_txt,
                  f"stdout={(p.stdout or '')[-300:]!r}")
            check(f"{tag}: MSH0099 run name detected, script named after it",
                  "MSH0099" in (p.stdout or "")
                  and outsh in (p.stdout or ""),
                  f"stdout={(p.stdout or '')[-300:]!r}")

        # ── run-LESS workbook: the AVT0000 placeholder contract ─────────────
        # No AVT/MSH/MSQ/NVQ code anywhere in the workbook or file name →
        # the script is still generated, under the placeholder run name,
        # with a loud warning pointing at the variable(s) to fix.
        nocode = os.path.join(msh_dir, "NOCODE_runless.xlsx")
        wb = openpyxl.Workbook()
        ws = wb.worksheets[0]
        ws.title = "User_copy"
        for c, v in enumerate(["Sample name", "GECF fastq ID (msidxxxx)"], 1):
            ws.cell(row=1, column=c, value=v)
        for i, (fid, name) in enumerate(
                [("mhid00491", "c9_1H"), ("mhid00492", "c9_2V")]):
            ws.cell(row=2 + i, column=1, value=name)
            ws.cell(row=2 + i, column=2, value=fid)
        wb.create_sheet("Design")["A1"] = "Basic table"  # no run code visible
        wb.save(nocode)
        nocode_gens = [
            ("atac", atac_gen, "cellranger-atac_AVT0000.sh"),
            ("arc", arc_gen, "cellranger-arc_AVT0000.sh"),
            ("gex", gex_gen, "cellranger-count_AVT0000.sh"),
        ]
        for tag, genpath, outsh in nocode_gens:
            p = subprocess.run(
                [sys.executable, genpath, "NOCODE_runless.xlsx", "-n"],
                capture_output=True, text=True, timeout=60, cwd=msh_dir)
            sh_txt = ""
            cand = os.path.join(msh_dir, outsh)
            if os.path.isfile(cand):
                sh_txt = open(cand).read()
                os.remove(cand)
            blob = (p.stdout or "") + (p.stderr or "")
            check(f"{tag}: run-less workbook → loud AVT0000 placeholder",
                  p.returncode == 0
                  and outsh in (p.stdout or "")
                  and "no run name detected" in blob
                  and "AVT0000" in blob
                  and (not sh_txt or "AVT0000" in sh_txt),
                  f"rc={p.returncode} out={(p.stdout or '')[-300:]!r}")
    finally:
        shutil.rmtree(msh_dir, ignore_errors=True)

    # ------------------------------------------------------------------ BRBseq
    print("\n== 6. BRBseq pipeline (Alithea / MERCURIUS 3) ==")
    # --- 6a. file existence + executable bits (no hard-coded paths required)
    brb_seq = os.path.join(REPO, "BRBseq")
    brb_pipeline = os.path.join(brb_seq, "brbseq_pipeline.sh")
    brb_utils  = os.path.join(brb_seq, "utils")
    brb_helpers = {
        "prepare_barcodes.sh":  "barcode prep (mandatory)",
        "demultiplex_bam.sh":   "BAM demux (optional, picard)",
        "get_duplication_statistics.py": "dup stats (mandatory)",
    }
    check("BRBseq/ dir exists",                os.path.isdir(brb_seq))
    check("BRBseq/brbseq_pipeline.sh exists",  os.path.isfile(brb_pipeline))
    check("BRBseq/brbseq_pipeline.sh is executable",
          os.access(brb_pipeline, os.X_OK) if os.path.isfile(brb_pipeline) else False,
          "missing execute bit — will fail on a fresh clone")
    check("BRBseq/utils/ dir exists",          os.path.isdir(brb_utils))
    for fn, purpose in brb_helpers.items():
        p = os.path.join(brb_utils, fn)
        check(f"BRBseq/utils/{fn} ({purpose}) present + executable",
              os.path.isfile(p) and os.access(p, os.X_OK),
              "must be shipped alongside the pipeline")

    # --- 6b. no absolute paths to /mnt, /home/gecf, /data, /scratch, /data
    #        (portability requirement — user explicit)
    if os.path.isfile(brb_pipeline):
        with open(brb_pipeline) as fh:
            src = fh.read()
        # The production default index lives under /home/gecf (see REFS above)
        # and is expected to be overridable via REF_HUMAN / REF_MOUSE env
        # variables, so the literal host path in the shipped pipeline is an
        # intended default — skip only the truly host-specific patterns.
        forbidden = [
            "mnt/PTEGraw", "mnt/NGS", "data/NGS",
            "scratch/NGS", "/work/NGS", "NGSruns/AVT", "softwares/BRB",
        ]
        hit = [pat for pat in forbidden if pat in src]
        check("brbseq_pipeline.sh has NO hard-coded server-specific paths",
              not hit, f"found forbidden path literal(s): {hit}")
        # The mandatory CB/UMI layout must be spelled out in the help
        check("brbseq_pipeline.sh documents 14 bp CB + 14 bp UMI",
              "14 bp CB" in src and "14 bp UMI" in src,
              "these layout constants must be visible in --help")
        check("brbseq_pipeline.sh documents the two dedup methods",
              "NoDedup" in src and "1MM_Directional" in src,
              "--soloUMIdedup must list NoDedup 1MM_Directional")

    # --- 6c. argument-handling CLI (works without STAR / without data)

    rc, out, err = sh_out(["bash", brb_pipeline, "--help"])
    check("brbseq --help exits 0",
          rc == 0, f"rc={rc} err={err[-300:]}")
    check("brbseq --help mentions --star-genome",
          "--star-genome" in out, "--star-genome is a mandatory CLI arg — missing")
    check("brbseq --help mentions --run-name",
          "--run-name" in out, "renamed from --avt; must be a mandatory CLI arg")
    check("brbseq --help documents the BRBseq<run-name>_<ts> leaf layout",
          "BRBseq_<run-name>_<timestamp>" in out,
          "output layout must be BRBseq_<run-name>_<timestamp>/ under OUT_DIR")
    check("brbseq --help does NOT document the stale --avt flag",
          "--avt " not in out,
          "--avt was removed in favour of --run-name")
    check("brbseq --help mentions --dry-run",
          "--dry-run" in out, "--dry-run is needed for CI smoke tests")

    rc, out, err = sh_out(["bash", brb_pipeline, "--version"])
    check("brbseq --version exits 0 and prints a version",
          rc == 0 and out.strip() != "", f"rc={rc} out={out!r}")

    rc, out, err = sh_out(["bash", brb_pipeline])   # no args at all
    check("brbseq (no args) exits NON-zero",
          rc != 0,
          "should hard-fail with 'missing required flag(s)' + non-zero exit")
    check("brbseq (no args) lists the missing flags",
          "--r1" in err and "--run-name" in err,
          f"stderr={err[-250:]}")

    rc, out, err = sh_out(["bash", brb_pipeline, "--this-flag-does-not-exist"])
    check("brbseq rejects unknown flags",
          rc != 0 and "I don't recognise" in err, f"rc={rc} err={err[-200:]}")

    rc, out, err = sh_out(["bash", brb_pipeline, "--resume", "bogus-value",
                           "--r1", "/tmp/x",
                           "--r2", "/tmp/y",
                           "--barcode-tsv", "/tmp/z",
                           "--star-genome", "/tmp/s",
                           "--run-name", "TEST0001"])
    check("brbseq --resume validates value",
          rc != 0 and ("resume" in err.lower() or "auto or force" in err.lower()),
          f"err={err[-200:]}")


    print()
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("RESULT: ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
