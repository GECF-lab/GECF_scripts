#!/usr/bin/env python3.12
"""
Deterministic smoke test for generate_cellranger_template_arc.py (Section B parity).

Covers (no real mounts, no interactive prompt, no network):
  [1] py_compile of the generator module
  [2] import + key symbols present
  [3] folder_safe() identity / safety cases
  [4] TSO helper body inlined: tso_percent() present, standalone runner excluded
  [5] parse_xlsx on a real fixture .xlsx -> (sample, fid, avt) triples + AVT fallback
  [6] build_script() end-to-end -> emitted .sh content assertions
  [7] bash -n syntax check on the emitted .sh

Run with:  python3.12 arc_smoke_test.py
Exit code 0 = all passed; 1 = failure (assertion message printed).
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

GEN = os.path.join(os.path.dirname(__file__),
                   "arc", "generate_cellranger_template_arc.py")

PASS = []
FAIL = []


def check(name, cond, detail=""):
    if cond:
        PASS.append(name)
        print("  ok   %s" % name)
    else:
        FAIL.append(name)
        print("  FAIL %s%s" % (name, ("  [%s]" % detail) if detail else ""))


def main():
    print("[1] py_compile")
    r = subprocess.run([sys.executable, "-m", "py_compile", GEN],
                       capture_output=True, text=True)
    check("py_compile rc==0", r.returncode == 0, (r.stderr or "").strip())

    print("[2] import module + symbols")
    spec = importlib.util.spec_from_file_location("arcgen", GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for sym in ("folder_safe", "_load_tso_helper_body", "_TSO_HELPER_BODY",
                "_row_avt", "parse_xlsx", "group_libraries", "build_script",
                "SCRIPT_VERSION"):
        check("has %s" % sym, hasattr(mod, sym))
    check("SCRIPT_VERSION set", bool(getattr(mod, "SCRIPT_VERSION", None)),
          str(getattr(mod, "SCRIPT_VERSION", None)))

    print("[3] folder_safe()  (must match GEX's folder_safe exactly)")
    check("Cèsar->Cesar", mod.folder_safe("Cèsar") == "Cesar", mod.folder_safe("Cèsar"))
    check("a.b/c d->a.b_c_d", mod.folder_safe("a.b/c d") == "a.b_c_d", mod.folder_safe("a.b/c d"))
    check("whitespace->''", mod.folder_safe("   ") == "")
    check("weird!@#->weird (runs collapse, trailing _ stripped)", mod.folder_safe("weird!@#") == "weird", mod.folder_safe("weird!@#"))
    check("a..b//c->a..b_c", mod.folder_safe("a..b//c") == "a..b_c", mod.folder_safe("a..b//c"))
    check("'  pad  '->pad", mod.folder_safe("  pad  ") == "pad", mod.folder_safe("  pad  "))
    check("user.name-x preserved", mod.folder_safe("user.name-x") == "user.name-x", mod.folder_safe("user.name-x"))
    check("日本語->'' (non-ASCII dropped)", mod.folder_safe("日本語") == "")
    check("x/y/z->x_y_z", mod.folder_safe("x/y/z") == "x_y_z", mod.folder_safe("x/y/z"))
    check("None->''", mod.folder_safe(None) == "")

    print("[4] TSO helper body")
    body = mod._TSO_HELPER_BODY
    check("tso_percent() defined", "tso_percent()" in body or "tso_percent (" in body)
    check("has tso_percent fn", "tso_percent" in body)
    check("standalone runner EXCLUDED", 'if [ "$#' not in body and 'tso_percent "$1"' not in body)
    check("body non-trivial", len(body) > 200, "len=%d" % len(body))

    print("[5] parse_xlsx on fixture")
    import openpyxl
    from openpyxl import Workbook
    tmpdir = tempfile.mkdtemp(prefix="arcsmt_")
    xlsx = os.path.join(tmpdir, "AVT0551_submission.xlsx")
    wb = Workbook()
    ws = wb.active
    ws.title = "sheet1"
    # header row (row 3) with the two required labels
    ws.cell(row=3, column=1, value="Sample name")
    ws.cell(row=3, column=2, value="GECF fastq ID (avid0001)")
    ws.cell(row=3, column=3, value="Run / AVT name")   # adjacency target for AVT
    # data rows (4..7)
    rows = [
        ("sol001_GEX",  "avid0001", "AVT0551"),   # per-row AVT in col 3
        ("sol001_ATAC", "avid0002", ""),          # AVT inherited via sheet fallback
        ("sol002_GEX",  "avid0003", ""),          # AVT inherited via sheet fallback
        ("sol003_GEX",  "avid0004", ""),          # AVT inherited via sheet fallback
    ]
    for i, (s, f, a) in enumerate(rows):
        rr = 4 + i
        ws.cell(row=rr, column=1, value=s)
        ws.cell(row=rr, column=2, value=f)
        ws.cell(row=rr, column=3, value=a)
    wb.save(xlsx)

    parsed, detected = mod.parse_xlsx(xlsx)
    check("parsed 4 rows", len(parsed) == len(rows), "len=%d" % len(parsed))
    check("detected AVT==AVT0551", (detected or "").upper() == "AVT0551", str(detected))
    # every row's avt resolved to AVT0551 (per-row OR sheet fallback)
    all_avt = all((r[2] or "").upper() == "AVT0551" for r in parsed)
    check("all rows avt==AVT0551", all_avt, str(parsed))
    # triple shape
    check("rows are (sample,fid,avt) triples",
          all(len(r) == 3 for r in parsed), str(parsed))

    print("[6] build_script end-to-end")
    from collections import OrderedDict
    groups = OrderedDict()
    # sample 1: GEX + ATAC (GEX-leg TSO% expected)
    groups["sol001"] = [
        ("sol001_GEX", "avid0001", "Gene Expression"),
        ("sol001_ATAC", "avid0002", "Chromatin Accessibility"),
    ]
    # sample 2: ATAC-only (TSO% skip expected)
    groups["sol002"] = [
        ("sol002_ATAC", "avid0003", "Chromatin Accessibility"),
    ]
    stub_resolve = lambda base, fid: "/fake/%s/%s/R2" % (base, fid)
    out = mod.build_script(
        reference_genome="/fake/ref/GRCa-mmm9",
        groups=groups,
        avt_run="AVT0551",
        user="Cèsar",                 # non-ASCII to exercise folder_safe in .sh
        run_dir="/fake/AVT0551",
        source_files=[xlsx],
        resolve=stub_resolve,
    )
    shpath = os.path.join(tmpdir, "cellranger-arc_smoke.sh")
    with open(shpath, "w") as fh:
        fh.write(out)

    check("scriptversion stamp matches SCRIPT_VERSION",
          ("scriptversion=" + mod.SCRIPT_VERSION) in out, mod.SCRIPT_VERSION)
    check("user= folder_safe applied (Cesar)", "user='Cesar'" in out,
          "expected user='Cesar'")
    check("TSO helper inlined (tso_percent)", "tso_percent" in out and "tso_percent()" in out)
    check("tso_report_parts mkdir", 'mkdir -p "${Server_folder}/tso_report_parts"' in out)
    check("GEX-leg TSO% tee (sol001)", 'tso_percent "${_tso_fq}"' in out
          and "tso_report_parts/01_sol001__avid0001.txt" in out,
          "01_sol001__avid0001.txt per-avid part name")
    check("ATAC-only TSO% skip (sol002)", "no GEX leg in this sample; TSO% step skipped (ATAC-only)" in out)
    check("single tso_report.txt assembly", '} > "${_tso_out}"' in out and "tso_report.txt" in out)
    check("tso_report rsync to archive", 'rsync -avh "${_tso_out}" "${archive_dir}/"' in out)
    check("merged summary_metrics.csv (Samplename)", "Samplename,${_sm_hdr}" in out
          and "summary_metrics.csv" in out)
    check("ref-existence check present",
          "Reference genome ${REFPATH} does not exist. Aborting." in out
          and 'if [[ ! -d "${REFPATH}" ]]' in out)
    check("status sweep _overall_rc", "_overall_rc=0" in out and 'exit "${_overall_rc}"' in out)
    check("dispatch serial/pool present", "arc_sample_1" in out and "arc_sample_2" in out)

    print("[7] bash -n on emitted .sh")
    rn = subprocess.run(["bash", "-n", shpath], capture_output=True, text=True)
    check("bash -n rc==0", rn.returncode == 0, (rn.stderr or "").strip()[:300])

    print()
    print("passed=%d  failed=%d" % (len(PASS), len(FAIL)))
    if FAIL:
        print("FAILURES:")
        for f in FAIL:
            print("  - %s" % f)
        return 1
    # keep the generated .sh path in output for manual review
    print("emitted script: %s" % shpath)
    return 0


if __name__ == "__main__":
    sys.exit(main())
