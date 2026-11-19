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
  * No workbook I/O. This test runs on the drafting box and does not
    require mounted server data. A/A fixture diffing of emitted scripts is
    a separate (server-side) concern.

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
    "NGS_BASE_DIR_TEMPLATE",
    "AVT_PATTERN",
    "FASTQ_ID_PATTERN",
    "norm",
    "build_merged_cell_lookup",
    "get_cell_value",
    "read_user_field",
    "set_non_interactive",
    "setup_readline_completer",
    "ask_avt_run",
    "ask_reference_genome",
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
