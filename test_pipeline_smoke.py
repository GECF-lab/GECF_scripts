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
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
PIPELINES_ROOT = os.path.join(REPO, "pipelines")

FILES = [
    "pipelines/common/intake.py",
    "pipelines/arc/generate_cellranger_template_arc.py",
    "pipelines/atac/generate_cellranger_template_atac.py",
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

# Emitted by the arc generator only (cellranger-arc count --localcores/
# --localmem). The atac generator routes resource tuning through the
# $CUSTOM_OPTIONS prompt instead, so these are arc-only atoms.
ARC_ONLY_ATOMS = [
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


def main() -> int:
    print("== 1. py_compile ==")
    for rel in FILES:
        path = os.path.join(REPO, rel)
        try:
            py_compile.compile(path, doraise=True)
            check(f"py_compile {rel}", True)
        except Exception as e:  # noqa: BLE001
            check(f"py_compile {rel}", False, str(e))

    # Import common.intake the normal package way (pipelines/ on sys.path)
    # so identity checks run against the ONE real module object that the
    # generators import via 'from common.intake import ...'.
    if PIPELINES_ROOT not in sys.path:
        sys.path.insert(0, PIPELINES_ROOT)
    import common.intake as common

    print("\n== 2. common.intake exports all shared atoms ==")
    for name in SHARED_ATOMS + ARC_ONLY_ATOMS:
        check(f"common.intake.{name} exists", hasattr(common, name))

    print("\n== 3. arc generator ==")
    arc = load_module_by_path(
        os.path.join(REPO, "pipelines/arc/generate_cellranger_template_arc.py"),
        "arc_gen_test",
    )
    for name in SHARED_ATOMS + ARC_ONLY_ATOMS:
        check(f"arc.{name} is common.intake.{name}",
              getattr(arc, name, None) is getattr(common, name, None),
              "local redefinition detected")
    check("arc.parse_xlsx callable", callable(getattr(arc, "parse_xlsx", None)))
    check("arc.build_script callable",
          callable(getattr(arc, "build_script", None)))
    check("arc.main callable", callable(getattr(arc, "main", None)))

    print("\n== 4. atac generator ==")
    atac = load_module_by_path(
        os.path.join(REPO, "pipelines/atac/generate_cellranger_template_atac.py"),
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
