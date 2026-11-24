# CellRanger pipelines — platform-agnostic discovery (roadmap)

Working list, agreed 2026-11-20. Goal: sample discovery that works without the
Excel being the only entry point AND without being welded to the production
AVT path or the `avid` leaf prefix, so new sequencers (MSH/MSQ/NVQ seen on
this server; more coming) are absorbed without code changes.

## Evidence gathered (this server, 2026-11-20)

| platform | run name      | demux output shape (verified) |
|----------|---------------|-------------------------------|
| AVT      | `AVT0565`     | `fastq_bases2fastq/<batch>/Samples/<group>/avid#####/<Name>_S?_R1_001.fastq.gz` |
| MSH      | `MSH0024`     | recent runs only `bcl/` + `fastq_BCL-Convert/` (no demuxed leaf confirmed) |
| MSQ      | `MSQ0350metr` | `MSQ0350_<date>/<group>/msidNNNNN/GFPlow_S1_R1_001.fastq.gz` |
| NVQ      | `NVQ0317`     | `NVQ0317_<date>/<group>/nvidNNNNN/1_A6_S1_R1_001.fastq.gz` |

Invariants (safe to rely on, verified on all three demuxed platforms):
- the leaf is a directory named `<prefix>id<digits>` containing the fastqs;
- fastq filenames keep `<samplename>_S?_R1|_R2_(|I1|I2)_001.fastq.gz`;
- therefore the sample name is always parseable from the R1 basename.

Workbook-side note (2026-11-23): the MSH0024 submission workbook names its
fastq IDs `mhid#####` (sheet header "GECF fastq ID (MHIDxxxx)"), so both
`mhid` and `msid` appear in the wild for MSH/MSQ.

Variables (never assume): platform root/prefix, run-name shape, number of
intermediate folders (`fastq_bases2fastq/`, `Samples/`, batch ts, group folder
names).

## Change list

1. **Generalize the ID token** — `common/intake.py::FASTQ_ID_PATTERN` (single
   source of truth, used by arc/gex/atac row gate and the add-sample path
   validation) from `^avid\d{4,6}$` to `^[a-z]*id\d{4,6}$` (case-insensitive;
   prefix = letters only — user confirmed no digits before the `id`).
   Accepts `avid/msid/nvid/xvid/…`; stays shell-safe for ARC's
   `fastq_{idx}_{fid}` variable names.
   **DONE 2026-11-20** (verified: unit accept/reject checks + arc smoke suite).

2. **Content-based discovery + anchor** — replace `find_fastq_path(main_avt,
   avid)` with `discover_libraries(anchor)`: walk an explicit directory
   (run folder / platform dir / anything), return every directory directly
   containing an R1+R2 (optional I1/I2) fastq set; name parsed from the R1
   basename; leaf recognized via the generalized ID token. Intermediate
   folders stop mattering.
   **DONE 2026-11-23** — as a variant: `find_fastq_path_with_pref` now walks
   the run tree (`run_base_dir(run)` for ANY of AVT/MSH/MSQ/NVQ), collects
   up to 2 hits (early exit), disambiguates duplicate index-attempt dirs via
   `choose_sample_dir`, and prunes BCL/basecalls/QC/raw-flowcell branches
   (AVT0551 absent-id walk: 113 s → 14 s).  Verified on real trees:
   AVT0594, AVT0551 (2 index attempts), MSH0023 (2 demux attempts), MSH0024,
   NVQ0312.  MSQ016metr has no per-sample leaves (barcode-amplicon project)
   → resolver returns None → FIXME_LOCATE (correct, out of CellRanger scope).
   NOT done: the standalone `discover_libraries(anchor)` API of item 3.

3. **Excel becomes an overlay, not a prerequisite** — intake: `anchor`
   required-or-prompted (production root pre-filled), Excel optional.
   With Excel: resolve each row against the discovery map, Excel stays
   authoritative for `--id`; unresolvable row → existing add-sample full-path
   prompt; Excel-name-vs-filename mismatch → warn, don't fail. Without Excel:
   the discovery map is the sample list, fed to the existing keep/drop menu
   (`ask_sample_selection`).

4. **Fail-loud preflight** — zero library-shaped dirs under the anchor →
   stop with one actionable line (expected leaf/filename contract + the
   anchor given). No guessing, no FIXME soup.
   **PARTIAL 2026-11-23** — unresolvable fastq id now warns + emits
   `FIXME_LOCATE_<id>.fastq` (per row, all three generators) and a missing
   base dir prints a "not on this machine" note; a hard stop-on-zero-libraries
   preflight line is NOT implemented (rows still get FIXME leaves so the
   draft script stays generatable).

5. **Output-side defaults overridable at generation time** —
   `SCRATCH_ROOT`, `NGSruns_archive`, `reference_genome` settable by flags at
   generation (production values stay the defaults) so a new server's script
   is born with the right defaults.

## Verification bar (per change)

- `arc_smoke_test.py` stays green.
- New fixtures from real trees: NVQ0317-style layout (no `Samples/` layer,
  `nvid` leaves) discovers the same libraries the AVT logic found; a
  no-Excel run over an AVT tree produces the same `cellranger-atac count`
  lines as the Excel path for the same samples.

## Deliberately NOT doing

- No platform-profile/config-file system (speculative until an unseen layout
  forces it).
- No change to per-sample params, TSO%, web-summary/merge,
  scratch→archive→rm flow. (Header/sheet-selection parsing IS in scope when a
  real non-AVT workbook demands it — see addendum.)
- A demux output that has no `*idNNNN` leaves is a demux-pipeline problem;
  item 4 reports it, it is not papered over here.

## Addendum (2026-11-23) — real MSH workbook

Checked `MSH0024 Bing_GECF_copy.xlsx`. Before the fix, atac fell back to the
GUIDELINES legend sheet as its sample source (emitting `PhiX_i7`,
`i5_backbone`, …); arc/gex would have taken wrong cells or failed. Fixed:
the header gate now matches any `GECF fastq ID` header cell (no `avid`
required); `guideline`-named sheets are excluded as data sources in all three
generators; atac treats QC `0` cells as blank (arc/gex parity); gex takes the
first `Sample name` column (arc/atac parity — the MSH Design sheet repeats
the phrase in a later column). New fixture: `test_pipeline_smoke.py` §5c
(synthetic MSH workbook through atac/arc/gex). All three now parse all 13
`mhid` samples from the real file with zero legend contamination.

### Sheet inventory of the real AVT workbooks (all 12 in `AVT_dossiers/`)

Classified every sheet of every workbook by presence of a real
`GECF fastq ID` + `Sample name` header row and by content:

**Sample-data sheets (valid fallbacks, in preference order):**
1. `User_x1x_copy` — the user's own table (always preferred by the
   generators);
2. `Design` / `User_copy` (MSH workbooks — same role);
3. `whole_severaluser` — multi-user shared run table;
4. `Sheet1` (atac workbooks where the user table was renamed);
5. `Option Index` — the index-assignment view.  It carries *real* sample
   rows (same names + fastq IDs as the user table, e.g. all 13 `mhid`
   samples in MSH0024) plus a couple of control rows with empty fastq ID
   (`PhiX`, `High_diversity_library`), which are skipped with a visible
   "row skipped" warning.  Kept as a data candidate: it is genuine data,
   and it only ever wins if no primary sheet parses — which no real
   workbook does (primary sheets rank first and have >= rows).

**Never data (excluded or unreachable):**
- `GUIDELINES` — excluded by name (all three generators).
- `Index table MM.YY` — the control-library reference (PhiX_i7, PhiX_i5,
  i5_backbone, Ad1noMX + ~12k index names).  Excluded by name in all three
  generators (user decision 2026-11-23: "we don't use this as a fallback");
  it has no `fastq ID` header anyway, so the exclusion also shields the
  atac positional fallback from its 12 897 rows.
- `Spike_1_copy` / `Spike_2_copy` — real rows, but spike controls belong to
  a different pipeline input; excluded by the pre-existing `spike` rule.
- `lane splitting (same pool)`, `Expected sample conc after dil`,
  `Paste Metrics Table MultiQC`, `nofilter`, `Firefly`,
  `Option Trimming Detail` — pool/QC/layout notes; none exposes the
  `fastq ID` + `Sample name` header row, so none can be picked.

Verification: full `test_pipeline_smoke.py` passes (including §5c, which now
also carries the `Option Index` and `Index table 08.26` legends); all three
generators on the real MSH0024 file emit exactly the 13 `mhid` samples,
zero legend names in the generated scripts; control rows with an empty
fastq ID raise one visible per-row warning (intentional — no silent
skipping of user errors).

**Placeholder run name (decision 2026-11-23, option 2):** in non-interactive
mode, when the detection chain (per-row AVT cell → sheet title → file name)
finds no run name at all, `ask_avt_run` still returns the `AVT0000`
placeholder (so a draft script can be generated and hand-edited, like the
FIXME_LOCATE_ leaf pattern) — but it now prints a loud WARNING telling the
user to edit `avt_run=` before running and how to avoid the situation
(name the file `AVT####…xlsx` or run interactively).  Behavior otherwise
unchanged: detected names still win, `intake.py` is the single patch point
for all three generators.  Regression check: `test_pipeline_smoke.py` §5c
now asserts the warning appears (fixture file carries no run name) for all
three generators, and the real `AVT0594 Sonia Cruciani.xlsx` file (name
intact) still resolves cleanly to `avt_run='AVT0594'` with no warning.

## Addendum (2026-11-24) — add-sample prompt generalized (arc ZZ_10XG fix)

User report: prompting to add a sheet-less sample in the arc generator
prefilled `…/AVT/AV240401/ZZ_10XG/arc` — the rsync *archive* destination
(`NGSruns_ARCHIVE_DIR`), where input samples never live. atac had the same
bug (worse: unmounted `/AA/NGSruns/…`). All three call sites now pass
`NGS_RUNS_ROOT` (`/mnt/PTEGraw/AA/NGSruns`) — user decision (2026-11-24):
an added sample may come from ANY sequencer (AVT/MSH/MSQ/NVQ/…), so the
prefill stops at the shared root level and never guesses a sequencer or
run; relative tails complete against it (`AVT/…`, `MSH/…`, `NVQ/…`),
full paths from anywhere else are accepted verbatim. The prompt itself is sequencer-neutral: "Full path to the sample's
FASTQ folder" (any `FASTQ_ID_PATTERN` basename — avid/mhid/nvid/…);
bare ids and non-id basenames are still rejected. Verified: piped
prompt-flow tests (prefill tail + foreign abs path + all rejection paths)
and full smoke suite ALL CHECKS PASSED. `NGSruns_ARCHIVE_DIR` is kept only
for the rsync destination line in generated scripts.

## Addendum (2026-11-23) — platform-agnostic run intake, items 2–4 wired in

Run-name + fastq-path intake is now platform-agnostic across **AVT, MSH, MSQ,
NVQ** (Visium generators deliberately untouched — they keep their local
copies; user instruction 2026-11-23 "don't touch visium"):

- `common/intake.py`: new `KNOWN_RUN_CODES` + `RUN_BASE_SUBPATH` (AVT keeps
  its `AV240401/` extra level) and `run_base_dir(run)` resolving run name →
  base folder for every code (prefix-matches `MSQ016` → `MSQ016metr`).
  `AVT_PATTERN` now matches all four codes; `ask_avt_run`'s prompt/warning
  are code-agnostic.  `NGS_BASE_DIR_TEMPLATE` stays for backward
  compatibility (gex/atac/arc no longer import it — all switch to
  `run_base_dir`; `test_pipeline_smoke.py` asserts the new atoms are the
  shared objects).
- Disambiguation: runs with duplicate index-attempt directories
  (AVT0551 `1_wrong_index…`/`2_corrected_indexes…`, MSH0023
  `1_original_wrong_index`/`2_corrected_index_given`) resolve via a
  "collect ≤2 hits + `choose_sample_dir`" walk instead of first hit wins;
  BCL/basecalls/QC/raw-flowcell branches are pruned (113 s → 14 s on the
  worst case, AVT0551).
- §5c of the smoke test now asserts the **detected** MSH run: a workbook
  whose sheet names `MSH0099` yields `cellranger-{atac,arc,count}_MSH0099.sh`
  for all three generators, and a run-LESS workbook (no AVT/MSH/MSQ/NVQ
  code anywhere) yields the loud `AVT0000` placeholder behavior.
  Both paths: ALL PASS, full suite green (2026-11-23).

Evidence table correction (same date): MSH runs DO carry demuxed leaves —
`MSH/MSH0024/fastq_BCL-Convert/<attempt>/mhid#####/` (and likewise MSH0023
with two attempts and a "corrected index given" copy).  NVQ confirmed with
`NVQ/NVQ0312/…/nvid06797/`.

## Addendum — 2026-11-24 (archive root moved off the sequencer folder)

- **New shared archive root** `/mnt/PTEGraw/AA/NGSruns/ZZ_10XG/` (created by
  the user), mirroring the internal layout of the former
  `AVT/AV240401/ZZ_10XG/` (`aggr/ arc/ atac/ count/ flex/ multi/
  spaceranger/`).
- All remaining writers now deploy there regardless of sequencer: gex
  (`/count`), atac (`/atac`), arc (`/arc`); README table + standard-archive
  rule updated to match.
- Provenance: generator versions bumped for the destination change — gex
  1.1.7, atac 1.3.5, arc 1.3.4.  Smoke test now asserts each generator's
  `NGSruns_ARCHIVE_DIR` sits under the shared root.
- **Legacy pipeline code removed from the repo (2026-11-24)** per user
  decision: `cellranger_pipelines/multi/` (CRmulti driver),
  `cellranger_pipelines/visium/` (run_visium.sh + spaceRanger
  generator), and the legacy arc bash driver
  `arc/CRarc_patched_full`.  Server-side copies under
  `/home/gecf/NGS_pipelines/CR_pipeline/` are untouched; README sections for
  the three were removed (a "legacy drivers removed" note explains where
  the server-side workflows still live).
- Historical results under `AVT/AV240401/ZZ_10XG/` stay put; new results
  from all sequencers land in the shared tree.
- Note: *input* (FASTQ) trees are unchanged — still per-sequencer
  (`RUN_BASE_SUBPATH` below `NGS_RUNS_ROOT`).  Only the output destination
  was generalized.  This is a new static default, not item 5 (per-run
  overridable output defaults), which remains open.
