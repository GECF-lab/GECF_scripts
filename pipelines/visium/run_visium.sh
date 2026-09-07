#!/usr/bin/env bash
#
# run_visium.sh — single entry point for the Visium HD SpaceRanger pipeline.
# Parse + validate every sample row, preview the spaceranger commands, ask
# for one confirmation, write the run script into the CURRENT directory
# (and optionally --launch it, follow the log, report PASS/FAIL).
# Details: bash run_visium.sh (no args).

set -uo pipefail

# ---------------------------------------------------------------------------
# Argument split: positional = workbook (output dir = current directory);
# the rest are forwards.
# --yes is consumed here; every other -*/token goes to the generator.
# Flags that take a value have their value consumed as well to avoid
# misclassifying it as a second positional (e.g. --species human would
# otherwise count "human" as the output directory).
# ---------------------------------------------------------------------------
ASSUME_YES=0
LAUNCH=0
POSITIONAL=()
GEN_ARGS=()

FLAGS_TAKING_VALUE=(
    "--species" "--chemistry" "--probe-root" "--slide" "--area"
    "--cytaimage" "--image" "--unknown-slide" "--loupe" "--fastqs"
    "--max-jobs" "--create-bam" "--spaceranger-bin"
)

takes_value() {
    local f="$1" want
    for want in "${FLAGS_TAKING_VALUE[@]}"; do
        [ "$f" = "$want" ] && return 0
    done
    return 1
}

add_positional() {  # append to the positional list (max 1 consumed later)
    if [ "${#POSITIONAL[@]}" -lt 1 ]; then
        POSITIONAL+=("$1")
    else
        echo "ERROR: unexpected positional argument: $1 (expected: just the workbook; the output directory is the current directory)" >&2
        exit 2
    fi
}

# The wrapper owns the output location (--output-dir is the second positional
# here, and it assigns --output-name itself). Strip any the caller passed
# through and mention it, rather than silently double-setting.
DROPPED_OUT=0
# Collect all args into an indexed array so we can look one ahead cleanly.
ARGS=( "$@" )
n=${#ARGS[@]}
i=0
while [ "$i" -lt "$n" ]; do
    a="${ARGS[$i]}"
    case "$a" in
        --output-dir|--output-name)
            DROPPED_OUT=1
            i=$((i+2)); continue
            ;;
        --yes|-y) ASSUME_YES=1; i=$((i+1)); continue ;;
        --launch) LAUNCH=1; i=$((i+1)); continue ;;
        -*)
            if takes_value "$a" && [ $((i+1)) -lt "$n" ]; then
                GEN_ARGS+=( "$a" "${ARGS[$((i+1))]}" )
                i=$((i+2))
            else
                GEN_ARGS+=( "$a" )
                i=$((i+1))
            fi
            ;;
        *) add_positional "$a"; i=$((i+1)) ;;
    esac
done
[ "$DROPPED_OUT" -eq 1 ] && echo "NOTE: ignoring --output-dir/--output-name — the wrapper sets the output location from its own arguments." >&2

usage() {
    cat >&2 <<'EOF'
Usage:
  bash run_visium.sh <workbook.xlsx> [run parameters ...]

The output directory is always the current directory — cd into the
folder where you want the run script and log written.

Wrapper options:
  --yes, -y    skip the y/n confirmation (non-interactive)
  --launch     also launch the run script after writing it (nohup) and
               follow the log until PASS/FAIL (default: create-only)

All other parameters go straight to the generator, e.g.
  --species, --image, --cytaimage, --slide, --area, --chemistry, --loupe
  (full list: python generate_spaceRanger_visium.py --help)

Example:
  cd /data/runs/AVT0598
  bash /path/to/run_visium.sh /data/samples.xlsx \
      --species human --cytaimage /data/images/cyt.tif \
      --image /data/images/20x.ome.tiff --slide V10J25-015 --area A1
EOF
    exit "${1:-2}"
}

if [ "${#POSITIONAL[@]}" -ne 1 ]; then
    echo "ERROR: need exactly one positional argument: <workbook.xlsx> (output goes to the current directory)" >&2
    usage 2
fi

WORKBOOK="${POSITIONAL[0]}"
OUT="$(pwd)"

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GENERATOR="$SCRIPT_DIR/generate_spaceRanger_visium.py"

if [ ! -f "$WORKBOOK" ]; then
    echo "ERROR: workbook not found: $WORKBOOK" >&2
    exit 2
fi
if [ ! -f "$GENERATOR" ]; then
    echo "ERROR: generator not found next to this script: $GENERATOR" >&2
    exit 2
fi

# Pick a python that can run the generator (openpyxl required).
pick_python() {
    local c
    for c in \
        /home/gecf/miniconda3/envs/test_env/bin/python \
        "${PYTHON:-}" \
        "$(command -v python3 || true)" \
        "$(command -v python || true)"; do
        [ -n "$c" ] || continue
        if "$c" -c "import openpyxl" >/dev/null 2>&1; then
            printf '%s\n' "$c"
            return 0
        fi
    done
    return 1
}

PY="$(pick_python)" || {
    echo "ERROR: no usable python 3 with openpyxl found. Install it:" >&2
    echo "       pip3 install openpyxl" >&2
    exit 2
}

mkdir -p "$OUT"
RUN_NAME="$(basename "${WORKBOOK%.*}")"
TS="$(date +%Y%m%d_%H%M%S)"
LAUNCH_LOG="$OUT/run_${TS}.log"
GEN_SCRIPT="$OUT/run_${TS}.sh"

echo
echo "======================================================================"
echo " Visium HD SpaceRanger run  —  $RUN_NAME"
echo " Workbook : $WORKBOOK"
echo " Output   : $OUT"
echo " Python   : $PY"
echo "======================================================================"
echo

# ---------------------------------------------------------------------------
# Step 1 — preview the exact commands (validates every sample row)
# ---------------------------------------------------------------------------
echo ">> 1/3  Checking the samples (validates workbook, images, references)"
set +e
PREVIEW_OUT="$( "$PY" "$GENERATOR" "$WORKBOOK" --output-dir "$OUT" \
    --output-name "run_${TS}" "${GEN_ARGS[@]+${GEN_ARGS[@]}}" \
    --print-commands 2>&1 )"
PREVIEW_RC=$?
[ -n "$PREVIEW_OUT" ] && echo "$PREVIEW_OUT"
if [ "$PREVIEW_RC" -ne 0 ]; then
    cat <<EOF

======================================================================
 FAIL — the workbook did not pass validation.
 Nothing was launched.

 Fix the highlighted row/column (or the run parameters) and run this
 command again.
======================================================================
EOF
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 2 — confirm
# ---------------------------------------------------------------------------
echo
echo ">> 2/3  Continue to write the run script?"
if [ "$ASSUME_YES" -eq 1 ]; then
    echo "   (--yes given: skipping confirmation)"
else
    printf "   Enter 'y' to continue, anything else to abort: "
    read -r answer
    case "$answer" in
        y|Y|yes|YES|Yes) : ;;
        *)
            echo "Aborted — nothing was written or launched."
            exit 0 ;;
    esac
fi

# ---------------------------------------------------------------------------
# Step 3 — write the generator script
# ---------------------------------------------------------------------------
echo
echo ">> 3/3  Writing the run script"
set +e
"$PY" "$GENERATOR" "$WORKBOOK" --output-dir "$OUT" \
    --output-name "run_${TS}" "${GEN_ARGS[@]+${GEN_ARGS[@]}}" > "$LAUNCH_LOG" 2>&1
GEN_RC=$?
if [ "$GEN_RC" -ne 0 ] || [ ! -f "$GEN_SCRIPT" ]; then
    cat <<EOF

======================================================================
 FAIL — the generator script could not be written.
 Last log lines:
$(tail -n 15 "$LAUNCH_LOG" 2>/dev/null)
======================================================================
EOF
    exit 1
fi
echo "   Runner written: $GEN_SCRIPT"

# ---------------------------------------------------------------------------
# Default: stop here. The run script is ready for you to launch manually.
# ---------------------------------------------------------------------------
if [ "$LAUNCH" -ne 1 ]; then
    echo
    echo "----------------------------------------------------------------------"
    echo " DONE — run script created; nothing was launched."
    echo "      Run script : $GEN_SCRIPT"
    echo "      Log        : $LAUNCH_LOG"
    echo
    echo "   When you're ready, launch it with:"
    echo "     nohup bash '$GEN_SCRIPT' >> '$LAUNCH_LOG' 2>&1 &"
    echo
    echo "   (or: bash run_visium.sh ... --launch to do it all at once)"
    echo "----------------------------------------------------------------------"
    exit 0
fi

# ---------------------------------------------------------------------------
# --launch — launch the run script, follow its log, report the verdict
# ---------------------------------------------------------------------------
echo
echo ">> Launching SpaceRanger"
set +e
nohup bash "$GEN_SCRIPT" >> "$LAUNCH_LOG" 2>&1 &
RUN_PID=$!
echo "   SpaceRanger launched (pid $RUN_PID) — following the log"
echo "   (to keep watching it later:  tail -f '$LAUNCH_LOG')"
echo

# ---------------------------------------------------------------------------
# Final status
# ---------------------------------------------------------------------------
echo ">> Final status"
# Stream the run log live while the run executes.
tail --lines=0 -f --pid="$RUN_PID" "$LAUNCH_LOG" || true
wait "$RUN_PID"
RUN_RC=$?

echo
echo "----------------------------------------------------------------------"
if [ "$RUN_RC" -eq 0 ]; then
    echo " PASS — all samples finished and the run deployed cleanly."
    echo "      Samples   : $OUT/<sample>/outs"
    echo "      Summaries : $OUT/Summaries"
    echo "      Log       : $LAUNCH_LOG"
    echo
    echo "   Open per-sample results in a browser:"
    echo "     $OUT/Summaries/web_summaries/"
else
    echo " FAIL — the run finished with errors (exit $RUN_RC)."
    echo "      Log       : $LAUNCH_LOG"
    echo
    echo "   Recent log lines:"
    tail -n 25 "$LAUNCH_LOG" 2>/dev/null
fi
echo "----------------------------------------------------------------------"

exit "$RUN_RC"
