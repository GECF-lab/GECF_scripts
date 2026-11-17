#!/usr/bin/env bash
set -euo pipefail
# setup_basecaller.sh - create the 'basecaller' conda env (if it doesn't
# already exist), then install/update dorado INSIDE that env.
# Usage: ./setup_basecaller.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="basecaller"
ENV_YML="${SCRIPT_DIR}/basecaller_env.yml"
DORADO_SH="${SCRIPT_DIR}/update_dorado.sh"

# --- locate conda ---
CONDA_BIN="$(command -v conda 2>/dev/null || true)"
for p in "$HOME/miniconda3/bin/conda" \
         "$HOME/anaconda3/bin/conda" \
         /home/gecf/miniconda3/bin/conda; do
  [[ -z "$CONDA_BIN" && -x "$p" ]] && CONDA_BIN="$p"
done
if [[ -z "$CONDA_BIN" ]]; then
  echo "ERROR: conda not found in PATH or in common install locations." >&2
  exit 1
fi
echo "Using conda: $CONDA_BIN"

# --- sanity checks ---
if [[ ! -f "$ENV_YML" ]]; then
  echo "ERROR: env file not found: $ENV_YML" >&2
  exit 1
fi
if [[ ! -f "$DORADO_SH" ]]; then
  echo "ERROR: dorado update script not found: $DORADO_SH" >&2
  exit 1
fi

# --- step 1: create the conda env if it doesn't already exist ---
if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Conda env '$ENV_NAME' already exists - skipping creation."
else
  echo "Creating conda env '$ENV_NAME' from $ENV_YML ..."
  "$CONDA_BIN" env create -f "$ENV_YML"
fi

# --- resolve the env prefix (works regardless of conda quirks) ---
ENV_PREFIX="$("$CONDA_BIN" env list | awk -v e="$ENV_NAME" '$1==e {print $2}')"
if [[ -z "$ENV_PREFIX" || ! -d "$ENV_PREFIX" ]]; then
  echo "ERROR: could not resolve prefix of conda env '$ENV_NAME'." >&2
  exit 1
fi
echo "Using env prefix: $ENV_PREFIX"

# --- step 2: install/update dorado inside the env ---
# 'conda activate' is not used on purpose: update_dorado.sh only needs
# CONDA_PREFIX, CONDA_DEFAULT_ENV and the env bin first in PATH, so we set
# them directly. The tarball is downloaded/removed in a temp dir, never in
# your current working directory.
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"
export CONDA_PREFIX="$ENV_PREFIX"
export CONDA_DEFAULT_ENV="$ENV_NAME"
export PATH="$ENV_PREFIX/bin:$PATH"
bash "$DORADO_SH"

echo "Done. Activate the env ('conda activate $ENV_NAME') - dorado is on its PATH."
