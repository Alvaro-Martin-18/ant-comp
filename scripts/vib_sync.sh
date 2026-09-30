#!/usr/bin/env bash
# Push the local working tree to the vib workspace. This is what makes a local edit
# to tools/ or activation/ live on the cluster without a commit or a push — the
# equivalent of Modal shipping your working tree into the image.
#
# Run from the repo root, before every submit:
#
#     bash scripts/vib_sync.sh
#
# `outputs/` and `.venv/` live inside the workspace and must never be synced.
# --exclude alone already protects them from --delete (verified), but the `protect`
# filters are kept because they also hold under --delete-excluded, which --exclude
# does not. NEVER add --delete-excluded: it deletes them outright, and outputs/ is
# every result you have.

set -euo pipefail

cd "$(dirname "$0")/.."

set -a
# shellcheck disable=SC1091
. ./.env
set +a

: "${SAPIA_VIB_HOST:?set SAPIA_VIB_HOST in .env}"
: "${SAPIA_VIB_WORKSPACE:?set SAPIA_VIB_WORKSPACE in .env}"

rsync -az --delete \
  --filter='protect outputs/' --filter='protect .venv/' \
  --exclude='.venv' --exclude='outputs' \
  --exclude='.git' --exclude='results_shard_0' \
  --exclude='__pycache__' --exclude='.ruff_cache' \
  -e "ssh -o BatchMode=yes -o ConnectTimeout=20 -x" \
  "$@" \
  ./ "$SAPIA_VIB_HOST:$SAPIA_VIB_WORKSPACE/"
