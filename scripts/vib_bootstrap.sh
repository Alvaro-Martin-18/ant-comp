#!/usr/bin/env bash
# Stand up (or repair) the vib workspace: the rsync'd copy of this repo on cluster
# storage that every `sapia` command runs from.
#
# Run it from the repo root, once per person, and again whenever uv.lock changes:
#
#     bash scripts/vib_bootstrap.sh
#
# Deliberately uses NOTHING but what the cluster already provides — the Miniconda3
# module's Python and stdlib `venv`/`pip`. `uv` is not required anywhere on vib (on
# this cluster it is a per-user install under one person's home, so it cannot be
# assumed), and the only uv artefact used here is `uv.lock`, read locally to pin the
# same prosapia commit the Modal side runs.
#
# The build runs inside a SLURM allocation, not on the login node: the login node is
# heavily throttled and a venv build there wedges partway through.

set -euo pipefail

cd "$(dirname "$0")/.."

set -a
# shellcheck disable=SC1091
. ./.env
set +a

: "${SAPIA_VIB_HOST:?set SAPIA_VIB_HOST in .env}"
: "${SAPIA_VIB_WORKSPACE:?set SAPIA_VIB_WORKSPACE in .env}"
: "${SAPIA_VIB_ACCOUNT:?set SAPIA_VIB_ACCOUNT in .env — srun/sbatch are rejected without -A}"

SSH=(ssh -o BatchMode=yes -o ConnectTimeout=20 -x)

# Pin prosapia to whatever uv.lock resolved, so vib and Modal run the same library.
COMMIT="$(grep -oE 'prosapia\?branch=[a-zA-Z0-9._-]+#[0-9a-f]{40}' uv.lock | head -1 | cut -d'#' -f2)"
[ -n "$COMMIT" ] || { echo "could not read the prosapia commit from uv.lock" >&2; exit 1; }
echo "prosapia pinned to ${COMMIT:0:8} (from uv.lock)"

echo "==> creating $SAPIA_VIB_WORKSPACE"
"${SSH[@]}" "$SAPIA_VIB_HOST" "mkdir -p '$SAPIA_VIB_WORKSPACE'"

echo "==> syncing the working tree"
bash scripts/vib_sync.sh

echo "==> building the venv on a compute node (several minutes, no output until done)"
"${SSH[@]}" "$SAPIA_VIB_HOST" "srun -A '$SAPIA_VIB_ACCOUNT' -c 8 -t 2:0:0 --job-name=venv_bootstrap bash -lc '
    set -e
    cd \"$SAPIA_VIB_WORKSPACE\"
    module load Miniconda3
    echo \"building with \$(python -V) on \$(hostname)\"
    rm -rf .venv
    python -m venv .venv
    .venv/bin/pip install --quiet --upgrade pip setuptools wheel
    .venv/bin/pip install \"prosapia[modal] @ git+https://github.com/jlmoraleshellin/prosapia@$COMMIT\"
'"

echo "==> verifying"
"${SSH[@]}" "$SAPIA_VIB_HOST" "cd '$SAPIA_VIB_WORKSPACE' && $SAPIA_VIB_ACTIVATE && sapia run --help >/dev/null && echo 'sapia OK'"

echo "done."
