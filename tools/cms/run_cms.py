#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) scoring the binder/target interface of each design
with contact molecular surface (CMS) and shape complementarity (SC), on the GPU.

CMS (Longxing Cao / Brian Coventry) is a distance-weighted interface area on the
target, ``sum(area * exp(-0.5 * d**2))`` over the target's buried surface dots: it
rewards an interface that is both large and tight, and is the interface metric the
Baker-lab binder pipelines filter on. SC (Lawrence & Colman 1993) comes off the same
molecular surface for free. The numbers come from cms-cuda, which matches the
original NumPy code to ~1e-13 and is 36-165x faster.

``action: update`` -- an interface score is a property of a design that already
exists, so it annotates the table in place. Run it on whichever structure column you
want scored (``-i``): ``rfdiffusion3_path`` for the designed pose, ``boltz_path`` /
``alphafold3_path`` for the predicted one, a relaxed structure for the refined one.
Label each with ``-l`` so the columns do not collide.

``default_input_column`` is the ``"not applicable"`` sentinel (as in usalign and
chainsel): there is no honest default structure column, and a wrong one fails
silently, so the builder refuses the run unless ``-i`` names a column the table has.

Several designs are packed into one task (``--designs-per-task``): the GPU work is
well under a second per design, while a container cold start and the first-call
kernel compilation take ~1 minute (measured on an L4), so a task per design would spend almost all its time
starting up.

Usage:
    sapia run cms outputs/20260928_100140_7ojg_binder \
        --table table1 \
        --input-column boltz_path \
        --binder-chains C --target-chains A,B \
        --dir-label pred
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path

# The sentinel default_input_column (see the module docstring).
NO_DEFAULT_COLUMN = "not applicable"
DEVICES = ("auto", "cuda", "cpu")


class CmsArgs(CommonArgs):
    binder_chains: str
    target_chains: str
    exclude_resnames: str
    no_sc: bool
    max_cms: bool
    device: str
    designs_per_task: int


def add_run_cms_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--binder-chains",
        type=str,
        required=True,
        help="The binder's chain IDs, comma-joined (e.g. 'C', or 'H,L' for an "
        "antibody). A chain named here but ABSENT from a design is an error for that "
        "design, never a smaller selection.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        required=True,
        help="The target's chain IDs, comma-joined (e.g. 'A,B'). CMS is reported on "
        "the target side by convention (cms_target); the binder side is cms_binder.",
    )
    parser.add_argument(
        "--exclude-resnames",
        type=str,
        default="",
        help="Residue names to drop from both sides, comma-joined (e.g. lipids or a "
        "ligand you do not want counted as interface). Waters, single-atom ions and "
        "hydrogens are always dropped. Default: none.",
    )
    parser.add_argument(
        "--no-sc",
        action="store_true",
        help="Skip shape complementarity. It shares the CMS surface, so it costs "
        "little; the default computes it.",
    )
    parser.add_argument(
        "--max-cms",
        action="store_true",
        help="Also compute each side's maximum possible CMS (its whole molecular "
        "surface) and the fraction in contact (cms_frac_target / cms_frac_binder). "
        "Off by default: two extra surface builds, and the fraction mostly matters "
        "for small-molecule targets.",
    )
    parser.add_argument(
        "--device",
        type=str,
        choices=DEVICES,
        default="auto",
        help="'auto' (default): cuda when the task has a GPU, cpu under -g 0. 'cuda' "
        "never falls back -- a task without a usable GPU records an error for every "
        "design instead of silently running ~100x slower.",
    )
    parser.add_argument(
        "--designs-per-task",
        type=int,
        default=100,
        help="Designs scored per task (default 100). On an L4 each task pays ~60 s "
        "of one-time kernel compilation, then ~0.6 s per design. Lower it for more parallelism.",
    )


def _split(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def build_cms_manifest(ctx: ManifestCtx[CmsArgs]) -> list[tuple[str, ...]]:
    binder = _split(ctx.args.binder_chains)
    target = _split(ctx.args.target_chains)
    if not binder or not target:
        raise ValueError(
            "--binder-chains and --target-chains must both name at least one chain."
        )
    overlap = sorted(set(binder) & set(target))
    if overlap:
        raise ValueError(
            f"chains {overlap} are in both --binder-chains and --target-chains; an "
            f"interface needs two disjoint sides."
        )
    if ctx.args.designs_per_task < 1:
        raise ValueError("--designs-per-task must be >= 1.")

    device = ctx.args.device
    if device == "auto":
        device = "cpu" if ctx.args.gpus_per_task == 0 else "cuda"
    elif device == "cpu":
        ctx.args.gpus_per_task = 0  # do not pay for a GPU the task will not use

    column = ctx.args.input_column
    if column == NO_DEFAULT_COLUMN or column not in ctx.df.columns:
        available = ", ".join(
            str(c) for c in ctx.df.columns if str(c).endswith("_path")
        )
        raise ValueError(
            f"cms has no default input column: pass -i/--input-column with the "
            f"structure column to score (got {column!r}, which table "
            f"'{ctx.args.table}' does not have). Structure columns available: "
            f"{available or '(none)'}."
        )

    ready = ctx.ready
    members: list[tuple[str, str]] = []
    for name in ready.index:
        name = cast(str, name)
        src = Path(str(ready.at[name, column]))
        if not src.exists():
            print(f"{name}: MISSING {src} (skipping)")
            continue
        # No CIF->PDB staging: Biopython reads either format.
        members.append((name, str(volume_path(src))))

    # One sub-manifest per task, and a top-level row per task pointing at it.
    tasks_dir = ctx.out_dir / "cms_tasks"
    tasks_dir.mkdir(parents=True, exist_ok=True)
    x = ctx.args.designs_per_task
    manifest_rows: list[tuple[str, ...]] = []
    for t, i in enumerate(range(0, len(members), x)):
        task_file = tasks_dir / f"task_{t}.tsv"
        with open(task_file, "w") as f:
            for name, src in members[i : i + x]:
                f.write(f"{name}\t{src}\n")
        manifest_rows.append(
            (
                str(volume_path(task_file)),
                ",".join(binder),
                ",".join(target),
                # May be empty, so never last.
                ",".join(_split(ctx.args.exclude_resnames)),
                "nosc" if ctx.args.no_sc else "sc",
                "max" if ctx.args.max_cms else "nomax",
                # Kept last: always non-empty.
                device,
            )
        )

    return manifest_rows
