#!/usr/bin/env python3
"""
Collect BindCraft2 campaign results into the (child) binder table.

``run_bindcraft2.py`` points each campaign's ``project_folder`` at
``<out_dir>/campaigns/<name>/``, where BindCraft2 writes three stage folders, each
with a table this collector can read:

    campaigns/<name>/1_Trajectories/!_Trajectories.csv        one row per attempt
    campaigns/<name>/1_Trajectories/<design>/<design>_trajectory.cif
    campaigns/<name>/2_Refolded/!_Refolded.csv       every scored candidate + outcome
    campaigns/<name>/2_Refolded/Complexes/<design>.cif
    campaigns/<name>/3_Ranked/!_Ranked.csv          the accepted designs, best-first
    campaigns/<name>/3_Ranked/<design>.cif

Which one to read is **decided by the run, not guessed here**: a `--trajectory-only`
run fills `1_Trajectories` and never accepts anything, so `3_Ranked` stays empty. The
run records that in the out_dir sidecar and `--stage auto` (the default) follows it,
the same contract that keeps `input_column` consistent between the two phases. Pass
`--stage` explicitly only to override.

* ``trajectories`` -- the hallucinated **backbones**, before any sequence redesign.
  The rows to hand to `atomium` / `proteinmpnn`. Note their `sequence` is the one the
  gradient happened to land on, which is exactly what you are about to replace.
* ``refolded`` -- every candidate the campaign scored, rejects included, with
  ``outcome`` and ``failed_filters``. What to read when a campaign accepted nothing.
* ``ranked`` -- only what the campaign accepted.

A trajectory row's ``terminated`` column names the stage the attempt stopped at;
blank means it ran to completion. They are collected rather than dropped, so filter
on it (``-f``) instead of assuming every row is a usable backbone.

Rows are keyed by BindCraft2's own ``design`` identity rather than by position:
``!_Ranked.csv`` is re-sorted by ``i_pDAE`` after every acceptance, so an index into
it is not stable across a re-collect, while the design name (campaign, modality,
length and recipe hash) is. The submitter sets ``campaign_name`` to the parent row's
name, so a design name already starts with its parent; the prefix is only added here
when a --set override has changed that.

Each accepted complex is converted to PDB (downstream tools consume PDB) and becomes
the row's ``<leaf>_path``; the source mmCIF stays available as ``<leaf>_cif_path``.
Note the structure is the **complex**, binder + target, not the binder alone.

Safe to re-run: rows are rebuilt from the CSV and the structures on disk.

Usage:
    sapia collect bindcraft2 outputs/RUN --table <the table the run reserved>
    sapia collect bindcraft2 outputs/RUN --table <table> --stage refolded --force
"""

import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple

import pandas as pd

from prosapia.core import (
    RUN_META_FILENAME,
    CollectArgs,
    CollectCtx,
    CollectEach,
    Collected,
    DesignCtx,
)
from prosapia.utils import ensure_pdb

# Layout run_bindcraft2.py imposes on the out_dir -- keep the two in step.
CAMPAIGNS_DIRNAME = "campaigns"

# --stage default: read what the run recorded in the sidecar instead of guessing.
STAGE_AUTO = "auto"


def _flat_structure(search_dir: Path, design: str) -> Path | None:
    """``<search_dir>/<design>.cif``, or the target-suffixed variant a multi-target
    campaign writes. The free-binder prediction (``_monomer``) is a different
    structure and is never returned."""
    exact = search_dir / f"{design}.cif"
    if exact.is_file():
        return exact
    matches = sorted(
        p for p in search_dir.glob(f"{design}*.cif") if "_monomer" not in p.name
    )
    return matches[0] if matches else None


def _ranked_structure(stage_dir: Path, design: str) -> Path | None:
    return _flat_structure(stage_dir, design)


def _refolded_structure(stage_dir: Path, design: str) -> Path | None:
    return _flat_structure(stage_dir / "Complexes", design)


def _trajectory_structure(stage_dir: Path, design: str) -> Path | None:
    """``1_Trajectories/<design>/<design>_trajectory.cif`` -- a per-design subfolder,
    unlike the two later stages, and the filename carries the design name twice
    (BindCraft2's ``trajectory_output_path`` prefixes it)."""
    design_dir = stage_dir / design
    if not design_dir.is_dir():
        return None
    return _flat_structure(design_dir, f"{design}_trajectory")


class Stage(NamedTuple):
    folder: str  # under project_folder
    table: str  # its CSV, inside that folder
    find: Callable[[Path, str], Path | None]  # (stage_dir, design) -> structure


STAGES: dict[str, Stage] = {
    "trajectories": Stage("1_Trajectories", "!_Trajectories.csv", _trajectory_structure),
    "refolded": Stage("2_Refolded", "!_Refolded.csv", _refolded_structure),
    "ranked": Stage("3_Ranked", "!_Ranked.csv", _ranked_structure),
}

# The column carrying the design's identity, and the one carrying its sequence.
DESIGN_COL = "design"
SEQUENCE_COL = "Binder_Sequence"
# 1_Trajectories only: the stage the attempt stopped at, blank when it completed.
TERMINATED_COL = "terminated"

# Collected by default: the metrics a binder campaign is actually read on. Everything
# else in the CSV is reachable with --metrics or --all-metrics rather than being
# poured into the table by default -- BindCraft2 writes ~60 columns per design.
CORE_METRICS = (
    # confidence
    "i_pDAE",
    "i_pTM",
    "i_pAE",
    "pLDDT",
    "pTM",
    "Unbound_Binder_pLDDT",
    # interface
    "Interface_Residues",
    "Interface_BuriedArea",
    "Hotspot_Contact_Fraction",
    # developability
    "Surface_Hydrophobicity",
    "Binder_Length",
    "Binder_Net_Charge",
    "Binder_Free_Cysteines",
    # provenance / triage. Which of these exist depends on the stage: `outcome` and
    # `failed_filters` on refolded, `rank` on ranked, `terminated` and `autotuned`
    # on trajectories. Absent ones are simply not collected.
    "rank",
    "hash",
    "trajectory",
    "outcome",
    "failed_filters",
    "terminated",
    "autotuned",
    "length",
)


class BindCraft2CollectArgs(CollectArgs):
    stage: str
    metrics: str
    all_metrics: bool


def add_collect_bindcraft2_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--stage",
        choices=(STAGE_AUTO, *STAGES),
        default=STAGE_AUTO,
        help="Which campaign stage to collect. 'auto' (the default) reads what the "
        "run recorded: 'trajectories' after a --trajectory-only run, 'ranked' "
        "otherwise. 'trajectories' takes the hallucinated backbones before any "
        "sequence redesign; 'refolded' takes every scored candidate, rejects "
        "included, with its `outcome` and `failed_filters`; 'ranked' takes only what "
        "the campaign accepted. Override with -l/--dir-label to keep two stages of "
        "one run in separate columns.",
    )
    parser.add_argument(
        "--metrics",
        type=str,
        default="",
        help="Comma-separated extra CSV columns to collect on top of the core set "
        "(e.g. 'Binder_pI,Binder_Helix_Fraction,Interface_W_Count'). Names are "
        "BindCraft2's own, as spelled in !_Ranked.csv.",
    )
    parser.add_argument(
        "--all-metrics",
        action="store_true",
        help="Collect every column in the stage's CSV instead of the core set. "
        "Wide: BindCraft2 writes ~60 measurements per design.",
    )


def _run_meta(out_dir: Path) -> dict:
    """The run's sidecar, or ``{}`` when absent or unreadable."""
    path = out_dir / RUN_META_FILENAME
    try:
        return json.loads(path.read_text()) if path.is_file() else {}
    except (OSError, json.JSONDecodeError):
        return {}


def resolve_stage(ctx: CollectCtx[BindCraft2CollectArgs], meta: dict) -> str:
    """The stage to collect: the flag, or what the run recorded in the sidecar.

    A run older than the sidecar key (or one whose settings set ``trajectory_only``
    through the verbatim ``--set``) reads as a full campaign, which is the safe
    default -- it looks in 3_Ranked and reports nothing rather than inventing rows.
    """
    if ctx.args.stage != STAGE_AUTO:
        return ctx.args.stage
    return "trajectories" if meta.get("trajectory_only") else "ranked"


def resolve_campaigns_root(ctx: CollectCtx[BindCraft2CollectArgs], meta: dict) -> Path:
    """Where the campaigns actually live.

    Normally this out_dir, but a ``--reuse-campaigns`` run reserves a table over an
    EARLIER run's campaigns so one campaign can be collected into two tables (its
    backbones into one, its accepted designs into another). The run records the path;
    collect follows it rather than assuming the two coincide.
    """
    recorded = meta.get("campaigns_root")
    return Path(recorded) if recorded else ctx.out_dir / CAMPAIGNS_DIRNAME


def _read_stage_table(campaign_dir: Path, stage: str) -> tuple[pd.DataFrame, Path]:
    """The stage's CSV as a frame, plus the stage dir its structures hang off.

    An absent stage folder or CSV yields an empty frame: a campaign that accepted
    nothing (or was killed before its first acceptance) is a normal outcome, not an
    error -- `2_Refolded` then still says why.
    """
    folder, csv_name, _ = STAGES[stage]
    stage_dir = campaign_dir / folder
    csv_path = stage_dir / csv_name
    if not csv_path.is_file():
        return pd.DataFrame(), stage_dir
    try:
        return pd.read_csv(csv_path), stage_dir
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        print(f"  unreadable {csv_path}: {e}")
        return pd.DataFrame(), stage_dir


def _row_data(row: pd.Series, columns: list[str]) -> dict[str, Any]:
    """The tool-specific columns for one design, as bare names (the driver prefixes)."""
    data: dict[str, Any] = {}
    if SEQUENCE_COL in row.index:
        # Named `sequence` like every other sequence-producing tool here, so a
        # downstream predictor takes `-i bindcraft2_sequence` and nothing else.
        data["sequence"] = row[SEQUENCE_COL]
    for col in columns:
        if col in row.index:
            data[col] = row[col]
    return data


def collect_bindcraft2(ctx: CollectCtx[BindCraft2CollectArgs]) -> CollectEach:
    """Per-campaign BindCraft2 collector. A create tool: the framework iterates the
    ready parents (table rows, or the root design groups the run recorded) and this
    mints one child row per design the campaign produced, carrying ``parent`` for
    lineage. A campaign with no results yields no rows. The framework stamps
    status/path/parent_name from each Collected."""
    meta = _run_meta(ctx.out_dir)
    stage = resolve_stage(ctx, meta)
    find_structure = STAGES[stage].find
    campaigns_root = resolve_campaigns_root(ctx, meta)
    run_dir = ctx.args.run_dir
    extra_metrics = [m.strip() for m in ctx.args.metrics.split(",") if m.strip()]
    how = "from --stage" if ctx.args.stage != STAGE_AUTO else "per the run's sidecar"
    print(f"Collecting '{stage}' designs ({how}) from campaigns in {campaigns_root}")

    def one(d: DesignCtx) -> Iterable[Collected]:
        campaign_dir = campaigns_root / d.name
        if not campaign_dir.is_dir():
            print(f"{d.name}: no campaign dir, skipping")
            return

        df, stage_dir = _read_stage_table(campaign_dir, stage)
        if df.empty or DESIGN_COL not in df.columns:
            print(f"{d.name}: no '{stage}' designs")
            return

        # archive_trajectories packs each design folder into <design>.zip, and the
        # structures then exist only inside it. Say that, rather than reporting a
        # campaign's worth of designs as having no structure.
        if stage == "trajectories" and any(stage_dir.glob("*.zip")):
            print(
                f"{d.name}: {stage_dir} holds archived trajectories — run "
                f"`bindcraft unarchive {campaign_dir}` before collecting"
            )

        columns = (
            [c for c in df.columns if c not in (DESIGN_COL, SEQUENCE_COL)]
            if ctx.args.all_metrics
            else [c for c in (*CORE_METRICS, *extra_metrics) if c in df.columns]
        )

        # Sorted by identity, not by rank: !_Ranked.csv is re-sorted after every
        # acceptance, so collecting in file order would be collecting in an order
        # that changes under a resumed campaign.
        n = 0
        ordered = df.sort_values(DESIGN_COL).set_index(DESIGN_COL)
        for design, row in ordered.iterrows():
            design = str(design)
            cif = find_structure(stage_dir, design)
            if cif is None:
                print(f"{d.name}: {design} has no structure under {stage_dir}, skipping")
                continue
            data = _row_data(row, columns)
            data["cif_path"] = str(cif)
            if stage == "trajectories":
                # `terminated` is blank when the attempt ran to completion, so its
                # NA means success -- the opposite of NA everywhere else in a
                # prosapia table. Carry the polarity explicitly so a --filter reads
                # `bindcraft2_completed == True` instead of testing for a blank.
                stopped_at = row.get(TERMINATED_COL)
                data["completed"] = bool(
                    pd.isna(stopped_at) or not str(stopped_at).strip()
                )
            # A design whose mmCIF won't parse still has its metrics, and one bad
            # file among hundreds shouldn't abort the collect: keep the row, point
            # it at the mmCIF, and say so in its status.
            try:
                path, status = ensure_pdb(cif, run_dir), "OK"
            except (ValueError, RuntimeError, OSError) as e:
                print(f"{d.name}: {design} mmCIF unreadable ({e})")
                path, status = cif, f"error: unreadable mmCIF ({type(e).__name__})"
            yield Collected(
                name=design if design.startswith(d.name) else f"{d.name}_{design}",
                parent=d.name,
                path=path,
                status=status,
                data=data,
            )
            n += 1
        print(f"{d.name}: OK ({n} {stage} design(s))")

    return one
