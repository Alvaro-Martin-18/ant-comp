#!/usr/bin/env python3
"""
Submit BindCraft2 binder-design campaigns, one array task per target.

BindCraft2 is campaign-driven, not design-driven: one ``bindcraft design
<settings.json>`` process takes a target, hallucinates binder backbones with AF2,
redesigns their sequences with ProteinMPNN, refolds and filters the candidates,
and keeps going until ``number_of_final_designs`` have been accepted or
``max_trajectories`` attempts are spent. So the unit of work here is **one campaign
per target**, not one task per design -- the designs only exist once the campaign
has run, which is why this is a ``create`` tool and why the row count is known only
at collect time.

This tool writes each campaign's settings JSON itself (from the flags below, merged
with ``--extra-settings``) and points ``project_folder`` at
``<out_dir>/campaigns/<name>/``. Everything BindCraft2 exposes that has no dedicated
flag is reachable through ``--extra-settings`` (a file, merged into every campaign)
or ``--set KEY=VALUE`` (verbatim, per run).

Hotspots, coldspots and binder lengths are authored in BindCraft2's own syntax, with
``{expr}`` placeholders resolved per-design against the table lineage (integers, bare
column names, and + - * // arithmetic; see resolve_expr):

    --hotspots 'A54,A56,A66-70'            # literal, BindCraft2 syntax untouched
    --hotspots 'A{epitope_start}-{epitope_end}'   # resolved per row

Usage:
    # child run: one campaign per target structure already in a table
    sapia run bindcraft2 outputs/RUN --table table0 -i pdb_path \\
        --hotspots 'A54,A56,A66,A115' --binder-lengths 60-100 --num-designs 10

    # root run: one campaign against a target not in any table yet
    sapia run bindcraft2 outputs/RUN \\
        --target-pdb targets/PDL1.pdb --chains A \\
        --hotspots 'A54,A56' --modality VHH --property humanize

    # root run against a target BindCraft2 ships
    sapia run bindcraft2 outputs/RUN --shipped-target hPDL1 --num-designs 10

    # backbones only: stop before BindCraft2's own ProteinMPNN, to redesign the
    # sequences with this workspace's tools instead
    sapia run bindcraft2 outputs/RUN --target-pdb targets/PDL1.pdb \\
        --trajectory-only --max-trajectories 40 --hotspots 'A54,A56'
"""

import json
import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, cast

import yaml

from prosapia.core import CommonArgs, ManifestCtx, build_tool_leaf
from prosapia.core.data_manager import LookupFn
from prosapia.core.executors import volume_path
from prosapia.utils import resolve_template

TOOL_NAME = "bindcraft2"

# Layout this tool imposes on its out_dir. collect_bindcraft2.py reads the same
# convention -- keep the two in step.
CAMPAIGNS_DIRNAME = "campaigns"
SETTINGS_DIRNAME = "settings"

# Design-property presets BindCraft2 ships under settings/property/, each turned on
# by its own ``--<name>`` flag on the bindcraft CLI. Kept as a list so --property
# validates against it instead of forwarding a typo that bindcraft would swallow as
# a stray path argument.
PROPERTY_PRESETS = (
    "bigbang",
    "disulfide_staple",
    "forced_targeting",
    "humanize",
    "initial_guess",
    "mixed_topology",
    "protease_stable",
    "termini_accessible",
    "termini_together",
)

# A {expr} placeholder island, resolved per-design up the table lineage. Meaningless
# in a root run (no table), so we reject it there.
_HAS_PLACEHOLDER = re.compile(r"\{[^}]*\}")


class SettingsConfigError(ValueError):
    """A run-wide settings misconfiguration that applies to every campaign (e.g. an
    --extra-settings key colliding with a dedicated flag). Unlike a per-row error it
    is not swallowed by the warn-and-skip loop -- it fails the whole submit up front."""


class BindCraft2Args(CommonArgs):
    trajectory_only: bool
    reuse_campaigns: str | None
    target_pdb: Path | None
    shipped_target: str | None
    chains: str | None
    hotspots: str | None
    coldspots: str | None
    binder_lengths: str | None
    num_designs: int | None
    max_trajectories: int | None
    modality: str | None
    property: list[str]
    core: str | None
    campaign_seed: int | None
    no_resume: bool
    extra_settings: str | None
    set: list[str]


def add_run_bindcraft2_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--trajectory-only",
        action="store_true",
        help="Stop each campaign after the AF2 hallucination stage: produce "
        "BACKBONES and nothing else, with no ProteinMPNN redesign, no refold, no "
        "filtering and no accepted design. Use it to hand the backbones to this "
        "workspace's own sequence designers (`atomium`, `proteinmpnn`) instead of "
        "BindCraft2's. Implies `save_design_trajectory` so the structures are "
        "written, and makes --max-trajectories the only budget (BindCraft2 defaults "
        "it to 100). `sapia collect` picks the matching stage up automatically.",
    )
    parser.add_argument(
        "--reuse-campaigns",
        type=str,
        default=None,
        metavar="TABLE[:LABEL]",
        help="Submit nothing; reserve a table over the campaigns an EARLIER "
        "bindcraft2 run already wrote, named by the table it collected into (and its "
        "-l label, if any). Use it to collect a second stage of one campaign into a "
        "second table -- e.g. the backbones into one and BindCraft2's own accepted "
        "designs into another, from the same GPU hours. Pass the SAME -t the "
        "original run used, so the design groups still line up, and pick the stage "
        "with `sapia collect --stage`. Incompatible with every target and campaign "
        "flag, since no campaign is run.",
    )
    parser.add_argument(
        "--target-pdb",
        type=Path,
        default=None,
        help="Single target structure to design binders against in a ROOT run (no "
        "--table): a PDB/mmCIF/FASTA not in any table yet. Only valid without "
        "--table (with a table, targets come from --input-column). The design group "
        "is named `<stem>_bc2`.",
    )
    parser.add_argument(
        "--shipped-target",
        type=str,
        default=None,
        help="Name of a target BindCraft2 ships (hPDL1, hPD1, mPDL1, hIL2R, hIL7RA, "
        "dynorphin_a; `bindcraft design --list-targets` is authoritative). ROOT runs "
        "only, and mutually exclusive with --target-pdb. The design group is named "
        "`<name>_bc2`.",
    )
    parser.add_argument(
        "--chains",
        type=str,
        default=None,
        help="Target chains to design against (per-target `chains`, e.g. 'A' or "
        "'A,B'). Omitted by default (BindCraft2 uses every chain in the file).",
    )
    parser.add_argument(
        "--hotspots",
        type=str,
        default=None,
        help="Target residues the binder should contact (per-target `hotspots`), in "
        "BindCraft2's own syntax: comma-separated residues and ranges, chain-prefixed "
        "(e.g. 'A54,A56,A66-70'). May embed {expr} placeholders resolved per-design "
        "up the lineage. Omitted by default (BindCraft2 picks the epitope itself).",
    )
    parser.add_argument(
        "--coldspots",
        type=str,
        default=None,
        help="Target regions to avoid contacting (per-target `coldspots`), same "
        "syntax as --hotspots. Omitted by default.",
    )
    parser.add_argument(
        "--binder-lengths",
        type=str,
        default=None,
        help="Binder size (`binder_lengths`): 'N' for one length, 'min-max' for a "
        "range drawn from per trajectory (e.g. '80' or '60-100'). May embed {expr} "
        "placeholders. Omitted by default (BindCraft2's own default, or the "
        "modality preset's).",
    )
    parser.add_argument(
        "--num-designs",
        type=int,
        default=None,
        help="Accepted designs to stop the campaign at (`number_of_final_designs`). "
        "This is the number of CHILD ROWS a campaign aims to produce, not a batch "
        "size -- BindCraft2 keeps spending trajectories until it has them. Omitted "
        "by default.",
    )
    parser.add_argument(
        "--max-trajectories",
        type=int,
        default=None,
        help="Design attempts to spend before giving up (`max_trajectories`). The "
        "real cost knob: a campaign runs until --num-designs are accepted OR this "
        "many attempts are spent. Omitted by default.",
    )
    parser.add_argument(
        "--modality",
        type=str,
        default=None,
        help="Binder format (-> `bindcraft design --modality`): binder, VHH, "
        "peptide, cyclic_peptide, ARP, scFv, Fab, large_binder, homo_oligomer, "
        "multidomain, induced_fit, fold_switch. Comma-separated to combine. "
        "Defaults to BindCraft2's own default (binder).",
    )
    parser.add_argument(
        "--property",
        action="append",
        default=[],
        choices=PROPERTY_PRESETS,
        metavar="NAME",
        help="Design-property preset to switch on (-> `bindcraft design --<name>`). "
        "Repeatable. One of: " + ", ".join(PROPERTY_PRESETS) + ".",
    )
    parser.add_argument(
        "--core",
        type=str,
        default=None,
        help="Core profile applied under every preset (-> `bindcraft design --core`), "
        "e.g. 'benchmark' for a reproducible run. Omitted by default.",
    )
    parser.add_argument(
        "--campaign-seed",
        type=int,
        default=None,
        help="Seed every trajectory is drawn from (`campaign_seed`). Set it with "
        "--core benchmark for a reproducible campaign. Omitted by default.",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Start each campaign from scratch instead of carrying on into a "
        "project_folder already written. By default this tool sets `resume: true`, "
        "so re-running a submit continues the campaigns it already started (the "
        "framework's own resume filter can't help here -- a create tool's status "
        "column lives in the child table, not the one it reads).",
    )
    parser.add_argument(
        "--extra-settings",
        type=str,
        default=None,
        help="Path to a YAML or JSON file: a mapping of (extra) BindCraft2 campaign "
        "settings merged into EVERY campaign's settings file (e.g. objective, "
        "aa_bias, min_iptm_final, save_design_trajectory). String values may embed "
        "{expr} placeholders resolved per-design up the lineage. "
        "`bindcraft design --list-settings` names every setting it accepts.",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Extra `bindcraft design --set` override, appended verbatim and applied "
        "over the generated settings file. Repeatable. Escape hatch for settings "
        "without a dedicated flag. Must contain no spaces (the task script "
        "word-splits these tokens); use --extra-settings for anything richer.",
    )


def load_extra_settings(extra_settings: str | None) -> dict[str, Any]:
    """Parse the --extra-settings YAML/JSON file into a mapping of campaign settings.

    Returns ``{}`` when unset or empty. YAML is a JSON superset, so ``yaml.safe_load``
    parses both. Raises FileNotFoundError for a missing path and ValueError if the
    top level isn't a mapping.
    """
    if not extra_settings:
        return {}
    path = Path(extra_settings)
    if not path.is_file():
        raise FileNotFoundError(f"--extra-settings file not found: {extra_settings}")
    data = yaml.safe_load(path.read_text())
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(
            f"--extra-settings must be a mapping of campaign settings, got "
            f"{type(data).__name__}"
        )
    return data


def _tree_has_placeholder(obj: Any) -> bool:
    """True if any string key/value anywhere in ``obj`` embeds a ``{expr}`` island.

    Only string leaves are inspected -- structural dict/list braces don't count.
    """
    if isinstance(obj, dict):
        return any(
            _tree_has_placeholder(k) or _tree_has_placeholder(v) for k, v in obj.items()
        )
    if isinstance(obj, list):
        return any(_tree_has_placeholder(v) for v in obj)
    if isinstance(obj, str):
        return bool(_HAS_PLACEHOLDER.search(obj))
    return False


def _resolve_tree(obj: Any, lookup: LookupFn, name: str) -> Any:
    """Recursively resolve ``{expr}`` placeholders in a parsed settings structure.

    Strings (and dict keys) pass through ``resolve_template``; dicts and lists are
    walked; other scalars (int/float/bool/None) are returned unchanged. So native
    YAML/JSON types survive except where a string embeds ``{expr}``.
    """
    if isinstance(obj, dict):
        return {
            resolve_template(str(k), lookup, name): _resolve_tree(v, lookup, name)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [_resolve_tree(v, lookup, name) for v in obj]
    if isinstance(obj, str):
        return resolve_template(obj, lookup, name)
    return obj


def _coerce_binder_lengths(spec: str) -> list[int]:
    """'80' -> [80]; '60-100' / '60,100' -> [60, 100] (BindCraft2's [min, max])."""
    parts = [p.strip() for p in re.split(r"[-,]", spec) if p.strip()]
    try:
        return [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(
            f"--binder-lengths {spec!r} is not a length or a 'min-max' range"
        ) from e


def resolve_reuse_root(spec: str, run_dir: Path) -> Path:
    """``TABLE[:LABEL]`` -> the campaigns dir of an earlier bindcraft2 run.

    A run's outputs live at ``run_dir/<the table it reserved>/<leaf>/campaigns``, so
    naming that table (and the ``-l`` label it used, if any) is enough to find them.
    """
    table, _, label = spec.partition(":")
    if not table:
        raise ValueError(
            f"--reuse-campaigns {spec!r} names no table; expected TABLE[:LABEL], "
            f"e.g. 'table1' or 'table1:traj'."
        )
    root = run_dir / table / build_tool_leaf(TOOL_NAME, label) / CAMPAIGNS_DIRNAME
    if not root.is_dir():
        raise FileNotFoundError(
            f"--reuse-campaigns {spec!r} points at {root}, which does not exist. "
            f"Name the table the earlier bindcraft2 run collected into, plus its "
            f"-l label after a colon if it used one."
        )
    return root


def _reject_campaign_flags(args: BindCraft2Args) -> None:
    """A reuse run submits nothing, so anything describing a campaign is a mistake."""
    named = [
        flag
        for flag, value in (
            ("--target-pdb", args.target_pdb),
            ("--shipped-target", args.shipped_target),
            ("--trajectory-only", args.trajectory_only),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
            ("--chains", args.chains),
            ("--binder-lengths", args.binder_lengths),
            ("--num-designs", args.num_designs),
            ("--max-trajectories", args.max_trajectories),
            ("--modality", args.modality),
            ("--core", args.core),
            ("--campaign-seed", args.campaign_seed),
            ("--extra-settings", args.extra_settings),
            ("--property", args.property),
            ("--set", args.set),
        )
        if value
    ]
    if named:
        raise ValueError(
            f"--reuse-campaigns runs no campaign, so {', '.join(named)} would have "
            f"no effect. Drop them; the campaign was already run with its own "
            f"settings, and only `sapia collect --stage` still applies."
        )


def trajectory_only_run(args: BindCraft2Args, extra_fields: dict[str, Any]) -> bool:
    """Whether this run stops at backbones.

    Reads the flag *and* --extra-settings, because a user may set ``trajectory_only``
    in the settings file instead. The verbatim ``--set`` tokens are not inspected --
    setting it that way leaves the sidecar saying otherwise, and ``sapia collect``
    will default to the wrong stage.
    """
    return bool(args.trajectory_only or extra_fields.get("trajectory_only"))


def _build_settings(
    name: str,
    target_path: Path | None,
    args: BindCraft2Args,
    lookup: LookupFn,
    campaign_dir: Path,
    extra_fields: dict[str, Any],
) -> dict[str, Any]:
    """Build one campaign's BindCraft2 settings mapping.

    ``target_path`` is the (absolute) target structure, or ``None`` when the target
    comes from ``--shipped-target`` or from --extra-settings. Raises ValueError on a
    per-row problem (an unresolvable {expr}) and the run-wide ``SettingsConfigError``
    on an extra-settings collision.
    """
    # Settings this tool merely defaults, so --extra-settings can still turn them
    # off. Unlike tool_fields these are not collision-checked.
    defaults: dict[str, Any] = {}

    # Bookkeeping this tool owns outright: the campaign's identity and where it
    # writes. Not negotiable via --extra-settings, hence checked for collisions.
    tool_fields: dict[str, Any] = {
        "campaign_name": name,
        "project_folder": str(campaign_dir),
        "resume": not args.no_resume,
    }

    # Keep the fold each trajectory ended on. Under --trajectory-only it is the
    # entire output; under a full campaign it is what lets the same GPU hours be
    # collected a second time as backbones (see --reuse-campaigns). One CIF per
    # trajectory is cheap, so it is a default rather than a second flag --
    # --extra-settings can still turn it off.
    defaults["save_design_trajectory"] = True
    if args.trajectory_only:
        tool_fields["trajectory_only"] = True

    if target_path is not None:
        target: dict[str, Any] = {"name": name, "target_path": str(target_path)}
        if args.chains is not None:
            target["chains"] = resolve_template(args.chains, lookup, name)
        if args.hotspots is not None:
            target["hotspots"] = resolve_template(args.hotspots, lookup, name)
        if args.coldspots is not None:
            target["coldspots"] = resolve_template(args.coldspots, lookup, name)
        tool_fields["targets"] = [target]
    elif args.shipped_target is not None:
        # A shipped target is named rather than described: its own preset carries the
        # path and chains, so per-target overrides would have nowhere to land.
        tool_fields["target"] = args.shipped_target
        for flag, value in (
            ("--chains", args.chains),
            ("--hotspots", args.hotspots),
            ("--coldspots", args.coldspots),
        ):
            if value is not None:
                raise SettingsConfigError(
                    f"{flag} describes a target file and cannot be combined with "
                    f"--shipped-target (its preset already names the epitope). Pass "
                    f"the structure with --target-pdb, or override the preset with "
                    f"--extra-settings."
                )

    if args.binder_lengths is not None:
        tool_fields["binder_lengths"] = _coerce_binder_lengths(
            resolve_template(args.binder_lengths, lookup, name)
        )
    if args.num_designs is not None:
        tool_fields["number_of_final_designs"] = args.num_designs
    if args.max_trajectories is not None:
        tool_fields["max_trajectories"] = args.max_trajectories
    if args.campaign_seed is not None:
        tool_fields["campaign_seed"] = args.campaign_seed

    resolved_extra = _resolve_tree(extra_fields, lookup, name)
    collisions = sorted(set(tool_fields) & set(resolved_extra))
    if collisions:
        raise SettingsConfigError(
            f"setting(s) {collisions} set by both a dedicated flag and "
            "--extra-settings; remove them from one source"
        )

    return {**defaults, **resolved_extra, **tool_fields}


def _cli_flags(args: BindCraft2Args) -> str:
    """The run-wide ``bindcraft design`` flags (identical for every campaign).

    Presets and --set live on the command line rather than in the settings file
    because that is the interface BindCraft2 documents for them: --modality/--core
    name preset files to layer under the campaign, and --set is applied over it.
    """
    flags: list[str] = []
    if args.core is not None:
        flags += ["--core", args.core]
    if args.modality is not None:
        flags += ["--modality", args.modality]
    for prop in args.property:
        flags.append("--" + prop.replace("_", "-"))
    for assignment in args.set:
        flags += ["--set", assignment]
    return " ".join(flags)


def _target_groups(ctx: ManifestCtx[BindCraft2Args]) -> list[tuple[str, Path | None]]:
    """The (name, target_path) groups this run designs against: one per ready table
    row for a child run, a single group for a root run."""
    args = ctx.args

    if args.table is not None:
        if args.target_pdb is not None or args.shipped_target is not None:
            raise ValueError(
                "--target-pdb/--shipped-target are only valid for a root run (no "
                "--table); with --table, targets come from the table's "
                "--input-column. Drop one of them."
            )
        groups: list[tuple[str, Path | None]] = []
        for name in ctx.ready.index:
            name = cast(str, name)
            target_path = volume_path(str(ctx.ready.at[name, args.input_column]))
            if not target_path.exists():
                print(f"{name}: MISSING {target_path} (skipping)")
                continue
            groups.append((name, target_path))
        return groups

    # Root run: {expr} placeholders resolve up a table lineage this run doesn't have.
    if any(
        s and _HAS_PLACEHOLDER.search(s)
        for s in (args.hotspots, args.coldspots, args.chains, args.binder_lengths)
    ):
        raise ValueError(
            "--hotspots/--coldspots/--chains/--binder-lengths contain a {expr} "
            "placeholder, but this is a root run (no --table) with no table lineage "
            "to resolve it against. Use literal values, or run with --table."
        )

    if args.target_pdb is not None and args.shipped_target is not None:
        raise ValueError(
            "--target-pdb and --shipped-target both name a target; pass exactly one."
        )
    if args.target_pdb is not None:
        target_path = volume_path(args.target_pdb)
        if not target_path.exists():
            raise FileNotFoundError(f"--target-pdb {target_path} does not exist.")
        return [(f"{target_path.stem}_bc2", target_path)]
    if args.shipped_target is not None:
        return [(f"{args.shipped_target}_bc2", None)]

    raise ValueError(
        "A root run (no --table) needs a target: pass --target-pdb <file> or "
        "--shipped-target <name>. With --table, targets come from --input-column."
    )


def _reuse_existing_campaigns(ctx: ManifestCtx[BindCraft2Args]) -> list[tuple[str, ...]]:
    """Reserve this run's table over an earlier run's campaigns, and submit nothing.

    The table and its out_dir are reserved by the driver before the manifest is
    built, so returning no rows still leaves a table for `sapia collect` to fill --
    it just fills it from the recorded campaigns rather than from this out_dir.
    """
    _reject_campaign_flags(ctx.args)
    root = resolve_reuse_root(ctx.args.reuse_campaigns or "", ctx.args.run_dir)
    ctx.write_meta(campaigns_root=str(volume_path(root)), trajectory_only=False)

    if ctx.args.table is None:
        # Root reuse: no parent table to iterate, so the campaign dirs on disk ARE
        # the design groups -- the same record a root campaign run would write.
        groups = sorted(p.name for p in root.iterdir() if p.is_dir())
        ctx.write_meta(root_designs=groups)
    else:
        groups = sorted(map(str, ctx.ready.index))

    print(f"Reusing {len(groups)} campaign(s) under {root}")
    print("Nothing to submit. Collect the stage you want, e.g.:")
    print("  sapia collect bindcraft2 <run_dir> -t <this table> --stage ranked")
    return []


def build_bindcraft2_manifest(
    ctx: ManifestCtx[BindCraft2Args],
) -> list[tuple[str, ...]]:
    if ctx.args.reuse_campaigns:
        return _reuse_existing_campaigns(ctx)

    extra_fields = load_extra_settings(ctx.args.extra_settings)

    if ctx.args.table is None and _tree_has_placeholder(extra_fields):
        raise ValueError(
            "--extra-settings contains a {expr} placeholder, but this is a root run "
            "(no --table) with no table lineage to resolve it against."
        )

    trajectory_only = trajectory_only_run(ctx.args, extra_fields)
    if trajectory_only and ctx.args.num_designs is not None:
        print(
            "--num-designs is meaningless with --trajectory-only: no design is ever "
            "accepted, so the budget is --max-trajectories alone (BindCraft2 "
            "defaults it to 100)."
        )
    # The stage `sapia collect` should read is decided here, not guessed there:
    # a trajectory-only run fills 1_Trajectories and leaves 3_Ranked empty, and a
    # collect that looked in the wrong place would report zero rows rather than a
    # mismatch. Same contract as input_column -- the run records, collect reads.
    ctx.write_meta(
        trajectory_only=trajectory_only,
        campaigns_root=str(volume_path(ctx.out_dir) / CAMPAIGNS_DIRNAME),
    )

    groups = _target_groups(ctx)
    if not groups:
        return []

    settings_dir = ctx.out_dir / SETTINGS_DIRNAME
    settings_dir.mkdir(parents=True, exist_ok=True)
    campaigns_root = volume_path(ctx.out_dir) / CAMPAIGNS_DIRNAME

    flags = _cli_flags(ctx.args)
    manifest_rows: list[tuple[str, ...]] = []
    submitted: list[str] = []
    for name, target_path in groups:
        try:
            settings = _build_settings(
                name,
                target_path,
                ctx.args,
                ctx.lookup,
                campaigns_root / name,
                extra_fields,
            )
        except SettingsConfigError:
            # A run-wide misconfiguration hits every row identically: fail fast
            # instead of silently skipping the entire table.
            raise
        except ValueError as e:
            # One bad row (e.g. an unresolvable {expr}) shouldn't sink the whole
            # array: warn and skip it.
            print(f"{name}: {e} (skipping)")
            continue
        settings_json = settings_dir / f"{name}.json"
        settings_json.write_text(json.dumps(settings, indent=2))
        manifest_rows.append((name, str(volume_path(settings_json)), flags))
        submitted.append(name)

    if ctx.args.table is None:
        # A root run has no parent table for collect to iterate; record the group
        # names so `sapia collect` can find their campaigns and rebuild their rows.
        ctx.write_meta(root_designs=submitted)

    return manifest_rows
