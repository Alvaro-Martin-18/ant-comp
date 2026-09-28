#!/usr/bin/env python3
"""
Submit an array (SLURM or Modal) scoring how each design fits a reference assembly.

Each task runs ringfit_worker.py on one design: it superimposes the design's
target chains (two adjacent protomers) onto the matching chains of a reference
full assembly, then measures (a) whether the binder straddles BOTH protomers
(buried surface area per target chain + bridge_ratio) and (b) whether it would
clash with the rest of the assembly or its lipid belt -- surface that is
artificially exposed when designing against an isolated protomer pair. Use
``sapia collect ringfit`` to merge the metrics back into the table.

The structure column is whatever ``--input-column`` says (default
``rfdiffusion3_path``); point it at ``boltz_path`` to score predictions instead
of backbones. Nothing about the tool assumes a particular column.

Usage:
    sapia run ringfit outputs/20260928_100140_7ojg_binder \
        --table table1 \
        --input-column boltz_path \
        --ref-structure inputs/7ojg_assembly.cif \
        --target-chains A,B --ref-target-chains A,B \
        --hotspots A35,B30 --hotspots-in-ref-numbering
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import cast

from prosapia.core import CommonArgs, ManifestCtx
from prosapia.core.executors import volume_path


class RingfitArgs(CommonArgs):
    ref_structure: str
    target_chains: str
    binder_chains: str
    ref_target_chains: str
    hotspots: str
    hotspots_in_ref_numbering: bool
    resnum_match: str
    clash_cutoff: float
    contact_cutoff: float
    lipid_resnames: str


def add_run_ringfit_args(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--ref-structure",
        type=str,
        required=True,
        help="Reference full assembly (PDB or mmCIF) the design is fitted into, "
        "e.g. the C11 ring. Must live on the runs volume under --executor modal.",
    )
    parser.add_argument(
        "--target-chains",
        type=str,
        default="A,B",
        help="The design's two target-protomer chains, comma-joined and ordered: "
        "the first is t1 and the second t2 in the bsa_t*/n_contact_res_t* columns. "
        "Default 'A,B'.",
    )
    parser.add_argument(
        "--binder-chains",
        type=str,
        default="auto",
        help="The design's binder chains, comma-joined. 'auto' (default) takes "
        "every protein chain that is not a target chain.",
    )
    parser.add_argument(
        "--ref-target-chains",
        type=str,
        default="A,B",
        help="The reference chains the design's --target-chains are superimposed "
        "onto, in the same order. Every other reference protein chain is what the "
        "binder is checked against. Default 'A,B'.",
    )
    parser.add_argument(
        "--hotspots",
        type=str,
        default="",
        help="Target residues the binder is meant to touch, e.g. 'A47,B42'. "
        "Given in the DESIGN's numbering and chain names unless "
        "--hotspots-in-ref-numbering is set. Empty (default) skips "
        "hotspot_recall/hotspot_hits (collected as NA).",
    )
    parser.add_argument(
        "--hotspots-in-ref-numbering",
        action="store_true",
        help="Read --hotspots in the REFERENCE's numbering and chain names "
        "(--ref-target-chains) and translate them to the design through the "
        "residue match, instead of taking them as design numbering. Set this "
        "whenever you copied hotspots off the reference structure: a generator "
        "such as rfdiffusion3 renumbers every chain from 1, so reference "
        "residue A35 is design residue A18 and a hotspot list in reference "
        "numbering silently scores the WRONG residues -- a meaningless "
        "hotspot_recall rather than an error. Check ringfit_resnum_offset "
        "(0 = the two numberings agree, so the flag makes no difference).",
    )
    parser.add_argument(
        "--resnum-match",
        type=str,
        choices=("auto", "ordinal", "resnum"),
        default="auto",
        help="How design residues are paired with reference residues for the "
        "superposition. 'ordinal' pairs the i-th residue of a design target "
        "chain with the i-th of its reference chain (right when the design was "
        "renumbered from 1); 'resnum' pairs equal residue numbers (right when "
        "the design kept the reference's numbering, and tolerant of gaps); "
        "'auto' (default) takes ordinal when the two chains hold the same "
        "number of CA atoms and resnum otherwise. Either way the match is "
        "checked by residue identity and reported as ringfit_seq_match_frac.",
    )
    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=2.5,
        help="Heavy-atom distance (A) below which a binder atom counts as clashing "
        "with the rest of the reference assembly. Default 2.5.",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=5.0,
        help="Heavy-atom distance (A) below which a target residue counts as "
        "contacted by the binder. Default 5.0.",
    )
    parser.add_argument(
        "--lipid-resnames",
        type=str,
        default="PLM,LPP,L8Z",
        help="Reference HETATM component ids marking the membrane/lipid belt, "
        "comma-joined. Empty disables the lipid check (collected as NA). "
        "Default 'PLM,LPP,L8Z'.",
    )


def build_ringfit_manifest(ctx: ManifestCtx[RingfitArgs]) -> list[tuple[str, ...]]:
    ctx.args.gpus_per_task = 0  # CPU-only tool

    ref = Path(ctx.args.ref_structure)
    if not ref.exists():
        raise FileNotFoundError(
            f"--ref-structure {ref} does not exist. Under --executor modal it must "
            f"be a path on the runs volume (the workstation's /runs)."
        )
    # Must resolve inside a task container, not just in the submitting shell.
    ref_path = str(volume_path(ref))

    ready = ctx.ready

    manifest_rows: list[tuple[str, ...]] = []
    for name in ready.index:
        name = cast(str, name)
        # No CIF->PDB staging: the worker reads either format with gemmi.
        design = Path(str(ready.at[name, ctx.args.input_column]))
        if not design.exists():
            print(f"{name}: MISSING {design} (skipping)")
            continue
        manifest_rows.append(
            (
                name,
                str(volume_path(design)),
                ref_path,
                ctx.args.target_chains,
                ctx.args.binder_chains,
                ctx.args.ref_target_chains,
                ctx.args.hotspots,
                ctx.args.lipid_resnames,
                # Kept last: always non-empty, so the manifest line never ends
                # on an empty field.
                str(ctx.args.clash_cutoff),
                str(ctx.args.contact_cutoff),
                "ref" if ctx.args.hotspots_in_ref_numbering else "design",
                ctx.args.resnum_match,
            )
        )

    return manifest_rows
