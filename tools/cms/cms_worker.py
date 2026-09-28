#!/usr/bin/env python3
"""
Contact molecular surface (CMS) and shape complementarity (SC) for a batch of complexes.

This is the per-array-task step of the cms tool. One task scores several designs,
because a design takes ~0.6 s on an L4 while a container cold start and the
first-call CUDA kernel compilation take ~1 minute: packing designs into a task
pays those costs once. Each design still gets its own calculator (as cms-cuda
recommends: GPU memory is pooled and reused between structures) and its own error
handling, so one bad structure never costs the rest of the batch.

The numbers come from cms-cuda (github.com/ullahsamee/cms-cuda), the CUDA port of
Brian Coventry's py_contact_ms / Longxing Cao's contact molecular surface, plus
Lawrence & Colman SC off the same surface. By convention CMS is reported on the
TARGET side (``cms_target``); the binder-side CMS (``cms_binder``) comes from the same
calculation at no extra cost.

The atom loader follows the rules of the cms-cuda tutorial (and its
examples/calculate_cms_sc.py): first model only, one altloc per atom, heavy atoms
only, waters and single-atom ions dropped, radii ONLY from ``get_radii_from_names``.
Any other choice changes the surface, so it is not configurable beyond
--exclude-resnames.

Error discipline (as in chainsel): **a requested chain that is absent is an error,
not a smaller selection.** A binder with its chain missing would otherwise score a
plausible, wrong CMS. An empty selection is an error too. A CMS of 0 is NOT an
error: two selections more than 8 A apart honestly score 0, a real (bad) result.

Writes, per design, ``<name>.tsv`` (one row: name, status, path, metrics...) and, on
success, ``<name>_per_residue.tsv`` (per-residue CMS on both sides -- the hotspots).
Errors are recorded as data (status 'error: ...').

Normally invoked per array task by the cms tool -- run ``sapia run cms`` rather than
calling this directly.

Usage (standalone):
    printf 'design_0\\tcomplex.pdb\\n' > task.tsv
    python tools/cms/cms_worker.py --task-file task.tsv \\
        --binder-chains B --target-chains A --out-dir out/ --device cpu
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

# Metric columns of the per-design TSV, in order (bare names: collect_cms.py hands
# them to the driver, which leaf-prefixes them to cms_<name>).
METRIC_COLUMNS = [
    "target",
    "binder",
    "sc",
    "sc_area",
    "sc_median_dist",
    "max_target",
    "max_binder",
    "frac_target",
    "frac_binder",
    "n_atoms_binder",
    "n_atoms_target",
    "n_radius0_binder",
    "n_radius0_target",
    "device",
    "seconds",
]
RESULT_COLUMNS = ["name", "status", "path", *METRIC_COLUMNS]
PER_RESIDUE_COLUMNS = ["side", "chain", "resnum", "resname", "cms"]

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}


def split_list(value: str) -> list[str]:
    """Comma-joined list -> stripped, non-empty tokens."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def load_model(path: Path):
    """First model of a PDB or mmCIF file (Biopython; auth chain IDs for mmCIF)."""
    from Bio.PDB import MMCIFParser, PDBParser

    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    is_cif = path.name.lower().endswith((".cif", ".mmcif"))
    parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(path))
    models = list(structure)
    if not models:
        raise ValueError(f"no model in {path.name}")
    return models[0]


def select_atoms(model, chains: list[str], exclude: set[str]):
    """Heavy atoms of the named chains -> (xyz, radii, info).

    info is one (chain, resnum+icode, resname) tuple per atom, for the per-residue
    table. Mirrors cms-cuda's ``load_atoms`` with its defaults (no water, no ions,
    no hydrogens), plus ``exclude`` residue names.
    """
    from cms_cuda import get_radii_from_names

    xyz, res_names, atom_names, info = [], [], [], []
    for chain in model:
        if chain.id not in chains:
            continue
        for residue in chain:  # Biopython keeps one altloc per atom
            resname = residue.get_resname().strip()
            hetflag = residue.id[0]
            if resname in exclude:
                continue
            if hetflag == "W" or resname in WATER_NAMES:
                continue
            if hetflag.startswith("H_") and len(residue) == 1:
                continue  # single-atom HETATM = ion
            for atom in residue:
                name = atom.get_name().strip()
                element = (atom.element or "").strip().upper()
                is_h = element in ("H", "D") or (
                    element in ("", "X")
                    and name.lstrip("0123456789").startswith(("H", "D"))
                )
                if is_h:
                    continue
                xyz.append(atom.coord)
                res_names.append(resname)
                atom_names.append(name)
                info.append(
                    (chain.id, f"{residue.id[1]}{residue.id[2].strip()}", resname)
                )

    xyz = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
    radii = get_radii_from_names(res_names, atom_names)
    return xyz, radii, info


def per_residue(side: str, info, per_atom) -> list[list]:
    """Sum a per-atom CMS array over residues, keeping file order."""
    sums: dict[tuple, float] = {}
    for key, value in zip(info, per_atom):
        sums[key] = sums.get(key, 0.0) + float(value)
    return [[side, *key, round(v, 4)] for key, v in sums.items()]


def make_calculator(device: str):
    if device == "cuda":
        from cms_cuda import MolecularSurfaceCalculatorGPU

        return MolecularSurfaceCalculatorGPU()
    from cms_cuda import MolecularSurfaceCalculator

    return MolecularSurfaceCalculator()


def score(
    src: Path,
    binder_chains: list[str],
    target_chains: list[str],
    exclude: set[str],
    do_sc: bool,
    do_max: bool,
    device: str,
) -> tuple[dict, list[list]]:
    """One complex -> (metrics, per-residue rows). Raises on anything unusable."""
    from cms_cuda import calculate_maximum_possible_contact_ms

    t0 = time.perf_counter()
    model = load_model(src)
    present = [chain.id for chain in model]
    for chain_id in binder_chains + target_chains:
        if chain_id not in present:
            raise ValueError(
                f"chain {chain_id} not in {src.name} (have: {','.join(present)})"
            )

    bx, br, binfo = select_atoms(model, binder_chains, exclude)
    tx, tr, tinfo = select_atoms(model, target_chains, exclude)
    if len(bx) == 0 or len(tx) == 0:
        raise ValueError(
            f"empty selection (binder {len(bx)} atoms, target {len(tx)} atoms)"
        )

    calc = make_calculator(device)
    calc.add_binder_and_target(bx, br, tx, tr)
    cms_target, per_atom_target = calc.CalcLoaded()
    cms_binder, per_atom_binder = calc.calc_contact_molecular_surface(
        target_side=False
    )
    metrics: dict = {
        "target": round(float(cms_target), 4),
        "binder": round(float(cms_binder), 4),
        "n_atoms_binder": len(bx),
        "n_atoms_target": len(tx),
        # Atoms with no entry in the radius table get radius 0 and are invisible.
        "n_radius0_binder": int((br == 0).sum()),
        "n_radius0_target": int((tr == 0).sum()),
        "device": device,
    }
    if do_sc:
        sc, sc_area, median_dist = calc.CalcLoadedSC()
        metrics.update(
            sc=round(float(sc), 4),
            sc_area=round(float(sc_area), 4),
            sc_median_dist=round(float(median_dist), 4),
        )
    if do_max:
        # A fresh calculator per side: max-CMS cannot share the CMS surface.
        max_t = float(calculate_maximum_possible_contact_ms(tx, tr, device=device)[0])
        max_b = float(calculate_maximum_possible_contact_ms(bx, br, device=device)[0])
        metrics.update(
            max_target=round(max_t, 4),
            max_binder=round(max_b, 4),
            frac_target=round(float(cms_target) / max_t, 4) if max_t else None,
            frac_binder=round(float(cms_binder) / max_b, 4) if max_b else None,
        )
    metrics["seconds"] = round(time.perf_counter() - t0, 3)

    rows = per_residue("target", tinfo, per_atom_target) + per_residue(
        "binder", binfo, per_atom_binder
    )
    return metrics, rows


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def write_result(out_dir: Path, name: str, status: str, path: str, metrics: dict):
    """One-row TSV. A metric that did not apply is written as NA (an empty cell),
    never as 0 -- a CMS of 0 is a real, alarming value (no interface)."""
    write_tsv(
        out_dir / f"{name}.tsv",
        RESULT_COLUMNS,
        [
            [name, status, path]
            + ["" if metrics.get(c) is None else str(metrics[c]) for c in METRIC_COLUMNS]
        ],
    )


class Args(argparse.Namespace):
    task_file: Path
    binder_chains: str
    target_chains: str
    exclude_resnames: str
    no_sc: bool
    max_cms: bool
    device: str
    out_dir: Path


def main() -> None:
    ap = argparse.ArgumentParser(description="CMS + SC for a batch of complexes.")
    ap.add_argument(
        "--task-file",
        type=Path,
        required=True,
        help="Tab-separated 'name<TAB>structure' rows, one per design.",
    )
    ap.add_argument("--binder-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--target-chains", required=True, help="Comma-joined chain IDs.")
    ap.add_argument("--exclude-resnames", default="", help="Comma-joined resnames.")
    ap.add_argument("--no-sc", action="store_true", help="Skip shape complementarity.")
    ap.add_argument("--max-cms", action="store_true", help="Also compute max CMS.")
    ap.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(namespace=Args())

    binder_chains = split_list(args.binder_chains)
    target_chains = split_list(args.target_chains)
    exclude = set(split_list(args.exclude_resnames))

    with open(args.task_file) as f:
        members = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    if args.device == "cuda":
        # Fail every design loudly rather than silently run 100x slower on the CPU.
        from cms_cuda import gpu_available

        if not gpu_available():
            for name, _src in members:
                write_result(
                    args.out_dir,
                    name,
                    "error: no usable CUDA GPU (CuPy missing, wrong CUDA build, or "
                    "no GPU visible); rerun with -g 0 for the CPU path",
                    "",
                    {},
                )
            print("ERROR: device cuda requested but no usable GPU", file=sys.stderr)
            return

    for name, src in members:
        per_res_path = args.out_dir / f"{name}_per_residue.tsv"
        try:
            metrics, rows = score(
                Path(src),
                binder_chains,
                target_chains,
                exclude,
                not args.no_sc,
                args.max_cms,
                args.device,
            )
            write_tsv(per_res_path, PER_RESIDUE_COLUMNS, rows)
            write_result(args.out_dir, name, "OK", str(per_res_path), metrics)
            print(
                f"{name}: cms_target={metrics['target']} cms_binder={metrics['binder']}"
                f" sc={metrics.get('sc')} ({metrics['seconds']} s, {args.device})"
            )
        except Exception as e:  # noqa: BLE001 - errors are recorded as data
            print(f"{name}: ERROR {e}", file=sys.stderr)
            per_res_path.unlink(missing_ok=True)
            write_result(args.out_dir, name, f"error: {e}", "", {})


if __name__ == "__main__":
    main()
