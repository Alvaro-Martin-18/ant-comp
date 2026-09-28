#!/usr/bin/env python3
"""
Score one design's fit into a reference oligomeric assembly (the "ring test").

This is the per-array-task step of the ringfit tool (and works standalone for a
single design). It answers two questions a binder designed against an isolated
pair of adjacent protomers cannot answer on its own:

1. *Does the binder really bridge BOTH protomers?* The design's target chains are
   two adjacent protomers; a binder that only grips one of them and brushes the
   other is not a seam binder. Measured as buried surface area (BSA) per target
   chain, and summarised by ``bridge_ratio = min(bsa_t1, bsa_t2) / max(...)``,
   which is 1.0 for a perfectly even straddle and ~0 for a one-protomer binder.
2. *Does the binder fit the real assembly?* Designing against an isolated dimer
   leaves the flanking faces of the two protomers artificially exposed, so a
   binder may happily occupy surface that is buried by the neighbouring
   protomers (or by the lipid belt) in the biological assembly. Measured by
   superimposing the design's target chains onto the corresponding chains of a
   reference full assembly and counting binder atoms that land on the *other*
   reference chains (``n_clash``, ``min_dist_ring``) or on its lipids
   (``lipid_clash``, ``min_dist_lipid``).

The superposition is a Kabsch fit on CA atoms of design target chain i against
reference target chain i. Residues are paired **ordinally** (i-th to i-th) when
the two chains have the same CA count, and by residue number otherwise
(``--resnum-match`` forces either). Ordinal is the default because a design's
target chains are routinely renumbered from 1 by the generator while the
reference keeps its crystallographic numbering (7OJG: 18-155) -- matching on
resnum would then align the wrong residues and yield a wrong-but-plausible fit.
The pairing is checked, not trusted: ``seq_match_frac`` is the fraction of
matched pairs whose one-letter residue identity agrees (1.0 when the pairing is
right), and ``resnum_offset`` the modal ``ref_resnum - design_resnum`` (0 when
the numbering agreed, 17 for the 7OJG case). A ``seq_match_frac`` below 0.9
warns on stderr and downgrades the status to ``warn: ...`` while still reporting
every metric. ``align_rmsd`` remains the trust metric for every ring/lipid
number below it.

BSA and contacts are computed on the DESIGN (they are properties of the design,
not of the reference), and are frame-invariant. SASA is a Shrake-Rupley
numerical estimate (golden-spiral sphere points, 1.4 A probe, Bondi radii) with
a cell-list neighbour search -- gemmi has no SASA and there is no scipy in these
images.

``t1``/``t2`` refer to the **order of --target-chains**, not to the literal chain
letters: with ``--target-chains A,B``, ``bsa_t1``/``n_contact_res_t1`` are chain A's.

Writes a one-row TSV (name, status, aligned_path, then the metrics) that
collect_ringfit.py merges back into the table, and the superposed design as
``<name>_onref.pdb`` so the fit can be inspected in the reference frame. Errors
are recorded as data (a status starting with 'error:') rather than only
crashing, so partial array runs still collect.

Normally invoked per array task by the ringfit tool -- run ``sapia run ringfit``
rather than calling this directly.

Usage (standalone, single design):
    python tools/ringfit/ringfit_worker.py \\
        --name design_0 --design design_0.pdb --ref 7ojg_assembly.cif \\
        --target-chains A,B --binder-chains auto --ref-target-chains A,B \\
        --hotspots A47,B42 --result-tsv design_0.tsv --out-pdb design_0_onref.pdb
"""

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

# Shrake-Rupley: enough points that the area of a lone atom is within ~0.1% of
# the analytical 4*pi*r^2, cheap enough for a few thousand atoms.
N_SPHERE_POINTS = 960
PROBE_RADIUS = 1.4
# Bondi van der Waals radii for the elements a protein/lipid model contains.
BONDI_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
DEFAULT_RADIUS = 1.70
MAX_RADIUS = max([*BONDI_RADII.values(), DEFAULT_RADIUS])
# Two atoms' solvent spheres can only overlap within this distance, so it bounds
# both the SASA neighbour search and the "can this atom's SASA change?" test.
SASA_INFLUENCE = 2 * (MAX_RADIUS + PROBE_RADIUS)
# Fixed cutoff for the lipid-belt check (the lipid flag is a yes/no, not a scan).
LIPID_CLASH_CUTOFF = 4.0
# How design residues are paired with reference residues for the superposition.
MATCH_MODES = ("auto", "ordinal", "resnum")
# Below this fraction of identity-agreeing matched pairs the pairing is suspect
# (a register shift, or simply the wrong chains): warn and downgrade the status.
SEQ_MATCH_MIN = 0.9

# Metric columns of the per-design TSV, in order (bare names: collect_ringfit.py
# hands them to the driver, which leaf-prefixes them to ringfit_<name>).
METRIC_COLUMNS = [
    "align_rmsd",
    "seq_match_frac",
    "resnum_offset",
    "n_clash",
    "n_clash_res",
    "min_dist_ring",
    "lipid_clash",
    "min_dist_lipid",
    "bsa_t1",
    "bsa_t2",
    "bsa_total",
    "bridge_ratio",
    "n_contact_res_t1",
    "n_contact_res_t2",
    "hotspot_recall",
    "hotspot_hits",
    "binder_len",
    "binder_chains",
]
RESULT_COLUMNS = ["name", "status", "aligned_path", *METRIC_COLUMNS]


# ---------------------------------------------------------------------------
# Atom bookkeeping
# ---------------------------------------------------------------------------


@dataclass
class AtomSet:
    """Heavy atoms of a residue selection, flattened for numpy.

    ``res_ids[k]`` indexes ``res_labels`` for atom ``k``, so a per-atom boolean
    mask collapses to a set of residues with one fancy-index.
    """

    coords: np.ndarray  # (N, 3) float
    radii: np.ndarray  # (N,) float, van der Waals (probe NOT added)
    res_ids: np.ndarray  # (N,) int, index into res_labels
    res_labels: list[str]  # e.g. "A47"

    def __len__(self) -> int:
        return int(self.coords.shape[0])

    @property
    def n_res(self) -> int:
        return len(self.res_labels)

    def residues_in(self, mask: np.ndarray) -> list[str]:
        """Labels of the residues with at least one atom selected by ``mask``."""
        ids = np.unique(self.res_ids[mask]) if len(self) else np.empty(0, dtype=int)
        return [self.res_labels[i] for i in ids]


def _empty_atom_set() -> AtomSet:
    return AtomSet(np.empty((0, 3)), np.empty(0), np.empty(0, dtype=int), [])


def _is_heavy(atom: gemmi.Atom) -> bool:
    """Non-hydrogen. Tested on the element name rather than
    ``Element.is_hydrogen``, which is a method in some gemmi versions and a
    property in others (0.7.x)."""
    return atom.element.name.upper() not in ("H", "D")


def _is_amino_acid(res: gemmi.Residue) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


def _atom_set(items: list[tuple[str, gemmi.Residue]]) -> AtomSet:
    """Flatten ``(label, residue)`` pairs into an AtomSet of their heavy atoms."""
    coords: list[list[float]] = []
    radii: list[float] = []
    res_ids: list[int] = []
    labels: list[str] = []
    for label, res in items:
        heavy = [a for a in res if _is_heavy(a)]
        if not heavy:
            continue
        idx = len(labels)
        labels.append(label)
        for atom in heavy:
            coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
            radii.append(BONDI_RADII.get(atom.element.name.upper(), DEFAULT_RADIUS))
            res_ids.append(idx)
    if not coords:
        return _empty_atom_set()
    return AtomSet(
        np.asarray(coords, dtype=float),
        np.asarray(radii, dtype=float),
        np.asarray(res_ids, dtype=int),
        labels,
    )


def _label(chain: gemmi.Chain, res: gemmi.Residue) -> str:
    """Residue key in the ``A47`` style used by --hotspots."""
    return f"{chain.name}{res.seqid.num}"


def protein_atoms(model: gemmi.Model, chain_names: list[str]) -> AtomSet:
    """Heavy atoms of the amino-acid residues of the named chains (in that order)."""
    items: list[tuple[str, gemmi.Residue]] = []
    for name in chain_names:
        for chain in model:
            if chain.name != name:
                continue
            items += [(_label(chain, r), r) for r in chain if _is_amino_acid(r)]
    return _atom_set(items)


def protein_chain_names(model: gemmi.Model) -> list[str]:
    """Names of every chain holding at least one amino-acid residue."""
    return [c.name for c in model if any(_is_amino_acid(r) for r in c)]


def het_atoms(model: gemmi.Model, resnames: set[str]) -> AtomSet:
    """Heavy atoms of every residue whose component id is in ``resnames``."""
    items = [
        (_label(chain, res), res)
        for chain in model
        for res in chain
        if res.name.upper() in resnames
    ]
    return _atom_set(items)


def one_letter(res: gemmi.Residue) -> str:
    """One-letter code of an amino acid, 'X' when gemmi does not know it."""
    info = gemmi.find_tabulated_residue(res.name)
    code = (info.one_letter_code if info else "").upper()
    return code if code.isalpha() else "X"


@dataclass
class CaResidues:
    """One chain's amino-acid CA atoms, in file order.

    ``resnums[i]``/``codes[i]``/``coords[i]`` describe the same residue, so a
    pairing can be expressed as index pairs regardless of how it was derived.
    """

    resnums: list[int]
    codes: list[str]
    coords: np.ndarray  # (N, 3) float

    def __len__(self) -> int:
        return len(self.resnums)


def ca_residues(model: gemmi.Model, chain_name: str) -> CaResidues:
    """CA atoms of one chain's amino acids, in the order the file lists them."""
    resnums: list[int] = []
    codes: list[str] = []
    coords: list[list[float]] = []
    for chain in model:
        if chain.name != chain_name:
            continue
        for res in chain:
            if not _is_amino_acid(res):
                continue
            for atom in res:
                if atom.name == "CA":
                    resnums.append(res.seqid.num)
                    codes.append(one_letter(res))
                    coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    break
    return CaResidues(
        resnums, codes, np.asarray(coords, dtype=float).reshape(len(coords), 3)
    )


def match_residues(
    design: CaResidues, ref: CaResidues, mode: str
) -> tuple[list[tuple[int, int]], str]:
    """Pair design residues with reference residues -> ``(index pairs, mode used)``.

    ``ordinal`` pairs the i-th residue of each chain (what a renumbered design
    needs); ``resnum`` intersects on residue number (tolerates gaps and a
    different order, but assumes the numbering agrees); ``auto`` takes ordinal
    when the two chains hold the same number of CA atoms and resnum otherwise.
    """
    if mode not in MATCH_MODES:
        raise ValueError(f"--resnum-match must be one of {MATCH_MODES}, got {mode!r}")
    if mode == "auto":
        mode = "ordinal" if len(design) == len(ref) and len(design) else "resnum"

    if mode == "ordinal":
        if len(design) != len(ref):
            raise ValueError(
                f"ordinal matching needs an equal number of CA atoms, got "
                f"{len(design)} (design) vs {len(ref)} (reference); use "
                f"--resnum-match auto or resnum"
            )
        return [(i, i) for i in range(len(design))], mode

    # Last occurrence wins on a duplicated residue number, as before.
    d_index = {num: i for i, num in enumerate(design.resnums)}
    r_index = {num: i for i, num in enumerate(ref.resnums)}
    shared = sorted(set(d_index) & set(r_index))
    return [(d_index[n], r_index[n]) for n in shared], mode


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def kabsch(p: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Optimal rigid superposition of ``p`` onto ``q`` (both (N, 3), paired).

    Returns ``(rot, trans, rmsd)`` with ``p @ rot.T + trans`` the moved ``p``.
    Reflections are excluded the usual way (flip the last singular vector when
    the naive rotation has determinant -1).
    """
    if p.shape != q.shape or p.shape[0] < 3:
        raise ValueError(f"need >=3 paired points, got {p.shape[0]}")
    p_mean, q_mean = p.mean(axis=0), q.mean(axis=0)
    pc, qc = p - p_mean, q - q_mean
    u, _, vt = np.linalg.svd(pc.T @ qc)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    trans = q_mean - rot @ p_mean
    rmsd = float(np.sqrt((((pc @ rot.T) - qc) ** 2).sum(axis=1).mean()))
    return rot, trans, rmsd


class _Grid:
    """Uniform cell list over a point set; ``near`` returns the 27-cell candidates.

    Cells are ``cell`` wide, so every point within ``cell`` of the query is in
    one of the 27 cells searched (the caller applies the exact distance test).
    """

    def __init__(self, coords: np.ndarray, cell: float) -> None:
        self.cell = cell
        self.coords = coords
        self.buckets: dict[tuple[int, int, int], list[int]] = {}
        if coords.size:
            keys = np.floor(coords / cell).astype(int)
            for i, key in enumerate(map(tuple, keys)):
                self.buckets.setdefault(key, []).append(i)  # type: ignore[arg-type]
        self.arrays = {k: np.asarray(v) for k, v in self.buckets.items()}

    def near(self, point: np.ndarray) -> np.ndarray:
        cx, cy, cz = (int(v) for v in np.floor(point / self.cell))
        found = [
            self.arrays[(cx + dx, cy + dy, cz + dz)]
            for dx in (-1, 0, 1)
            for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)
            if (cx + dx, cy + dy, cz + dz) in self.arrays
        ]
        return np.concatenate(found) if found else np.empty(0, dtype=int)


def _sphere_points(n: int) -> np.ndarray:
    """``n`` near-uniform points on the unit sphere (golden-spiral / Fibonacci)."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0**0.5) * i
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1
    )


SPHERE = _sphere_points(N_SPHERE_POINTS)


def sasa(
    coords: np.ndarray, radii: np.ndarray, eval_idx: np.ndarray | None = None
) -> np.ndarray:
    """Shrake-Rupley SASA (A^2) per atom, occluded by *all* of ``coords``.

    ``eval_idx`` restricts which atoms are evaluated (not which occlude), so a
    buried-area difference only has to touch the atoms that can actually change.
    """
    n_atoms = coords.shape[0]
    idx = np.arange(n_atoms) if eval_idx is None else np.asarray(eval_idx, dtype=int)
    if n_atoms == 0 or idx.size == 0:
        return np.zeros(idx.size)

    r = radii + PROBE_RADIUS
    grid = _Grid(coords, SASA_INFLUENCE)
    out = np.zeros(idx.size)
    for k, i in enumerate(idx):
        cand = grid.near(coords[i])
        cand = cand[cand != i]
        if cand.size:
            delta = coords[cand] - coords[i]
            dist = np.sqrt(np.einsum("ij,ij->i", delta, delta))
            keep = dist < r[i] + r[cand]
            cand, dist = cand[keep], dist[keep]
            cand = cand[np.argsort(dist)]  # nearest first: occludes the most
        points = coords[i] + r[i] * SPHERE
        for j in cand:
            delta = points - coords[j]
            open_mask = np.einsum("ij,ij->i", delta, delta) > r[j] ** 2
            if not open_mask.all():
                points = points[open_mask]
                if points.shape[0] == 0:
                    break
        out[k] = 4.0 * np.pi * r[i] ** 2 * points.shape[0] / N_SPHERE_POINTS
    return out


def pair_stats(
    a: np.ndarray, b: np.ndarray, cutoff: float, chunk: int = 256
) -> tuple[float | None, np.ndarray]:
    """``(min distance a<->b, mask of a-atoms within cutoff of any b-atom)``.

    The minimum is None when either set is empty. Chunked over ``a`` so the
    pairwise block stays small for a few thousand atoms on each side.
    """
    mask = np.zeros(a.shape[0], dtype=bool)
    if a.size == 0 or b.size == 0:
        return None, mask
    best = np.inf
    for start in range(0, a.shape[0], chunk):
        block = a[start : start + chunk]
        delta = block[:, None, :] - b[None, :, :]
        dist = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
        nearest = dist.min(axis=1)
        best = min(best, float(nearest.min()))
        mask[start : start + chunk] = nearest < cutoff
    return best, mask


def buried_area(a: AtomSet, b: AtomSet) -> float:
    """SASA(a) + SASA(b) - SASA(a+b), in A^2 (0 for non-touching sets).

    Only atoms within ``SASA_INFLUENCE`` of the other set are evaluated -- every
    other atom contributes the same SASA on both sides of the subtraction, so
    the value is the full difference, just without computing the cancelling
    terms (the design is ~370 residues and this runs six times).
    """
    if len(a) == 0 or len(b) == 0:
        return 0.0
    _, near_a = pair_stats(a.coords, b.coords, SASA_INFLUENCE)
    _, near_b = pair_stats(b.coords, a.coords, SASA_INFLUENCE)
    idx_a, idx_b = np.flatnonzero(near_a), np.flatnonzero(near_b)
    if idx_a.size == 0 or idx_b.size == 0:
        return 0.0

    free = sasa(a.coords, a.radii, idx_a).sum() + sasa(b.coords, b.radii, idx_b).sum()
    coords = np.concatenate([a.coords, b.coords])
    radii = np.concatenate([a.radii, b.radii])
    complexed = sasa(coords, radii, np.concatenate([idx_a, idx_b + len(a)])).sum()
    return max(float(free - complexed), 0.0)


# ---------------------------------------------------------------------------
# Per-design computation
# ---------------------------------------------------------------------------


def split_list(value: str) -> list[str]:
    """Comma-joined CLI list -> items, dropping blanks."""
    return [v.strip() for v in value.split(",") if v.strip()]


def parse_hotspots(value: str) -> list[str]:
    """``'A47, B42'`` -> ``['A47', 'B42']`` (chain letters + residue number)."""
    out: list[str] = []
    for token in split_list(value):
        m = re.fullmatch(r"([A-Za-z]{1,4})\s*(-?\d+)[A-Za-z]?", token)
        if not m:
            raise ValueError(f"cannot parse hotspot {token!r} (expected e.g. 'A47')")
        out.append(f"{m.group(1)}{int(m.group(2))}")
    return out


def load_structure(path: Path) -> gemmi.Structure:
    """Read a PDB or mmCIF file into a single-model, single-conformer structure."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")
    return st


def _fmt(value: object) -> str:
    """TSV cell: NA (None) becomes empty, floats get 3 decimals."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def compute(args: "Args") -> tuple[dict[str, object], gemmi.Structure]:
    """Every ringfit metric for one design, plus the design moved into the
    reference frame. Raises on anything unusable (the caller records the message
    as data)."""
    target_chains = split_list(args.target_chains)
    ref_target_chains = split_list(args.ref_target_chains)
    if len(target_chains) != 2:
        raise ValueError(
            f"--target-chains needs exactly 2 chains (the two protomers of the "
            f"seam), got {target_chains}"
        )
    if len(ref_target_chains) != len(target_chains):
        raise ValueError(
            f"--ref-target-chains {ref_target_chains} must have the same length as "
            f"--target-chains {target_chains}"
        )
    hotspots = parse_hotspots(args.hotspots)
    if args.hotspots_numbering not in ("design", "ref"):
        raise ValueError(
            f"--hotspots-numbering must be 'design' or 'ref', got "
            f"{args.hotspots_numbering!r}"
        )
    # Hotspots are matched against the DESIGN's residue labels; in 'ref' they are
    # given in the reference's numbering and translated once the match is known.
    hotspots_in_ref_numbering = args.hotspots_numbering == "ref"
    lipid_resnames = {r.upper() for r in split_list(args.lipid_resnames)}

    design = load_structure(args.design)
    ref = load_structure(args.ref)
    ref.remove_waters()
    d_model, r_model = design[0], ref[0]

    d_chains = protein_chain_names(d_model)
    missing = [c for c in target_chains if c not in d_chains]
    if missing:
        raise ValueError(
            f"target chain(s) {','.join(missing)} absent from the design "
            f"(protein chains: {','.join(d_chains) or 'none'})"
        )
    if args.binder_chains.strip().lower() == "auto":
        binder_chains = [c for c in d_chains if c not in target_chains]
    else:
        binder_chains = split_list(args.binder_chains)
        absent = [c for c in binder_chains if c not in d_chains]
        if absent:
            raise ValueError(
                f"binder chain(s) {','.join(absent)} absent from the design "
                f"(protein chains: {','.join(d_chains) or 'none'})"
            )
    if not binder_chains:
        raise ValueError(
            f"no binder chain in the design: every protein chain "
            f"({','.join(d_chains) or 'none'}) is a target chain"
        )
    both = [c for c in binder_chains if c in target_chains]
    if both:
        raise ValueError(
            f"chain(s) {','.join(both)} given as both target and binder chains"
        )

    r_chains = protein_chain_names(r_model)
    absent = [c for c in ref_target_chains if c not in r_chains]
    if absent:
        raise ValueError(
            f"reference target chain(s) {','.join(absent)} absent from the reference "
            f"(protein chains: {','.join(r_chains) or 'none'})"
        )

    # 1. Superimpose the design's target chains onto the reference's. Residues
    #    are paired ordinally when the chains are the same length (a design
    #    renumbered from 1 still matches) and by residue number otherwise; the
    #    pairing is then checked residue identity by residue identity.
    moving: list[np.ndarray] = []
    fixed: list[np.ndarray] = []
    offsets: list[int] = []
    n_pairs = n_same = 0
    ref_to_design: dict[str, str] = {}
    for d_chain, r_chain in zip(target_chains, ref_target_chains):
        d_ca = ca_residues(d_model, d_chain)
        r_ca = ca_residues(r_model, r_chain)
        pairs, used_mode = match_residues(d_ca, r_ca, args.resnum_match)
        if not pairs:
            raise ValueError(
                f"no residue matched between design chain {d_chain} "
                f"({len(d_ca)} CA) and reference chain {r_chain} ({len(r_ca)} CA) "
                f"under {used_mode} matching"
            )
        same = sum(d_ca.codes[i] == r_ca.codes[j] for i, j in pairs)
        frac = same / len(pairs)
        if frac < SEQ_MATCH_MIN:
            print(
                f"WARNING: design chain {d_chain} <-> reference chain {r_chain}: "
                f"only {frac:.2f} of the {len(pairs)} matched residues agree in "
                f"identity ({used_mode} matching). The superposition -- and every "
                f"metric under it -- is probably wrong: check the chain pairing "
                f"and try --resnum-match ordinal/resnum.",
                file=sys.stderr,
            )
        n_pairs += len(pairs)
        n_same += same
        moving += [d_ca.coords[i] for i, _ in pairs]
        fixed += [r_ca.coords[j] for _, j in pairs]
        offsets += [r_ca.resnums[j] - d_ca.resnums[i] for i, j in pairs]
        ref_to_design.update(
            {
                f"{r_chain}{r_ca.resnums[j]}": f"{d_chain}{d_ca.resnums[i]}"
                for i, j in pairs
            }
        )
    rot, trans, align_rmsd = kabsch(np.asarray(moving), np.asarray(fixed))
    seq_match_frac = n_same / n_pairs
    resnum_offset = int(np.bincount(np.asarray(offsets) - min(offsets)).argmax()) + min(
        offsets
    )

    if hotspots_in_ref_numbering:
        converted: list[str] = []
        for h in hotspots:
            if h not in ref_to_design:
                print(
                    f"WARNING: --hotspots-in-ref-numbering: reference residue {h} is "
                    f"not among the matched target residues; left as is",
                    file=sys.stderr,
                )
            converted.append(ref_to_design.get(h, h))
        hotspots = converted

    transform = gemmi.Transform()
    transform.mat.fromlist(rot.tolist())
    transform.vec.fromlist(trans.tolist())
    d_model.transform_pos_and_adp(transform)

    binder = protein_atoms(d_model, binder_chains)
    t1 = protein_atoms(d_model, target_chains[:1])
    t2 = protein_atoms(d_model, target_chains[1:2])
    targets = protein_atoms(d_model, target_chains)
    if len(binder) == 0:
        raise ValueError(f"binder chain(s) {','.join(binder_chains)} have no atoms")

    # 2. Ring clashes: the binder, now in the reference frame, against every
    #    reference protein chain that is not one of the two it was fitted on.
    other_chains = [c for c in r_chains if c not in ref_target_chains]
    ring = protein_atoms(r_model, other_chains)
    if len(ring) == 0:
        print(
            f"WARNING: the reference has no protein chain outside "
            f"{','.join(ref_target_chains)}; the ring check is vacuous",
            file=sys.stderr,
        )
    min_dist_ring, clash_mask = pair_stats(
        binder.coords, ring.coords, args.clash_cutoff
    )
    n_clash = int(clash_mask.sum())
    n_clash_res = len(binder.residues_in(clash_mask))

    # 3. Lipid belt: same idea, against the reference's membrane HETATMs.
    if lipid_resnames:
        lipids = het_atoms(r_model, lipid_resnames)
        if len(lipids) == 0:
            print(
                f"WARNING: no HETATM with comp id in "
                f"{','.join(sorted(lipid_resnames))} in the reference",
                file=sys.stderr,
            )
        min_dist_lipid, lipid_mask = pair_stats(
            binder.coords, lipids.coords, LIPID_CLASH_CUTOFF
        )
        lipid_clash: int | None = int(lipid_mask.sum())
    else:
        min_dist_lipid, lipid_clash = None, None

    # 4. Buried surface area, on the design: does the binder straddle the seam?
    bsa_t1 = buried_area(t1, binder)
    bsa_t2 = buried_area(t2, binder)
    bsa_total = buried_area(targets, binder)
    hi = max(bsa_t1, bsa_t2)
    bridge_ratio = min(bsa_t1, bsa_t2) / hi if hi > 0 else None

    # 5. Contacts, per target chain and against the hotspot selection.
    _, mask_t1 = pair_stats(t1.coords, binder.coords, args.contact_cutoff)
    _, mask_t2 = pair_stats(t2.coords, binder.coords, args.contact_cutoff)
    _, mask_targets = pair_stats(targets.coords, binder.coords, args.contact_cutoff)
    contacted = set(targets.residues_in(mask_targets))

    if hotspots:
        unknown = [h for h in hotspots if h not in set(targets.res_labels)]
        if unknown:
            print(
                f"WARNING: hotspot(s) {','.join(unknown)} not found in the design's "
                f"target chains; counted as not contacted",
                file=sys.stderr,
            )
        hits = [h for h in hotspots if h in contacted]
        hotspot_recall: float | None = len(hits) / len(hotspots)
        hotspot_hits: str | None = ",".join(hits)
    else:
        hotspot_recall, hotspot_hits = None, None

    metrics = {
        "align_rmsd": align_rmsd,
        "seq_match_frac": seq_match_frac,
        "resnum_offset": resnum_offset,
        "n_clash": n_clash,
        "n_clash_res": n_clash_res,
        "min_dist_ring": min_dist_ring,
        "lipid_clash": lipid_clash,
        "min_dist_lipid": min_dist_lipid,
        "bsa_t1": bsa_t1,
        "bsa_t2": bsa_t2,
        "bsa_total": bsa_total,
        "bridge_ratio": bridge_ratio,
        "n_contact_res_t1": len(t1.residues_in(mask_t1)),
        "n_contact_res_t2": len(t2.residues_in(mask_t2)),
        "hotspot_recall": hotspot_recall,
        "hotspot_hits": hotspot_hits,
        "binder_len": binder.n_res,
        "binder_chains": ",".join(binder_chains),
    }
    return metrics, design


def write_result(
    result_tsv: Path, name: str, status: str, aligned_path: str, metrics: dict
) -> None:
    result_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(result_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(RESULT_COLUMNS)
        writer.writerow(
            [name, status, aligned_path]
            + [_fmt(metrics.get(c)) for c in METRIC_COLUMNS]
        )


class Args(argparse.Namespace):
    name: str
    design: Path
    ref: Path
    target_chains: str
    binder_chains: str
    ref_target_chains: str
    hotspots: str
    hotspots_numbering: str
    resnum_match: str
    clash_cutoff: float
    contact_cutoff: float
    lipid_resnames: str
    out_pdb: Path
    result_tsv: Path


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Score one design's fit into a reference oligomeric assembly."
    )
    ap.add_argument("--name", required=True)
    ap.add_argument(
        "--design", type=Path, required=True, help="Design structure (PDB or mmCIF)."
    )
    ap.add_argument(
        "--ref",
        type=Path,
        required=True,
        help="Reference full assembly (PDB or mmCIF).",
    )
    ap.add_argument("--target-chains", default="A,B")
    ap.add_argument("--binder-chains", default="auto")
    ap.add_argument("--ref-target-chains", default="A,B")
    ap.add_argument("--hotspots", default="")
    ap.add_argument(
        "--hotspots-numbering",
        choices=("design", "ref"),
        default="design",
        help="Numbering --hotspots are written in: 'design' (default) or 'ref', "
        "which translates them through the design<->reference residue match.",
    )
    ap.add_argument(
        "--resnum-match",
        choices=MATCH_MODES,
        default="auto",
        help="How design residues are paired with reference residues for the "
        "superposition: ordinal (i-th to i-th), resnum (by residue number), or "
        "auto (ordinal when the chains are equally long). Default auto.",
    )
    ap.add_argument("--clash-cutoff", type=float, default=2.5)
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--lipid-resnames", default="PLM,LPP,L8Z")
    ap.add_argument(
        "--out-pdb", type=Path, required=True, help="Superposed design to write."
    )
    ap.add_argument(
        "--result-tsv", type=Path, required=True, help="Per-design result TSV to write."
    )
    args = ap.parse_args(namespace=Args())

    try:
        metrics, superposed = compute(args)
        # Keep the design in the reference frame, so the fit can be eyeballed
        # against the full assembly in PyMOL/ChimeraX.
        args.out_pdb.parent.mkdir(parents=True, exist_ok=True)
        superposed.setup_entities()
        superposed.write_pdb(str(args.out_pdb))
        seq_match = float(metrics["seq_match_frac"])  # type: ignore[arg-type]
        print(
            f"{args.name}: align_rmsd {metrics['align_rmsd']:.3f} A, "
            f"seq_match_frac {seq_match:.3f}, resnum_offset "
            f"{metrics['resnum_offset']}, n_clash {metrics['n_clash']}, bsa "
            f"{metrics['bsa_t1']:.0f}/{metrics['bsa_t2']:.0f} A^2, bridge_ratio "
            f"{_fmt(metrics['bridge_ratio'])}"
        )
        # A suspect residue match is reported, never swallowed: every metric is
        # still written, but the status says the superposition is not trustworthy.
        status = (
            "OK"
            if seq_match >= SEQ_MATCH_MIN
            else f"warn: seq_match_frac {seq_match:.2f}"
        )
        write_result(args.result_tsv, args.name, status, str(args.out_pdb), metrics)
    except Exception as e:  # noqa: BLE001 - errors are recorded as data
        print(f"{args.name}: ERROR {e}", file=sys.stderr)
        write_result(args.result_tsv, args.name, f"error: {e}", "", {})


if __name__ == "__main__":
    main()
