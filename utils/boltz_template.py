#!/usr/bin/env python3
"""
Build a forced-template mmCIF for Boltz -- and the ``templates:`` YAML that goes
with it -- from a subset of chains of a reference structure.

NOT a sapia tool (no spec.py, so `discover()` never picks it up). It is a
workstation-side utility, run once per template while staging inputs. The
workstation mounts the user tools dir at its own absolute path and exports
``$PROSAPIA_TOOLS_DIR``, so:

    sapia modal-shell --cmd 'python $PROSAPIA_TOOLS_DIR/_utils/boltz_template.py \\
        --fetch 7OJG --out /runs/inputs/7ojg_tmpl_kabc.cif \\
        --chains K,A,B,C --rename-to A,B,C,D --residues K:19-59+107-155 \\
        --yaml /runs/inputs/7ojg_tmpl_kabc.yaml'

Write outputs under ``/runs/`` -- every ``modal-shell`` call is a FRESH
container, so a file left in /tmp does not exist for the next one.

It exists because the template CIF is the one input in this workbench whose
failure modes are SILENT and EXPENSIVE, and because writing it ad hoc gets the
same four things wrong every time (see `.claude/skills/boltz/SKILL.md`):

1. **Boltz names template chains by ``label_asym_id``, not the auth chain ID.**
   A CIF written by gemmi's own ``setup_entities()`` labels subchains ``Axp``,
   ``Bxp``, ... and ``template_id: [A, ...]`` then dies with *"Template chain A
   is not one of the protein chains"*. Here every chain is given
   ``label_asym_id == auth_asym_id`` by construction.

2. **``_entity_poly_seq`` must equal the MODELLED residue count of every
   entity.** Boltz's ``parse_polymer`` uses it as the token list, and a declared
   residue absent from the model still consumes a token index -- so a
   too-long declaration shifts every residue index downstream and nothing
   errors. The classic way in: residues cloned from a source file keep their
   source ``subchain`` label, ``setup_entities()`` merges two chains into one
   entity, and ``full_sequence`` is taken from whichever came first (a real
   case: ``_entity_poly_seq`` = 49 for 41-residue chains). Here
   ``setup_entities()`` is never called -- one entity is built per output chain
   with ``full_sequence`` taken from the residues actually written, so the two
   cannot disagree.

3. **Boltz falls back to a sequence SEARCH unless ``chain_id`` and
   ``template_id`` are both given and equal length**, which with sequence-
   identical chains (any homo-oligomer) can assign the wrong template chain and
   destroy the geometry without a word. The emitted YAML always writes both.

4. Ligands, waters and a second model quietly become extra entities. Only the
   first model is used and only tabulated amino/nucleic acids are kept
   (``--keep-het`` to opt out).

Every write is followed by ``verify()``, which RE-READS the written file as raw
mmCIF categories -- not through gemmi's structure view, which would paper over
exactly the mismatch being looked for -- and counts ``_entity_poly_seq.mon_id``
rows per entity against the modelled residues of its subchains. A failed check
exits 1 and says which chain and which count. Run it on a CIF from anywhere with
``--check <file>``.

What it deliberately does NOT do: decide which chains to template. Do not
template the binder chain, and remember a forced template STEERS, it does not
constrain -- superpose the predicted target chains back onto this file
afterwards and drop the rows where the target did not land.

Modes
-----
build   --out with --chains (and --in or --fetch)
check   --check <file.cif> [--expect-chains A,B,C]  -- verify only, writes nothing
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import gemmi

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"


# --------------------------------------------------------------------------
# argument parsing
# --------------------------------------------------------------------------


def split_list(value: str) -> list[str]:
    """'A, B ,C' -> ['A', 'B', 'C']; '' -> []."""
    return [tok.strip() for tok in value.split(",") if tok.strip()]


def parse_residues(value: str) -> dict[str, list[tuple[int, int]]]:
    """'K:19-59+107-155,A:36-149' -> {'K': [(19,59), (107,155)], 'A': [(36,149)]}.

    Keyed by the SOURCE chain ID, applied before any rename. Bounds are
    inclusive and compared against the source residue numbers, so they read the
    same as the numbering in the reference file (and in the paper).
    """
    ranges: dict[str, list[tuple[int, int]]] = {}
    for spec in split_list(value):
        if spec.count(":") != 1:
            raise ValueError(
                f"--residues entry {spec!r} must be '<chain>:<first>-<last>[+<first>-<last>...]'"
            )
        chain, segments = (tok.strip() for tok in spec.split(":"))
        parsed: list[tuple[int, int]] = []
        for seg in segments.split("+"):
            first_str, _, last_str = seg.strip().partition("-")
            if not first_str or not last_str:
                raise ValueError(f"--residues segment {seg!r} must be '<first>-<last>'")
            first, last = int(first_str), int(last_str)
            if last < first:
                raise ValueError(f"--residues segment {seg!r} ends before it starts")
            parsed.append((first, last))
        ranges.setdefault(chain, []).extend(parsed)
    return ranges


def resolve_chains(chains: str, rename_to: str) -> list[tuple[str, str]]:
    """[(source_id, output_id), ...] in the order given; rename is positional."""
    kept = split_list(chains)
    renamed = split_list(rename_to)
    if not kept:
        raise ValueError("--chains is empty")
    if renamed and len(renamed) != len(kept):
        raise ValueError(
            f"--rename-to has {len(renamed)} IDs but --chains has {len(kept)}"
        )
    out_ids = renamed or kept
    if len(set(out_ids)) != len(out_ids):
        raise ValueError(f"output chain IDs are not unique: {','.join(out_ids)}")
    return list(zip(kept, out_ids))


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------


def fetch_rcsb(pdb_id: str, dest: Path) -> Path:
    """Download an entry from RCSB to ``dest``, uncompressed.

    Done here rather than with ``modal volume put`` because this network
    intercepts TLS to Modal's blob storage and the upload fails; fetching from
    inside the workstation sidesteps it entirely.
    """
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as resp, gzip.GzipFile(
        fileobj=resp
    ) as gz, dest.open("wb") as out:
        shutil.copyfileobj(gz, out)
    print(f"fetched {url} -> {dest} ({dest.stat().st_size} bytes)")
    return dest


def load_structure(path: Path) -> gemmi.Structure:
    """Read a PDB or mmCIF file, single conformer, no hydrogens."""
    if not path.exists():
        raise FileNotFoundError(f"structure missing: {path}")
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    if len(st) == 0:
        raise ValueError(f"no model in {path}")
    return st


def is_polymer_residue(res: gemmi.Residue) -> bool:
    """True for a tabulated amino-acid or nucleic-acid residue (MSE included)."""
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and (info.is_amino_acid() or info.is_nucleic_acid()))


def one_letter(res: gemmi.Residue) -> str:
    info = gemmi.find_tabulated_residue(res.name)
    code = info.one_letter_code if info else ""
    return code.upper() if code else "X"


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------


@dataclass
class ChainReport:
    source: str
    out: str
    n_res: int
    first: str
    last: str
    sequence: str


@dataclass
class BuildResult:
    structure: gemmi.Structure
    chains: list[ChainReport] = field(default_factory=list)


def select(
    st: gemmi.Structure,
    groups: list[tuple[str, str]],
    ranges: dict[str, list[tuple[int, int]]],
    strip_het: bool,
    renumber: bool,
    model_index: int,
) -> BuildResult:
    """Copy the requested chains into a fresh structure, one entity per chain.

    A requested chain that is not in the file, or a residue range that selects
    nothing, is an ERROR -- never a shorter template. A silently partial
    template is well-formed, folds, and is wrong.
    """
    if model_index >= len(st):
        raise ValueError(f"--model {model_index} but the file has {len(st)} model(s)")
    model = st[model_index]
    present = {chain.name: chain for chain in model}
    unknown = [src for src, _ in groups if src not in present]
    if unknown:
        raise ValueError(
            f"chain(s) {','.join(unknown)} not in the structure "
            f"(have: {','.join(sorted(present))})"
        )

    out_st = gemmi.Structure()
    out_st.name = st.name or "template"
    out_st.spacegroup_hm = "P 1"
    out_st.cell = gemmi.UnitCell()
    out_model = gemmi.Model("1")
    result = BuildResult(structure=out_st)

    for src, out_id in groups:
        segments = ranges.get(src)
        out_chain = gemmi.Chain(out_id)
        for res in present[src]:
            if strip_het and not is_polymer_residue(res):
                continue
            num = res.seqid.num
            if segments and not any(lo <= num <= hi for lo, hi in segments):
                continue
            copy = res.clone()
            # The whole point: no source subchain label survives into the new
            # file, so nothing can merge two chains into one entity behind our
            # back. label_asym_id is set explicitly to the auth chain ID.
            copy.subchain = out_id
            copy.het_flag = "A"
            copy.label_seq = None
            if renumber:
                copy.seqid = gemmi.SeqId(len(out_chain) + 1, " ")
            out_chain.add_residue(copy)
        if len(out_chain) == 0:
            raise ValueError(
                f"chain {src} selected no residues"
                + (f" in range {ranges_str(segments)}" if segments else "")
            )
        if segments:
            want = sum(hi - lo + 1 for lo, hi in segments)
            if len(out_chain) != want:
                print(
                    f"note: chain {src} -> {out_id}: {len(out_chain)} residues for a "
                    f"range spanning {want} -- the source has gaps or unmodelled "
                    f"residues there. The template declares what is MODELLED, which "
                    f"is correct; check the gap is where you think it is.",
                    file=sys.stderr,
                )
        out_model.add_chain(out_chain)
        result.chains.append(
            ChainReport(
                source=src,
                out=out_id,
                n_res=len(out_chain),
                first=str(out_chain[0].seqid.num),
                last=str(out_chain[len(out_chain) - 1].seqid.num),
                sequence="".join(one_letter(r) for r in out_chain),
            )
        )

    out_st.add_model(out_model)
    set_entities(out_st)
    return result


def ranges_str(segments: list[tuple[int, int]]) -> str:
    return "+".join(f"{lo}-{hi}" for lo, hi in segments)


def set_entities(st: gemmi.Structure) -> None:
    """One polymer entity per chain, ``full_sequence`` = the modelled residues.

    ``setup_entities()`` is deliberately NOT used. It would relabel subchains
    (``Axp``) and deduplicate sequence-identical chains into a shared entity
    whose ``full_sequence`` comes from one of them -- the two ways this file
    goes wrong. Duplicating identical entities is redundant mmCIF, not invalid,
    and it makes the per-chain token count impossible to get wrong.
    """
    st.entities.clear()
    for chain in st[0]:
        ent = gemmi.Entity(chain.name)
        ent.entity_type = gemmi.EntityType.Polymer
        ent.polymer_type = gemmi.PolymerType.PeptideL
        ent.subchains = [chain.name]
        ent.full_sequence = [res.name for res in chain]
        st.entities.append(ent)
    # label_seq_id then runs 1..N per chain, matching the token indices Boltz
    # builds from _entity_poly_seq.
    st.assign_label_seq_id(True)


def write_cif(st: gemmi.Structure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    doc = st.make_mmcif_document()
    doc.write_file(str(out))


# --------------------------------------------------------------------------
# verify -- reads the written file back as raw mmCIF
# --------------------------------------------------------------------------


def verify(path: Path, expect_chains: list[str]) -> tuple[bool, list[str]]:
    """Check the invariants Boltz depends on. Returns (ok, report lines).

    Reads ``_atom_site``, ``_entity_poly_seq`` and ``_struct_asym`` directly:
    going through ``gemmi.read_structure`` would re-derive what is being
    checked and pass on a file Boltz rejects.
    """
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    lines: list[str] = []
    failures: list[str] = []

    atoms = block.find(
        "_atom_site.",
        ["label_asym_id", "auth_asym_id", "label_entity_id", "label_seq_id", "label_comp_id"],
    )
    if not len(atoms):
        return False, ["FAIL: no _atom_site rows"]

    # modelled residues per label_asym_id, in file order
    residues: dict[str, list[tuple[str, str]]] = {}
    label_to_auth: dict[str, set[str]] = {}
    label_to_entity: dict[str, set[str]] = {}
    for row in atoms:
        label_asym, auth_asym, entity, label_seq, comp = (row[i] for i in range(5))
        label_to_auth.setdefault(label_asym, set()).add(auth_asym)
        label_to_entity.setdefault(label_asym, set()).add(entity)
        seen = residues.setdefault(label_asym, [])
        if not seen or seen[-1] != (label_seq, comp):
            seen.append((label_seq, comp))

    # 1. label_asym_id == auth_asym_id
    mismatched = {
        lab: sorted(auth) for lab, auth in label_to_auth.items() if auth != {lab}
    }
    if mismatched:
        failures.append(
            "label_asym_id != auth_asym_id for "
            + "; ".join(f"{lab} (auth {','.join(a)})" for lab, a in mismatched.items())
            + " -- Boltz addresses template chains by label_asym_id, so template_id "
            "must name these, or the run dies with 'not one of the protein chains'"
        )
    else:
        lines.append(f"OK  label_asym_id == auth_asym_id for all {len(residues)} chains")

    # 2. _entity_poly_seq rows per entity == modelled residues of its subchains
    poly = block.find("_entity_poly_seq.", ["entity_id", "num", "mon_id"])
    declared: dict[str, list[str]] = {}
    for row in poly:
        declared.setdefault(row[0], []).append(row[2])
    if not declared:
        failures.append(
            "no _entity_poly_seq rows at all -- Boltz's parse_polymer builds its "
            "token list from this category, so there is nothing to template"
        )
    for label_asym, res_list in sorted(residues.items()):
        entities = label_to_entity[label_asym]
        if len(entities) != 1:
            failures.append(
                f"chain {label_asym} spans entities {','.join(sorted(entities))}"
            )
            continue
        entity = entities.pop()
        mon_ids = declared.get(entity, [])
        modelled = [comp for _, comp in res_list]
        if len(mon_ids) != len(modelled):
            failures.append(
                f"chain {label_asym} (entity {entity}): _entity_poly_seq declares "
                f"{len(mon_ids)} residues but {len(modelled)} are modelled -- every "
                "declared-but-absent residue still consumes a token index, so all "
                "downstream residue indices shift silently. A shared entity between "
                "two sequence-identical chains is the usual cause"
            )
        elif mon_ids != modelled:
            first = next(
                i for i, (a, b) in enumerate(zip(mon_ids, modelled)) if a != b
            )
            failures.append(
                f"chain {label_asym} (entity {entity}): _entity_poly_seq and the "
                f"modelled residues differ at position {first + 1} "
                f"({mon_ids[first]} declared, {modelled[first]} modelled)"
            )
        else:
            lines.append(
                f"OK  chain {label_asym}: _entity_poly_seq = {len(mon_ids)} residues, "
                "identical to what is modelled"
            )

    # 3. label_seq_id contiguous from 1 (the token index Boltz builds)
    for label_asym, res_list in sorted(residues.items()):
        nums = [r[0] for r in res_list]
        if nums != [str(i + 1) for i in range(len(nums))]:
            failures.append(
                f"chain {label_asym}: label_seq_id is not 1..{len(nums)} "
                f"(starts {','.join(nums[:4])})"
            )

    # 4. the chains the YAML will name are actually there
    missing = [c for c in expect_chains if c not in residues]
    if missing:
        failures.append(
            f"template_id names {','.join(missing)} but the CIF has "
            f"{','.join(sorted(residues))}"
        )
    elif expect_chains:
        lines.append(f"OK  template chains present: {','.join(expect_chains)}")

    for f in failures:
        lines.append(f"FAIL {f}")
    return not failures, lines


# --------------------------------------------------------------------------
# YAML
# --------------------------------------------------------------------------


def template_yaml(
    cif: str, chain_ids: list[str], template_ids: list[str], force: bool, threshold: float
) -> str:
    """The ``templates:`` block, at column 0 -- boltz splices it verbatim.

    ``chain_id`` (prediction chains) and ``template_id`` (chains in this CIF)
    are ALWAYS both written and equal length: that is what makes Boltz zip them
    positionally instead of falling back to a sequence search.
    """
    if len(chain_ids) != len(template_ids):
        raise ValueError(
            f"--predict-chains has {len(chain_ids)} IDs but the template has "
            f"{len(template_ids)} chains -- they must zip 1:1"
        )
    body = [
        "templates:",
        f"  - cif: {cif}",
        f"    chain_id: [{', '.join(chain_ids)}]",
        f"    template_id: [{', '.join(template_ids)}]",
    ]
    if force:
        body.append("    force: true")
        body.append(f"    threshold: {threshold}")
    return "\n".join(body) + "\n"


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


class Args(argparse.Namespace):
    pass


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Build a Boltz forced-template mmCIF and its templates: YAML.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  # 8-chain template from two protomers of a fetched assembly\n"
            "  boltz_template.py --fetch 7OJG --out /runs/inputs/tmpl.cif \\\n"
            "      --chains K,A,B,C --residues K:19-59+107-155 \\\n"
            "      --yaml /runs/inputs/tmpl.yaml\n\n"
            "  # verify a CIF written by something else\n"
            "  boltz_template.py --check /runs/inputs/tmpl.cif --expect-chains A,B,C,D\n"
        ),
    )
    ap.add_argument("--in", dest="src", type=Path, help="Source PDB/mmCIF (or --fetch).")
    ap.add_argument("--fetch", type=str, default="", help="PDB ID to download from RCSB.")
    ap.add_argument(
        "--fetch-to",
        type=Path,
        help="Where to keep the fetched entry (default: <out>'s dir, <ID>.cif).",
    )
    ap.add_argument("--out", type=Path, help="Template mmCIF to write.")
    ap.add_argument("--chains", type=str, default="", help="Source chain IDs, in order.")
    ap.add_argument(
        "--rename-to", type=str, default="", help="Output chain IDs, positional."
    )
    ap.add_argument(
        "--residues",
        type=str,
        default="",
        help="Per-SOURCE-chain residue ranges, e.g. 'K:19-59+107-155,A:36-149'.",
    )
    ap.add_argument(
        "--predict-chains",
        type=str,
        default="",
        help="chain_id in the YAML (chains of the PREDICTION). Default: the "
        "template's own chain IDs.",
    )
    ap.add_argument("--yaml", dest="yaml_out", type=Path, help="Write the YAML block here.")
    ap.add_argument(
        "--cif-path",
        type=str,
        default="",
        help="Path to name in the YAML's cif: field (default: --out, which is "
        "already a workstation path).",
    )
    ap.add_argument(
        "--threshold", type=float, default=2.0, help="force threshold (default 2.0)."
    )
    ap.add_argument(
        "--no-force",
        action="store_true",
        help="Omit force/threshold (template as a hint, not forced).",
    )
    ap.add_argument("--renumber", action="store_true", help="Renumber each chain 1..N.")
    ap.add_argument(
        "--keep-het",
        action="store_true",
        help="Keep ligands/ions/waters (default: polymer residues only).",
    )
    ap.add_argument("--model", type=int, default=0, help="Model index (default 0).")
    ap.add_argument("--check", type=Path, help="Verify an existing CIF and exit.")
    ap.add_argument(
        "--expect-chains", type=str, default="", help="With --check: chains that must exist."
    )
    return ap


def main() -> None:
    args = build_parser().parse_args(namespace=Args())

    if args.check:
        if not args.check.exists():
            # Each modal-shell call is a FRESH container: anything written to
            # /tmp by an earlier call is gone. Templates belong on the Volume.
            raise SystemExit(
                f"no such file: {args.check} -- if a previous call wrote it to /tmp, "
                "note every `sapia modal-shell` call is a new container; write "
                "templates under /runs/"
            )
        ok, lines = verify(args.check, split_list(args.expect_chains))
        print(f"--- {args.check}")
        print("\n".join(lines))
        raise SystemExit(0 if ok else 1)

    if not args.out:
        raise SystemExit("--out is required (or use --check)")
    if not args.src and not args.fetch:
        raise SystemExit("one of --in / --fetch is required")

    src = args.src
    if args.fetch:
        dest = args.fetch_to or args.out.parent / f"{args.fetch.upper()}.cif"
        src = fetch_rcsb(args.fetch, dest)

    groups = resolve_chains(args.chains, args.rename_to)
    built = select(
        load_structure(src),
        groups,
        parse_residues(args.residues),
        strip_het=not args.keep_het,
        renumber=args.renumber,
        model_index=args.model,
    )
    write_cif(built.structure, args.out)

    out_ids = [c.out for c in built.chains]
    predict = split_list(args.predict_chains) or out_ids
    block = template_yaml(
        cif=args.cif_path or str(args.out),
        chain_ids=predict,
        template_ids=out_ids,
        force=not args.no_force,
        threshold=args.threshold,
    )

    print(f"--- {args.out}")
    for c in built.chains:
        seq = c.sequence if len(c.sequence) <= 60 else c.sequence[:57] + "..."
        print(
            f"  {c.source} -> {c.out}: {c.n_res} res, {c.first}-{c.last}  {seq}"
        )
    ok, lines = verify(args.out, out_ids)
    print("\n".join(f"  {line}" for line in lines))

    if args.yaml_out:
        args.yaml_out.parent.mkdir(parents=True, exist_ok=True)
        args.yaml_out.write_text(block)
        print(f"--- {args.yaml_out}")
    else:
        print("--- templates: block (pass to boltz --template-yaml)")
    print(block, end="")

    if not ok:
        print(
            "\nNOT WRITTEN AS USABLE: the checks above failed. Fix before spending.",
            file=sys.stderr,
        )
        raise SystemExit(1)


if __name__ == "__main__":
    main()
