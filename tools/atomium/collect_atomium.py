"""Collect AtomiUM FASTA outputs into a table.

Scans the output directory for the batched subgroup dirs, parses the FASTAs and
writes one child row per designed sequence.

Lineage is derived from the directory structure:

    <output_dir>/grp_<g>/seqs/<fasta_stem>.fa

  * grp_<g>     = a batched subgroup run by atomium.sh (several designs sharing
    params, run in one atomium.py call).
  * fasta_stem  = the staged input filename AtomiUM processed. The submitter
    symlinks each input as ``<design_name>.pdb``, so the stem IS the immediate
    parent-table row -- stamped as ``parent_name``.

AtomiUM's FASTA is NOT ProteinMPNN's, so this does not share proteinmpnn's parser
(verified against ``atomium.py``, which writes the file directly):

    >_init_sequence
    MKT...                                     <- the native input sequence
    >gen_sample_0, temperature=0.1, seq_rec=0.43
    MRE...

There is no ``score`` or ``global_score`` -- AtomiUM reports only sequence recovery
against the input, so a design cannot be ranked by model likelihood here the way an
mpnn ``score`` column allows. Two consequences worth knowing when reading the table:

  * ``seq_rec`` is similarity to the INPUT sequence, not a quality score. On a de
    novo backbone whose input sequence is poly-glycine it is close to meaningless;
    it only carries information when redesigning a real sequence.
  * ``sample`` is NOT unique per row. AtomiUM numbers samples ``n % num_seq_per_target``
    and walks temperatures with ``n // num_seq_per_target``, so with several
    ``--sampling-temp`` values the same ``sample`` id recurs once per temperature.
    The row key therefore uses the running entry index (``<parent>_a<i>``), and
    (``sample``, ``temperature``) together identify a draw.

Each row carries only ``parent_name``; ancestor values are resolved on demand by
walking the lineage, so nothing is propagated here.

Usage:
    sapia collect atomium outputs/<run> --table <the table the run reserved>
"""

from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import pandas as pd

from prosapia.core import (
    Collected,
    CollectCtx,
    CollectEach,
    DesignCtx,
)

# The header AtomiUM gives the native sequence it echoes first in every FASTA.
INIT_HEADER = "_init_sequence"

# Trailing ``key=value`` fields on a generated header.
HEADER_FIELDS: Dict[str, type] = {
    "temperature": float,
    "seq_rec": float,
}


def parse_atomium_header(header: str) -> Dict[str, Any]:
    """Parse ``gen_sample_0, temperature=0.1, seq_rec=0.43`` into columns.

    The sample id is the ``gen_sample_<k>`` prefix rather than a ``key=value`` pair.
    Unparseable fields become NA rather than raising: a partially written header
    should still collect its sequence.
    """
    out: Dict[str, Any] = {"sample": pd.NA}
    out.update({k: pd.NA for k in HEADER_FIELDS})
    for part in header.split(","):
        part = part.strip()
        if part.startswith("gen_sample_"):
            try:
                out["sample"] = int(part[len("gen_sample_") :])
            except ValueError:
                pass
            continue
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key, value = key.strip(), value.strip()
        if key in HEADER_FIELDS:
            try:
                out[key] = HEADER_FIELDS[key](value)
            except ValueError:
                out[key] = pd.NA
    return out


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    header = ""
    seq_lines: List[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header:
        entries.append((header, "".join(seq_lines)))
    return entries


def collect_atomium(ctx: CollectCtx) -> CollectEach:
    """Per-parent AtomiUM collector. atomium is a create tool: the framework iterates
    the ready parents and this mints one child row (``<parent>_a<i>``) per sampled
    sequence, carrying ``parent`` for lineage. The framework stamps status/path/
    parent_name from each Collected; child rows are discovered on disk, so re-running
    rebuilds them (idempotent)."""
    # Each grp_<g>/ output dir holds seqs/<design>.fa, where the staged input was
    # symlinked as <design>.pdb -- so the FASTA stem IS the parent-table row name.
    fasta_by_parent: Dict[str, Path] = {}
    for subdir in sorted(p for p in ctx.out_dir.iterdir() if p.is_dir()):
        if subdir.name in ("atomium_logs", "atomium_tasks"):
            continue
        seqs_dir = subdir / "seqs"
        if not seqs_dir.is_dir():
            continue
        for fasta_path in seqs_dir.glob("*.fa"):
            fasta_by_parent[fasta_path.stem] = fasta_path

    def one(d: DesignCtx) -> Iterable[Collected]:
        fasta_path = fasta_by_parent.get(d.name)
        if fasta_path is None:
            print(f"{d.name}: no FASTA found (skipping)")
            return

        for i, (header, sequence) in enumerate(parse_fasta(fasta_path)):
            # The first entry is AtomiUM's echo of the native input sequence. Matched
            # by its header rather than by position, so a future extra preamble entry
            # cannot silently shift which row is dropped.
            if header.strip() == INIT_HEADER:
                continue

            # Bare column names; the driver leaf-prefixes them.
            data: Dict[str, Any] = {"iteration": i, "sequence": sequence}
            data.update(parse_atomium_header(header))
            yield Collected(
                name=f"{d.name}_a{i}",
                parent=d.name,
                path=fasta_path,
                data=data,
            )

    return one
