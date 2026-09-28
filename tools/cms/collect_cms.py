#!/usr/bin/env python3
"""
Collect cms results into the table.

Reads the per-design <name>.tsv files that cms_worker.py wrote under
<run_dir>/<table>/cms[_<label>]/ and merges them back as <prefix>_status /
<prefix>_path (the per-residue CMS table: side, chain, resnum, resname, cms -- the
interface hotspots) plus one <prefix>_<field> column per metric:

    target          CMS on the target side, A^2 (the conventional CMS number)
    binder          CMS on the binder side, A^2
    sc              Lawrence & Colman shape complementarity (0-1; ~0.5-0.75 good)
    sc_area         trimmed interface area of both sides, A^2
    sc_median_dist  median interface separation, A
    max_target      max possible CMS of the target alone   (--max-cms only)
    max_binder      max possible CMS of the binder alone   (--max-cms only)
    frac_target     target / max_target                    (--max-cms only)
    frac_binder     binder / max_binder                    (--max-cms only)
    n_atoms_binder  heavy atoms scored on each side
    n_atoms_target
    n_radius0_*     atoms with no radius-table entry (invisible to the surface)
    device          cuda or cpu -- what actually ran
    seconds         wall time of the design's calculation

NA (empty) rather than 0 marks a field that does not apply: the design failed, no TSV
was written, or the metric was switched off. A CMS of 0 is a real value -- no atoms
within 8 A of the partner -- and is not a failure.

Usage:
    sapia collect cms outputs/20260928_100140_7ojg_binder --table table1 -l pred
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (cms_<name>).
FLOAT_COLUMNS = [
    "target",
    "binder",
    "sc",
    "sc_area",
    "sc_median_dist",
    "max_target",
    "max_binder",
    "frac_target",
    "frac_binder",
    "seconds",
]
INT_COLUMNS = ["n_atoms_binder", "n_atoms_target", "n_radius0_binder", "n_radius0_target"]
STR_COLUMNS = ["device"]
RESULT_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    data = _empty()
    for col in FLOAT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = float(value)
    for col in INT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value))
    for col in STR_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    return data


def collect_cms(ctx: CollectCtx) -> CollectEach:
    """Per-design cms collector. Variants -- another structure column or chain split
    off the same table -- are distinguished via --dir-label, matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"
        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype={"device": str})
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row = result_df.iloc[0]
        status = str(row["status"])
        if status != "OK":
            yield Collected(status=status, path="", data=_empty())
            return

        path = _present(row, "path")
        yield Collected(
            status=status,
            path="" if path is None else str(path),
            data=_fields(row),
        )

    return one
