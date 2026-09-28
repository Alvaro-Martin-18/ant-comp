#!/usr/bin/env python3
"""
Collect mkcomplex results into the table.

Scans <run_dir>/<table>/mkcomplex/ for the per-design TSV files written by
mkcomplex.sh and merges the assembled complex back into the table as
<prefix>_status plus one column per field:

    sequence            the assembled multi-chain string, e.g. 'SEQA/SEQB/DESIGN'
                        -- feed it to a predictor with
                        `-i mkcomplex_sequence` so it folds the COMPLEX
    n_chains            number of chains in it (separator-delimited segments)
    chain_lens          per-chain residue counts, comma-joined, in chain order
    design_chain_index  0-based index of the design's own chain, i.e. which
                        chain is the binder (chain letter chr(ord('A') + index)
                        once a predictor maps segments onto letters)
    total_len           residues summed over all chains

<prefix>_path is always empty: mkcomplex writes no structure, only a sequence.
A design whose sequence was empty, already multi-chain, or whose fixed chains
carried a non-standard residue gets an 'error: ...' status and empty columns --
the task still exited 0, so the rest of the array collects normally.

Usage:
    sapia collect mkcomplex outputs/20260928_100140_7ojg_binder --table table1
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (mkcomplex_<name>).
INT_COLUMNS = ["n_chains", "design_chain_index", "total_len"]
STR_COLUMNS = ["sequence", "chain_lens"]
RESULT_COLUMNS = STR_COLUMNS + INT_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in RESULT_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _fields(row: pd.Series) -> dict[str, Any]:
    """One result row -> collected columns, NA-safe (a blank cell is a field that
    was not produced, not a zero)."""
    data = _empty()
    for col in STR_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = str(value)
    for col in INT_COLUMNS:
        if (value := _present(row, col)) is not None:
            data[col] = int(float(value))
    return data


def collect_mkcomplex(ctx: CollectCtx) -> CollectEach:
    """Per-design mkcomplex collector. The framework iterates ready designs and stamps
    status (keyed by the tool leaf); this reads one design's one-row <name>.tsv and adds
    the assembled sequence columns. Variants are distinguished via --dir-label, matching
    the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"

        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t", dtype={"sequence": str})
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row_data = result_df.iloc[0]
        status = str(row_data["status"])

        if status != "OK":
            yield Collected(status=status, path="", data=_empty())
            return

        # No structure is produced: the result is the sequence column, not a path.
        yield Collected(status=status, path="", data=_fields(row_data))

    return one
