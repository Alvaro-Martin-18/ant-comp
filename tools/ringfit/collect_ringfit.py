#!/usr/bin/env python3
"""
Collect ringfit results into the table.

Scans <run_dir>/<table>/ringfit/ for the per-design TSV files written by
ringfit.sh, and merges the ring-fit metrics back into the table as
<prefix>_status / <prefix>_path (the design superposed onto the reference) plus
one <prefix>_<metric> column per metric:

    align_rmsd                  CA RMSD (A) of the target-chain superposition --
                                the trust metric for every ring/lipid number
    seq_match_frac              fraction of the matched design/reference residue
                                pairs whose identity agrees: 1.0 = the pairing is
                                right, low = a register shift or the wrong chains
                                (below 0.9 the status becomes 'warn: ...', and the
                                metrics are still reported)
    resnum_offset               modal ref_resnum - design_resnum over the matched
                                pairs: 0 = the numbering agreed, 17 = the design
                                was renumbered from 1 against a 18-155 reference
    n_clash, n_clash_res        binder heavy atoms / residues inside the clash
                                cutoff of the reference's other protein chains
    min_dist_ring               closest binder-to-other-chain heavy-atom distance
    lipid_clash, min_dist_lipid same, against the reference's lipid HETATMs
    bsa_t1, bsa_t2, bsa_total   buried surface area (A^2) against target chain 1,
                                target chain 2, and both together (t1/t2 follow
                                the order of --target-chains at run time)
    bridge_ratio                min(bsa_t1, bsa_t2) / max(...): 1.0 = an even
                                straddle of the seam, ~0 = a one-protomer binder
    n_contact_res_t1/_t2        target residues contacted per protomer
    hotspot_recall, hotspot_hits fraction / list of --hotspots residues contacted
    binder_len, binder_chains   size and identity of the binder chain(s)

NA (empty) rather than 0 marks a metric that was not applicable: no --hotspots,
no --lipid-resnames, or no reference chains/lipids to measure against.

Usage:
    sapia collect ringfit outputs/20260928_100140_7ojg_binder --table table1
"""

from collections.abc import Iterable
from typing import Any

import pandas as pd
from prosapia.core import CollectCtx, CollectEach, Collected, DesignCtx

# Bare column names; the driver leaf-prefixes them (ringfit_<name>).
FLOAT_COLUMNS = [
    "align_rmsd",
    "seq_match_frac",
    "min_dist_ring",
    "min_dist_lipid",
    "bsa_t1",
    "bsa_t2",
    "bsa_total",
    "bridge_ratio",
    "hotspot_recall",
]
INT_COLUMNS = [
    "resnum_offset",
    "n_clash",
    "n_clash_res",
    "lipid_clash",
    "n_contact_res_t1",
    "n_contact_res_t2",
    "binder_len",
]
STR_COLUMNS = ["hotspot_hits", "binder_chains"]
METRIC_COLUMNS = FLOAT_COLUMNS + INT_COLUMNS + STR_COLUMNS


def _empty() -> dict[str, Any]:
    return {c: None for c in METRIC_COLUMNS}


def _present(row: pd.Series, col: str) -> Any | None:
    """The cell's value, or None when the column is absent, NaN or blank."""
    value = row.get(col)
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    return value


def _metrics(row: pd.Series) -> dict[str, Any]:
    """One result row -> collected columns, NA-safe (a missing or blank cell is
    a metric that did not apply, not a zero)."""
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


def collect_ringfit(ctx: CollectCtx) -> CollectEach:
    """Per-design ringfit collector. The framework iterates ready designs and stamps
    status/path (keyed by the tool leaf); this reads one design's one-row <name>.tsv
    and adds the ring-fit metric columns. Variants are distinguished via --dir-label,
    matching the output dir."""

    def one(d: DesignCtx) -> Iterable[Collected]:
        tsv_path = ctx.out_dir / f"{d.name}.tsv"

        if not tsv_path.is_file():
            yield Collected(status="missing", path="", data=_empty())
            return

        result_df = pd.read_csv(tsv_path, sep="\t")
        if result_df.empty:
            yield Collected(status="error: empty tsv", path="", data=_empty())
            return

        row_data = result_df.iloc[0]
        status = str(row_data["status"])

        # 'warn: ...' means the metrics were computed but the design/reference
        # residue match looks wrong (low seq_match_frac) -- keep the numbers so
        # the suspicion can be judged from the table, while the status keeps the
        # row out of any `== "OK"` selection.
        if status != "OK" and not status.startswith("warn:"):
            yield Collected(status=status, path="", data=_empty())
            return

        aligned_path = row_data.get("aligned_path")
        yield Collected(
            status=status,
            path=""
            if aligned_path is None or pd.isna(aligned_path)
            else str(aligned_path),
            data=_metrics(row_data),
        )

    return one
