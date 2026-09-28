---
name: usalign
description: How to run the usalign structure-comparison tool on Modal — the --col-a/--col-b pair instead of an input column, lineage resolution for self-consistency, the prefix that must match at collect, and the fact that failures exit 0. Load before composing a usalign run.
---

# usalign

Compares two structures per design with [USalign](https://github.com/pylelab/USalign) and
writes TM-score / RMSD columns. **`action: update`** — annotates the table it reads, in
place. `-t` is required.

There is no `docs/tools/usalign.md`; the source is the reference
(`src/prosapia/tools/usalign/` in `../prosapia`), and `sapia run usalign --help` is
authoritative for flags. Everything below was read from that source.

## It does not use `--input-column`

Its `default_input_column` is the literal string `"not applicable"`. The manifest builder
selects rows on **`--col-a`** instead (rows where that column is non-empty), so `-i` is
meaningless here. You must pass:

- **`--col-a`** — the column holding structure A. Required.
- **`--col-b`** *or* **`--ref`** — structure B. Exactly one; passing both, or neither,
  raises at submit time.

**`--col-b` is resolved up the lineage.** For each row it takes that row's own value if
set, otherwise the nearest ancestor's (via `lookup`). That is what makes self-consistency
possible: a child row's prediction compared against the backbone its parent holds, with no
column copied forward.

`--ref` is a single fixed structure for every comparison. It is resolved to an absolute
`/runs` path, so **the file must already be on the Volume**.

## Self-consistency: the canonical run

The real test of a design — does the predicted structure of the designed sequence fold
back onto the backbone it came from? `table1` holds the sequences and their Boltz
predictions; `table0` holds the parent backbones.

```bash
sapia run     usalign <run_dir> -t table1 --col-a boltz_path --col-b rfdiffusion3_path
sapia collect usalign <run_dir> -t table1 --col-a boltz_path --col-b rfdiffusion3_path
```

`boltz_path` lives on `table1`; `rfdiffusion3_path` does not — it is found on `table0`
through the parent link. Columns land as `usalign_boltz_vs_rfdiffusion3_*`.

## The prefix is the contract between run and collect

Every column and the results dir are keyed by a **prefix**, derived from the flags:
`<col_a minus "_path">_vs_<col_b minus "_path">` (with `--ref`, the ref file's stem), or
whatever `--output-prefix` you pass.

**Collect recomputes it from its own flags**, so it must be given the same ones — that is
why `--col-a`/`--col-b` reappear on the collect line above. If they disagree, collect
fails with `Results dir not found: .../<prefix>`. Passing `--output-prefix <name>` on both
phases is the safest form when the comparison is anything non-obvious.

This prefix is also what lets one table hold several comparisons side by side
(`boltz_vs_rfdiffusion3`, `boltz_vs_openfold3`, …) without `-l/--dir-label`.

## Flags

| Flag | Default | Note |
| --- | --- | --- |
| `--col-a` | — | **Required.** Column for structure A; also selects which rows run. |
| `--col-b` | none | Column for structure B, resolved up the lineage. Mutually exclusive with `--ref`. |
| `--ref` | none | One fixed structure for all comparisons. Must be on the Volume. |
| `--output-prefix` | `<a>_vs_<b>` | Names the results dir and every column. Must match at collect. |
| `--mm` | `1` | USalign `-mm`: 0 monomer, 1 multimer, 2 chain-to-complex, 3 circular permutation, 4 >2 structures, 5 fully non-sequential, 6 semi-non-sequential. **For single-chain designs, `--mm 0` is the matching mode** — the default is multimer. |
| `--ter` | `0` | USalign `-ter`: 0 all chains in all models, 1 all chains of the first model, 2 first chain only, 3 first chain split at TER. |

Resources: **1 CPU, 4 GiB, 30 min timeout, no GPU.** The manifest builder forces
`gpus_per_task = 0` itself, so **you do not need `-g 0`** — passing it is harmless but
redundant.

The image compiles USalign from source (`git clone` + `g++`), so the *first* submit builds
it inside the call — give that Bash call a long timeout like any other first run. After
that it is cached.

Structures are converted CIF→PDB at **manifest-build time**, in the workstation, cached
under `<run_dir>/.cif_to_pdb/`. So a large table makes the submit call itself slow; the
tasks are trivial.

## Gotchas

- **A failed comparison still exits 0.** The task script catches a missing input, a
  USalign crash and unparsable output, writes `ERROR: <reason>` into that design's TSV, and
  exits `0`. So **`n_tasks` `.exit` files all `0` does not mean the step succeeded** — for
  this tool the `.exit` loop only proves the tasks ran. The truth is the
  `usalign_<prefix>_status` column after collect. Check it before reporting success.
- **Rows can be dropped from the manifest silently.** If `--col-b` resolves to nothing on a
  row *and* its ancestors, that design is skipped. `Submitting N designs` may be well under
  the table's row count — compare the two and say so if they differ.
- **`missing` at collect means no TSV on disk** — usually one of those skipped rows, since
  collect iterates every row in the table (the "not applicable" input column matches no
  column, so no row is filtered out).
- **It writes no `usalign_status` / `usalign_path`.** Status is per-comparison
  (`usalign_<prefix>_status`), which suppresses the framework's usual leaf status stamp.
  Two consequences: don't look for `usalign_status`, and collect never skips already-done
  rows, so re-collecting re-reads every TSV. That is idempotent and safe.
- **`--force` on the run** re-submits every design, printing
  `Re-running USalign for all designs (including those with status OK).` Don't use it to
  fix a failed comparison until you have read the reason out of the status column.

## What it collects

Columns are leaf-prefixed *and* prefix-keyed: `usalign_<prefix>_<field>` (with `-l`, the
leaf becomes `usalign_<label>`).

| Column | Meaning |
| --- | --- |
| `_status` | `OK`, `ERROR: <reason>`, or `missing`. The only success signal. |
| `_sup_path` | The superposition PDB USalign wrote. |
| `_TM1`, `_TM2` | TM-score normalized by structure A and by structure B. Differ when the two have different lengths. |
| `_RMSD` | RMSD over the aligned residues, in Å. |
| `_ID1`, `_ID2`, `_IDali` | Sequence identity, normalized by A, by B, and over the alignment. |
| `_L1`, `_L2`, `_Lali` | Lengths of A, of B, and of the alignment. |

On disk: `<run_dir>/<table>/usalign/<prefix>/<name>.tsv` plus `<name>.pdb` (the
superposition). *(The collector's docstring names an older layout — trust this one.)*

Reading the result: for de-novo self-consistency the usual bar is **RMSD < 2 Å**, with
TM-score above ~0.9 saying the predicted fold matches the designed one. A design with
confident Boltz metrics but a poor RMSD back to its backbone is a failed design, not a
failed prediction — report both numbers together.
