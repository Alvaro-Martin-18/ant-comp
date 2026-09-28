---
name: proteinmpnn
description: How to run the proteinmpnn sequence-design tool on Modal — flags, the input-column trap when coming from rfdiffusion3, symmetry, and what it collects. Load before composing a proteinmpnn run.
---

# proteinmpnn

Designs sequences for input backbones. **`action: create`** — mints a child table, one row
per designed sequence.

Full reference: `docs/tools/proteinmpnn.md` in the prosapia repo, and
`sapia run proteinmpnn --help`.

## Verified invocation

```bash
sapia run proteinmpnn <run_dir> -t table0 -i rfdiffusion3_path --num-seq-per-target 2
```

5 backbones → 10 sequences in `table1`, one task on an L4, ~1.5 min including a first-time
image build.

## The input-column trap

`default_input_column` is **`rfdiffusion_path`** — RFdiffusion, *not* rfdiffusion3. Coming
from an rfd3 table you **must** pass `-i rfdiffusion3_path`, or the run finds no ready rows
and submits nothing. Check the parent table's columns before composing the command.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--num-seq-per-target` | — | Sequences sampled per input backbone. |
| `--sampling-temp` | — | Lower = more conservative. |
| `--designs-per-task` | 10 | Designs bin-packed per task. Designs with identical params share one batched call. Lower for more parallelism. |
| `--chains-to-design` | all | Chain mini-language: `A:C,E` → `A B C E`. |
| `--fixed-positions` | none | Keep positions fixed. `/` separates chains, `,` fragments, `start:end` inclusive 1-indexed, `{expr}` resolves up the lineage. E.g. `9:23/10,11,18:20`. |
| `--tied-positions` | none | Tie positions across chains. Same language. Mutually exclusive with `--symmetry`. |
| `--symmetry` | none | Homo-oligomer convenience: `auto` or an integer chain count. Ties all chains, and lets `--fixed-positions` describe a single asymmetric unit. **A plain number, never a point group** (`12`, not `C12`). |
| `--bias-aa` | none | `'D:1.39 E:1.39 K:1.39'`. |
| `--set` | none | Raw flag passed through, repeatable: `--set '--omit_AAs C'`. |
| `--seed`, `--batch-size` | — | |

Default Modal resources: **L4**, 8 CPU, 8 GiB, 1 h timeout. Weights ship inside the image,
so there is no cache Volume and no download step.

## What it collects

Child rows named `<parent>_f1`, `<parent>_f2`, … each linked to its parent backbone.

Columns (leaf-prefixed `proteinmpnn_`): `sequence`, `score` (lower is better),
`global_score`, `seq_recovery`, `T` (sampling temp), `sample`, `iteration`, plus `_status`
and `_path` (the FASTA).

`proteinmpnn_sequence` is what the structure predictors consume — it is Boltz's default
input column, so a Boltz run straight after needs no `-i`.
