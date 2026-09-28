---
name: boltz
description: How to run the boltz structure-prediction tool on Modal — flags, its weight cache Volume, timing, and the confidence columns it writes. Load before composing a boltz run.
---

# boltz

Predicts structures for designed sequences. **`action: update`** — annotates the table it
reads, in place. It does **not** create a table, and `-t` is required.

There is no dedicated doc page; `sapia run boltz --help` is authoritative.

IMPORTANT: prosapia's bundled boltz is a wrapper of the original. Check the github repo for all information: https://github.com/jwohlwend/boltz/tree/main/docs

## Verified invocation

```bash
sapia run boltz <run_dir> -t table1
```

10 sequences predicted in ~16 min on one task, including the first-time weight download.
Later runs reuse the cache and start much faster.

`default_input_column` is `proteinmpnn_sequence`, so straight after a ProteinMPNN step no
`-i` is needed. From any other source, pass `-i <sequence column>`.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--shard-size` | 10 | YAML inputs per shard directory / task. |
| `--devices` | 1 | GPUs per task; also sets `--gpus-per-task` to match. |
| `--use-msa-server` | off | Adds MSA information. Slower, and calls an external server — don't enable it without being asked. |
| `--chains` | all | Chain mini-language, e.g. `A:D`. Letters beyond the sequence's chain count are dropped. |
| `--positions` | full | Position mini-language: `/` maps groups onto chains, `start:end` inclusive 1-indexed, open ends allowed (`10:`, `:50`), `{expr}` islands resolve up the lineage. |
| `--template-yaml` | none | File spliced verbatim as a `templates:` block into every input YAML. |

Default Modal resources: **A10**, 24 CPU, 64 GiB, 8 h timeout.

## The weight cache

Weights and the CCD download on first use into the Volume named by `SAPIA_MODAL_VOLUME_BOLTZ_CACHE`, mounted at `/boltz_cache`. It holds `boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols/` and `mols.tar`. Budget extra time on the very first run of a fresh cache; never delete that Volume to "clean up".

## What it collects

Columns added to the **same** table (leaf-prefixed `boltz_`): `confidence_score`, `ptm`, `iptm`, `ligand_iptm`, `protein_iptm`, `complex_plddt`, `complex_iplddt`, `complex_pde`, `complex_ipde`, plus `_status` and `_path` (the predicted structure).

Reading them: `confidence_score` and `complex_plddt` around 0.9+ is a confident prediction for a de-novo monomer; `ptm` tracks global fold confidence. High confidence means the predictor believes the fold — it is **not** proof the sequence folds to the backbone it was designed for. For that, compare the prediction back to its parent backbone with `usalign`.

##  Important notes

- Forced templates steer; they don't constrain fully (e.g.: `force: true, threshold: 2.0`). Use the flag anyway but always verify independently and exclude rows where the target didn't land.
