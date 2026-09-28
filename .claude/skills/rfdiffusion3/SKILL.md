---
name: rfdiffusion3
description: How to run the rfdiffusion3 (rfd3) backbone-generation tool on Modal — flags, the root vs child run, what it collects, and its checkpoint Volume. Load before composing an rfdiffusion3 run.
---

# rfdiffusion3

Generates protein backbones. **`action: create`** — always mints a new table.

Full reference: `docs/tools/rfdiffusion3.md` in the prosapia repo, and
`sapia run rfdiffusion3 --help` (authoritative for flags).

## Two shapes of run

- **Root run (no `-t`)** — starts a fresh lineage in `table0`. One design group, named
  `denovo_diff` (pure de novo) or `<pdb-stem>_diff` (with `--input-pdb`). Rows collect as
  `<group>_0`, `<group>_1`, …
- **Child run (`-t <table>`)** — one design group per ready row, inputs from
  `-i/--input-column`. Mints a child table.

`{expr}` placeholders in contigs/length resolve up the lineage, so they need `-t`. A root
run must use literal values.

## Verified de-novo invocation

```bash
sapia run rfdiffusion3 <run_dir> --length 80-120 --num-designs 5
```

No `-t`, no contig, no input. Produced 5 backbones in `table0` in ~3 min on an A100
(≈53 s of that was inference; the rest was container start and model load).

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--length` | none | `N` or `min-max`. One subunit's length when symmetric. |
| `--contigs` | none | rfd3 contig syntax: `A40-60` motif, bare `30` designed, `/0` chain break. |
| `--num-designs` | 1 | Designs per input key (`diffusion_batch_size`). Model loads once for all. |
| `--shard-size` | 10 | Design keys packed per task; they run **in series** in one process. |
| `--symmetry` | none | `C5`, `D4`, or `auto`. Write **one asymmetric unit's** contig — the sampler replicates it. There is no `--replicate`. |
| `--input-pdb` | none | Root-only, for diffusing a structure not yet in a table. |
| `--num-timesteps` | 200 | |
| `--step-scale` | 1.5 | Higher = less diverse, more designable. |
| `--partial-t` | none | Partial diffusion noise, ~5–15 Å. |
| `--extra-spec` | none | YAML/JSON of extra rfd3 spec fields merged into every design. |

Default Modal resources: **A10**, 8 CPU, 32 GiB, 4 h timeout. Override with `--modal-gpu`.

## Gotchas

- **`--length min-max` does not vary length within a batch.** rfd3 draws one length per
  batch, so `--num-designs 5` gives 5 backbones of the *same* length. For a spread, use
  several batches (`--num-designs 1 --set n_batches=5`) or several design keys — at the
  cost of loading the model more often.
- **Checkpoints live on a Volume** (`SAPIA_MODAL_VOLUME_RFD3_CKPT`, mounted at
  `/checkpoints`) and are **not installed automatically**. If a run fails for a missing
  checkpoint, the Volume is empty — say so rather than retrying. Populating it is a
  one-off `foundry install rfd3 --checkpoint-dir /checkpoints` from a container of the
  tool's image (~2.5 GiB, ~3 min).
- Expect harmless stderr noise: an `atomworks` warning about an unset env var, and a log
  line showing Modal's internal `/__modal/volumes/...` path. Files still land under `/runs`.

## What it collects

Columns (leaf-prefixed `rfdiffusion3_`): `iteration`, `rfd3_batch`, `rfd3_model`,
`rfd3_ca_rmsd_to_input` (NaN for de novo), per-chain length, plus `_status` and `_path`.

`rfdiffusion3_path` points at a **PDB** that collect converts from rfd3's `.cif.gz`, under
`<run_dir>/.cif_to_pdb/`. That is the column to feed downstream (e.g. ProteinMPNN's
`-i rfdiffusion3_path`).
