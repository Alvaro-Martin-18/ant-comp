---
name: rfdiffusion3
description: How to run the rfdiffusion3 (rfd3) backbone-generation tool on Modal — flags, the root vs child run, what it collects, and its checkpoint Volume. Load before composing an rfdiffusion3 run.
---

# rfdiffusion3

Generates protein backbones. **`action: create`** — always mints a new table.

Full reference: `docs/tools/rfdiffusion3.md` in the prosapia repo, and
`sapia run rfdiffusion3 --help` (authoritative for flags).

IMPORTANT: prosapia's bundled rfdiffusion3 is a wrapper of the original. Check the github repo for all information: https://github.com/RosettaCommons/foundry/tree/production/models/rfd3/docs

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

## Binder design: hotspots and the two-target contig

There is no binder flag. A binder job is a **root run** that holds the target chains as
motif and appends one designed chain, with hotspots supplied through `--extra-spec`.

```bash
sapia run rfdiffusion3 <run_dir> \
    --input-pdb target/7ojg_AB.pdb \
    --contigs 'A18-155,/0,B18-155,/0,70-100' \
    --extra-spec target/hotspots.yaml \
    --num-designs 2 --set n_batches=4 --modal-gpu A100
```

Verified: 8 backbones, 3 chains each (`A` 138 / `B` 138 / **`C` = the binder**), ~4 min on
an A100 for a 276-residue motif plus a ~90-residue binder.

### `select_hotspots` — the syntax

Not documented anywhere in prosapia (`grep hotspot docs/` returns nothing), and **classic
RFdiffusion's `ppi.hotspot_res=[A30,A33]` list form does not carry over.** In rfd3 the
field is typed `Optional[InputSelection]` on `DesignInputSpecification`, and
`InputSelection.from_any` accepts only **a contig-style string, a bool, or a dict** — a
list raises `ValueError: Cannot convert <class 'list'> to InputSelection`.

```yaml
# hotspots.yaml — a whole-residue selection (ALL atoms of each residue)
select_hotspots: "A47,A123,B42,B155"     # ranges work too: "A47-49"
infer_ori_strategy: hotspots             # places the origin token 12 Å out
                                         # along the outward normal from the hotspot COM
```

```yaml
# atom-level form: the dict picks which atoms carry the annotation
select_hotspots:
  A47: TIP        # sidechain tip atoms
  B42: BKBN       # backbone
  B155: [CA, CB]  # explicit atom names
```

That string-vs-dict choice is exactly what the field's "atom-level or token-level"
docstring refers to.

**Specify few hotspots.** rfd3 was trained with hotspots present in 75% of PPI examples,
showing only a random subset of up to **20%** of the true hotspot atoms (ground truth =
target atoms within 4.5 Å of the binder). A dense patch is off-distribution. It is real
conditioning, not a no-op: the shipped checkpoint carries trained
`token_initializer.…is_atom_level_hotspot.weight` tensors.

Constraints worth knowing:

- Hotspots need `input` set in the same spec, and must lie **inside the contig's motif
  ranges** — annotations on residues outside the contig never enter the built structure.
- Chain/residue ids are the **input file's own numbering**, not renumbered.
- prosapia resolves `{expr}` in every `--extra-spec` string, so avoid literal braces.
- Setting `contig`/`length`/`input`/`symmetry`/`partial_t` in both a flag and
  `--extra-spec` is a `SpecConfigError`. The task script prevalidates, so a bad spec fails
  fast and cheaply — let it, rather than building a probe container to check.

### Contig traps for a two-chain target

- **Use `/0` for the chain break.** `elif part == "/0"` is the only literal the parser
  recognises. The colon form `A18-20:B18-20` does **not** error — it silently matches the
  prefix only and **drops chain B**, designing against half the target.
- **Do not pass `--length` with a binder contig.** It is the *total*, and
  `length_min -= num_motif_residues`, so `--length 400` against a 276-residue motif demands
  a 124-residue binder and fails validation. Let the contig's `70-100` govern.
- **The designed length is drawn once per spec build**, so `--num-designs 8` gives 8
  backbones of one length. `--set n_batches=K` gives K fresh draws
  (total = `n_batches × num_designs`); the value lands in `extra["sampled_contig"]`.
- rfd3 **renumbers every output chain from 1**, so a target numbered 18–155 comes back as
  1–138. Any downstream step that names target residues needs that offset.

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
