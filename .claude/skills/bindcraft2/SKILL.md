---
name: bindcraft2
description: How to run the custom bindcraft2 tool on Modal — BindCraft2 binder-design campaigns, where one task is a whole campaign rather than one design. Covers the campaign cost model (--num-designs vs --max-trajectories), root vs child runs, hotspots with {expr}, --trajectory-only for handing backbones to atomium/proteinmpnn instead of BC2's own MPNN, --reuse-campaigns for collecting one campaign into two tables (how to compare BC2's own designs against atomium's on the same backbones), the table/dir labels that fork needs and the silent table collision when you omit them, the AlphaFold-parameter cache Volume, the three collect stages, and why its own confidence scores are not independent validation. Load before composing a bindcraft2 run or reading its columns.
---

# bindcraft2

**Custom tool** (lives in `tools/bindcraft2/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per binder design a campaign produced.

Upstream: [PacesaLab/BindCraft2](https://github.com/PacesaLab/BindCraft2), pinned to tag
`v1.0.3` in `modal_image.py`. Its `docs/source/reference.md` is authoritative for
settings, `docs/source/outputs.md` for columns.

> **Not yet run end to end.** Everything below is read off the upstream source and this
> tool's own smoke test (manifest + collector, checked locally). The first real submit
> should be one small campaign, not a fan-out — see *First run* below.

## The thing to understand first: a task is a campaign

Every other design tool here maps one design to one unit of work. BindCraft2 does not.
One `bindcraft design` process takes **one target** and loops: hallucinate a backbone
with AF2, redesign its sequence with ProteinMPNN, refold it, filter it, keep or discard —
until `--num-designs` have been **accepted** or `--max-trajectories` attempts are spent.

Consequences that shape how you drive it:

- **One manifest row = one campaign = one target.** `-C/--max-concurrent` throttles
  targets, not designs.
- **The row count is unknown until collect.** A campaign can accept fewer designs than
  you asked for, or none at all. That is a normal outcome, not a failure.
- **`--num-designs` is a stopping condition, not a batch size.** The real cost knob is
  `--max-trajectories`: a campaign runs until one of the two fires.
- **It already contains the rfd3 → mpnn → predict → filter chain.** Do not compose those
  tools around it; compose *independent* checks after it (see the last section).

## Two shapes of run

- **Root run (no `-t`)** — starts a fresh lineage in `table0`, one campaign. The target is
  `--target-pdb <file>` (design group `<stem>_bc2`) or `--shipped-target <name>` (group
  `<name>_bc2`), for a target BindCraft2 ships: `hPDL1`, `hPD1`, `mPDL1`, `hIL2R`,
  `hIL7RA`, `dynorphin_a`.
- **Child run (`-t <table>`)** — one campaign per ready row, targets from
  `-i/--input-column` (default `pdb_path` — set it to whatever column actually holds your
  target, e.g. `mkcomplex_path`, `chainsel_path`). Mints a child table.

`{expr}` placeholders resolve up the lineage, so they need `-t`. A root run must use
literal values, and says so rather than resolving to nonsense.

## Invocation

```bash
# root: one campaign against a target not in any table yet
sapia run bindcraft2 <run_dir> \
    --target-pdb targets/PDL1.pdb --chains A \
    --hotspots 'A54,A56,A66,A115' \
    --binder-lengths 60-100 --num-designs 10 --max-trajectories 200

# child: one campaign per target row, epitope read from the table
sapia run bindcraft2 <run_dir> -t table0 -i pdb_path \
    --hotspots 'A{epitope_start}-{epitope_end}' --num-designs 5

sapia collect bindcraft2 <run_dir> -t <the table the run reserved>
```

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--trajectory-only` | off | Stop at **backbones** — no MPNN, no refold, no filter, nothing accepted. The handoff to this workspace's own designers; see below. |
| `--reuse-campaigns` | none | `TABLE[:LABEL]`. Submit nothing; reserve a second table over an earlier run's campaigns, to collect a second stage of the same GPU hours. See below. |
| `--max-trajectories` | BC2's own | **The cost knob.** Attempts spent before giving up. Set it or a hard target can burn the whole timeout. With `--trajectory-only` it is the *only* budget (BC2 defaults it to 100). |
| `--num-designs` | BC2's own | Accepted designs to stop at (`number_of_final_designs`). A target, not a guarantee. |
| `--hotspots` | none | Target residues the binder should contact, BC2 syntax (`A54,A56,A66-70`). Takes `{expr}`. Omit to let BC2 pick the epitope. |
| `--coldspots` | none | Regions to avoid, same syntax. |
| `--chains` | all | Target chains to design against (`A`, `A,B`). |
| `--binder-lengths` | modality's | `80` or `60-100`. Takes `{expr}`. |
| `--modality` | `binder` | `binder`, `VHH`, `peptide`, `cyclic_peptide`, `ARP`, `scFv`, `Fab`, `large_binder`, `homo_oligomer`, `multidomain`, `induced_fit`, `fold_switch`. Comma-separated to combine. |
| `--property` | none | Repeatable preset: `humanize`, `protease_stable`, `disulfide_staple`, `forced_targeting`, `initial_guess`, `mixed_topology`, `termini_together`, `termini_accessible`, `bigbang`. Validated at submit time. |
| `--core` | none | Core profile under every preset; `benchmark` for a reproducible run. Pair with `--campaign-seed`. |
| `--extra-settings` | none | YAML/JSON merged into every campaign's settings (`objective`, `aa_bias`, `min_iptm_final`, `save_design_trajectory`, …). Values take `{expr}`. |
| `--set KEY=VALUE` | none | Verbatim `bindcraft design --set`, applied over the generated file. **No spaces** — the task script word-splits it. Use `--extra-settings` for anything richer. |
| `--no-resume` | off | See below. |

`bindcraft design --list-settings` names every setting `--extra-settings` and `--set`
accept. This tool writes the settings file itself; a key set by both a flag and
`--extra-settings` is refused at submit time rather than silently resolved.

Default Modal resources: **A100**, 8 CPU, 32 GiB, **12 h** timeout. Size the timeout to
the campaign with `-T`, and remember a killed container leaves no `.exit` file.

### `resume` is on by default

This tool writes `"resume": true` into every settings file, so re-submitting continues
the campaigns already started rather than refusing or restarting. It has to: the
framework's own resume filter reads `<leaf>_status`, which for a `create` tool lives in
the *child* table, not the one the run reads — so `sapia run` cannot skip finished
campaigns on its own. `--no-resume` forces a clean start.

## Backbones only: `--trajectory-only`

BC2's ProteinMPNN is **not** swappable by configuration — it is a JAX/Haiku
reimplementation inside the package, `jit`-compiled, fed in-memory arrays. `mpnn_model`
and `mpnn_variant` only choose among `.npz` checkpoints of that same architecture
(`neutral` / `negative` / `positive`). Nothing reaches a PyTorch model.

What BC2 *does* ship is a clean seam. With `--trajectory-only` it builds the backbone
and stops — `mpnn_model` is never constructed, no design is ever accepted — so you can
redesign the sequences with `atomium` or `proteinmpnn` and finish the campaign with this
workspace's own tools:

```bash
sapia run bindcraft2 <run_dir> -t table0 -l traj \
    --trajectory-only --max-trajectories 40 --hotspots 'A54,A56'
sapia collect bindcraft2 <run_dir> -t <the table the run reserved> -l traj
```

`sapia collect` needs no `--stage`: the run records `trajectory_only` in the out_dir
sidecar and collect follows it, the same contract that keeps `input_column` consistent
across the two phases.

### Two things to get right

- **Filter on `bindcraft2_completed`.** `!_Trajectories.csv` has one row per attempt,
  *terminated ones included*. Those are collected, not silently dropped — a terminated
  attempt that wrote no structure yields no row, but one that died late still can. BC2's
  own `terminated` column is blank when the attempt *succeeded*, which inverts prosapia's
  "blank means not applicable", so the collector adds an explicit boolean. Filter on it
  before spending sequence design on junk backbones.
- **Read the chain IDs before composing the atomium run.** The collected structure is the
  **complex**; the binder chain id comes from BC2's `binder_chain` setting and the target
  chains keep theirs. Look at the first collected PDB, then pass that chain to atomium's
  `--chains-to-design`. Do *not* `chainsel` the binder out first here — you want it
  redesigned in the target's context. (`chainsel` comes later, before `usalign`/`cms`,
  which do want the binder alone.)

Table naming and labels for this fork are in *Two tables from one campaign* below.

### What you give up

BC2's redesign call is not a plain MPNN pass. `mark_redesign_residues` in `MPNN_stage.py`
holds interface contacts fixed when `redesign_interface` is set, ties chains for
oligomers, masks fold-conditioning scaffold and inter-domain linkers, and applies the
campaign's `aa_bias`. A vanilla `atomium` run reproduces none of that unless you rebuild
it with `--fixed-positions`. That is a real difference in *what is being designed*, not
just in which model designs it — decide it deliberately rather than discovering it in the
results.

You also lose the whole filter battery, which is the point: you are replacing it with
`boltz` + `usalign` + `cms`, which are independent of the AF2 that built the backbone.

## Two tables from one campaign, and how they are named

One run reserves **one** table, and its output dir is derived from it
(`run_dir/<the table it reserved>/bindcraft2[_label]/`). So a second
`sapia collect -t <other table>` looks in a directory nothing ever wrote to and finds
nothing — even though the campaign it wants is right there, since one campaign folder
holds all three stages.

That bites the moment you want to **compare BindCraft2's own designs against atomium's
on the same backbones**. Both branches have to come out of one campaign, or they are
built on different trajectories and the comparison is confounded.

`--reuse-campaigns` closes it: a run that submits nothing, reserves a table, and records
in its sidecar that its campaigns live in an earlier run's out_dir. Collect follows the
record instead of assuming the two coincide.

```bash
# 1. one FULL campaign, reserving table1. Collect its backbones.
sapia run bindcraft2 R -t table0 -l traj --hotspots 'A54,A56' --num-designs 10
sapia collect bindcraft2 R -t table1 -l traj --stage trajectories

# 2. the SAME campaign's accepted designs, into table1_bc2. No GPU.
sapia run bindcraft2 R -t table0 -l bc2 --table-label bc2 --reuse-campaigns table1:traj
sapia collect bindcraft2 R -t table1_bc2 -l bc2 --stage ranked

# 3. the atomium branch off the backbones, reserving table2
sapia run atomium R -t table1 -i bindcraft2_traj_path --chains-to-design <binder chain>
```

Pass the **same `-t` the original run used** (here `table0`). The collector maps each
ready design to a campaign folder by name, and those folders are named after the
*original* run's design groups — reusing with `-t table1` would hunt for campaigns named
after backbones and find none. A root campaign needs no `-t` at all; the reuse run reads
the group names off disk.

### The two labels are different things, and one of them is mandatory

| Flag | Names | Omit it and… |
| --- | --- | --- |
| `-l/--dir-label` | the out_dir leaf and the column prefix (`bindcraft2_traj_path`) | both tables' columns are `bindcraft2_*`, with no way to tell a backbone from a validated complex |
| `--table-label` | the table itself (`table1_bc2` instead of `table1`) | **the second run silently reserves the first run's table** |

That second one is the trap, and it is silent. `derive_new_table` names a child
`table<gen>[_<label>]` from the parent's generation, so two unlabelled runs off `table0`
both derive `table1` — and `register_table` treats identical lineage as idempotent rather
than as a conflict. No warning; the two runs just share a table. **Give every fork a
`--table-label`.** Both labels have to match between `run` and `collect`.

Note the resulting name: a child of `table0` labelled `bc2` is `table1_bc2`, not `table2`.
`table2` is gen 2 — what the atomium run off `table1` reserves.

### The lineage is a fork, joined on `hash`

```
table0  targets
  ├ table1      bindcraft2 backbones (-l traj) ──→ table2 atomium ──→ boltz / usalign / cms
  └ table1_bc2  bindcraft2 ranked    (-l bc2)      BindCraft2's own answer
```

Each branch is an ordinary generation — `create` mints a child, `lookup` walks up, so
`{expr}` and ancestor columns still resolve from the bottom. Nothing about lineage needs
bending; there is simply one more generation than a plain campaign has.

The two bindcraft2 tables are **siblings**, not ancestor and descendant, so pair them
with BindCraft2's own join key rather than with lineage: `bindcraft2_traj_hash` on a
backbone equals `bindcraft2_bc2_hash` on every accepted design that came from it. One
backbone can yield several (`_seq0`, `_seq1`, … up to `kept_sequences`), so it is a
one-to-many join.

Compare like with like: `table1_bc2` carries BC2's own AF2 numbers, the atomium branch
carries Boltz's. The honest comparison is Boltz-on-both — predict `table1_bc2`'s
sequences too (`-i bindcraft2_bc2_sequence`) rather than reading BC2's `i_pTM` against
Boltz's.

## The AlphaFold parameters are a Volume, and they are not warm

~5.3 GB, downloaded on first use into `$XDG_CACHE_HOME/bindcraft` — mounted at
`/bindcraft_cache` from the `bindcraft-cache` Volume
(`SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE`). ProteinMPNN's weights ship inside the package, so
these are the only download.

**Warm it with one campaign before fanning out.** Submit N targets cold and N containers
each pull the same 5.3 GB. `BINDCRAFT_AF2_PARAMS` in `.env` overrides the location if the
weights already live somewhere shared.

Per the workspace convention the Volume is named without the `sapia-` prefix so it can be
shared. **Never delete it to tidy up** — it re-downloads, slowly.

## First run

The image builds inside the submit call: Ubuntu + jax cuda13 wheels + an editable install
of the repo. Allow **15+ minutes** before assuming a first submit is stuck; later runs are
seconds. Then:

1. One target, `--max-trajectories 10` or so, and watch it finish. This warms the weights
   Volume and proves the GPU path.
2. Check `campaigns/<name>/1_Trajectories/!_Trajectories.csv` appears. If jax fell back to
   the CPU the campaign still *runs*, ~100× slower, and says so only in a buried warning —
   that is the failure mode the image's `ldconfig` step exists to prevent.
3. Only then fan out.

Results are written incrementally, so you can read `3_Ranked/!_Ranked.csv` while a
campaign is still going.

## What it collects

Child rows keyed by BindCraft2's own design identity (campaign, modality, length, recipe
hash, `_seq<n>`), each linked to its parent target. **Not** by rank: `!_Ranked.csv` is
re-sorted by `i_pDAE` after every acceptance, so a positional key would not survive a
re-collect.

`--stage` picks what to collect, and defaults to **`auto`** — what the run recorded in
its sidecar, so you normally never pass it:

- **`trajectories`** (auto after `--trajectory-only`) — `1_Trajectories/!_Trajectories.csv`
  plus `1_Trajectories/<design>/<design>_trajectory.cif`, the hallucinated backbones.
  Adds `terminated`, `autotuned`, `length` and the derived `completed` boolean.
- **`refolded`** — `2_Refolded/!_Refolded.csv`, every scored candidate including rejects,
  with `outcome` and `failed_filters`. **This is the one to reach for when a campaign
  accepted nothing** and the question is why. Collecting it *as well as* another stage
  needs a `--reuse-campaigns` run, not just a second `-l` — see below.
- **`ranked`** (auto otherwise) — `3_Ranked/!_Ranked.csv`, only what the campaign accepted.

All three carry the same metric battery: BC2 scores a trajectory on the same filters it
runs at refold, so the columns are comparable across stages.

If `archive_trajectories` was on, the per-design folders are zipped and the collector says
so — run `bindcraft unarchive <campaign folder>` first.

Columns (leaf-prefixed `bindcraft2_`): `sequence`, `i_pDAE`, `i_pTM`, `i_pAE`, `pLDDT`,
`pTM`, `Unbound_Binder_pLDDT`, `Interface_Residues`, `Interface_BuriedArea`,
`Hotspot_Contact_Fraction`, `Surface_Hydrophobicity`, `Binder_Length`, `Binder_Net_Charge`,
`Binder_Free_Cysteines`, `rank`, `hash`, `trajectory`, `outcome`, `failed_filters`, plus
`cif_path`, `_status` and `_path`. `--metrics a,b,c` adds more; `--all-metrics` takes
every column (~60).

Two things about `bindcraft2_path`:

- It is the **complex** — binder *and* target — converted to PDB. `usalign` and `cms`
  want the binder alone, so `chainsel` it out first for those. Sequence redesign
  (`atomium`, `proteinmpnn`) is the exception: keep the complex and name the binder chain,
  so the target context is there.
- A design whose mmCIF will not parse still collects, with its metrics and a
  `_status` of `error: unreadable mmCIF`, pointing at the `.cif`. One bad file does not
  abort the collect.

`bindcraft2_sequence` is what a predictor consumes — not Boltz's default column, so a
Boltz run after bindcraft2 needs `-i bindcraft2_sequence`.

## Its own scores are not independent validation

`i_pDAE`, `i_pTM` and pLDDT here come from the **same AF2 that designed the binder**. They
are the objective the campaign optimised against, so a high value is partly a statement
that the optimiser succeeded, not that the binder is real. BindCraft2's internal filters
are a floor, not evidence.

The honest follow-ups are the ones AF2 did not see:

- **A different predictor** on the complex — `boltz` or `alphafold3` with
  `-i bindcraft2_sequence`. Agreement between two architectures means something; AF2
  agreeing with itself does not.
- **`usalign`** back to the predicted pose, to confirm the independent prediction puts the
  binder in the same place rather than merely folding it.
- **`pyrosetta`** for an energy that is not a confidence score at all.

Read `bindcraft2_rank` as position in *that* ranking and nothing more.
