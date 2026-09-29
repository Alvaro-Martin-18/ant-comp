# What this project is

A **protein-design workspace** built on [`prosapia`](https://github.com/jlmoraleshellin/prosapia)
(CLI: `sapia`), running on one of **two backends**: **Modal** (a container per task) or the
**VIB DataCore** (SLURM array jobs over ssh). No design software is installed locally —
every tool runs remotely, and the data stays there too, on a Modal Volume or on cluster
group storage.

The intended way of working: **Claude drives the campaign, the backend does the compute.**
Claude never holds the data; it submits steps, waits for them, collects them, and reads the
tables.

The science is identical either way — same CLI, same tables, same lineage, same skills.
What differs is the execution machinery, and one thing that will change plans: **the six
custom tools in `tools/` (`atomium`, `bindcraft2`, `chainsel`, `cms`, `mkcomplex`,
`ringfit`) exist on Modal only.** They ship Modal images, not cluster activation scripts, so
vib registers the 11 built-ins and nothing else (verified with `sapia run --help` on both).
A chain needing `mkcomplex`/`chainsel`/`cms` must run on Modal. **A single `run_dir` lives
on one backend** — the Volume and the cluster filesystem are separate worlds.

## prosapia in one page

A **workbench, not a pipeline.** There is no fixed order of steps — there is a shared
tabular data format and tools that consume and produce it. You compose a workflow
dynamically, forking and back-tracking as the science demands.

Two principles:

1. **A design is a row; a generation of designs is a table.** Each row is keyed by a design
   `name`; each tool contributes columns (a structure path, a sequence, a pLDDT, an RMSD).
2. **When a protein diverges in sequence or structure it is a child, not the same protein**
   — so it needs a new table.

Everything happens inside a **`run_dir`**, minted by `sapia new_run` and never by a tool:

```
run_dir/
├── _registry.tsv      catalog of every table + its lineage (parent, gen)
├── table0.tsv         one row per design, keyed by `name`
├── table1.tsv         a child table (another generation)
├── .manifests/        transient per-run manifests
└── table0/            per-table, per-tool outputs on disk
    └── rfdiffusion3/
```

**Two phases per tool.** `sapia run` fans tasks out (each writes only its own files);
`sapia collect` folds those files into the table. This keeps table writes safe under heavy
parallelism and repeatable on reruns.

**Two kinds of tool.** The `action` decides where output lands, and **you never name the
output table** — the driver derives it:

- **`create`** (rfdiffusion3, proteinmpnn) mints a **child table**, `gen+1`, and links each
  new row to its parent. Produces new entities.
- **`update`** (boltz, alphafold3, usalign) annotates the **same table in place**, adding
  columns. Derives a property of designs that already exist.

Columns are leaf-prefixed per tool (`boltz_ptm`, `proteinmpnn_score`), and
`<leaf>_status == "OK"` marks a design that actually succeeded. Reruns skip already-`OK`
rows unless `--force`.

## How work gets done here: the agents

Claude Code subagents **cannot spawn further subagents**, so the setup is two layers, not
three:

```
main session = thinker                    (claude --agent thinker; Opus)
   ├─ modal-orchestrator subagent         (Sonnet; Bash/Read/Skill)   → Modal
   └─ vib-orchestrator subagent           (Sonnet; Bash/Read/Skill)   → VIB DataCore
         └─ reads .claude/skills/<tool>/SKILL.md on demand
```

- **`thinker`** (`.claude/agents/thinker.md`) — owns the scientific problem: goals, what to
  try next, reading result tables, what to keep. **Never runs `sapia`, `modal` or `ssh`
  itself.** Delegates intent ("20 backbones, length 90–110") and requires the run_dir,
  table, row count and failures back. **Asks the user which backend** at the start of a
  session if they haven't said, then uses that one orchestrator throughout.
- **`modal-orchestrator`** (`.claude/agents/modal-orchestrator.md`) — executes on Modal:
  the workstation, tool images, the `.exit` wait loop.
- **`vib-orchestrator`** (`.claude/agents/vib-orchestrator.md`) — executes on the VIB
  DataCore: ssh to the login node, `sbatch`, partitions, the `squeue`/`sacct` wait loop.
- Both orchestrators report back and **stop**; neither chains into the next tool on its own.
- **Per-tool skills** (`.claude/skills/<tool>/SKILL.md`) — flags,
  verified invocations, collected columns and the specific traps of each tool. Loaded by
  the orchestrator instead of re-reading the full docs. **They were written against
  Modal**: the tool flags, input/output columns and scientific traps hold everywhere, but
  anything about Volumes, images, `modal-shell`, `--gpu-type` or `.exit` files does not
  apply on vib.

Start a campaign with `claude --agent thinker`.

## Execution model A: Modal

`run_dir`s live **only on a Modal Volume**, mounted at `/runs`. `sapia` runs next to it in
a small container, the **workstation**, not on this machine:

```bash
sapia modal-shell                       # interactive shell; cwd is /runs
sapia modal-shell --cmd '<command>'     # one command, exits with its code
```

Non-negotiables, each learned the hard way:

- **Never run `sapia` outside the workstation.** Locally there is no `/runs`, so the run
  fails or, worse, writes paths no task container can resolve.
- **Run `sapia` from `/runs`** (the workstation's cwd) and use the relative `run_dir` that
  `new_run` printed. Never `cd` into a run_dir — stored paths become relative to it and
  break for every later tool.
- **`sapia run` is detached.** It returns once tasks are queued. The `.exit` files below are
  the only reliable completion signal.
- **A tool's first run builds its Modal image inside the submit call** — Boltz took over
  10 minutes. Allow 15+ minutes before assuming a submit is stuck. Later runs are seconds.
- **`-g 0` for CPU-only tools.** `--gpus-per-task` defaults to 1.

### Task status: the `.exit` loop

`<run_dir>/<table>/<leaf>/<script>_logs/` holds, per task, `<script>_<id>.out`,
`.err`, `.exit` (the exit code — written even on failure, `255` if the wrapper itself
failed), plus `<script>_modal.json` = `{"app_id", "n_tasks"}`.

| State | How to tell |
| --- | --- |
| done | `n_tasks` `.exit` files, all `0` |
| failed | an `.exit` that isn't `0` → read the matching `.err` |
| running | `.exit` files missing, app still running |
| killed | `.exit` files missing, app `stopped` (timeout/OOM) → `modal app logs <app_id>` |

The loop is **submit → poll `.exit` → check codes → collect**. Poll every 30–60 s (each
workstation call costs 5–10 s of cold start).

`modal` CLI commands (`app list`, `app logs`, `volume ls`) run **locally**, not in the
workstation. Set `NO_COLOR=1` before parsing their output — ANSI codes break JSON parsing.

## Execution model B: the VIB DataCore (SLURM)

No containers here. `sapia` lives in a shared env on the cluster, runs on the **login
node**, and submits SLURM array jobs; each tool's environment comes from an activation
script. Everything goes over one ssh hop, configured by four variables in **this repo's
`.env`** (read locally by the orchestrator, never by prosapia): `SAPIA_VIB_HOST`,
`SAPIA_VIB_ENV_DIR`, `SAPIA_VIB_ACTIVATE`, `SAPIA_VIB_PROJECT_DIR`.

```bash
set -a; . ./.env; set +a
timeout 120 ssh -o BatchMode=yes -o ConnectTimeout=20 -x "$SAPIA_VIB_HOST" \
  "cd $SAPIA_VIB_ENV_DIR && $SAPIA_VIB_ACTIVATE && <command>"
```

Non-negotiables, each learned the hard way:

- **Always `cd` to `$SAPIA_VIB_ENV_DIR` first.** A SLURM task starts in the submit
  directory and its prelude sources `.env` relative to it. Submit from anywhere else and
  the CLI still works, but every task dies with
  `sapia: set SAPIA_ACTIVATE_<TOOL> in your .env`.
- **Run_dirs are absolute**, minted with `--base "$SAPIA_VIB_PROJECT_DIR/outputs"`. Group
  storage is shared by login and compute nodes, so an absolute path resolves wherever a
  task lands. Still never `cd` into a run_dir.
- **GPU jobs need `--partitions`.** The default partition `gp_64C_128T_512GB` has no GPUs,
  so a GPU task submitted without one never runs. Check `sinfo -o '%P %G'` and `sinfo -s`;
  don't pick a partition unilaterally for a large batch. `--gpu-type` is Modal-only.
- **Stay under half a partition's GPUs.** Each holds 4–16 in total across its nodes;
  `--max-gpu-fraction` defaults to `0.5` and caps concurrency accordingly. Leave it there
  unless the caller has a reason — the capacity is the lab's, not ours.
- **Prefer a `_co_pi` partition** when one exists for that GPU; it is this group's
  entitlement. Plain `gpu_h100_*` / `gpu_b300_*` belong to others even when idle. a100,
  l40s, `gpu_ds` and `gpu_short` have no `_co_pi` form and are fine as they are.
- **Every `sapia run` needs `-a $SAPIA_VIB_ACCOUNT`.** prosapia omits `--account` from the
  `sbatch` line when unset, so the failure surfaces at the scheduler, not in `sapia`. If
  the variable is empty, ask the user — never submit without it or guess a value.
- **Certificates expire.** Auth is a short-lived SSH cert from a Smallstep CA via an Azure
  AD browser sign-in. If an ssh call hangs or says `Permission denied`, **stop and ask the
  user to sign in** — never retry in a loop. Always wrap ssh in `timeout`.
- **The login node and filesystem are shared with the lab.** `new_run`, `run`, `collect`
  and read-only inspection there are fine; never run a tool's compute on it. Never touch
  another user's directories, the shared `envs/`, or the cluster's `.env` / `activation/`.
- **`-g 0` for CPU-only tools**, same as Modal.

### Task status: the SLURM loop

**There are no `.exit` files.** The scheduler is the source of truth; the log dir
(`<script>_<jobid>_<taskidx>.out` / `.err`) is the evidence. Record every job ID —
`--partitions` or >1000 tasks produce several arrays.

```bash
squeue -h -j <jobid> -o '%i %T %R'
sacct -n -P -X -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

| State | Meaning |
| --- | --- |
| `PENDING` | queued; `%R` gives the reason. `ReqNodeNotAvail` / `PartitionConfig` will never start — stop and report |
| `RUNNING` | keep polling |
| `COMPLETED` `0:0` | done — collect only when **every** array task is here |
| `FAILED` | read that task's `.err` |
| `TIMEOUT`, `OUT_OF_MEMORY`, `CANCELLED` | killed; report `Elapsed`/`MaxRSS` |

The loop is **submit → poll `squeue`/`sacct` → check states → collect**. Poll every
60–120 s; queue waits run from minutes to hours, and each call is a fresh ssh connection.

A task that fails with **empty `.out` and `.err`** died before its first statement — look
at the activation script and the prelude, not the tool.

## Local configuration

`prosapia` is installed **from GitHub, branch `dev`**, pinned to a commit in `uv.lock`:

```toml
[tool.uv.sources]
prosapia = { git = "https://github.com/jlmoraleshellin/prosapia", branch = "dev" }
```

Consequences, both of which matter:

- **The install is not editable.** Edits in a local `../prosapia` checkout have **no effect**
  on what runs — neither locally nor in the workstation image, which mounts the installed
  package from `.venv/lib/python3.13/site-packages/prosapia/`. To pick up new commits on
  `dev`: `uv lock --upgrade-package prosapia && uv sync`. To develop the library, switch the
  source back to `{ path = "../prosapia", editable = true }` for the duration.
- **`docs/` is not shipped.** The wheel carries `src/prosapia/` only. The library *source* is
  therefore readable at `.venv/lib/python3.13/site-packages/prosapia/`, but the prose docs
  that several skills cite (`docs/running-on-modal.md`, `docs/lineage-and-tables.md`, …) exist
  only in a checkout of the repo. Keeping a sibling `../prosapia` clone on `dev` is still
  worthwhile for that reason alone — it is just no longer the dependency.

This pin governs **Modal only**. The cluster env is installed and updated separately by
whoever maintains it; a `uv.lock` bump here does not move it, and the two can drift.

`dev` keeps moving. Read the source before trusting a detail, and prefer
`sapia run <tool> --help` (in the workstation, or on vib) over memory.

`.env` (no secrets; Modal auth lives in `~/.modal.toml`, vib auth in a short-lived SSH
cert). Note `SAPIA_EXECUTOR` is **not** set here: the Modal workstation exports
`SAPIA_EXECUTOR=modal` itself, and the cluster's own `.env` sets `slurm`.

| Variable | Value here | Meaning |
| --- | --- | --- |
| `SAPIA_MODAL_RUNS_VOLUME` | `sapia-runs` | The runs Volume. **Required.** Created on first use. |
| `SAPIA_MODAL_VOLUME_RFD3_CKPT` | `rfd3-checkpoints` | rfd3 checkpoints, mounted at `/checkpoints`. |
| `SAPIA_MODAL_VOLUME_BOLTZ_CACHE` | `boltz-cache` | Boltz weights + CCD, mounted at `/boltz_cache`. |
| `SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE` | `bindcraft-cache` | BindCraft2's AlphaFold params (~5.3 GB). |
| `SAPIA_VIB_HOST` | `vib` | ssh alias for the DataCore login node (`~/.ssh/config`). |
| `SAPIA_VIB_ENV_DIR` | `/data/groups/csb/…/prosapia-workstation-dev` | Where prosapia is installed on the cluster. Every remote command runs from here. |
| `SAPIA_VIB_ACTIVATE` | `source .venv/bin/activate` | Puts `sapia` on PATH, after cd-ing there. |
| `SAPIA_VIB_PROJECT_DIR` | **empty — per user** | Where run_dirs go. Blank by design; each person sets their own before the first vib run. |
| `SAPIA_VIB_ACCOUNT` | **empty — per user** | SLURM account for `-a`. Blank by design; required on every vib run. |

The weight Volumes were deliberately renamed **without** a `sapia-` prefix so teammates who
don't use prosapia can share them. Note that prosapia's own defaults are still
`sapia-rfd3-checkpoints` / `sapia-boltz-cache`, and `get_named_volume` uses
`create_if_missing=True` — so a teammate missing these `.env` lines silently gets a new
**empty** Volume rather than an error. Keep these lines when copying the project.

**Never delete a weight Volume to tidy up.** `rfd3-checkpoints` is 2.5 GiB and must be
populated by hand (see below); `boltz-cache` re-downloads on first use.

## Validated pipeline

The full chain has been run end to end **on Modal**. Reference run: `outputs/20260927_211428_rfd3_denovo`
on `sapia-runs-test` (the previous test Volume; the current one is `sapia-runs`). The
timings below are Modal's; on vib, add the queue wait and drop the image-build time.
**The chain has not yet been validated end to end on vib** — treat a first run there as a
test, and start small.

| Step | Command | Result |
| --- | --- | --- |
| rfd3 de novo | `sapia run rfdiffusion3 <run_dir> --length 80-120 --num-designs 5` | 5 backbones → `table0`, ~3 min |
| ProteinMPNN | `sapia run proteinmpnn <run_dir> -t table0 -i rfdiffusion3_path --num-seq-per-target 2` | 10 sequences → `table1`, ~1.5 min |
| Boltz | `sapia run boltz <run_dir> -t table1` | 10 predictions annotated onto `table1`, ~16 min (incl. first weight download) |
| PyRosetta | `sapia run pyrosetta <run_dir> -t table1` | FastRelax (1 cycle) + ref2015 metrics onto `table1`, ~1.5 min for 2 designs (`outputs/20260928_113443_pyrosetta_test` on `sapia-runs`) |

All 10 Boltz predictions came back confident (confidence 0.91–0.97, pLDDT 0.93–0.97). Note
that high confidence is **not** proof the sequence folds to its designed backbone — the real
test is self-consistency: compare each prediction back to its parent backbone with
`usalign` (an `update` tool). That step has not been run yet.

## Known gaps

- **Tool weights are installed by hand.** There is no `sapia` command for it. rfd3's
  checkpoint was populated with a one-off `foundry install rfd3 --checkpoint-dir /checkpoints`
  from a container of the tool's image. A `setup()` hook in each `modal_image.py` plus a
  `sapia modal-setup <tool>` verb would let an agent do this itself.
- **`--length min-max` does not vary length within a batch.** rfd3 draws one length per
  batch, so `--num-designs 5` gave 5 backbones of the *same* length (103 aa). For a spread,
  use several batches (`--set n_batches=5`) or several design keys.
- **ProteinMPNN's `default_input_column` is `rfdiffusion_path`** — RFdiffusion, not
  rfdiffusion3. Coming from an rfd3 table you must pass `-i rfdiffusion3_path` or the run
  silently submits nothing.
- **Untested:** polling from inside a single long-running workstation container (unclear
  whether its Volume mount refreshes to show task commits). The orchestrator therefore uses
  repeated short `modal-shell` calls, which is what was actually verified.
- **The custom tools are Modal-only** (see the top of this file). Giving one to vib means
  writing an activation script there and a `SAPIA_ACTIVATE_<NAME>` entry in the cluster's
  `.env` — which lives in the shared env dir and is not ours to edit unilaterally.
  `tool-creator` scaffolds a `modal_image.py`, not an activation script.
- **Nothing syncs the two backends.** No command moves a `run_dir` between the Modal Volume
  and cluster storage; a campaign started on one finishes on that one.
- **Untested on vib:** the full chain, and every custom tool by definition.

## Network note (imec)

This network intercepts TLS to Modal's blob storage, so `modal volume get` / `put` can fail
with a certificate error while `modal volume ls` works. Because data stays on the Volume and
`sapia` runs in the workstation, normal work is unaffected — only local up/downloads are.
Do not work around it by weakening TLS verification; move the data through the workstation
or raise it with IT.
