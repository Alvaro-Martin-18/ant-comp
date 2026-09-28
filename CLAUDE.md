# What this project is

A **protein-design workspace** built on [`prosapia`](https://github.com/jlmoraleshellin/prosapia)
(CLI: `sapia`), running entirely on **Modal**. No design software is installed locally —
every tool runs in its own Modal container, and all data lives on a Modal Volume.

The intended way of working: **Claude drives the campaign, Modal does the compute.** Claude
never holds the data; it submits steps, waits for them, collects them, and reads the tables.

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
main session = thinker              (claude --agent thinker; Opus)
   └─ modal-orchestrator subagent   (Sonnet; Bash/Read/Skill)
         └─ reads .claude/skills/<tool>/SKILL.md on demand
```

- **`thinker`** (`.claude/agents/thinker.md`) — owns the scientific problem: goals, what to
  try next, reading result tables, what to keep. **Never runs `sapia` itself.** Delegates
  intent ("20 backbones, length 90–110") and requires the run_dir, table, row count and
  failures back.
- **`modal-orchestrator`** (`.claude/agents/modal-orchestrator.md`) — the only thing that
  executes. Knows the Modal mechanics and the wait loop. Reports back and stops; it does
  not chain into the next tool on its own.
- **Per-tool skills** (`.claude/skills/<tool>/SKILL.md`) — flags,
  verified invocations, collected columns and the specific traps of each tool. Loaded by
  the orchestrator instead of re-reading the full docs.

Start a campaign with `claude --agent thinker`.

## The Modal execution model

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

## Local configuration

`prosapia` is installed **editable from a sibling checkout**, so edits to the library take
effect immediately:

```toml
[tool.uv.sources]
prosapia = { path = "../prosapia", editable = true }
```

Currently `../prosapia` on branch `new-scheduler`. The Modal executor and the workstation
are new work on that branch — expect it to keep moving, and read the code before trusting
a detail.

`.env` (no secrets; Modal auth lives in `~/.modal.toml`):

| Variable | Value here | Meaning |
| --- | --- | --- |
| `SAPIA_EXECUTOR` | `modal` | Default executor for `sapia run`. |
| `SAPIA_MODAL_RUNS_VOLUME` | `sapia-runs` | The runs Volume. **Required.** Created on first use. |
| `SAPIA_MODAL_VOLUME_RFD3_CKPT` | `rfd3-checkpoints` | rfd3 checkpoints, mounted at `/checkpoints`. |
| `SAPIA_MODAL_VOLUME_BOLTZ_CACHE` | `boltz-cache` | Boltz weights + CCD, mounted at `/boltz_cache`. |

The weight Volumes were deliberately renamed **without** a `sapia-` prefix so teammates who
don't use prosapia can share them. Note that prosapia's own defaults are still
`sapia-rfd3-checkpoints` / `sapia-boltz-cache`, and `get_named_volume` uses
`create_if_missing=True` — so a teammate missing these `.env` lines silently gets a new
**empty** Volume rather than an error. Keep these lines when copying the project.

**Never delete a weight Volume to tidy up.** `rfd3-checkpoints` is 2.5 GiB and must be
populated by hand (see below); `boltz-cache` re-downloads on first use.

## Validated pipeline

The full chain has been run end to end on Modal. Reference run: `outputs/20260927_211428_rfd3_denovo`
on `sapia-runs-test` (the previous test Volume; the current one is `sapia-runs`).

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

## Network note (imec)

This network intercepts TLS to Modal's blob storage, so `modal volume get` / `put` can fail
with a certificate error while `modal volume ls` works. Because data stays on the Volume and
`sapia` runs in the workstation, normal work is unaffected — only local up/downloads are.
Do not work around it by weakening TLS verification; move the data through the workstation
or raise it with IT.
