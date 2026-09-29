# ant-comp — a protein-design workspace driven by Claude Code

This repository is not a program you run. It is a **workspace**: a set of instructions, agents and custom tools that let you sit in a terminal, describe a protein-design problem in plain language, and have Claude Code plan it, run it on [Modal](https://modal.com) or in the VIB DataCore, and read the results back to you.

No design software is installed on your laptop. Every tool — RFdiffusion3, ProteinMPNN, Boltz, PyRosetta, BindCraft2 and the custom ones in `tools/` — runs elsewhere: either in its own container on **Modal**, or as a **SLURM job on the VIB DataCore cluster**. The data stays where it ran, on a Modal Volume or on group storage. Your machine only holds this repo, the [`prosapia`](https://github.com/jlmoraleshellin/prosapia) library, and Claude Code.

You start a session with:

```bash
claude --agent thinker
```

and then talk to it about the science.

---

## 1. The idea in one page

Everything rests on [`prosapia`](https://github.com/jlmoraleshellin/prosapia) (CLI: `sapia`), a **workbench, not a pipeline**. There is no fixed order of steps — there is a shared tabular data format and tools that consume and produce it. You compose a workflow as the science demands, forking and back-tracking freely.

Two principles carry the whole model:

1. **A design is a row; a generation of designs is a table.** Each row is keyed by a design `name`, and each tool contributes columns — a structure path, a sequence, a pLDDT, an RMSD.
2. **When a protein diverges in sequence or structure it is a child, not the same protein** — so it needs a new table.

Work happens inside a `run_dir` — on the Modal Volume, or on cluster group storage, depending on where you are running:

```
run_dir/
├── _registry.tsv      catalog of every table + its lineage (parent, gen)
├── table0.tsv         one row per design, keyed by `name`
├── table1.tsv         a child table (another generation)
└── table0/            per-table, per-tool outputs on disk
    └── rfdiffusion3/
```

Every tool has an **action** that decides where its output lands, and you never name the output table — the driver derives it:

- **`create`** (rfdiffusion3, proteinmpnn, atomium, bindcraft2) mints a **child table**, generation + 1, linking each new row to its parent. It produces new entities.
- **`update`** (boltz, usalign, pyrosetta, cms, chainsel, mkcomplex, ringfit) annotates the **same table in place**, adding columns. It derives a property of designs that already exist.

Columns are leaf-prefixed per tool (`boltz_ptm`, `proteinmpnn_score`), and `<leaf>_status == "OK"` marks a design that actually succeeded. A blank is "not applicable", never zero.

The corollary you should internalise before your first session: **if a measurement produces one value per design, it belongs in the table as a column** — not in an ad-hoc analysis script. Columns can be filtered on, carry a status, survive into child tables through lineage, and are the campaign's audit trail. A script's output is a file nobody else can see, including the next session.

---

## 2. How Claude is wired here

Claude Code subagents cannot spawn further subagents, so the setup is deliberately two layers:

```
you  ──talk about the science──▶  thinker              (main session, Opus)
                                     │ delegates one step at a time
                        ┌────────────┴────────────┐
                        ▼                         ▼
              modal-orchestrator          vib-orchestrator     (subagents, Sonnet)
                        │                         │
                        ▼                         ▼
              Modal containers            SLURM array jobs
              + the runs Volume           on the VIB DataCore
```

| Piece | File | What it owns |
| --- | --- | --- |
| **`thinker`** | `.claude/agents/thinker.md` | The scientific problem: the goal, what to try next, what the numbers mean, what to keep. **Never runs `sapia` itself.** It delegates intent ("20 backbones, length 90–110") and demands the run_dir, table, row count and failures back. |
| **`modal-orchestrator`** | `.claude/agents/modal-orchestrator.md` | Executes on Modal. Knows the workstation, the tool images, and the `.exit` polling loop. |
| **`vib-orchestrator`** | `.claude/agents/vib-orchestrator.md` | Executes on the VIB DataCore over ssh. Knows the login node, `sbatch` flags, partitions, and the `squeue`/`sacct` polling loop. |
| **`tool-creator`** | `.claude/agents/tool-creator.md` | Builds a new prosapia tool when no existing one answers a measurement. |
| **Per-tool skills** | `.claude/skills/<tool>/SKILL.md` | Flags, verified invocations, collected columns, and the specific traps of each tool. Loaded on demand instead of re-reading full docs. |
| **`CLAUDE.md`** | repo root | Read automatically at every session start. The project's standing context. |

Both orchestrators report back and **stop** — neither chains into the next tool on its own. The division of labour is the point: the thinker is expensive and thinks about biology; the orchestrator is cheap and thinks about exit codes. Keeping the polling loop out of the main session is what keeps a long campaign affordable and the reasoning uncluttered.

### 2.1 Two backends, one workbench

The science is identical either way — same `sapia` CLI, same tables, same lineage, same skills. Only the execution machinery differs, and the orchestrator hides it. **The thinker will ask which backend you want at the start of a session** if you don't say; you can also state it up front ("run this on vib").

| | **Modal** | **VIB DataCore** |
| --- | --- | --- |
| Agent | `modal-orchestrator` | `vib-orchestrator` |
| How it runs | one container per task, image built from the tool's `modal_image.py` | SLURM array job; environment comes from an activation script per tool |
| Where data lives | the `sapia-runs` Volume, at `/runs` | group storage under `$SAPIA_VIB_PROJECT_DIR`, absolute paths |
| Getting a shell | `sapia modal-shell --cmd '…'` | `ssh` to the login node, from the env dir |
| Completion signal | `.exit` files per task | `squeue` / `sacct` job states — **there are no `.exit` files** |
| Waiting | no queue; first run of a tool builds its image (10–15 min) | queue waits from minutes to hours; images are already there |
| GPUs | `-g N`, plus `--gpu-type` | `-g N` **and `--partitions`** — the default partition has none |
| Cost | billed per second to the Modal workspace | free at point of use, but nodes are shared with the whole lab |
| **Tools available** | **all 17** — built-ins *and* the six customs in `tools/` | **the 11 built-ins only** |

That last row is the one that will actually change your plans. `atomium`, `bindcraft2`, `chainsel`, `cms`, `mkcomplex` and `ringfit` live in this repo and ship as Modal images, not activation scripts, so they are **not registered on vib** unless someone sets `PROSAPIA_TOOLS_DIR` there to a copy of `tools/`. In practice that means a binder campaign needing `mkcomplex` → `chainsel` → `cms` has to run on Modal, while a de-novo chain of `rfdiffusion3` → `proteinmpnn` → `boltz` → `usalign` → `pyrosetta` runs happily on either. Verified with `sapia run --help` on both.

Nothing stops you moving between them **between** campaigns, but a single `run_dir` lives on one backend — the volume and the cluster filesystem are separate worlds, and prosapia will not sync them for you.

---

## 3. Setup

### 3.1 What you need

- **Python 3.13** and [`uv`](https://docs.astral.sh/uv/).
- **Claude Code** (`npm i -g @anthropic-ai/claude-code`), version 2.1 or newer.
- A checkout of **`prosapia`** as a *sibling* directory of this repo (optional, but see §3.2).
- **At least one backend**, and you can set up both:
  - **Modal** — an account with access to the workspace holding the shared Volumes. Billed per second, so know who is paying before you launch a campaign.
  - **VIB DataCore** — an account on the cluster, ssh access to the login node, and membership of the group that owns the storage. Free at point of use; the nodes are shared with the lab.

### 3.2 Clone both repos side by side

This repo is all you strictly need. `prosapia` is pulled straight from GitHub (branch `dev`) by `uv`, pinned to an exact commit in `uv.lock`, so everyone runs the same library version:

```bash
mkdir -p ~/projects && cd ~/projects
git clone git@github.com:jlmoraleshellin/ant-comp.git
cd ant-comp
```

Cloning `prosapia` beside it is still worth doing: the installed package carries the library's *source* but not its `docs/`, and several skills cite those pages. Neither orchestrator has web access, so without a local checkout they cannot read them.

```bash
cd ~/projects
git clone git@github.com:jlmoraleshellin/prosapia.git
cd prosapia && git checkout dev
```

Target layout:

```
~/projects/
├── prosapia/     (branch: dev)     ← optional, for docs and source reading
└── ant-comp/     (branch: main)    ← you work here
```

Note that this clone is **reference material only**. Editing it changes nothing about what runs — see [§6.3](#63-updating-prosapia) for how the dependency actually moves.

If you use HTTPS rather than SSH keys, swap in `https://github.com/jlmoraleshellin/ant-comp.git` and the same for `prosapia`.

### 3.3 Install

```bash
uv sync                 # fetches prosapia[modal] from GitHub at the commit in uv.lock
```

`uv sync` needs read access to the `prosapia` repo. If it fails on authentication, that is a GitHub access problem, not a Python one — ask for access to the repo.

Everything below is per backend. Do the one you intend to use; doing both is fine.

### 3.4 Modal backend

```bash
uv run modal setup      # writes ~/.modal.toml — Modal auth lives there, never in this repo
```

`.env` holds no secrets, only Volume names, so it is committed. Keep these lines when you copy the project:

| Variable | Value here | Meaning |
| --- | --- | --- |
| `SAPIA_MODAL_RUNS_VOLUME` | `sapia-runs` | The runs Volume. **Required.** |
| `SAPIA_MODAL_VOLUME_RFD3_CKPT` | `rfd3-checkpoints` | RFdiffusion3 checkpoints, mounted at `/checkpoints`. |
| `SAPIA_MODAL_VOLUME_BOLTZ_CACHE` | `boltz-cache` | Boltz weights + CCD, mounted at `/boltz_cache`. |
| `SAPIA_MODAL_VOLUME_BINDCRAFT_CACHE` | `bindcraft-cache` | BindCraft2's AlphaFold parameters (~5.3 GB). |

The weight Volumes were renamed **without** a `sapia-` prefix so colleagues who don't use prosapia can share them. prosapia's own defaults are still `sapia-rfd3-checkpoints` / `sapia-boltz-cache`, and volume lookup uses `create_if_missing=True` — so if you drop these lines you silently get a **new empty Volume** instead of an error, and your first run fails for a missing checkpoint. This is the single most common setup mistake.

**Never delete a weight Volume to tidy up.** `rfd3-checkpoints` is 2.5 GiB and has to be populated by hand (`foundry install rfd3 --checkpoint-dir /checkpoints` from a container of the tool's image, ~3 min); `boltz-cache` re-downloads on first use, costing you ten minutes.

Check it works:

```bash
uv run sapia modal-shell --cmd 'sapia run --help'
```

(Prefix `sapia` with `uv run` unless you have activated `.venv` yourself — Claude does the same.) That opens the **workstation** — a small Modal container with the runs Volume mounted at `/runs` — and lists every registered tool, built-in and custom. If you see the tool list, you are set up. The first call takes ~10 s of cold start; that is normal and it is why the orchestrator polls at 30–60 s intervals rather than continuously.

### 3.5 VIB DataCore backend

Here `sapia` is **not** installed on your machine at all for this purpose — it lives in a shared env on the cluster, and the orchestrator drives it over ssh. Nothing about the cluster is containerised: each tool's environment comes from an activation script, and the compute happens on whatever node your SLURM job lands on.

**1. An ssh host alias.** Every command goes through one ssh hop, so put the cluster in `~/.ssh/config` under the name you will use in `.env`:

```
Host vib
    HostName <login-node-hostname>
    User <your-username>
```

Confirm it works before involving Claude — `ssh vib true` should return silently.

**2. Authentication is a browser sign-in, and it expires.** The cluster issues short-lived SSH certificates from a Smallstep CA via a Microsoft (Azure AD) login. While a certificate is valid every call just works; when it lapses, the next `ssh` opens a sign-in page in your browser and waits for you. The orchestrator is told never to loop on this — **if it reports that ssh hung or returned `Permission denied`, that is your cue to sign in**, either in the browser window it opened or by running `! ssh vib true` yourself in the Claude session.

**3. Five variables in `.env`**, read locally by the orchestrator, never by prosapia:

| Variable | Value here | Meaning |
| --- | --- | --- |
| `SAPIA_VIB_HOST` | `vib` | The ssh alias from step 1. |
| `SAPIA_VIB_ENV_DIR` | `/data/groups/csb/.../prosapia-workstation-dev` | Where prosapia is installed on the cluster. **Every remote command runs from here.** |
| `SAPIA_VIB_ACTIVATE` | `source .venv/bin/activate` | Run after cd-ing there, to put `sapia` on PATH. |
| `SAPIA_VIB_PROJECT_DIR` | **empty — set this yourself** | Where your run_dirs and outputs go. |
| `SAPIA_VIB_ACCOUNT` | **empty — set this yourself** | Your SLURM account, passed as `-a` on every run. |

**The last two ship blank on purpose — they are yours, not the lab's, and you must fill both in before your first vib run.** `SAPIA_VIB_PROJECT_DIR` is the directory on group storage where your run_dirs are minted; without it the orchestrator has nowhere to put anything.

`SAPIA_VIB_ACCOUNT` is the subtler one. **Every `sapia run` on vib has to carry `-a <account>`**, and prosapia leaves `--account` off the `sbatch` line entirely when it is unset — so forgetting it doesn't produce a helpful error, it produces a job the scheduler rejects or charges to the wrong place. The orchestrator is instructed to read this variable and to stop and ask you if it is empty, rather than submit without one. Set it to your own account; don't commit your value, since the next person's differs.

The env dir requirement is not arbitrary. A SLURM task starts in the directory its job was submitted from, and the prelude sources `.env` relative to that directory — so submitting from anywhere else gives you tasks that die with `sapia: set SAPIA_ACTIVATE_<TOOL> in your .env`. The cluster's own `.env` and `activation/` belong to whoever maintains that shared env; **do not edit them**.

**4. Check it works:**

```bash
set -a; . ./.env; set +a
ssh "$SAPIA_VIB_HOST" "cd $SAPIA_VIB_ENV_DIR && $SAPIA_VIB_ACTIVATE && sapia run --help"
```

You should get the tool list — **11 built-ins, and none of this repo's six custom tools** (see §2.1). If a tool you need is missing, say so rather than improvising; adding it means giving it an activation script on the cluster, not a Modal image.

---

## 4. Your first session

```bash
cd ant-comp
claude --agent thinker
```

The thinker reads `CLAUDE.md`, loads the `all-tools` skill (the catalog of what exists), and starts by **asking you about the design decisions it foresees** — including **which backend to run on**, if you haven't said. Expect questions, not immediate action; that is the agent working correctly. Answer them, because they are the decisions that determine whether the campaign means anything.

You can settle the backend in your opening line — "…and run it on vib" — and it will pick the matching orchestrator. If your chain needs a custom tool, it will (or should) tell you that Modal is the only option.

A first conversation might run like this:

> **You:** I want a de-novo binder against the periplasmic domain of PDB 7OJG.
>
> **Thinker:** *(inspects the assembly, reports that 7OJG is an 11-mer ring, identifies the membrane belt, proposes three candidate epitopes with approach vectors, and asks you to choose.)*
>
> **You:** Go with epitope 2. Twenty backbones to start.
>
> **Thinker:** *(delegates to the orchestrator; comes back with the run_dir, the table, 20 rows, 0 failures, and a proposal for the next step.)*

### What a good instruction looks like

Talk to the thinker in **intent and constraints**, not commands. It translates; the orchestrator executes. Compare:

| Say this | Not this |
| --- | --- |
| "20 backbones, length 90–110, no symmetry" | "run `sapia run rfdiffusion3 …`" |
| "4 sequences per backbone at sampling temp 0.2" | "use `--num-seq-per-target 4`" |
| "predict a random subset of 10 first, so we don't burn GPU on a bad batch" | "run boltz on everything" |
| "which of these actually fold back to their own backbone?" | "compute the RMSDs" |

If you *do* know the flag you want, say it — the thinker will pass it through. The point is that you are never required to.

### What to ask for back

The thinker will hand you numbers. Push on them. Useful questions: *what is the spread, not just the best?* *how many rows are `OK`?* *what did the failures fail on?* *does this number survive the geometry gate?* The agent is instructed to say plainly when a batch is bad, and to distrust large interfaces before the target RMSD has been checked — hold it to that.

### Ending a campaign

Ask the thinker to **write the campaign report** into `campaigns/` before you close the session. Skills accumulate tool knowledge; nothing accumulates scientific knowledge unless it is written down. See `campaigns/20260928_7ojg_slyb_binder.md` for the shape: the target's real geometry, each decision and why, the results, the traps found with their signatures, and the gaps left open. The next session — yours or a colleague's — starts from that file instead of re-deriving it.

---

## 5. What is in the repo

```
CLAUDE.md                  standing project context, read at every session start
.claude/agents/            thinker, modal-orchestrator, vib-orchestrator, tool-creator
.claude/skills/            one SKILL.md per tool + the cross-cutting ones
tools/                     custom prosapia tools (auto-discovered from ./tools)
campaigns/                 campaign reports — the scientific record
utils/                     one-off helpers (boltz_template.py)
.env                       Volume names and cluster paths, no secrets
```

### The custom tools

Beyond prosapia's built-ins, this workspace ships the following — **all of them Modal-only**, since they carry `modal_image.py` rather than cluster activation scripts:

| Tool | Answers |
| --- | --- |
| `cms` | How much real interface is there, and does it fit? Contact molecular surface + shape complementarity, on GPU. |
| `chainsel` | Pull the binder chain out of a complex; merge split protomers back into one chain. |
| `mkcomplex` | Put the fixed target chains back around a designed binder sequence, so the predictor folds the complex and not the binder alone. |
| `ringfit` | Does the binder straddle two protomers of an oligomeric target, and would it clash with the rest of the assembly or its lipid belt? |
| `bindcraft2` | A whole binder campaign as one step — AF2 hallucination + MPNN + refold + filter, looping until enough designs are accepted. |
| `atomium` | ProteinMPNN-like sequence design with a private noise-conditioned model. **Requires access to two Modal Secrets** (`github-token`, `github-username`) to build its image; without them the build fails. |

### The cross-cutting skills

- **`all-tools`** — the catalog. Which tool answers which question, what it returns, and what is registered but not actually runnable. The thinker loads this first, every time.
- **`binder-campaign`** — how to sequence and judge a de-novo binder campaign: the order the gates must be applied in, and what each metric is blind to.
- **`prosapia`** — the CLI contract, tables, lineage, labels, and the traps that make a run silently do nothing.
- **`authoring-a-tool` / `editing-a-tool`** — for extending the workbench.

---

## 6. Working with the repo day to day

Nothing scientific is stored in git here. The run data — structures, tables, logs — lives on the Modal Volume, and the repo holds only text: agent definitions, skills, tool source, campaign reports. That makes collaboration cheap, but it also means **we are all editing the same handful of files**, so a little git discipline saves a lot of annoyance.

### 6.1 Start every session with a pull

```bash
cd ~/projects/ant-comp
git pull
uv sync                        # only if pyproject.toml or uv.lock changed
```

That is the whole routine. `uv sync` reinstalls `prosapia` only when someone has committed a new pin in `uv.lock` — it does not chase the `dev` branch on its own, which is deliberate: your environment does not change under you mid-campaign.

### 6.2 Pulling when you have local edits

This is the common case: you tweaked a `SKILL.md` mid-campaign, and now you want everyone else's changes too. Git refuses to pull over dirty files, so park your work, pull, and put it back:

```bash
git stash              # park your uncommitted edits
git pull               # fast-forward to the shared state
git stash pop          # replay your edits on top
```

`git stash pop` applies the stash **and deletes it**. If you would rather keep a safety copy until you are sure the result is right, use `git stash apply` instead and drop it yourself afterwards:

```bash
git stash apply        # replay, but keep the stash
# ...check everything looks right...
git stash drop         # now discard it
```

Two things worth knowing before you need them:

- **If `pop` hits a conflict**, the stash is *not* deleted — git leaves it in place and marks the conflicting files. Resolve the markers, `git add` them, then `git stash drop` manually.
- **`git stash list`** shows what you have parked, and **`git stash show -p stash@{0}`** prints the diff. A forgotten stash is the easiest way to lose an afternoon of skill edits, so check the list if something you wrote seems to have vanished.

### 6.3 Updating prosapia

The library is pinned to a commit in `uv.lock`, so it moves only when someone deliberately moves it. To take the current tip of `dev`:

```bash
uv lock --upgrade-package prosapia
uv sync
```

Then **test before committing the new pin** — `dev` is an active branch and can break things under you:

```bash
uv run sapia run --help                              # tools still register?
uv run sapia modal-shell --cmd 'sapia run --help'    # workstation image still builds?
```

(Only the Modal side is pinned by this repo's lock. The cluster install is a separate env that someone else updates — a `uv.lock` bump here does not move it, and the two can drift.)

If both pass, commit `uv.lock` so everyone lands on the same commit. If something broke, `git checkout uv.lock && uv sync` puts you back.

**To work on the library itself**, point the dependency at a local checkout for the duration:

```toml
[tool.uv.sources]
prosapia = { path = "../prosapia", editable = true }
```

Run `uv sync` after the edit, and your changes take effect immediately, locally *and* in the workstation image. Switch the line back to the `git`/`branch` form before committing — do not push an editable path, since it resolves to a directory nobody else has.

### 6.4 Sharing your changes back

Improvements to skills, tools and campaign reports are the whole point of the shared repo — a trap you wrote down is a trap nobody else hits. Work on a branch and open a PR:

```bash
git checkout -b skills/boltz-msa-note
git add .claude/skills/boltz/SKILL.md
git commit -m "skills(boltz): note the empty-MSA trap on natural targets"
git push -u origin skills/boltz-msa-note
```

Keep out of commits: anything under a run_dir, loose `.cif` / `.pdb` / `.json` dumps pulled off the Volume, and `.venv`. If you changed `.env`, only commit it when a Volume name genuinely changed for everyone — not when you pointed at a scratch Volume of your own.

---

## 7. Things that will bite you

These were all learned the hard way. They are in `CLAUDE.md` and the agent files too, so Claude knows them — this section is so that *you* do.

**Whichever backend you use**

- **Never run `sapia` on your own machine.** Neither `/runs` nor the cluster filesystem exists locally, so the run either fails or, worse, writes paths no task can resolve. It goes through `sapia modal-shell`, or over ssh. (The `modal` CLI itself — `app list`, `app logs`, `volume ls` — is the exception; that is local.)
- **Never `cd` into a run_dir.** Stored paths become relative to it and break for every later tool.
- **`sapia run` is detached.** It returns as soon as tasks are queued; queued is not finished. The orchestrator waits and checks states before collecting.
- **`No designs to submit.` exits 0.** A run that queues nothing looks like a success. It is almost always a wrong input column — ProteinMPNN's default is `rfdiffusion_path`, and coming from an rfd3 table you must pass `-i rfdiffusion3_path`.
- **`Submitting N designs` often isn't N designs.** Several tools bin-pack, so N counts manifest rows. The honest number is the count of staged input files.
- **`--length min-max` does not vary length within a batch.** rfd3 draws one length per batch, so `--num-designs 5` gives five backbones of the *same* length. For a spread, use several batches.
- **High confidence is not proof.** A confident Boltz prediction does not mean the sequence folds to its designed backbone. The real test is self-consistency: predict, then compare back to the parent backbone with `usalign`. For binders, fold and pose are separate questions and a design can pass one at 1 Å while failing the other by 90 Å.
- **The broken predictions produce the best-looking numbers.** A huge interface is evidence of a collapsed target until the target RMSD says otherwise. Gate on geometry before reading any interface number.

**Modal only**

- **A tool's first run builds its image inside the submit call.** Boltz took over 10 minutes. Allow 15+ minutes before assuming a submit is stuck. Later runs are seconds.
- **`.exit` files are the completion signal**, under `<run_dir>/<table>/<leaf>/<script>_logs/` — written even on failure, `255` if the wrapper itself died.
- Run_dirs are **relative**, resolved against `/runs`.

**VIB DataCore only**

- **GPU jobs need `--partitions`.** The default partition `gp_64C_128T_512GB` has no GPUs at all, so a GPU task submitted without one simply never runs. Ask the orchestrator to report `sinfo` load and choose a partition deliberately for anything large.
- **Don't take more than half a partition's GPUs.** Each GPU partition holds 4–16 in total across all its nodes, and lab etiquette is to keep your concurrent tasks under 50% of that. `--max-gpu-fraction` defaults to `0.5` and enforces it for you — raising it is a decision to take other people's capacity, so have a reason.
- **Prefer the `_co_pi` partitions** where one exists for the GPU you want; those are this group's entitlement. The plain `gpu_h100_*` / `gpu_b300_*` partitions are somebody else's even when `sinfo` shows them idle. Where there is no `_co_pi` equivalent (a100, l40s, `gpu_ds`, `gpu_short`), the ordinary partition is fine.
- **Every run needs `-a <account>`** — see `SAPIA_VIB_ACCOUNT` in §3.5. A missing account fails at the scheduler, not in `sapia`.
- **There are no `.exit` files.** SLURM state is the truth: `squeue` while it runs, `sacct` afterwards. `COMPLETED 0:0` on *every* array task is the bar.
- **Run_dirs are absolute**, minted under `$SAPIA_VIB_PROJECT_DIR` with `--base`. That is correct, not a workaround: group storage is shared by login and compute nodes, so an absolute path resolves wherever the task lands.
- **Your ssh certificate expires.** When it does, the next call opens a browser sign-in and blocks. The orchestrator will stop and tell you rather than retry — sign in, then let it continue.
- **The login node is shared with the lab.** `new_run`, `run`, `collect` and read-only inspection there are fine; actual compute never is. Likewise, a big GPU array occupies nodes everyone else needs — flag the size before submitting, not after.
- **A task that fails with empty `.out` and `.err`** died before its first statement. Look at the activation script and the prelude, not the tool.
- **Exit code 0 does not mean it worked.** Some task scripts (`usalign`, `pyrosetta`) catch their own errors and still exit 0. `<leaf>_status` after collect is the proof.

### Network note (imec)

The imec network intercepts TLS to Modal's blob storage, so `modal volume get` / `put` can fail with a certificate error while `modal volume ls` works fine. Because data stays on the Volume and `sapia` runs in the workstation, normal work is unaffected — only local up- and downloads are. Do not work around it by weakening TLS verification; move data through the workstation, or raise it with IT.

---

## 8. Extending the workspace

**To add a measurement**, ask the thinker for it in scientific terms. If no tool covers it, it will route to `tool-creator`, which checks first whether an existing collector already writes those columns (building a duplicate is worse than building nothing), then scaffolds the five files a prosapia tool needs: `spec.py`, `run_<name>.py`, `collect_<name>.py`, the `.sh` task script, and `modal_image.py`. Tools in `./tools` are discovered automatically.

A tool built that way is **Modal-only**, because `modal_image.py` is what gives it an environment. Making it run on the cluster as well is a separate job: an activation script on vib that puts the binary on PATH, and a `SAPIA_ACTIVATE_<NAME>` entry in the cluster's `.env` — which lives in the shared env dir and is not yours to edit unilaterally. Worth agreeing with whoever maintains that env before promising a tool on both backends.

**To teach the workspace something you learned**, edit the relevant `SKILL.md` — that is what they are for. A trap you hit once and wrote down is a trap the next session skips. Measured facts with their numbers are worth far more than general advice; the existing skills are written that way on purpose.

**To change how the agents behave**, keep the division of labour intact whatever you do: the thinker must not execute, and the orchestrator must not decide. That constraint is what keeps a long campaign readable and affordable.

**If you want your own variant, copy the agent rather than editing the shared one.** `thinker.md`, `modal-orchestrator.md` and `vib-orchestrator.md` are shared files that everybody pulls. Rewriting them to suit your project means a merge conflict for your colleagues on every pull, and it quietly imposes your preferences on their campaigns. Give your version its own name instead:

```bash
cp .claude/agents/thinker.md .claude/agents/thinker-membrane.md
```

Then **edit the `name:` field in the frontmatter to match the filename** — that field, not the filename, is what `--agent` resolves:

```yaml
---
name: thinker-membrane
description: Protein-design lead for membrane-protein binder campaigns.
model: opus
---
```

and start your sessions with it:

```bash
claude --agent thinker-membrane
```

A new file merges cleanly forever, and you can still pull improvements to the original. Pick a name that says what makes it different (`thinker-membrane`, `thinker-nanobody`) rather than `thinker2` or your initials — a colleague reading `.claude/agents/` should be able to tell which one they want.

The same applies to skills: a variant belongs in a new `.claude/skills/<name>/SKILL.md`, not in an edit to the shared one. **Corrections and newly-discovered traps are the exception** — those go straight into the shared skill, because everyone needs them.

---

## 9. A note on cost and honesty

Two habits make the difference between a campaign that produces a result and one that produces a bill — or, on the cluster, a queue full of someone else's stalled jobs. Modal charges you by the second; the DataCore charges you in goodwill, which is harder to refund.

**Test on a subset first.** Predicting ten designs to see whether a batch is worth predicting is nearly always the right call. The thinker knows to do this; ask for it anyway.

**Make the agent show you the failures.** Ask for row counts and `<leaf>_status` breakdowns, not just the top hits. A table where 40 of 100 rows are blank is a different result from one where all 100 succeeded and 60 scored poorly, and only one of those two is a statement about your designs.
