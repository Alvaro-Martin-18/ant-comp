---
name: atomium
description: How to run the custom atomium tool on Modal — AtomiUM, a private ProteinMPNN-like sequence designer with noise-conditioned weights. Covers --model-noise and how to pick it, the multi-temperature sampling count, the private-repo image built with Modal Secrets, and why its FASTA carries no score column (so designs cannot be ranked the way proteinmpnn's score allows). Load before composing an atomium run or interpreting its columns.
---

# atomium

**Custom tool** (lives in `tools/atomium/`, not bundled with prosapia).
**`action: create`** — mints a child table, one row per designed sequence.

AtomiUM is a private PyG re-implementation of ProteinMPNN-style sequence design, with
**noise-conditioned weights**: one checkpoint per training noise level. It ships
ProteinMPNN's own helper scripts and takes the same jsonl inputs, so this tool is
deliberately the `proteinmpnn` tool with the same two mini-languages, the same
signature-grouping/bin-packing, and the same lineage contract. **Everything you know
about composing a proteinmpnn run transfers** — read that skill for the mini-languages.

Repo: `AndreiSokolovskii/develop_atomium`, branch `pure_wo_jit`, **private**.

## Verified invocation

```bash
sapia run atomium <run_dir> -t table0 --num-seq-per-target 2
sapia collect atomium <run_dir> -t <the table the run reserved>
```

`default_input_column` is **`rfdiffusion3_path`** — unlike the bundled proteinmpnn,
whose default is the older `rfdiffusion_path` and silently submits nothing after an
rfd3 run. Coming from a different parent, pass `-i` explicitly.

## Flags that matter

| Flag | Default | Note |
| --- | --- | --- |
| `--model-noise` | `n05` | **The one flag that is not in proteinmpnn.** Which weight file to design with: `n00`…`n07` (training sigma 0.0 Å…0.7 Å), plus the `n05_v2` retrain. Validated at submit time, so a typo fails fast instead of per container. |
| `--num-seq-per-target` | 2 | Sequences per backbone **per temperature**. |
| `--sampling-temp` | `0.1` | One or more space-separated temperatures: `'0.1 0.2'`. See the count rule below. |
| `--seed` | 37 | **`0` means pick a random seed**, it does not mean seed zero. |
| `--designs-per-task` | 10 | Designs bin-packed per task; identical-param designs share one batched call. |
| `--chains-to-design`, `--fixed-positions`, `--tied-positions`, `--symmetry`, `--bias-aa`, `--set` | | Identical to proteinmpnn, including the guardrail that positions require a chain list. |

Default Modal resources: **L4**, 8 CPU, 16 GiB, 1 h timeout. Weights ship inside the
image — no cache Volume, no download step.

### Picking `--model-noise`

The suffix is the coordinate noise the checkpoint was **trained to denoise**, so it is a
statement about how much you trust the input geometry:

- **Low (`n00`–`n02`)** — trusts the backbone as given. For crystal structures or
  already-relaxed models.
- **Mid (`n03`–`n05`)** — the useful default range for **de novo backbones**, which carry
  generator-specific geometry quirks. `n03` is AtomiUM's own default; this tool defaults
  to `n05`.
- **High (`n06`–`n07`)** — tolerates rough or coarsely sampled backbones.

Noise level is a *cheap* axis to scan: run the same parent table twice with different
`--table-label`s and compare downstream self-consistency, not the sequences themselves.

```bash
sapia run atomium <run_dir> -t table0 --table-label n03 --model-noise n03
sapia run atomium <run_dir> -t table0 --table-label n05 --model-noise n05
```

## The sampling count is a product

`BATCH_COPIES = num_seq_per_target * len(temperatures)` — each temperature gets a **full**
`--num-seq-per-target` set. So `--num-seq-per-target 4 --sampling-temp '0.1 0.2 0.3'` is
**12 sequences per backbone**, not 4. Easy way to accidentally triple a run.

There is **no `--batch-size`**: `atomium.py` accepts the flag but never reads it. Do not
pass it through `--set` expecting an effect.

## The image is built from a private repo with Modal Secrets

`modal_image.py` clones the repo at **build** time using the `github-token` and
`github-username` Modal Secrets, then deletes `.git` so the credentialed remote is not
baked into the image.

**This has to happen at build time.** prosapia's Modal executor hardcodes the task
secret to the run's `.env` (`executors/modal.py`) and exposes no per-tool `secrets()`
hook, so a named Secret is simply not reachable from a running task. Cloning in the
task (as an ad-hoc Modal script would) also re-clones a private repo on every container.

Consequences:

- **You need access to those two Secrets** in the Modal workspace, or the image build
  fails. They are workspace-level, not personal.
- The checkout is **pinned to a commit** in `modal_image.py` (`ATOMIUM_COMMIT`). The
  branch is under active development, so bump that constant deliberately — a rebuild
  will not drift on its own.
- First run pays a **long image build inside the submit call** (torch 2.11 + cu128 and
  ~470 MB of weights). Allow 15+ minutes before assuming a first submit is stuck; later
  runs are seconds.

`torch_cluster` is a **required** dependency and easy to drop: nothing imports it
directly, but `model_lib.py` calls `torch_geometric.nn.radius_graph`, which is a
torch-cluster binding. Without it every task dies at model construction.

## What it collects

Child rows named `<parent>_a1`, `<parent>_a2`, … each linked to its parent backbone.
(The `_a` suffix distinguishes them from proteinmpnn's `_f` rows.)

Columns (leaf-prefixed `atomium_`): `sequence`, `sample`, `temperature`, `seq_rec`,
`iteration`, plus `_status` and `_path` (the FASTA).

`atomium_sequence` is what the structure predictors consume. Note it is **not** Boltz's
default input column (that is `proteinmpnn_sequence`), so a Boltz run after atomium
needs `-i atomium_sequence` — or `-i mkcomplex_sequence` on a binder table, where
mkcomplex still belongs in between.

### There is no score column — this changes how you rank

AtomiUM's FASTA reports only sequence recovery. **There is no `score` or `global_score`**,
so the model-likelihood ranking that `proteinmpnn_score` supports **does not exist here**.

- **`seq_rec` is not a quality metric.** It is similarity to the *input* sequence. On a de
  novo backbone whose input is effectively poly-glycine it is close to meaningless; it
  only carries information when redesigning a real sequence.
- **`sample` is not unique per row.** AtomiUM numbers samples `n % num_seq_per_target` and
  walks temperatures with `n // num_seq_per_target`, so with several temperatures the same
  `sample` id recurs once per temperature. Use the row name, or the
  (`sample`, `temperature`) pair, to identify a draw.

Rank atomium designs by what comes **after** them — predictor confidence and, above all,
self-consistency (`usalign` back to the parent backbone). That is the honest test for any
sequence designer, and here it is the only one available.
