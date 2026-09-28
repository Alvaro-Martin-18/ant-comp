---
name: chainsel
description: How to run the custom chainsel tool on Modal — extracting a named subset of chains (typically the binder chain) from a structure column into its own file, and merging split-protomer chains back into one chain per protomer. Covers --chains, --rename-to, --merge-groups, the numbering contract (--renumber / --renumber-from), the absent-chain contract, the dir-label rule, and the columns it collects. Load before composing a chainsel run, before any binder self-consistency comparison, or before running ringfit on a prediction whose target protomers came back split across two chains.
---

# chainsel

**Custom tool** (lives in `tools/chainsel/`, not bundled with prosapia).
**`action: update`** — a sub-structure is a *property* of a design that already exists,
so it annotates the table it reads, in place, and mints no child. `-t` is required.
**CPU-only**; the manifest builder forces `gpus_per_task = 0`, so you do not pass `-g 0`.

Extracts the chains you name from each design's structure into a new file — optionally
**merging** several source chains into one output chain — and collects its path as
`chainsel_path`, an ordinary structure column that `usalign`, `pyrosetta`, `ringfit` or
a viewer then consumes like any other.

## Why it exists

Every structure comparison in this workbench is **whole-file**, and a binder complex is
mostly target:

- **`usalign --mm 1`** (the default) aligns the whole multi-chain file. On a
  target+binder complex the score is dominated by the large, fixed target and reads
  ~0.95–1.0 **whatever the binder did** — including when the binder is a completely
  different fold.
- **`usalign --ter 2`** isolates only the **first** chain. In every complex we build
  (`mkcomplex` puts the fixed target chains first) the binder is the **last** chain, so
  `--ter 2` compares target to target.

So before chainsel there was **no way to compare a binder to a binder**. That matters
because a binder campaign has two separate questions, and one number cannot answer both:

| Question | What to compare | How |
| --- | --- | --- |
| **Did it fold?** | binder chain alone, prediction vs. design | `chainsel` both, then `usalign --mm 0` |
| **Did it stay?** | the pose in the target's frame | whole complex: `usalign --mm 1`, or `ringfit` |

A binder can fold perfectly and sit in the wrong place, or hold its pose in a prediction
that folded into a different topology. Measure them separately.

### The second reason: split protomers

A target trimmed of a membrane belt is **two chains per protomer**. SlyB here is
residues 19–59 (41 aa) and 107–155 (49 aa), with a real 13.4 Å gap between the
segments, so a predictor is handed two entities per protomer and returns a complex
where protomer 1 is chains `A` + `C` and protomer 2 is chains `B` + `D`.

`ringfit` takes **exactly two** `--target-chains` — one per protomer — and there is no
syntax for "t1 is the union of A and C". On a split-protomer file the per-protomer
metrics `bsa_t1`, `bsa_t2` and `bridge_ratio` (does the binder straddle the seam?)
cannot be computed *at all*. `--merge-groups "A+C:A,B+D:B"` writes the same coordinates
back as one chain per protomer, and `ringfit --target-chains A,B` works again.

## Invocation

```bash
sapia run chainsel <run_dir> \
    --table table1 \
    --input-column boltz_path \
    --chains E --rename-to A \
    --dir-label pred_binder
sapia collect chainsel <run_dir> --table table1 --dir-label pred_binder
```

| Flag | Default | Notes |
| --- | --- | --- |
| `-i/--input-column` | **effectively required** | `default_input_column` is the sentinel `"not applicable"` (the `usalign` convention). No column is a defensible default — the tool is equally at home on `rfdiffusion3_path`, `boltz_path`, `alphafold3_path` — so the builder **raises** at submit time unless you pass a column that the table actually has. |
| `--chains` | `""` | Chain IDs to KEEP, comma-joined — **one output chain each**. **The output order follows the order given**, so `--chains E,A` writes E first. **Required unless `--merge-groups` is given.** |
| `--rename-to` | `""` | New IDs, same count and order as `--chains`. Empty keeps the originals. Use it to normalise a binder to one fixed letter across tables whose lettering differs. Two kept chains may **not** be renamed onto the same ID — that guard stays a hard error, because there it would be an accidental merge. |
| `--merge-groups` | `""` | **Merge several source chains into one output chain.** `"<chain>[+<chain>...]:<out_id>"`, groups comma-joined: `"A+C:A,B+D:B"`. Replaces `--chains`/`--rename-to` (it states both). The only way to get two source chains into one output chain. |
| `--out-format` | `pdb` | `pdb` or `cif`. `pdb` is what USalign and Rosetta want. With `pdb`, an output chain ID longer than one character is refused at submit time (PDB has a single chain-ID column). |
| `--keep-het` | off | Default **strips** HETATM — ligands, ions, waters. `--keep-het` keeps them. Modified polymer residues (MSE and friends) are kept either way: the test is "tabulated amino/nucleic acid", not the HETATM flag. |
| `--renumber` | off | Renumber each **output** chain 1..N in written order, insertion codes dropped — **continuous across a merge seam** (`A+C:A` with 41+49 residues comes out 1..90). Off preserves the source numbering. |
| `--renumber-from` | `""` | Restore auth numbering **per source chain**: `"<chain>:<first>[-<last>]"`, comma-joined — `"A:19-59,C:107-155"`. Mutually exclusive with `--renumber`. |
| `-l/--dir-label` | `""` | **In practice required** — see the traps. |

Reads PDB or mmCIF transparently (gemmi). Alternative conformations and hydrogens are
dropped on read; only the **first model** is used.

## Merging chains (`--merge-groups`)

```
--merge-groups "A+C:A,B+D:B"
                └┬┘ └┬┘
   source chains ─┘   └─ the output chain ID
```

Merging is **never** implicit. `--rename-to` still refuses to map two kept chains onto
one ID; only `--merge-groups` merges, and it says so in the flag.

**Residue order inside a merged chain is the order the group writes its sources** —
all of A's residues, then all of C's — and is **never re-sorted by residue number**. A
predictor numbers each chain 1..N independently, so sorting by number would interleave
the two segments catastrophically.

### The numbering contract — read this before every merge

Exactly one of three regimes applies, and the output always carries one of them:

| Flags | What comes out |
| --- | --- |
| neither (default) | Every residue keeps its **source** number and insertion code. A merged chain whose segments both start at 1 therefore **collides** → that design is an **error**, not a file with repeated numbers. |
| `--renumber` | Each **output** chain is numbered **1..N in written order**, insertion codes dropped. `A+C:A` with a 41- and a 49-residue segment → **1..90 continuous across the seam** (not 1..41 then 1..49). |
| `--renumber-from "A:19-59,C:107-155"` | Each **source** chain's segment is numbered **consecutively from its own first number**, insertion codes dropped. Output chain A carries **19..59 then 107..155** — the target's real auth numbering, which is what `ringfit --resnum-match resnum` pairs against the reference. |

One sentence, for the record: **`--merge-groups` writes each group's sources in the
order given, and the residue numbers in the output are the source numbers unchanged
(default, an error if they collide inside a merged chain), or 1..N per output chain
with `--renumber`, or `<first>`, `<first>+1`, … per source chain with
`--renumber-from`.**

Two things `--renumber-from` does *not* do: it does not reproduce numbering **gaps
inside** a segment (it writes `first`, `first+1`, … in file order), and it does not
guess — every kept source chain must be listed. Give the `-<last>` whenever you know
it: the residue count actually written is checked against the range, so a segment of
the wrong length is an error rather than a silently shifted numbering.

### Guards, all at submit time (the run refuses; nothing is queued)

| Combination | Result |
| --- | --- |
| neither `--chains` nor `--merge-groups` | error — chainsel needs a selection |
| `--chains` **and** `--merge-groups` | error — they are alternatives |
| `--rename-to` **and** `--merge-groups` | error — a group's output ID is the part after `:` |
| `--renumber` **and** `--renumber-from` | error — both set the numbering, pick one |
| a source chain in two groups, or twice in one | error — each source chain is written once |
| two groups with the same output ID | error — write them as one group (`A+C:A`) |
| `--rename-to` mapping two chains onto one ID | error (unchanged) — say `--merge-groups` if you mean it |
| a malformed group (`A+C`, `A+C:`, `A:B:C`) | error naming the entry |
| `--renumber-from` naming a chain this run does not keep | error naming it |
| `--renumber-from` missing a kept source chain | error — a half-renumbered file mixes two regimes |
| `--renumber-from` with `end < start` or a non-integer | error naming the entry |
| `--out-format pdb` with an output ID longer than 1 char | error (unchanged) |

Per design, recorded as data in `chainsel_status` (not a crash, no file written):

- **a chain named in a group but absent** → `error: chain C not in <file> (have: A,B)`;
- a segment length that contradicts its `--renumber-from` range →
  `error: chain C of <file> has 48 residues but --renumber-from says 107-155 (49)`;
- a merged chain whose residue numbers collide →
  `error: merged chain A of <file> (A+C) has duplicate residue numbers (1,2,...)`.

## Traps

**An absent chain is an error, never a shorter file** — including a chain named inside
a merge group. This is the whole point of the
tool's error discipline. If `--chains E` and the structure has no chain E, that design
comes back `chainsel_status = "error: chain E not in <file> (have: A,B,C)"`, with an
empty path and NA counts, and **no file is written**. A silently partial extraction
would be a perfectly well-formed PDB that every downstream TM-score, RMSD, BSA and
energy would then be computed on — wrong, plausible, undetectable.

**But a chain that is present and WRONG is not caught.** Extracting chain `A` when the
binder is chain `E` succeeds: you get a clean file of the target. Nothing errors.
Guard it by reading `chainsel_n_res` against the binder length you designed, and
`chainsel_chains` against what you asked for. If `n_res` is the target's length, you
extracted the target.

**Chain lettering differs between tools.** rfdiffusion3 and a structure predictor need
not agree on which letter the binder gets: the predictor assigns letters in the order of
the `/`-separated entities `mkcomplex` built, so the binder is the **last** letter there,
while the designed backbone may call it something else. Look at an actual file before
choosing `--chains` — do not assume. `--rename-to` then normalises the two extractions
to the same letter.

**`-l/--dir-label` is effectively required.** The tool leaf (`chainsel`) names both the
output dir and every column, so a second chainsel run on the same table with a different
chain set would overwrite the first's `chainsel_path`. Label each extraction
(`-l design_binder`, `-l pred_binder`) and the columns become
`chainsel_design_binder_*` / `chainsel_pred_binder_*`. **The same `-l` must be passed to
`collect`** — it resolves the same output dir.

**Failures still exit 0.** Errors are recorded as data, so `n_tasks` `.exit` files all
`0` does not mean the step succeeded. The truth is the `chainsel_status` column after
collect.

**Rows can be dropped from the manifest silently.** A design whose input cell points at
a file that does not exist is printed as `MISSING ... (skipping)` and left out.
`Submitting N designs` under the table's row count means exactly that — compare them.

**`--renumber` changes what "residue 12" means.** USalign aligns by sequence, so it does
not care; anything that pairs residues **by number** across two files does
(`ringfit --resnum-match resnum`, a hotspot list, a mutation table). Renumber only when
the consumer needs matching numbering, and never when you will look residues up by their
original numbers. On a merged chain `--renumber` runs **across the seam**, so the second
segment's residue 1 becomes residue 42 — use `--renumber-from` instead whenever the
downstream tool matches by number against a reference.

**A merge lowers `n_chains` but not `n_res`.** `--merge-groups "A+C:A,B+D:B"` on a
4-chain prediction collects `chainsel_n_chains = 2` and `chainsel_chains = "A,B"`,
while `chainsel_n_res` still counts every residue written — 41 + 49 = **90 on chain
A**. If `n_res` dropped, a segment went missing; if `n_chains` did not drop, the merge
did not happen.

## Columns collected (`chainsel_` prefix, or `chainsel_<label>_` with `-l`)

| Column | Meaning |
| --- | --- |
| `_path` | **The extracted structure.** The point of the tool — feed it to the next tool as `--col-a` / `--col-b` / `-i`. |
| `_chains` | The chain IDs **actually written**, comma-joined, in output order (AFTER `--rename-to` or a group's `:<out_id>`). This is what you pass to `ringfit --target-chains`. |
| `_n_chains` | How many **output** chains were written = the number of `--chains` IDs, or the number of **groups**. A merge makes it smaller than the source chain count. |
| `_n_res` | Residues written, summed over the output chains. A merge does **not** change it — a merged chain contributes all its segments (41 + 49 = 90). Your sanity check against the expected length. |
| `_n_atoms` | Heavy atoms written, summed the same way (hydrogens and alt-confs are already gone). |
| `_status` | `OK`, `error: <reason>`, or `missing`. The only success signal. |

**NA (empty) means not applicable, not zero** — the extraction failed, or no TSV was
written. A blank `_n_chains` is not "zero chains", and a blank `_path` is not "nothing
to extract".

On disk: `<run_dir>/<table>/chainsel[_<label>]/<name>.{pdb,cif}` plus `<name>.tsv`.

## Worked example: did the binder fold?

Setup: `table0` holds rfd3 designed complexes (`rfdiffusion3_path`, target chains `A,B`
+ binder chain `C`); `table1` holds the ProteinMPNN sequences, their `mkcomplex`
complex sequence, and the Boltz prediction of that complex (`boltz_path`, chains `A,B`
target + `C` binder — **check yours**).

```bash
# 1. binder chain out of the DESIGNED complex — on table0, where rfdiffusion3_path lives
sapia run     chainsel <run_dir> -t table0 -i rfdiffusion3_path \
                  --chains C --rename-to A -l design_binder
sapia collect chainsel <run_dir> -t table0 -l design_binder

# 2. binder chain out of the PREDICTED complex — on table1, where boltz_path lives
sapia run     chainsel <run_dir> -t table1 -i boltz_path \
                  --chains C --rename-to A -l pred_binder
sapia collect chainsel <run_dir> -t table1 -l pred_binder

# 3. binder vs binder, monomer mode
sapia run     usalign <run_dir> -t table1 \
                  --col-a chainsel_pred_binder_path \
                  --col-b chainsel_design_binder_path \
                  --mm 0 --output-prefix binder_fold
sapia collect usalign <run_dir> -t table1 \
                  --col-a chainsel_pred_binder_path \
                  --col-b chainsel_design_binder_path \
                  --output-prefix binder_fold
```

Three things make that work:

- **chainsel runs on the table that holds the column.** It reads `ctx.df` — the table it
  is updating — and does **not** walk the lineage. `rfdiffusion3_path` is on `table0`, so
  the first extraction targets `table0`.
- **usalign's `--col-b` *does* walk the lineage.** From a `table1` row it finds
  `chainsel_design_binder_path` on the parent `table0` row. Nothing is copied forward.
- **`--mm 0` is monomer mode**, right for a one-chain file. `--mm 1` (usalign's default)
  would be wrong here in the opposite direction: it expects multiple chains.

Read the result as `usalign_binder_fold_TM1` / `_RMSD`: RMSD < 2 Å with TM > 0.9 means
the sequence really folds to the binder you designed. **Run the whole-complex comparison
too** (`usalign --mm 1` on `boltz_path` vs `rfdiffusion3_path`, or `ringfit`) and report
the pair — a good binder fold with a bad pose is a binder that folded and then went
somewhere else, and it is a failed design either way.

Verify before scaling: extract one design, check `chainsel_n_res` equals the binder
length and `chainsel_chains` is what you asked for, and open `chainsel_path` in a viewer.

## Worked example: a split protomer, so ringfit can score the seam

Setup: SlyB is trimmed of its membrane belt, so each target protomer is two segments —
**19–59 (41 aa)** and **107–155 (49 aa)** — with a real 13.4 Å gap. Boltz was given
four target entities plus the binder, so `boltz_path` holds **five** chains: `A`+`C` =
protomer 1, `B`+`D` = protomer 2, `E` = binder (**check yours** — the predictor letters
them in the order `mkcomplex` built them). `ringfit` needs exactly two target chains, so
as it stands `bsa_t1`, `bsa_t2` and `bridge_ratio` cannot be computed at all.

```bash
# 1. merge each protomer into one chain, restoring the auth numbering,
#    and carry the binder through untouched as its own (degenerate) group.
sapia run     chainsel <run_dir> -t table1 -i boltz_path \
                  --merge-groups 'A+C:A,B+D:B,E:E' \
                  --renumber-from 'A:19-59,C:107-155,B:19-59,D:107-155' \
                  -l protomers
sapia collect chainsel <run_dir> -t table1 -l protomers

# expected on every row:
#   chainsel_protomers_n_chains = 3          (A, B, E — not 5)
#   chainsel_protomers_chains   = "A,B,E"
#   chainsel_protomers_n_res    = 90 + 90 + <binder length>

# 2. now ringfit has one chain per protomer, and the numbering the reference uses
sapia run     ringfit <run_dir> -t table1 \
                  -i chainsel_protomers_path \
                  --ref-structure inputs/slyb_assembly.cif \
                  --target-chains A,B --ref-target-chains A,B \
                  --resnum-match resnum
sapia collect ringfit <run_dir> -t table1
```

Why each piece:

- **`--merge-groups`, not `--rename-to`.** `--rename-to A,A,B,B` is a hard error by
  design — an accidental collapse must never be silent. The group spelling is the
  explicit request.
- **The binder needs its own group.** `--merge-groups` replaces `--chains` entirely, so
  any chain you do not name is **dropped**. `E:E` keeps the binder (ringfit's
  `--binder-chains auto` then picks it up as the non-target chain).
- **`--renumber-from`, not `--renumber`.** Boltz numbers each chain 1..N, so both
  segments start at 1 and the default (preserve) would collide — the run would come back
  `error: merged chain A ... has duplicate residue numbers`. `--renumber` would fix the
  collision but give 1..90, which only `ringfit --resnum-match ordinal` could use, and
  ordinal is wrong here: the reference chain also carries residues 60–106, so the i-th
  residues do not correspond. `--renumber-from` restores 19–59 / 107–155, and
  `--resnum-match resnum` pairs equal numbers and tolerates the gap.
- **The `-<last>` is a checksum.** `C:107-155` asserts 49 residues; if Boltz returned a
  different count that design errors instead of shifting the numbering by a few.
- **Order is the group's order.** `A+C:A` writes 19–59 then 107–155. Writing `C+A:A`
  would write 107–155 first — legal, but then the chain's numbers run backwards at the
  seam. Write the segments in sequence order.

Sanity-check one design before scaling: `chainsel_protomers_n_chains == 3`,
`chainsel_protomers_n_res == 180 + binder`, then after ringfit read
`ringfit_resnum_offset` (should be **0** — the merged chain now speaks the reference's
numbering) and `ringfit_seq_match_frac` (should be ~1.0). A non-zero offset or a poor
match means the merge or the ranges were wrong, and `bridge_ratio` would be
wrong-but-plausible.

> **Not yet run.** The tool was authored and its worker exercised locally on synthetic
> structures; it has never been submitted to Modal. Treat the first run as a shakedown —
> two designs, then read `chainsel_status` and `chainsel_n_res` before scaling.
