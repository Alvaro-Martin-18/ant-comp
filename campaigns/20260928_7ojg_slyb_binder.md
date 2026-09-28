# Campaign: composite-epitope binder against SlyB (7OJG) protomers A/B

**Date:** 2026-09-28 · **run_dir:** `outputs/20260928_133724_7ojg_patch1` · **recon:** `outputs/20260928_131320_slyb_recon`
**Brief:** design a binder contacting *both* subunits A and B of PDB 7OJG.
**Outcome: no validated binder.** 4–5 designs pass every *structural* gate — they fold, contact
both protomers at the designed epitope, and are physically compatible with the real 11-mer —
but Boltz's own interface confidence rejects all of them (best binder–target ipTM **0.357**,
best interface PAE **19.5 Å**, and **every** design *loses* pLDDT when docked). See Result 5.
The structural machinery is now validated end to end and the failure is localised and
measurable, which the previous attempt at this brief was not.

**Confound to resolve before trusting Result 5:** every prediction ran with `msa: empty`
(`--use-msa-server` was never enabled). Correct for the de-novo binder, **wrong for the
target** — SlyB has a deep MSA. Some of the uniformly poor interface confidence may be the
setup rather than the designs.

---

## 1. The target

7OJG is **not a hetero-complex**. It is undecameric SlyB — a **C11 homo-oligomeric
outer-membrane lipoprotein ring**, 11 identical protomers (chains A–K), cryo-EM 3.40 Å.
"Subunits A and B" are two adjacent identical protomers, so the brief is a **composite
epitope straddling the A/B seam**.

Per protomer: residues **18–155 modelled** (138 aa), no internal gaps, `label_seq_id ==
auth_seq_id`. Residues 1–17 are the unmodelled lipoprotein signal peptide. Ligands per
protomer: PLM ×3, GOL ×1, L8Z ×1 (Ra-LPS), LPP ×1 — 66 het residues in the assembly.

### Ring geometry

| quantity | value |
| --- | --- |
| symmetry axis | verified two ways: inertia tensor, and the A→B superposition rotation of **32.702°** (360/11 = 32.727°), `\|dot\| = 0.999998` |
| ring height (CA) | 73.86 Å |
| outer radius | 53.18 Å (protein, all atoms 54.93) |
| pore radius | 24.87 Å |
| lipid slab | z ≈ **+9 to +50 Å** (from the ligands themselves) |

### The membrane belt — residues 62–106

Identified three independent ways that agree:

- **direct lipid contact** (< 4.5 Å of PLM/LPP/L8Z/GOL): contiguous core **62–106**, axial z **+8.5 to +41.9 Å**
- **hydrophobicity**: 59–107 is 0.755 hydrophobic vs 0.390 / 0.396 either side; Gly fraction 0.306
- **seven GxxxG motifs**, all inside the belt — a textbook Gly-zipper

Also **Cys18**, the S-diacylglyceryl lipoprotein anchor, contacting GOL and PLM.

### Key geometric facts that shaped every later decision

1. **The belt is not at the widest radius.** Lipid-facing core sits at r ≈ 25–38 Å while the
   *soluble* skirt flares to r ≈ 53 Å at z < 0. The bindable surface is the most radially
   accessible surface on the particle and points away from the membrane.
2. **The exposed face points down, not out.** Most high-SASA skirt residues (rel_sasa up to
   0.98) are exposed **axially along −z**; the ring is flat-bottomed. A radial `sc_out` test
   is the wrong approach-vector criterion below z ≈ −10.
3. **The ring is too narrow to bind one seam in isolation.** Arc per protomer at skirt radius
   is only **17–27 Å**, and every candidate epitope patch is **18–20 Å across**. Any binder
   large enough to straddle A/B necessarily comes within a few Å of the next seam. This is
   11-fold geometry, not an artefact.

---

## 2. Target definition

**Trim: delete 60–106 and Cys18; keep 19–59 + 107–155** (90 of 138 residues per protomer).

Gated on a check run *before* anything was built: does the soluble body survive losing the belt?

| contact set (chain A, heavy atoms < 5 Å) | atom pairs | residue pairs |
| --- | --- | --- |
| **19–59 ↔ 107–155 (both retained)** | **633** | **97** |
| 19–59 ↔ belt 60–106 | 63 | 6 |
| 107–155 ↔ belt | 59 | 7 |

The two retained segments are **hydrogen-bonded partners in the same β-sheet** — 18 backbone
H-bonds in reciprocal antiparallel ladder pairs (V46↔E115, T44↔R117, S41↔L143, V51↔G111,
I53↔T109, N39↔S145). Belt contacts are only the covalent neighbours at the cut (54–59, 107–108).

Compactness ratio (Rg / 2.2·N^0.38): **trimmed unit 1.14** vs **intact protomer 1.64**. The
belt is what makes SlyB elongated; removing it leaves a *better*-behaved globule.

**Padding: K, A, B, C — four protomers, 360 residues.** K and C exist solely to cap the
artificial faces that cutting A+B out of a ring creates. See §5, result 1: this was not a
precaution, it was load-bearing.

**Chain break: mandatory.** CA(59)–CA(107) = **13.44 Å**. Any tool given one fused 90-residue
sequence per protomer will try to bond residues 13.4 Å apart — silently, by distorting the
domain. Every protomer must be presented as **two chains** (41 aa + 49 aa) to any predictor.

---

## 3. Epitope

Candidate patches were found by SASA in the 4-mer context (rel_sasa ≥ 0.25, axial z < +5 Å,
within 14 Å of the partner chain), then clustered. *Single-linkage clustering percolated* —
one cluster of 73 residues, 49 Å across — because the soluble skirt is one continuous exposed
surface; seed-centred 11 Å footprint patches with non-max suppression were used instead.

**Chosen: patch 1** — best balance of the candidate set.

| | |
| --- | --- |
| chain A | T123, N146, G147, S148, Q149 |
| chain B | Q36, V37, N39, N146, G147 |
| balance `min/max` per-chain area | **0.488** (323 Å² on A vs 338 Å² on B) |
| centroid | z = −18.50, r = 30.49 |
| clearance to belt plane (z = +8.5) | **19.68 Å** |
| hotspots given to rfd3 | **A146, A148, B37, B39** |

**Two better-balanced-looking patches were rejected on physical grounds.** Patch 2 (balance
0.410, 787 Å²) sits at centroid z = −1.0 with an 18 Å footprint — a binder there reaches
*into* the membrane — and 2.4 Å from the 59/107 cut, i.e. partly binding a surface that exists
only because 47 residues were deleted. That is precisely the error that sank the previous
campaign, and it was the second-best-scoring option.

---

## 4. The pipeline as run

```
/runs/inputs/                     7ojg_ref_assembly.cif    full 11-mer + lipids (ringfit ref)
                                  7ojg_soluble_KABC.pdb    trimmed target, 4 × 90 res
                                  7ojg_tmpl_kabc.cif/.yaml 8-chain forced template  ✅
                                  7ojg_tmpl_ab.cif/.yaml   4-chain template         ✗ dead branch
                                  patch1_hotspots.yaml     A146,A148,B37,B39
```

| # | table | tool · leaf | action | rows | filter |
| --- | --- | --- | --- | --- | --- |
| 1 | → **table0** | `rfdiffusion3` · `binder_patch1_v1` | create | **8** | — |
| 2 | table0 | `chainsel` · `design_binder` | update | 8 | — |
| 3 | table0 → **table1** | `proteinmpnn` · `mpnn_t02` | create | **32** | — |
| 4 | table1 | `boltz` · `boltz_binder_alone` | update | 32 | — |
| 5 | table1 | `usalign` · prefix `binder_fold` | update | 32 | — |
| 6 | table1 | `mkcomplex` · `mkc_ab` / `mkc_kabc` | update | 32 each | — |
| 7–9 | table1 | `boltz`/`chainsel`/`ringfit` · `*_ab` + `*_kabc` | update | **2** each | `pilot2_filter.py` |
| 10–12 | table1 | `boltz`/`chainsel`/`ringfit` · `cofold_kabc`, `prot_kabc`, `rf_kabc` | update | **+25 → 27** | `fold_pass_filter.py` |
| 13 | table1 | `cms` · `cms_ab` / `cms_all` | update | **9** each | `landed9_filter.py` |
| 14 | table1 | `pyrosetta` · `pyr` | update | **9** | same |

**Filters**

| filter | rule | kept |
| --- | --- | --- |
| `pilot2_filter.py` | two named designs | 2 / 32 |
| `filters/fold_pass_filter.py` | `binder_alone_complex_plddt ≥ 0.80` **and** `binder_fold_TM1 ≥ 0.80` | **27** / 32 |
| `filters/landed9_filter.py` | `cofold_kabc_confidence_score ≥ 0.67`, cross-checked against measured target RMSD ≤ 2.0 Å (sets identical) | **9** / 27 |

**Funnel:** 8 backbones → 32 sequences → 27 fold → **9 targets landed** → 4–5 shortlisted.

**Cost:** one A100 (~5.5 min), one L4, 32 small monomer predictions, ~15 A10-minutes of
co-folding, everything else CPU.

### rfd3 invocation

```
--contigs 'A19-59,A107-155,/0,B19-59,B107-155,/0,K19-59,K107-155,/0,C19-59,C107-155,/0,70-90'
--extra-spec patch1_hotspots.yaml --num-designs 1 --set n_batches=8 --modal-gpu A100
```

A **bare comma keeps two motif segments in the same output chain** and runs numbering straight
through; `/0` is the only chain-break literal. Output chains were renumbered 1..N, mapping
verified by Kabsch over all 24 permutations: **(A,B,K,C) at 0.050 Å vs ≥16.0 Å for every
alternative**. Binder = chain E.

The 13.4 Å belt gap **survived geometrically** (CA41–CA42 = 13.407 Å) but **not in numbering**
(41→42 contiguous) — which is why the split-chain treatment in `mkcomplex` was mandatory.

---

## 5. Results

### Result 1 — deleting the decoy is not interchangeable with detecting it

The same design, same sequence, under two target definitions:

| | AB (K, C absent) | KABC (complete) |
| --- | --- | --- |
| `bsa_t1` / `bsa_t2` | 1238 / **0** Å² | 696 / 482 Å² |
| `bridge_ratio` | **0.000** | **0.693** |
| `hotspot_recall` | 0.00 | **0.75** |
| `n_clash` | **312** over 55 of 86 residues | 2 |
| `min_dist_ring` | **0.278 Å** | 1.34 Å |

diff_6_f1 under AB did not drift — it **docked into the void where protomer K belongs**: 1027
heavy-atom pairs under 2.5 Å against chain K specifically, interpenetrating to 0.28 Å, 64% of
its residues involved. It abandoned the seam entirely to bury 1238 Å² in a protomer-shaped hole
that exists only because AB omits it.

**AB scored 0/2 at the designed epitope** (`hotspot_recall` 0.00 even on its clean row) while
producing one impossible pose. KABC: 2/2 with real hotspot recall. *A surface you removed by
cutting the target is still, to the predictor, a surface.*

### Result 2 — the best-looking numbers came from the broken predictions

Target landing was **bimodal with an empty band**: 9 rows at 1.162–1.242 Å, nothing between,
18 rows from 4.73 to 13.24 Å. And:

| name | bsa_total | hotspot_recall | bridge | **target RMSD** |
| --- | --- | --- | --- | --- |
| diff_3_f1 | 2783 | 0.50 | **0.954** | **5.20** ✗ |
| diff_3_f3 | 2739 | **1.00** | 0.410 | **6.95** ✗ |
| diff_7_f1 | 2248 | **1.00** | 0.565 | **8.44** ✗ |
| diff_1_f4 | 2664 | **1.00** | 0.180 | **7.79** ✗ |
| diff_5_f3 | 2523 | **1.00** | 0.277 | **7.43** ✗ |

**Every design with perfect hotspot recall and a 2200–2800 Å² interface had a collapsed
target** — binders being engulfed by a ring folding around them (`n_clash` 32, 374, 388).
Ranked on interface size or hotspot recall without a geometric gate, the five worst designs
would have been selected as the five best.

### Result 3 — whole-complex confidence is a target-assembly detector, not a binder score

| metric | Pearson r vs target RMSD |
| --- | --- |
| `confidence_score` | **−0.896** |
| `complex_plddt` | −0.872 |
| `iptm` | −0.839 |

`confidence_score ≥ 0.67` classified **all 27 correctly** (landed 0.706–0.736, distorted
0.542–0.634). This is mechanistically expected: the complex has 9 chains = 36 interchain pairs,
of which **35 are target–target and template-forced** and exactly **one** is binder–target. The
collected `iptm` is therefore overwhelmingly a statement about template reproduction.

Excellent gate. Useless as a binder score. **The per-pair value that would be a binder score
(`pair_chains_iptm`) is written by Boltz and is not collected by prosapia — see §7.**

### Result 4 — fold and pose are independent

diff_6_f1 had the **best self-consistency in the campaign** (TM 0.971, RMSD 0.58 Å) and, under
the wrong target definition, a physically impossible pose. Fold quality carries no information
about pose quality.

### Result 5 — the predictor does not believe any of these interfaces

`pair_chains_iptm` and `chains_ptm` are in Boltz's confidence JSON (9×9, chain index positional
from the input YAML, index 8 = chain I = binder). The matrix is **not symmetric** (max asymmetry
0.079); values below are symmetrised means. Also read: per-residue pLDDT from
`B_iso_or_equiv` (cross-checked against `plddt_*.npz` to 7e-5) and PAE from `pae_*.npz`.

| | best in set | typical "predicted binder" |
| --- | --- | --- |
| binder–target ipTM | **0.357** (diff_4_f4) | ~0.5–0.6 |
| interface PAE to the seam | **19.5 Å** | < 10 Å |
| binder pLDDT, alone → in complex | **−7.5 to −23.3** (mean −13.7) | should *increase* |

**Not one of 27 designs gains confidence from being docked.** The binder's own fold is confident
everywhere (`chains_ptm[I]` 0.87–0.99, intra-binder PAE 1.5–6.4 Å); only the placement is not.

**`confidence_score` is anti-informative once past the landing gate.** Spearman vs binder ipTM:
**−0.250** among the 9 landed, versus +0.679 across all 27 — i.e. it separates "target landed"
from "target didn't" and carries no binder signal thereafter. **diff_2_f4** is the clean
demonstration: highest `confidence_score` of all 27 (0.736) and highest target pLDDT (81.5),
yet near-bottom binder interface confidence (0.158).

**Per-pair ipTM did not deliver the better ranking it promised.** Spearman +0.883 against the
collected `iptm`; the top two are unchanged; it is *worse* at predicting `bridge_ratio`
(+0.27 vs +0.50) and marginally worse for `hotspot_recall` (+0.867 vs +0.913). The
seam-vs-neighbour split is rank-identical (rho = +1.000; seam/neigh ratio a uniform 1.36 ± 0.15).
**Binder pLDDT is the one orthogonal axis** (rho +0.05 with iptm) and re-ranks strongly.

Boltz interface confidence and geometric straddling are close to **independent** here:
diff_2_f4 has `bridge_ratio` 0.847 with the second-worst binder ipTM.

*Next refinement, not done:* `pair_chains_iptm` is whole-chain vs whole-chain, so a binder
touching a small patch of a 41-residue chain is diluted the same way the global number is. An
**interface-restricted PAE** (contacting residue pairs only) is computable from the npz already
on the Volume and should separate these better than anything above.

### Shortlist (all gated on target RMSD ≤ 2 Å)

| design | bridge | hotspots | Δ CMS | cms_ab | sc | if_dG | packstat | relax_ca_rmsd |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **diff_6_f3** | **0.864** | 2/4 | **0.03** | 240 | 0.370 | −13.1 | 0.590 | 1.30 |
| **diff_4_f1** | **0.869** | **3/4** | 75.6 | 224 | 0.496 | −11.8 | **0.597** | 1.11 |
| **diff_6_f1** | 0.693 | **3/4** | 147.6 | 265 | 0.500 | **−28.4** | 0.557 | 1.29 |
| **diff_4_f4** | 0.353 | **3/4** | 42.2 | **297** | **0.532** | −23.9 | 0.542 | **1.01** |
| diff_4_f2 | 0.376 | 3/4 | 21.1 | 250 | 0.501 | −22.4 | 0.544 | 1.68 |

All are 28–32 Å clear of lipid. Δ CMS = `cms_all_target − cms_ab_target`, i.e. interface with
the *neighbouring* protomers K/C rather than the designed A/B seam.

**The two-target CMS comparison earned its keep**: diff_6_f4 puts only **50.8 of its 253.8 Å²**
on the designed seam — its interface is essentially all with K and C. Its `sc` of 0.7886 is the
set maximum and is an artefact of a 46 Å² trimmed area. diff_2_f4 and diff_6_f3 have Δ ≈ 0.00,
entirely on-target.

**Honest verdict on quality:** `sc` 0.36–0.53 against a well-packed reference of ~0.5–0.75;
`packstat` below 0.6 for 7 of 9; `if_dG` −11.8 to −28.4 REU. These pass every structural gate
and are **not** strong binders.

### Sequence composition — a problem no structural metric caught

| | designs | typical soluble protein |
| --- | --- | --- |
| % charged | median 10.3, **min 4.5** | ~25 |
| net charge | **−4 to +3, median −1** | — |
| % aromatic | median 3.35 | ~8 |
| % TSVIG | **median 57.1, max 70.5** | ~30 |

Seven of eight backbones are all-β (the epitope is a β-sandwich skirt, so rfd3 answered with β);
ProteinMPNN is at its worst there. Boltz folds these at pLDDT 0.96 regardless — self-consistency
cannot see solubility.

A prediction made in advance and **falsified**: interface `buried_unsat` was expected to be
high. Total counts are 25–44, but the interface-attributable share (`if_delta_unsat`) is only
**4–16, median 8**. The composition problem is real but manifests as *packing quality*
(`sc`, `packstat`), not unsatisfied interface polars.

---

## 6. Traps found, with their signatures

| # | trap | signature | status |
| --- | --- | --- | --- |
| 1 | `usalign.sh` merged stderr via `2>&1`, then `grep -v '^#' \| head -1` grabbed USalign's `-o` warning | `ERROR: expected 9 metrics, got -1` on **every** row, **every `.exit` = 0** | fixed, `199a103` |
| 2 | `GROUPS` is a bash special variable; `GROUPS=$(…)` fails rc=1 under `set -e` | `.exit` ≠ 0 with **`.out` and `.err` both 0 bytes**, nothing in `modal app logs` | fixed; recorded in `authoring-a-tool` |
| 3 | manifest not keyed by leaf → concurrent same-tool runs overwrote each other | all `.exit` = 0, both collects report full row counts, 7/32 rows silently wrong | fixed upstream (now `{table}_{leaf}_manifest.txt`) |
| 4 | gemmi `subchain` labels carried from source → `setup_entities()` merged protomer halves | **`_entity_poly_seq` = 49 for 41-residue chains**; every other check passed | caught pre-spend |
| 5 | Boltz names template chains by `label_asym_id`, not auth chain ID | `Template chain A … is not one of the protein chains!` | caught pre-spend |
| 6 | Boltz template falls back to **sequence search** unless `chain_id` *and* `template_id` are both given and equal length | silent mis-assignment; all four SEG1 chains are sequence-identical | avoided by construction |
| 7 | `build_boltz_manifest` never clears `boltz_shards/` | stale YAMLs re-predict rows you meant to skip | caught by input-file count |
| 8 | `"Submitting N designs"` counts **shards** for boltz, and **parameter groups** for proteinmpnn | reads as 1 when 8 or 32 designs ran | verify by counting `*_inputs/` files |
| 9 | `name` is the DataFrame **index**, not a column | `KeyError: 'name'` in a filter module | — |
| 10 | merged chains both numbered from 1 → `ringfit --resnum-match resnum` builds `{resnum: index}` | would silently drop half a protomer | now a hard error in `chainsel` |

---

## 7. Known gaps

1. **Binder-specific confidence is not collected** — `pair_chains_iptm` and `chains_ptm` are in
   Boltz's confidence JSON but prosapia lifts only the scalar summaries. Extracted by script
   (`analysis_binder_conf/`) and reported in Result 5; **still not a collected column**, so it
   is not available to a filter. Worth adding to the boltz collector.
2. **No MSA was used anywhere** (`msa: empty`, `--use-msa-server` never enabled). Correct for a
   de-novo binder; **wrong for the target**, which is a natural protein with a deep MSA. This
   plausibly depresses both target assembly (only 9/27 landed) and interface confidence
   (Result 5). **Untested and it is the cheapest high-value control available.**
3. **No per-chain relax RMSD.** `pyrosetta_relax_ca_rmsd` is whole-pose over 435–450 residues,
   dominated by the eight excised target chains and their artificial termini. 8 of 9 rows exceed
   1 Å under *constrained* relax, and it cannot be determined whether that is the target settling
   or the binder shifting. Every `if_dG` here carries that asterisk.
4. **No fold check inside the complex.** Binder fold was measured only in isolation. Extracting
   chain I from the co-fold and running `usalign --mm 0` against the designed backbone would
   measure the fold in the context that matters. Pure CPU, not yet run.
5. **rfd3 cannot be steered toward helix.** No SS field on `DesignInputSpecification`
   (`extra="forbid"`), the helix/sheet annotations are training-only, and the released
   `rfd3_net.yaml` `token_1d_features` has no SS channel. Oversample and select is the only route.
   `is_non_loopy` is wired but does not distinguish helix from sheet.
6. **Absolute Rosetta scores are not meaningful here** — artificial termini at every belt cut add
   terminal charges that do not exist in the real protein. Within-batch ranking only.

---

## 8. Round 2

Ordered by information gained per unit cost. **The first item is a control, not a design step** —
Result 5 currently cannot distinguish "these binders are bad" from "the prediction setup was
handicapped", and everything downstream depends on which it is.

1. **Re-run the co-fold on the 9 landed designs with an MSA for the target** (`--use-msa-server`;
   the binder is de novo and stays MSA-free). ~9 predictions. If interface ipTM and interface
   PAE improve materially, the current verdict is an artefact of my setup and the shortlist is
   live again. If they do not, the designs are genuinely not believed and no amount of sequence
   redesign on these backbones will fix it. **Do this before spending anything else.**
2. **Interface-restricted PAE** on the existing predictions (gap 1 note in Result 5). Pure CPU,
   files already on the Volume, and the most likely source of real discrimination among the 9.
3. **The two cheap checks still owed**: extract binder chain I from the co-fold and `usalign
   --mm 0` it against the designed backbone (fold *in complex*, gap 4); and per-chain relax RMSD
   so `if_dG` stops carrying an asterisk (gap 3).
4. **ProteinMPNN with soluble weights + AA bias, on the diff_4 / diff_6 backbones.** Zero code:
   `--set '--use_soluble_model'`, `--bias-aa` to suppress T/S/V and enrich D/E/K and aromatics,
   `--set '--omit_AAs C'`, `--sampling-temp 0.25–0.3` (a single value — a space-separated list
   is word-split by bash and breaks argparse). Attacks the packing weakness (`sc` 0.36–0.53,
   `packstat` < 0.6) and leaves the validated fold/pose machinery untouched.
5. **diff_4 and diff_6 are the productive backbones** — 4/4 and 3/4 targets landed; between them
   they supply every shortlisted design. Oversample from the same epitope.
6. **Drive a larger interface at design time** — current interfaces are 14–22 residues and
   508–837 Å² as designed, at or below the low end for a minibinder. More hotspots, or a longer
   binder against a patch at larger ring radius where the arc per protomer is widest (patch 3 at
   r = 42.9, patch 5 at r = 47.7).

**Add to the boltz collector** regardless: `pair_chains_iptm` and `chains_ptm`. Without them the
binder-interface confidence cannot be used in a `--filter`, which is where it belongs.
