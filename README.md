# docking-pilot-multistate

Does side-chain conformational heterogeneity seen in crystals (qFit-at-scale multiconformer models) change ligand
docking pose recovery? Pilot on non-homologous apo/holo pairs, before any packer comparison.

Design (decided 2026-10-28): 2×2 rigid-receptor docking, **holo / apo backbone × single / multiple side-chain states**,
side chains taken from the qFit-at-scale ground truth (no packer). Single = the qFit major state; multiple = up to
K = 8 receptors of observed chi1/chi2 states grafted on the same major backbone. AutoDock Vina 1.2 (standard
Meeko/PDB2PQR preparation). Main cohort (whole-pocket backbone RMSD ≤ 0.75 Å) and a separate backbone-change stratum.
The follow-up added GNINA 1.1 in two roles, as a rescorer of the Vina pose pools and as a separate search arm.
Not done: flexible side-chain (flexres) comparator arm; packer (single-state B vs multistate C) ensembles.

**The study ran to completion and closed on 2026-10-29.** Its registered endpoints were not met; the oracle
(candidate-pool) gain was significant three times. Read `followup_v1/budget50/RESULTS.md` for the verdict and
`review_20261104_v1/REVIEW.md` for the audit that corrected three readings of it.

## Result as it stands

Apo multi − single, positive favours multiple states. Cases are averaged within homology component and the component
means are bootstrapped; the 13-case pilot bootstraps cases directly, its representatives being one per component.
The 50-call row is the registered confirmatory test.

| run | cohort | calls / arm | top-1 RMSD (primary) | oracle RMSD |
|---|---|---:|---|---|
| `docking_v1` | 13 cases / 13 components | 8 | +0.565 [−0.093, +1.434] — not met | +0.453 [−0.014, +1.130] |
| `expand_v1` | 128 cases / 39 components | 8 | −0.084 [−0.383, +0.240] — not met | **+0.364 [+0.129, +0.657]** |
| `followup_v1/budget50` | 128 cases / 39 components | 50 | −0.017, Bonferroni [−0.306, +0.318] — not met | **+0.326 [+0.123, +0.580]** |

Multiple states put near-native poses into the candidate pool more often (60 → 68 of 128 cases at 50 calls); the
scoring function does not then select them (top-1 success 18 → 18). The gain appears only on apo backbones; on holo
backbones the oracle contrast is 0. The 13-case pilot's +0.565 did not survive expansion: at 128 cases the large
positive and large negative cases balance and the median is 0.

Two positives that appeared at 8 calls — Vina raw top-5 (+0.517) and GNINA rescoring (+0.541) — both vanished at
50 calls, as the frozen plans predicted a budget-dependent metric would. GNINA is nonetheless much the better
selector in absolute terms (apo top-1 success 18 → 40 of the same pose pools), but it lifts both arms equally.

## Pipeline (folders in order)

| step | folder | what | status |
|---|---|---|---|
| 1 | `run_v1`–`run_v3`, `selection_v3` | candidate apo/holo pairs, homology vs train/validation (MMseqs2), geometry QC, strict selection | done (Codex) |
| 2 | `filter_review_v1` | filter attrition and relaxation sensitivity (443 cases) | done (Codex) |
| 3 | `final_review_v1` | fixed input: 79 cases / 23 homology components, deposited mmCIF, metadata | done (Codex) |
| 4 | `final_selection_v1` | ligand roles, biological assembly vs crystal contacts, cofactor/partner checks, apo-minor/holo-major residues with real PDB numbers | done |
| 5 | `state_audit_v1` | qFit pocket side-chain states (all pocket chains, chi1/chi2), whole-pocket apo pose refit | done |
| 6 | `design2x2_v1` | 2×2 representatives (13 main + 4 backbone-change), receptor ensembles by local-frame graft | done |
| 7 | `docking_v1` | Vina environment, ligand/receptor preparation, 544 dockings, first result on 13 cases | done 2026-10-28 |
| 8 | `expand_v1` | cohort expanded by releasing the packer-homology exclusion (no packer is used here), 572 tasks, 143 analysed cases | done 2026-10-29 |
| 9 | `followup_v1/budget50` | registered confirmatory test at 50 calls per arm, co-primaries P1 (top-1) and P2 (distinct top-5) | done 2026-10-29 |
| 10 | `followup_v1/rescore` | Vinardo and GNINA CNNscore reselection of the same saved pose pools, at 8 and 50 calls | done 2026-10-29 |
| 11 | `followup_v1/gnina_dock` | experiment C: GNINA does the search, 8 calls, 572 tasks on a10/a40 | done 2026-10-29 |
| 12 | `review_20261104_v1` | independent audit of all of the above: tables recomputed, three readings corrected, manuscript paragraphs drafted | done 2026-11-04 |

Steps 1–11 were run on this cluster. Step 12 was produced on the desktop machine, reading this tree over SSH at
commit `8b6814b`, and committed here afterwards; `review_20261104_v1/CLUSTER_COPY_NOTE.md` records the two evidence
archives it leaves out and why.

Each folder has `RESULTS.md` (decisions, counts, limitations), `code/` (scripts, Slurm files) and, where written,
`files_sha256.json`. Large regenerable files are not tracked (`.gitignore`); their hashes are in `files_sha256.json`
or in the consuming step's recorded `inputs_sha256`. Docked poses (`**/docked/`) are untracked and live only on this
filesystem: 4.1 GB in all, 2.5 GB of it `expand_v1`.

## What the audit corrected (2026-11-04)

The measured outcomes stand; three descriptions of them do not.

- Experiment C ran GNINA's **default rescore**, not CNN-guided search. The claim that one CNN both generated and
  selected the poses, and that the sampling–scoring mismatch was therefore removed, does not hold.
- The 13-case pilot and the 8-call and 50-call runs are **nested** cohort and budget extensions, not three
  independent replications.
- The reported first-near-native rank is `argmin(RMSD)` after score sorting, i.e. the rank of the single closest
  pose, not of the first usable one. The audit recomputes it against the conditional random-permutation null.

Still open: whether observed minor-state **identity** matters beyond ordinary conformational diversity. That needs
the matched rotamer-prior control designed in `review_20261104_v1/NEXT_EXPERIMENTS.md` (~6,400 apo calls), which has
not been run. `followup_v1/PLAN.md` argues that a zero gain from measured states caps the packer's ceiling at zero;
the audit rejects that, because the qFit ensemble holds only states observed above 0.10 occupancy in chi1 wells and
is not a physical upper bound. That disagreement is unresolved and matters for the parent manuscript's wording.

## Dependencies outside this repository

- Parent research repository `sidechain-multistate` (this folder lives at `experiments/docking_pilot_20261028/` in it),
  commit `0f26be84d9b9fab59ccc7c85551a677afdb3bf63`. Used: `scripts/105_joint_rotamers.py` (wells),
  `scripts/180_apo_holo.py` (apo/holo source rows, index convention), `scripts/177/178` (additive/cofactor lists),
  `results/apo_holo/residues.csv.gz`, `reference/frozen/data/qfit_at_scale/replicate_cache/`.
- qFit-at-scale models: `/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip` (`qfit_PDBs/<pdb>_qFit_010.pdb`).
- Environments: `/gpfs/jaeohshin/work/envs/multistate-analysis` (gemmi 0.7.5, Bio; audits),
  `/gpfs/jaeohshin/work/envs/dock-vina` (Vina 1.2.7, Meeko 0.8.0, RDKit 2026.03.1, OpenBabel 3.2.1, PDB2PQR 3.8.0,
  PROPKA 3.5.1; export in `docking_v1/dock-vina_env.txt`),
  `/gpfs/jaeohshin/work/envs/gnina-bin/gnina` (GNINA v1.1 `master:e4cb380+`, built Dec 18 2023; hash and `--help`
  output in `review_20261104_v1/gnina_runtime_evidence.txt`).
- Public RCSB files (mmCIF, wwPDB validation reports), fetched by the scripts.

## Rules kept throughout

- Metrics, thresholds, strata, cases and bootstrap units were fixed in `ANALYSIS_PLAN.md` / `PLAN.md` **before** the
  results existed, and were not changed afterwards. This is what caught the two 8-call positives above as budget
  effects rather than reporting them as findings.
- The apo-minor/holo-major chi1 signal and docking results are never used to select cases or representatives;
  signal-positive cases are reported separately.
- Source residue indices (`i`) are 0-based qFit polymer positions, not PDB numbers; tables show author numbers.
- No docking, inference or training without explicit approval; earlier results are not overwritten. Each run writes
  into its own folder, and superseded attempts (`run_v1`/`run_v2`, `receptors_v0_superseded`) are kept as they were.
- Scope ends at pose RMSD. Virtual screening and enrichment were ruled out in advance (`followup_v1/PLAN.md`) so that
  a negative result could not be widened into a search for a positive one.
