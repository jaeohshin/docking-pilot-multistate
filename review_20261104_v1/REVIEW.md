# Docking application review — 2026-11-04

**All new analyses in this folder are exploratory.** The original cohort, exclusions, success threshold, endpoints and verdict are retained. Frozen result files and manuscripts were not edited. No new docking, model inference, training, installation, commit or push was performed.

The application supports a bounded result: under this receptor construction and matched docking-call budget, observed qFit alternatives improve the best pose available in the Vina candidate pool, but improvement of the final selected pose remains unconfirmed. It does not yet establish that observed minor-state identity contributes beyond ordinary conformational diversity. This is useful application evidence for the parent paper's distinction between agreement/coverage and discrimination, but it is not downstream validation of the occupancy-trained packer.

## Provenance and actual study locations

- Local parent: `/data/work/sidechain-multistate`, HEAD `193c6c5`. It has untracked `experiments/` and `paper/FOLLOWUP_BINDING_SITES_2026-10-25.md`; preserved.
- No separate completed docking checkout was found in the searched local work directories. The local original experiment folder contains preparation/selection material.
- Existing completed checkout: Syntax, `/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028`, origin `git@github.com:jaeohshin/docking-pilot-multistate.git`.
- Docking HEAD is exactly `8b6814b97fe28e75f38bf0f7ddd923d645401e4d`; no subsequent changes and the checkout was clean. This is the reviewed commit, not a new project.
- Remote parent HEAD: `1c139d0`. Latest manuscript read: `paper/short_20261104/main.tex` and `si.tex`; the parent has existing untracked manuscript build files, a review folder and job scripts, all preserved.
- No applicable AGENTS.md was found in either parent/docking checkout or their checked ancestors. The AGENTS instructions supplied in this conversation were followed. Broad home-tree searches were not used as evidence of absence everywhere.
- Evidence archives, selected-coordinate input hashes, exact-score evidence, manuscript snapshot and SHA256 manifests accompany this review. The 73.7 MB selected-file readback was limited to existing logs/build manifests and coordinates of the ten oracle switch cases; its compressed transfer is approximately 13 MB. No qFit archive or docking environment was installed/downloaded.

## Frozen results and aggregation audit

The analysed main cohort is **128 cases in 39 homology components**. There are 15 analysed backbone-change cases in 7 components, kept separate. `cohort_final.tsv` labels 142 main and 16 backbone-change cases eligible; 14 + 1 did not reach the analysed complete-case cohort. Every missing case has a recorded preparation/build error (`eligible_not_analysed.tsv`), including inconsistent arm preparation and chain/Meeko failures. This limits generalisation to successfully prepared cases; no exclusions were changed here.

`expand_v1/ANALYSIS_PLAN.md` specifies: calculate each paired improvement, average cases within each component, then average component means with equal component weights. The bootstrap resamples these component means with replacement. It is **not** a case-weighted mean with a cluster-robust CI. The implementation matches that estimand. The initial 13-case pilot uses case bootstrap because its representatives are distinct components.

Independent table recalculation reproduces frozen means and bootstrap intervals for the pilot, expanded Vina run, Vina50 co-primaries, GNINA docking alternatives and both rescoring selectors, with the original RNG sequence and 10,000 draws. Newly added subset/rank/interaction intervals use 20,000 component draws and seed 20261104. They are unadjusted exploratory intervals. See `checks.tsv`, `EXTENDED_AUDIT.json` and `TABLE_AUDIT.json`.

Improvement is single RMSD minus multi RMSD; positive favours multi:

| Frozen main contrast | Equal-component mean (Å) | Interval | Interpretation |
|---|---:|---|---|
| Vina50 top-1, P1 | −0.0172 | Bonferroni 97.5% [−0.3060, +0.3179] | Registered improvement criterion not met |
| Vina50 distinct top-5, P2 | +0.1916 | Bonferroni 97.5% [−0.1435, +0.5447] | Registered improvement criterion not met |
| Vina50 oracle | +0.3262 | 95% [+0.1234, +0.5803] | Secondary candidate-pool coverage signal |
| Vina50 → CNNscore rescore top-1 |  +0.1262 | 95% [−0.1854, +0.4887] | Exploratory improvement unresolved |
| GNINA8 default-rescore top-1 | −0.2804 | 95% [−0.9313, +0.2735] | Exploratory improvement unresolved |

The full precision frozen numbers above are read from the saved JSON files; see the audit-generated tables for additional endpoints. The P1/P2 verdict is based on each corrected interval's lower bound exceeding zero. Failure to meet that rule is not equivalence, absence of any possible benefit, or proof that receptor representation is adequate.

**Two averages must not be mixed.** Vina50 apo oracle case-weighted arm means are 2.8386 → 2.5790 Å, a case-weighted paired gain of **0.2597 Å**. The registered equal-component gain is **0.3262 Å**. Likewise, case-weighted top-1 gain is −0.0368 Å, whereas the registered equal-component gain is −0.0172 Å. Both are valid summaries of different populations. Absolute success counts are case counts, not component-weighted rates. `subgroup_effects.tsv` gives both.

## GNINA experiment C: actual mode

The current binary at the recorded execution path reports:

- `gnina v1.1 master:e4cb380+ Built Dec 18 2023`.
- SHA256 `cd5f12d88f63bef2f637ef183beb24f1620ba97944d88173b76fbd373e46d84f`.
- `--cnn_scoring`: **rescore (default)**; refinement and all are separate modes.
- Default output sort: CNNscore.

Evidence: `gnina_runtime_evidence.txt`; frozen `followup_v1/gnina_dock/code/run_gnina_dock.py:69–88` and `code/dock.sbatch`. The launcher passes receptor, converted SDF ligand, box, exhaustiveness 32, modes 9, seed and 8 CPUs. It passes **no** `--cnn_scoring`, `--scoring`, `--cnn`, or `--config`. Logged settings and receptor/seed assignments agree. All 572 C case-arm logs contain eight successful calls.

Thus C is **GNINA's empirical Vina-style search/minimisation, followed by CNN rescoring/selection**, using the default CNN ensemble. The v1.1 source uses the empirical Gaussian/repulsion/hydrophobic/hydrogen-bond/torsional terms when no alternative `--scoring` is specified. CNN-guided refinement or search was not requested. This agrees with the [official v1.1 README](https://github.com/gnina/gnina/tree/v1.1#usage) and [v1.1 main.cpp](https://github.com/gnina/gnina/blob/v1.1/gninasrc/main/main.cpp) (default CNN option around line 1110; empirical default terms around 1439).

The retained wrapper does not save the complete subprocess argv or successful subprocess stdout/stderr. Therefore exact historical stdout and a historical binary hash are unavailable; the binary/version/help above were read now from the recorded path. The frozen command-building code, sbatch, call logs and output score fields consistently support default rescore. No evidence supports a hidden CNN-guided mode. `gnina_runtime_evidence.txt` also records this limitation through the available wrapper/log structure.

Consequences:

- Withdraw “the same CNN generated and selected poses” and “sampling–scoring mismatch was removed/ruled out.”
- Keep the measured C outcomes, relabelled **GNINA default rescore docking**.
- C changes more than a scoring function: GNINA's search implementation, ligand SDF/torsion treatment and candidate retention also differ from Vina1.2/Meeko docking. Equal exhaustiveness and calls do not prove identical search work or isolate CNN-guided search.
- A learned scorer can respond differently to major/minor geometry or holo-like states after training exposure. Pairing removes case-level variation; it does **not** guarantee cancellation of training-set bias. Exposure was not audited here. Absolute and paired CNN contrasts require this qualification.

## Correct first-near-native rank and conditional random baseline

Frozen `analyse_budget50.py:80` computes `argmin(RMSD)+1` after score sorting: the rank of **one closest pose**. Its overall median also includes cases with no RMSD ≤2 Å pose. Consequently the reported 238/216 ranks do not measure the first usable pose, and comparing rank/N to 0.5 cannot diagnose random selection. The analogous GNINA `rank_best` has the same issue.

This review uses all saved poses, sorts by the appropriate score, and reports the first RMSD ≤2 Å rank only when G>0. Missing-good cases are separate. For N candidates and G good candidates, the exact random-permutation distribution is

`Pr(R > r) = choose(N-r,G) / choose(N,G)`, with `E[R]=(N+1)/(G+1)`.

We report the observed rank, that case-specific expectation, their ratio, and the null mid-CDF `Pr(R<observed)+0.5 Pr(R=observed)`. Under uniform random permutations its expectation is 0.5; smaller values mean earlier retrieval. Mid-CDF values are averaged within component then across components conditional on G>0. The toy exact-enumeration checks verify the formula. This avoids interpreting first-success ranks against the wrong null.

| Candidate pool / selector | Apo arm | Has good / no good | First-good rank median | Mean observed rank / mean random expected rank | Component mean null mid-CDF [95%] |
|---|---|---|---:|---|---|
| Vina50 / Vina | single | 60 / 68 | 62.5 | 120.10 / 47.04 | 0.690 [0.567, 0.805] |
| Vina50 / Vina | multi | 68 / 60 | 40.5 | 91.37 / 38.13 | 0.710 [0.596, 0.819] |
| Vina50 / CNNscore | single | 60 / 68 | 1 | 29.93 / 47.04 | 0.294 [0.173, 0.423] |
| Vina50 / CNNscore | multi | 68 / 60 | 1 | 18.18 / 38.13 | 0.249 [0.152, 0.357] |
| GNINA8 / CNNscore | single | 81 / 47 | 1 | 8.46 / 9.97 | 0.327 [0.238, 0.425] |
| GNINA8 / CNNscore | multi | 83 / 45 | 1 | 7.96 / 10.07 | 0.308 [0.231, 0.391] |

Vina's conditional first-good retrieval is late on average under this null, but it still selects 18 good apo poses per arm, compared with random-first-pose expected totals of about 8.22/7.78 on those saved lists. These statements coexist: retrieval is heterogeneous, with early successes and severe burial in other cases. Neither “almost random” nor a claim that Vina has no pose information is justified. CNNscore retrieves good poses earlier on these fixed pools, subject to exposure and candidate-list limitations. This does not establish a multi-minus-single benefit.

**Duplicates and precision.** These are raw candidate lists with repeated modes across calls. The null permutes individual saved candidates; N/G accounts for the multiplicity of good candidates, but this is not inference over independent geometries, nor a uniform draw over distinct pose clusters. We do not equate repeated successful calls with independent biological replication. Vina score ties use deterministic call/pose order, with tie-averaged expectation and min/max ranks also exported. Original GNINA SDF score precision was recovered to remove export-induced ties; the 26 RMSD rows exported as exactly 2.000 Å were recomputed from saved coordinates, without moving them. No main success-generation count changed. The rescoring outputs retain their saved 4-decimal scores, so residual tie ranges are explicit. See `precision_evidence.tsv`, `first_good_rank_per_case.tsv`, and `first_good_rank_summary.tsv` (all four arms and alternative selectors).

Distinct top-5 was already defined before the 50-call analysis as greedy score-order selection with symmetry-aware in-place RMSD >1 Å from every retained pose. It yields 36/128 → 40/128 successes, whereas raw top-5 yields 20 → 20 and top-10 raw yields 22 → 23. This is direct evidence that duplicates materially affect top-N retrieval. The frozen distinct metrics/counts are retained, and their geometry calculation was independently rerun for **all ten switch cases, both apo arms**, reproducing top-5 RMSDs exactly (`distinct_top5_geometry_check.tsv`). A full 128-case distinct-rank null was not computed: that requires cross-pose geometry for the entire cohort, not just the exported RMSD-to-reference table.

## Success transitions and actual apo ensemble size

| Vina50 apo endpoint | Both fail | Single fail → multi success | Single success → multi fail | Both success |
|---|---:|---:|---:|---:|
| Top-1 | 108 | **2** | **2** | 16 |
| Distinct top-5 | 86 | **6** | **2** | 34 |
| Oracle | 59 | **9** | **1** | 59 |

These exactly verify the earlier transition claims. The nine gains span **six** components, not nine independent proteins. Case IDs, components, ligand instances and winning calls/receptors are in the accompanying tables.

The old `has_multiple_states` flag is **apo OR holo** (`analyse_budget50.py:148`); its “states only” subset is not an apo-specific subset. Actual apo prepared receptor counts give:

| Exploratory subset | Cases / components | Apo top-1 gain [95%] | Distinct top-5 gain [95%] | Oracle gain [95%] |
|---|---|---|---|---|
| Apo receptor count =1 | 21 / 12 | 0 [0,0] | 0 [0,0] | 0 [0,0] |
| Apo receptor count ≥2 | 107 / 35 | −0.0226 [−0.3113,+0.3022] | +0.2161 [−0.1144,+0.5502] | +0.3713 [+0.1469,+0.6455] |

All 21 one-receptor control arms have exactly identical score/RMSD tables by `(call,pose)` keys. Their zero difference is a designed computational control, not equivalence evidence for alternative states. The component counts overlap across subsets, so these are not independent subgroup experiments. The 128-case primary analysis remains unchanged.

## Apo–holo interaction, not a comparison of separate significance labels

For each case, calculate `(apo single−multi) − (holo single−multi)`, average within component, then bootstrap components:

| Endpoint | Exploratory interaction (Å) | 95% CI |
|---|---:|---|
| Top-1 | −0.1203 | [−0.4792,+0.2477] |
| Distinct top-5 | +0.3095 | [+0.0017,+0.6228] |
| Oracle | +0.3345 | [+0.1317,+0.5722] |

The **oracle** interaction supports a larger candidate-pool gain on apo in this dataset. It does not establish apo-specific improvement of final top-1 selection. The distinct-top-5 interaction is borderline, unadjusted and exploratory; it does not rescue P2. A non-significant holo effect is not proof that the holo effect is zero. No new endpoint or interaction was promoted to a primary test.

## Dependence and overinterpretation corrections

| Verified issue | Evidence | Required interpretation |
|---|---|---|
| “Three independent replications / all three significant” | All 13 initial cases are included in expanded main; first 8 calls of Vina50 match expanded poses by all 39,676 keys, scores and RMSDs. Initial `docking_v1/RESULTS.md` reports oracle **+0.453 [−0.014,+1.130]**, not +0.565. +0.565 is initial top-1. | Nested cohort/budget extensions, not independent validation. Initial oracle CI includes zero. Exploratory fresh bootstrap gives +0.4529 [−0.0077,+1.1516], also crossing zero. |
| “No effect” / all packers' ceiling | CIs include meaningful positive and negative values; no equivalence margin/test. `build_receptors.py:249–297` retains at most 8 independent occupancy-sampled state combinations, fixes major backbone and filters clashes. | No confirmed selected-pose improvement in this protocol. qFit ensemble is a restricted observed-state benchmark, not an absolute physical or packer upper bound. |
| “Receptor representation is not limiting / only scoring limits performance” | Oracle gains coexist with 60 multi cases without any good Vina50 pose; no matched generic-diversity control; C did not remove mismatch. | Generation and selection both remain incomplete. Selection loses usable poses, but its exclusive causal responsibility is not established. |
| “Observed minor states specifically enable binding conformations” | Some switch structures support steric/holo-like interpretation; many do not; arbitrary diverse states were not compared. | Examples support plausibility. Specific minor-state information versus generic diversity remains unresolved. |
| “Budget effect explains why 8-call positives disappeared” | Nested 8/50 endpoints differ; no mechanistic test of candidate-count explanation. | Describe budget dependence and loss of the earlier positive interval, without treating a non-significant later result as proof that the earlier effect was false. |
| “Paired comparison cancels CNN leakage” | No training-exposure audit; arm geometry differs. | Pairing alone cannot cancel arm-dependent exposure bias. |
| “Apo significant / holo not → different effects” | Direct interaction above. | Oracle interaction supports that limited claim; top-1 interaction does not. |
| Reported arm mean difference equals registered gain | Arm means are case weighted; contrasts average component means. | Label both estimands explicitly; do not subtract case means to reproduce component contrast. |
| Registration heading “50-fold search” | Calls 8→50. | 50 calls per arm is **6.25×** the previous calls, not 50×. Equal calls are an algorithmic budget, not equal measured compute across receptor/engine settings. |

Registration timing is documented by the frozen plans; this audit verifies the saved design/code/results consistency, not an external immutable preregistration timestamp. The review does not rewrite those plans retrospectively.

## Connection to the actual latest parent manuscript

The 2026-11-04 short manuscript's central result is a small occupancy-supervision gain in pooled minor-well discrimination, confirmed on its independent corrected-model holdout, with narrower secondary identification claims. Its Discussion explicitly states that downstream gains in docking/design/binding-associated conformations remain untested. Those packer experiments and this docking experiment must not be conflated: docking here used observed qFit states, **no trained packer**, and relaxed packer-training homology exclusions when expanding the cohort.

A defensible integration is a short application/feasibility subsection (or supplementary application result): improving the **availability** of near-native poses does not guarantee their **retrieval**, paralleling the paper's distinction between agreement/coverage and identification. It provides concrete evidence that downstream evaluation needs generation and selection endpoints. It cannot establish that the independently measured +0.0183 holdout AUC gain improves docking, or that a learned packer has reached a usefulness ceiling. If integrated, replace the blanket “downstream docking remains untested” with a precise distinction: a restricted observed-state docking application was assessed; transfer from the learned packer remains untested.

`MANUSCRIPT_DRAFT.md` contains ready-to-edit English Results/Discussion paragraphs. `STRUCTURE_CASES.md` explains all ten switch cases and their remaining limits. `NEXT_EXPERIMENTS.md` evaluates further work without launching it.

## What was and was not verified

- 1,716 case-arm logs checked: successful calls, seed sequence, round-robin receptor assignment, requested runs/modes/exhaustiveness and counts; no failed logged calls. Vina50 initial eight-call pose table matches the expanded run. Vina50 rescoring uses the same keyed pose table.
- Main numeric tests and success transitions reproduced. Component and original cohort assignment checked; source sequence homology graph was not recomputed with MMseqs, and biological eligibility/RSCC was not re-adjudicated for every expanded ligand.
- Ten switch cases: all saved Vina poses read again, reference RMSDs and distinct top-5 independently recalculated; all-major receptor heavy atoms and identical-seed outputs agree, and changed-receptor backbones remain equal.
- Original GNINA score precision and threshold-boundary RMSDs recovered. Rounded exported GNINA scores can change the selected **RMSD** by up to 0.024 Å in eight arm/endpoint entries (see `rounded_score_tie_differences.tsv`); these are export/tie differences, not evidence that frozen full-precision analysis is wrong. No top-1 success label in those entries changes. Original frozen summaries remain the primary numerical record.
- Static local side-chain comparison uses explicitly mapped residue keys and local backbone fits; its atom-name RMSD does not symmetry-correct equivalent side-chain atom labels and should be treated descriptively. Joint-well identity is an additional coarser check; a bin match is not full-conformation identity.
- Missing evidence: an exposure audit for GNINA training data; joint-state probabilities beyond independent construction; an equal-budget/equal-receptor generic rotamer control; independent ensemble/seed replication; whole-cohort distinct-cluster ranking. These gaps constrain interpretation, not the preserved primary verdict.
