# Docking exploratory review, 2026-11-04

Start with [REVIEW.md](REVIEW.md). The original primary cohort and verdict are retained. All new analyses are exploratory. No new docking/model inference/training, environment installation, commit or push was performed.

## Main deliverables

- [REVIEW.md](REVIEW.md): verified problems, evidence/code locations, original numeric audit, corrected ranks, subgroups, interactions, parent manuscript connection and remaining limits.
- [STRUCTURE_CASES.md](STRUCTURE_CASES.md): all nine oracle gains and one loss, residue/state interpretations, seed evidence and missing evidence.
- [MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md): English Results/Discussion and a precise edit to the latest parent manuscript's downstream-application limitation.
- [NEXT_EXPERIMENTS.md](NEXT_EXPERIMENTS.md): matched rotamer-prior control, interpretation, budget and priorities; design only.
- `success_transitions.tsv`, `frozen_main_contrasts.tsv`, `absolute_results.tsv`, `subgroup_effects.tsv`, `interaction_effects.tsv`, `component_effects.tsv`: estimands/counts, with case and component means labelled.
- `first_good_rank_per_case.tsv` and `first_good_rank_summary.tsv`: first ≤2 Å retrieval, G=0 cases, N/G-aware random baseline and tie diagnostics for Vina, CNN/Vinardo rescoring, and GNINA selectors.
- `structure_case_summary.tsv`, `changed_residues.tsv`, `receptor_call_results.tsv`, `matched_seed_calls.tsv`: all ten requested cases, all their receptors/calls, and best-combination residues.
- `ensemble_inventory.tsv`, `eligible_not_analysed.tsv`, `initial_overlap.tsv`, `leave_one_component.tsv`: preparation/cohort/ensemble/dependence audit.
- `TABLE_AUDIT.json`, `EXTENDED_AUDIT.json`, `LOG_AUDIT.json`, `structure_checks.tsv`, `distinct_top5_geometry_check.tsv`, `checks.tsv`: validation evidence.
- `precision_evidence.tsv`, `threshold_boundary.tsv`, `rounded_score_tie_differences.tsv`: export precision and threshold-boundary audit. Frozen metrics retain full-precision originals.
- `gnina_runtime_evidence.txt`, `PROVENANCE.json`, `SHA256SUMS.json`: mode/version/help and provenance.

## Figures and coordinates

- `generation_selection.png`: main Vina50 success counts and gains/losses.
- `first_good_rank.png`: observed first-good rank against case-specific random expectation; raw repeated poses retained.
- `structure_examples.png`: coordinate projections of Met237 steric relief and Ile621 holo-like alternative. Distances are 3D.
- `structures/<case_id>/`: original built PDBs, saved oracle pose SDFs, reference ligand, build mapping and optional PyMOL commands. No viewer installation is required to read the tabular evidence.

All three PNGs were visually inspected. Structures are illustrative; no density map was retrieved or validated.

## Reproduce from the existing completed repository

Every script writes to its supplied review output folder and reads the source root. Do not invoke frozen original analysis scripts: those overwrite frozen outputs.

On the existing Syntax checkout, use its already installed `dock-vina` Python for geometry; use a Python with NumPy/Pandas/SciPy/Matplotlib for table analysis. The local environments used here are recorded in `PROVENANCE.json`.

```bash
python reanalyse.py --source /path/to/docking-pilot-root --out /path/to/new-review-output
python audit_logs.py --source /path/to/docking-pilot-root --out /path/to/new-review-output
python extended_audit.py --source /path/to/docking-pilot-root --out /path/to/new-review-output
python structure_cases.py --source /path/to/docking-pilot-root --out /path/to/new-review-output
```

`reanalyse.py` expects `precision_evidence.tsv` copied into the output directory before running; that file contains original GNINA SDF properties and coordinate-based boundary RMSDs already collected in this review. It does not run a learned scorer. `make_structure_figure.py` reads the `structures/` folder beside itself.

For an offline reproduction, make a temporary source-evidence directory and extract, in this order, `evidence/docking_source_8b6814b.tar.gz` and `evidence/selected_existing.tar.gz` into it. These are **evidence snapshots**, not a new checkout. Then pass that directory as `--source`. Original tracked `inputs` symlinks refer to the cluster; the review scripts use actual selected files under `expand_v1/inputs` and do not execute those symlinks. The snapshot has all inputs needed for the review's table/log/ten-case analyses, but not the full docking dataset/environment needed to repeat docking. Do not run archived docking launchers.

`collect_existing.py` and `collect_precision.py` document the read-only collection from the discovered Syntax path. Their products are already included; they need not be rerun. They require the existing cluster data and environment and never start a docking engine or learned model. `COLLECTED_MANIFEST.json` records the selected files' original relative paths, sizes and SHA256 hashes.

`evidence/parent_manuscript_1c139d0.tar.gz` and the extracted TeX copies preserve the latest manuscript text that was read. Original manuscript sources were not edited.

## Statistical conventions

Frozen reproduction uses the original sequential RNG and 10,000 draws. New exploratory confidence intervals use 20,000 equal-component bootstrap draws, seed 20261104. Primary P1/P2 retain their original Bonferroni 97.5% intervals. Additional intervals are unadjusted and are not new primary tests. The rank null conditions on the actual saved candidate list; repeated poses/calls are not independent biological replicates.
