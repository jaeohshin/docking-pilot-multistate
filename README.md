# docking-pilot-multistate

Does side-chain conformational heterogeneity seen in crystals (qFit-at-scale multiconformer models) change ligand
docking pose recovery? Pilot on non-homologous apo/holo pairs, before any packer comparison.

Design (decided 2026-10-28): 2×2 rigid-receptor docking, **holo / apo backbone × single / multiple side-chain states**,
side chains taken from the qFit-at-scale ground truth (no packer). Single = the qFit major state; multiple = up to
K = 8 receptors of observed chi1/chi2 states grafted on the same major backbone. AutoDock Vina 1.2 (standard
Meeko/PDB2PQR preparation). Main cohort (whole-pocket backbone RMSD ≤ 0.75 Å) and a separate backbone-change stratum.
Later, not yet: flexible side-chain (flexres) comparator arm; packer (single-state B vs multistate C) ensembles.

## Pipeline (folders in order)

| step | folder | what | status |
|---|---|---|---|
| 1 | `run_v1`–`run_v3`, `selection_v3` | candidate apo/holo pairs, homology vs train/validation (MMseqs2), geometry QC, strict selection | done (Codex) |
| 2 | `filter_review_v1` | filter attrition and relaxation sensitivity (443 cases) | done (Codex) |
| 3 | `final_review_v1` | fixed input: 79 cases / 23 homology components, deposited mmCIF, metadata | done (Codex) |
| 4 | `final_selection_v1` | ligand roles, biological assembly vs crystal contacts, cofactor/partner checks, apo-minor/holo-major residues with real PDB numbers | done |
| 5 | `state_audit_v1` | qFit pocket side-chain states (all pocket chains, chi1/chi2), whole-pocket apo pose refit | done |
| 6 | `design2x2_v1` | 2×2 representatives (13 main + 4 backbone-change), receptor ensembles | in progress |
| 7 | `docking_v1` | Vina environment, ligand/receptor preparation, docking (not run) | in progress |

Each folder has `RESULTS.md` (decisions, counts, limitations), `code/` (scripts, Slurm files) and, where written,
`files_sha256.json`. Large regenerable files are not tracked (`.gitignore`); their hashes are in `files_sha256.json`
or in the consuming step's recorded `inputs_sha256`.

## Dependencies outside this repository

- Parent research repository `sidechain-multistate` (this folder lives at `experiments/docking_pilot_20261028/` in it),
  commit `0f26be84d9b9fab59ccc7c85551a677afdb3bf63`. Used: `scripts/105_joint_rotamers.py` (wells),
  `scripts/180_apo_holo.py` (apo/holo source rows, index convention), `scripts/177/178` (additive/cofactor lists),
  `results/apo_holo/residues.csv.gz`, `reference/frozen/data/qfit_at_scale/replicate_cache/`.
- qFit-at-scale models: `/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip` (`qfit_PDBs/<pdb>_qFit_010.pdb`).
- Environments: `/gpfs/jaeohshin/work/envs/multistate-analysis` (gemmi 0.7.5, Bio; audits),
  `/gpfs/jaeohshin/work/envs/dock-vina` (Vina 1.2.7, Meeko 0.8.0, RDKit 2026.03.1, OpenBabel 3.2.1, PDB2PQR 3.8.0,
  PROPKA 3.5.1; export in `docking_v1/dock-vina_env.txt`).
- Public RCSB files (mmCIF, wwPDB validation reports), fetched by the scripts.

## Rules kept throughout

- The apo-minor/holo-major chi1 signal and docking results are never used to select cases or representatives;
  signal-positive cases are reported separately.
- Source residue indices (`i`) are 0-based qFit polymer positions, not PDB numbers; tables show author numbers.
- No docking, inference or training without explicit approval; earlier results are not overwritten.
