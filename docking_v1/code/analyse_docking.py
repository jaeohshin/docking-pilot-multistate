"""Score the 2x2 docking pilot exactly as ANALYSIS_PLAN.md fixed it, before any result was produced.

Primary metric: the symmetry-aware heavy-atom RMSD of the best-scoring pose over an arm's 8 calls, against the
reference ligand in that backbone's own frame (holo: deposited pose; apo: the whole-pocket refit pose). RMSD is
computed in place (rdMolAlign.CalcRMS), never after alignment. Success is <= 2 A.
Secondary: the lowest RMSD any pose reached (oracle), the best score, and the success rate.

Primary contrast: multiple minus single on the apo backbone, as a paired per-case improvement
(single RMSD - multiple RMSD, so positive means the multiple-state receptor helped), with a percentile
bootstrap over cases. Strata (main cohort, backbone-change) are reported separately and never pooled.
Cases whose multiple arm has one receptor stay in as zero-difference controls; a subset restricted to cases
that actually have multiple states is reported alongside.
"""
import csv, json, sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
SUCCESS, BOOT, SEED = 2.0, 10000, 20261028


def scored_poses(pdbqt, ref):
    """[(vina score, RMSD to ref)] for every pose in one Vina output; RMSD is in place, symmetry aware."""
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    rd = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0])
    scores = [float(pose.score) for pose in pm]
    assert rd.GetNumConformers() == len(scores), f'{pdbqt}: {rd.GetNumConformers()} poses, {len(scores)} scores'
    return [(s, rdMolAlign.CalcRMS(rd, ref, prbId=i)) for i, s in enumerate(scores)]


def arm_result(case_dir, arm, ref):
    log = json.loads((case_dir/arm/'docking_log.json').read_text())
    scored = []
    for call in log['calls']:
        scored += scored_poses(case_dir/arm/call['out'], ref)
    assert scored, f'no poses for {case_dir.name} {arm}'
    top = min(scored, key=lambda x: x[0])                 # best Vina score over all calls
    rmsds = [r for _, r in scored]
    return dict(poses=len(scored), receptors=len(log['receptors']), best_score=round(top[0], 3),
                rmsd_top=round(top[1], 3), rmsd_oracle=round(min(rmsds), 3),
                success_top=int(top[1] <= SUCCESS), success_any=int(min(rmsds) <= SUCCESS))


def boot(values, rng):
    v = np.asarray(values, float)
    if len(v) < 2:
        return [float('nan'), float('nan')]
    draws = v[rng.integers(0, len(v), size=(BOOT, len(v)))].mean(1)
    return [round(float(np.percentile(draws, 2.5)), 4), round(round(float(np.percentile(draws, 97.5)), 4), 4)]


def main():
    reps = {r['case_id']: r for r in csv.DictReader(open(PILOT/'design2x2_v1/representatives.tsv'), delimiter='\t')}
    cases = [c['case_id'] for c in json.loads((HERE/'inputs_manifest.json').read_text())['cases']]
    rows = []
    for cid in cases:
        d, inp = HERE/'docked'/cid, HERE/'inputs'/cid
        refs = {s: Chem.RemoveHs(Chem.MolFromMolFile(str(inp/f'ligand_ref_{s}.sdf'), removeHs=False)) for s in ('holo', 'apo')}
        row = dict(case_id=cid, stratum=reps[cid]['stratum'], apo=reps[cid]['apo'], holo=reps[cid]['holo'],
                   ligand=reps[cid]['ligand'], ligand_name=reps[cid]['ligand_name'][:30])
        for arm in ARMS:
            for k, v in arm_result(d, arm, refs[arm.split('_')[0]]).items():
                row[f'{arm}_{k}'] = v
        for side in ('holo', 'apo'):
            row[f'{side}_gain_top'] = round(row[f'{side}_single_rmsd_top'] - row[f'{side}_multi_rmsd_top'], 3)
            row[f'{side}_gain_oracle'] = round(row[f'{side}_single_rmsd_oracle'] - row[f'{side}_multi_rmsd_oracle'], 3)
        row['has_multiple_states'] = int(row['apo_multi_receptors'] > 1 or row['holo_multi_receptors'] > 1)
        rows.append(row)
        print(row['case_id'][:8], row['stratum'][:4], 'apo single/multi %.2f/%.2f' %
              (row['apo_single_rmsd_top'], row['apo_multi_rmsd_top']), flush=True)

    with open(HERE/'docking_results.tsv', 'w', newline='') as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)

    rng = np.random.default_rng(SEED)
    summary = {}
    for stratum in ['main', 'backbone_change']:
        sel = [r for r in rows if r['stratum'] == stratum]
        if not sel:
            continue
        s = dict(cases=len(sel), with_multiple_states=sum(r['has_multiple_states'] for r in sel))
        for arm in ARMS:
            s[arm] = dict(mean_rmsd_top=round(float(np.mean([r[f'{arm}_rmsd_top'] for r in sel])), 3),
                          median_rmsd_top=round(float(np.median([r[f'{arm}_rmsd_top'] for r in sel])), 3),
                          success_top=sum(r[f'{arm}_success_top'] for r in sel),
                          mean_rmsd_oracle=round(float(np.mean([r[f'{arm}_rmsd_oracle'] for r in sel])), 3))
        for name, vals in [('apo_multi_minus_single', [r['apo_gain_top'] for r in sel]),
                           ('holo_multi_minus_single', [r['holo_gain_top'] for r in sel]),
                           ('apo_minus_holo_single', [r['apo_single_rmsd_top'] - r['holo_single_rmsd_top'] for r in sel]),
                           ('interaction', [r['apo_gain_top'] - r['holo_gain_top'] for r in sel])]:
            s[name] = dict(mean=round(float(np.mean(vals)), 4), ci95=boot(vals, rng), n=len(vals))
        sub = [r for r in sel if r['has_multiple_states']]
        if sub:
            s['apo_multi_minus_single_states_only'] = dict(
                mean=round(float(np.mean([r['apo_gain_top'] for r in sub])), 4),
                ci95=boot([r['apo_gain_top'] for r in sub], rng), n=len(sub))
        summary[stratum] = s
    summary['note'] = ('Positive gain means the multiple-state receptor gave a lower RMSD. Primary: '
                       'apo_multi_minus_single in stratum main, decided by the CI lower bound above zero.')
    (HERE/'docking_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
