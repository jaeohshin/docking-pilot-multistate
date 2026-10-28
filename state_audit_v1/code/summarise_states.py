"""Join the selection calls, the whole-pocket refit and the qFit pocket states into one table, and count cases and
components usable for the 2x2 design (holo/apo backbone x single/multiple side-chain states).

Whole-pocket check: the fixed input rule was pocket backbone RMSD <= 0.75 A, measured on the receptor chain only.
Re-measured over every pocket chain (refit_pose.py); eligible cases above 0.75 A are moved to a backbone-change
stratum (not excluded). Multi-state: a pocket residue with >= 2 occupied wells (chi1, or joint chi1/chi2).
holo arm informative = holo pocket has >= 1 multi-state residue; apo arm likewise; a case without any is still a
valid zero-difference control for that backbone.
"""
import csv, json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
Q = {o['case_id']: o for o in json.loads((HERE/'qfit_state_audit.json').read_text())['cases']}
F = {o['case_id']: o for o in json.loads((HERE/'pose_refit.json').read_text())['cases']}
REV = {r['case_id']: r for r in csv.DictReader(open(PILOT/'final_selection_v1/case_review_all79.tsv'), delimiter='\t')}
rows = []
for cid, q in Q.items():
    f, r = F[cid], REV[cid]
    h, a = q['holo_states'], f['apo_states']
    status = r['status']
    if status == 'accept' and f['refit_rmsd'] > 0.75:
        status = 'backbone_change'
    rows.append(dict(case_id=cid, homology_component=q['homology_component'], apo=q['apo'], holo=q['holo'], ligand=q['ligand'],
                     selection_status=r['status'], status_2x2=status, pilot_representative=int(r['pilot_representative']),
                     whole_pocket_rmsd=round(f['refit_rmsd'], 3), receptor_chain_rmsd=r['pocket_backbone_rmsd_deposited'],
                     fitted_residues=f['refit_residues'], ligand_shift_vs_first_transfer=round(f['ligand_shift_max'], 3),
                     holo_pocket=h['pocket_residues'], holo_multi_chi1=h['multi_chi1'], holo_multi_chi12=h['multi_chi12'],
                     apo_pocket=a['pocket_residues'], apo_multi_chi1=a['multi_chi1'], apo_multi_chi12=a['multi_chi12'],
                     holo_multi_residues=';'.join(h['multi_chi12_residues']), apo_multi_residues=';'.join(a['multi_chi12_residues']),
                     holo_max_backbone_altloc_spread=h['max_backbone_spread'], apo_max_backbone_altloc_spread=a['max_backbone_spread'],
                     holo_chi12_combinations=h['chi12_combinations'], apo_chi12_combinations=a['chi12_combinations'],
                     apo_crystal_dependence_refit=f['apo_contacts']['crystal_dependence'],
                     minor_state_residues=r['minor_state_residues']))
rows.sort(key=lambda x: (x['homology_component'], x['status_2x2'], -x['pilot_representative'], x['case_id']))
with open(HERE/'state_summary.tsv', 'w', newline='') as h:
    w = csv.DictWriter(h, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)


def count(sel, key):
    hs = [r for r in sel if r[f'holo_multi_{key}']]
    ap = [r for r in sel if r[f'apo_multi_{key}']]
    both = [r for r in sel if r[f'holo_multi_{key}'] and r[f'apo_multi_{key}']]
    comp = lambda s: len({r['homology_component'] for r in s})
    return dict(cases=len(sel), components=comp(sel), holo_multi=[len(hs), comp(hs)], apo_multi=[len(ap), comp(ap)],
                both=[len(both), comp(both)], either=[len({r['case_id'] for r in hs + ap}), comp(hs + ap)])


ok = [r for r in rows if r['status_2x2'] == 'accept']
reps = [r for r in ok if r['pilot_representative']]
summary = dict(status_2x2=Counter(r['status_2x2'] for r in rows),
               backbone_change_cases=[(r['case_id'], r['holo'], r['ligand'], r['whole_pocket_rmsd']) for r in rows if r['status_2x2'] == 'backbone_change'],
               eligible_chi1=count(ok, 'chi1'), eligible_chi12=count(ok, 'chi12'),
               representatives_still_eligible=len(reps), representatives_chi1=count(reps, 'chi1'), representatives_chi12=count(reps, 'chi12'),
               components_lost_to_backbone_change=sorted({r['homology_component'] for r in rows if r['status_2x2'] == 'backbone_change'}
                                                         - {r['homology_component'] for r in ok}),
               backbone_altloc_spread_gt_0p5=dict(holo=sum(r['holo_max_backbone_altloc_spread'] > 0.5 for r in ok),
                                                  apo=sum(r['apo_max_backbone_altloc_spread'] > 0.5 for r in ok), of=len(ok)))
(HERE/'state_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
print(json.dumps(summary, indent=1))
for r in ok:
    print('R' if r['pilot_representative'] else ' ', r['homology_component'][1:7], r['case_id'][:8], r['holo'], r['ligand'], 'rmsd', r['whole_pocket_rmsd'],
          'H', r['holo_multi_chi1'], r['holo_multi_chi12'], 'A', r['apo_multi_chi1'], r['apo_multi_chi12'])
