"""Representatives for the 2x2 design (user decision 2026-10-28): one per homology component, separately for the main
cohort (whole-pocket RMSD <= 0.75 A) and the backbone_change stratum (> 0.75 A).

Order within a component: (not [holo and apo pockets each have >= 1 chi1/chi2 multi-state residue],
then the final_selection_v1 order: pose not clear [RSCC >= 0.90 and mean occupancy >= 0.90], environment tier
max(holo, apo crystal dependence) with the apo side from the whole-pocket refit, worse resolution, case_id).
Multi-state presence is the treatment variable, not an outcome; docking results and the apo-minor/holo-major signal
are not used. Components without a qualifying case keep their best case as a zero-difference control.
"""
import csv, json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
S = {r['case_id']: r for r in csv.DictReader(open(PILOT/'state_audit_v1/state_summary.tsv'), delimiter='\t')}
V = {r['case_id']: r for r in csv.DictReader(open(PILOT/'final_selection_v1/case_review_all79.tsv'), delimiter='\t')}
TIER = {'none': 0, 'minor': 1, 'major': 2}


def key(cid):
    s, v = S[cid], V[cid]
    treated = int(s['holo_multi_chi12']) >= 1 and int(s['apo_multi_chi12']) >= 1
    rscc = float(v['rscc']) if v['rscc'] not in ('', 'None') else None
    clear = rscc is not None and rscc >= 0.90 and float(v['ligand_mean_occupancy']) >= 0.90
    tier = max(TIER[v['holo_crystal_dependence']], TIER[s['apo_crystal_dependence_refit']])
    return (not treated, not clear, tier, max(float(v['res_apo']), float(v['res_holo'])), cid)


out = []
for stratum, status in (('main', 'accept'), ('backbone_change', 'backbone_change')):
    comp = defaultdict(list)
    for cid, s in S.items():
        if s['status_2x2'] == status:
            comp[s['homology_component']].append(cid)
    for k, ids in sorted(comp.items()):
        best = min(ids, key=key)
        s, v = S[best], V[best]
        out.append(dict(stratum=stratum, homology_component=k, case_id=best, apo=s['apo'], holo=s['holo'], ligand=s['ligand'],
                        ligand_name=v['ligand_name'], ligand_role=v['ligand_role'], treated=int(not key(best)[0]),
                        holo_multi_chi12=s['holo_multi_chi12'], apo_multi_chi12=s['apo_multi_chi12'],
                        whole_pocket_rmsd=s['whole_pocket_rmsd'], rscc=v['rscc'], res_apo=v['res_apo'], res_holo=v['res_holo'],
                        holo_crystal_dependence=v['holo_crystal_dependence'], apo_crystal_dependence=s['apo_crystal_dependence_refit'],
                        previous_representative=int(v['pilot_representative']), candidates_in_component=len(ids)))
with open(HERE/'representatives.tsv', 'w', newline='') as h:
    w = csv.DictWriter(h, list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)
for r in out:
    print(r['stratum'][:4], r['homology_component'][1:9], r['case_id'], r['holo'], r['ligand'], 'treated', r['treated'],
          'H/A', r['holo_multi_chi12'], r['apo_multi_chi12'], 'prev rep', r['previous_representative'], 'n', r['candidates_in_component'])
