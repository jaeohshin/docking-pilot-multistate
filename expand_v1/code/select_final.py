"""ANALYSIS_PLAN.md의 포함·제외 기준을 확대 cohort에 적용한다. 전부 구조에서 나오는 값이며 도킹 결과와 무관하다."""
import json, sys, csv
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
sys.path.insert(0, str(HERE/'code'))
from ligand_roles_auto import role

C = {c['case_id']: c for c in json.loads((HERE/'candidates.json').read_text())}
A = {o['case_id']: o for o in json.loads((HERE/'assembly_contact_audit.json').read_text())['cases']}
F = {o['case_id']: o for o in json.loads((HERE/'pose_refit.json').read_text())['cases']}
S = {o['case_id']: o for o in json.loads((HERE/'qfit_state_audit.json').read_text())['cases']}

rows = []
for cid, c in C.items():
    o, f = A.get(cid), F.get(cid)
    if o is None:
        rows.append(dict(case_id=cid, status='exclude', reason='assembly audit failed')); continue
    if f is None:
        rows.append(dict(case_id=cid, status='exclude', reason='apo pose refit failed')); continue
    h = o['holo_contacts'][o['holo_assembly_used']]
    dep = [d for d in o['holo_neighbour_dependence'] if not d['present_at_apo_site']]
    miss = [p for p in o['holo_partner_chains'] if not p['apo_in_assembly_same_sequence_copies']]
    rscc = float(o['validation']['rscc']) if o['validation'].get('rscc') else None
    whole = f['refit_rmsd']
    if role(c['ligand_id']) == 'exclude':
        st, why = 'exclude', 'buffer or additive'
    elif h['crystal_dependence'] == 'major':
        st, why = 'exclude', f"crystal-context dependent ({h['outside_share']:.2f} of contact atoms)"
    elif dep:
        st, why = 'exclude', 'cofactor or co-ligand absent at the apo site'
    elif miss:
        st, why = 'exclude', 'partner chain absent at the apo site'
    elif rscc is not None and rscc < 0.80:
        st, why = 'exclude', f'reference pose unclear (RSCC {rscc:.2f})'
    elif whole > 0.75:
        st, why = 'backbone_change', f'whole-pocket RMSD {whole:.2f} A'
    else:
        st, why = 'accept', 'eligible'
    hs, as_ = S[cid]['holo_states'], S[cid]['apo_states']
    rows.append(dict(case_id=cid, status=st, reason=why, apo=c['apo'], holo=c['holo'], ligand=c['ligand_id'],
                     ligand_name=c['ccd']['name'][:40] if c.get('ccd') else '', source_group=c['group'],
                     res_apo=c['res_apo'], res_holo=c['res_holo'], rscc=rscc,
                     whole_pocket_rmsd=round(whole, 3), ligand_heavy_atoms=c['ligand_heavy_atoms'],
                     holo_multi_chi12=hs['multi_chi12'], apo_multi_chi12=as_['multi_chi12'],
                     holo_crystal_dependence=h['crystal_dependence'],
                     apo_crystal_dependence=f['apo_contacts']['crystal_dependence']))

with open(HERE/'cohort_final.tsv', 'w', newline='') as h:
    w = csv.DictWriter(h, sorted({k for r in rows for k in r}), delimiter='\t'); w.writeheader(); w.writerows(rows)
import collections
by = collections.Counter(r['status'] for r in rows)
acc = [r for r in rows if r['status'] == 'accept']
summary = dict(total=len(rows), status=dict(by),
               accept_cases=len(acc), accept_source_groups=len({r['source_group'] for r in acc}),
               accept_with_multiple_states=sum(1 for r in acc if r['holo_multi_chi12'] or r['apo_multi_chi12']),
               exclusion_reasons=dict(collections.Counter(r['reason'] for r in rows if r['status'] == 'exclude')))
(HERE/'cohort_final.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
