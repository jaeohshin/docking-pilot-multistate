"""확대 cohort 분석. ANALYSIS_PLAN.md(도킹 결과 생성 전 고정)대로만 집계한다.

1차와 다른 점은 하나뿐이다: bootstrap 단위가 사례가 아니라 **상동성 component**다.
같은 단백질의 여러 사례는 독립이 아니므로 사례 단위로 재표집하면 CI가 거짓으로 좁아진다.
component를 복원 추출하고, 뽑힌 component의 사례 평균을 먼저 낸 뒤 그 평균들의 평균을 통계량으로 쓴다.
"""
import csv, json, sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
HERE = Path(__file__).resolve().parents[1]
ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
SUCCESS, BOOT, SEED = 2.0, 10000, 20261028


def scored(pdbqt, ref):
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    rd = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0])
    s = [float(p.score) for p in pm]
    assert rd.GetNumConformers() == len(s), pdbqt
    return [(x, rdMolAlign.CalcRMS(rd, ref, prbId=i)) for i, x in enumerate(s)]


def arm_result(d, arm, ref):
    log = json.loads((d/arm/'docking_log.json').read_text())
    out = []
    for call in log['calls']:
        if call['returncode'] == 0:
            out += scored(d/arm/call['out'], ref)
    if not out:
        return None
    top = min(out, key=lambda x: x[0])
    r = [x for _, x in out]
    return dict(poses=len(out), receptors=len(log['receptors']), best_score=round(top[0], 3),
                rmsd_top=round(top[1], 3), rmsd_oracle=round(min(r), 3),
                success_top=int(top[1] <= SUCCESS), success_any=int(min(r) <= SUCCESS))


def cluster_boot(pairs, rng):
    """pairs: [(component, value)]. component를 복원 추출하고 component 평균들의 평균을 낸다."""
    comps = sorted({c for c, _ in pairs})
    by = {c: np.array([v for cc, v in pairs if cc == c], float) for c in comps}
    means = np.array([by[c].mean() for c in comps])
    if len(comps) < 2:
        return [float('nan')] * 2
    draws = means[rng.integers(0, len(comps), size=(BOOT, len(comps)))].mean(1)
    return [round(float(np.percentile(draws, 2.5)), 4), round(float(np.percentile(draws, 97.5)), 4)]


def main():
    comp = {c['case_id']: c['homology_component'] for c in json.loads((HERE/'candidates.json').read_text())}
    sel = {r['case_id']: r for r in csv.DictReader(open(HERE/'cohort_final.tsv'), delimiter='\t')}
    rows = []
    for d in sorted((HERE/'docked').iterdir()):
        cid = d.name
        if cid not in sel or not all((d/a/'docking_log.json').exists() for a in ARMS):
            continue
        inp = HERE/'inputs'/cid
        refs = {s: Chem.RemoveHs(Chem.MolFromMolFile(str(inp/f'ligand_ref_{s}.sdf'), removeHs=False)) for s in ('holo', 'apo')}
        if any(v is None for v in refs.values()):
            continue
        res = {a: arm_result(d, a, refs[a.split('_')[0]]) for a in ARMS}
        if any(v is None for v in res.values()):
            continue
        r = dict(case_id=cid, component=comp[cid], stratum=sel[cid]['status'], apo=sel[cid]['apo'],
                 holo=sel[cid]['holo'], ligand=sel[cid]['ligand'])
        for a in ARMS:
            for k, v in res[a].items():
                r[f'{a}_{k}'] = v
        for side in ('holo', 'apo'):
            r[f'{side}_gain_top'] = round(r[f'{side}_single_rmsd_top'] - r[f'{side}_multi_rmsd_top'], 3)
            r[f'{side}_gain_oracle'] = round(r[f'{side}_single_rmsd_oracle'] - r[f'{side}_multi_rmsd_oracle'], 3)
        r['has_multiple_states'] = int(r['apo_multi_receptors'] > 1 or r['holo_multi_receptors'] > 1)
        rows.append(r)
        if len(rows) % 20 == 0:
            print(' scored', len(rows), flush=True)
    assert rows, 'no complete cases'
    with open(HERE/'docking_results.tsv', 'w', newline='') as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)

    rng = np.random.default_rng(SEED)
    summary = {}
    for stratum in ['accept', 'backbone_change']:
        s = [r for r in rows if r['stratum'] == stratum]
        if not s:
            continue
        o = dict(cases=len(s), components=len({r['component'] for r in s}),
                 with_multiple_states=sum(r['has_multiple_states'] for r in s))
        for a in ARMS:
            o[a] = dict(mean_rmsd_top=round(float(np.mean([r[f'{a}_rmsd_top'] for r in s])), 3),
                        median_rmsd_top=round(float(np.median([r[f'{a}_rmsd_top'] for r in s])), 3),
                        success_top=sum(r[f'{a}_success_top'] for r in s),
                        mean_rmsd_oracle=round(float(np.mean([r[f'{a}_rmsd_oracle'] for r in s])), 3))
        for name, key in [('apo_multi_minus_single', 'apo_gain_top'), ('holo_multi_minus_single', 'holo_gain_top'),
                          ('apo_oracle', 'apo_gain_oracle'), ('holo_oracle', 'holo_gain_oracle')]:
            pairs = [(r['component'], r[key]) for r in s]
            cm = [np.mean([v for c, v in pairs if c == comp_]) for comp_ in sorted({c for c, _ in pairs})]
            o[name] = dict(mean_of_component_means=round(float(np.mean(cm)), 4),
                           ci95=cluster_boot(pairs, rng), components=len(cm), cases=len(pairs))
        sub = [r for r in s if r['has_multiple_states']]
        if sub:
            p = [(r['component'], r['apo_gain_top']) for r in sub]
            cm = [np.mean([v for c, v in p if c == q]) for q in sorted({c for c, _ in p})]
            o['apo_multi_minus_single_states_only'] = dict(
                mean_of_component_means=round(float(np.mean(cm)), 4), ci95=cluster_boot(p, rng),
                components=len(cm), cases=len(p))
        summary[stratum] = o
    summary['note'] = ('Positive gain means the multiple-state receptor gave a lower RMSD. Bootstrap unit is the '
                       'homology component, not the case. Primary: apo_multi_minus_single in stratum accept, '
                       'decided by the CI lower bound above zero.')
    (HERE/'docking_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
