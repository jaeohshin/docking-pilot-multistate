"""Experiment C analysis, fixed by gnina_dock/PLAN.md before any gnina-search pose was scored.

Everything that defines the comparison is inherited from the frozen analyses: symmetry-aware in-place
heavy-atom RMSD against that backbone's own reference ligand, the paired per-case contrast
(single RMSD minus multi RMSD, so positive means the multi-state receptor helped), the cluster bootstrap
over homology components, success at 2 A, and the top-5-distinct rule with its 1 A separation.

Two things differ from the Vina analyses, both because gnina writes SDF with its own scores:
  - poses are read from SDF, not pdbqt, so Meeko is not involved on the pose side;
  - the selector is CNNscore, gnina's pose-ranking quantity (PLAN.md). CNNaffinity and minimizedAffinity
    are carried through as alternative selectors and reported beside it, never chosen after the fact.

C is a registered exploratory experiment, not a third confirmatory test: the confirmatory family stays
P1 and P2 of budget50. Its purpose is to remove the sampling-scoring mismatch of the rescoring
experiment, and it is read for that purpose only.
"""
import csv, gzip, json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog('rdApp.*')
HERE = Path(__file__).resolve().parents[1]
ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
SUCCESS, BOOT, SEED = 2.0, 10000, 20261028
TOPN, DISTINCT, CURVE = 5, 1.0, 10
# every selector is written so that LOWER IS BETTER, matching the Vina analyses
SELECTORS = {'cnnscore': ('CNNscore', -1.0), 'cnnaffinity': ('CNNaffinity', -1.0),
             'minaffinity': ('minimizedAffinity', 1.0)}
PRIMARY = 'cnnscore'


def read_sdf(path, ref):
    """[(scores dict, rmsd, conf_id, mol)] for one gnina output file, in file order."""
    out, mols = [], []
    for m in Chem.SDMolSupplier(str(path), removeHs=False):
        if m is None:
            continue
        # RemoveAllHs, not the supplier's removeHs: the supplier keeps a hydrogen that defines
        # double-bond stereochemistry, which left one amidine ligand with 36 heavy atoms against the
        # reference's 35 and made CalcRMS fail to match. Checked on the 1440 poses that already
        # matched: every RMSD is identical to the supplier-stripped value, so this only adds the
        # case that could not be read before.
        mols.append(Chem.RemoveAllHs(m))
    if not mols:
        return []
    # gnina writes each pose as its own record; merge them into one mol so CalcRMS can use conformer ids
    base = Chem.Mol(mols[0])
    base.RemoveAllConformers()
    for m in mols:
        base.AddConformer(m.GetConformer(), assignId=True)
    for i, m in enumerate(mols):
        sc = {k: float(m.GetProp(p)) * s for k, (p, s) in SELECTORS.items() if m.HasProp(p)}
        out.append((sc, rdMolAlign.CalcRMS(base, ref, prbId=i), i, base))
    return out


def pick_distinct(order, n, thresh):
    kept = []
    for cand in order:
        if all(rdMolAlign.CalcRMS(cand[3], k[3], prbId=cand[2], refId=k[2]) > thresh for k in kept):
            kept.append(cand)
            if len(kept) == n:
                break
    return kept


def arm_result(d, arm, ref, sink, cid):
    log = json.loads((d/arm/'docking_log.json').read_text())
    out, per_call = [], []
    for call in log['calls']:
        if call['returncode']:
            continue
        ps = read_sdf(d/arm/call['out'], ref)
        if not ps:
            continue
        per_call.append(min(p[1] for p in ps))
        for sc, rm, i, _ in ps:
            sink.append((cid, arm, call['run'], i, round(rm, 3)) +
                        tuple(round(sc.get(k, float('nan')), 4) for k in SELECTORS))
        out += ps
    if not out:
        return None
    r = np.array([x[1] for x in out])
    res = dict(poses=len(out), calls=len(per_call), receptors=len(log['receptors']),
               rmsd_oracle=round(float(r.min()), 3), success_any=int(r.min() <= SUCCESS),
               frac_good=round(float((r <= SUCCESS).mean()), 4),
               calls_good=sum(1 for x in per_call if x <= SUCCESS))
    for name in SELECTORS:
        usable = [x for x in out if name in x[0]]
        if not usable:
            continue
        order = sorted(usable, key=lambda x: x[0][name])
        dist = pick_distinct(order, TOPN, DISTINCT)
        res[f'{name}_rmsd_top'] = round(order[0][1], 3)
        res[f'{name}_rmsd_top5'] = round(float(min(x[1] for x in dist)), 3)
        res[f'{name}_rmsd_top5_raw'] = round(float(min(x[1] for x in order[:TOPN])), 3)
        res[f'{name}_n_distinct'] = len(dist)
        res[f'{name}_success_top'] = int(order[0][1] <= SUCCESS)
        res[f'{name}_success_top5'] = int(min(x[1] for x in dist) <= SUCCESS)
        res[f'{name}_rank_best'] = int(np.argmin([x[1] for x in order])) + 1
        for n in range(1, CURVE + 1):
            res[f'{name}_top{n}'] = int(min(x[1] for x in order[:n]) <= SUCCESS)
    return res


def cluster_boot(pairs, rng, lo=2.5, hi=97.5):
    comps = sorted({c for c, _ in pairs})
    by = {c: np.array([v for cc, v in pairs if cc == c], float) for c in comps}
    means = np.array([by[c].mean() for c in comps])
    if len(comps) < 2:
        return [float('nan')] * 2
    draws = means[rng.integers(0, len(comps), size=(BOOT, len(comps)))].mean(1)
    return [round(float(np.percentile(draws, lo)), 4), round(float(np.percentile(draws, hi)), 4)]


def contrast(sel, key, rng):
    pairs = [(r['component'], r[key]) for r in sel if key in r]
    if not pairs:
        return None
    cm = [np.mean([v for c, v in pairs if c == q]) for q in sorted({c for c, _ in pairs})]
    return dict(mean_of_component_means=round(float(np.mean(cm)), 4), ci95=cluster_boot(pairs, rng),
                components=len(cm), cases=len(pairs))


def main():
    comp = {c['case_id']: c['homology_component'] for c in json.loads((HERE/'candidates.json').read_text())}
    sel = {r['case_id']: r for r in csv.DictReader(open(HERE/'cohort_final.tsv'), delimiter='\t')}
    rows, sink = [], []
    for d in sorted((HERE/'docked').iterdir()):
        cid = d.name
        if cid not in sel or not all((d/x/'docking_log.json').exists() for x in ARMS):
            continue
        inp = HERE/'inputs'/cid
        refs = {s: Chem.RemoveHs(Chem.MolFromMolFile(str(inp/f'ligand_ref_{s}.sdf'), removeHs=False)) for s in ('holo', 'apo')}
        if any(v is None for v in refs.values()):
            continue
        mine = []
        res = {x: arm_result(d, x, refs[x.split('_')[0]], mine, cid) for x in ARMS}
        if any(v is None for v in res.values()):
            continue
        sink += mine
        r = dict(case_id=cid, component=comp[cid], stratum=sel[cid]['status'], apo=sel[cid]['apo'],
                 holo=sel[cid]['holo'], ligand=sel[cid]['ligand'])
        for x in ARMS:
            for k, v in res[x].items():
                r[f'{x}_{k}'] = v
        for side in ('holo', 'apo'):
            r[f'{side}_gain_oracle'] = round(r[f'{side}_single_rmsd_oracle'] - r[f'{side}_multi_rmsd_oracle'], 3)
            for name in SELECTORS:
                for lbl in ('top', 'top5', 'top5_raw'):
                    a, b = f'{side}_single_{name}_rmsd_{lbl}', f'{side}_multi_{name}_rmsd_{lbl}'
                    if a in r and b in r:
                        r[f'{side}_{name}_gain_{lbl}'] = round(r[a] - r[b], 3)
        r['has_multiple_states'] = int(r['apo_multi_receptors'] > 1 or r['holo_multi_receptors'] > 1)
        rows.append(r)
        if len(rows) % 20 == 0:
            print(' scored', len(rows), flush=True)
    assert rows, 'no complete cases'
    with open(HERE/'docking_results.tsv', 'w', newline='') as h:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(h, keys, delimiter='\t', restval=''); w.writeheader(); w.writerows(rows)
    with gzip.open(HERE/'pose_table.tsv.gz', 'wt', newline='') as h:
        w = csv.writer(h, delimiter='\t')
        w.writerow(['case_id', 'arm', 'call', 'pose', 'rmsd'] + list(SELECTORS))
        w.writerows(sink)
    print('poses written', len(sink))

    rng = np.random.default_rng(SEED)
    summary = {}
    for stratum in ['accept', 'backbone_change']:
        s = [r for r in rows if r['stratum'] == stratum]
        if not s:
            continue
        o = dict(cases=len(s), components=len({r['component'] for r in s}),
                 with_multiple_states=sum(r['has_multiple_states'] for r in s),
                 primary_selector=PRIMARY)
        for x in ARMS:
            e = dict(mean_poses=round(float(np.mean([r[f'{x}_poses'] for r in s])), 1),
                     success_any=sum(r[f'{x}_success_any'] for r in s),
                     mean_rmsd_oracle=round(float(np.mean([r[f'{x}_rmsd_oracle'] for r in s])), 3),
                     mean_frac_good=round(float(np.mean([r[f'{x}_frac_good'] for r in s])), 4))
            for name in SELECTORS:
                k = f'{x}_{name}_rmsd_top'
                if k not in s[0]:
                    continue
                e[name] = dict(success_top=sum(r[f'{x}_{name}_success_top'] for r in s),
                               success_top5=sum(r[f'{x}_{name}_success_top5'] for r in s),
                               mean_rmsd_top=round(float(np.mean([r[k] for r in s])), 3),
                               median_rmsd_top=round(float(np.median([r[k] for r in s])), 3),
                               median_rank_best=float(np.median([r[f'{x}_{name}_rank_best'] for r in s])),
                               success_curve={f'top{n}': sum(r[f'{x}_{name}_top{n}'] for r in s) for n in range(1, CURVE + 1)})
            o[x] = e
        for side in ('holo', 'apo'):
            o[f'{side}_oracle'] = contrast(s, f'{side}_gain_oracle', rng)
            for name in SELECTORS:
                for lbl in ('top', 'top5', 'top5_raw'):
                    c = contrast(s, f'{side}_{name}_gain_{lbl}', rng)
                    if c:
                        o[f'{side}_multi_minus_single_{name}_{lbl}'] = c
        summary[stratum] = o
    summary['note'] = (
        'Experiment C: gnina both searches and selects, so there is no sampling-scoring mismatch. '
        'Positive gain means the multiple-state receptor did better. Bootstrap unit is the homology '
        'component. The primary selector is CNNscore, gnina\'s pose-ranking quantity, fixed before any '
        'pose was scored; CNNaffinity and minimizedAffinity are alternatives reported beside it. '
        'EXPLORATORY: the confirmatory family is budget50 P1 and P2, and a positive result here does not '
        'change that. gnina may have our holo complexes in its training set, and letting it steer the '
        'search can amplify that, so absolute success rates are assumed optimistic and only the paired '
        'contrast is used. The pose set differs from the Vina runs, so oracle and success counts are not '
        'compared across experiments.')
    (HERE/'docking_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
