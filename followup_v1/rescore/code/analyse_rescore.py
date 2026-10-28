"""Judge each scoring function on the same fixed pose set, the way followup_v1/PLAN.md fixed it.

Only the selection step differs between scoring functions. The poses, and the RMSD of every pose, come from
one pose table that was written once (pose_table.py), so no geometry is recomputed here and nothing about the
comparison can drift between functions. For each (case, arm) and each scoring function we re-rank that arm's
own poses by that function's score, take the new rank-1 pose, and read its RMSD.

Every number here is EXPLORATORY by the plan's multiplicity rule (§다중 검정 4): the single confirmatory
family lives in budget50. A positive result here is reported as one cell of the full table, never as "the"
result, and the absolute top-1 success rate of each function is printed beside every contrast -- a function
that cannot pick poses at all says nothing about the receptor question.

Top-5 is reported without the distinctness filter that budget50 applies, because that filter needs the pose
geometry and this script deliberately holds only scores and RMSDs. It is labelled top5_raw throughout.
"""
import argparse, csv, gzip, json
from collections import defaultdict
from pathlib import Path

import numpy as np

ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
SUCCESS, BOOT, SEED, TOPN = 2.0, 10000, 20261028, 5


def read_pose_table(path):
    """{(case, arm): [(call, pose, vina_score, rmsd)]}"""
    t = defaultdict(list)
    with gzip.open(path, 'rt') as h:
        for r in csv.DictReader(h, delimiter='\t'):
            t[(r['case_id'], r['arm'])].append((int(r['call']), int(r['pose']),
                                                float(r['vina_score']), float(r['rmsd'])))
    return t


def read_scores(d):
    """{(case, arm): {(call, pose): score}} from one scoring function's output directory."""
    s = {}
    for f in sorted(Path(d).glob('*__*.tsv')):
        cid, arm = f.stem.split('__')
        with open(f) as h:
            s[(cid, arm)] = {(int(r['call']), int(r['pose'])): float(r['score']) for r in csv.DictReader(h, delimiter='\t')}
    return s


def cluster_boot(pairs, rng, lo=2.5, hi=97.5):
    comps = sorted({c for c, _ in pairs})
    by = {c: np.array([v for cc, v in pairs if cc == c], float) for c in comps}
    means = np.array([by[c].mean() for c in comps])
    if len(comps) < 2:
        return [float('nan')] * 2
    draws = means[rng.integers(0, len(comps), size=(BOOT, len(comps)))].mean(1)
    return [round(float(np.percentile(draws, lo)), 4), round(float(np.percentile(draws, hi)), 4)]


def contrast(pairs, rng):
    cm = [np.mean([v for c, v in pairs if c == q]) for q in sorted({c for c, _ in pairs})]
    return dict(mean_of_component_means=round(float(np.mean(cm)), 4), ci95=cluster_boot(pairs, rng),
                components=len(cm), cases=len(pairs))


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--run', required=True, help='the docking run the poses came from')
    a.add_argument('--pose-table', required=True)
    a.add_argument('--scores', nargs='+', required=True, help='name=dir, e.g. vinardo=out/vinardo')
    a.add_argument('--out', required=True)
    a = a.parse_args()
    run = Path(a.run)
    comp = {c['case_id']: c['homology_component'] for c in json.loads((run/'candidates.json').read_text())}
    sel = {r['case_id']: r for r in csv.DictReader(open(run/'cohort_final.tsv'), delimiter='\t')}
    table = read_pose_table(a.pose_table)
    funcs = {'vina': None}
    for spec in a.scores:
        name, d = spec.split('=', 1)
        funcs[name] = read_scores(d)

    rows = []
    for cid in sorted({c for c, _ in table}):
        if cid not in sel:
            continue
        if not all((cid, arm) in table for arm in ARMS):
            continue
        r = dict(case_id=cid, component=comp[cid], stratum=sel[cid]['status'])
        ok = True
        for name, sc in funcs.items():
            for arm in ARMS:
                poses = table[(cid, arm)]
                if name == 'vina':
                    key = {(c, p): v for c, p, v, _ in poses}
                else:
                    key = sc.get((cid, arm))
                    if not key or len(key) != len(poses):
                        ok = False
                        break
                order = sorted(poses, key=lambda x: key[(x[0], x[1])])
                r[f'{name}_{arm}_rmsd_top'] = round(order[0][3], 3)
                r[f'{name}_{arm}_rmsd_top5_raw'] = round(min(x[3] for x in order[:TOPN]), 3)
                r[f'{name}_{arm}_success_top'] = int(order[0][3] <= SUCCESS)
                r[f'{name}_{arm}_success_top5_raw'] = int(min(x[3] for x in order[:TOPN]) <= SUCCESS)
            if not ok:
                break
            for side in ('holo', 'apo'):
                r[f'{name}_{side}_gain_top'] = round(r[f'{name}_{side}_single_rmsd_top'] - r[f'{name}_{side}_multi_rmsd_top'], 3)
                r[f'{name}_{side}_gain_top5_raw'] = round(r[f'{name}_{side}_single_rmsd_top5_raw'] - r[f'{name}_{side}_multi_rmsd_top5_raw'], 3)
        # the oracle is a property of the pose set, so it is identical for every scoring function
        for arm in ARMS:
            r[f'{arm}_rmsd_oracle'] = round(min(x[3] for x in table[(cid, arm)]), 3)
            r[f'{arm}_poses'] = len(table[(cid, arm)])
        if ok:
            rows.append(r)
    assert rows, 'no case had scores from every function'
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with open(out/'rescore_results.tsv', 'w', newline='') as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)

    rng = np.random.default_rng(SEED)
    summary = {}
    for stratum in ['accept', 'backbone_change']:
        s = [r for r in rows if r['stratum'] == stratum]
        if not s:
            continue
        o = dict(cases=len(s), components=len({r['component'] for r in s}),
                 mean_poses_per_arm=round(float(np.mean([r['apo_single_poses'] for r in s])), 1),
                 oracle_success_any_apo_single=sum(int(r['apo_single_rmsd_oracle'] <= SUCCESS) for r in s),
                 oracle_success_any_apo_multi=sum(int(r['apo_multi_rmsd_oracle'] <= SUCCESS) for r in s))
        for name in funcs:
            e = {}
            for arm in ARMS:
                e[arm] = dict(success_top=sum(r[f'{name}_{arm}_success_top'] for r in s),
                              success_top5_raw=sum(r[f'{name}_{arm}_success_top5_raw'] for r in s),
                              mean_rmsd_top=round(float(np.mean([r[f'{name}_{arm}_rmsd_top'] for r in s])), 3),
                              median_rmsd_top=round(float(np.median([r[f'{name}_{arm}_rmsd_top'] for r in s])), 3))
            for side in ('holo', 'apo'):
                for lbl, key in [('top1', f'{name}_{side}_gain_top'), ('top5_raw', f'{name}_{side}_gain_top5_raw')]:
                    e[f'{side}_multi_minus_single_{lbl}'] = contrast([(r['component'], r[key]) for r in s], rng)
            o[name] = e
        summary[stratum] = o
    summary['note'] = (
        'EXPLORATORY. Only the selection step differs between scoring functions; the pose set and every RMSD '
        'are fixed and shared. Positive gain means the multiple-state receptor did better. Bootstrap unit is '
        'the homology component. Read each function\'s absolute success_top first: a function that picks poses '
        'worse than vina cannot settle the receptor question whatever its contrast says. gnina was trained on '
        'PDBbind-like complexes and our holo structures may be in that training set, so its absolute numbers '
        'are assumed optimistic; the paired contrast is unaffected because the leak enters both arms equally. '
        'top5_raw has no distinctness filter -- budget50 applies that.')
    (out/'rescore_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
