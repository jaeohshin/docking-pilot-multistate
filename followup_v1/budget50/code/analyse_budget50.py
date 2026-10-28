"""50-call analysis, fixed by followup_v1/PLAN.md before any 50-call pose existed.

Everything that decides the confirmatory test is inherited unchanged from expand_v1/code/analyse_expand.py:
the primary metric (symmetry-aware in-place heavy-atom RMSD of the best-Vina-scoring pose), the primary
contrast (apo multiple minus single as a paired per-case improvement), the cluster bootstrap over homology
components, and the rule that the CI lower bound must exceed zero.

There are two co-primary tests (PLAN.md, declared before any 50-call pose existed): P1 on the rank-1 pose,
which is the only condition comparable with the closed 8-call result, and P2 on the best of the top 5
*distinct* poses, which is far less sensitive to a single scoring mistake. Distinct means each selected pose
is more than DISTINCT A from every pose already selected, walking the score order; without that, the top 5 of
450 poses can be five copies of one pose and P2 collapses onto P1. Two co-primaries are judged at a
Bonferroni-corrected 97.5% two-sided interval, and the uncorrected 95% one is reported next to it.

What this file adds are the sampling metrics PLAN.md declared, which separate generating the right pose from
selecting it:
  frac_good   fraction of all poses within 2 A of the reference  -- generation density
  calls_good  calls (out of --runs) that produced at least one such pose -- reproducibility of generation
  rank_best   1-based rank of the closest pose when all poses are sorted by Vina score -- the size of the
              selection failure. rank_best == 1 means scoring picked the best pose it had; a large rank means
              the right answer was generated and buried.
It also writes pose_table.tsv.gz (case, arm, call, pose, score, rmsd) so the rescoring experiment can swap in
another scoring function without re-deriving any RMSD.
"""
import csv, gzip, json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
HERE = Path(__file__).resolve().parents[1]
ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
SUCCESS, BOOT, SEED = 2.0, 10000, 20261028
TOPN, DISTINCT, CURVE = 5, 1.0, 10        # top-5 distinct poses, 1 A apart; success curve to N=10


def scored(pdbqt, ref):
    """[(vina score, rmsd to ref, conformer id, mol)] for one Vina output file."""
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    rd = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0])
    s = [float(p.score) for p in pm]
    assert rd.GetNumConformers() == len(s), pdbqt
    return [(x, rdMolAlign.CalcRMS(rd, ref, prbId=i), i, rd) for i, x in enumerate(s)]


def pick_distinct(order, n, thresh):
    """Walk the score order and keep a pose only if it is > thresh from every pose already kept.

    `order` is [(score, rmsd, conf_id, mol)] sorted by score. Distinctness is the same symmetry-aware in-place
    RMSD used everywhere else, so "distinct" means the same thing as the metric does.
    """
    kept = []
    for cand in order:
        if all(rdMolAlign.CalcRMS(cand[3], k[3], prbId=cand[2], refId=k[2]) > thresh for k in kept):
            kept.append(cand)
            if len(kept) == n:
                break
    return kept


def arm_result(d, arm, ref, sink, cid):
    """Per-arm summary, and every pose appended to `sink` for the rescoring experiment."""
    log = json.loads((d/arm/'docking_log.json').read_text())
    out, per_call = [], []
    for call in log['calls']:
        if call['returncode']:
            continue
        ps = scored(d/arm/call['out'], ref)
        per_call.append(min(p[1] for p in ps))
        for sc, rm, i, _ in ps:
            sink.append((cid, arm, call['run'], i, round(sc, 3), round(rm, 3)))
        out += ps
    if not out:
        return None
    order = sorted(out, key=lambda x: x[0])                    # ascending Vina score: best first
    r = np.array([x[1] for x in out])
    rank_best = int(np.argmin([x[1] for x in order])) + 1
    raw5 = min(x[1] for x in order[:TOPN])                     # no dedupe: the secondary version
    dist = pick_distinct(order, TOPN, DISTINCT)
    top5 = min(x[1] for x in dist)
    curve = {f'top{n}': int(min(x[1] for x in order[:n]) <= SUCCESS) for n in range(1, CURVE + 1)}
    return dict(poses=len(out), calls=len(per_call), receptors=len(log['receptors']),
                best_score=round(order[0][0], 3), rmsd_top=round(order[0][1], 3),
                rmsd_top5=round(float(top5), 3), rmsd_top5_raw=round(float(raw5), 3), n_distinct=len(dist),
                rmsd_oracle=round(float(r.min()), 3),
                success_top=int(order[0][1] <= SUCCESS), success_top5=int(top5 <= SUCCESS),
                success_any=int(r.min() <= SUCCESS),
                frac_good=round(float((r <= SUCCESS).mean()), 4),
                calls_good=sum(1 for x in per_call if x <= SUCCESS), rank_best=rank_best, **curve)


def cluster_boot(pairs, rng, lo=2.5, hi=97.5):
    """pairs: [(component, value)]. Resample components, average within one first, then across."""
    comps = sorted({c for c, _ in pairs})
    by = {c: np.array([v for cc, v in pairs if cc == c], float) for c in comps}
    means = np.array([by[c].mean() for c in comps])
    if len(comps) < 2:
        return [float('nan')] * 2
    draws = means[rng.integers(0, len(comps), size=(BOOT, len(comps)))].mean(1)
    return [round(float(np.percentile(draws, lo)), 4), round(float(np.percentile(draws, hi)), 4)]


def contrast(sel, key, rng, coprimary=False):
    """One paired contrast. A co-primary also carries the Bonferroni-corrected interval it is judged on."""
    pairs = [(r['component'], r[key]) for r in sel]
    cm = [np.mean([v for c, v in pairs if c == q]) for q in sorted({c for c, _ in pairs})]
    o = dict(mean_of_component_means=round(float(np.mean(cm)), 4), ci95=cluster_boot(pairs, rng),
             components=len(cm), cases=len(pairs))
    if coprimary:
        # two co-primary tests, family-wise 5%: judge on the 97.5% interval, report the 95% one beside it
        o['ci975_bonferroni'] = cluster_boot(pairs, rng, 1.25, 98.75)
        o['decided_on'] = 'ci975_bonferroni'
    return o


def main():
    comp = {c['case_id']: c['homology_component'] for c in json.loads((HERE/'candidates.json').read_text())}
    sel = {r['case_id']: r for r in csv.DictReader(open(HERE/'cohort_final.tsv'), delimiter='\t')}
    rows, sink = [], []
    for d in sorted((HERE/'docked').iterdir()):
        cid = d.name
        if cid not in sel or not all((d/a/'docking_log.json').exists() for a in ARMS):
            continue
        inp = HERE/'inputs'/cid
        refs = {s: Chem.RemoveHs(Chem.MolFromMolFile(str(inp/f'ligand_ref_{s}.sdf'), removeHs=False)) for s in ('holo', 'apo')}
        if any(v is None for v in refs.values()):
            continue
        mine = []
        res = {a: arm_result(d, a, refs[a.split('_')[0]], mine, cid) for a in ARMS}
        if any(v is None for v in res.values()):
            continue
        sink += mine
        r = dict(case_id=cid, component=comp[cid], stratum=sel[cid]['status'], apo=sel[cid]['apo'],
                 holo=sel[cid]['holo'], ligand=sel[cid]['ligand'])
        for a in ARMS:
            for k, v in res[a].items():
                r[f'{a}_{k}'] = v
        for side in ('holo', 'apo'):
            r[f'{side}_gain_top'] = round(r[f'{side}_single_rmsd_top'] - r[f'{side}_multi_rmsd_top'], 3)
            r[f'{side}_gain_top5'] = round(r[f'{side}_single_rmsd_top5'] - r[f'{side}_multi_rmsd_top5'], 3)
            r[f'{side}_gain_top5_raw'] = round(r[f'{side}_single_rmsd_top5_raw'] - r[f'{side}_multi_rmsd_top5_raw'], 3)
            r[f'{side}_gain_oracle'] = round(r[f'{side}_single_rmsd_oracle'] - r[f'{side}_multi_rmsd_oracle'], 3)
            r[f'{side}_gain_frac_good'] = round(r[f'{side}_multi_frac_good'] - r[f'{side}_single_frac_good'], 4)
            r[f'{side}_gain_calls_good'] = r[f'{side}_multi_calls_good'] - r[f'{side}_single_calls_good']
        r['has_multiple_states'] = int(r['apo_multi_receptors'] > 1 or r['holo_multi_receptors'] > 1)
        rows.append(r)
        if len(rows) % 20 == 0:
            print(' scored', len(rows), flush=True)
    assert rows, 'no complete cases'
    with open(HERE/'docking_results.tsv', 'w', newline='') as h:
        w = csv.DictWriter(h, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
    with gzip.open(HERE/'pose_table.tsv.gz', 'wt', newline='') as h:
        w = csv.writer(h, delimiter='\t'); w.writerow(['case_id', 'arm', 'call', 'pose', 'vina_score', 'rmsd'])
        w.writerows(sink)
    print('poses written', len(sink))

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
                        mean_rmsd_top5=round(float(np.mean([r[f'{a}_rmsd_top5'] for r in s])), 3),
                        median_rmsd_top5=round(float(np.median([r[f'{a}_rmsd_top5'] for r in s])), 3),
                        success_top=sum(r[f'{a}_success_top'] for r in s),
                        success_top5=sum(r[f'{a}_success_top5'] for r in s),
                        success_any=sum(r[f'{a}_success_any'] for r in s),
                        mean_n_distinct=round(float(np.mean([r[f'{a}_n_distinct'] for r in s])), 2),
                        success_curve={f'top{n}': sum(r[f'{a}_top{n}'] for r in s) for n in range(1, CURVE + 1)},
                        mean_rmsd_oracle=round(float(np.mean([r[f'{a}_rmsd_oracle'] for r in s])), 3),
                        mean_frac_good=round(float(np.mean([r[f'{a}_frac_good'] for r in s])), 4),
                        mean_calls_good=round(float(np.mean([r[f'{a}_calls_good'] for r in s])), 2),
                        median_rank_best=float(np.median([r[f'{a}_rank_best'] for r in s])),
                        mean_poses=round(float(np.mean([r[f'{a}_poses'] for r in s])), 1))
        # buried: the right pose was generated but Vina did not put it first
        for a in ARMS:
            gen = [r for r in s if r[f'{a}_success_any']]
            o[a]['generated_good'] = len(gen)
            o[a]['buried'] = sum(1 for r in gen if not r[f'{a}_success_top'])
        # the two co-primary tests; everything after them is secondary or exploratory
        o['P1_apo_multi_minus_single_top1'] = contrast(s, 'apo_gain_top', rng, coprimary=True)
        o['P2_apo_multi_minus_single_top5'] = contrast(s, 'apo_gain_top5', rng, coprimary=True)
        for name, key in [('apo_top5_raw', 'apo_gain_top5_raw'),
                          ('holo_multi_minus_single_top1', 'holo_gain_top'),
                          ('holo_multi_minus_single_top5', 'holo_gain_top5'),
                          ('apo_oracle', 'apo_gain_oracle'), ('holo_oracle', 'holo_gain_oracle'),
                          ('apo_frac_good', 'apo_gain_frac_good'), ('holo_frac_good', 'holo_gain_frac_good'),
                          ('apo_calls_good', 'apo_gain_calls_good'), ('holo_calls_good', 'holo_gain_calls_good')]:
            o[name] = contrast(s, key, rng)
        sub = [r for r in s if r['has_multiple_states']]
        if sub:
            o['apo_top1_states_only'] = contrast(sub, 'apo_gain_top', rng)
            o['apo_top5_states_only'] = contrast(sub, 'apo_gain_top5', rng)
        summary[stratum] = o
    summary['note'] = (
        'Positive gain means the multiple-state receptor did better. Bootstrap unit is the homology component. '
        'Two CO-PRIMARY tests in stratum accept (followup_v1/PLAN.md, fixed before any 50-call pose existed): '
        'P1 on the rank-1 pose and P2 on the best of the top-5 distinct poses (>1 A apart). Each is decided on '
        'ci975_bonferroni, whose lower bound must exceed zero, keeping the family-wise error at 5%; ci95 is '
        'printed beside it for transparency and is not the decision rule. Everything else is secondary or '
        'exploratory and is reported whatever it shows. Oracle, top5 and success counts depend on the call '
        'budget, so they are never compared across the 8-call and 50-call runs.')
    (HERE/'docking_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
