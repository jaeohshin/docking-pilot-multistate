"""Vina 1.2 docking for the 2x2 pilot: holo/apo backbone x single/multiple side-chain states.

Search budget is matched across arms (fixed before any result is seen):
  every arm gets --runs docking calls for a case. The single arm runs its one receptor with --runs different
  seeds; the multi arm spreads the same --runs calls over its receptors round-robin, so a multi arm with one
  receptor is exactly the single arm with the same seeds. Seeds are SEED0 + 0..runs-1, so a larger --runs is a
  strict superset of a smaller one: --runs 50 repeats expand_v1's 8 calls byte for byte and adds 42 more. Box, exhaustiveness and number of modes are
  identical everywhere; the box comes from the reference ligand of that backbone (inputs_manifest.json).

Each call writes its poses; scoring and RMSD are not computed here (see analyse_docking.py), so this script
cannot be tuned against the outcome. One case-arm per task: `--case <id> --arm holo_single|holo_multi|apo_single|apo_multi`.
"""
import argparse, json, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ENV = Path('/gpfs/jaeohshin/work/envs/dock-vina/bin')
EXHAUSTIVENESS, MODES, SEED0 = 32, 9, 20261028


def manifest():
    return {c['case_id']: c for c in json.loads((HERE/'inputs_manifest.json').read_text())['cases']}


def receptors(case, side, kind):
    """PDBQT paths for one arm, in a fixed order."""
    files = case[side]['receptors']
    if kind == 'single':
        return [files['single_0']]
    return [files[k] for k in sorted(files) if k.startswith('multi_')] or [files['single_0']]


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--case', required=True)
    a.add_argument('--arm', required=True, choices=['holo_single', 'holo_multi', 'apo_single', 'apo_multi'])
    a.add_argument('--cpu', type=int, default=8)
    a.add_argument('--runs', type=int, default=8)   # 8 keeps expand_v1 behaviour byte for byte
    a = a.parse_args()
    side, kind = a.arm.split('_')
    case = manifest()[a.case]
    inp, out = HERE/'inputs'/a.case, HERE/'docked'/a.case/a.arm
    out.mkdir(parents=True, exist_ok=True)
    recs = receptors(case, side, kind)
    box = case[side]
    log = dict(case=a.case, arm=a.arm, receptors=recs, box_center=box['box_center'], box_size=box['box_size'],
               exhaustiveness=EXHAUSTIVENESS, runs=a.runs, modes=MODES, calls=[])
    for k in range(a.runs):
        rec = recs[k % len(recs)]                     # round-robin: the budget is shared, not multiplied
        tag = f'run{k:02d}_{Path(rec).stem}'
        cmd = [ENV/'vina', '--receptor', inp/rec, '--ligand', inp/'ligand_dock.pdbqt',
               '--center_x', box['box_center'][0], '--center_y', box['box_center'][1], '--center_z', box['box_center'][2],
               '--size_x', box['box_size'][0], '--size_y', box['box_size'][1], '--size_z', box['box_size'][2],
               '--exhaustiveness', EXHAUSTIVENESS, '--num_modes', MODES, '--seed', SEED0 + k,
               '--cpu', a.cpu, '--out', out/f'{tag}.pdbqt']
        t = time.time()
        r = subprocess.run([str(x) for x in cmd], capture_output=True, text=True)
        log['calls'].append(dict(run=k, receptor=rec, seed=SEED0 + k, out=f'{tag}.pdbqt',
                                 returncode=r.returncode, seconds=round(time.time() - t, 1),
                                 stderr=r.stderr.strip()[-300:] if r.returncode else ''))
        print(a.case, a.arm, k, rec, f'{time.time()-t:.1f}s rc={r.returncode}', flush=True)
        if r.returncode:
            print(r.stderr[-800:], file=sys.stderr, flush=True)
    (out/'docking_log.json').write_text(json.dumps(log, indent=1) + '\n')
    bad = [c for c in log['calls'] if c['returncode']]
    print('done', a.case, a.arm, 'failed', len(bad))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
