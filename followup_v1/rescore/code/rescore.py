"""Re-score poses that are already on disk with a different scoring function. No re-docking.

The pose set is fixed, identical to what Vina produced, and identical between arms in the sense that each arm
keeps its own poses. Only the *selection* step changes: for each arm we re-rank that arm's poses by the new
score and take the new rank-1 pose. RMSD is never recomputed here -- it is read from the pose table that
analyse_budget50.py (or analyse_rescore_from_expand.py) wrote, so the geometry side of the analysis cannot
drift between scoring functions.

Scoring functions:
  vinardo  smina --scoring vinardo --score_only            (CPU, empirical)
  gnina    gnina --score_only --cnn_scoring rescore        (CUDA, learned 3D CNN; use CNNaffinity)

One (case, arm) per task so this parallelises the same way the docking did.
Writes <out>/<case>__<arm>.tsv with one line per pose: call, pose, new score. Nothing else.
"""
import argparse, csv, json, re, subprocess, tempfile
from pathlib import Path

from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
SMINA = Path('/gpfs/jaeohshin/work/envs/rescore/bin/smina')
GNINA = Path('/gpfs/jaeohshin/work/envs/gnina-bin/gnina')


def to_sdf(pdbqt, dst):
    """Vina's multi-MODEL pdbqt is not readable as a ligand by smina/gnina, so hand them an SDF instead.

    The conversion goes through Meeko, the same path the RMSD is derived through, so conformer k of the SDF
    is pose k of the pdbqt and the new scores line up with the pose table without any reordering.
    """
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    rd = RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0]
    w = Chem.SDWriter(str(dst))
    for i in range(rd.GetNumConformers()):
        w.write(rd, confId=i)
    w.close()
    return rd.GetNumConformers()


def run_vinardo(receptor, ligand):
    r = subprocess.run([str(SMINA), '--receptor', str(receptor), '--ligand', str(ligand),
                        '--scoring', 'vinardo', '--score_only'], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-600:])
    return [float(x) for x in re.findall(r'^Affinity:\s+(-?\d+\.?\d*)', r.stdout, re.M)]


def run_gnina(receptor, ligand):
    """gnina's pose-ranking quantity is CNNscore; CNNaffinity is its affinity prediction.

    Both are recorded. CNNscore is the primary selector (gnina_dock/PLAN.md), CNNaffinity is kept because
    the first pass of this experiment selected on it and that result stands as reported. Both are
    higher-is-better, so both are negated to keep "lower score is better" everywhere downstream.
    """
    r = subprocess.run([str(GNINA), '--receptor', str(receptor), '--ligand', str(ligand),
                        '--score_only', '--cnn_scoring', 'rescore'], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-600:])
    sc = [float(x) for x in re.findall(r'^CNNscore:\s+(-?\d+\.?\d*)', r.stdout, re.M)]
    aff = [float(x) for x in re.findall(r'^CNNaffinity:\s+(-?\d+\.?\d*)', r.stdout, re.M)]
    if len(sc) != len(aff):
        raise RuntimeError(f'{len(sc)} CNNscore vs {len(aff)} CNNaffinity')
    return [-x for x in sc], [-x for x in aff]


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--run', required=True, help='directory holding inputs/ and docked/ (budget50 or expand_v1)')
    a.add_argument('--case', required=True)
    a.add_argument('--arm', required=True, choices=['holo_single', 'holo_multi', 'apo_single', 'apo_multi'])
    a.add_argument('--scoring', required=True, choices=['vinardo', 'gnina'])
    a.add_argument('--out', required=True)
    a = a.parse_args()
    run = Path(a.run)
    armdir = run/'docked'/a.case/a.arm
    log = json.loads((armdir/'docking_log.json').read_text())
    fn = run_vinardo if a.scoring == 'vinardo' else run_gnina
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        sdf = Path(tmp)/'poses.sdf'
        for call in log['calls']:
            if call['returncode']:
                continue
            n = to_sdf(armdir/call['out'], sdf)
            got = fn(run/'inputs'/a.case/call['receptor'], sdf)
            cols = got if isinstance(got, tuple) else (got,)
            for c in cols:
                if len(c) != n:
                    raise RuntimeError(f"{armdir/call['out']}: {n} poses, {len(c)} {a.scoring} scores")
            rows += [(call['run'], i) + tuple(round(c[i], 4) for c in cols) for i in range(n)]
    head = ['call', 'pose', 'score'] + (['score_alt'] if a.scoring == 'gnina' else [])
    with open(out/f'{a.case}__{a.arm}.tsv', 'w', newline='') as h:
        w = csv.writer(h, delimiter='\t'); w.writerow(head); w.writerows(rows)
    print('rescored', a.case, a.arm, a.scoring, len(rows), 'poses', flush=True)


if __name__ == '__main__':
    main()
