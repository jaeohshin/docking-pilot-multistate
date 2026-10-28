"""Dump (case, arm, call, pose, vina_score, rmsd) for every pose of a finished run, read-only.

The RMSD is derived exactly as the frozen analyses derive it: symmetry aware, in place (rdMolAlign.CalcRMS,
no alignment), against the reference ligand of that arm's own backbone. Writing it once here means every
scoring function in the rescoring experiment is judged against identical geometry, and that the frozen run
directory is never modified -- the table goes wherever --out says.
"""
import argparse, csv, gzip, json
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
ARMS = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--run', required=True)
    a.add_argument('--out', required=True)
    a = a.parse_args()
    run = Path(a.run)
    sel = {r['case_id']: r for r in csv.DictReader(open(run/'cohort_final.tsv'), delimiter='\t')}
    rows, n = [], 0
    for d in sorted((run/'docked').iterdir()):
        cid = d.name
        if cid not in sel or not all((d/x/'docking_log.json').exists() for x in ARMS):
            continue
        inp = run/'inputs'/cid
        refs = {s: Chem.RemoveHs(Chem.MolFromMolFile(str(inp/f'ligand_ref_{s}.sdf'), removeHs=False)) for s in ('holo', 'apo')}
        if any(v is None for v in refs.values()):
            continue
        for arm in ARMS:
            ref = refs[arm.split('_')[0]]
            for call in json.loads((d/arm/'docking_log.json').read_text())['calls']:
                if call['returncode']:
                    continue
                pm = PDBQTMolecule.from_file(str(d/arm/call['out']), skip_typing=True)
                rd = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0])
                sc = [float(p.score) for p in pm]
                assert rd.GetNumConformers() == len(sc), d/arm/call['out']
                rows += [(cid, arm, call['run'], i, round(s, 3),
                          round(rdMolAlign.CalcRMS(rd, ref, prbId=i), 3)) for i, s in enumerate(sc)]
        n += 1
        if n % 20 == 0:
            print(' cases', n, 'poses', len(rows), flush=True)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, 'wt', newline='') as h:
        w = csv.writer(h, delimiter='\t')
        w.writerow(['case_id', 'arm', 'call', 'pose', 'vina_score', 'rmsd'])
        w.writerows(rows)
    print('cases', n, 'poses', len(rows), '->', out)


if __name__ == '__main__':
    main()
