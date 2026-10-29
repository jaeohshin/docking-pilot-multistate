"""gnina 1.1 docking for the 2x2 pilot: holo/apo backbone x single/multiple side-chain states.

This is experiment C (gnina_dock/PLAN.md): gnina both generates and selects, so the sampling-scoring
mismatch of the rescoring experiment is gone. Receptors, ligands and boxes are expand_v1's, reused
unchanged, and the per-call settings match it exactly (exhaustiveness 32, num_modes 9, 8 calls, seeds
SEED0 + 0..7), so this run, expand_v1 and the rescoring sit side by side as a three-way comparison in
which only the search and selection function differ.

Search budget is matched across arms, as in every arm of this pilot: the single arm runs its one receptor
with RUNS seeds, the multi arm spreads the same RUNS calls over its receptors round-robin, so a multi arm
with one receptor is exactly the single arm with the same seeds. A multi arm therefore cannot win by
searching more.

Each call writes its poses and gnina's own scores; no RMSD is computed here (see analyse_gnina.py), so
this script cannot be tuned against the outcome. One case-arm per task.
"""
import argparse, json, subprocess, sys, tempfile, time
from pathlib import Path

from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog('rdApp.*')
HERE = Path(__file__).resolve().parents[1]
GNINA = Path('/gpfs/jaeohshin/work/envs/gnina-bin/gnina')
EXHAUSTIVENESS, RUNS, MODES, SEED0 = 32, 8, 9, 20261028


def manifest():
    return {c['case_id']: c for c in json.loads((HERE/'inputs_manifest.json').read_text())['cases']}


def ligand_sdf(pdbqt, dst):
    """Hand gnina the ligand as SDF, carrying over Vina's exact starting conformer.

    gnina 1.1 descends from Vina 1.1.2 and rejects Meeko's macrocycle pseudo-atom types (CG0, G0 ...),
    which 9 of our 154 ligands carry. Converting through Meeko -- the same path the RMSD is derived
    through -- keeps the starting geometry identical to the Vina run for every case and keeps the input
    format uniform across cases rather than mixing pdbqt and sdf. gnina then does its own rotatable-bond
    detection from the SDF, which is part of letting it run the search end to end.
    """
    pm = PDBQTMolecule.from_file(str(pdbqt), skip_typing=True)
    rd = RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0]
    assert rd.GetNumConformers() == 1, f'{pdbqt}: {rd.GetNumConformers()} starting conformers'
    w = Chem.SDWriter(str(dst)); w.write(rd, confId=0); w.close()
    return dst


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
    a = a.parse_args()
    side, kind = a.arm.split('_')
    case = manifest()[a.case]
    inp, out = HERE/'inputs'/a.case, HERE/'docked'/a.case/a.arm
    out.mkdir(parents=True, exist_ok=True)
    recs = receptors(case, side, kind)
    box = case[side]
    log = dict(case=a.case, arm=a.arm, engine='gnina-1.1', receptors=recs, box_center=box['box_center'],
               box_size=box['box_size'], exhaustiveness=EXHAUSTIVENESS, runs=RUNS, modes=MODES,
               ligand_format='sdf from ligand_dock.pdbqt via meeko', calls=[])
    with tempfile.TemporaryDirectory() as tmp:
        lig = ligand_sdf(inp/'ligand_dock.pdbqt', Path(tmp)/'ligand.sdf')
        for k in range(RUNS):
            rec = recs[k % len(recs)]                 # round-robin: the budget is shared, not multiplied
            tag = f'run{k:02d}_{Path(rec).stem}'
            cmd = [GNINA, '--receptor', inp/rec, '--ligand', lig,
                   '--center_x', box['box_center'][0], '--center_y', box['box_center'][1], '--center_z', box['box_center'][2],
                   '--size_x', box['box_size'][0], '--size_y', box['box_size'][1], '--size_z', box['box_size'][2],
                   '--exhaustiveness', EXHAUSTIVENESS, '--num_modes', MODES, '--seed', SEED0 + k,
                   '--cpu', a.cpu, '--out', out/f'{tag}.sdf']
            t = time.time()
            r = subprocess.run([str(x) for x in cmd], capture_output=True, text=True)
            log['calls'].append(dict(run=k, receptor=rec, seed=SEED0 + k, out=f'{tag}.sdf',
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
