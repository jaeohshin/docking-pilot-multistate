"""Checks of the built receptors: (1) single = qFit major altloc coordinates for ASU copies; (2) every multi receptor has
the single receptor's backbone and non-pocket atoms exactly, differing only in pocket side chains; (3) reference ligand
in the holo frame vs the ligand in the holo qFit model; (4) apo qFit vs deposited apo receptor-chain CA (no fit)."""
import json, sys
from pathlib import Path
import gemmi, numpy as np
HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE/'code'))
import build_receptors as B
M = json.loads((HERE/'build_manifest.json').read_text())
C = {c['case_id']: c for c in json.loads((B.PILOT/'final_review_v1/candidates79.json').read_text())}
def atoms(p):
    st = gemmi.read_structure(str(p)); out = {}
    for ch in st[0]:
        for r in ch:
            for a in r:
                out[(ch.name, r.seqid.num, r.seqid.icode, r.name, a.name)] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out
rep = []
for case in M['cases']:
    cid = case['case_id']; c = C[cid]; d = HERE/'receptors'/cid; row = dict(case_id=cid)
    for side in ('holo', 'apo'):
        single = atoms(d/f'{side}_single_0.pdb')
        diffs, bbmax = [], 0.0
        for p in sorted(d.glob(f'{side}_multi_*.pdb')):
            m = atoms(p)
            assert set(m) == set(single), (cid, p.name, 'atom set differs')
            moved = {k[:4] for k in m if np.linalg.norm(m[k] - single[k]) > 1e-3}
            bb = max([np.linalg.norm(m[k] - single[k]) for k in m if k[4] in ('N', 'CA', 'C', 'O')] + [0])
            bbmax = max(bbmax, bb); diffs.append(len(moved))
        row[f'{side}_multi_files'] = len(diffs); row[f'{side}_moved_residues'] = diffs; row[f'{side}_backbone_max_change'] = bbmax
        # single vs qFit major altloc (first chain copy = ASU identity in all cases here)
        st = B.S.qfit(c[side]); maj = {}
        for ch in st[0]:
            for r in ch.get_polymer():
                for a in B.A.select_atoms(r)[0]:
                    maj[(r.seqid.num, r.seqid.icode.strip(), r.name, a.name, ch.name)] = np.array([a.pos.x, a.pos.y, a.pos.z])
        rch = c[f'{side}_chain']
        pdb_chain = case[side]['chain_names_in_pdb'].get(f'{rch}0')
        dev = [np.linalg.norm(x - maj[(k[1], k[2].strip(), k[3], k[4], rch)]) for k, x in single.items()
               if k[0] == pdb_chain and (k[1], k[2].strip(), k[3], k[4], rch) in maj]
        row[f'{side}_single_vs_qfit_major_max'] = max(dev) if dev else None
        row[f'{side}_single_clashes'] = case[side]['single_clashes']
    hq = B.S.qfit(c['holo']); comp, lch, lnum = c['ligand_id'].split(':')
    lq = [r for ch in hq[0] if ch.name == lch for r in ch if r.name == comp and r.seqid.num == int(lnum)][0]
    la = {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) for a in B.A.select_atoms(lq)[0]}
    ref = atoms(d/'ligand_ref_holo.pdb')
    row['ligand_ref_vs_qfit_holo_rmsd'] = float(np.sqrt(np.mean([np.sum((x - la[k[4]]) ** 2) for k, x in ref.items() if k[4] in la])))
    aq = B.S.qfit(c['apo']); ad = B.A.load(c['apo'])
    ca = lambda st, ch: {(r.seqid.num, r.seqid.icode.strip()): np.array([a.pos.x, a.pos.y, a.pos.z]) for x in st[0] if x.name == ch for r in x for a in r if a.name == 'CA'}
    q, dpos = ca(aq, c['apo_chain']), ca(ad, c['apo_deposited_chain'])
    row['apo_qfit_vs_deposited_ca_rmsd_nofit'] = float(np.sqrt(np.mean([np.sum((q[k] - dpos[k]) ** 2) for k in q if k in dpos])))
    rep.append(row)
    print(row)
(HERE/'verify_receptors.json').write_text(json.dumps(rep, indent=1, default=float) + '\n')
