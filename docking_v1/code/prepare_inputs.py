"""Vina inputs for the 2x2 pilot (design2x2_v1 receptors), standard Vina 1.2 / Meeko workflow. No docking here.

ligand    reference heavy atoms (design2x2_v1/receptors/<case>/ligand_ref_{holo,apo}.pdb) + CCD SMILES: bond orders
          from the template (AssignBondOrdersFromTemplate), explicit H (addCoords) -> reference SDF. Docking input:
          one RDKit ETKDG conformer (seed 20261028) + MMFF, prepared with Meeko -> PDBQT (starts away from the pose).
receptor  PDB2PQR (--noopt --nodebump) + PROPKA at pH 7.4 (AMBER names) on each backbone's single receptor; the titration states found
          there (HIS/ASP/GLU/LYS/CYS names) are written into every multi receptor of that backbone, which is then
          protonated with --titration-state-method none, so the arms of one backbone differ only in side-chain
          placement. Heavy atoms must come out unchanged (checked). Meeko mk_prepare_receptor -> rigid PDBQT.
box       centre = reference ligand centroid in the backbone's frame; edge per axis = max(extent + 10 A, 22.5 A);
          one box per backbone, shared by its arms.
Writes docking_v1/inputs/<case>/... and inputs_manifest.json.
"""
import json, hashlib, subprocess, sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
ENV = Path('/gpfs/jaeohshin/work/envs/dock-vina/bin')
SEED, PH = 20261028, 7.4
TITRATABLE = {'HIS', 'HID', 'HIE', 'HIP', 'ASP', 'ASH', 'GLU', 'GLH', 'LYS', 'LYN', 'CYS', 'CYX', 'CYM'}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode:
        raise RuntimeError(f'{Path(str(cmd[0])).name} failed: ' + ' | '.join(l.strip() for l in (r.stderr + r.stdout).splitlines()
                                                                         if l.strip() and 'WARNING' not in l)[-600:])
    return r


def smiles(ccd):
    pref = [d for d in ccd['descriptors'] if d[0] == 'SMILES_CANONICAL' and 'CACTVS' in d[1]] or \
           [d for d in ccd['descriptors'] if d[0].startswith('SMILES')]
    return pref[0][2].strip('"')


def reference_mol(pdb, smi):
    raw = Chem.MolFromPDBFile(str(pdb), removeHs=True, sanitize=False)
    tmpl = Chem.MolFromSmiles(smi)
    tmpl = Chem.RemoveHs(tmpl)
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    Chem.SanitizeMol(mol)
    return Chem.AddHs(mol, addCoords=True)


def docking_ligand(ref, out_pdbqt):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.Mol(ref); m.RemoveAllConformers()
    ps = AllChem.ETKDGv3(); ps.randomSeed = SEED
    assert AllChem.EmbedMolecule(m, ps) == 0, 'embedding failed'
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    setups = MoleculePreparation().prepare(m)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    assert ok, err
    Path(out_pdbqt).write_text(s)
    return m


def residue_names(pdb):
    st = gemmi.read_structure(str(pdb))
    return {(ch.name, r.seqid.num, r.seqid.icode): r.name for ch in st[0] for r in ch}


def heavy(pdb):
    st = gemmi.read_structure(str(pdb)); st.remove_hydrogens()
    # keyed without the residue name: PDB2PQR renames titratable residues (HIS -> HID ...)
    return {(ch.name, r.seqid.num, r.seqid.icode, a.name): np.array([a.pos.x, a.pos.y, a.pos.z]) for ch in st[0] for r in ch for a in r}


def all_atoms(pdb):
    """{(chain, num, icode, atom, resname-without-protonation-variant): xyz} including hydrogens."""
    canon = {'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'ASH': 'ASP', 'GLH': 'GLU', 'LYN': 'LYS', 'CYX': 'CYS', 'CYM': 'CYS'}
    st = gemmi.read_structure(str(pdb))
    return {(ch.name, r.seqid.num, r.seqid.icode, a.name, canon.get(r.name, r.name)): np.array([a.pos.x, a.pos.y, a.pos.z])
            for ch in st[0] for r in ch for a in r}


# Freely rotatable hydrogen groups: (anchor atom, pivot atom, hydrogens turned together about anchor->pivot).
ROTATABLE = {
    'SER': [('CB', 'OG', ['HG'])], 'CYS': [('CB', 'SG', ['HG'])], 'TYR': [('CZ', 'OH', ['HH'])],
    'THR': [('CB', 'OG1', ['HG1']), ('CB', 'CG2', ['HG21', 'HG22', 'HG23'])],
    'ASH': [('CG', 'OD2', ['HD2'])], 'GLH': [('CD', 'OE2', ['HE2'])],
    'LYS': [('CD', 'NZ', ['HZ1', 'HZ2', 'HZ3'])], 'LYN': [('CD', 'NZ', ['HZ2', 'HZ3'])],
    'ALA': [('CA', 'CB', ['HB1', 'HB2', 'HB3'])],
    'VAL': [('CB', 'CG1', ['HG11', 'HG12', 'HG13']), ('CB', 'CG2', ['HG21', 'HG22', 'HG23'])],
    'LEU': [('CG', 'CD1', ['HD11', 'HD12', 'HD13']), ('CG', 'CD2', ['HD21', 'HD22', 'HD23'])],
    'ILE': [('CB', 'CG2', ['HG21', 'HG22', 'HG23']), ('CG1', 'CD1', ['HD11', 'HD12', 'HD13'])],
    'MET': [('SD', 'CE', ['HE1', 'HE2', 'HE3'])],
}
PLANAR = {('TYR', 'HH'), ('ASH', 'HD2'), ('GLH', 'HE2')}   # conjugated: only the two in-plane positions


def rotate(p, a, b, ang):
    """Rotate point p about the axis a->b by ang radians (Rodrigues)."""
    k = (b - a) / np.linalg.norm(b - a); v = p - b
    return b + v * np.cos(ang) + np.cross(k, v) * np.sin(ang) + k * np.dot(k, v) * (1 - np.cos(ang))


def relax_rotatable_h(pdb):
    """With PDB2PQR optimisation off, rotatable hydrogens keep a default torsion and can overlap a neighbour,
    which also makes Meeko perceive a spurious bond (e.g. a valine methyl 1.36 A from a cysteine sulfur).
    Each rotatable group whose clearance is below 1.6 A is turned as a rigid unit about its own bond to the
    torsion with the largest clearance. Heavy atoms are never moved. Returns the groups changed."""
    st = gemmi.read_structure(str(pdb)); st.setup_entities()
    ns = gemmi.NeighborSearch(st[0], st.cell, 4).populate()

    def clearance(ch, r, points, own=()):
        """Smallest clearance of `points` to any atom that is not one of this group's own atoms.
        Atoms of the same residue count too, except the rotating hydrogens and the pivot they hang from:
        PDB2PQR can place an NH3 hydrogen 0.89 A from its own CE, which Meeko then reads as a bond."""
        best = 9.0
        for x in points:
            for m in ns.find_atoms(gemmi.Position(*x), '\0', radius=3.0):
                c = m.to_cra(st[0])
                same = (c.chain.name, c.residue.seqid.num) == (ch.name, r.seqid.num)
                if same and c.atom.name in own:
                    continue
                best = min(best, c.atom.pos.dist(gemmi.Position(*x)) - (0.0 if c.atom.element.name == 'H' else 0.35))
        return best

    changed = []
    for ch in st[0]:
        for r in ch:
            for anchor, pivot, names in ROTATABLE.get(r.name, []):
                A_, B_ = r.find_atom(anchor, '*'), r.find_atom(pivot, '*')
                hs = [r.find_atom(n, '*') for n in names]
                hs = [(n, h) for n, h in zip(names, hs) if h is not None]
                if not (A_ and B_ and hs):
                    continue
                a, b = np.array(A_.pos.tolist()), np.array(B_.pos.tolist())
                x0 = [np.array(h.pos.tolist()) for _, h in hs]
                own = set(names) | {pivot}
                if clearance(ch, r, x0, own) >= 1.6:
                    continue
                planar = (r.name, hs[0][0]) in PLANAR
                angles = [0.0, np.pi] if planar else np.radians(np.arange(0, 360, 10))
                best = max(angles, key=lambda t: clearance(ch, r, [rotate(x, a, b, t) for x in x0], own))
                if abs(best) < 1e-9:
                    continue
                for (_, h), x in zip(hs, x0):
                    h.pos = gemmi.Position(*rotate(x, a, b, best))
                changed.append(f'{ch.name}{r.seqid.num}:{r.name}:{pivot}')
    if changed:
        st.write_pdb(str(pdb))
    return changed


def pdb2pqr(src, dst, names=None):
    """Protonate src -> dst (PDB with H, AMBER names). With `names`, titration states are imposed, not predicted."""
    tmp = dst.with_suffix('.in.pdb')
    if names:
        st = gemmi.read_structure(str(src))
        for ch in st[0]:
            for r in ch:
                k = (ch.name, r.seqid.num, r.seqid.icode)
                if k in names and r.name in TITRATABLE | {'HIS'}:
                    r.name = names[k]
        st.write_pdb(str(tmp))
    else:
        tmp.write_text(Path(src).read_text())
    # --noopt: no hydrogen-bond optimisation, which flips ASN/GLN/HIS heavy atoms; heavy atoms must stay as built
    # --nodebump: debumping moved existing heavy atoms (e.g. a grafted ARG by 0.9 A); hydrogens are checked below
    args = [ENV/'pdb2pqr30', '--ff', 'AMBER', '--ffout', 'AMBER', '--keep-chain', '--noopt', '--nodebump', '--pdb-output', dst]
    args += ['--titration-state-method', 'propka', '--with-ph', PH] if not names else []
    run(args + [tmp, dst.with_suffix('.pqr')])
    before, after = heavy(src), heavy(dst)
    common = set(before) & set(after)
    shift = {k: float(np.linalg.norm(before[k] - after[k])) for k in common}
    moved = max(shift.values())
    assert moved < 1e-3, f'PDB2PQR moved heavy atoms by {moved:.3f} A in {src}'
    rotated = relax_rotatable_h(dst)
    st = gemmi.read_structure(str(dst)); st.setup_entities()
    ns = gemmi.NeighborSearch(st[0], st.cell, 3).populate()
    close = []
    for ch in st[0]:
        for r in ch:
            for a in r:
                if a.element.name != 'H':
                    continue
                for m in ns.find_atoms(a.pos, '\0', radius=1.2):
                    c = m.to_cra(st[0])
                    lim = 1.0 if c.atom.element.name == 'H' else 1.2
                    if (c.chain.name, c.residue.seqid.num) != (ch.name, r.seqid.num) and a.pos.dist(c.atom.pos) < lim:
                        close.append(f'{ch.name}{r.seqid.num}:{a.name}-{c.chain.name}{c.residue.seqid.num}:{c.atom.name}')
    # recorded, not asserted: what must hold is that heavy atoms are unchanged and that Meeko accepts the receptor.
    return dict(heavy_in=len(before), heavy_out=len(after), heavy_common=len(common), max_heavy_shift=moved,
                rotated_hydrogens=rotated, close_hydrogens=close[:10],
                moved_residues=sorted({f'{k[0]}:{k[1]}{k[2].strip()}' for k, v in shift.items() if v > 1e-3}))


def receptor_pdbqt(pdb_h, stem, center, size):
    # Meeko resolves HIS (incl. N/C-terminal variants) from the hydrogens present; HID/HIE/HIP names are not
    # ambiguous keys in its templates, so a terminal HID would fail. Other PDB2PQR names (ASH, GLH, LYN, CYX) are.
    meeko_in = Path(str(stem) + '.meeko.pdb')
    st = gemmi.read_structure(str(pdb_h))
    for ch in st[0]:
        for r in ch:
            if r.name in ('HID', 'HIE', 'HIP'):
                r.name = 'HIS'
    st.write_pdb(str(meeko_in))
    # unmatched residues are deleted only if their centre lies > 8 A (Vina cutoff) outside the box; otherwise an error
    r = run([ENV/'mk_prepare_receptor.py', '--read_pdb', meeko_in, '-o', stem, '-p', '-v',
             '--box_center', *center, '--box_size', *size, '--delete_bad_res_from_box_radius', 8])
    deleted = [l.strip() for l in (r.stdout + r.stderr).splitlines() if 'delet' in l.lower() or 'ignor' in l.lower()]
    return Path(f'{stem}.pdbqt'), deleted


def main():
    C = {c['case_id']: c for c in json.loads((PILOT/'final_review_v1/candidates79.json').read_text())}
    B = json.loads((PILOT/'design2x2_v1/build_manifest.json').read_text())
    only = set(sys.argv[1:])
    manifest, failures = [], []
    for case in B['cases']:
        cid = case['case_id']
        if only and cid not in only:
            continue
        try:
            manifest.append(prepare_case(case, C))
        except Exception as e:                       # collect every failure in one pass
            failures.append(dict(case_id=cid, error=repr(e))); print('FAILED', cid, repr(e)[:700], flush=True)
    (HERE/'inputs_manifest.json').write_text(json.dumps(dict(cases=manifest, failures=failures, seed=SEED, pH=PH, code_sha256=sha(__file__)), indent=1) + '\n')
    print('prepared', len(manifest), 'failed', len(failures))


def prepare_case(case, C):
    if True:
        cid = case['case_id']
        B_case = case
        c = C[cid]; src = PILOT/'design2x2_v1/receptors'/cid; out = HERE/'inputs'/cid; out.mkdir(parents=True, exist_ok=True)
        smi = smiles(c['ccd']); rec = dict(case_id=cid, smiles=smi)
        for side in ('holo', 'apo'):
            ref = reference_mol(src/f'ligand_ref_{side}.pdb', smi)
            Chem.MolToMolFile(ref, str(out/f'ligand_ref_{side}.sdf'))
            X = ref.GetConformer().GetPositions()[[a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum() > 1]]
            center = [round(float(v), 3) for v in X.mean(0)]
            size = [round(float(max(e + 10.0, 22.5)), 3) for e in (X.max(0) - X.min(0))]
            single = src/f'{side}_single_0.pdb'
            h_single = out/f'{side}_single_0.H.pdb'
            chk = {'single_0': pdb2pqr(single, h_single)}
            names = {k: v for k, v in residue_names(h_single).items() if v in TITRATABLE}
            f0, d0 = receptor_pdbqt(h_single, out/f'{side}_single_0', center, size)
            files, deleted = {'single_0': f0.name}, {'single_0': d0}
            for p in sorted(src.glob(f'{side}_multi_*.pdb')):
                tag = p.stem.split('_', 1)[1]
                h = out/f'{side}_{tag}.H.pdb'
                chk[tag] = pdb2pqr(p, h, names)
                f1, d1 = receptor_pdbqt(h, out/f'{side}_{tag}', center, size)
                files[tag], deleted[tag] = f1.name, d1
            rec[side] = dict(box_center=center, box_size=size, receptors=files, meeko_deleted=deleted, protonation_check=chk,
                             titration_states=sorted({v for v in names.values()}))
            # arms of one backbone must agree outside the pocket after protonation (heavy atoms and hydrogens)
            names_pdb = B_case[side]['chain_names_in_pdb']        # copy label (qFit chain + copy index) -> PDB chain
            pocket = {(names_pdb[l.split(':')[0]], l.split(':')[1]) for l in B_case[side]['pocket_residues']}
            ref_at = all_atoms(h_single)
            worst = 0.0
            for h in sorted(out.glob(f'{side}_multi_*.H.pdb')):
                at = all_atoms(h)
                assert set(at) == set(ref_at), (cid, h.name, 'atom set differs from single after protonation')
                for k, x in at.items():
                    rot = {n for _, _, ns_ in ROTATABLE.get(k[4], []) for n in ns_}
                    if (k[0], f'{k[4]}{k[1]}{k[2].strip()}') in pocket or k[3] in rot:
                        continue                         # pocket atoms, and hydroxyl H turned per receptor
                    worst = max(worst, float(np.linalg.norm(x - ref_at[k])))
            assert len({json.dumps(v) for v in deleted.values()}) == 1, (cid, side, 'Meeko deleted different residues across arms', deleted)
            rec[side]['max_nonpocket_difference_vs_single'] = worst
            assert worst < 1e-3, (cid, side, 'non-pocket atoms differ between arms', worst)
        docking_ligand(reference_mol(src/'ligand_ref_holo.pdb', smi), out/'ligand_dock.pdbqt')
        print(cid, {s: (len(rec[s]['receptors']), max(x['max_heavy_shift'] for x in rec[s]['protonation_check'].values())) for s in ('holo', 'apo')}, flush=True)
        return rec


if __name__ == '__main__':
    main()
