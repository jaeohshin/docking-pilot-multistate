"""Pocket side-chain state audit of the qFit holo and apo models for the 79 candidates (2x2 docking design:
holo/apo backbone x single/multiple side-chain states, no packer). Data only: no docking, inference or training.

Pocket = residues of the ligand's assembly copies with a heavy atom within 4.5 A of the ligand reference pose
(holo: deposited pose; apo: pose moved into the deposited apo frame), taken from final_selection_v1's
assembly_contact_audit.json and mapped to the qFit model by chain and author number (name checked).
Per pocket residue and qFit model: altloc states with occupancy >= 0.10 (105_joint_rotamers FLOOR), chi1 wells
(well3) and joint chi1/chi2 wells (chi2_well), occupancies renormalised, as the frozen analyses. A residue is
multi-state when >= 2 wells are occupied. Backbone spread = max distance between altloc copies of N/CA/C/O.
Checks: chi1 wells of the source receptor chain equal the replicate cache's; qFit cell/space group equal the
deposited ones; every pocket residue and receptor chain is present in the qFit model; ligand/HETs kept in holo.
"""
import csv, glob, hashlib, importlib.util, json, sys, zipfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
ROOT = PILOT.parents[1]
ZIP = Path('/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip')
spec = importlib.util.spec_from_file_location('j105', ROOT/'scripts/105_joint_rotamers.py')
j105 = importlib.util.module_from_spec(spec); sys.modules['j105'] = j105; spec.loader.exec_module(j105)
FLOOR = j105.FLOOR
CHI = {
    'ARG': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')], 'ASN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')],
    'ASP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'OD1')], 'CYS': [('N', 'CA', 'CB', 'SG')],
    'GLN': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')], 'GLU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'HIS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'ND1')], 'ILE': [('N', 'CA', 'CB', 'CG1'), ('CA', 'CB', 'CG1', 'CD1')],
    'LEU': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')], 'LYS': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')],
    'MET': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'SD')], 'PHE': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'PRO': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD')], 'SER': [('N', 'CA', 'CB', 'OG')], 'THR': [('N', 'CA', 'CB', 'OG1')],
    'TRP': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')], 'TYR': [('N', 'CA', 'CB', 'CG'), ('CA', 'CB', 'CG', 'CD1')],
    'VAL': [('N', 'CA', 'CB', 'CG1')]}
BB = ('N', 'CA', 'C', 'O')
Z = zipfile.ZipFile(ZIP)
_models = {}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def qfit(pid):
    if pid not in _models:
        st = gemmi.read_pdb_string(Z.read(f'qfit_PDBs/{pid}_qFit_010.pdb').decode('utf-8', 'replace'))
        st.setup_entities(); st.remove_hydrogens(); _models[pid] = st
    return _models[pid]


def dihedral(p):
    return float(np.degrees(gemmi.calculate_dihedral(*[gemmi.Position(*x) for x in p])))


def residue_states(res):
    """[(altloc, occupancy, chi list)] of the altloc states >= FLOOR; shared atoms count for every altloc."""
    alts = sorted({a.altloc for a in res if a.altloc not in ('\x00', ' ')}) or ['']
    out = []
    for alt in alts:
        atoms = {a.name: a for a in res if a.altloc in ('\x00', ' ')}
        atoms.update({a.name: a for a in res if a.altloc == alt})
        chis = []
        for q in CHI.get(res.name, []):
            if not all(n in atoms for n in q):
                break
            chis.append(dihedral([[atoms[n].pos.x, atoms[n].pos.y, atoms[n].pos.z] for n in q]))
        own = [a.occ for a in res if a.altloc == alt]
        occ = float(np.mean(own)) if own else 1.0
        if alt == '' or occ >= FLOOR - j105.TIE:
            out.append((alt, occ, chis))
    return out


def wells(resname, states):
    """chi1 and joint chi1/chi2 well occupancies, renormalised over the kept states with a chi1."""
    s = [(o, c) for _, o, c in states if c]
    if not s:
        return None, None
    tot = sum(o for o, _ in s) or 1.0
    w1, w12 = defaultdict(float), defaultdict(float)
    for o, c in s:
        w1[j105.well3(c[0])] += o / tot
        if len(c) > 1 and j105.chi2_class(resname):
            w12[j105.well3(c[0]) + j105.chi2_well(resname, c[1])] += o / tot
        else:
            w12[j105.well3(c[0])] += o / tot
    keep = lambda w: {k: round(v, 3) for k, v in w.items() if v >= FLOOR - j105.TIE}
    return keep(w1), keep(w12)


def backbone_spread(res):
    d = 0.0
    for n in BB:
        P = [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.name == n]
        for a, b in combinations(P, 2):
            d = max(d, float(np.linalg.norm(np.subtract(a, b))))
    return d


def find(st, chain, num, icode, name):
    hits = [r for ch in st[0] if ch.name == chain for r in ch if r.seqid.num == num and r.seqid.icode.strip() == icode]
    hits = [r for r in hits if r.name == name]
    return hits[0] if len(hits) == 1 else None


def parse_label(lbl):
    """'A:HIS4' / 'A:SER127B' -> (chain, name, num, icode)."""
    ch, rest = lbl.split(':')
    name, tail = rest[:3], rest[3:]
    icode = tail[-1] if tail and tail[-1].isalpha() else ''
    return ch, name, int(tail[:-1] if icode else tail), icode


def pocket(records, assembly, side_chain_map):
    """Unique pocket residues (qFit chain, name, num, icode) from in-assembly polymer contacts; symmetry copies of a
    chain map to the same qFit residues."""
    out = {}
    for r in records:
        if r['kind'] != 'Polymer' or not r['contact_atoms'] or not r['in_assembly'].get(assembly):
            continue
        for lbl in r['contact_residues']:
            ch, name, num, ic = parse_label(lbl)
            out[(side_chain_map.get(ch, ch), name, num, ic)] = r['copy']
    return out


def audit_side(pid, records, assembly, chain_map):
    st = qfit(pid)
    res_rows, missing = [], []
    for (ch, name, num, ic), copy in sorted(pocket(records, assembly, chain_map).items()):
        r = find(st, ch, num, ic, name)
        if r is None:
            missing.append(f'{ch}:{name}{num}{ic}'); continue
        states = residue_states(r)
        w1, w12 = wells(name, states)
        res_rows.append(dict(residue=f'{ch}:{name}{num}{ic}', copy=copy, altlocs=len(states),
                             chi1_wells=w1, chi12_wells=w12, chi1_states=len(w1) if w1 else 0,
                             chi12_states=len(w12) if w12 else 0, backbone_spread=round(backbone_spread(r), 3)))
    return res_rows, missing


def side_summary(rows):
    sc = [r for r in rows if r['chi1_wells']]
    return dict(pocket_residues=len(rows), with_chi=len(sc), multi_chi1=sum(r['chi1_states'] >= 2 for r in sc),
                multi_chi12=sum(r['chi12_states'] >= 2 for r in sc),
                multi_chi1_residues=[r['residue'] for r in sc if r['chi1_states'] >= 2],
                multi_chi12_residues=[r['residue'] for r in sc if r['chi12_states'] >= 2],
                max_backbone_spread=max([r['backbone_spread'] for r in rows] + [0.0]),
                residues_backbone_spread_gt_0p5=sum(r['backbone_spread'] > 0.5 for r in rows),
                chi12_combinations=int(np.prod([max(1, r['chi12_states']) for r in sc])) if sc else 1)


def cell_match(qst, dep_path):
    d = gemmi.read_structure(str(dep_path))
    a = np.array(qst.cell.parameters); b = np.array(d.cell.parameters)
    return bool(qst.spacegroup_hm.replace(' ', '') == d.spacegroup_hm.replace(' ', '') and np.all(np.abs(a - b) <= np.maximum(0.005 * b, 0.05)))


def main():
    C = {c['case_id']: c for c in json.loads((PILOT/'final_review_v1/candidates79.json').read_text())}
    AUD = {o['case_id']: o for o in json.loads((PILOT/'final_selection_v1/assembly_contact_audit.json').read_text())['cases']}
    REV = {r['case_id']: r for r in csv.DictReader(open(PILOT/'final_selection_v1/case_review_all79.tsv'), delimiter='\t')}
    cache = {}
    need = {c['holo'] for c in C.values()} | {c['apo'] for c in C.values()}
    for f in sorted(glob.glob(str(ROOT/'reference/frozen/data/qfit_at_scale/replicate_cache/*.jsonl'))):
        for line in open(f):
            if line[12:16] in need:
                e = json.loads(line); cache[e['pdb_id']] = e
    out, cache_checks = [], []
    for cid, c in C.items():
        o = AUD[cid]
        # deposited -> qFit chain names: the receptor chains have recorded aliases; other chains keep their name
        hmap = {c['holo_deposited_chain']: c['holo_chain']}
        amap = {c['apo_deposited_chain']: c['apo_chain']}
        hrows, hmiss = audit_side(c['holo'], o['holo_records'], o['holo_assembly_used'], hmap)
        arows, amiss = audit_side(c['apo'], o['apo_records'], o['apo_assembly_used'], amap)
        # chi1 wells of the receptor chain against the replicate cache (independent implementation check)
        for side, pid, ch, rows in (('holo', c['holo'], c['holo_chain'], hrows), ('apo', c['apo'], c['apo_chain'], arows)):
            st = qfit(pid)
            poly = [p for p in st[0] if p.name == ch and len(p.get_polymer())][0].get_polymer()
            idx = {(r.seqid.num, r.seqid.icode.strip()): i for i, r in enumerate(poly)}
            states = {s['i']: s for s in cache[pid]['chains'][ch]['states']}
            for r in rows:
                rc, name, num, ic = parse_label(r['residue'])
                if rc != ch or not r['chi1_wells']:
                    continue
                i = idx[(num, ic)]
                if i in states:
                    w = states[i]['wells']; tot = sum(w.values()) or 1
                    ref = {k for k, v in w.items() if v / tot >= FLOOR - 1e-9}
                    cache_checks.append(dict(case_id=cid, side=side, residue=r['residue'], match=ref == set(r['chi1_wells'])))
        hq, aq = qfit(c['holo']), qfit(c['apo'])
        comp, lch, lnum = c['ligand_id'].split(':')
        lig_in_qfit = any(r.name == comp and r.seqid.num == int(lnum) for ch in hq[0] if ch.name == c['ligand_id'].split(':')[1] for r in ch)
        hets = sorted({x for n in o['holo_nonpolymer_neighbours'] for x in n['components']})
        hq_hets = {r.name for ch in hq[0] for r in ch if r.het_flag == 'H'}
        h, a = side_summary(hrows), side_summary(arows)
        out.append(dict(case_id=cid, homology_component=o['homology_component'], apo=c['apo'], holo=c['holo'], ligand=c['ligand_id'],
                        selection_status=REV[cid]['status'], selection_reason=REV[cid]['reason'],
                        pilot_representative=int(REV[cid]['pilot_representative']),
                        holo_states=h, apo_states=a, holo_missing_in_qfit=hmiss, apo_missing_in_qfit=amiss,
                        holo_cell_match=cell_match(hq, PILOT/'final_review_v1/cifs'/f"{c['holo']}.cif.gz"),
                        apo_cell_match=cell_match(aq, PILOT/'final_review_v1/cifs'/f"{c['apo']}.cif.gz"),
                        ligand_in_holo_qfit=lig_in_qfit, holo_het_neighbours=hets,
                        holo_het_neighbours_missing_in_qfit=sorted(set(hets) - hq_hets),
                        holo_residues=hrows, apo_residues=arows))
        print(cid, 'holo multi chi1/chi12', h['multi_chi1'], h['multi_chi12'], 'apo', a['multi_chi1'], a['multi_chi12'], flush=True)
    res = dict(cases=out, cache_check=dict(n=len(cache_checks), mismatches=[x for x in cache_checks if not x['match']]),
               code_sha256=sha(__file__), inputs_sha256={p: sha(PILOT/p) for p in ['final_review_v1/candidates79.json',
               'final_selection_v1/assembly_contact_audit.json', 'final_selection_v1/case_review_all79.tsv']},
               qfit_zip=str(ZIP), gemmi=gemmi.__version__)
    (HERE/'qfit_state_audit.json').write_text(json.dumps(res, indent=1) + '\n')
    print('cache check', len(cache_checks), 'mismatches', len(res['cache_check']['mismatches']))


if __name__ == '__main__':
    main()
