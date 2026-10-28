"""Biological-assembly and crystal-contact audit of the 79 fixed docking-pilot candidates (deposited mmCIF).

For each case (apo, holo, receptor chain, ligand instance):
  holo  every polymer/non-polymer subchain copy under every space-group operator and nearby lattice translation is
        placed explicitly; copies with a heavy atom within 6 A of the ligand reference pose are kept. A copy is
        'asu' (identity operator, zero translation) or 'symmetry'. Each copy is tested for membership of every
        deposited assembly, relative to the ligand's own copy in that assembly (g_lig^-1 g_copy == crystal transform,
        checked on the copy's atoms, tolerance 1 A).
  apo   the ligand pose is moved into the deposited apo frame by a pocket backbone fit through `residue_map`
        (auth keys of the receptor chain), and the same decomposition is made relative to the apo receptor chain.
Also recorded: ligand struct_conn records, instances of the same component, validation-report RSCC/RSR, non-water
HET neighbours (cofactors, buffers) and whether the same component sits at the apo site.

Contact = heavy atom pair <= 4.5 A (the source contact definition). Fixed rules, set before looking at results:
  crystal_dependence  none: no non-assembly polymer atom within 4.5 A; minor: non-assembly share of the ligand's
                      protein contact atoms < 0.15 and of ligand atoms < 0.25; major otherwise.
  cofactor_dependence a non-polymer, non-water neighbour of >= 6 heavy atoms (or a listed cofactor) contacting
                      >= 10% of ligand atoms at 4.0 A that has no same-component copy within 3 A at the apo site.
No docking, inference or training. Input files are read only; output goes to this folder.
"""
import gzip, hashlib, json, sys, xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
ROOT = PILOT.parents[1]
CIFS = PILOT/'final_review_v1/cifs'
CANDS = PILOT/'final_review_v1/candidates79.json'
RECORDS = PILOT/'filter_review_v1/review_records.json'
VAL = HERE/'validation_xml'
CONTACT, SHELL, COF_CONTACT = 4.5, 6.0, 4.0
SOLVENT = {'HOH', 'DOD'}
sys.path.insert(0, str(ROOT/'scripts'))
COFACTORS = set("""HEM HEC HEA HEB HDD FAD FMN FDA NAD NAP NDP NAI NAJ NMN ADP ATP ANP ACP AMP GDP GTP GNP GSP GCP
SAM SAH SFG COA ACO CAA 01A PLP PMP TPP TDP BTN CLA BCL CHL SF4 FES F3S MGD PQQ TOP UMP UDP UTP CTP CDP C5P B12
COB H4B THF MTX FOL RET GSH""".split())  # 178_pocket_variation.COFACTORS + GSH


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def xyz(atoms): return np.array([[a.pos.x, a.pos.y, a.pos.z] for a in atoms], float).reshape(-1, 3)


def select_atoms(res):
    """Heavy atoms of the highest-mean-occupancy altloc (same rule as the geometry audit)."""
    ats = [a for a in res if not a.is_hydrogen()]
    alts = sorted({a.altloc for a in ats if a.altloc not in ('\x00', ' ', '')})
    alt = sorted(alts, key=lambda x: (-np.mean([a.occ for a in ats if a.altloc == x]), x))[0] if alts else '\x00'
    return [a for a in ats if a.altloc in (alt, '\x00', ' ', '')], alt


def load(pid):
    st = gemmi.read_structure(str(CIFS/f'{pid}.cif.gz'))
    st.setup_entities(); st.remove_hydrogens()
    return st


def subchains(st):
    """{subchain: dict(chain, kind, residues, atoms xyz, atom->residue index)} for model 1 (all altlocs, heavy atoms)."""
    out = {}
    for ch in st[0]:
        for r in ch:
            if r.name in SOLVENT:
                continue
            d = out.setdefault(r.subchain, dict(chain=ch.name, kind=r.entity_type.name, residues=[], atoms=[], rix=[], names=[]))
            d['residues'].append(r)
            for a in r:
                if not a.is_hydrogen():
                    d['atoms'].append([a.pos.x, a.pos.y, a.pos.z]); d['rix'].append(len(d['residues']) - 1); d['names'].append(a.name)
    for d in out.values():
        d['atoms'] = np.array(d['atoms'], float).reshape(-1, 3); d['rix'] = np.array(d['rix'], int)
    return out


def crystal_ops(st):
    """Cartesian (R, t, op index) of the space-group operators; translations are added per copy."""
    sg = gemmi.find_spacegroup_by_name(st.spacegroup_hm) if st.spacegroup_hm else None
    ops = list(sg.operations()) if sg else [gemmi.Op('x,y,z')]
    Fm, Fv = np.array(st.cell.frac.mat.tolist()), np.array(st.cell.frac.vec.tolist())
    Om, Ov = np.array(st.cell.orth.mat.tolist()), np.array(st.cell.orth.vec.tolist())
    out = []
    for k, op in enumerate(ops):
        S = np.array(op.float_seitz())
        out.append((S[:3, :3], S[:3, 3]))
    return out, (Fm, Fv, Om, Ov)


def copies_near(sub, L, ops, cellm, radius=SHELL):
    """All (subchain, op, shift) copies with an atom within radius of ligand atoms L. Returns list of dicts with the
    Cartesian transform and transformed atoms."""
    Fm, Fv, Om, Ov = cellm
    cen = L.mean(0)
    reach = np.sqrt(((L - cen) ** 2).sum(1)).max() + radius + 2.0
    cf = Fm @ cen + Fv
    found = []
    for name, d in sub.items():
        X = d['atoms']
        if not len(X):
            continue
        Xf = X @ Fm.T + Fv
        for k, (R, t) in enumerate(ops):
            Y = Xf @ R.T + t
            # per atom: lattice shift putting it nearest the ligand centre; only shifts of atoms that can reach the pose
            S = np.round(cf - Y)
            near = np.sqrt((((Y + S - cf) @ Om.T) ** 2).sum(1)) <= reach
            base = np.unique(S[near], axis=0) if near.any() else np.zeros((0, 3))
            cand = {tuple(int(v) for v in b + np.array(dv) - 1) for b in base for dv in np.ndindex(3, 3, 3)}
            if k == 0:
                cand.add((0, 0, 0))
            for s in cand:
                Yc = (Y + np.array(s)) @ Om.T + Ov
                close = np.sqrt(((Yc - cen) ** 2).sum(1)) <= reach
                dmin = np.full(len(Yc), np.inf)
                if close.any():
                    dmin[close] = np.sqrt(((Yc[close][:, None] - L[None]) ** 2).sum(-1)).min(1)
                if dmin.min() <= radius:
                    Rc = Om @ R @ Fm
                    tc = Om @ (R @ Fv + t + np.array(s)) + Ov
                    found.append(dict(subchain=name, op=k, shift=list(s), R=Rc, t=tc, X=Yc, dmin=dmin))
    # the same placement can be reached twice (equivalent op/shift); keep one
    uniq, seen = [], set()
    for f in found:
        key = (f['subchain'], tuple(np.round(f['X'][:3].ravel(), 2)))
        if key not in seen:
            seen.add(key); uniq.append(f)
    return uniq


def assembly_index(st, sub):
    """{assembly name: {subchain: [Cartesian (R, t)]}}; auth-chain-only generators are expanded to subchains."""
    out = {}
    for a in st.assemblies:
        m = defaultdict(list)
        for g in a.generators:
            names = list(g.subchains) or [s for s, d in sub.items() if d['chain'] in set(g.chains)]
            for o in g.operators:
                R, t = np.array(o.transform.mat.tolist()), np.array(o.transform.vec.tolist())
                for s in names:
                    m[s].append((R, t))
        out[a.name] = dict(members=m, details=a.oligomeric_details, author=bool(a.author_determined),
                           software=bool(a.software_determined))
    return out


def in_assembly(asm, anchor, copy, X0):
    """Is this copy of a subchain in the same assembly instance as the anchor's identity (ASU) copy?"""
    m = asm['members']
    if anchor not in m or copy['subchain'] not in m:
        return False
    P = X0[: min(len(X0), 60)]
    target = copy['X'][: len(P)]
    for Ra, ta in m[anchor]:
        Rinv = Ra.T
        for Rs, ts in m[copy['subchain']]:
            Y = ((P @ Rs.T + ts) - ta) @ Rinv.T   # g_anchor^-1 g_copy applied to ASU coordinates
            if np.abs(Y - target).max() < 1.0:
                return True
    return False


def residue_label(d, i):
    r = d['residues'][i]
    return f"{d['chain']}:{r.name}{r.seqid.num}{r.seqid.icode.strip()}"


def environment(st, L, anchor, tag):
    """Decompose everything within 6 A of pose L in structure st relative to the anchor subchain's ASU copy."""
    sub = subchains(st)
    ops, cellm = crystal_ops(st)
    asms = assembly_index(st, sub)
    copies = copies_near(sub, L, ops, cellm)
    recs = []
    for c in copies:
        d = sub[c['subchain']]
        within = c['dmin'] <= CONTACT
        lig_atoms = np.sqrt(((L[:, None] - c['X'][None]) ** 2).sum(-1)).min(1) <= CONTACT
        lig_atoms_cof = np.sqrt(((L[:, None] - c['X'][None]) ** 2).sum(-1)).min(1) <= COF_CONTACT
        res_contact = sorted({int(i) for i in d['rix'][within]})
        bb = np.array([n in ('N', 'CA', 'C', 'O') for n in d['names']])
        recs.append(dict(subchain=c['subchain'], chain=d['chain'], kind=d['kind'], op=c['op'], shift=c['shift'], heavy_atoms=len(d['atoms']),
                         copy='asu' if c['op'] == 0 and c['shift'] == [0, 0, 0] else 'symmetry',
                         components=sorted({d['residues'][i].name for i in res_contact}) if d['kind'] != 'Polymer' else [],
                         min_dist=float(c['dmin'].min()), contact_atoms=int(within.sum()),
                         ligand_atoms_contacted=int(lig_atoms.sum()), ligand_atoms_within_4A=int(lig_atoms_cof.sum()),
                         backbone_atoms_within_2p5=int(((c['dmin'] <= 2.5) & bb).sum()),
                         atoms_within_2p5=int((c['dmin'] <= 2.5).sum()),
                         contact_residues=[residue_label(d, i) for i in res_contact],
                         in_assembly={n: in_assembly(a, anchor, c, d['atoms']) for n, a in asms.items()}))
    return recs, {n: dict(details=a['details'], author=a['author'], software=a['software'], 
                          anchor_present=anchor in a['members']) for n, a in asms.items()}


def summarise(recs, nlig, asm_name):
    """Polymer contacts in vs outside the chosen assembly."""
    poly = [r for r in recs if r['kind'] == 'Polymer' and r['contact_atoms']]
    inside = [r for r in poly if r['in_assembly'].get(asm_name)]
    outside = [r for r in poly if not r['in_assembly'].get(asm_name)]
    tot = sum(r['contact_atoms'] for r in poly)
    out_atoms = sum(r['contact_atoms'] for r in outside)
    lig_out = max([r['ligand_atoms_contacted'] for r in outside] + [0])
    share = out_atoms / tot if tot else 0.0
    dep = 'none' if not outside else ('minor' if share < 0.15 and lig_out / nlig < 0.25 else 'major')
    return dict(assembly=asm_name, polymer_contact_atoms=tot, outside_contact_atoms=out_atoms, outside_share=share,
                outside_max_ligand_atom_fraction=lig_out / nlig, crystal_dependence=dep,
                inside_copies=[f"{r['chain']}/{r['subchain']}:{r['copy']}:op{r['op']}{tuple(r['shift'])}" for r in inside],
                outside_copies=[f"{r['chain']}/{r['subchain']}:{r['copy']}:op{r['op']}{tuple(r['shift'])}:{r['contact_atoms']}at" for r in outside],
                partner_chains_inside=sorted({r['chain'] for r in inside}))


def validation(pid):
    p = VAL/f'{pid}_validation.xml.gz'
    if not p.exists():
        return {}
    out = {}
    for _, el in ET.iterparse(gzip.open(p)):
        if el.tag == 'ModelledSubgroup' and el.get('resname'):
            out[(el.get('chain'), el.get('resname'), el.get('resnum'), (el.get('icode') or '').strip(), el.get('altcode') or '')] = \
                {k: el.get(k) for k in ('rscc', 'rsr', 'rsrz', 'ligRSRZ', 'avgoccu', 'ligand_num_clashes', 'ligand_clashes_outlier') if el.get(k) is not None}
        el.clear() if el.tag == 'ModelledSubgroup' else None
    return out


def fit(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, V = np.linalg.svd((P - pc).T @ (Q - qc))
    R = U @ np.diag([1, 1, np.linalg.det(U @ V)]) @ V
    return R, qc - pc @ R, float(np.sqrt((((P @ R + qc - pc @ R) - Q) ** 2).sum(1).mean()))


ALIGN = None


def entity_sequences(st):
    """{polymer subchain: one-letter SEQRES of its entity}."""
    out = {}
    for e in st.entities:
        if e.entity_type == gemmi.EntityType.Polymer:
            seq = gemmi.one_letter_code([gemmi.Entity.first_mon(m) for m in e.full_sequence])
            for s in e.subchains:
                out[s] = seq
    return out


def same_protein(a, b):
    global ALIGN
    if a == b:
        return True
    if ALIGN is None:
        from Bio import Align
        ALIGN = Align.PairwiseAligner(mode='global', match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-.1)
    aln = ALIGN.align(a, b)[0]
    ident = sum(a[i0 + k] == b[j0 + k] for (i0, i1), (j0, j1) in zip(*aln.aligned) for k in range(i1 - i0))
    return ident / max(1, min(len(a), len(b))) >= 0.95


def residue_by_key(st, chain, key):
    num, icode, name = key
    hits = [r for ch in st[0] if ch.name == chain for r in ch if r.seqid.num == num and r.seqid.icode.strip() == icode
            and r.name == name and r.entity_type == gemmi.EntityType.Polymer]
    return hits[0] if len(hits) == 1 else None


def audit(c, comp_of):
    hs, aps = load(c['holo']), load(c['apo'])
    comp, _, lnum = c['ligand_id'].split(':')
    lch = c['ligand_deposited_chain']
    ligs = [r for ch in hs[0] if ch.name == lch for r in ch if r.name == comp and r.seqid.num == int(lnum)
            and r.seqid.icode.strip() == c['ligand_icode']]
    assert len(ligs) == 1, 'ligand instance not unique in deposited holo'
    lig = ligs[0]
    la, alt = select_atoms(lig)
    L = xyz(la)
    out = dict(case_id=c['case_id'], apo=c['apo'], holo=c['holo'], ligand=c['ligand_id'], ligand_subchain=lig.subchain,
               ligand_altloc=alt.strip('\x00'), ligand_heavy_atoms=len(la), homology_component=comp_of[c['case_id']])
    # same-component instances and struct_conn
    out['same_component_instances'] = [f'{ch.name}:{r.seqid.num}' for ch in hs[0] for r in ch if r.name == comp]
    out['struct_conn'] = [dict(id=x.name, type=x.type.name, partner1=f'{x.partner1.chain_name}:{x.partner1.res_id.name}{x.partner1.res_id.seqid.num}:{x.partner1.atom_name}',
                               partner2=f'{x.partner2.chain_name}:{x.partner2.res_id.name}{x.partner2.res_id.seqid.num}:{x.partner2.atom_name}',
                               distance=x.reported_distance)
                          for x in hs.connections for p in (x.partner1, x.partner2)
                          if p.chain_name == lch and p.res_id.name == comp and p.res_id.seqid.num == int(lnum)]
    v = validation(c['holo'])
    out['validation'] = next((val for k, val in v.items() if k[0] == lch and k[1] == comp and k[2] == lnum), {})
    pocket_b = [a.b_iso for ch in hs[0] if ch.name == c['holo_deposited_chain'] for r in ch
                if r.entity_type == gemmi.EntityType.Polymer for a in r
                if not a.is_hydrogen() and np.sqrt(((L - [a.pos.x, a.pos.y, a.pos.z]) ** 2).sum(1)).min() <= CONTACT]
    out['ligand_mean_b'] = float(np.mean([a.b_iso for a in la])); out['pocket_mean_b'] = float(np.mean(pocket_b)) if pocket_b else None
    out['ligand_mean_occupancy'] = float(np.mean([a.occ for a in la]))
    # holo environment relative to the ligand's own subchain
    hrec, hasm = environment(hs, L, lig.subchain, 'holo')
    out['holo_assemblies'] = hasm
    receptor_sub = [r['subchain'] for r in hrec if r['chain'] == c['holo_deposited_chain'] and r['kind'] == 'Polymer' and r['copy'] == 'asu']
    out['holo_receptor_subchain'] = receptor_sub[0] if receptor_sub else None
    assert receptor_sub, 'receptor chain ASU copy does not contact the ligand'
    first = sorted(hasm, key=lambda n: (not hasm[n]['anchor_present'], n))[0]
    out['holo_assembly_used'] = first
    out['holo_receptor_in_assembly'] = any(r['in_assembly'].get(first) for r in hrec if r['subchain'] == receptor_sub[0] and r['copy'] == 'asu')
    out['holo_contacts'] = {n: summarise(hrec, len(L), n) for n in hasm}
    out['holo_records'] = hrec
    # non-polymer neighbours in the holo (cofactors, buffers, other ligands), any copy within the assembly or ASU
    nb = []
    for r in hrec:
        if r['kind'] == 'Polymer' or r['subchain'] == lig.subchain and r['copy'] == 'asu' or not r['contact_atoms']:
            continue
        nb.append(dict(components=r['components'], chain=r['chain'], subchain=r['subchain'], copy=r['copy'],
                       in_assembly=r['in_assembly'].get(first), heavy_atoms=r['heavy_atoms'], contact_atoms=r['contact_atoms'],
                       ligand_atom_fraction_4A=r['ligand_atoms_within_4A'] / len(L), min_dist=r['min_dist']))
    out['holo_nonpolymer_neighbours'] = nb
    # move the pose into the deposited apo frame through the receptor chain's pocket backbone
    P, Q, used = [], [], 0
    for m in c['residue_map']:
        rh = residue_by_key(hs, c['holo_deposited_chain'], tuple(m['holo_key']))
        ra = residue_by_key(aps, c['apo_deposited_chain'], tuple(m['apo_key']))
        if rh is None or ra is None:
            continue
        ah, _ = select_atoms(rh); aa, _ = select_atoms(ra)
        ah = {a.name: a for a in ah}; aa = {a.name: a for a in aa}
        if all(n in ah and n in aa for n in ('N', 'CA', 'C', 'O')):
            used += 1
            for n in ('N', 'CA', 'C', 'O'):
                P.append([ah[n].pos.x, ah[n].pos.y, ah[n].pos.z]); Q.append([aa[n].pos.x, aa[n].pos.y, aa[n].pos.z])
    R, T, rmsd = fit(np.array(P), np.array(Q))
    X = L @ R + T
    Xq = L @ np.array(c['rotation_holo_to_apo']) + np.array(c['translation_holo_to_apo'])
    out.update(apo_fit_residues=used, apo_fit_mapped=len(c['residue_map']), deposited_pocket_backbone_rmsd=rmsd,
               qfit_vs_deposited_transfer_rmsd=float(np.sqrt(((X - Xq) ** 2).sum(1).mean())))
    apo_anchor = [s for s, d in subchains(aps).items() if d['chain'] == c['apo_deposited_chain'] and d['kind'] == 'Polymer']
    assert len(apo_anchor) == 1, 'apo receptor polymer subchain not unique'
    arec, aasm = environment(aps, X, apo_anchor[0], 'apo')
    afirst = sorted(aasm, key=lambda n: (not aasm[n]['anchor_present'], n))[0]
    out['apo_assemblies'] = aasm; out['apo_assembly_used'] = afirst
    out['apo_contacts'] = {n: summarise(arec, len(X), n) for n in aasm}
    out['apo_records'] = arec
    out['apo_nonpolymer_at_site'] = [dict(components=r['components'], chain=r['chain'], copy=r['copy'], min_dist=r['min_dist'],
                                          contact_atoms=r['contact_atoms']) for r in arec if r['kind'] != 'Polymer' and r['contact_atoms']]
    # cofactor dependence: holo neighbour with no same component at the apo site
    # same component in the apo near the transferred pose (within 3 A of a ligand atom)
    apo_comps = {x for r in arec if r['kind'] != 'Polymer' and r['min_dist'] <= 3.0 for x in r['components']}
    dep = []
    for n in nb:
        listed = bool(set(n['components']) & COFACTORS)
        if n['ligand_atom_fraction_4A'] >= 0.10 and (listed or n['heavy_atoms'] >= 6):
            dep.append(dict(components=n['components'], chain=n['chain'], copy=n['copy'], in_assembly=n['in_assembly'],
                            ligand_atom_fraction_4A=n['ligand_atom_fraction_4A'], heavy_atoms=n['heavy_atoms'],
                            listed_cofactor=listed, present_at_apo_site=bool(set(n['components']) & apo_comps)))
    out['holo_neighbour_dependence'] = dep
    # partner chains needed in the holo assembly and their presence at the apo site
    # partner identity by deposited entity full sequence (SEQRES), identity >= 0.95 over the shorter sequence
    hs_seq, as_seq = entity_sequences(hs), entity_sequences(aps)
    partners = []
    for r in hrec:
        if r['kind'] != 'Polymer' or not r['contact_atoms'] or not r['in_assembly'].get(first) or (r['subchain'] == receptor_sub[0] and r['copy'] == 'asu'):
            continue
        apo_hits = [x for x in arec if x['kind'] == 'Polymer' and x['contact_atoms'] and x['in_assembly'].get(afirst)
                    and not (x['subchain'] == apo_anchor[0] and x['copy'] == 'asu') and same_protein(as_seq[x['subchain']], hs_seq[r['subchain']])]
        partners.append(dict(holo_copy=f"{r['chain']}/{r['subchain']}:{r['copy']}:op{r['op']}{tuple(r['shift'])}", contact_atoms=r['contact_atoms'],
                             ligand_atoms_contacted=r['ligand_atoms_contacted'],
                             same_sequence_as_receptor=same_protein(hs_seq[r['subchain']], hs_seq[receptor_sub[0]]),
                             apo_in_assembly_same_sequence_copies=[f"{x['chain']}/{x['subchain']}:{x['copy']}:{x['contact_atoms']}at" for x in apo_hits]))
    out['holo_partner_chains'] = partners
    return out


def main():
    C = json.loads(CANDS.read_text())
    comp_of = {r['case_id']: r['homology_component'] for r in json.loads(RECORDS.read_text())}
    res, errors = [], []
    for k, c in enumerate(C):
        try:
            res.append(audit(c, comp_of))
        except Exception as e:
            errors.append(dict(case_id=c['case_id'], error=repr(e))); print('ERROR', c['case_id'], repr(e), flush=True)
        print(k, c['case_id'], flush=True)
    inputs = {str(p.relative_to(PILOT)): sha(p) for p in [CANDS, RECORDS, *sorted(CIFS.glob('*.cif.gz')), *sorted(VAL.glob('*.xml.gz'))]}
    (HERE/'assembly_contact_audit.json').write_text(json.dumps(dict(cases=res, errors=errors, inputs_sha256=inputs,
                                                                    code_sha256=sha(__file__), gemmi=gemmi.__version__), indent=1, default=float) + '\n')
    print('done', len(res), 'errors', len(errors))


if __name__ == '__main__':
    main()
