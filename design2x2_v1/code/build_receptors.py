"""Build the 2x2 receptors of the representatives from qFit-at-scale models (no packer, no docking).

  arms   holo_single, holo_multi, apo_single, apo_multi on the qFit holo / apo model of each case
  chains whole in-assembly polymer copies within 6 A of the ligand reference pose (deposited assembly membership from
         final_selection_v1; symmetry copies made with the deposited operators applied to qFit coordinates). All
         HETATM and waters removed.
  backbone per chain, the altloc path of highest occupancy with intact peptide junctions (backbone_path); side chains
         grafted on it. single every residue's side chain in its major altloc (highest mean occupancy, shared atoms kept); background (non-pocket)
         side chains that clash are switched to their next altloc on the major backbone, once per backbone.
  pocket residues of these copies with any heavy atom (any altloc) within 4.5 A of the pose (holo: deposited pose;
         apo: whole-pocket refit pose). Each copy instance is a separate residue.
  states per pocket residue the occupied chi1/chi2 wells (>= 0.10, 105_joint_rotamers); a well is represented by its
         highest-occupancy altloc; the major state is the well of largest summed occupancy.
  graft  a state's side chain is moved onto the major backbone by the rigid fit of its own N/CA/C onto the major
         N/CA/C (local frame), keeping its refined internal geometry; the major backbone is never changed.
  single arm: major state on every pocket residue. multi arm: K = 8 receptors; all combinations if <= 8, else the
         all-major combination plus combinations drawn with per-residue well occupancies (seed 20261028), rejecting
         any receptor with a heavy-atom contact < 2.2 A between a grafted side chain and any other residue. Fewer than K valid unique combinations are kept as found.
Writes receptors/<case>/<arm>_<k>.pdb, ligand_ref_{holo,apo}.pdb and build_manifest.json.
"""
import csv, hashlib, importlib.util, json, sys
from itertools import product
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
K, SEED, CLASH, POCKET = 8, 20261028, 2.2, 4.5
BB = ('N', 'CA', 'C', 'O')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m)
    return m


A = load('aac', PILOT/'final_selection_v1/code/audit_assembly_contacts.py')
S = load('aqs', PILOT/'state_audit_v1/code/audit_qfit_states.py')
j105 = S.j105


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pos(a): return np.array([a.pos.x, a.pos.y, a.pos.z])


def altloc_atoms(res, alt):
    """{name: atom} for one altloc (shared atoms included)."""
    d = {a.name: a for a in res if not a.is_hydrogen() and a.altloc in ('\x00', ' ')}
    d.update({a.name: a for a in res if not a.is_hydrogen() and a.altloc == alt})
    return d


def major_alt(res):
    return A.select_atoms(res)[1]


def wells_with_rep(res):
    """[(well key, summed occupancy, representative altloc)] of the occupied chi1/chi2 wells, major first."""
    states = S.residue_states(res)
    s = [(alt, o, c) for alt, o, c in states if c]
    if not s:
        return []
    tot = sum(o for _, o, _ in s) or 1.0
    acc = {}
    for alt, o, c in s:
        k = j105.well3(c[0]) + (j105.chi2_well(res.name, c[1]) if len(c) > 1 and j105.chi2_class(res.name) else '')
        w = acc.setdefault(k, [0.0, None, -1.0])
        w[0] += o / tot
        if o > w[2]:
            w[1], w[2] = alt, o
    out = [(k, v[0], v[1]) for k, v in acc.items() if v[0] >= j105.FLOOR - j105.TIE]
    return sorted(out, key=lambda x: -x[1])


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, V = np.linalg.svd((P - pc).T @ (Q - qc))
    R = U @ np.diag([1, 1, np.linalg.det(U @ V)]) @ V
    return R, qc - pc @ R


def residue_coords(res, state_alt, Rc, tc, bb_alt=None):
    """[(name, element, xyz)] of the residue: backbone of altloc bb_alt (default: its major altloc), side chain of
    state_alt (default: its major altloc) grafted onto that backbone by the N/CA/C local frame."""
    bbk = major_alt(res) if bb_alt is None else bb_alt
    bba = altloc_atoms(res, bbk)
    out = [(n, bba[n].element.name, pos(bba[n])) for n in BB if n in bba]
    src = altloc_atoms(res, major_alt(res) if state_alt is None else state_alt)
    sc = [n for n in src if n not in BB]
    if all(n in src and n in bba for n in ('N', 'CA', 'C')):
        R, t = kabsch(np.array([pos(src[n]) for n in ('N', 'CA', 'C')]), np.array([pos(bba[n]) for n in ('N', 'CA', 'C')]))
    else:
        R, t = np.eye(3), np.zeros(3)
    out += [(n, src[n].element.name, pos(src[n]) @ R + t) for n in sc]
    return [(n, e, x @ Rc.T + tc) for n, e, x in out]


def backbone_path(poly):
    """Backbone altloc per residue along a chain: highest summed occupancy subject to intact peptide junctions
    (C(i)-N(i+1) within 1.2-1.5 A for sequence neighbours), by dynamic programming; a junction that no altloc pair
    can close costs the same for every choice. Per-residue majors can come from conformers that do not join."""
    opts = []
    for r in poly:
        alts = sorted({a.altloc for a in r if a.altloc not in ('\x00', ' ') and a.name in BB}) or [major_alt(r)]
        occ = {x: (np.mean([a.occ for a in r if a.altloc == x and a.name in BB]) if x not in ('\x00', ' ', '') else 1.0) for x in alts}
        opts.append([(x, occ[x], altloc_atoms(r, x)) for x in alts])
    def broken(i, a, b):
        if poly[i + 1].seqid.num - poly[i].seqid.num != 1 or 'C' not in a or 'N' not in b:
            return 0.0
        return 0.0 if 1.2 <= np.linalg.norm(pos(a['C']) - pos(b['N'])) <= 1.5 else 1.0
    cost = [[-o for _, o, _ in opts[0]]]; back = [[None] * len(opts[0])]
    for i in range(1, len(poly)):
        row, brow = [], []
        for x, o, at in opts[i]:
            c = [cost[i - 1][j] + 100.0 * broken(i - 1, opts[i - 1][j][2], at) for j in range(len(opts[i - 1]))]
            j = int(np.argmin(c)); row.append(c[j] - o); brow.append(j)
        cost.append(row); back.append(brow)
    j = int(np.argmin(cost[-1])); path = [None] * len(poly)
    for i in range(len(poly) - 1, -1, -1):
        path[i] = opts[i][j][0]; j = back[i][j] if back[i][j] is not None else 0
    return path


def chis(name, coords):
    d = {n: x for n, _, x in coords}
    out = []
    for q in S.CHI.get(name, []):
        if not all(n in d for n in q):
            break
        out.append(float(np.degrees(gemmi.calculate_dihedral(*[gemmi.Position(*d[n]) for n in q]))))
    return out


def deposited_copies(o, side, st_dep):
    """Cartesian (chain, R, t, copy label) of in-assembly polymer copies within 6 A, from the deposited audit."""
    ops, cellm = A.crystal_ops(st_dep)
    Fm, Fv, Om, Ov = cellm
    used = o[f'{side}_assembly_used']
    out, seen = [], set()
    for r in o[f'{side}_records']:
        if r['kind'] != 'Polymer' or not r['in_assembly'].get(used):
            continue
        Rf, tf = ops[r['op']]
        Rc = Om @ Rf @ Fm
        tc = Om @ (Rf @ Fv + tf + np.array(r['shift'])) + Ov
        key = (r['chain'], r['op'], tuple(r['shift']))
        if key not in seen:
            seen.add(key); out.append((r['chain'], Rc, tc, r['copy']))
    return out


def build_side(pid, copies, L, outdir, tag, rng, chain_map):
    st = S.qfit(pid)
    residues = []          # (copy index, chain label, residue, Rc, tc)
    paths = {}
    for k, (dch, Rc, tc, kind) in enumerate(copies):
        qch = chain_map.get(dch, dch)
        polys = [ch.get_polymer() for ch in st[0] if ch.name == qch and len(ch.get_polymer())]
        assert len(polys) == 1, (pid, qch)
        if qch not in paths:
            paths[qch] = backbone_path(polys[0])
        for res, bb in zip(polys[0], paths[qch]):
            residues.append((k, qch, res, Rc, tc, bb))
    # pocket: any altloc heavy atom within 4.5 A of the pose
    pocket = []
    for idx, (k, qch, res, Rc, tc, _) in enumerate(residues):
        X = np.array([pos(a) for a in res if not a.is_hydrogen()]) @ Rc.T + tc
        if np.sqrt(((X[:, None] - L[None]) ** 2).sum(-1)).min() <= POCKET:
            pocket.append(idx)
    states = {i: wells_with_rep(residues[i][2]) for i in pocket}
    multi = [i for i in pocket if len(states[i]) >= 2]
    rc = lambda i, alt: residue_coords(residues[i][2], alt, residues[i][3], residues[i][4], residues[i][5])
    base = {i: rc(i, None) for i in range(len(residues))}

    def coords_for(combo):
        c = dict(base)
        for i in pocket:
            st_list = states[i]
            alt = st_list[combo.get(i, 0)][2] if st_list else None
            c[i] = rc(i, alt)
        return c

    def clashes(c, changed):
        """Heavy-atom contacts < CLASH between grafted side chains and any other residue."""
        allX, owner, isbb = [], [], []
        for i, at in c.items():
            for n, _, x in at:
                allX.append(x); owner.append(i); isbb.append(n in BB)
        allX, owner, isbb = np.array(allX), np.array(owner), np.array(isbb)
        worst = []
        for i in changed:
            sc = np.array([x for n, _, x in c[i] if n not in BB])
            if not len(sc):
                continue
            valid = owner != i   # neighbours included: side chain to neighbour backbone is > 2.2 A in normal geometry
            d = np.sqrt(((sc[:, None] - allX[valid][None]) ** 2).sum(-1))
            if (d < CLASH).any():
                worst.append(i)
        return worst

    # Background (non-pocket) side chains: per-residue major altlocs need not coexist in the qFit model, so a clashing
    # background side chain is switched to its next altloc (grafted on the major backbone), lowest major occupancy
    # first, until no background-involved clash remains or no alternative is left. Fixed once per backbone, so every
    # arm of the backbone shares it.
    pocket_set = set(pocket)
    def alts_by_occ(res):
        a = sorted({x.altloc for x in res if x.altloc not in ('\x00', ' ')})
        return sorted(a, key=lambda z: (-np.mean([x.occ for x in res if x.altloc == z]), z))
    switched, tried = {}, {}
    for _ in range(500):
        bad = [i for i in clashes(coords_for({}), range(len(residues))) if i not in pocket_set]
        cand = [i for i in bad if len(alts_by_occ(residues[i][2])) > 1 + tried.get(i, 0)]
        if not cand:
            break
        i = min(cand, key=lambda j: (np.mean([x.occ for x in residues[j][2] if x.altloc == major_alt(residues[j][2])] or [1.0]), j))
        tried[i] = tried.get(i, 0) + 1
        alt = alts_by_occ(residues[i][2])[tried[i]]
        switched[i] = alt
        base[i] = rc(i, alt)
    unresolved = clashes(coords_for({}), range(len(residues)))
    # disulfides are bonds, not clashes
    def ss_only(i):
        c = coords_for({})
        sg = [x for n, _, x in c[i] if n == 'SG']
        if not sg:
            return False
        others = [x for j, at in c.items() if j != i for n, _, x in at]
        near = [x for x in others if np.linalg.norm(x - sg[0]) < CLASH]
        return all(any(n == 'SG' and np.linalg.norm(x - y) < 1e-6 for j, at in c.items() if j != i for n, _, y in at) for x in near)
    unresolved = [i for i in unresolved if not (residues[i][2].name == 'CYS' and ss_only(i))]
    # background residues still clashing: side chain beyond CB removed (rebuilt later by PDB2PQR), in every arm
    stripped = [i for i in unresolved if i not in pocket_set and residues[i][2].name not in ('GLY', 'ALA')]
    for i in stripped:
        base[i] = [(n, e, x) for n, e, x in base[i] if n in BB or n == 'CB']
    unresolved = [i for i in unresolved if i not in stripped]

    # chi deviation of every grafted state relative to its source altloc
    dev = 0.0
    for i in pocket:
        for _, _, alt in states[i]:
            src = altloc_atoms(residues[i][2], alt)
            ref = chis(residues[i][2].name, [(n, '', pos(a)) for n, a in src.items()])
            got = chis(residues[i][2].name, residue_coords(residues[i][2], alt, np.eye(3), np.zeros(3), residues[i][5]))
            for x, y in zip(ref, got):
                dev = max(dev, abs((x - y + 180) % 360 - 180))
    combos = [dict()]                                   # all-major
    space = [range(len(states[i])) for i in multi]
    total = int(np.prod([len(s) for s in space])) if multi else 1
    if multi and total <= K:
        combos = [dict(zip(multi, c)) for c in product(*space)]
    elif multi:
        seen, tries = {tuple([0] * len(multi))}, 0
        while len(combos) < K and tries < 5000:
            tries += 1
            pick = tuple(int(rng.choice(len(states[i]), p=np.array([w[1] for w in states[i]]) / sum(w[1] for w in states[i]))) for i in multi)
            if pick in seen:
                continue
            seen.add(pick)
            combos.append(dict(zip(multi, pick)))
    kept, rejected = [], 0
    for combo in combos:
        c = coords_for(combo)
        bad = clashes(c, [i for i, s in combo.items() if s != 0])
        if bad and combo:
            rejected += 1; continue
        kept.append((combo, c))
    names = {}
    def write(c, path):
        s2 = gemmi.Structure(); m = gemmi.Model('1')
        chains = {}
        for i, at in c.items():
            k, qch, res = residues[i][:3]
            label = f'{qch}{k}'
            if label not in chains:
                letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
                chains[label] = gemmi.Chain(letters[len(chains)]); names[label] = chains[label].name
            r = gemmi.Residue(); r.name = res.name; r.seqid = res.seqid; r.het_flag = 'A'
            for n, e, x in at:
                a = gemmi.Atom(); a.name = n; a.element = gemmi.Element(e); a.pos = gemmi.Position(*x); a.occ = 1.0; a.b_iso = 20.0
                r.add_atom(a)
            chains[label].add_residue(r)
        for ch in chains.values():
            m.add_chain(ch)
        s2.add_model(m); s2.setup_entities(); s2.write_pdb(str(path))
    write(coords_for({}), outdir/f'{tag}_single_0.pdb')
    for n, (combo, c) in enumerate(kept):
        write(c, outdir/f'{tag}_multi_{n}.pdb')
    label = lambda i: f"{residues[i][1]}{residues[i][0]}:{residues[i][2].name}{residues[i][2].seqid.num}{residues[i][2].seqid.icode.strip()}"
    single_bad = [label(i) for i in clashes(coords_for({}), pocket)]   # the deposited-major state, reported only
    return dict(chains=[f'{qch}{k}' for k, (dch, _, _, kind) in enumerate(copies) for qch in [chain_map.get(dch, dch)]],
                chain_names_in_pdb=names, pocket_residues=[label(i) for i in pocket],
                multi_state_residues={label(i): [(w[0], round(w[1], 3), w[2]) for w in states[i]] for i in multi},
                combinations_total=total, multi_receptors=len(kept), rejected_for_clash=rejected,
                combos=[{label(i): states[i][s][0] for i, s in combo.items()} for combo, _ in kept],
                max_graft_chi_deviation_deg=round(dev, 2), single_clashes=single_bad,
                background_switched={label(i): a for i, a in switched.items()},
                background_stripped=[label(i) for i in stripped], residual_clashes=[label(i) for i in unresolved])


def write_ligand(lig_atoms, X, path):
    s2 = gemmi.Structure(); m = gemmi.Model('1'); ch = gemmi.Chain('L'); r = gemmi.Residue()
    r.name = lig_atoms[0][0]; r.seqid = gemmi.SeqId(1, ' '); r.het_flag = 'H'
    for (_, n, e), x in zip(lig_atoms, X):
        a = gemmi.Atom(); a.name = n; a.element = gemmi.Element(e); a.pos = gemmi.Position(*x); a.occ = 1.0
        r.add_atom(a)
    ch.add_residue(r); m.add_chain(ch); s2.add_model(m); s2.write_pdb(str(path))


def main():
    C = {c['case_id']: c for c in json.loads((PILOT/'final_review_v1/candidates79.json').read_text())}
    O = {o['case_id']: o for o in json.loads((PILOT/'final_selection_v1/assembly_contact_audit.json').read_text())['cases']}
    F = {o['case_id']: o for o in json.loads((PILOT/'state_audit_v1/pose_refit.json').read_text())['cases']}
    reps = list(csv.DictReader(open(HERE/'representatives.tsv'), delimiter='\t'))
    manifest = []
    for rep in reps:
        cid = rep['case_id']; c, o, f = C[cid], O[cid], F[cid]
        rng = np.random.default_rng(SEED)
        outdir = HERE/'receptors'/cid; outdir.mkdir(parents=True, exist_ok=True)
        hs, aps = A.load(c['holo']), A.load(c['apo'])
        comp, _, lnum = c['ligand_id'].split(':')
        lig = [r for ch in hs[0] if ch.name == c['ligand_deposited_chain'] for r in ch if r.name == comp and r.seqid.num == int(lnum)][0]
        la, _ = A.select_atoms(lig)
        L = A.xyz(la); X = np.array(f['apo_pose'])
        atoms = [(comp, a.name, a.element.name) for a in la]
        write_ligand(atoms, L, outdir/'ligand_ref_holo.pdb'); write_ligand(atoms, X, outdir/'ligand_ref_apo.pdb')
        # qFit receptor chains equal deposited names except the recorded receptor aliases
        hmap = {c['holo_deposited_chain']: c['holo_chain']}; amap = {c['apo_deposited_chain']: c['apo_chain']}
        # apo copies: re-derive the in-assembly copies around the refit pose
        arec, aasm = A.environment(aps, X, [s for s, d in A.subchains(aps).items() if d['chain'] == c['apo_deposited_chain'] and d['kind'] == 'Polymer'][0], 'apo')
        o_apo = dict(apo_records=arec, apo_assembly_used=sorted(aasm, key=lambda n: (not aasm[n]['anchor_present'], n))[0])
        holo = build_side(c['holo'], deposited_copies(o, 'holo', hs), L, outdir, 'holo', rng, hmap)
        apo = build_side(c['apo'], deposited_copies(o_apo, 'apo', aps), X, outdir, 'apo', rng, amap)
        manifest.append(dict(rep, holo=holo, apo=apo))
        print(cid, rep['stratum'], 'holo multi', len(holo['multi_state_residues']), 'recs', holo['multi_receptors'], 'rej', holo['rejected_for_clash'],
              '| apo multi', len(apo['multi_state_residues']), 'recs', apo['multi_receptors'], 'rej', apo['rejected_for_clash'],
              '| chi dev', holo['max_graft_chi_deviation_deg'], apo['max_graft_chi_deviation_deg'], flush=True)
    (HERE/'build_manifest.json').write_text(json.dumps(dict(cases=manifest, K=K, seed=SEED, clash=CLASH, pocket=POCKET,
                                                            code_sha256=sha(__file__)), indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
