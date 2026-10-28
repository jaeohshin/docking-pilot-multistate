"""Re-place the ligand reference pose in the deposited apo frame with every pocket chain, then redo the apo pocket
and its qFit states.

The first transfer (final_selection_v1) fitted only the receptor chain's mapped residues; for interface sites that
can be 4-7 residues of the minor side. Here every in-assembly holo polymer copy within 6 A of the ligand is paired
with the in-assembly apo copy of the same protein lying nearest to it under the current transfer; residues are
matched by a global alignment of the two copies' observed sequences (numbering may differ between entries), and
N/CA/C/O of all holo shell residues (6 A) and their apo partners are fitted. Two passes (start: the first transfer).
Writes pose_refit.json: per case the fit, the ligand displacement from the first transfer, and the apo pocket with
its qFit states under the new pose. No docking, inference or training.
"""
import hashlib, importlib.util, json, sys
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m)
    return m


A = load('aac', PILOT/'final_selection_v1/code/audit_assembly_contacts.py')
S = load('aqs', HERE/'code/audit_qfit_states.py')
from Bio import Align  # noqa: E402
ALN = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-3, extend_gap_score=-.5)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def copies(st, L, anchor, radius):
    """In-assembly polymer copies within radius of L (first assembly holding the anchor)."""
    sub = A.subchains(st)
    ops, cellm = A.crystal_ops(st)
    asms = A.assembly_index(st, sub)
    name = sorted(asms, key=lambda n: (anchor not in asms[n]['members'], n))[0]
    out = []
    for c in A.copies_near(sub, L, ops, cellm, radius):
        d = sub[c['subchain']]
        if d['kind'] == 'Polymer' and A.in_assembly(asms[name], anchor, c, d['atoms']):
            out.append(dict(subchain=c['subchain'], chain=d['chain'], residues=d['residues'], R=c['R'], t=c['t'],
                            dmin=c['dmin'], rix=d['rix']))
    return out


def backbone(res, R, t):
    atoms, _ = A.select_atoms(res)
    return {a.name: np.array([a.pos.x, a.pos.y, a.pos.z]) @ R.T + t for a in atoms if a.name in ('N', 'CA', 'C', 'O')}


def index_map(a, b):
    aln = ALN.align(a, b)[0]
    return {i0 + k: j0 + k for (i0, i1), (j0, j1) in zip(*aln.aligned) for k in range(i1 - i0) if a[i0 + k] == b[j0 + k]}


def one_letter(residues): return gemmi.one_letter_code([r.name for r in residues])


MIN_COPY = 3  # copies with fewer shell residues are recorded but not fitted (a 1-2 residue contact can mis-pair)


def collect(hcopies, acopies, hseq, aseq, R, T, receptor):
    """Holo/apo backbone pairs over the holo shell residues of every paired copy, under transfer (R, T).
    The holo receptor chain's ASU copy is paired with the apo receptor chain's ASU copy (the source pairing)."""
    P, Q, info = [], [], []
    for h in hcopies:
        shell = sorted({int(i) for i in h['rix'][h['dmin'] <= A.SHELL]})
        if len(shell) < MIN_COPY:
            info.append(dict(holo_copy=h['chain'], shell_residues=len(shell), fitted_residues=0, note='below MIN_COPY')); continue
        hb = {i: backbone(h['residues'][i], h['R'], h['t']) for i in shell}
        hca = np.array([v['CA'] @ R + T for v in hb.values() if 'CA' in v])
        if not len(hca):
            continue
        is_receptor = h['subchain'] == receptor[0] and np.allclose(h['R'], np.eye(3)) and np.allclose(h['t'], 0)
        best = None
        for a in acopies:
            if is_receptor and not (a['subchain'] == receptor[1] and np.allclose(a['R'], np.eye(3)) and np.allclose(a['t'], 0)):
                continue
            if not A.same_protein(aseq[a['subchain']], hseq[h['subchain']]):
                continue
            aca = np.array([[x.pos.x, x.pos.y, x.pos.z] for r in a['residues'] for x in r if x.name == 'CA']) @ a['R'].T + a['t']
            d = float(np.mean(np.sqrt(((hca[:, None] - aca[None]) ** 2).sum(-1)).min(1)))
            if best is None or d < best[0]:
                best = (d, a)
        if best is None:
            info.append(dict(holo_copy=h['chain'], paired=None)); continue
        a = best[1]
        m = index_map(one_letter(h['residues']), one_letter(a['residues']))
        n, start = 0, len(P)
        for i in shell:
            if i not in m:
                continue
            ab = backbone(a['residues'][m[i]], a['R'], a['t'])
            common = [k for k in ('N', 'CA', 'C', 'O') if k in hb[i] and k in ab]
            if len(common) == 4:
                n += 1
                P += [hb[i][k] for k in common]; Q += [ab[k] for k in common]
        info.append(dict(holo_copy=h['chain'], apo_copy=a['chain'], receptor_pair=bool(is_receptor), pre_fit_mean_ca_gap=round(best[0], 3),
                         shell_residues=len(shell), fitted_residues=n, rows=[start, len(P)]))
    return np.array(P), np.array(Q), info


def case(c, o):
    hs, aps = A.load(c['holo']), A.load(c['apo'])
    comp, _, lnum = c['ligand_id'].split(':')
    lig = [r for ch in hs[0] if ch.name == c['ligand_deposited_chain'] for r in ch if r.name == comp and r.seqid.num == int(lnum)
           and r.seqid.icode.strip() == c['ligand_icode']][0]
    la, _ = A.select_atoms(lig)
    L = A.xyz(la)
    # first transfer, as in final_selection_v1 (receptor chain residue_map)
    P0, Q0 = [], []
    for m in c['residue_map']:
        rh = A.residue_by_key(hs, c['holo_deposited_chain'], tuple(m['holo_key']))
        ra = A.residue_by_key(aps, c['apo_deposited_chain'], tuple(m['apo_key']))
        ah = {a.name: a for a in A.select_atoms(rh)[0]}; aa = {a.name: a for a in A.select_atoms(ra)[0]}
        for n in ('N', 'CA', 'C', 'O'):
            P0.append([ah[n].pos.x, ah[n].pos.y, ah[n].pos.z]); Q0.append([aa[n].pos.x, aa[n].pos.y, aa[n].pos.z])
    R, T, rmsd0 = A.fit(np.array(P0), np.array(Q0))
    X0 = L @ R + T
    anchor = [s for s, d in A.subchains(aps).items() if d['chain'] == c['apo_deposited_chain'] and d['kind'] == 'Polymer'][0]
    hseq, aseq = A.entity_sequences(hs), A.entity_sequences(aps)
    hcopies = copies(hs, L, lig.subchain, A.SHELL)
    acopies = copies(aps, X0, anchor, 12.0)
    hrec = [s for s, d in A.subchains(hs).items() if d['chain'] == c['holo_deposited_chain'] and d['kind'] == 'Polymer'][0]
    for _ in range(2):
        P, Q, info = collect(hcopies, acopies, hseq, aseq, R, T, (hrec, anchor))
        R, T, rmsd = A.fit(P, Q)
    for x in info:                               # per-copy backbone RMSD under the joint fit
        if x.get('rows') and x['rows'][1] > x['rows'][0]:
            a0, a1 = x.pop('rows')
            x['rmsd_in_joint_fit'] = float(np.sqrt((((P[a0:a1] @ R + T) - Q[a0:a1]) ** 2).sum(1).mean()))
        x.pop('rows', None)
    X = L @ R + T
    arec, aasm = A.environment(aps, X, anchor, 'apo')
    afirst = sorted(aasm, key=lambda n: (not aasm[n]['anchor_present'], n))[0]
    rows, miss = S.audit_side(c['apo'], arec, afirst, {c['apo_deposited_chain']: c['apo_chain']})
    return dict(case_id=c['case_id'], first_fit_residues=len(c['residue_map']), first_fit_rmsd=rmsd0,
                refit_residues=int(len(P) // 4), refit_rmsd=rmsd, chains=info,
                ligand_shift_mean=float(np.sqrt(((X - X0) ** 2).sum(1)).mean()), ligand_shift_max=float(np.sqrt(((X - X0) ** 2).sum(1)).max()),
                apo_pose=X.tolist(), apo_contacts=A.summarise(arec, len(X), afirst), apo_states=S.side_summary(rows),
                apo_residues=rows, apo_missing_in_qfit=miss)


def main():
    C = json.loads((PILOT/'final_review_v1/candidates79.json').read_text())
    O = {o['case_id']: o for o in json.loads((PILOT/'final_selection_v1/assembly_contact_audit.json').read_text())['cases']}
    out, err = [], []
    for c in C:
        try:
            r = case(c, O[c['case_id']]); out.append(r)
            print(c['case_id'], 'fit', r['first_fit_residues'], '->', r['refit_residues'], 'rmsd %.2f' % r['refit_rmsd'],
                  'shift %.2f/%.2f' % (r['ligand_shift_mean'], r['ligand_shift_max']), flush=True)
        except Exception as e:
            err.append(dict(case_id=c['case_id'], error=repr(e))); print('ERROR', c['case_id'], repr(e), flush=True)
    (HERE/'pose_refit.json').write_text(json.dumps(dict(cases=out, errors=err, code_sha256=sha(__file__)), indent=1, default=float) + '\n')
    print('done', len(out), 'errors', len(err))


if __name__ == '__main__':
    main()
