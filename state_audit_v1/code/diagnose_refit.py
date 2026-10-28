"""Per-chain diagnosis of high-RMSD multi-chain refits: fit each paired holo/apo copy alone and report the other
copies' backbone RMSD under that fit (quaternary difference at the site vs. a pairing error)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import refit_pose as RP
A = RP.A
PILOT = RP.PILOT
C = {c['case_id']: c for c in json.loads((PILOT/'final_review_v1/candidates79.json').read_text())}
ids = sys.argv[1:]
for cid in ids:
    c = C[cid]
    hs, aps = A.load(c['holo']), A.load(c['apo'])
    comp, _, lnum = c['ligand_id'].split(':')
    lig = [r for ch in hs[0] if ch.name == c['ligand_deposited_chain'] for r in ch if r.name == comp and r.seqid.num == int(lnum)][0]
    L = A.xyz(A.select_atoms(lig)[0])
    hseq, aseq = A.entity_sequences(hs), A.entity_sequences(aps)
    hc = RP.copies(hs, L, lig.subchain, A.SHELL)
    anchor = [s for s, d in A.subchains(aps).items() if d['chain'] == c['apo_deposited_chain'] and d['kind'] == 'Polymer'][0]
    # pair every holo copy with every same-protein apo copy near the site and record per-pair fits
    # every apo ASU polymer subchain (identity copy); frames differ between entries, so no distance filter
    ac = [dict(subchain=k, chain=d['chain'], residues=d['residues'], R=np.eye(3), t=np.zeros(3))
          for k, d in A.subchains(aps).items() if d['kind'] == 'Polymer']
    print('==', cid, c['holo'], c['ligand_id'], 'holo copies', [h['chain'] for h in hc], 'apo copies', [a['chain'] for a in ac])
    for h in hc:
        shell = sorted({int(i) for i in h['rix'][h['dmin'] <= A.SHELL]})
        res = []
        for a in ac:
            if not A.same_protein(aseq[a['subchain']], hseq[h['subchain']]):
                continue
            m = RP.index_map(RP.one_letter(h['residues']), RP.one_letter(a['residues']))
            P, Q = [], []
            for i in shell:
                if i in m:
                    hb = RP.backbone(h['residues'][i], h['R'], h['t']); ab = RP.backbone(a['residues'][m[i]], a['R'], a['t'])
                    if all(k in hb and k in ab for k in ('N', 'CA', 'C', 'O')):
                        P += [hb[k] for k in ('N', 'CA', 'C', 'O')]; Q += [ab[k] for k in ('N', 'CA', 'C', 'O')]
            if len(P) >= 8:
                R, T, r = A.fit(np.array(P), np.array(Q))
                res.append((round(r, 2), a['chain'], len(P) // 4, round(float(np.linalg.norm((L @ R + T).mean(0) - L.mean(0))), 1)))
        print('  holo', h['chain'], 'shell', len(shell), 'best apo pairings (rmsd, apo chain, residues, ligand centroid move A):', sorted(res)[:3])
