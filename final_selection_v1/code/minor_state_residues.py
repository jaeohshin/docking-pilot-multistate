"""apo-minor chi1 == holo-major chi1 residues of the 79 candidates, with real residue names and author numbers.

Source definition (scripts/180_apo_holo.py): per holo qFit chain, `i` is the 0-based index in the chain's polymer
(gemmi get_polymer order, equal to the replicate-cache index); wells occupied at normalised occupancy >= 0.1;
selected = holo major well != apo major well and holo major in the apo occupied wells. Rows count for a case
when the case's ligand is in the row's `ligands` (holo-defined 4.5 A contact).
Here every case is recomputed from the replicate cache (chi1 wells) and the source rows, the 0-based index is
mapped to the qFit/deposited author residue through the case's `residue_map`, and the residue names are checked.
Reads only; writes minor_state_residues.tsv and minor_state_check.json in this folder.
"""
import csv, glob, gzip, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
ROOT = PILOT.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
OCC = 0.1
C = json.loads((PILOT/'final_review_v1/candidates79.json').read_text())
need = {c['holo'] for c in C} | {c['apo'] for c in C}
E = {}
for f in sorted(glob.glob(str(ROOT/'reference/frozen/data/qfit_at_scale/replicate_cache/*.jsonl'))):
    for line in open(f):
        pid = line[12:16]
        if pid in need:
            e = json.loads(line)
            E[e['pdb_id']] = e
rows = list(csv.DictReader(gzip.open(ROOT/'results/apo_holo/residues.csv.gz', 'rt')))
# the source aligner and index map, imported from the frozen source script
import importlib.util
spec = importlib.util.spec_from_file_location('s180', ROOT/'scripts/180_apo_holo.py'); s180 = importlib.util.module_from_spec(spec)
sys.modules['s180'] = s180; spec.loader.exec_module(s180)


import zipfile
import gemmi
Z = zipfile.ZipFile('/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip')
_models = {}


def qfit_key(pid, chain, idx):
    """(author number, icode, resname) of the idx-th (0-based) polymer residue of a qFit chain (source indexing)."""
    if pid not in _models:
        st = gemmi.read_pdb_string(Z.read(f'qfit_PDBs/{pid}_qFit_010.pdb').decode('utf-8', 'replace'))
        st.setup_entities(); st.remove_hydrogens(); _models[pid] = st
    polys = [ch.get_polymer() for ch in _models[pid][0] if ch.name == chain and len(ch.get_polymer())]
    assert len(polys) == 1
    r = polys[0][int(idx)]
    return (r.seqid.num, r.seqid.icode.strip(), r.name)


def deposited_chain(pid, qchain, key):
    """Deposited auth chain holding this author residue (same chain name expected; checked on the mmCIF)."""
    st = gemmi.read_structure(str(PILOT/'final_review_v1/cifs'/f'{pid}.cif.gz'))
    hits = [ch.name for ch in st[0] for r in ch if r.seqid.num == key[0] and r.seqid.icode.strip() == key[1] and r.name == key[2]
            and r.entity_type == gemmi.EntityType.Polymer]
    assert qchain in hits, (pid, qchain, key, hits)
    return qchain


def occupied(w):
    tot = sum(w.values()) or 1.0
    return {k: v / tot for k, v in w.items() if v / tot >= OCC}


out, checks = [], []
for c in C:
    eh, ea = E[c['holo']], E[c['apo']]
    hseq = eh['chains'][c['holo_chain']]['seq']
    # the source's apo chain: best alignment over all apo chains (180_apo_holo.work); may differ from c['apo_chain']
    src_ac, imap = max([(k, s180.index_map(hseq, v['seq'])) for k, v in ea['chains'].items()], key=lambda t: len(t[1]))
    hs = {s['i']: s for s in eh['chains'][c['holo_chain']]['states']}
    as_ = {s['i']: s for s in ea['chains'][src_ac]['states']}
    src = [r for r in rows if r['holo'] == c['holo'] and r['apo'] == c['apo'] and r['chain'] == c['holo_chain']
           and c['ligand_id'] in r['ligands'].split(';')]
    rmap = {m['holo_index']: m for m in c['residue_map']}
    recomputed = []
    for r in src:
        i = int(r['i']); j = imap.get(i)
        # same skips as 180_apo_holo.compare (every source row passed them, so a skip here is an inconsistency)
        assert j is not None and i in hs and j in as_ and as_[j]['resname'] == hs[i]['resname'], (c['case_id'], i, j)
        H, A = occupied(hs[i]['wells']), occupied(as_[j]['wells'])
        h = max(H, key=H.get); am = max(A, key=A.get)
        cls = 'same' if h == am else 'selected' if h in A else 'new'
        assert cls == r['cls'], (c['case_id'], i, cls, r['cls'])
        if cls != 'selected':
            continue
        m = rmap[i]
        hk = m['holo_key']
        ak = qfit_key(c['apo'], src_ac, j)          # author key of the source's apo residue, from the qFit model
        if src_ac == c['apo_chain']:
            assert m['apo_index'] == j and tuple(m['apo_key']) == ak, (c['case_id'], i, 'apo index/key differ')
        assert qfit_key(c['holo'], c['holo_chain'], i) == tuple(hk), (c['case_id'], i, 'holo key differs')
        assert hs[i]['resname'] == as_[j]['resname'] == hk[2] == ak[2] == r['resname']
        apo_dep = c['apo_deposited_chain'] if src_ac == c['apo_chain'] else deposited_chain(c['apo'], src_ac, ak)
        recomputed.append(i)
        out.append(dict(case_id=c['case_id'], homology_component=None, apo=c['apo'], holo=c['holo'], ligand=c['ligand_id'],
                        source_index_0based=i, apo_index_0based=j, resname=r['resname'],
                        holo_residue=f"{c['holo_deposited_chain']}:{r['resname']}{m['holo_key'][0]}{m['holo_key'][1]}",
                        apo_residue=f"{apo_dep}:{r['resname']}{ak[0]}{ak[1]}",
                        holo_qfit_chain=c['holo_chain'], apo_qfit_chain=src_ac, candidate_apo_chain=c['apo_chain'],
                        holo_chi1_wells=json.dumps({k: round(v, 3) for k, v in sorted(H.items())}),
                        apo_chi1_wells=json.dumps({k: round(v, 3) for k, v in sorted(A.items())}),
                        holo_major=h, apo_major=am, apo_minor_matching_holo_major=round(A[h], 3)))
    checks.append(dict(case_id=c['case_id'], source_apo_chain=src_ac, candidate_apo_chain=c['apo_chain'], source_rows=len(src), recomputed_selected=sorted(recomputed),
                       candidate_selected=sorted(int(x) for x in c['selected_residues']),
                       match=sorted(recomputed) == sorted(int(x) for x in c['selected_residues'])))
comp = {r['case_id']: r['homology_component'] for r in json.loads((PILOT/'filter_review_v1/review_records.json').read_text())}
for r in out:
    r['homology_component'] = comp[r['case_id']]
with open(HERE/'minor_state_residues.tsv', 'w', newline='') as h:
    w = csv.DictWriter(h, list(out[0]), delimiter='\t'); w.writeheader(); w.writerows(out)
(HERE/'minor_state_check.json').write_text(json.dumps(dict(cases=len(checks), all_match=all(x['match'] for x in checks),
                                                          positive_cases=sorted({r['case_id'] for r in out}), checks=checks), indent=1) + '\n')
print('cases', len(checks), 'all match', all(x['match'] for x in checks), 'positive', sorted({r['case_id'] for r in out}))
for r in out:
    print(r['case_id'], r['holo'], r['ligand'], 'i=', r['source_index_0based'], r['holo_residue'], '| apo', r['apo_residue'], r['holo_chi1_wells'], r['apo_chi1_wells'])
