import json,csv,glob,gzip,zipfile,hashlib,time,urllib.request,sys
from pathlib import Path
from collections import Counter
import gemmi
import numpy as np
from Bio import Align

ROOT=Path('/gpfs/jaeohshin/work/sidechain-multistate')
OUT=Path(__file__).resolve().parents[1]
CACHE=OUT.parent/'run_v1/download_cache'; CACHE.mkdir(exist_ok=True)
# clean_cache.json은 상동성 정제 집합(333 구조)만 담는다. 확대 집합은 전체 replicate cache에서 읽는다.
_need={x for c in json.loads((OUT.parent/'run_v3/homology_candidates.json').read_text()) for x in (c['apo'],c['holo'])}
E={}
for _f in sorted(glob.glob('/gpfs/jaeohshin/work/sidechain-multistate/reference/frozen/data/qfit_at_scale/replicate_cache/*.jsonl')):
    for _line in open(_f):
        if _line[12:16] in _need:
            _e=json.loads(_line); E[_e['pdb_id']]=_e
assert len(E)==len(_need), f'replicate cache missing {len(_need)-len(E)} structures'
CS=json.loads((OUT.parent/'run_v3/homology_candidates.json').read_text())  # C: packer용 상동성 차단 해제 (이 실험은 packer를 쓰지 않음)
ALIGN=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
Z=zipfile.ZipFile('/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip')
BB=('N','CA','C','O')
METALS=set('LI NA K RB CS MG CA SR BA MN FE CO NI CU ZN CD HG AL GA PB YB GD'.split())
MODEL={}; CCD={}; DEP={}

def xyz(a): return np.array([a.pos.x,a.pos.y,a.pos.z])
def heavy(res): return [a for a in res if not a.is_hydrogen()]
def key(res): return (res.seqid.num,res.seqid.icode.strip(),res.name)
def dist(atoms,L):
    if not atoms: return float('inf')
    X=np.array([xyz(a) for a in atoms]); return float(np.sqrt(((X[:,None]-L[None])**2).sum(-1).min()))
def select_atoms(res):
    ats=heavy(res); alts=sorted({a.altloc for a in ats if a.altloc not in ('\x00',' ','')})
    if alts:
        alts.sort(key=lambda alt:(-np.mean([a.occ for a in ats if a.altloc==alt]),alt))
        alt=alts[0]
    else: alt='\x00'
    selected=[a for a in ats if a.altloc in (alt,'\x00',' ','')]
    names=[a.name for a in selected]
    if len(set(names))!=len(names): raise ValueError('duplicate atom name in selected conformer')
    return {a.name:a for a in selected},alt
def model(pid):
    if pid not in MODEL:
        content=Z.read('qfit_PDBs/'+pid+'_qFit_010.pdb').decode()
        st=gemmi.read_pdb_string(content); st.setup_entities(); st.remove_hydrogens(); st.setup_cell_images()
        MODEL[pid]=st
    return MODEL[pid]
def polymer(st,ch):
    matches=[list(c.get_polymer()) for c in st[0] if c.name==ch and len(c.get_polymer())]
    if len(matches)!=1: raise ValueError('ambiguous polymer chain')
    p=matches[0]
    if len({key(x)[:2] for x in p})!=len(p): raise ValueError('microheterogeneous polymer residue')
    return p
def sequence(poly): return ''.join(gemmi.find_tabulated_residue(x.name).one_letter_code.upper() for x in poly)
def fetch(url,path):
    if not path.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url,timeout=30) as f: data=f.read()
                if not data: raise ValueError('empty download')
                path.write_bytes(data); break
            except Exception:
                if attempt==2: raise
                time.sleep(1+attempt)
    return path
def ccd(cid):
    if cid not in CCD:
        p=fetch('https://files.rcsb.org/ligands/download/'+cid+'.cif',CACHE/(cid+'.cif'))
        b=gemmi.cif.read(str(p)).sole_block()
        atomrows=[list(x) for x in b.find('_chem_comp_atom.',['atom_id','type_symbol'])]
        descriptors=[list(x) for x in b.find('_pdbx_chem_comp_descriptor.',['type','program','descriptor'])]
        def val(k):
            v=b.find_value(k)
            return gemmi.cif.as_string(v) if v else ''
        CCD[cid]={'id':cid,'type':val('_chem_comp.type'),'name':val('_chem_comp.name'),
                  'parent':val('_chem_comp.mon_nstd_parent_comp_id'),'atoms':atomrows,'descriptors':descriptors,
                  'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return CCD[cid]
def deposition(pid):
    if pid not in DEP:
        existing=ROOT/'work/binding_site/cif'/(pid+'.cif.gz')
        p=existing if existing.exists() else fetch('https://files.rcsb.org/download/'+pid+'.cif.gz',CACHE/(pid+'.cif.gz'))
        st=gemmi.read_structure(str(p)); st.setup_entities(); st.remove_hydrogens(); st.setup_cell_images()
        DEP[pid]=(st,p)
    return DEP[pid]
def fit(P,Q):
    pc=P.mean(0); qc=Q.mean(0)
    U,S,V=np.linalg.svd((P-pc).T@(Q-qc))
    R=U@np.diag([1,1,np.linalg.det(U@V)])@V
    T=qc-pc@R
    return R,T,float(np.sqrt(((P@R+T-Q)**2).sum(1).mean()))

# Sensitivity audit: retain the original flags, but do not skip downstream checks for these design restrictions.
SOFT={'heavy_atom_count_outside_10_50','ligand_occupancy_below_0_8','ccd_not_nonpolymer','pocket_backbone_rmsd_above_0_75',
      'apo_multichain_pocket','holo_multichain_pocket','apo_crystal_mate_pocket','holo_crystal_mate_pocket',
      'apo_nearby_metal','holo_nearby_metal','apo_nearby_other_heterocomponent','holo_nearby_other_heterocomponent'}
def hard_failure(reasons): return any(x not in SOFT for x in reasons)

def audit(c):
    d=dict(c); reasons=[]; d['geometry_reasons']=reasons
    hs=model(c['holo']); aps=model(c['apo'])
    hp=polymer(hs,c['holo_chain']); ap=polymer(aps,c['apo_chain'])
    hseq=sequence(hp); aseq=sequence(ap)
    if hseq!=E[c['holo']]['chains'][c['holo_chain']]['seq'] or aseq!=E[c['apo']]['chains'][c['apo_chain']]['seq']:
        reasons.append('qfit_cache_sequence_mismatch'); return d
    comp,lch,lnum=c['ligand_id'].split(':')
    ligs=[res for ch in hs[0] if ch.name==lch for res in ch if res.name==comp and res.seqid.num==int(lnum)]
    if len(ligs)!=1:
        reasons.append('ambiguous_ligand_instance'); return d
    lig=ligs[0]; la,lalt=select_atoms(lig); L=np.array([xyz(a) for a in la.values()])
    d.update(ligand_icode=lig.seqid.icode.strip(),ligand_altloc=lalt,ligand_heavy_atoms=len(la),
             ligand_min_occupancy=min(a.occ for a in la.values()))
    if not 10<=len(la)<=50: reasons.append('heavy_atom_count_outside_10_50')
    if d['ligand_min_occupancy']<.8-1e-5: reasons.append('ligand_occupancy_below_0_8')
    if hard_failure(reasons): return d
    chem=ccd(comp); d['ccd']=chem
    if chem['type'].upper()!='NON-POLYMER': reasons.append('ccd_not_nonpolymer')
    expected={gemmi.cif.as_string(n) for n,e in chem['atoms'] if e not in ('H','D')}
    if set(la)!=expected: reasons.append('ligand_heavy_atom_names_incomplete_or_extra')
    if hard_failure(reasons): return d
    aln=ALIGN.align(hseq,aseq)[0]; mapping={}
    for (h0,h1),(a0,a1) in zip(*aln.aligned):
        for n in range(h1-h0):
            if hseq[h0+n]==aseq[a0+n]: mapping[int(h0+n)]=int(a0+n)
    shell=[i for i,x in enumerate(hp) if dist(heavy(x),L)<=6.5]
    contact=[i for i in shell if dist(heavy(hp[i]),L)<=4.5]
    P=[]; Q=[]; mapped=[]; bbs_h={}; bbs_a={}
    for i in shell:
        if i not in mapping: continue
        j=mapping[i]
        ha,_=select_atoms(hp[i]); aa,_=select_atoms(ap[j])
        if not all(n in ha and n in aa for n in BB): continue
        mapped.append(i); bbs_h[i]=ha; bbs_a[j]=aa
        for n in BB: P.append(xyz(ha[n])); Q.append(xyz(aa[n]))
    d.update(pocket_residues=len(shell),contact_residues=len(contact),mapped_pocket_residues=len(mapped),
             pocket_coverage=len(mapped)/len(shell) if shell else 0)
    if not contact or len(mapped)<3: reasons.append('insufficient_pocket_atoms')
    if set(contact)-set(mapped): reasons.append('contact_backbone_mapping_or_atoms_missing')
    if d['pocket_coverage']<.95: reasons.append('pocket_coverage_below_95pct')
    if hard_failure(reasons): return d
    P=np.array(P); Q=np.array(Q); R,T,rmsd=fit(P,Q)
    d.update(pocket_backbone_rmsd=rmsd,pocket_ca_rmsd=float(np.sqrt(((P[1::4]@R+T-Q[1::4])**2).sum(1).mean())),
             rotation_holo_to_apo=R.tolist(),translation_holo_to_apo=T.tolist())
    if rmsd>.75: reasons.append('pocket_backbone_rmsd_above_0_75')
    # Check peptide junctions touching the shell, without synthesizing missing atoms.
    for tag,pol,idxs in [('holo',hp,set(shell)),('apo',ap,{mapping[i] for i in mapped})]:
        for i in sorted({j for k in idxs for j in (k-1,k) if 0<=j<len(pol)-1}):
            x,_=select_atoms(pol[i]); y,_=select_atoms(pol[i+1])
            if 'C' not in x or 'N' not in y:
                reasons.append(tag+'_pocket_junction_missing_atoms'); break
            dcn=float(np.linalg.norm(xyz(x['C'])-xyz(y['N'])))
            if not .9<=dcn<=2.0:
                reasons.append(tag+'_pocket_junction_gap_or_distortion'); break
    d['residue_map']=[{'holo_index':i,'apo_index':mapping[i],'holo_key':key(hp[i]),'apo_key':key(ap[mapping[i]])} for i in mapped]
    X=L@R+T
    d['box_size']=(np.maximum(X.max(0)-X.min(0)+10,20)).tolist()
    d['box_center']=((X.max(0)+X.min(0))/2).tolist()
    if max(d['box_size'])>30: reasons.append('box_edge_above_30')
    if hard_failure(reasons): return d
    # Exclude actual extra-chain contacts, crystal mates and metal/cofactor dependence.
    for tag,st,ref,chosen_chain in [('holo',hs,L,c['holo_chain']),('apo',aps,X,c['apo_chain'])]:
        for ch in st[0]:
            for res in ch:
                atoms=heavy(res)
                if tag=='holo' and ch.name==lch and key(res)==key(lig): continue
                mind=dist(atoms,ref)
                if mind<=4.5 and res.entity_type==gemmi.EntityType.Polymer and ch.name!=chosen_chain:
                    reasons.append(tag+'_multichain_pocket')
                if mind<=6 and any(a.element.name.upper() in METALS for a in atoms): reasons.append(tag+'_nearby_metal')
                if mind<=6 and res.het_flag=='H' and res.name not in ('HOH','DOD') and len(atoms)>=6:
                    reasons.append(tag+'_nearby_other_heterocomponent')
        ns=gemmi.NeighborSearch(st[0],st.cell,4.5).populate(include_h=False)
        mate=False
        for pos in ref:
            gp=gemmi.Position(*pos)
            for mark in ns.find_atoms(gp,'\x00',radius=4.5):
                cra=mark.to_cra(st[0])
                if cra.residue.entity_type!=gemmi.EntityType.Polymer: continue
                if mark.image_idx!=0 or np.linalg.norm(xyz(cra.atom)-pos)>4.5+1e-5:
                    mate=True; break
            if mate: break
        if mate: reasons.append(tag+'_crystal_mate_pocket')
    if hard_failure(reasons):
        d['geometry_reasons']=sorted(set(reasons)); return d
    # Deposited mmCIF identifies polymeric ligands and explicit covalent/metal links.
    dep,p=deposition(c['holo']); d['deposited_cif_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    dep_lch=c.get('ligand_deposited_chain',lch)
    dl=[res for ch in dep[0] if ch.name==dep_lch for res in ch if key(res)==key(lig)]
    if len(dl)!=1: reasons.append('deposited_ligand_instance_mismatch')
    elif dl[0].entity_type!=gemmi.EntityType.NonPolymer: reasons.append('deposited_ligand_not_nonpolymer')
    else:
        da,_=select_atoms(dl[0])
        if set(da)!=set(la): reasons.append('deposited_ligand_atoms_mismatch')
        else:
            delta=float(np.sqrt(np.mean([np.sum((xyz(da[n])-xyz(la[n]))**2) for n in la])))
            d['qfit_vs_deposited_ligand_rmsd']=delta
            if delta>.25: reasons.append('qfit_deposited_ligand_frame_mismatch')
    for conn in dep.connections:
        for addr in (conn.partner1,conn.partner2):
            if addr.chain_name==dep_lch and addr.res_id.name==comp and addr.res_id.seqid.num==int(lnum) and addr.res_id.seqid.icode.strip()==lig.seqid.icode.strip():
                if conn.type in (gemmi.ConnectionType.Covale,gemmi.ConnectionType.MetalC): reasons.append('deposited_covalent_or_metal_link')
    # Conservative geometry catch for unannotated covalent contacts.
    for ch in dep[0]:
        for res in ch:
            if res.entity_type!=gemmi.EntityType.Polymer: continue
            if dist(heavy(res),L)<1.8: reasons.append('possible_unannotated_covalent_contact')
    d['geometry_reasons']=sorted(set(reasons))
    return d

# Independent rigid-transform check for the coordinate alignment implementation.
rng=np.random.default_rng(7); A=rng.normal(size=(20,3)); U,_,V=np.linalg.svd(rng.normal(size=(3,3))); RR=U@V
if np.linalg.det(RR)<0: RR[:,0]*=-1
B=A@RR+np.array([3,-8,2]); rot,tr,err=fit(A,B)
assert err<1e-10 and np.linalg.det(rot)>.999999
results=[]
for i,c in enumerate(CS):
    try: d=audit(c)
    except Exception as e: d=dict(c,geometry_reasons=['audit_exception'],audit_exception=repr(e))
    results.append(d)
    if i%25==0: print('geometry',i,'/',len(CS),'pass',sum(not x['geometry_reasons'] for x in results),'exceptions',sum('audit_exception' in x for x in results),flush=True)
    (OUT/'sensitivity_candidates.json').write_text(json.dumps(results,indent=1)+'\n')
passed=[x for x in results if not x['geometry_reasons']]
summary={'status':'sensitivity_flags_complete_rdkit_pending','cases':len(results),'passed_cases':len(passed),
         'passed_groups':len({x['group'] for x in passed}),'reasons':dict(Counter(y for x in results for y in x['geometry_reasons'])),
         'exceptions':[{k:x.get(k) for k in ('case_id','apo','holo','audit_exception')} for x in results if 'audit_exception' in x]}
(OUT/'sensitivity_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)

