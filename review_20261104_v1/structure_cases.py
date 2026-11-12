#!/usr/bin/env python3
"""Exploratory geometry of all 10 preselected oracle switches; existing coordinates only.
No minimization, docking, or model inference. Requires existing dock environment with Meeko/RDKit.
"""
import argparse,json,re,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule,RDKitMolCreate
RDLogger.DisableLog('rdApp.*')
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();S=a.source;O=a.out;D=O/'structures';D.mkdir(exist_ok=True,parents=True)
b=pd.read_csv(S/'followup_v1/budget50/docking_results.tsv',sep='\t');b=b[b.stratum=='accept'];cases=b[b.apo_single_success_any!=b.apo_multi_success_any]
t=pd.read_csv(S/'followup_v1/budget50/pose_table.tsv.gz',sep='\t');cand={x['case_id']:x for x in json.load(open(S/'expand_v1/candidates.json'))};state={x['case_id']:x for x in json.load(open(S/'expand_v1/qfit_state_audit.json'))['cases']}
checks=[];result=[];resrows=[];recrows=[];callrows=[];dedup=[]
def atoms(p,pdbqt=False):
 d={}
 for l in p.read_text().splitlines():
  if not l.startswith(('ATOM  ','HETATM')):continue
  elem=l.split()[-1] if pdbqt else l[76:78].strip()
  if elem in ['H','HD','HS'] or (not elem and l[12:16].strip().startswith('H')):continue
  k=(l[21],l[22:27].strip(),l[17:20],l[12:16].strip());d[k]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 return d
def contacts(at,L):
 x=np.array(list(at.values()));d=np.linalg.norm(x[:,None]-L[None,:],axis=-1)
 return dict(min_heavy_distance=float(d.min()),pairs_lt2p2=int((d<2.2).sum()),ligand_atoms_lt2p2=int((d<2.2).any(axis=0).sum()))
def aligned_rms(x,y):
 bb=['N','CA','C'];names=sorted((set(x)&set(y))-{'N','CA','C','O'})
 if not names or not all(n in x and n in y for n in bb):return None
 P=np.array([x[n] for n in bb]);Q=np.array([y[n] for n in bb]);pc=P.mean(0);qc=Q.mean(0);u,_,v=np.linalg.svd((P-pc).T@(Q-qc));R=u@np.diag([1,1,np.linalg.det(u@v)])@v
 return float(np.sqrt(np.mean(np.sum((((np.array([x[n] for n in names])-pc)@R+qc)-np.array([y[n] for n in names]))**2,axis=1))))
def residue(d,ch,num,name):return {k[3]:v for k,v in d.items() if k[:3]==(ch,str(num),name)}
for row in cases.itertuples():
 cid=row.case_id;out=D/cid;out.mkdir(exist_ok=True);inp=S/'expand_v1/inputs'/cid;rec=S/'expand_v1/receptors'/cid;build=json.load(open(rec/'build.json'));c=cand[cid]
 ref=Chem.RemoveHs(Chem.MolFromMolFile(str(inp/'ligand_ref_apo.sdf'),removeHs=False));L=ref.GetConformer().GetPositions();maj=atoms(inp/'apo_single_0.pdbqt',True)
 armdata={};maxerr=0
 for arm in ['apo_single','apo_multi']:
  base=S/'followup_v1/budget50/docked'/cid/arm;log=json.load(open(base/'docking_log.json'));records=[]
  for call in log['calls']:
   f=base/call['out'];pm=PDBQTMolecule.from_file(str(f),skip_typing=True);mol=Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm,only_cluster_leads=False)[0]);scores=[float(p.score) for p in pm]
   saved=t[(t.case_id==cid)&(t.arm==arm)&(t.call==call['run'])].sort_values('pose');assert len(saved)==len(scores)==mol.GetNumConformers()
   rms=[]
   for i,score in enumerate(scores):
    rm=rdMolAlign.CalcRMS(mol,ref,prbId=i);err=abs(round(rm,3)-saved.iloc[i].rmsd);maxerr=max(maxerr,err);assert err<.00101,(cid,arm,call['run'],i,err)
    records.append(dict(call=call['run'],pose=i,score=score,rmsd=rm,mol=mol,receptor=call['receptor'],seed=call['seed']));rms.append(rm)
   callrows.append(dict(case_id=cid,arm=arm,call=call['run'],seed=call['seed'],receptor=call['receptor'],oracle=min(rms),success=int(min(rms)<=2),good_poses=int(np.sum(np.array(rms)<=2))))
  order=sorted(records,key=lambda z:z['score']);kept=[]
  for r in order:
   if all(rdMolAlign.CalcRMS(r['mol'],k['mol'],prbId=r['pose'],refId=k['pose'])>1 for k in kept):kept.append(r)
   if len(kept)==5:break
  calc=min(r['rmsd'] for r in kept);expected=getattr(row,arm+'_rmsd_top5');assert abs(round(calc,3)-expected)<.00101
  dedup.append(dict(case_id=cid,arm=arm,recomputed_top5=round(calc,3),frozen_top5=expected,raw_top5=round(min(r['rmsd'] for r in order[:5]),3),first_good_rank=next((i+1 for i,r in enumerate(order) if r['rmsd']<=2),None),n_kept=len(kept),kept_call_pose=';'.join(f"{r['call']}:{r['pose']}" for r in kept)))
  armdata[arm]=records
  best=min(records,key=lambda r:r['rmsd']);w=Chem.SDWriter(str(out/(arm+'_oracle.sdf')));w.write(best['mol'],confId=best['pose']);w.close()
 best=min(armdata['apo_multi'],key=lambda r:r['rmsd']);ri=int(re.search(r'multi_(\d+)',best['receptor']).group(1));combo=build['apo']['combos'][ri];changes={res:well for res,well in combo.items() if well!=build['apo']['multi_state_residues'][res][0][0]}
 winning=atoms(inp/best['receptor'],True);bp=best['mol'].GetConformer(best['pose']).GetPositions()
 sm=contacts(maj,bp);mm=contacts(winning,bp);sr=contacts(maj,L);mr=contacts(winning,L)
 # Save only actual PDB geometries required for viewing, preserving original source files.
 for fn in ['apo_single_0.pdb',f'apo_multi_{ri}.pdb','holo_single_0.pdb','ligand_ref_apo.pdb']:(out/fn).write_bytes((rec/fn).read_bytes())
 (out/'build.json').write_bytes((rec/'build.json').read_bytes())
 # Matched-seed all-major checks (coordinates and scores in existing outputs).
 single={r['call']:[] for r in armdata['apo_single']}
 for r in armdata['apo_single']:single[r['call']].append(r)
 eq=[]
 for r in armdata['apo_multi']:
  if r['receptor']=='apo_multi_0.pdbqt':
   ss=single[r['call']][r['pose']];eq.append(np.array_equal(ss['mol'].GetConformer(ss['pose']).GetPositions(),r['mol'].GetConformer(r['pose']).GetPositions()) and ss['score']==r['score'])
 major_equal=set(maj)==set(atoms(inp/'apo_multi_0.pdbqt',True)) and all(np.array_equal(v,atoms(inp/'apo_multi_0.pdbqt',True)[k]) for k,v in maj.items())
 checks.append(dict(case_id=cid,pose_table_max_rounding_difference=maxerr,all_major_same_atoms=major_equal,all_major_same_seed_outputs_equal=all(eq),compared_poses=len(eq),backbone_equal=all(k in winning and np.array_equal(v,winning[k]) for k,v in maj.items() if k[3] in ['N','CA','C','O'])))
 cr=pd.DataFrame([q for q in callrows if q['case_id']==cid]);multi=cr[cr.arm=='apo_multi'];singlecalls=cr[cr.arm=='apo_single'].set_index('call');successrecs=[]
 for rn,rr in multi.groupby('receptor'):
  rnum=int(re.search(r'multi_(\d+)',rn).group(1));comb=build['apo']['combos'][rnum];chg={r:w for r,w in comb.items() if w!=build['apo']['multi_state_residues'][r][0][0]}
  recrows.append(dict(case_id=cid,receptor=rn,calls=len(rr),successful_calls=int(rr.success.sum()),good_poses=int(rr.good_poses.sum()),best_rmsd=rr.oracle.min(),paired_single_successful_calls=int(singlecalls.loc[rr.call,'success'].sum()),changed_states=json.dumps(chg,sort_keys=True)))
  if rr.success.sum():successrecs.append(rn)
 ap=atoms(rec/'apo_single_0.pdb');wp=atoms(rec/f'apo_multi_{ri}.pdb');hp=atoms(rec/'holo_single_0.pdb')
 for rn,newwell in changes.items():
  label,rest=rn.split(':');mat=re.fullmatch(r'([A-Z]{3})(-?\d+)([A-Za-z]?)',rest);name,num,icode=mat.groups();numkey=num+icode;ch=build['apo']['chain_names_in_pdb'][label];aa=residue(ap,ch,numkey,name);ww=residue(wp,ch,numkey,name);sts=build['apo']['multi_state_residues'][rn];new=[x for x in sts if x[0]==newwell][0]
  # Explicit source residue_map only for the primary receptor chain/copy; other copies not inferred.
  mappings=[m for m in c.get('residue_map',[]) if m['apo_key']==[int(num),icode,name]]
  hc=build['holo']['chain_names_in_pdb'].get(c['holo_chain']+'0');hh={};hk=None
  if label==c['apo_chain']+'0' and hc is not None and len(mappings)==1:
   hk=mappings[0]['holo_key'];hh=residue(hp,hc,str(hk[0])+hk[1],hk[2])
  oldr=aligned_rms(aa,hh) if hh else None;newr=aligned_rms(ww,hh) if hh else None
  def minres(z,ll):return float(np.linalg.norm(np.array([v for n,v in z.items() if n not in ['N','CA','C','O']])[:,None]-ll[None,:],axis=-1).min())
  qh=[r for r in state[cid]['holo_residues'] if r['residue']==f"{c['holo_chain']}:{hk[2]}{hk[0]}{hk[1]}"] if hk else []
  wells=qh[0].get('chi12_wells') if len(qh)==1 else None
  holowell=max(wells,key=wells.get) if wells else None
  resrows.append(dict(case_id=cid,receptor=best['receptor'],residue_source=rn,pdb_chain=ch,author_residue=numkey,resname=name,major_well=sts[0][0],major_occupancy=sts[0][1],major_altloc=sts[0][2],chosen_well=newwell,chosen_occupancy=new[1],chosen_altloc=new[2],holo_residue_key=json.dumps(hk),holo_major_well=holowell,chosen_matches_holo_well=(newwell==holowell) if holowell else None,major_matches_holo_well=(sts[0][0]==holowell) if holowell else None,major_local_sc_rmsd_to_holo=oldr,chosen_local_sc_rmsd_to_holo=newr,major_min_distance_to_winning_ligand=minres(aa,bp),chosen_min_distance_to_winning_ligand=minres(ww,bp),major_min_distance_to_reference=minres(aa,L),chosen_min_distance_to_reference=minres(ww,L)))
 result.append(dict(case_id=cid,apo=row.apo,holo=row.holo,ligand=row.ligand,component=row.component,transition='gain' if row.apo_multi_success_any else 'loss',single_oracle=row.apo_single_rmsd_oracle,multi_oracle=row.apo_multi_rmsd_oracle,best_multi_call=best['call'],best_multi_pose=best['pose'],best_multi_seed=best['seed'],best_multi_receptor=best['receptor'],best_multi_score=best['score'],changed_residues=json.dumps(changes,sort_keys=True),successful_multi_receptors=';'.join(successrecs),single_successful_calls=int(cr[cr.arm=='apo_single'].success.sum()),multi_successful_calls=int(multi.success.sum()),winning_pose_clashes_major=sm['pairs_lt2p2'],winning_pose_clashes_selected=mm['pairs_lt2p2'],winning_pose_min_distance_major=sm['min_heavy_distance'],winning_pose_min_distance_selected=mm['min_heavy_distance'],reference_clashes_major=sr['pairs_lt2p2'],reference_clashes_selected=mr['pairs_lt2p2'],note='Static geometry and fixed-seed evidence; no single-residue causality or independent-seed replication.'))
 print(cid,len(changes),'changes','good calls',result[-1]['multi_successful_calls'],'clashes',sm['pairs_lt2p2'],mm['pairs_lt2p2'],flush=True)
for x,n in [(result,'structure_case_summary.tsv'),(resrows,'changed_residues.tsv'),(recrows,'receptor_call_results.tsv'),(callrows,'matched_seed_calls.tsv'),(dedup,'distinct_top5_geometry_check.tsv'),(checks,'structure_checks.tsv')]:pd.DataFrame(x).to_csv(O/n,sep='\t',index=False)
print('done',len(cases),'cases')
