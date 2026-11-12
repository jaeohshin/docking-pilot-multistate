"""Static coordinate illustration of two exploratory cases; projection never used for distances."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
O=Path(__file__).resolve().parent
CASES=[('28cfc7bf72ea90b1','MET','237','A','MET237: steric relief'),('88c5041cdb0042d0','ILE','621','A','ILE621: holo-like alternative')]
def atoms(p):
 out={}
 for l in p.read_text().splitlines():
  if l.startswith(('ATOM  ','HETATM')) and l[76:78].strip()!='H':out[(l[21],l[22:27].strip(),l[17:20],l[12:16].strip())]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 return out
def res(d,ch,num,name):return {k[-1]:v for k,v in d.items() if k[:3]==(ch,num,name)}
def fit(x,y):
 P=np.array([x[n] for n in ['N','CA','C']]);Q=np.array([y[n] for n in ['N','CA','C']]);pc=P.mean(0);qc=Q.mean(0);u,_,v=np.linalg.svd((P-pc).T@(Q-qc));R=u@np.diag([1,1,np.linalg.det(u@v)])@v
 return {k:(z-pc)@R+qc for k,z in x.items()}
fig,axes=plt.subplots(1,2,figsize=(11,5))
for ax,(cid,name,num,ch,title) in zip(axes,CASES):
 d=O/'structures'/cid;b=json.load(open(d/'build.json'));p=sorted(d.glob('apo_multi_*.pdb'))[0];x=res(atoms(d/'apo_single_0.pdb'),ch,num,name);y=res(atoms(p),ch,num,name);h=fit(res(atoms(d/'holo_single_0.pdb'),ch,num,name),x)
 m=Chem.RemoveHs(Chem.MolFromMolFile(str(d/'apo_multi_oracle.sdf'),removeHs=False));L=m.GetConformer().GetPositions();sc=[k for k in x if k not in ['N','CA','C','O']];allx=np.concatenate([L,np.array(list(x.values())),np.array(list(y.values()))]);center=allx.mean(0);_,_,v=np.linalg.svd(allx-center);rot=v[:2].T;pr=lambda z:(np.array(z)-center)@rot
 ll=pr(L)
 for bond in m.GetBonds():
  i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx();ax.plot(ll[[i,j],0],ll[[i,j],1],c='#555555',lw=1,zorder=1)
 ax.scatter(ll[:,0],ll[:,1],s=22,c='#555555',label='Oracle ligand',zorder=2)
 for zz,color,label,ls in [(x,'#377eb8','Apo major','-'),(y,'#e41a1c','Selected alternative','-'),(h,'#4daf4a','Holo major (local fit)','--')]:
  keys=list(zz);X=np.array([zz[k] for k in keys]);q=pr(X)
  for i in range(len(keys)):
   for j in range(i):
    if .7<np.linalg.norm(X[i]-X[j])<1.95:ax.plot(q[[i,j],0],q[[i,j],1],c=color,lw=1.5,ls=ls,alpha=.75)
  ax.scatter(q[:,0],q[:,1],s=15,c=color,label=label)
 for k in sc:
  if k in ['SG','SD','CD1','CG2']:
   z=pr(x[k]);ax.annotate(k,z,xytext=(3,4),textcoords='offset points',fontsize=8,color='#377eb8')
 old=np.linalg.norm(np.array([x[k] for k in sc])[:,None]-L[None,:],axis=-1).min();new=np.linalg.norm(np.array([y[k] for k in sc])[:,None]-L[None,:],axis=-1).min()
 ax.set_title(title+'\n'+cid+'\n'+f'3D side-chain–ligand minimum: {old:.2f} → {new:.2f} Å',fontsize=10);ax.set_aspect('equal');ax.set_xlabel('Coordinate projection 1 (Å)');ax.set_ylabel('Coordinate projection 2 (Å)');ax.legend(fontsize=7,loc='best')
fig.suptitle('Exploratory structure examples: saved coordinates, no minimization',fontsize=12);fig.tight_layout();fig.savefig(O/'structure_examples.png',dpi=180)
