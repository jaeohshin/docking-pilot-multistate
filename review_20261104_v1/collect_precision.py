"""Read existing SDF properties and boundary RMSDs; no engine or model calls."""
import sys,csv,gzip,json,re,importlib.util
from pathlib import Path
R=Path('/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028')
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule,RDKitMolCreate
RDLogger.DisableLog('rdApp.*')
w=csv.writer(sys.stdout,delimiter='\t');w.writerow(['engine','case_id','arm','call','pose','cnnscore','cnnaffinity','minaffinity','rmsd_exact_boundary'])
# Full precision properties can be read without recomputing coordinates/RMSD.
for f in sorted((R/'followup_v1/gnina_dock/docked').glob('*/*/docking_log.json')):
 cid=f.parent.parent.name;arm=f.parent.name
 for call in json.load(open(f))['calls']:
  text=(f.parent/call['out']).read_text()
  for i,rec in enumerate(text.split('$$$$')):
   if not rec.strip():continue
   props={k:float(v) for k,v in re.findall(r'>\s*<([^>]+)>\s*\n\s*([-+\deE.]+)',rec)}
   w.writerow(['GNINA8',cid,arm,call['run'],i,-props['CNNscore'],-props['CNNaffinity'],props['minimizedAffinity'],''])
for run,eng in [('followup_v1/budget50','Vina50'),('followup_v1/gnina_dock','GNINA8')]:
 with gzip.open(R/run/'pose_table.tsv.gz','rt') as h:bound=[r for r in csv.DictReader(h,delimiter='\t') if float(r['rmsd'])==2]
 for r in bound:
  cid,arm=r['case_id'],r['arm'];calli=int(r['call']);i=int(r['pose']);d=R/run/'docked'/cid/arm;call=json.load(open(d/'docking_log.json'))['calls'][calli];ref=Chem.RemoveAllHs(Chem.MolFromMolFile(str(R/'expand_v1/inputs'/cid/('ligand_ref_'+arm.split('_')[0]+'.sdf')),removeHs=False))
  if eng=='Vina50':
   pm=PDBQTMolecule.from_file(str(d/call['out']),skip_typing=True);mol=Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm,only_cluster_leads=False)[0]);rm=rdMolAlign.CalcRMS(mol,ref,prbId=i)
  else:
   mol=Chem.RemoveAllHs(list(Chem.SDMolSupplier(str(d/call['out']),removeHs=False))[i]);rm=rdMolAlign.CalcRMS(mol,ref)
  w.writerow([eng,cid,arm,calli,i,'','','',rm])
