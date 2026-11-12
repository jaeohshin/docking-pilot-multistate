"""Read-only collector run on Syntax: streams only existing selected files as a tar.gz."""
import csv,io,json,tarfile,sys,hashlib
from pathlib import Path
root=Path('/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028')
rows=list(csv.DictReader(open(root/'followup_v1/budget50/docking_results.tsv'),delimiter='\t'))
ids=[r['case_id'] for r in rows if r['stratum']=='accept' and r['apo_single_success_any']!=r['apo_multi_success_any']]
files=[]
for r in rows:
 c=r['case_id']
 files.extend((root/'expand_v1/receptors'/c).glob('build.json'))
 for run in ['expand_v1','followup_v1/budget50','followup_v1/gnina_dock']:
  files.extend((root/run/'docked'/c).glob('*/docking_log.json'))
for c in ids:
 for run in ['followup_v1/budget50']:
  for arm in ['apo_single','apo_multi']:
   files.extend((root/run/'docked'/c/arm).glob('*.pdbqt'))
 files.extend((root/'expand_v1/receptors'/c).glob('*.pdb'))
 for pattern in ['apo*.pdbqt','holo_single_0.pdbqt','ligand*.sdf','*.json']:
  files.extend((root/'expand_v1/inputs'/c).glob(pattern))
manifest=[]
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for p in sorted(set(files)):
  tar.add(p,arcname=str(p.relative_to(root)),recursive=False)
  manifest.append({'path':str(p.relative_to(root)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 raw=json.dumps({'selected_cases':ids,'files':manifest},indent=2).encode()
 ti=tarfile.TarInfo('COLLECTED_MANIFEST.json');ti.size=len(raw);tar.addfile(ti,io.BytesIO(raw))
print('files',len(manifest),'bytes',sum(x['bytes'] for x in manifest),file=sys.stderr)
