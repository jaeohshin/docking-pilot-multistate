import json,csv,hashlib,random
from pathlib import Path
from collections import Counter,defaultdict
import gemmi
from rdkit import Chem,rdBase
from rdkit.Chem import rdMolDescriptors

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'selection_v3'; OUT.mkdir(exist_ok=True)
records=json.loads((ROOT/'remote_results_v3/geometry_v3_candidates.json').read_text())
hom=json.loads((ROOT/'remote_results_v3/homology.json').read_text())
assert len(records)==hom['remaining']['ligand_chain_cases']
assert len({c['case_id'] for c in records})==len(records)

# Connectivity components include all homology-clean sequences, including geometry failures.
parent={}
def root(x):
    parent.setdefault(x,x)
    if parent[x]!=x: parent[x]=root(parent[x])
    return parent[x]
for h in hom['edges']:
    assert float(h[2])>=.3 and max(float(h[6]),float(h[7]))>=.5 and float(h[4])<=.001
    a,b=root(h[0]),root(h[1]); parent[max(a,b)]=min(a,b)
for c in records: c['homology_component']=root(c['query_id'])

def chemistry(c):
    reasons=[]
    if c['geometry_reasons']: return ['not_evaluated_geometry_failed']
    cd=c['ccd']; smiles=[]
    for typ,program,s in cd['descriptors']:
        typ,program,s=[gemmi.cif.as_string(x) for x in (typ,program,s)]
        if typ.startswith('SMILES'): smiles.append((typ!='SMILES_CANONICAL',program,s))
    molecule=None
    for _,prog,s in sorted(smiles):
        m=Chem.MolFromSmiles(s)
        if m is not None:
            molecule=Chem.RemoveHs(m); c['smiles']=s; c['smiles_program']=prog; break
    if molecule is None: return ['ccd_smiles_unreadable']
    c['rdkit_heavy_atoms']=molecule.GetNumHeavyAtoms()
    c['rotatable_bonds']=rdMolDescriptors.CalcNumRotatableBonds(molecule,strict=True)
    c['formal_charge']=Chem.GetFormalCharge(molecule)
    c['molecular_formula']=rdMolDescriptors.CalcMolFormula(molecule)
    c['ligand_name']=cd['name']
    if c['rdkit_heavy_atoms']!=c['ligand_heavy_atoms']: reasons.append('ccd_smiles_structure_heavy_count_mismatch')
    if not 10<=c['rdkit_heavy_atoms']<=50: reasons.append('rdkit_size_outside_10_50')
    if c['rotatable_bonds']>12: reasons.append('rotatable_bonds_above_12')
    if len(Chem.GetMolFrags(molecule))!=1: reasons.append('disconnected_ligand_graph')
    allowed={1,5,6,7,8,9,15,16,17,35,53}
    if not any(a.GetAtomicNum()==6 for a in molecule.GetAtoms()): reasons.append('no_carbon')
    if any(a.GetAtomicNum() not in allowed for a in molecule.GetAtoms()): reasons.append('metal_or_unsupported_element')
    # Supplement the source's explicit cofactor exclusions for common glutathione species.
    if cd['id'] in {'GSH','GSS','GSF','GTT'}: reasons.append('glutathione_cofactor_family')
    return reasons

for c in records:
    c['chemistry_reasons']=chemistry(c)
    c['eligible']=not c['geometry_reasons'] and not c['chemistry_reasons']
eligible=[c for c in records if c['eligible']]
bycomp=defaultdict(list)
for c in eligible: bycomp[c['homology_component']].append(c)
reps=[]
for comp,cs in sorted(bycomp.items()):
    # Fixed quality tie-breaker; do not sort by RMSD, observed minor state, or docking outcomes.
    rep=min(cs,key=lambda c:(max(float(c['res_apo']),float(c['res_holo'])),float(c['res_apo'])+float(c['res_holo']),c['apo'],c['holo'],c['holo_chain'],c['ligand_id'],c['apo_chain']))
    reps.append(rep)
rng=random.Random(20261028)
rng.shuffle(reps)
# Keep development structures separate; never silently relax eligibility to reach 24.
n_dev=3 if len(reps)>=8 else 1 if len(reps)>=3 else 0
dev=reps[:n_dev]; pilot=reps[n_dev:n_dev+24]; reserve=reps[n_dev+24:]
for c in records: c['allocation']='not_selected'
for label,subset in [('development',dev),('pilot',pilot),('reserve',reserve)]:
    for c in subset: c['allocation']=label
fields=['case_id','allocation','homology_component','group','apo','apo_chain','holo','holo_chain','ligand_id','ligand_icode','ligand_name','molecular_formula','res_apo','res_holo','pocket_backbone_rmsd','pocket_ca_rmsd','pocket_coverage','ligand_heavy_atoms','rotatable_bonds','formal_charge','ligand_min_occupancy','qfit_vs_deposited_ligand_rmsd','apo_mc_count','selected_count','changed_count','smiles']
def table(name,rs):
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
        for c in rs:
            d=dict(c,apo_mc_count=len(c['apo_mc_residues']),selected_count=len(c['selected_residues']),changed_count=len(c['changed_residues']))
            w.writerow(d)
table('eligible_cases.tsv',eligible);table('pilot_manifest.tsv',pilot);table('development_manifest.tsv',dev);table('reserve_manifest.tsv',reserve)
with (OUT/'exclusions.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['case_id','apo','holo','holo_chain','ligand_id','geometry_reasons','chemistry_reasons','exception'])
    for c in records:
        if not c['eligible']:w.writerow([c['case_id'],c['apo'],c['holo'],c['holo_chain'],c['ligand_id'],';'.join(c['geometry_reasons']),';'.join(c['chemistry_reasons']),c.get('audit_exception','')])
summary={'status':'selection_complete_under_defined_filters','rdkit_version':rdBase.rdkitVersion,'seed':20261028,'homology_screened_cases':len(records),'geometry_pass':sum(not c['geometry_reasons'] for c in records),'eligible_cases':len(eligible),'eligible_pairs':len({(c['apo'],c['holo']) for c in eligible}),'eligible_identical_sequence_groups':len({c['group'] for c in eligible}),'eligible_nonhomologous_components':len(bycomp),'development_cases':len(dev),'pilot_cases':len(pilot),'reserve_cases':len(reserve),'geometry_exclusion_reasons':dict(Counter(x for c in records for x in c['geometry_reasons'])),'chemistry_exclusion_reasons':dict(Counter(x for c in records if not c['geometry_reasons'] for x in c['chemistry_reasons'])),'pilot_apo_multiconformer_cases':sum(bool(c['apo_mc_residues']) for c in pilot),'pilot_selected_state_cases':sum(bool(c['selected_residues']) for c in pilot),'planned_pilot_docking_calls':len(pilot)*5*8*3,'inputs_sha256':{str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'remote_results_v3/geometry_v3_candidates.json',ROOT/'remote_results_v3/homology.json')}}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'all_records.json').write_text(json.dumps(records,indent=1)+'\n')
print(json.dumps(summary,indent=2))
