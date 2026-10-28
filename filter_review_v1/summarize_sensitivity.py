import json, csv, hashlib
from pathlib import Path
from collections import Counter
import gemmi
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

OUT=Path(__file__).resolve().parent
BASE=OUT.parent if (OUT.parent/'selection_v3').exists() else Path('/data/work/sidechain-multistate/experiments/docking_pilot_20261028')
rows=json.loads((OUT/'sensitivity_candidates.json').read_text())
old=json.loads((BASE/'selection_v3/all_records.json').read_text())
hom=json.loads((BASE/'remote_results_v3/homology.json').read_text())
assert len(rows)==443 and len({r['case_id'] for r in rows})==443
assert {r['case_id'] for r in rows}=={r['case_id'] for r in old}
parent={r['query_id']:r['query_id'] for r in rows}
def root(x):
    parent.setdefault(x,x)
    if parent[x]!=x: parent[x]=root(parent[x])
    return parent[x]
for h in hom['edges']:
    a,b=root(h[0]),root(h[1]);parent[max(a,b)]=min(a,b)
ENV={f'{side}_{flag}' for side in ('apo','holo') for flag in ('multichain_pocket','crystal_mate_pocket','nearby_metal','nearby_other_heterocomponent')}
DESIGN={'heavy_atom_count_outside_10_50','ligand_occupancy_below_0_8','ccd_not_nonpolymer','pocket_backbone_rmsd_above_0_75'}
SOFT=ENV|DESIGN
def chemistry(c):
    if set(c['geometry_reasons'])-SOFT: return ['not_evaluated_hard_geometry_failure']
    assert 'deposited_cif_sha256' in c and 'qfit_vs_deposited_ligand_rmsd' in c
    smiles=[]
    for descriptor in c['ccd']['descriptors']:
        typ,prog,s=map(gemmi.cif.as_string,descriptor)
        if typ.startswith('SMILES'): smiles.append((typ!='SMILES_CANONICAL',prog,s))
    m=None
    for _,prog,s in sorted(smiles):
        candidate=Chem.MolFromSmiles(s)
        if candidate is not None:
            m=Chem.RemoveHs(candidate);c['smiles']=s;break
    if m is None:return ['ccd_smiles_unreadable']
    c['rdkit_heavy_atoms']=m.GetNumHeavyAtoms()
    c['rotatable_bonds']=rdMolDescriptors.CalcNumRotatableBonds(m,strict=True)
    c['formal_charge']=Chem.GetFormalCharge(m)
    c['ligand_name']=c['ccd']['name']
    reasons=[]
    if m.GetNumHeavyAtoms()!=c['ligand_heavy_atoms']:reasons.append('ccd_smiles_structure_heavy_count_mismatch')
    if c['rotatable_bonds']>12:reasons.append('rotatable_bonds_above_12')
    if len(Chem.GetMolFrags(m))!=1:reasons.append('disconnected_ligand_graph')
    if not any(a.GetAtomicNum()==6 for a in m.GetAtoms()):reasons.append('no_carbon')
    if any(a.GetAtomicNum() not in {1,5,6,7,8,9,15,16,17,35,53} for a in m.GetAtoms()):reasons.append('metal_or_unsupported_element')
    if c['ccd']['id'] in {'GSH','GSS','GSF','GTT'}:reasons.append('glutathione_cofactor_family')
    return reasons
for c in rows:
    c['homology_component']=root(c['query_id'])
    c['chemistry_reasons']=chemistry(c)
def count(rs):
    return dict(cases=len(rs),pairs=len({(r['apo'],r['holo']) for r in rs}),source_groups=len({r['group'] for r in rs}),homology_components=len({root(r['query_id']) for r in rs}),apo_mc_cases=sum(bool(r['apo_mc_residues']) for r in rs),selected_state_cases=sum(bool(r['selected_residues']) for r in rs))
def select(env=(),occ=.8,size=10,rmsd=.75,ccd=True,exclude_additives=False):
    return [r for r in rows if not r['chemistry_reasons'] and not(set(r['geometry_reasons'])-SOFT) and not(set(r['geometry_reasons'])&(ENV-set(env))) and r['ligand_min_occupancy']>=occ and size<=r['ligand_heavy_atoms']<=50 and r['pocket_backbone_rmsd']<=rmsd and (not ccd or r['ccd']['type'].upper()=='NON-POLYMER') and (not exclude_additives or r['ccd']['id'] not in {'B3P','BTB','9JE'})]
scenarios={
    'strict':{},
    'crystal_only':{'env':{x for x in ENV if 'crystal' in x}},
    'other_chain_only':{'env':{x for x in ENV if 'multichain' in x}},
    'other_hetero_only':{'env':{x for x in ENV if 'heterocomponent' in x}},
    'metal_proximity_only':{'env':{x for x in ENV if 'metal' in x}},
    'occupancy_only':{'occ':.5},
    'size_only':{'size':6},
    'rmsd_only':{'rmsd':1.5},
    'environment_no_metal':{'env':{x for x in ENV if 'metal' not in x}},
    'environment_all':{'env':ENV},
    'environment_occ':{'env':ENV,'occ':.5},
    'environment_occ_size':{'env':ENV,'occ':.5,'size':6},
    'environment_occ_size_rmsd':{'env':ENV,'occ':.5,'size':6,'rmsd':1.5},
    'environment_occ_size_rmsd_ccdtype':{'env':ENV,'occ':.5,'size':6,'rmsd':1.5,'ccd':False},
    'simple_occ_size_rmsd':{'occ':.5,'size':6,'rmsd':1.5},
    'no_metal_occ_size_rmsd':{'env':{x for x in ENV if 'metal' not in x},'occ':.5,'size':6,'rmsd':1.5},
    'no_metal_occ_size':{'env':{x for x in ENV if 'metal' not in x},'occ':.5,'size':6},
    'expanded_exclude_candidate_additives':{'env':ENV,'occ':.5,'size':6,'rmsd':1.5,'exclude_additives':True},
}
summary={name:count(select(**args)) for name,args in scenarios.items()}
assert {r['case_id'] for r in select()}=={r['case_id'] for r in old if r['eligible']}
selected_ids={n:{r['case_id'] for r in select(**a)} for n,a in scenarios.items()}
for c in rows:c['review_scenarios']=[n for n,ids in selected_ids.items() if c['case_id'] in ids]
stages=[('mapping',{'qfit_cache_sequence_mismatch'}),('size_or_occupancy',{'heavy_atom_count_outside_10_50','ligand_occupancy_below_0_8'}),('ccd_type_or_atom_completeness',{'ccd_not_nonpolymer','ligand_heavy_atom_names_incomplete_or_extra'}),('pocket_mapping',{'pocket_coverage_below_95pct','insufficient_pocket_atoms','contact_backbone_mapping_or_atoms_missing'}),('backbone_or_junction_or_box',{'pocket_backbone_rmsd_above_0_75','apo_pocket_junction_gap_or_distortion','holo_pocket_junction_gap_or_distortion','box_edge_above_30'}),('environment',ENV)]
remaining=old.copy();attrition=[]
for name,reasons in stages:
    drop=[r for r in remaining if set(r['geometry_reasons'])&reasons]
    remaining=[r for r in remaining if not(set(r['geometry_reasons'])&reasons)]
    attrition.append(dict(stage=name,excluded_cases=len(drop),remaining=count(remaining)))
drop=[r for r in remaining if r['geometry_reasons'] or r['chemistry_reasons']]
remaining=[r for r in remaining if not r['geometry_reasons'] and not r['chemistry_reasons']]
attrition.append(dict(stage='deposition_covalence_chemistry',excluded_cases=len(drop),remaining=count(remaining)))
assert len(remaining)==12
result=dict(status='sensitivity_audit_complete_preparation_review_pending',rdkit_version=rdBase.rdkitVersion,initial=count(rows),sequential_original_attrition=attrition,scenarios=summary,chemistry_reasons=dict(Counter(x for r in rows for x in r['chemistry_reasons'])),inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (OUT/'sensitivity_candidates.json',BASE/'selection_v3/all_records.json',BASE/'remote_results_v3/homology.json')})
(OUT/'review_summary.json').write_text(json.dumps(result,indent=2)+'\n')
(OUT/'review_records.json').write_text(json.dumps(rows,indent=1)+'\n')
with (OUT/'review_candidates.tsv').open('w') as f:
    fields=['case_id','homology_component','apo','apo_chain','holo','holo_chain','ligand_id','ligand_name','ligand_heavy_atoms','ligand_min_occupancy','pocket_backbone_rmsd','apo_mc_count','selected_count','environment_flags','review_scenarios']
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
    for r in rows:
        if r['review_scenarios']:w.writerow(dict(r,apo_mc_count=len(r['apo_mc_residues']),selected_count=len(r['selected_residues']),environment_flags=';'.join(sorted(set(r['geometry_reasons'])&ENV)),review_scenarios=';'.join(r['review_scenarios'])))
print(json.dumps(result,indent=2))
