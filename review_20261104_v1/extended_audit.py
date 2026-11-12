import argparse,json
from pathlib import Path
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();S=a.source;O=a.out
b=pd.read_csv(S/'followup_v1/budget50/docking_results.tsv',sep='\t');e=pd.read_csv(S/'expand_v1/docking_results.tsv',sep='\t');p=pd.read_csv(S/'docking_v1/docking_results.tsv',sep='\t');checks=[]
def boot(d,k,rng):
 cm=np.array([np.mean(x[k].to_numpy()) for _,x in d.groupby('component',sort=True)]);x=cm[rng.integers(0,len(cm),(10000,len(cm)))].mean(1);return round(float(cm.mean()),4),[round(float(z),4) for z in np.percentile(x,[2.5,97.5])]
rng=np.random.default_rng(20261028);summary=json.load(open(S/'expand_v1/docking_summary.json'))
for st in ['accept','backbone_change']:
 d=e[e.stratum==st]
 for name,k in [('apo_multi_minus_single','apo_gain_top'),('holo_multi_minus_single','holo_gain_top'),('apo_oracle','apo_gain_oracle'),('holo_oracle','holo_gain_oracle')]:
  m,ci=boot(d,k,rng);s=summary[st][name];assert m==s['mean_of_component_means'] and ci==s['ci95'],(st,name if 'name' in locals() else '',m,ci,s);checks.append(dict(check=st+' '+name,passed=True,mean=m,ci=ci))
 d=d[d.has_multiple_states==1];m,ci=boot(d,'apo_gain_top',rng);s=summary[st]['apo_multi_minus_single_states_only'];assert m==s['mean_of_component_means'] and ci==s['ci95'],(st,name if 'name' in locals() else '',m,ci,s);checks.append(dict(check=st+' states_only',passed=True,mean=m,ci=ci))
# Complete-case flow and group heterogeneity; no exclusions changed.
cm=b[b.stratum=='accept'].groupby('component').agg(cases=('case_id','size'),apo_top=('apo_gain_top','mean'),apo_top5=('apo_gain_top5','mean'),apo_oracle=('apo_gain_oracle','mean'),holo_oracle=('holo_gain_oracle','mean'))
cm['oracle_interaction']=cm.apo_oracle-cm.holo_oracle;cm.to_csv(O/'component_effects.tsv',sep='\t')
loo=[]
for comp in cm.index:
 d=cm.drop(index=comp);loo.append(dict(omitted_component=comp,component_mean_oracle=d.apo_oracle.mean(),component_mean_top1=d.apo_top.mean(),component_mean_top5=d.apo_top5.mean()))
pd.DataFrame(loo).to_csv(O/'leave_one_component.tsv',sep='\t',index=False)
# Actual combination number and collision attrition, all prepared analysed cases.
builds=[]
for cid in b.case_id:
 j=json.load(open(S/'expand_v1/receptors'/cid/'build.json'))
 for side in ['apo','holo']:
  x=j[side];builds.append(dict(case_id=cid,stratum=b.set_index('case_id').loc[cid,'stratum'],side=side,combinations_total=x['combinations_total'],kept_receptors=x['multi_receptors'],rejected_for_clash=x['rejected_for_clash'],multi_state_residues=len(x['multi_state_residues'])))
pd.DataFrame(builds).to_csv(O/'ensemble_inventory.tsv',sep='\t',index=False)
# Registered P1 CI exactly reproduced from initial rows, original RNG sequence.
rng=np.random.default_rng(20261028);s=json.load(open(S/'docking_v1/docking_summary.json'))
for st in ['main','backbone_change']:
 d=p[p.stratum==st]
 vals=[('apo_multi_minus_single',d.apo_gain_top),('holo_multi_minus_single',d.holo_gain_top),('apo_minus_holo_single',d.apo_single_rmsd_top-d.holo_single_rmsd_top),('interaction',d.apo_gain_top-d.holo_gain_top),('apo_multi_minus_single_states_only',d[d.has_multiple_states==1].apo_gain_top)]
 for name,v in vals:
  v=v.to_numpy();ci=[round(float(z),4) for z in np.percentile(v[rng.integers(0,len(v),(10000,len(v)))].mean(1),[2.5,97.5])];saved=s[st][name];assert round(v.mean(),4)==saved['mean'] and ci==saved['ci95'];checks.append(dict(check='pilot '+st+' '+name,passed=True,mean=round(v.mean(),4),ci=ci))
(O/'EXTENDED_AUDIT.json').write_text(json.dumps(dict(checks=checks,passed=True,leave_one_component_oracle_range=[min(x['component_mean_oracle'] for x in loo),max(x['component_mean_oracle'] for x in loo)],note='Pilot oracle interval is reported in RESULTS.md but was not exported by frozen analysis; exploratory rebootstrap provided separately.'),indent=2))
print('pass',len(checks),'LOO oracle range',min(x['component_mean_oracle'] for x in loo),max(x['component_mean_oracle'] for x in loo))
