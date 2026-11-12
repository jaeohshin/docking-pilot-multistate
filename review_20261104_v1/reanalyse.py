#!/usr/bin/env python3
"""Exploratory audit of frozen docking tables; no docking/inference and no source writes.
Usage: python reanalyse.py --source ../source --out .
"""
import argparse,json,hashlib,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
a=ap.parse_args();S=a.source;O=a.out;O.mkdir(exist_ok=True,parents=True)
ARMS=['holo_single','holo_multi','apo_single','apo_multi'];SEED=20261028
read=lambda p:pd.read_csv(S/p,sep='\t')
load=lambda p:json.loads((S/p).read_text())
save=lambda x,n:pd.DataFrame(x).to_csv(O/n,sep='\t',index=False)
B=read('followup_v1/budget50/docking_results.tsv');C=read('followup_v1/gnina_dock/docking_results.tsv');E=read('expand_v1/docking_results.tsv');P=read('docking_v1/docking_results.tsv')
R=read('followup_v1/rescore/out/b50/rescore_results.tsv')
V=read('followup_v1/budget50/pose_table.tsv.gz');G=read('followup_v1/gnina_dock/pose_table.tsv.gz')
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
 assert ok,(name,detail)
def component(d,key,rng,n=10000,q=(2.5,97.5)):
 cm=d.groupby('component',sort=True)[key].mean().values
 draws=cm[rng.integers(0,len(cm),(n,len(cm)))].mean(axis=1)
 return float(cm.mean()),np.percentile(draws,q).round(4).tolist()
def contrast_check(d,key,saved,rng,coprimary=False):
 m,ci=component(d,key,rng)
 check(key+' mean',round(m,4)==saved['mean_of_component_means'],(m,saved))
 check(key+' ci95',ci==saved['ci95'],ci)
 if coprimary:
  _,ci=component(d,key,rng,q=(1.25,98.75));check(key+' ci975',ci==saved['ci975_bonferroni'],ci)
# Independent groupby + vectorized implementation, same RNG consumption as frozen analyses.
rng=np.random.default_rng(SEED);sb=load('followup_v1/budget50/docking_summary.json')
for st in ['accept','backbone_change']:
 d=B[B.stratum==st];s=sb[st]
 for name,key,co in [('P1_apo_multi_minus_single_top1','apo_gain_top',True),('P2_apo_multi_minus_single_top5','apo_gain_top5',True),('apo_top5_raw','apo_gain_top5_raw',False),('holo_multi_minus_single_top1','holo_gain_top',False),('holo_multi_minus_single_top5','holo_gain_top5',False),('apo_oracle','apo_gain_oracle',False),('holo_oracle','holo_gain_oracle',False),('apo_frac_good','apo_gain_frac_good',False),('holo_frac_good','holo_gain_frac_good',False),('apo_calls_good','apo_gain_calls_good',False),('holo_calls_good','holo_gain_calls_good',False)]:
  contrast_check(d,key,s[name],rng,co)
 sub=d[d.has_multiple_states==1]
 for name,key in [('apo_top1_states_only','apo_gain_top'),('apo_top5_states_only','apo_gain_top5')]:contrast_check(sub,key,s[name],rng)
rng=np.random.default_rng(SEED);sc=load('followup_v1/gnina_dock/docking_summary.json')
for st in ['accept','backbone_change']:
 d=C[C.stratum==st]
 for side in ['holo','apo']:
  contrast_check(d,side+'_gain_oracle',sc[st][side+'_oracle'],rng)
  for sel in ['cnnscore','cnnaffinity','minaffinity']:
   for met in ['top','top5','top5_raw']:contrast_check(d,f'{side}_{sel}_gain_{met}',sc[st][f'{side}_multi_minus_single_{sel}_{met}'],rng)
rng=np.random.default_rng(SEED);sr=load('followup_v1/rescore/out/b50/rescore_summary.json')
for st in ['accept','backbone_change']:
 d=R[R.stratum==st]
 for sel in ['vina','vinardo','gnina']:
  for side in ['holo','apo']:
   for lbl,met in [('top1','top'),('top5_raw','top5_raw')]:contrast_check(d,f'{sel}_{side}_gain_{met}',sr[st][sel][f'{side}_multi_minus_single_{lbl}'],rng)
# Key integrity and reconstruction from pose tables.
key=['case_id','arm','call','pose']
for name,t,d,runs,selectors in [('vina50',V,B,50,['vina_score']),('gnina8',G,C,8,['cnnscore','cnnaffinity','minaffinity'])]:
 check(name+' unique pose keys',not t.duplicated(key).any())
 check(name+' complete case sets',set(t.case_id)==set(d.case_id))
 check(name+' finite scores/RMSD',np.isfinite(t[['rmsd']+selectors].values).all())
 for (cid,arm),g in t.groupby(['case_id','arm'],sort=False):
  row=d.set_index('case_id').loc[cid]
  check(name+' '+cid+' '+arm+' N/calls',len(g)==row[arm+'_poses'] and g.call.nunique()==runs)
  check(name+' '+cid+' '+arm+' oracle',abs(g.rmsd.min()-row[arm+'_rmsd_oracle'])<1e-9)
# Threshold rounding audit: frozen RMSDs are rounded to .001 A.
threshold=V[V.rmsd==2].assign(engine='vina50');threshold=pd.concat([threshold,G[G.rmsd==2].assign(engine='gnina8')]);save(threshold,'threshold_boundary.tsv')
# Preserve original confirmatory results; supplementary subsets and interactions are all exploratory.
main=B[B.stratum=='accept'].copy();check('main 128/39',len(main)==128 and main.component.nunique()==39)
cohort=read('followup_v1/budget50/cohort_final.tsv');cand={c['case_id']:c['homology_component'] for c in load('followup_v1/budget50/candidates.json')}
check('main cohort membership subset',set(main.case_id)<=set(cohort[cohort.status=='accept'].case_id))
missing=cohort[(cohort.status.isin(['accept','backbone_change']))&(~cohort.case_id.isin(B.case_id))].copy()
failed={}
for f in (S/'expand_v1').glob('build_prep_shard*.json'):
 for r in json.loads(f.read_text()).get('failed',[]):failed[r['case_id']]=r['error']
missing['recorded_preparation_failure']=missing.case_id.map(failed);save(missing,'eligible_not_analysed.tsv')
check('component keys',all(cand[r.case_id]==r.component for r in B.itertuples()))
check('experiments share analysed cohorts',set(B.case_id)==set(C.case_id)==set(E.case_id)==set(R.case_id))
summary=[];trans=[];effects=[];subrows=[]
for side in ['apo','holo']:
 for met,succ in [('top','success_top'),('top5','success_top5'),('oracle','success_any')]:
  x=main[f'{side}_single_{succ}'];y=main[f'{side}_multi_{succ}'];trans.append(dict(side=side,metric=met,both_fail=int(((x==0)&(y==0)).sum()),gain=int(((x==0)&(y==1)).sum()),loss=int(((x==1)&(y==0)).sum()),both_success=int(((x==1)&(y==1)).sum()),single_success=int(x.sum()),multi_success=int(y.sum())))
for label,d in [('all128',main),('apo_one',main[main.apo_multi_receptors==1]),('apo_multiple',main[main.apo_multi_receptors>=2])]:
 for side in ['apo','holo']:
  for met in ['top','top5','top5_raw','oracle']:
   k=f'{side}_gain_{met}';m,ci=component(d,k,np.random.default_rng(20261104),n=20000)
   subrows.append(dict(subset=label,side=side,metric=met,cases=len(d),components=d.component.nunique(),case_mean=d[k].mean(),component_mean=m,ci95_low=ci[0],ci95_high=ci[1],single_case_mean=d[f'{side}_single_rmsd_{met}'].mean(),multi_case_mean=d[f'{side}_multi_rmsd_{met}'].mean()))
 for met in ['top','top5','oracle']:
  d=d.copy();d['interaction']=d['apo_gain_'+met]-d['holo_gain_'+met];m,ci=component(d,'interaction',np.random.default_rng(20261104),n=20000)
  effects.append(dict(subset=label,metric=met,definition='apo gain minus holo gain',cases=len(d),components=d.component.nunique(),component_mean=m,ci95_low=ci[0],ci95_high=ci[1]))
for name,d,pre in [('Vina50',main,''),('GNINA8',C[C.stratum=='accept'],'cnnscore_')]:
 for arm in ARMS:
  summary.append(dict(engine=name,arm=arm,cases=len(d),top1=int(d[arm+'_'+pre+'success_top'].sum()),distinct_top5=int(d[arm+'_'+pre+'success_top5'].sum()),oracle=int(d[arm+'_success_any'].sum()),top1_case_mean=d[arm+'_'+pre+'rmsd_top'].mean(),oracle_case_mean=d[arm+'_rmsd_oracle'].mean()))
save(trans,'success_transitions.tsv');save(subrows,'subgroup_effects.tsv');save(effects,'interaction_effects.tsv');save(summary,'absolute_results.tsv')
# No independence assertion for nested cohorts: exact case IDs and biological key overlaps.
p=P[P.stratum=='main'];overlap=p.merge(main[['case_id','component']],on='case_id',how='left',validate='one_to_one');save(overlap[['case_id','apo','holo','ligand','component']],'initial_overlap.tsv')
check('pilot overlap',overlap.component.notna().all(),len(overlap))
for met in ['top','oracle']:
 v=p['apo_gain_'+met].values;draw=v[np.random.default_rng(20261104).integers(0,len(v),(20000,len(v)))].mean(1)
 effects.append(dict(subset='pilot13',metric=met,definition='apo gain',cases=len(p),components=len(p),component_mean=v.mean(),ci95_low=np.percentile(draw,2.5),ci95_high=np.percentile(draw,97.5)))
save(effects,'interaction_effects.tsv')
# Restore exact labels only at the .001 A export boundary, and full SDF score precision.
precision=pd.read_csv(O/'precision_evidence.tsv',sep='\t')
key=['case_id','arm','call','pose']
for engine,table in [('Vina50',V),('GNINA8',G)]:
 z=precision[(precision.engine==engine)&precision.rmsd_exact_boundary.notna()]
 for r in z.itertuples():
  mask=(table.case_id==r.case_id)&(table.arm==r.arm)&(table.call==r.call)&(table.pose==r.pose)
  check(engine+' boundary key '+str((r.case_id,r.arm,r.call,r.pose)),mask.sum()==1)
  table.loc[mask,'rmsd']=r.rmsd_exact_boundary
z=precision[(precision.engine=='GNINA8')&precision.cnnscore.notna()]
check('GNINA full precision score key equality',set(map(tuple,z[key].values))==set(map(tuple,G[key].values)))
G=G.drop(columns=['cnnscore','cnnaffinity','minaffinity']).merge(z[key+['cnnscore','cnnaffinity','minaffinity']],on=key,validate='one_to_one',how='left')
# Score merge uses explicit four-part keys, never row position.
vs=V.copy()
for sel,dirname,col in [('gnina','gnina','score'),('gnina_affinity','gnina','score_alt'),('vinardo','vinardo','score')]:
 ss=[]
 for f in sorted((S/'followup_v1/rescore/out/b50'/dirname).glob('*__*.tsv')):
  cid,arm=f.stem.split('__');z=pd.read_csv(f,sep='\t');z['case_id']=cid;z['arm']=arm;ss.append(z[key+[col]].rename(columns={col:sel}))
 s=pd.concat(ss);check(sel+' score keys unique',not s.duplicated(key).any());check(sel+' score key equality',set(map(tuple,s[key].values))==set(map(tuple,V[key].values)))
 vs=vs.merge(s,on=key,how='left',validate='one_to_one',sort=False)
# Rank of FIRST <=2 A, conditional on G>0. Exact random-without-replacement null.
# P(R<=r)=1-C(N-r,G)/C(N,G); scipy hypergeom.sf(0,N,G,r).
# Mid-CDF = P(R<observed)+.5 P(R=observed); random expectation .5, small is better.
def ranks(table,engine,selectors):
 rows=[]
 for (cid,arm),g0 in table.groupby(['case_id','arm'],sort=True):
  if cid not in set(main.case_id):continue
  for sel in selectors:
   g=g0.sort_values([sel,'call','pose'],kind='stable');good=(g.rmsd<=2).to_numpy();N=len(g);G0=int(good.sum());row=dict(engine=engine,selector=sel,case_id=cid,component=cand[cid],arm=arm,N=N,G=G0,good_exists=G0>0,top1_success=int(good[0]))
   if G0:
    idx=int(np.flatnonzero(good)[0]);r=idx+1;score=g.iloc[idx][sel];tie=g[sel]==score;T=int(tie.sum());gt=int((tie&(g.rmsd<=2)).sum());before=int((g[sel]<score).sum());E0=(N+1)/(G0+1)
    F=hypergeom.sf(0,N,G0,r);Fp=hypergeom.sf(0,N,G0,r-1)
    row.update(first_good_rank=r,random_expected_rank=E0,rank_over_random_expectation=r/E0,random_midcdf=(F+Fp)/2,random_prob_rank_le_observed=F,first_good_tie_expected_rank=before+(T+1)/(gt+1),first_good_tie_min=before+1,first_good_tie_max=before+T-gt+1,random_top1_probability=G0/N)
   rows.append(row)
 return rows
rank=pd.DataFrame(ranks(vs,'Vina50',['vina_score','gnina','gnina_affinity','vinardo'])+ranks(G,'GNINA8',['cnnscore','cnnaffinity','minaffinity']));save(rank,'first_good_rank_per_case.tsv')
rs=[]
for (eng,sel,arm),d in rank.groupby(['engine','selector','arm']):
 g=d[d.good_exists];m,ci=component(g,'random_midcdf',np.random.default_rng(20261104),n=20000)
 rs.append(dict(engine=eng,selector=sel,arm=arm,cases=len(d),generated_good=len(g),no_good=len(d)-len(g),components_good=g.component.nunique(),selected_good=int(g.top1_success.sum()),first_rank_median=g.first_good_rank.median(),first_rank_mean=g.first_good_rank.mean(),random_expected_rank_mean=g.random_expected_rank.mean(),random_expected_rank_median=g.random_expected_rank.median(),ratio_median=g.rank_over_random_expectation.median(),mean_midcdf_case=g.random_midcdf.mean(),mean_midcdf_component=m,midcdf_component_ci_low=ci[0],midcdf_component_ci_high=ci[1],expected_random_top1_success=g.random_top1_probability.sum(),tie_expected_rank_median=g.first_good_tie_expected_rank.median(),cases_with_tie_range=int((g.first_good_tie_max>g.first_good_tie_min).sum())))
save(rs,'first_good_rank_summary.tsv')
# Explicit random permutation verification for small toy sets (all combinations exact).
for N,ng in [(5,1),(5,2),(6,5),(5,5)]:
 mins=np.array([min(c)+1 for c in itertools.combinations(range(N),ng)])
 check('rank null expectation '+str((N,ng)),abs(mins.mean()-(N+1)/(ng+1))<1e-12)
 q=.5*(hypergeom.sf(0,N,ng,mins)+hypergeom.sf(0,N,ng,mins-1));check('midCDF null '+str((N,ng)),abs(q.mean()-.5)<1e-12)
# Observed selected-pose success against random first pose; no random independent-pose significance test.
fig,ax=plt.subplots(1,2,figsize=(10,4))
for i,side in enumerate(['apo','holo']):
 for label,sel,color in [('Vina','vina_score','#377eb8'),('CNN rescore','gnina','#e41a1c')]:
  d=rank[(rank.engine=='Vina50')&(rank.selector==sel)&(rank.arm==side+'_multi')&rank.good_exists]
  ax[i].scatter(d.random_expected_rank,d.first_good_rank,label=label,alpha=.65,s=18,c=color)
 ax[i].plot([1,450],[1,450],'k--',lw=1,label='Observed = random expectation');ax[i].set(xscale='log',yscale='log',xlabel='Random first-good rank: (N+1)/(G+1)',ylabel='Observed first-good rank',title=side+' multi, conditional on G > 0');ax[i].legend(fontsize=7)
fig.suptitle('Exploratory: raw pose lists; duplicates retained');fig.tight_layout();fig.savefig(O/'first_good_rank.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.5,4));labels=['Top-1','Distinct top-5','Oracle'];x=np.arange(3)
t=pd.DataFrame(trans).query("side=='apo'")
ax.bar(x-.18,t.single_success,.36,label='Single',color='#377eb8');ax.bar(x+.18,t.multi_success,.36,label='Multiple',color='#e41a1c');ax.set(xticks=x,xticklabels=labels,ylabel='Successful cases / 128',ylim=(0,80),title='Vina, 50 calls: generation and selection');ax.legend()
for i,r in enumerate(t.itertuples()):ax.text(i,max(r.single_success,r.multi_success)+2,f'+{r.gain} / −{r.loss}',ha='center')
fig.tight_layout();fig.savefig(O/'generation_selection.png',dpi=180);plt.close(fig)
save(checks,'checks.tsv')
(O/'TABLE_AUDIT.json').write_text(json.dumps(dict(exploratory=True,checks=len(checks),all_passed=all(x['passed'] for x in checks),threshold_boundary_rows=len(threshold),pilot_cases_in_main=int(overlap.component.notna().sum()),apo_one=int((main.apo_multi_receptors==1).sum()),apo_multiple=int((main.apo_multi_receptors>=2).sum()),note='New CIs 20,000 component resamples, seed 20261104; frozen CI exact reproduction uses 10,000 and original sequential RNG. Raw-pose null is conditional on saved candidate lists, not a claim of independent biological samples.'),indent=2))
print(pd.DataFrame(rs).query("arm in ['apo_single','apo_multi'] and selector in ['vina_score','gnina','cnnscore']").to_string(index=False));print(pd.DataFrame(effects).to_string(index=False));print('all',len(checks),'checks passed')
