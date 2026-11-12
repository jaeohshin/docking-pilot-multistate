import json,hashlib,argparse
from pathlib import Path
import pandas as pd,numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();S=a.source;O=a.out
B=pd.read_csv(S/'followup_v1/budget50/docking_results.tsv',sep='\t');V=pd.read_csv(S/'followup_v1/budget50/pose_table.tsv.gz',sep='\t');E=pd.read_csv(S/'followup_v1/rescore/out/expand8/pose_table.tsv.gz',sep='\t');R=pd.read_csv(S/'followup_v1/rescore/out/b50/pose_table.tsv.gz',sep='\t')
K=['case_id','arm','call','pose'];assert not E.duplicated(K).any();v=V[V.call<8].set_index(K).sort_index();e=E.set_index(K).sort_index();assert v.index.equals(e.index);assert np.array_equal(v[['vina_score','rmsd']],e[['vina_score','rmsd']]);assert V.set_index(K).sort_index().equals(R.set_index(K).sort_index())
rows=[]
for run,n in [('expand_v1',8),('followup_v1/budget50',50),('followup_v1/gnina_dock',8)]:
 for cid in B.case_id:
  for arm in ['apo_single','apo_multi','holo_single','holo_multi']:
   log=json.load(open(S/run/'docked'/cid/arm/'docking_log.json'));calls=log['calls'];assert len(calls)==n;assert [c['run'] for c in calls]==list(range(n));assert [c['seed'] for c in calls]==[20261028+k for k in range(n)];assert all(c['returncode']==0 for c in calls);assert log['exhaustiveness']==32 and log['modes']==9
   assert all(c['receptor']==log['receptors'][k%len(log['receptors'])] for k,c in enumerate(calls))
   row=B.set_index('case_id').loc[cid];assert len(log['receptors'])==row[arm+'_receptors']
   rows.append(dict(run=run,case_id=cid,arm=arm,calls=n,returncode_failures=0,receptors=len(log['receptors']),seconds=sum(c['seconds'] for c in calls)))
# Full frozen table top-1/raw top-5 reconstruction, including score rounding/tie effects.
errs=[]
for run,pose,results,sel,pre in [('Vina50',V,B,'vina_score',''),('GNINA8',pd.read_csv(S/'followup_v1/gnina_dock/pose_table.tsv.gz',sep='\t'),pd.read_csv(S/'followup_v1/gnina_dock/docking_results.tsv',sep='\t'),'cnnscore','cnnscore_')]:
 res=results.set_index('case_id')
 for (cid,arm),d in pose.groupby(['case_id','arm'],sort=False):
  d=d.sort_values([sel,'call','pose'],kind='stable')
  for met,calc in [('top',d.iloc[0].rmsd),('top5_raw',d.iloc[:5].rmsd.min())]:
   orig=res.loc[cid,arm+'_'+pre+'rmsd_'+met]
   if calc!=orig:errs.append(dict(engine=run,case_id=cid,arm=arm,metric=met,from_rounded_pose_table=calc,frozen_full_precision=orig))
# Single apo receptor control, strict key equality of scores/RMSD between arms.
controls=[]
for cid in B[(B.stratum=='accept')&(B.apo_multi_receptors==1)].case_id:
 x=V[(V.case_id==cid)&(V.arm=='apo_single')].set_index(['call','pose'])[['vina_score','rmsd']].sort_index();y=V[(V.case_id==cid)&(V.arm=='apo_multi')].set_index(['call','pose'])[['vina_score','rmsd']].sort_index();assert x.equals(y);controls.append(cid)
pd.DataFrame(rows).to_csv(O/'log_audit.tsv',sep='\t',index=False);pd.DataFrame(errs,columns=['engine','case_id','arm','metric','from_rounded_pose_table','frozen_full_precision']).to_csv(O/'rounded_score_tie_differences.tsv',sep='\t',index=False)
(O/'LOG_AUDIT.json').write_text(json.dumps(dict(logs_checked=len(rows),failed_calls=0,first8_pose_rows_identical=len(v),b50_rescore_pose_table_identical=True,apo_one_receptor_controls_exact=len(controls),rounded_score_tie_differences=errs,total_logged_seconds_by_run=pd.DataFrame(rows).groupby('run').seconds.sum().to_dict()),indent=2))
print('logs',len(rows),'first8',len(v),'score-tie differences',errs)
