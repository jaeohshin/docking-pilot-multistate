"""재채점 task: (사례, arm, 점수함수). pose가 이미 있는 것만."""
import argparse, csv
from pathlib import Path
a = argparse.ArgumentParser(); a.add_argument('--run', required=True); a.add_argument('--tag', required=True)
a.add_argument('--scoring', nargs='+', default=['vinardo', 'gnina']); a = a.parse_args()
H, run = Path(__file__).resolve().parents[1], Path(a.run)
sel = {r['case_id'] for r in csv.DictReader(open(run/'cohort_final.tsv'), delimiter='\t')
       if r['status'] in ('accept', 'backbone_change')}
arms = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
t = [f'{s} {c} {m}' for c in sorted(sel) for m in arms for s in a.scoring
     if (run/'docked'/c/m/'docking_log.json').exists()]
(H/f'tasks_{a.tag}.txt').write_text('\n'.join(t) + '\n')
print('tasks', len(t))
