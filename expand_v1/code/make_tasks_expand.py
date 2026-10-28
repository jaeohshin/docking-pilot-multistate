"""입력 준비가 끝난 사례 × arm 4개를 도킹 task 목록으로 쓴다."""
import json, csv
from pathlib import Path
H = Path(__file__).resolve().parents[1]
rows = [r for r in csv.DictReader(open(H/'cohort_final.tsv'), delimiter='\t') if r['status'] in ('accept','backbone_change')]
ready = [r['case_id'] for r in rows if (H/'inputs'/r['case_id']/'ligand_dock.pdbqt').exists()]
arms = ['holo_single','holo_multi','apo_single','apo_multi']
(H/'docked_tasks.txt').write_text('\n'.join(f'{c} {a}' for c in ready for a in arms) + '\n')
print('cases in cohort', len(rows), '| prepared', len(ready), '| tasks', len(ready)*len(arms))
