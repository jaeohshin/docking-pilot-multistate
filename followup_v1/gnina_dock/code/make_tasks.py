"""expand_v1에서 입력 준비가 끝난 사례 x arm 4개. 사례 목록은 expand_v1과 동일하다."""
import csv
from pathlib import Path
H = Path(__file__).resolve().parents[1]
rows = [r for r in csv.DictReader(open(H/'cohort_final.tsv'), delimiter='\t') if r['status'] in ('accept', 'backbone_change')]
ready = [r['case_id'] for r in rows if (H/'inputs'/r['case_id']/'ligand_dock.pdbqt').exists()]
arms = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
(H/'docked_tasks.txt').write_text('\n'.join(f'{c} {a}' for c in ready for a in arms) + '\n')
print('cohort', len(rows), '| prepared', len(ready), '| tasks', len(ready)*len(arms))
