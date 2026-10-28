"""One line per (case, arm) for the docking array, in a fixed order."""
import json
from pathlib import Path
H = Path(__file__).resolve().parents[1]
cases = [c['case_id'] for c in json.loads((H/'inputs_manifest.json').read_text())['cases']]
arms = ['holo_single', 'holo_multi', 'apo_single', 'apo_multi']
(H/'docked_tasks.txt').write_text('\n'.join(f'{c} {a}' for c in cases for a in arms) + '\n')
print(len(cases), 'cases x', len(arms), 'arms =', len(cases)*len(arms), 'tasks')
