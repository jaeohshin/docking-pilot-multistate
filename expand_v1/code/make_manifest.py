"""준비된 입력 파일에서 inputs_manifest.json을 재구성한다(run_docking.py가 읽는 형식).

build_and_prep.py가 prepare_case의 반환값을 버렸으므로, 디스크에 남은 산출물에서 되살린다:
수용체 목록은 inputs/<case>/<side>_<tag>.pdbqt, box는 Meeko가 쓴 <side>_single_0.box.txt에서 읽는다.
"""
import json, re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
cases, skipped = [], []
for d in sorted((HERE/'inputs').iterdir()):
    if not d.is_dir():
        continue
    rec = dict(case_id=d.name)
    ok = True
    for side in ('holo', 'apo'):
        box = d/f'{side}_single_0.box.txt'
        pdbqts = sorted(p.name for p in d.glob(f'{side}_*.pdbqt'))
        if not box.exists() or not pdbqts:
            ok = False; break
        v = dict(re.findall(r'(\w+)\s*=\s*(-?[\d.]+)', box.read_text()))
        receptors = {p[len(side) + 1:-len('.pdbqt')]: p for p in pdbqts}
        if 'single_0' not in receptors:
            ok = False; break
        rec[side] = dict(box_center=[float(v['center_x']), float(v['center_y']), float(v['center_z'])],
                         box_size=[float(v['size_x']), float(v['size_y']), float(v['size_z'])],
                         receptors=receptors)
    (cases if ok else skipped).append(rec if ok else d.name)
(HERE/'inputs_manifest.json').write_text(json.dumps(dict(cases=cases, incomplete=skipped), indent=1) + '\n')
print('cases', len(cases), '| incomplete', len(skipped))
import collections
print('arm별 수용체 수', collections.Counter(len(c[s]['receptors']) for c in cases for s in ('holo', 'apo')))
