"""확대 cohort의 상동성 component. 분석의 bootstrap 단위이며 도킹 결과와 무관하게 서열로만 계산한다.

기존 설계와 같은 조건: identity >= 0.30, max(query coverage, target coverage) >= 0.50, E <= 0.001,
sensitivity 7.5, 양방향. hit가 없는 query도 자기 자신만의 component로 남긴다(고립 query 누락 방지).
"""
import json, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
MMSEQS = '/gpfs/deepfold/others/mmseqs/bin/mmseqs'
TMP = HERE/'mmseqs_tmp'; TMP.mkdir(exist_ok=True)

cases = json.loads((HERE/'candidates.json').read_text())
want = {c['query_id'] for c in cases}
fa_in, keep, name = (HERE.parent/'run_v3/candidate_full_seqres.fasta').read_text(), [], None
for line in fa_in.splitlines():
    if line.startswith('>'):
        name = line[1:].strip(); continue
    if name in want:
        keep.append(f'>{name}\n{line}')
fa = TMP/'cohort.fasta'; fa.write_text('\n'.join(keep) + '\n')
print('sequences', len(keep), flush=True)

out = TMP/'hits.tsv'
subprocess.run([MMSEQS, 'easy-search', str(fa), str(fa), str(out), str(TMP/'work'),
                '-s', '7.5', '-e', '0.001', '--max-seqs', '4000', '--threads', '4',
                '--format-output', 'query,target,fident,qcov,tcov,evalue'], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

parent = {q: q for q in want}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[max(ra, rb)] = min(ra, rb)

edges = 0
for line in out.read_text().splitlines():
    q, t, fid, qcov, tcov, ev = line.split('\t')
    if q == t or q not in parent or t not in parent:
        continue
    if float(fid) >= 0.30 and max(float(qcov), float(tcov)) >= 0.50 and float(ev) <= 0.001:
        union(q, t); edges += 1

comp = {q: find(q) for q in want}
sizes = {}
for q, r in comp.items():
    sizes[r] = sizes.get(r, 0) + 1
for c in cases:
    c['homology_component'] = comp[c['query_id']]
(HERE/'candidates.json').write_text(json.dumps(cases, indent=1) + '\n')
summary = dict(queries=len(want), hit_edges=edges, components=len(set(comp.values())),
               largest_component=max(sizes.values()), singletons=sum(1 for v in sizes.values() if v == 1),
               cases=len(cases), cases_per_component_max=max(
                   sum(1 for c in cases if c['homology_component'] == r) for r in set(comp.values())))
(HERE/'components.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
