"""확대 cohort(233건)를 기존 단계 코드로 그대로 통과시킨다.

기존 스크립트는 79건 전용 경로가 코드에 박혀 있다. 그 파일들을 고치면 이미 끝난 결과의 출처가 흐려지므로,
여기서는 **모듈로 불러와 경로 상수만 바꾼 뒤** 각자의 main()을 호출한다. 계산 코드는 한 줄도 바뀌지 않는다.

단계마다 출력 파일이 있으면 건너뛴다(중단 후 이어받기). 79건과 달리 대표를 고르지 않고 **전 사례를 도킹**하므로
상동성 component는 이 단계에서 필요 없다(분석에서 bootstrap 단위로만 쓰인다).
"""
import importlib.util, json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent
CANDS = HERE/'candidates.json'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def stage(out, label, fn):
    if (HERE/out).exists():
        print(f'skip {label} ({out} exists)', flush=True)
        return
    t = time.time()
    print(f'--- {label}', flush=True)
    fn()
    print(f'--- {label} done in {time.time()-t:.0f}s', flush=True)


def assembly():
    m = load('aac', PILOT/'final_selection_v1/code/audit_assembly_contacts.py')
    m.CANDS, m.CIFS, m.HERE = CANDS, HERE/'cifs', HERE
    m.VAL = HERE/'validation_xml'                      # 비어 있으면 RSCC 없음으로 기록된다
    cases = json.loads(CANDS.read_text())
    m.json = json                                       # (모듈 네임스페이스 유지용 no-op)
    comp = {c['case_id']: None for c in cases}
    res, errors = [], []
    for k, c in enumerate(cases):
        try:
            res.append(m.audit(c, comp))
        except Exception as e:
            errors.append(dict(case_id=c['case_id'], error=repr(e)))
            print('ERROR', c['case_id'], repr(e)[:200], flush=True)
        if k % 25 == 0:
            print(' assembly', k, '/', len(cases), 'errors', len(errors), flush=True)
    (HERE/'assembly_contact_audit.json').write_text(json.dumps(
        dict(cases=res, errors=errors, code_sha256=m.sha(PILOT/'final_selection_v1/code/audit_assembly_contacts.py'),
             gemmi=m.gemmi.__version__), indent=1, default=float) + '\n')
    print('assembly: ok', len(res), 'errors', len(errors))


def states():
    m = load('aqs', PILOT/'state_audit_v1/code/audit_qfit_states.py')
    m.HERE, m.PILOT = HERE, PILOT
    sys.argv = ['audit_qfit_states.py']
    # 입력 경로만 바꾼 main: 79건 파일 대신 확대 cohort와 이 폴더의 감사 결과를 읽는다
    C = {c['case_id']: c for c in json.loads(CANDS.read_text())}
    AUD = {o['case_id']: o for o in json.loads((HERE/'assembly_contact_audit.json').read_text())['cases']}
    out, cache, need = [], {}, {c[k] for c in C.values() for k in ('apo', 'holo')}
    import glob
    for f in sorted(glob.glob(str(PILOT.parents[1]/'reference/frozen/data/qfit_at_scale/replicate_cache/*.jsonl'))):
        for line in open(f):
            if line[12:16] in need:
                e = json.loads(line); cache[e['pdb_id']] = e
    for k, (cid, c) in enumerate(C.items()):
        if cid not in AUD:
            continue
        o = AUD[cid]
        hrows, hmiss = m.audit_side(c['holo'], o['holo_records'], o['holo_assembly_used'], {c['holo_deposited_chain']: c['holo_chain']})
        arows, amiss = m.audit_side(c['apo'], o['apo_records'], o['apo_assembly_used'], {c['apo_deposited_chain']: c['apo_chain']})
        out.append(dict(case_id=cid, apo=c['apo'], holo=c['holo'], ligand=c['ligand_id'],
                        holo_states=m.side_summary(hrows), apo_states=m.side_summary(arows),
                        holo_missing_in_qfit=hmiss, apo_missing_in_qfit=amiss,
                        holo_residues=hrows, apo_residues=arows))
        if k % 25 == 0:
            print(' states', k, '/', len(C), flush=True)
    (HERE/'qfit_state_audit.json').write_text(json.dumps(dict(cases=out), indent=1, default=float) + '\n')
    print('states: ok', len(out))


def refit():
    m = load('rfp', PILOT/'state_audit_v1/code/refit_pose.py')
    m.HERE, m.PILOT = HERE, PILOT
    m.A.load.__globals__['CIFS'] = HERE/'cifs'
    C = json.loads(CANDS.read_text())
    O = {o['case_id']: o for o in json.loads((HERE/'assembly_contact_audit.json').read_text())['cases']}
    out, err = [], []
    for k, c in enumerate(C):
        if c['case_id'] not in O:
            continue
        try:
            out.append(m.case(c, O[c['case_id']]))
        except Exception as e:
            err.append(dict(case_id=c['case_id'], error=repr(e)))
            print('ERROR refit', c['case_id'], repr(e)[:160], flush=True)
        if k % 25 == 0:
            print(' refit', k, '/', len(C), 'errors', len(err), flush=True)
    (HERE/'pose_refit.json').write_text(json.dumps(dict(cases=out, errors=err), indent=1, default=float) + '\n')
    print('refit: ok', len(out), 'errors', len(err))


if __name__ == '__main__':
    stage('assembly_contact_audit.json', 'assembly + crystal contacts', assembly)
    stage('qfit_state_audit.json', 'qFit pocket states', states)
    stage('pose_refit.json', 'whole-pocket apo pose refit', refit)
    print('pipeline stages complete')
