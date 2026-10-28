"""확대 cohort의 수용체 구축 + Vina 입력 준비. 기존 코드를 모듈로 불러와 경로만 바꾼다(계산 코드는 그대로).

--shard i/n 으로 사례를 나눠 array로 돌린다. 이미 만들어진 사례는 건너뛴다.
"""
import argparse, importlib.util, json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
PILOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--shards', type=int, default=1)
    a = ap.parse_args()

    import csv
    rows = [r for r in csv.DictReader(open(HERE/'cohort_final.tsv'), delimiter='\t')
            if r['status'] in ('accept', 'backbone_change')]
    mine = [r for k, r in enumerate(rows) if k % a.shards == a.shard]
    print(f'shard {a.shard}/{a.shards}: {len(mine)} of {len(rows)} cases', flush=True)

    C = {c['case_id']: c for c in json.loads((HERE/'candidates.json').read_text())}
    A = {o['case_id']: o for o in json.loads((HERE/'assembly_contact_audit.json').read_text())['cases']}
    F = {o['case_id']: o for o in json.loads((HERE/'pose_refit.json').read_text())['cases']}

    B = load('brc', PILOT/'design2x2_v1/code/build_receptors.py')
    B.HERE, B.PILOT = HERE, PILOT
    B.A.load.__globals__['CIFS'] = HERE/'cifs'
    P = load('pin', PILOT/'docking_v1/code/prepare_inputs.py')
    # prepare_case는 src를 PILOT/'design2x2_v1/receptors'로 계산한다. PILOT을 여기로 돌리고
    # design2x2_v1/receptors를 이 폴더의 receptors로 링크해 같은 코드를 그대로 쓴다.
    P.HERE, P.PILOT = HERE, HERE

    # 확대 집합에는 Meeko가 template를 못 맞추는 잔기가 있는 수용체가 나온다(예: 접촉이 너무 가까워
    # 탄소 원자가가 5가 되는 LYS). 142건을 손으로 고칠 수 없으므로, 엄격 호출이 실패하면 해당 잔기를
    # 제거하고 무엇을 제거했는지 기록한다. prepare_case가 arm 사이에 제거 집합이 같은지 이미 단언하므로
    # 비교 가능성은 유지된다. 제거된 잔기가 pocket이면 그 사례는 하류에서 버린다.
    _strict = P.receptor_pdbqt

    def receptor_pdbqt(pdb_h, stem, center, size):
        try:
            return _strict(pdb_h, stem, center, size)
        except RuntimeError as e:
            meeko_in = Path(str(stem) + '.meeko.pdb')
            r = P.run([P.ENV/'mk_prepare_receptor.py', '--read_pdb', meeko_in, '-o', stem, '-p', '-v',
                       '--box_center', *center, '--box_size', *size, '-a'])
            dropped = [l.strip() for l in (r.stdout + r.stderr).splitlines()
                       if 'delet' in l.lower() or 'ignor' in l.lower() or 'removed' in l.lower()]
            return Path(f'{stem}.pdbqt'), dropped + [f'allow_bad_res after: {str(e)[:120]}']

    P.receptor_pdbqt = receptor_pdbqt

    import numpy as np
    built, prepared, failed = [], [], []
    for k, r in enumerate(mine):
        cid = r['case_id']
        try:
            if not (HERE/'receptors'/cid/'apo_single_0.pdb').exists():
                c, o, f = C[cid], A[cid], F[cid]
                rng = np.random.default_rng(B.SEED)
                outdir = HERE/'receptors'/cid; outdir.mkdir(parents=True, exist_ok=True)
                hs, aps = B.A.load(c['holo']), B.A.load(c['apo'])
                comp, _, lnum = c['ligand_id'].split(':')
                lig = [x for ch in hs[0] if ch.name == c['ligand_deposited_chain'] for x in ch
                       if x.name == comp and x.seqid.num == int(lnum)][0]
                la, _ = B.A.select_atoms(lig)
                L = B.A.xyz(la); X = np.array(f['apo_pose'])
                atoms = [(comp, x.name, x.element.name) for x in la]
                B.write_ligand(atoms, L, outdir/'ligand_ref_holo.pdb')
                B.write_ligand(atoms, X, outdir/'ligand_ref_apo.pdb')
                anchor = [s for s, d in B.A.subchains(aps).items()
                          if d['chain'] == c['apo_deposited_chain'] and d['kind'] == 'Polymer'][0]
                arec, aasm = B.A.environment(aps, X, anchor, 'apo')
                o_apo = dict(apo_records=arec,
                             apo_assembly_used=sorted(aasm, key=lambda n: (not aasm[n]['anchor_present'], n))[0])
                holo = B.build_side(c['holo'], B.deposited_copies(o, 'holo', hs), L, outdir, 'holo', rng,
                                    {c['holo_deposited_chain']: c['holo_chain']})
                apo = B.build_side(c['apo'], B.deposited_copies(o_apo, 'apo', aps), X, outdir, 'apo', rng,
                                   {c['apo_deposited_chain']: c['apo_chain']})
                (outdir/'build.json').write_text(json.dumps(dict(case_id=cid, holo=holo, apo=apo), indent=1, default=str) + '\n')
                built.append(cid)
            if not (HERE/'inputs'/cid/'ligand_dock.pdbqt').exists():
                case = json.loads((HERE/'receptors'/cid/'build.json').read_text())
                P.prepare_case(case, C)
                prepared.append(cid)
        except Exception as e:
            failed.append(dict(case_id=cid, error=repr(e)[:300]))
            print('FAILED', cid, repr(e)[:200], flush=True)
        if k % 10 == 0:
            print(f' {k}/{len(mine)} built {len(built)} prepped {len(prepared)} failed {len(failed)}', flush=True)
    (HERE/f'build_prep_shard{a.shard}.json').write_text(json.dumps(
        dict(shard=a.shard, cases=len(mine), built=built, prepared=prepared, failed=failed), indent=1) + '\n')
    print('shard done: built', len(built), 'prepared', len(prepared), 'failed', len(failed))


if __name__ == '__main__':
    main()
