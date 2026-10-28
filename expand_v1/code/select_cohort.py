"""확대 cohort 고정: 기하 QC 결과 + 자동 리간드 역할 -> candidates.json (79건 파일과 같은 형식).

하류 코드(assembly 감사, 상태 감사, 수용체 구축)가 그대로 읽을 수 있도록 필드 이름을 맞춘다.
homology_component는 아직 없다(확대 집합의 MMseqs all-vs-all 미실행). 도킹에는 필요 없고 분석에서만 쓰므로
`homology_component=None`으로 두고, 분석 전에 채운다.
"""
import json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE/'code'))
from ligand_roles_auto import role

SOFT = {'heavy_atom_count_outside_10_50', 'ligand_occupancy_below_0_8', 'ccd_not_nonpolymer',
        'pocket_backbone_rmsd_above_0_75', 'apo_multichain_pocket', 'holo_multichain_pocket',
        'apo_crystal_mate_pocket', 'holo_crystal_mate_pocket', 'apo_nearby_metal', 'holo_nearby_metal',
        'apo_nearby_other_heterocomponent', 'holo_nearby_other_heterocomponent'}
METAL = {'apo_nearby_metal', 'holo_nearby_metal'}


def keep(c):
    """79건 입력과 같은 조건: 금속 근접 제외, 점유율 >= 0.5, heavy 6-50, pocket RMSD <= 0.75."""
    r = set(c['geometry_reasons'])
    return (not (r - SOFT) and not (METAL & r)
            and c.get('ligand_min_occupancy', 0) >= 0.5 - 1e-6
            and 6 <= c.get('ligand_heavy_atoms', 0) <= 50
            and c.get('pocket_backbone_rmsd', 9) <= 0.75)


def main():
    R = json.loads((HERE/'sensitivity_candidates.json').read_text())
    geom = [c for c in R if keep(c)]
    buffers = [c for c in geom if role(c['ligand_id']) == 'exclude']
    sel = [dict(c, homology_component=None) for c in geom if role(c['ligand_id']) == 'accept']
    (HERE/'candidates.json').write_text(json.dumps(sel, indent=1) + '\n')
    summary = dict(audited=len(R), passed_geometry=len(geom), excluded_as_buffer=len(buffers),
                   cohort=len(sel), source_groups=len({c['group'] for c in sel}),
                   structures=len({x for c in sel for x in (c['apo'], c['holo'])}),
                   buffer_ccds=sorted({c['ligand_id'].split(':')[0] for c in buffers}),
                   note='homology_component is not assigned yet; it is only needed for the analysis bootstrap unit.')
    (HERE/'cohort_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
