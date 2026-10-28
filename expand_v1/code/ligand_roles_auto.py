"""확대 집합(C)의 리간드 역할 자동 판정. 79건 수작업 판정에 대해 검증한 규칙이다.

79건은 사람이 문헌·entry 제목·좌표를 보고 판정했다. 수천 건에는 쓸 수 없으므로 규칙으로 대체하되,
**같은 79건에 적용해 수작업과 비교해 성능을 재고 그 수치를 함께 보고한다.**

규칙: 결정화 완충제·동결보호제·세제에 해당하는 CCD만 제외한다. 나머지 판정은 하류의 기존 자동 필터
(결정 접촉 의존, 보조인자 부재, partner 사슬 부재, RSCC)에 맡긴다.

검토했다가 버린 규칙: "한 entry에 같은 CCD가 여러 개면 접촉이 가장 많은 것만 남긴다".
4량체의 네 사슬에 같은 자리로 붙은 경우(2WYW TCL A/B/C/D, 8TN6 RPB A/B/C)를 중복으로 오인해
수작업 accept 18건을 잘못 버렸다. 자리 동일성을 접촉 잔기로 식별해도 수용 사슬이 다른 경우를 가르지 못했다.

79건 검증 (`validate.py`):
- 거짓 제외 **0건**. 수작업 accept 62건을 하나도 잘못 버리지 않는다.
- 놓침 8건. 전부 "부차 자리의 여분 분자" 또는 "문헌 미확인"이며, 완충제는 4/4 모두 잡는다.
- **대표 선정은 13개 중 12개가 동일**하다. 유일한 차이는 5NY7의 다른 니코틴아미드 분자(NCA:A:302 대 303)다.

놓친 사례가 섞이는 것의 영향: 주 대비는 사례 내 쌍 비교(같은 사례의 single 대 multi)다. 의미 없는 리간드는
양쪽 arm에 똑같이 들어가므로 차이가 0에 가까워진다. 즉 **효과를 희석시키지 편향시키지 않는다.**
"""

# 원천 ADDITIVES(scripts/177)에 없던 결정화 완충제·동결보호제·세제. B3P/BTB/FLC는 79건 수작업에서 실제로 걸렸다.
BUFFERS = {
    'B3P', 'BTB', 'FLC', 'NHE', 'CXS', 'TAR', 'TCE',
    'PGE', 'PG4', 'P6G', '1PE', '2PE', '7PE', '12P', '15P', 'XPE', 'PE4', 'PE8', 'P33',
    'MPD', 'MRD', 'TRS', 'EPE', 'MES', 'DMS', 'IMD', 'CIT', 'TLA', 'BCT',
    'SIN', 'MLI', 'MLA', 'HEZ', 'PEO', 'BU3', 'PGO', 'PGR', 'DIO',
    'SPD', 'SPM', 'MOH', 'IPA', 'ETA', 'ACN', 'DTT', 'BME', 'OLC', 'LDA', 'C8E',
}


def role(ligand_id):
    """'CCD:chain:num' -> 'accept' | 'exclude'."""
    return 'exclude' if ligand_id.split(':')[0] in BUFFERS else 'accept'
