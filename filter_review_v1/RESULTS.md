# 도킹 pilot 필터 탈락 원인 및 완화 민감도 감사

2026-10-28 KST. 사용자 요청: 필터가 지나치게 엄격한지, 어디에서 많이 탈락했는지, 연구 질문을 보존하면서 어떤 조건을 완화할 수 있는지 확인.

## 결론

앞선 '12건·2개 비상동성 성분만 사용 가능'은 당시의 엄격한 **준비 편의 조건** 아래의 결과다. 데이터 전체에 유용한 도킹 사례가 2개뿐이라는 뜻으로 해석하면 안 된다. 특히 다른 사슬·결정 대칭체·다른 HET 분자·금속의 근접을 일괄 탈락시킨 조건이 과도했다.

학습/검증 상동성 차단과 좌표·화학 QC를 그대로 유지하고, 환경 조건을 자동 탈락에서 준비 검토로 바꾸면 67건·22개 비상동성 성분이 남는다. 점유율≥0.5, heavy atoms 6–50, 포켓 backbone RMSD≤1.5 Å를 적용한 확장안은 87건·25개 성분이다. 이는 **후속 구조 준비와 리간드 역할 검토를 할 후보 수**이며, 생물학적으로 유효한 최종 도킹 타깃 수는 아니다.

## 연구 목적과 이번 작업의 범위

질문은 같은 apo backbone에서 single-state 학습 B와 multistate 학습 C의 sidechain ensemble이 리간드 pose 회복에 차이를 주는지다. 구조 정확도 전반의 비교로 넓히지 않는다. 기존 apo–holo 결과에서 strict 대표 2개 모두 apo minor-state 신호가 없었던 점은 모델 효과 부재를 뜻하지 않는다.

고정: 원천 apo–holo 후보, train/validation split, 상동성 조건, 원자·잔기 매핑, 기존 구조 파일, 도킹 결과를 사용하지 않는 선택 원칙.
변경: 일부 설계/준비 필터를 통계적 층화 또는 수동 준비 검토 대상으로 전환. 기존 조기 탈락 때문에 생략된 구조 검사를 실제로 추가 수행했다.

실행: Slurm job 438440, CPU 2개, 2분38초, peak RSS 약704 MiB, COMPLETED. 443건 재검사, 예외 0건. 로컬 RDKit 2025.09.3 화학 검사 완료. 새 학습·모델 추론·도킹은 실행하지 않았다.

## 1. 어디에서 가장 많이 탈락했는가

먼저 상동성 검사 전 단계에서 비보조인자·해상도≤2.5 Å 후보는 3,527 ligand-chain cases, 1,118 source sequence groups였다. train/validation 동일 그룹·직접 중복을 차단하면 1,209건/385그룹, 추가 full-SEQRES 상동성 검사 후 443건/136그룹이다. 이 upstream 감소에는 동일 그룹 차단과 서열 상동성 차단이 함께 들어간다. 구조 필터와 구분해서 읽어야 한다.

443건에 대한 기존 실행 순서의 서로 겹치지 않는 탈락 수:

| 순서 | 기존 조건 | 이 단계 탈락 | 남은 cases | 남은 source groups |
|---|---|---:|---:|---:|
| 시작 | 상동성 검사 통과 | — | 443 | 136 |
| 1 | qFit/cache 서열 매핑 불일치 | 23 | 420 | 133 |
| 2 | heavy atoms 10–50, 최소 점유율≥0.8 | 142 | 278 | 101 |
| 3 | CCD NON-POLYMER 및 리간드 원자 완전성 | 69 | 209 | 85 |
| 4 | 포켓 backbone 매핑/원자/coverage | 21 | 188 | 79 |
| 5 | RMSD≤0.75 Å, peptide 연결, box 조건 | 41 | 147 | 68 |
| 6 | 다른 사슬·결정 대칭체·금속·다른 HET 근접 | 123 | 24 | 16 |
| 7 | deposited 구조 일치/공유결합 등 최종 QC | 12 | 12 | 7 |

절대 탈락 수는 크기·점유율 142건이 가장 크다. 그러나 마지막 병목은 환경 조건으로, 그 단계 입력 147건의 83.7%인 123건을 제거했다. 구조 준비가 복잡하다는 이유와 연구적으로 부적절하다는 이유를 섞었던 것이다.

기존 개별 플래그는 서로 겹친다. 크기 탈락 77건 중 73건은 heavy atom<10, 4건은 >50. 점유율 탈락 78건 중 71건은 0.5–0.8, 7건은 <0.5. 이 수를 더해 142건으로 해석하면 안 된다. 조기 탈락 이후 조건은 검사되지 않았으므로 기존 플래그를 지우기만 해서는 회복 후보 수를 알 수 없다.

## 2. 후속 QC까지 수행한 완화 결과

case는 apo–holo–receptor-chain–ligand-instance 조합이다. 같은 PDB pair도 여러 case를 만든다. 성분 수는 모든 homology-clean query와 후보 간 상동성 edge로 만든 연결 성분 수이며, 생물학적 독립성을 절대 보장하는 표현은 아니다.

| 누적 완화안 | cases | PDB pairs | source groups | 비상동성 성분 | apo MC 관측 cases | apo minor→holo major 관측 cases |
|---|---:|---:|---:|---:|---:|---:|
| 기존 strict | 12 | 8 | 7 | 2 | 6 | 0 |
| 환경 일괄 제외 해제 | 67 | 39 | 33 | 22 | 21 | 3 |
| 위 + 점유율≥0.5 | 75 | 44 | 37 | 24 | 22 | 4 |
| 위 + heavy atoms 6–50 | 83 | 48 | 39 | 25 | 23 | 4 |
| 위 + RMSD≤1.5 Å | 87 | 50 | 41 | 25 | 24 | 4 |
| 위 + CCD type 제한도 완화 | 89 | 51 | 42 | 26 | 24 | 4 |

위 표는 리간드 직접 공유결합/금속 연결, 원자 불완전성, 좌표 불일치 등 현재 hard QC를 계속 적용한 결과다. 마지막 CCD 완화는 탐색용이고 우선 권장하지 않는다. deposited entity가 NonPolymer인지 별도로 검증하며, 당류의 CCD 분류만으로 비공유결합 리간드를 부적절하게 제외할 수 있다는 점을 살펴본 것이다.

단일 조건만 완화하면: 결정 대칭체 25건/9성분, 다른 사슬 27건/8성분, 다른 HET 23건/7성분, 금속 근접 12건/2성분, 점유율 13건/3성분, 크기 14건/3성분, RMSD 12건/2성분. 조건의 중첩 때문에 개별 증가분을 더할 수 없다. RMSD 완화 단독은 0건 회복이며, 환경·점유율·크기 완화 이후에야 4건이 추가된다.

금속 근접은 계속 제외하고 다른 환경 조건만 검토로 돌리면 66건/21성분. 여기에 점유율≥0.5와 크기6–50을 적용하면 **79건/23성분**. RMSD≤1.5까지 적용하면 83건/23성분이다. 따라서 첫 pilot에서 금속 준비를 별도로 미뤄도 후보 확보의 대부분은 유지된다.

또한 환경 자동 제외를 유지한 채 크기·점유율·RMSD만 완화하면 16건/4성분뿐이다. **우선 바꿔야 하는 것은 환경 일괄 제외다.**

## 3. 어떤 기준을 어떻게 바꿀 것인가

### 환경 조건: 제외에서 준비 적합성 검토로

- 다른 사슬이 접촉한다고 실패 사례가 되는 것은 아니다. biological assembly와 실제 결합 부위 구성인지 확인하고, B/C 비교에서 같은 사슬과 환경을 유지한다. packer 입력에서 추가 사슬이 어떻게 처리되는지도 고정한다.
- 결정 대칭체의 근접만으로 제외하지 않는다. 결정 packing이 pose를 직접 안정화하는 경우와 주변 접촉을 구분한다. 용액상 receptor로 재현하려는 질문에서 필수 결정 접촉이 있는 사례는 별도 층으로 두거나 제외한다. 임의로 symmetry mate를 모두 넣어서 통과시키지 않는다.
- 다른 HET 분자는 buffer, cofactor, 다른 리간드, 공유 결합 성분을 구분한다. 보존/제거 원칙을 B/C에서 동일하게 적용한다. holo의 주변 유기 분자를 근거 없이 apo receptor에 옮기지 않는다.
- 6 Å 내 금속은 직접 배위와 다르다. 이번 표에서도 명시적 ligand–metal 연결은 계속 제외했다. 첫 pilot은 금속 근접 case를 별도 대기시켜도 된다. 직접 zinc 배위를 연구하려면 별도 점수화/준비가 필요하다. [Vina 공식 zinc docking 문서](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html).

### 포켓 backbone RMSD: 비교의 통제 조건보다 난이도 층화

B와 C는 각 case에서 **동일한 apo backbone**을 사용하므로, apo–holo RMSD가 0.75 Å보다 커져도 paired 비교의 backbone 고정은 유지된다. 0.75 Å hard cutoff는 이 비교의 필수조건이 아니다. 우선 ≤0.75 Å와 0.75–1.5 Å를 나누고, 더 큰 변화는 stress-test 후보로 남기는 편이 합리적이다. 단 큰 backbone 이동으로 두 방법 모두 실패하면 sidechain 방법의 무용성을 주장할 수 없다. native-holo redocking과 native-apo docking 대조군으로 해석한다.

### 점유율/크기: 범위를 확대하되 reference pose 품질 확인

0.8→0.5는 검토 범위를 넓히기 위한 운영 기준이지 보편적인 품질 기준이 아니다. 전자밀도, 대체 ligand pose, 원자별 점유율과 원자 완전성을 확인한다. 단일 pose가 불분명하면 복수 reference 또는 별도 층으로 처리한다. 높은 점유율만 남기는 것이 관심 있는 이질적 구조를 줄일 가능성은 있으나, 그것이 실제 선택 편향의 원인인지는 아직 검증하지 않았다.

6–9 heavy-atom fragment도 배제할 과학적 필연성은 없다. nicotinamide, catechol 등의 후보가 추가된다. fragment와 큰 ligand는 층화해서 보고하고 symmetry-aware pose RMSD와 reference ambiguity를 확인한다. 단, 작고 유기분자라는 이유만으로 생물학적 ligand라고 분류하지 않는다.

### 그대로 유지하거나 먼저 원인을 해결할 조건

- train/validation 직접 중복 차단, 기존 reciprocal homology 조건(identity≥30%, max coverage≥50%, E≤0.001)을 이번 분석에서는 유지했다.
- 원자 누락·잘못된 residue mapping·끊어진 pocket backbone은 수치 완화로 통과시키지 않는다. mapping 오류23건은 원자료 무효가 아니라 해결 대기 상태다.
- 비공유결합 docking scope에서는 공유결합 및 명시적 금속배위 ligand를 별도 처리한다.
- qFit와 deposited ligand 좌표 차이>0.25 Å는 이번에는 보수적으로 유지했다. 그러나 그 차이가 실제 refinement, altloc 선택, atom-name symmetry, 좌표계 중 무엇 때문인지 해결해야 하며, 그 자체가 영구적인 생물학적 부적격 판정은 아니다. 이 감사에서는 원자 이름 기반·무중첩 RMSD를 썼다. 이를 느슨하게 하려면 protein frame/ligand atom correspondence를 먼저 확인한다.
- homology 30% 차단은 '새 단백질 계열에 대한 일반화'의 조건이다. '기존 계열의 새로운 복합체에서도 유용한가'를 별도 질문으로 두면 seen-family 평가를 추가할 수 있지만, 분리해서 보고하고 exact evaluation leakage는 계속 금지한다. 지금은 환경 필터만 바꿔도 후보가 늘어 이 기준을 먼저 낮출 이유는 작다.

## 4. 후보를 과장하지 않기 위한 추가 확인

확장 87건에는 buffer/additive 가능성이 있는 B3P, BTB, 9JE가 포함된다. 화학 이름은 [B3P CCD](https://www.rcsb.org/ligand/B3P), [BTB CCD](https://www.rcsb.org/ligand/BTB)와 일치한다. 해당 복합체에서 생물학적으로 의미 있는 결합인지 원논문/실험 조건 확인은 미완료다. 이 3개 component ID의 ligand 사례를 예비 제외하면 82건/23성분이다. 이는 완전한 ligand-role 분류를 마친 수치가 아니다. 한 성분에 다른 ligand가 있으면 그 성분은 유지되므로 25−3으로 계산하지 않는다.

apo minor→holo major 관측은 확장안 4건: 3hre→3huo PNN:A:302, 6qo5→6qo9 BTB:B:401, 8ftp→8tav ZKR:A:303, 5hih→5hol APR:A:201. 이 중 BTB는 additive 검토 대상이다. 이 관측은 원자료의 holo-defined 4.5 Å contact/χ1 정의에 한정되며 전체 pocket χ1–χ4의 다중 상태나 동역학 부재/존재를 확정하지 않는다.

발견을 보여줄 mechanism subset으로 관측 양성 case를 별도 제시할 수 있지만, 이것만 골라 전체 docking 성능 개선을 주장하면 선택 편향이다. 일반 성능 cohort는 docking 결과와 minor-state 양성 여부에 의존하지 않는 기준으로 정하고, 관측 양성 subset 분석은 따로 명시한다.

## 5. 권장 다음 단계와 계산 비용

우선 79건/23성분(환경 검토 허용, 금속 근접 제외, 점유율≥0.5, 크기6–50, RMSD≤0.75 Å)을 대상으로 ligand 역할과 assembly/결정 접촉을 검토한다. 필요하면 RMSD 0.75–1.5 Å 4건을 별도 난이도 층으로 추가한다. 목표24개에 맞추기 위해 임계값을 사후 조정하지 않는다. 12–20개 수준의 독립 성분이 확정돼도 pilot의 실행 가능성과 paired 효과 범위를 확인하는 데 사용할 수 있다. 검출력 보장은 하지 않는다.

이 단계는 CPU 구조 준비와 문헌 검토이며 새 GPU 학습이 필요 없다. 이후 실제 docking은 기존 5조건×8 receptor/search slots×3 반복이면 case당120회, 12–20 representatives면 1,440–2,400회다. 실제 시간은 소규모 throughput 측정 전 확정하지 않는다. 오늘은 후보 감사만 실행했다.

## 재현성과 검증

- `audit_sensitivity.py`: 원 strict geometry 코드에서 설계 제한을 플래그로 유지하되, 그 제한만 걸린 사례도 후속 QC를 실행. hard QC 실패는 계속 중단한다.
- `sensitivity_candidates.json`, `sensitivity_summary.json`: 원시 검사 결과. summary의 passed는 원 strict 조건 통과 수를 의미한다.
- `summarize_sensitivity.py`: RDKit 화학 검사, 조건별 재집계, 원 strict12건 case_id 동일성 assert.
- `review_summary.json`: 단계별 탈락·완화안별 수·입력 SHA256.
- `review_records.json`, `review_candidates.tsv`: 후보별 플래그, 적용 가능한 검토 시나리오. 최종 pilot manifest가 아니다.
- 모든443 case_id의 일치/유일성 검사, 기존 strict12개 재현, 예외0개, rigid-transform 합성 검사 통과.
- 연결 성분 계산은 hit edge 없는 query도 포함한다. 전체443건은78성분이다. 이전 중간 집계72는 isolated query6개를 빠뜨린 수치였으며 여기서 정정한다. 최종 strict2성분 결과는 바뀌지 않는다.
- remote source root: `/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/filter_review_v1`.
- local persistent root: `/data/work/sidechain-multistate/experiments/docking_pilot_20261028/filter_review_v1`.

관측 사실: 완화 후 검사 가능한 후보와 비상동성 성분이 늘었다. 아직 알 수 없는 것: 이들 중 최종 receptor 준비가 가능한 수, C가 B보다 도킹에 유용한지, 효과 크기와 일반화 범위. 이 감사만으로 모델의 우위를 주장하지 않는다.
