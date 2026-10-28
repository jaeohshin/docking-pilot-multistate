# Apo–holo docking pilot 적격성 검사 — 2026-10-28

## 결론

추가 상동성, 리간드 화학 분류, 포켓 backbone RMSD 및 단순 포켓 조건 검사를 완료했다.
보수적인 자동 검사에서 **12개 ligand–chain 사례 / 8개 apo–holo PDB 쌍 / 7개 원본 서열 그룹**이 통과했다.
이들끼리 상동성을 묶으면 **2개 독립 계열**이며, 계열당 대표 한 개를 선정했다.
따라서 최초의 24-target pilot은 현재 후보 풀과 이 조건으로 충족되지 않는다.
선정된 두 사례는 **도킹 절차 점검용**이다. 논문 발견의 일반적 효용을 검증하는 평가셋으로는 부족하다.

학습, 모델 추론, docking은 실행하지 않았다. 기존 checkpoint와 실험 결과는 변경하지 않았다.
CPU 검사와 제한된 RCSB 구조/CCD 조회만 실행했고, 분석 결과는 별도 experiments 경로에 저장했다.

## 연구 목적과 판단 범위

검증할 가설은 같은 apo backbone과 도킹 예산에서 occupancy 학습 C가 single-state B보다 ligand pose recovery를 개선하는가이다.
기존 C−B minor-well AUC 향상은 이 가설의 출발점이며, docking 개선을 이미 보여준 결과는 아니다.
이번 작업은 데이터 적격성 검사다. 작은 최종 표본은 모델 성능의 음성 결과가 아니다.
qFit의 apo χ1 minor state가 관측되지 않는다고 실제 대체 상태가 없다고 단정하지 않는다.

## 단계별 결과

| 단계 | ligand–chain 사례 | apo–holo PDB 쌍 | 원 분석의 서열 그룹 |
|---|---:|---:|---:|
| 비보조인자, 양쪽 해상도 ≤2.5 Å | 3,527 | 1,691 | 1,118 |
| train/validation 동일 그룹 및 PDB 직접 중복 제외 | 1,209 | 598 | 385 |
| 전체 서열 검증 및 양방향 추가 상동성 제외 | 443 | 197 | 136 |
| 화학·포켓·결정 접촉 조건 통과 | 12 | 8 | 7 |
| 통과 사례끼리 상동성 군집화 후 대표 선정 | 2 | 2 | 2개 비상동 계열 |

원 분석은 동일서열 그룹당 최고 해상도 apo 1개와 최대 holo 3개를 선택한 자료다.
따라서 이는 모든 qFit/PDB apo–holo 조합의 최대 가용량이 아니다.
동일 ligand가 여러 chain과 접촉하는 사례는 중간 집계에서 별도 행이다. 이를 독립 단백질 수로 해석하지 않는다.

초기 43개 사례는 full SEQRES를 이름으로 찾지 못했다. 예를 들어 qFit chain A가 mmCIF auth chain AA에 해당했다.
공식 mmCIF의 관측 서열을 qFit cache와 정렬해 대응을 확인하고, polymer entity의 전체 서열로 복원했다.
전체 1,209개 사례의 paired full SEQRES를 확인했다. 392개 고유 receptor 서열 중 256개가 추가 상동성 검사에서 제외됐고 136개가 남았다.
새로 복원된 사례 중 2개가 상동성 검사를 더 통과했지만, 최종 구조 조건은 통과하지 못했다.

## 선정된 대표 사례

| 용도 | Apo chain | Holo chain | ligand instance | 포켓 backbone RMSD | 포켓 coverage | ligand heavy atoms | rotatable bonds |
|---|---|---|---|---:|---:|---:|---:|
| 기술 점검 | 7k4t A | 7k4z A | 6NT:A:402 | 0.267 Å | 100% | 12 | 1 |
| 기술 점검 | 8tn1 A | 8tnb A | IQM:A:201 | 0.253 Å | 100% | 22 | 4 |

7k4z의 deposited entity 이름은 Endo-1,4-beta-xylanase이며, 다른 통과 단백질들은 Kemp Eliminase HG3/변이체로 기재돼 있다.
이 계열 사이 직접 서열 일치도는 약 94.9–99.3%이므로 여러 독립 표적으로 세지 않았다.
8tnb/8tnc의 deposited entity는 De novo designed protein이다.
따라서 일반적인 자연 단백질 약물 표적군을 대표하지 않는다.

대표 선정은 계열당 양쪽 해상도의 최댓값, 합, PDB/chain/ligand ID 순으로 결정했다.
포켓 RMSD의 세부 값, minor-state 유무, docking 결과로 대표를 유리하게 고르지 않았다.

12개 전체 통과 사례 중 6개에는 원 분석의 4.5 Å holo-contact 잔기에 해당하는 apo χ1 다중상태 잔기가 있다.
그러나 모든 통과 사례에서 holo major가 apo minor에 해당하는 잔기는 0개다.
선정된 두 대표 사례에서는 해당 contact 잔기의 apo χ1 다중상태도 0개다.
이는 단일 χ1 기준/점유율 0.1/접촉 범위에 관한 관측이며, 모든 χ/포켓 전체의 동역학 부재를 뜻하지 않는다.

동일 계열의 보조 사례 5rg6 A → 5rgc A, 6NT:A:401에는 apo χ1 다중상태 접촉 잔기 2개가 있다.
`same_family_diagnostic.json`에 따로 남겼다. 샘플러의 대체 상태 생성 확인용으로 사용할 수 있으나 세 번째 독립 표적으로 세면 안 된다.
이 보조 사례도 observed apo-minor → holo-major 사례는 아니다.

## 검사 방법

### 상동성

- 실제 v2 run의 frozen train/validation full SEQRES를 사용했다. manifest와 FASTA를 snapshot SHA256과 대조했다.
- MMseqs2 reciprocal search: identity ≥30%, query 또는 target coverage ≥50%, E ≤1e−3, sensitivity 7.5.
- literal full-sequence equality도 별도로 제외했다. 후보끼리 같은 조건의 검색 후 connected component로 묶었다.
- 미세한 서열 차이를 서로 다른 독립 표적으로 세지 않았다. 관측된 검색 hit 부재는 모든 진화적 상동성 부재의 증명은 아니다.
- 사슬 alias 복원은 같은 PDB 내 관측 서열 정렬의 동일 잔기 수가 두 관측 서열 길이의 최댓값 대비 ≥98%인 경우로 제한했다.
  최상위 대응들이 서로 다른 전체 서열을 주면 보류하는 규칙을 사용했다. 관측 서열을 full SEQRES 대신 사용하지 않았다.

### 리간드

- 원 목록의 cofactor 제거 후 CCD type NON-POLYMER, heavy atoms 10–50, 완전한 CCD atom-name 대응, 주요 conformer 모든 원자의 occupancy ≥0.8.
- RDKit 2025.09.3으로 CCD SMILES를 읽어 heavy-atom 수 일치, 단일 연결 성분, 지원 유기 원소, rotatable bonds ≤12를 검사했다.
- CCD polymer/saccharide 형식, covalent/metal 연결, protein과의 1.8 Å 미만 접촉을 보류/제외했다.
- 마지막 12개는 RDKit 화학 검사에서도 모두 통과했다. 이것이 binding activity나 protonation 상태의 정답을 보증하지는 않는다.

### 좌표와 포켓

- qFit 원본 구조와 replicate cache의 chain 관측 서열 일치를 확인했다.
- residue별 최고 평균 occupancy altloc과 공통 atoms를 선택하되 같은 residue에서 다른 altloc atoms를 임의 혼합하지 않았다.
- holo ligand 기준 4.5 Å contact shell, 6.5 Å pocket shell. backbone N, CA, C, O를 apo/holo 동일 잔기로 매핑했다.
- contact backbone 결손 없음, shell coverage ≥95%, 공통 backbone 원자의 local Kabsch RMSD ≤0.75 Å.
- peptide junction 0.9–2.0 Å, box edge ≤30 Å. Kabsch 구현은 독립 rigid transform 합성 검사를 통과했다.
- 다른 polymer chain 또는 crystal image가 ligand 4.5 Å 안에 있는 사례를 첫 단일 사슬 pilot에서 보류했다.
- ligand 6 Å 안의 금속, 그리고 물 이외 6-heavy-atom 이상 다른 heterocomponent가 있는 사례도 보류했다.
  후자는 해당 성분의 필수성이 입증됐다는 뜻이 아니라, 첫 pilot에서 복합적인 preparation을 피하는 보수적 규칙이다.
- deposited mmCIF에서 ligand nonpolymer 분류와 covalent/metal links를 확인하고,
  qFit/deposited ligand 좌표 차이가 ≤0.25 Å인지 검사했다.
- 통과 구조의 crystal cell과 space group 존재를 별도 확인했다.

## 탈락과 보류의 해석

수정 아미노산/당, 작은 fragment, 낮은 ligand occupancy, 불완전한 ligand atom 목록,
금속·다중 사슬·crystal mate·다른 heterocomponent가 관여하는 포켓이 상당수였다.
backbone RMSD >0.75 Å는 34개 사례에서 기록됐다.
qFit와 cache의 관측 서열이 맞지 않는 23개 사례는 안전하게 보류했다. 생물학적으로 부적격하다고 판정한 것은 아니다.

`exclusions.tsv`에 사례별 사유가 있다. 한 사례에 여러 사유가 있을 수 있고,
앞 단계에서 탈락하면 불필요한 후속 계산은 하지 않았으므로 이유별 숫자를 합산하거나 완전한 독립 원인분해로 읽으면 안 된다.
특히 단일 사슬/금속 없는 포켓 등의 제한은 초기 pilot의 구현 단순화를 위한 선택이다.
이 제한으로 탈락한 구조가 docking 자체에 쓸 수 없다는 의미가 아니다.

## 실험 계획 갱신

현재 선정 2개는 B/C/B-temperature/apo/holo × receptor 또는 반복 search slot 8개 × docking seed 3개로 총 240회 도킹인 기술 점검에 사용할 수 있다.
추론과 도킹은 아직 실행하지 않았다. 시간은 측정하지 않았으므로 실제 ETA를 제시하지 않는다.
이 두 구조로 모델 우위, 일반적 pose-recovery 향상, virtual-screening 향상 또는 효과 부재를 결론내리지 않는다.
독립 개발/검증 분할을 만들 만큼 표본이 없으므로 별도 development manifest는 비어 있다.

논문의 발견이 실제 도킹에 도움이 되는지 검증하려면 다음 설계 결정이 필요하다:
1. 다른 apo/holo 또는 holo/holo 교차 도킹 자료로 학습 비상동 표적을 확대하거나,
2. 적절한 assembly·보조인자 준비를 지원해 보류한 복합 포켓을 다시 검토한다.
어느 경우에도 기준 변경을 사전에 기록하고 기존 strict 결과와 구분해야 한다. 이번 작업에서 자동 확장하지 않았다.

## 파일과 재현

- `pilot_manifest.tsv`: 최종 대표 2개.
- `eligible_cases.tsv`: 전체 통과 12개; 서로 상동인 사례 포함.
- `exclusions.tsv`: 추가 상동성 통과 후 443개 중 탈락/보류 사유.
- `selected_geometry.json`: 선정 사례의 residue map, holo→apo rotation/translation, box.
- `all_records.json`, `summary.json`: 전체 검사 결과와 집계.
- `same_family_diagnostic.json`: 동일 계열의 보조 다중상태 사례; 독립 표적 아님.

클러스터: `/gpfs/jaeohshin/work/sidechain-multistate/experiments/docking_pilot_20261028/`.
`run_v1`의 초깃값 검사와 `run_v3`의 최종 chain-alias 복원/상동성 결과를 보존했다.
`run_v2`는 alias 처리 전 full SEQRES 복원 시도이며, 최종 후보 선정에는 run_v3를 사용했다.
run_v1의 첫 geometry job은 numpy int64 JSON 직렬화 오류로 실패했고, int 변환을 명시한 geometry_v2로 수정했다.
최종 geometry_v3는 입력/사슬 대응이 같은 완료 사례만 재사용하고, 새로 통과한 사례를 추가 검사했다. 최종 처리 예외는 0개다.
원격에서 사용한 라이브러리는 Gemmi 0.7.5, Biopython 및 numpy; RDKit 화학 검사는 데스크톱에서 실행했다.
Slurm CPU jobs: 438408, 438424(실패), 438425, 438426, 438430, 438431. GPU는 요청하지 않았다.

초기의 `PLAN.txt`, `AUDIT.json`, `REMOTE_PRESCREEN_20261028.txt`는 접근 불가/중간 단계 기록이다.
현재 판단은 이 문서와 `selection_v3/summary.json`이 우선한다.

공식 구조/CCD: https://files.rcsb.org/download/ 와 https://files.rcsb.org/ligands/download/.
결정 이웃 검사 API: https://gemmi.readthedocs.io/en/latest/analysis.html.

## 실행 비용과 리비전

원격 저장소 HEAD: `0f26be84d9b9fab59ccc7c85551a677afdb3bf63`.
Slurm 작업 6개의 실제 CPU 사용 시간 합은 약 14분 17초이며, 작업별 경과 시간 합은 6분 7초다.
이는 데이터 확인·SSH·분석 코드 작성 시간을 포함한 전체 작업 시간이 아니다.
최종 작업 상태는 모두 완료이며, 보존한 첫 geometry 직렬화 오류 작업만 FAILED다.
