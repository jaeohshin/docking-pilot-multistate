# 2×2 docking 설계: 대표와 수용체 구축 (design2x2_v1)

2026-10-28 KST. 설계(사용자 결정): holo/apo backbone × 단일/다중 sidechain 상태. qFit-at-scale ground truth를 쓰고 packer는 쓰지 않는다.
도킹·모델 추론·학습은 실행하지 않았다. commit하지 않았다. 입력은 `final_selection_v1`, `state_audit_v1`이며 읽기만 했다.

## 확정한 규칙

- 단일 상태: 같은 qFit 모델의 major 상태다. pocket 잔기는 점유율 합이 가장 큰 χ1/χ2 well(원자료의 major 정의)이고, 그 well에서 점유율이 가장 높은 altloc으로 나타낸다. pocket 밖 잔기는 평균 점유율이 가장 높은 altloc이다.
  - pocket 잔기 416개 중 11개(2.6%)는 점유율 최고 altloc이 다른 well에 있다. 일부는 0.5/0.5 동률이고, 동률은 altloc 순서로 고정했다.
  - 이 11개 때문에 단일 수용체가 '점유율 최고 altloc 모델'과 일부 잔기에서 다르다(최대 4.6 Å).
- 다중 상태: major backbone은 고정하고 sidechain만 바꾼다.
  - 각 상태(well 대표 altloc)의 N/CA/C를 major N/CA/C에 강체 정합해 sidechain을 옮겼다(local-frame graft). refinement된 결합 길이·각은 그대로이고, χ 변화는 최대 3.9°다.
  - 처음에는 'χ를 표준 기하로 재구축'하는 방식을 제안했지만, 기하를 덜 바꾸는 이 방식으로 바꿨다.
  - 단일 arm도 같은 함수로 만들었다(major 상태는 자기 backbone이라 원좌표와 같다).
- receptor 예산 K = 8.
  - 조합 수 ≤8이면 전부 쓴다. 넘으면 모든 잔기가 major인 조합을 넣고, 잔기별 well 점유율로 독립 추출한다(seed 20261028).
  - 다른 잔기와 heavy-atom 거리 <2.2 Å인 조합은 버렸다(인접 잔기 backbone 제외). qFit altloc의 잔기 간 결합 관계는 알 수 없어 독립 추출로 생기는 충돌을 이렇게 걸렀다.
  - multi_0은 항상 모든 잔기가 major인 조합이라 단일 수용체와 같다.
- 수용체 구성: 리간드 6 Å 이내의 assembly 내 polymer copy 전체 사슬이다. 대칭 copy는 deposited 연산자를 qFit 좌표에 적용했다. HETATM과 물은 모두 뺐다.
- 탐색량 맞춤(도킹 단계 규칙, 결과 전 고정): 각 arm의 도킹 호출은 8회로 같다.
  - 단일 arm: 수용체 1개에 seed 8개
  - 다중 arm: 수용체 n개에 8회를 고르게 나눈다(round-robin, n=8이면 각 1회)
  - 같은 exhaustiveness와 box를 쓴다(box는 reference 리간드 기준이며 holo/apo가 각자 자기 좌표계).
  - 주 지표: 8회 중 docking score 1위 pose의 symmetry-aware RMSD. 보조 지표: 최저 RMSD(oracle).

## 대표 (`representatives.tsv`, `code/select_reps.py`)

- component 안 순서: holo와 apo pocket 모두에 다중 상태 잔기 ≥1개(처리 변수), 그다음 final_selection_v1 순서(pose 명확성, 환경 등급[apo는 전체 pocket 재정합 기준], 해상도, case_id).
- 결과·minor→major 신호는 쓰지 않았다.

| 층 | component | 대표 | 처리(H/A 다중 잔기) | 비고 |
|---|---|---|---|---|
| main | amidase | 5NXZ→5NY7 NCA:A:302 | 0/0 | 차이 0 대조 |
| main | class A β-lactamase | 3HRE→3HUO PNN:B:301 | 2/3 | 대표 변경(4PM7 → 3HUO) |
| main | Ubc9 | 5F6E→5F6V 5VL:A:201 | 2/1 | |
| main | Kemp eliminase | 8USG→8USH 6NT:A:401 | 4/3 | |
| main | CaM–SK2 | 4J9Y→4G27 PHU:R:1006 | 5/5 | |
| main | T. brucei macrodomain | 5FSV→5FSY AR6:A:1262 | 6/7 | |
| main | nucleoside phosphorylase | 3UAV→3UAW ADN:A:236 | 2/3 | |
| main | HMGR | 6EEU→6EEV MEV:A:501 | 2/0 | holo만 다중 |
| main | CoV macrodomain | 5HIH→5HOL APR:A:201 | 7/12 | minor→major 양성(선정에 쓰지 않음; 같은 component의 8SH8은 0/0) |
| main | PARPi binder | 8TN1→8TNB IQM:A:201 | 1/2 | 대표 변경(IQM:C → IQM:A) |
| main | GSTO3S | 6F43→6F66 CUT:A:300 | 3/3 | 대표 변경(CUT:B → CUT:A) |
| main | IRED | 5A9R→5A9T YLG:A:1291 | 3/3 | |
| main | EvaD | 1OFN→1OI6 TMP:A:1204 | 1/4 | |
| backbone_change | OMPDC | 3G1Y→3G22 U5P:A:229 | 1/4 | 0.90 Å |
| backbone_change | UPP | 3DPS→2PGA ANU:B:7016 | 1/0 | 4.49 Å, apo 자리 대응 불확실 |
| backbone_change | MTAN | 6AYM→6AYQ TDI:B:301 | 1/1 | 2.19 Å |
| backbone_change | SARS-CoV-2 Mac1 | 7TX4→7TX5 APR:A:201 | 4/3 | 0.76 Å |

## 구축 결과 (`receptors/<case>/`, `build_manifest.json`, `verify_receptors.json`)

- 17건 × {holo, apo} × {single 1개, multi 1–8개}, `ligand_ref_holo.pdb`/`ligand_ref_apo.pdb`(reference pose, heavy atom).
- 다중 수용체 수(holo/apo): 5NY7 1/1, 3HUO 4/8, 5F6V 6/2, 8USH 5/4, 4G27 3/8, 5FSY 8/8, 3UAW 4/8, 6EEV 4/1, 5HOL 7/4, 8TNB 2/8, 6F66 8/8, 5A9T 8/6, 1OI6 2/8, 3G22 1/8, 2PGA 2/1, 6AYQ 2/2, 7TX5 8/8.
  - 충돌로 버린 조합: 0–5개(5HOL apo 4, 4G27 holo 5).
  - 3G22 holo는 major 외 상태가 모두 충돌로 빠져 holo 다중 = 단일이다.
- 검증:
  - 모든 다중 수용체의 backbone이 단일과 같다(최대 0.000 Å). 바뀐 원자는 pocket sidechain뿐이다.
  - reference 리간드와 qFit holo 리간드의 차이는 ≤0.25 Å다. qFit apo와 deposited apo 수용 사슬 CA의 차이(무정합)는 ≤0.40 Å다. 좌표계는 일치한다.
  - 단일 수용체 자체의 충돌: 5A9T apo(Asn241, Tyr179), 1OI6 holo(B:Phe19). 원 qFit 상태의 조합에서 나온 것으로, 기록만 했다.

## 남은 일(도킹 전)

- 리간드 준비: bond order, protonation. 클러스터 env에 RDKit이 없어 데스크탑 RDKit이나 별도 env가 필요하다. 수용체 protonation과 atom typing도 필요하다.
- 도킹 엔진과 exhaustiveness, box 크기 확정. 소규모 throughput 측정은 승인 후에 한다.
- 대표 17건 × 4 arm × 8회 = 544 docking 호출이다(deposited holo sanity check 별도).
