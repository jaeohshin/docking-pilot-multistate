# qFit pocket 다중 상태 감사 — 2×2 docking 설계 (state_audit_v1)

2026-10-28 KST. 설계(사용자 결정): holo/apo backbone × 단일/다중 sidechain 상태, packer 없이 qFit-at-scale ground truth를 사용한다(packer B/C 비교는 이후).
단일 상태 = 같은 qFit 모델의 major conformer만 남긴 구조(제안, 미확정). 도킹·모델 추론·학습은 실행하지 않았다. commit하지 않았다.
입력: `final_review_v1/candidates79.json`, `final_selection_v1/{assembly_contact_audit.json, case_review_all79.tsv}`, qFit zip. 이전 결과는 수정하지 않았다.

## 1. 전체 pocket 기준 apo–holo backbone 재정합 (`code/refit_pose.py`)

이전 이식(final_selection_v1)과 입력 필터의 pocket RMSD는 수용 사슬(qFit `holo_chain`)의 매핑 잔기만 썼다. interface 부위에서는 접촉이 적은 쪽 사슬의 4–7잔기에 그쳤다.
이번에는 리간드 6 Å 이내의 모든 assembly 내 holo 사슬 copy를 apo의 같은 단백질 copy와 짝지어 맞췄다. 사슬마다 관측 서열을 정렬했다(번호 체계가 달라도 처리). shell 잔기 <3개인 copy는 정합에서 뺐다(1잔기 접촉이 오짝을 만든 사례 확인). 수용 사슬의 짝은 원자료(apo `apo_chain`)로 고정했다.

- 적격 51건 중 10건이 전체 pocket RMSD >0.75 Å였다: 3G22 U5P 0.90/0.94, 2PGA ANU 2.43–4.49, 6AYO/6AYQ C1Y·TDI 2.17–2.26, 7TX5 APR 0.76.
  - 6AYQ는 주 사슬 단독 정합에서도 apo 두 사슬 모두와 2.25 Å 차이가 난다. 실제 pocket 구조 변화다.
  - 2PGA는 apo dimer의 두 자리가 0.8 Å / 2.4–2.5 Å로 다르다. 원자료 사슬 기준 자리는 2.4–4.5 Å다. 짝 대응에 불확실성이 남는다.
- 입력 규칙(≤0.75 Å)을 전체 pocket에 올바르게 적용하면, 이 10건은 탈락이 아니라 `backbone_change` 층으로 옮긴다. 이에 따라 OMPDC(Q31b7)와 MTAN(Q8365) component의 적격 사례가 없어진다.
- 나머지 적격 사례의 리간드 이식 위치 변화는 ≤0.45 Å(4G27)였다. 대부분 <0.1 Å다.

## 2. pocket sidechain 상태 (`code/audit_qfit_states.py`)

pocket = assembly 내 polymer copy 중 리간드 reference pose(holo는 deposited pose, apo는 위 재정합 pose)와 4.5 Å 이내인 잔기(모든 사슬)를 qFit 모델 잔기로 대응시킨 것이다.
상태 = 점유율 ≥0.10인 altloc의 χ1 well, χ1/χ2 joint well이다(`scripts/105_joint_rotamers.py`의 정의). 다중 상태 = well이 2개 이상 점유된 잔기다.
검사: 수용 사슬 χ1 well 1,299잔기가 replicate cache와 전부 일치했다. 모든 pocket 잔기와 리간드·주변 HET이 qFit 모델에 있다. 예외로 7TX5 holo qFit의 cell이 deposited와 다르다.

| 집합 | cases / component | holo 다중(≥1잔기) | apo 다중 | 둘 다 | 어느 한쪽 |
|---|---|---|---|---|---|
| 적격 41 (χ1χ2) | 41 / 13 | 32 / 12 | 29 / 12 | 26 / 11 | 35 / 13 |
| 적격 41 (χ1만) | 41 / 13 | 24 / 11 | 21 / 9 | 19 / 7 | 26 / 11 |
| 기존 대표 중 적격 13 (χ1χ2) | 13 / 13 | 11 | 8 | 8 | 11 |
| 기존 대표 중 적격 13 (χ1만) | 13 / 13 | 9 | 7 | 7 | 9 |

- 기존 대표 중 다중 상태가 어느 쪽에도 없는 사례: 5NY7 NCA(amidase), 4PM7 CE3(CTX-M-14). 6EEV MEV와 6F66 CUT:B, 8TNB IQM:C는 holo에만 있다.
  - amidase component는 적격 3건 모두 다중 상태가 거의 없다(5NYB apo χ1χ2 1잔기).
  - CTX-M component에는 다중 상태가 있는 사례(6C79, 3HUO PNN:B, 3HVF)가 있다.
- 다중 상태가 없는 사례는 그 backbone에서 단일 = 다중이므로 차이 0의 대조로 쓸 수 있다.
- 주의: qFit은 pocket backbone에도 altloc을 붙인다. 적격 41건 중 holo 31건, apo 38건의 pocket에서 altloc 간 backbone 원자 거리가 최대 >0.5 Å이다. 다중 상태 receptor를 altloc 전체 conformer로 만들지(backbone이 조금씩 다름), major backbone 위 sidechain만 바꿀지 결정이 필요하다.
- χ1χ2 조합 수는 1부터 139,968(5HIH apo)까지다. receptor 개수 예산과 조합 방법을 결과 전에 고정해야 한다.

## 산출물

`qfit_state_audit.json`(사례·잔기별 상태, cache 대조), `pose_refit.json`(재정합 pose와 사슬별 RMSD, apo 상태), `state_summary.tsv`/`.json`(2×2 판정·집계).
코드: `code/audit_qfit_states.py`, `code/refit_pose.py`, `code/diagnose_refit.py`(진단용), `code/summarise_states.py`. Slurm: 438518, 438520(CPU, 예외 0).

## 결정 (사용자, 2026-10-28)

1. 단일 상태 = 같은 qFit 모델의 major conformer.
2. 다중 상태 receptor = major backbone 고정, sidechain만 상태별로 바꾼다. 각 상태의 χ 각도를 major backbone 위에 표준 결합 기하로 다시 세운다(altloc 좌표를 그대로 붙이지 않음). 단일 arm도 같은 절차로 다시 세워 구축 오차를 맞춘다. 충돌을 검사한다. 상태와 함께 움직이는 작은 backbone 변화(backrub 등)는 반영하지 않는다(한계).
3. 전체 pocket RMSD >0.75 Å인 10건(3G22 U5P×2, 2PGA ANU×3, 6AYO C1Y×2, 6AYQ TDI×2, 7TX5 APR)은 backbone_change 층으로 포함한다. 주 cohort와 분리해 보고한다. 2PGA는 apo 자리 대응이 불확실하므로(원자료 사슬 기준 2.4–4.5 Å, 다른 대칭 자리 0.8 Å), 이 층 안에서도 그 점을 명시한다.

## 남은 결정

- 대표 재선정: 새 설계에서 '다중 상태 존재(처리 변수)'를 대표 순서에 넣을지. 결과·minor→major 신호는 계속 쓰지 않는다.
- 다중 상태 조합 방법과 receptor 예산, 단일 arm과 탐색량을 맞추는 방법.
