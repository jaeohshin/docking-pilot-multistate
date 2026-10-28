# 2×2 도킹 분석 계획 (사전 고정, 2026-10-28 밤)

**도킹 결과를 보기 전에 고정한다.** 이후 수정은 무엇을·왜·결과를 보기 전인지 후인지와 함께 적는다.

## 묻는 것

> 결정 구조에서 관측된 **여러 sidechain 상태**를 수용체에 넣으면, 단일 상태보다 리간드 pose를 더 잘 회복하는가?

설계: **holo/apo backbone × single/multiple sidechain 상태**. packer는 쓰지 않고 qFit-at-scale 관측 상태를 그대로 쓴다.

## 비교 (사전 고정)

| 대비 | 묻는 것 |
|---|---|
| **multi − single (apo)** | **주 대비.** 실제 쓰임새에 해당 |
| multi − single (holo) | 상한선. 리간드에 맞는 backbone에서의 효과 |
| apo − holo (single) | backbone이 주는 난이도 |
| 상호작용 | 다중 상태의 이득이 apo에서 더 큰가 (가설의 핵심) |

## 탐색량 맞춤

arm마다 도킹 호출 **8회**로 동일. exhaustiveness 32, num_modes 9, box 동일.
- single: 수용체 1개 × seed 8개
- multi: 수용체 n개에 8회를 round-robin. n=1이면 single과 정확히 같아진다.

즉 다중 상태 arm이 **더 많이 탐색해서 이기는** 일은 구조적으로 불가능하다.

## 지표

- **주 지표: 8회 중 docking score 1위 pose의 symmetry-aware heavy-atom RMSD** (정합 없음, RDKit `CalcRMS`).
- 성공 기준 **≤ 2 Å** (표준).
- 보조: 8회 중 최저 RMSD(oracle), 상위 pose score, 성공률.
- 사례별 쌍 비교(같은 사례의 multi − single)를 쓰고, 사례 간 bootstrap으로 CI를 낸다.

## 판정

- **주 판정: apo에서 multi − single의 평균 RMSD 개선이 95% CI 하한 > 0.**
- 대표 17건은 homology component당 1건이라 독립이지만 n이 작다. **유의하지 않아도 "효과 없음"의 증명이 아니다**라고 적는다.
- 층은 분리해 보고한다: 주 cohort 13건, backbone 변화 층 4건(2PGA는 apo 자리 대응 불확실도 함께).

## 층과 사례

- 주 cohort 13건: apo–holo 전체 pocket backbone RMSD ≤ 0.75 Å
- backbone 변화 층 4건: 3G22(0.90 Å), 2PGA(4.49 Å), 6AYQ(2.19 Å), 7TX5(0.76 Å). 주 결과와 섞지 않는다.
- 다중 상태가 없어 multi = single인 사례는 **차이 0의 대조**로 그대로 둔다(5NY7 양쪽, 3G22 holo 등).
  이들을 제외하면 평균이 부풀려지므로 제외하지 않는다. 다만 "상태가 있는 사례만" 부분집합도 함께 보고한다.

## 결과와 무관하게 보고할 것

- 상태 조합은 K=8까지만 쓰고, 잔기별 well 점유율로 독립 추출했다. 잔기 간 상관은 반영하지 않는다.
- major backbone에 sidechain만 이식했으므로 상태와 함께 움직이는 작은 backbone 변화(backrub 등)는 빠진다.
- 수소는 PDB2PQR(--noopt --nodebump)로 붙이고 회전 가능한 것만 충돌 해소로 돌렸다. His 토토머는 PROPKA 결정을 따른다.
- 리간드는 결정 pose에서 시작하지 않도록 무작위 conformer에서 출발한다.
- 대표 선정에 minor→major 신호와 도킹 결과를 쓰지 않았다(`final_selection_v1`).

## 하지 않는 것

- 결과를 보고 지표·임계값·층·사례를 바꾸지 않는다.
- flexres arm과 packer B/C 비교는 이번 범위 밖이다.
