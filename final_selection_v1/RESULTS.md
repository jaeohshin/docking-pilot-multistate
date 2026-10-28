# Docking pilot 최종 후보 선정 (final_selection_v1)

2026-10-28 KST. 입력: `final_review_v1/candidates79.json`(79건/23 component, 고정). 도킹·모델 추론·학습은 실행하지 않았다.
기존 폴더(`final_review_v1`, `filter_review_v1`, `selection_v3`, `run_v*`)는 읽기만 했다. commit/push하지 않았다.

## 결론

| 집합 | cases | 비상동성 component | minor→major 양성 cases | 양성 component | 양성 잔기 |
|---|---:|---:|---:|---:|---:|
| 입력 79건 | 79 | 23 | 4 | 4 | 8 |
| 최종 적격(accept) | 51 | 15 | 1 | 1 | 3 |
| pilot 대표(component당 1) | 15 | 15 | 1 | 1 | 3 |
| 대표 민감도(holo-only 환경 등급) | 15 | 15 | 0 | 0 | 0 |

- 판정: accept 51, exclude 16, defer 12. 적격이 없는 component 8개(B3P, BTB, SKM, ZKR, 9JE/NPO, MXX, TCL, CFF 계열).
- 적격 중 minor-state 양성은 MERS-CoV macrodomain 5HIH→5HOL APR:A:201 한 건뿐이다. 이 사례가 대표가 된 것은 apo 측 결정 접촉을 환경 등급에 넣은 규칙 때문이다. holo 측만 쓰면 같은 component에서 8SH6→8SH8 AR6:A:201이 대표가 되고, 대표 중 양성은 0건이다. 두 결과를 모두 보고한다. 대표 선정에는 minor-state 신호를 쓰지 않았다(join은 선정 후).
- 나머지 양성 3건은 일반 cohort에서 빠진다: 3HUO PNN:A:302(촉매부위 밖 2차 결합, 보류), 6QO9 BTB:B:401(완충제, 제외), 8TAV ZKR:A:303(결정 접촉 위 비공유 여분 분자, 제외). mechanism subset(`mechanism_subset.tsv`)에 상태와 함께 따로 남겼다.

## minor-state 잔기 — 실제 PDB 번호

원자료 `i`는 holo qFit 사슬 polymer의 0-based 위치다(`scripts/180_apo_holo.py`, `enumerate(ch.get_polymer())`). PDB 번호가 아니다.
79건 전부를 replicate cache의 χ1 well에서 다시 계산했고, 4건·8잔기가 후보 기록과 일치했다(`minor_state_check.json`).
번호는 qFit 모델의 해당 위치 잔기에서 읽고, deposited mmCIF에서 같은 이름·번호를 확인했다.

| case | 판정 | apo → holo, ligand | 잔기 | holo 번호 | apo 번호 | 원자료 i | holo χ1 well | apo χ1 well |
|---|---|---|---|---|---|---:|---|---|
| fe02a159def13cce | **선정·대표** | 5HIH → 5HOL APR:A:201 | Ser | A:Ser127 | A:Ser127 | 127 | m .88, p .12 | m .25, **p .75** |
| | | | Ile | A:Ile130 | A:Ile130 | 130 | m .60, p .40 | m .38, **p .62** |
| | | | Val | A:Val154 | A:Val154 | 154 | m 1.0 | m .25, **p .51**, t .24 |
| fe0893e9d4e1f863 | 보류 | 3HRE → 3HUO PNN:A:302 | Gln | A:Gln188 | A:Gln188 | 161 | t 1.0 | **m .51**, t .49 |
| 475fd70298075f1c | 제외 | 6QO5 → 6QO9 BTB:B:401 | Val | B:Val75 | B:Val75 | 74 | t 1.0 | **p .51**, t .49 |
| | | | Leu | B:Leu77 | B:Leu77 | 76 | t 1.0 | **m .51**, t .49 |
| | | | Lys | B:Lys83 | B:Lys83 | 82 | t 1.0 | **m .51**, t .49 |
| 8c7d7558c0c71f80 | 제외 | 8FTP → 8TAV ZKR:A:303 | Lys | A:Lys70 | A:Lys70 | 65 | m .52, t .48 | t .51, m .49 (major t) |

굵게 표시한 것이 apo major well. 5HOL 3잔기는 원자료 i와 PDB 번호가 우연히 같다. 3HUO는 i=161 → Gln188, 8TAV는 i=65 → Lys70.
6QO9 BTB 사례에서 원자료는 apo 6QO5의 B 사슬을 썼다(후보 기록의 `apo_chain`은 A). 원자료 규칙(정렬이 가장 긴 apo 사슬)대로 B 사슬 번호를 표시했다.
주의: 3HUO·6QO9·8TAV의 apo well은 0.49/0.51처럼 거의 반반이라 'minor'가 경계값이다. χ1만의 정의이며 χ2–4, 전체 pose, 동역학적 전이나 conformational selection을 입증하지 않는다.

## 판정 규칙 (코드 순서대로)

1. 리간드 역할(`code/ligand_roles.py`, 사례별 수동 판정·근거 수준 L/S/U). 기질·생성물·유사체·TS analogue·억제제·검증된 fragment가 기능 부위에 있으면 accept. 완충제·결정화 첨가물·비기능 부위 여분 분자는 exclude. 근거가 없으면 defer.
2. holo 결정 접촉: 리간드가 속한 deposited assembly 밖 polymer copy가 리간드 단백질 접촉 원자의 ≥15% 또는 리간드 원자의 ≥25%를 4.5 Å 이내에서 접촉하면 major → exclude(crystal-context 층). 감사 코드에 결과 전에 고정했다.
3. 보조인자·공동 리간드: listed cofactor 또는 ≥6 heavy atom HET이 리간드 원자의 ≥10%와 4 Å 이내이고 apo 해당 위치에 같은 성분이 없으면 defer. holo HET을 apo에 자동 이식하지 않았다.
4. assembly 내 partner 사슬이 리간드와 접촉하는데 apo에 같은 단백질(SEQRES 동일성 ≥95%)의 assembly 내 copy가 그 자리에 없으면 defer. 해당 사례는 없었다.
5. 검증 보고서 RSCC < 0.80이면 defer.
6. 대표: (pose 명확성[RSCC ≥0.90, 평균 점유율 ≥0.90] 아님, 환경 등급 max(holo, apo 결정 의존) none<minor<major, apo/holo 중 나쁜 해상도, case_id) 최소. apo 측을 등급에 넣는 규칙은 감사 요약을 본 뒤, 대표 목록을 만들기 전에 정했다. holo-only 등급 결과를 민감도로 함께 보고한다(Q92b4 macrodomain: 5HOL ↔ 8SH8, Qf444 EvaD: 같은 TMP:A:1204 리간드의 수용 사슬 A ↔ B).

## 조립체·결정 접촉 감사 방법 (`code/audit_assembly_contacts.py`)

- deposited mmCIF(83개, `final_review_v1/cifs`), gemmi 0.7.5. 모든 space-group 연산 × 근방 격자 이동으로 모든 subchain copy를 명시적으로 생성했다. image_idx에 의존하지 않았다. ASU copy와 symmetry copy를 구분하고, 각 copy가 deposited assembly에서 리간드(apo는 수용 사슬) copy와 같은 assembly instance에 속하는지 판정했다(g_anchor⁻¹·g_copy = 결정 변환, 원자 좌표 1 Å 허용).
- apo: holo 리간드 pose를 `residue_map`의 author key로 deposited apo pocket backbone에 맞춰 옮겼다(79건 모두 매핑 잔기 전부 사용, deposited RMSD 0.04–0.72 Å). qFit 기반 변환과의 리간드 위치 차이는 77건에서 ≤0.34 Å였다. 예외는 6AYQ TDI:B:301/수용 A(0.98 Å)와 6AYO C1Y:B:301/수용 A(0.45 Å)로, 둘 다 4잔기 fit이고 대표가 아니다. interface 부위는 이후 reference 정합에 두 사슬 pocket을 함께 써야 한다.
- Slurm 438495, CPU 2, 예외 0. 검증 보고서 46개는 RCSB 공개 파일(`validation_xml/`)이다. 7TX5는 RSCC가 없다(중성자/X선 구조).

## 주요 제외·보류 근거

- 완충제·첨가물: 4MNK B3P, 6QO9 BTB, 6AYO BTB:B:305, 8OE5 FLC(침전제 1.02 M citrate가 BlaC 활성부위에 결합).
- 여분 분자: 5NY7 NCA:A:303/304, 3HUO PNN:A:304, PNK:B:305, 8TAV ZKR:A:302/303. 8TAV에서 공유결합 활성부위 adduct는 ZKR:A:301(B–Ser93 OG, B–His220 NE2)이고 후보 302/303은 비공유 여분 분자다.
- 결정 접촉 major(역할은 유효): 2IYS SKM(리간드 원자 42%), 5XH2 NPO(60%), 7JME CMP(27%), 1OI6 TMP:B:1203(38%, 2건).
- 보조인자 부재: 3GAM MXX(FAD, 리간드 원자 82%), 2WYW TCL(NAD, 65%), 6F68 4EU(GSH, 18%). 공동 리간드: 3G6M CFF 두 분자는 문헌상 둘 다 groove(subsite −1~+1)에 함께 쌓여 서로 접촉한다(64%/21%). FAD/NAD 등을 B/C 공통으로 이식하는 receptor recipe를 사용자 결정으로 허용하면 이 7건/4 component가 cofactor 층으로 추가될 수 있다.
- 역할 미확인: 3HUO PNN:A:302(촉매부위 밖, 기능 근거 없음), 7YKO 9JE(미출판).

## 대표 15개의 receptor recipe 요약

`pilot_manifest.tsv`의 `apo_receptor_chains`가 apo deposited 사슬(필요하면 assembly 대칭 copy)을 준다. 모든 HET과 물은 제거한다. holo 주변 HET(예: sulfate, GSH)은 `holo_het_neighbours_not_in_receptor`에 기록했다. interface 대표: 3G22 U5P(A+B), 6AYQ TDI(A+B), 4G27 PHU(B+R, CaM–SK2), 1OI6 TMP(A+B), 3UAW ADN·6EEV MEV·5A9T YLG(assembly 대칭 copy 포함). protonation, 원자 typing, 엔진 준비는 아직 검증하지 않았다.

## 산출물

`case_review_all79.tsv`(79건 전체), `eligible_cases.tsv`(51), `pilot_manifest.tsv`(15), `mechanism_subset.tsv`(양성 4), `minor_state_residues.tsv`, `minor_state_check.json`, `assembly_contact_audit.json`(사례별 copy·assembly·접촉 원자료와 입력 SHA256), `selection_summary.json`, `files_sha256.json`.
재현: `sbatch code/audit.sbatch` → `sbatch code/minor.sbatch` → `python code/select_pilot.py`(multistate-analysis env). `code/fetch_validation.py`는 검증 보고서를 받는다(이미 있으면 건너뜀).

## 한계

- 리간드 역할은 entry 제목·primary citation·좌표 부위로 판정했다. 원논문 본문까지 모두 확인하지는 않았다(근거 수준 열 참조). 8SH8은 미출판이다.
- 결정 접촉 임계값(15%/25%)은 사전 고정한 운영 기준이고, packing의 인과 의존을 입증하지 않는다. 경계 사례: 7JME 0.27, 8TAV ZKR:A:302 0.25.
- 대표 선정 규칙의 apo 등급 부분은 감사 결과를 본 뒤 정했다(위 민감도).
- 적격 51건 중 6NT 11건, PARPi binder 8건, 뉴클레오사이드 phosphorylase 8건처럼 한 component에 몰려 있다. 대표 15개는 component당 1건이다.

## 79건 전체 판정표

| component | case_id | apo>holo | ligand | 역할 | 판정 | 이유 | 대표 |
|---|---|---|---|---|---|---|---|
| 06d8f0a1 | 7ef4c23f14dbfb0c | 4mnj>4mnk A | B3P:A:407 2-[3-(2-HYDROXY-1,1-DIHYDROX | buffer | 제외 | ligand role: buffer |  |
| 071e2319 | 97ea259b19814c3e | 5nxz>5ny7 A | NCA:A:302 NICOTINAMIDE | substrate | 선정 | 적격 | 대표 |
| 071e2319 | 60e537d07177c7f6 | 5nxz>5nyb A | 9EE:A:301 adipamide | substrate | 선정 | 적격 | 예비 |
| 071e2319 | 84a6cc4cec17be6f | 5nxz>4izs A | BMD:A:304 BUTYRAMIDE | substrate | 선정 | 적격 | 예비 |
| 071e2319 | 1df758496479cfac | 5nxz>5ny7 A | NCA:A:304 NICOTINAMIDE | extra copy, surface | 제외 | ligand role: extra copy, surface |  |
| 071e2319 | 223bcc6ed7796669 | 5nxz>5ny7 A | NCA:A:303 NICOTINAMIDE | extra copy, secondary site | 제외 | ligand role: extra copy, secondary site |  |
| 0d511715 | 0c2914f1336e5905 | 4pm8>4pm7 A | CE3:A:301 (6R,7R)-3-(acetyloxymethyl)- | substrate | 선정 | 적격 | 대표 |
| 0d511715 | 01100dcbf8cafc35 | 3hre>3hvf B | PNK:B:302 (2R,4S)-2-{(R)-carboxy[(phen | product | 선정 | 적격 | 예비 |
| 0d511715 | 129663a7f774bc71 | 4pm6>4pm5 A | CE3:A:301 (6R,7R)-3-(acetyloxymethyl)- | substrate | 선정 | 적격 | 예비 |
| 0d511715 | 709de94af951aa44 | 6c7a>6c79 A | CE3:A:301 (6R,7R)-3-(acetyloxymethyl)- | substrate | 선정 | 적격 | 예비 |
| 0d511715 | 8b520349f920f468 | 3hre>3huo B | PNN:B:301 PENICILLIN G | substrate | 선정 | 적격 | 예비 |
| 0d511715 | fe0893e9d4e1f863 | 3hre>3huo A | PNN:A:302 PENICILLIN G | extra copy, secondary site | 보류 | ligand role: extra copy, secondary site |  |
| 0d511715 | af5b3fa499286445 | 3hre>3huo B | PNK:B:305 (2R,4S)-2-{(R)-carboxy[(phen | extra copy, surface | 제외 | ligand role: extra copy, surface |  |
| 0d511715 | c34cfe62e60dc71d | 8oe1>8oe5 C | FLC:C:301 CITRATE ANION | crystallisation additive | 제외 | ligand role: crystallisation additive |  |
| 0d511715 | cd901008b90ffbd4 | 3hre>3huo A | PNN:A:304 PENICILLIN G | extra copy, surface | 제외 | ligand role: extra copy, surface |  |
| 0f715600 | 35e1ea353db1a647 | 5f6e>5f6v A | 5VL:A:201 2-(2-hydroxyphenyl)phenol | validated allosteric fragment | 선정 | 적격 | 대표 |
| 18854647 | 3cb6f9e6ad8d8368 | 4qod>3gam A | MXX:B:231 5,8-dimethoxy-1,4-dimethylqu | inhibitor | 보류 | needs absent FAD (0.82 of ligand atoms <= 4 A) |  |
| 18854647 | 4bc0799f7dadd332 | 4qod>3gam B | MXX:B:231 5,8-dimethoxy-1,4-dimethylqu | inhibitor | 보류 | needs absent FAD (0.82 of ligand atoms <= 4 A) |  |
| 18854647 | 75c6cb1919e3410c | 4qod>3gam B | MXX:A:231 5,8-dimethoxy-1,4-dimethylqu | inhibitor | 보류 | needs absent FAD (0.82 of ligand atoms <= 4 A) |  |
| 244ca05b | a55370d81109dc3b | 8usg>8ush A | 6NT:A:401 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 대표 |
| 244ca05b | 28cfc7bf72ea90b1 | 5rg6>5rgc A | 6NT:A:401 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | 2f558d042f16f4f7 | 7k4t>7k4z A | 6NT:A:402 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | 85dc1e77a1cdb21d | 5rg9>5rgf B | 6NT:B:500 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | ac718cff8957d4c8 | 5rg6>5rgc B | 6NT:B:401 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | ca403edcbb1150e0 | 8fme>8fos A | 6NT:A:502 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | e71ec1febb7d05e2 | 5rg9>5rgf A | 6NT:A:401 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | f6bc7e3de3c6c8f1 | 8usg>8ush B | 6NT:B:401 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 244ca05b | fd5e7dc6163d1398 | 5rg8>5rge A | 6NT:A:500 6-NITROBENZOTRIAZOLE | TS analogue | 선정 | 적격 | 예비 |
| 2cf9871b | 169e40ddaed42efc | 2wyu>2wyw B | TCL:B:1260 TRICLOSAN | inhibitor | 보류 | needs absent NAD (0.65 of ligand atoms <= 4 A) |  |
| 2cf9871b | ac4609eed42760d7 | 2wyu>2wyw A | TCL:A:1261 TRICLOSAN | inhibitor | 보류 | needs absent NAD (0.65 of ligand atoms <= 4 A) |  |
| 2cf9871b | b18fd2e3f40c1574 | 2wyu>2wyw D | TCL:D:1260 TRICLOSAN | inhibitor | 보류 | needs absent NAD (0.65 of ligand atoms <= 4 A) |  |
| 2cf9871b | cdb554d7c5589c57 | 2wyu>2wyw C | TCL:C:1260 TRICLOSAN | inhibitor | 보류 | needs absent NAD (0.65 of ligand atoms <= 4 A) |  |
| 2e938254 | 4e2004374c441952 | 4j9y>4g27 B | PHU:R:1006 1-phenylurea | weak modulator / pocket fragment | 선정 | 적격 | 대표 |
| 31b796f7 | 2248260e300938dd | 3g1y>3g22 B | U5P:A:229 URIDINE-5'-MONOPHOSPHATE | product | 선정 | 적격 | 대표 |
| 31b796f7 | ce8f0e29cd1dca59 | 3g1y>3g22 A | U5P:B:229 URIDINE-5'-MONOPHOSPHATE | product | 선정 | 적격 | 예비 |
| 403d6c2f | 172d7da9c683ce21 | 5fsv>5fsy A | AR6:A:1262 [(2R,3S,4R,5R)-5-(6-AMINOPUR | substrate/product | 선정 | 적격 | 대표 |
| 4357fad1 | 274d968ed2de6c64 | 3uav>3uaw A | ADN:A:236 ADENOSINE | substrate | 선정 | 적격 | 대표 |
| 4357fad1 | 3501c3c7d181fc13 | 3dps>2pga D | ANU:B:7016 2,2'-Anhydro-(1-beta-D-arabi | inhibitor | 선정 | 적격 | 예비 |
| 4357fad1 | 9d217113b12c58ef | 3uav>2ac7 A | ADN:B:1215 ADENOSINE | substrate | 선정 | 적격 | 예비 |
| 4357fad1 | ccec5c8aea863078 | 3dps>2pga B | ANU:D:7013 2,2'-Anhydro-(1-beta-D-arabi | inhibitor | 선정 | 적격 | 예비 |
| 4357fad1 | d044697c9b260b63 | 3uav>2ac7 B | ADN:B:1215 ADENOSINE | substrate | 선정 | 적격 | 예비 |
| 4357fad1 | db219b9e4c93bafe | 3dps>2pga A | ANU:F:7014 2,2'-Anhydro-(1-beta-D-arabi | inhibitor | 선정 | 적격 | 예비 |
| 4357fad1 | ea2aa0effa40a210 | 3uav>2ac7 B | ADN:A:1216 ADENOSINE | substrate | 선정 | 적격 | 예비 |
| 4357fad1 | f103150d4318d464 | 3uav>3uax A | NOS:A:236 INOSINE | substrate | 선정 | 적격 | 예비 |
| 4819d46c | 2a61e70b2504eba1 | 6eeu>6eev A | MEV:A:501 (R)-MEVALONATE | product | 선정 | 적격 | 대표 |
| 5908f5eb | 475fd70298075f1c | 6qo5>6qo9 B | BTB:B:401 2-[BIS-(2-HYDROXY-ETHYL)-AMI | buffer | 제외 | ligand role: buffer |  |
| 5908f5eb | ad802dab0efb57b5 | 6qo5>6qo9 A | BTB:B:401 2-[BIS-(2-HYDROXY-ETHYL)-AMI | buffer | 제외 | ligand role: buffer |  |
| 664a4fbc | 827d8b35f26be3bf | 2iyt>2iys A | SKM:A:201 (3R,4S,5R)-3,4,5-TRIHYDROXYC | substrate | 제외 | crystal-context dependent: 0.13 of protein contact atoms and up to 0.42 of ligand atoms from non-assembly copies |  |
| 806f11f1 | 3d14b6f21e8b903c | 8ftp>8tav A | ZKR:A:302 [5-(trifluoromethyl)thiophen | extra copy, secondary site | 제외 | ligand role: extra copy, secondary site |  |
| 806f11f1 | 8c7d7558c0c71f80 | 8ftp>8tav A | ZKR:A:303 [5-(trifluoromethyl)thiophen | extra copy, secondary site | 제외 | ligand role: extra copy, secondary site |  |
| 836548e7 | 6a0fb1f5b4c558d2 | 6aym>6ayq B | TDI:A:301 (3R,4S)-1-[(4-AMINO-5H-PYRRO | TS analogue | 선정 | 적격 | 대표 |
| 836548e7 | ae557f538c35d329 | 6aym>6ayo B | C1Y:A:301 (3R,4S)-1-[(4-amino-5H-pyrro | TS analogue | 선정 | 적격 | 예비 |
| 836548e7 | b3a33fcc7a7495b2 | 6aym>6ayo A | C1Y:B:301 (3R,4S)-1-[(4-amino-5H-pyrro | TS analogue | 선정 | 적격 | 예비 |
| 836548e7 | c58ec64e1ac75a51 | 6aym>6ayq A | TDI:B:301 (3R,4S)-1-[(4-AMINO-5H-PYRRO | TS analogue | 선정 | 적격 | 예비 |
| 836548e7 | d965c9e63efef324 | 6aym>6ayo B | BTB:B:305 2-[BIS-(2-HYDROXY-ETHYL)-AMI | buffer | 제외 | ligand role: buffer |  |
| 92b4b8ac | fe02a159def13cce | 5hih>5hol A | APR:A:201 ADENOSINE-5-DIPHOSPHORIBOSE | substrate/product | 선정 | 적격 | 대표 |
| 92b4b8ac | 49e6f53ebc0b9814 | 8sh6>8sh8 A | AR6:A:201 [(2R,3S,4R,5R)-5-(6-AMINOPUR | substrate/product | 선정 | 적격 | 예비 |
| 92b4b8ac | 95ffd089bd47d70a | 7tx4>7tx5 A | APR:A:201 ADENOSINE-5-DIPHOSPHORIBOSE | substrate/product | 선정 | 적격 | 예비 |
| 92b4b8ac | ce056305803760d9 | 6wey>7jme A | CMP:A:401 ADENOSINE-3',5'-CYCLIC-MONOP | ligand (soaked) | 제외 | crystal-context dependent: 0.11 of protein contact atoms and up to 0.27 of ligand atoms from non-assembly copies |  |
| 94771a75 | cc74cee67640cfef | 7ykp>7yko A | 9JE:A:302 pentane-1,5-diol | unverified | 보류 | ligand role: unverified |  |
| 94771a75 | e3f6b2caa033ba2b | 5xfz>5xh2 A | NPO:A:305 P-NITROPHENOL | product | 제외 | crystal-context dependent: 0.35 of protein contact atoms and up to 0.60 of ligand atoms from non-assembly copies |  |
| b658135f | 02232eaff63e6f6a | 8tn1>8tnb C | IQM:C:201 5-fluoro-2-{4-[(methylamino) | designed-binder ligand | 선정 | 적격 | 대표 |
| b658135f | 2554c2e878337eda | 8tn1>8tn6 B | RPB:B:201 Rucaparib | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | 3defdfd8f936876e | 8tn1>8tnb A | IQM:A:201 5-fluoro-2-{4-[(methylamino) | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | 3e4042944ebec322 | 8tn1>8tn6 A | RPB:A:201 Rucaparib | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | 69bd5aa33dfdd4f3 | 8tn1>8tnb B | IQM:B:201 5-fluoro-2-{4-[(methylamino) | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | 7d46addf82382500 | 8tn1>8tnc C | 3JD:C:201 2-{4-[(3S)-piperidin-3-yl]ph | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | ab769351ccfcd8b9 | 8tn1>8tn6 C | RPB:C:201 Rucaparib | designed-binder ligand | 선정 | 적격 | 예비 |
| b658135f | eb0590317573e173 | 8tn1>8tnc A | 3JD:A:201 2-{4-[(3S)-piperidin-3-yl]ph | designed-binder ligand | 선정 | 적격 | 예비 |
| c50f788b | 1aa414d8b85e516d | 6f43>6f66 B | CUT:B:300 [2,4-bis(oxidanyl)phenyl]-ph | substrate analogue | 선정 | 적격 | 대표 |
| c50f788b | d1ef3ccf4bf68009 | 6f43>6f66 A | CUT:A:300 [2,4-bis(oxidanyl)phenyl]-ph | substrate analogue | 선정 | 적격 | 예비 |
| c50f788b | ff4d22e8ce445a8b | 6f43>6f68 B | 4EU:B:302 (2,4-dihydroxyphenyl)(4-hydr | substrate analogue | 보류 | needs absent GSH (0.18 of ligand atoms <= 4 A) |  |
| c73cfcfa | 8b5a09ff83941c82 | 5a9r>5a9t A | YLG:A:1291 (1R)-1-methyl-1,2,3,4-tetrah | product | 선정 | 적격 | 대표 |
| ee1463f6 | 28df252a153e6e16 | 3g6l>3g6m A | CFF:A:427 CAFFEINE | inhibitor | 보류 | needs absent CFF (0.21 of ligand atoms <= 4 A) |  |
| ee1463f6 | 7968c75e320af942 | 3g6l>3g6m A | CFF:A:1 CAFFEINE | inhibitor | 보류 | needs absent CFF (0.64 of ligand atoms <= 4 A) |  |
| f44483b2 | 8ec051785e59eb5b | 1ofn>1oi6 A | TMP:A:1204 THYMIDINE-5'-PHOSPHATE | substrate fragment | 선정 | 적격 | 대표 |
| f44483b2 | 2a6070389276a7a4 | 1ofn>1oi6 B | TMP:A:1204 THYMIDINE-5'-PHOSPHATE | substrate fragment | 선정 | 적격 | 예비 |
| f44483b2 | 06f8ca4f58af2891 | 1ofn>1oi6 A | TMP:B:1203 THYMIDINE-5'-PHOSPHATE | substrate fragment | 제외 | crystal-context dependent: 0.15 of protein contact atoms and up to 0.38 of ligand atoms from non-assembly copies |  |
| f44483b2 | 3f291c9437c562a0 | 1ofn>1oi6 B | TMP:B:1203 THYMIDINE-5'-PHOSPHATE | substrate fragment | 제외 | crystal-context dependent: 0.15 of protein contact atoms and up to 0.38 of ligand atoms from non-assembly copies |  |
