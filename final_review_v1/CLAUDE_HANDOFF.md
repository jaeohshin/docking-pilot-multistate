# Claude 인수인계: sidechain-multistate docking pilot 최종 후보 검토

작성: 2026-10-28. **상태: 79건에 대한 최종 생물학적 검토 미완료. 최종 pilot 수 미확정.**

## 사용자 요청과 작업 범위

사용자는 자신의 sidechain-multistate 연구에서 발견한 다중 sidechain 상태 정보가 실제 docking pose 회복에 도움이 되는지 평가하려 한다. 기존의 광범위한 backbone×sidechain packer 비교 아이디어는 잠시 보류했다.

현재 목적은 같은 apo backbone에서 single-state 학습 B와 multistate 학습 C의 sidechain ensemble을 비교할 docking pilot 후보를 선정하는 것이다. Virtual screening은 후속 단계다.

사용자의 누적 요청:
1. 완화된 필터를 통과한 79건의 리간드가 실제로 의미 있는 기질·생성물·억제제·기능성 fragment인지 확인한다.
2. biological assembly와 crystal packing 접촉을 검토하고, B/C에서 동일한 apo backbone·주변 환경으로 평가할 수 있는 최종 케이스를 선정한다.
3. 선정·제외·보류 이유를 표로 정리한다.
4. 최종 선정 사례 중 apo minor χ1 state가 holo major χ1 state와 일치하는 사례 수를 case/독립 서열군별로 집계하고 **잔기 이름과 실제 PDB 번호**도 표시한다.
5. 일반 성능 cohort는 이 minor-state 신호나 도킹 결과를 이용해 골라서는 안 된다. 신호 양성 사례는 선정 후 분석하거나 별도 mechanism subset으로 명시한다.

사용자는 지금 이 작업을 Claude에서 이어가도록 요청했다. 여기서 임의로 도킹·모델 추론·학습을 시작하지 말고, 먼저 최종 후보 선정과 재현 가능한 근거를 마무리한다. 이전 결과를 덮어쓰지 않고 새 출력 디렉터리를 사용한다. commit/push하지 않는다. 사용자에게는 간결한 한국어로 설명하고, 상세 근거는 experiments 문서에 남긴다.

## 위치와 접근

- 연구 저장소/원자료 클러스터: `/gpfs/jaeohshin/work/sidechain-multistate`
- SSH: `ssh jaeohshin@syntax.kias.re.kr` (Codex에서는 BatchMode=yes, ConnectTimeout=10으로 접속 성공)
- 데스크탑 저장소: `/data/work/sidechain-multistate`
- 실험 공통 상대 경로: `experiments/docking_pilot_20261028/`
- 이번 이어받기 묶음: `experiments/docking_pilot_20261028/final_review_v1/`
- 이전 필터 민감도 검토: `experiments/docking_pilot_20261028/filter_review_v1/`
- 이전 strict 선정: `experiments/docking_pilot_20261028/selection_v3/`

**`final_review_v1`은 최종 검토 작업 폴더 이름이며, 선정이 끝났다는 뜻이 아니다.**

Codex의 작업 cwd는 우연히 `/data/work/pythia_ppi`였지만 이번 연구와 무관하다. 그 저장소를 수정하지 않았으며 앞으로도 이번 작업을 그곳에서 수행하지 않는다. 별도 PreFold도 무관하다.

## 가장 먼저 읽을 파일

1. 이 문서.
2. `filter_review_v1/RESULTS.md`: 필터 정의·탈락·완화의 과학적 근거와 한계.
3. `filter_review_v1/review_summary.json`, `summarize_sensitivity.py`.
4. `final_review_v1/candidates79.json`: 이번 검토의 고정 입력 79건.
5. `final_review_v1/entry_metadata.json`: 83개 apo/holo 구조에서 추출한 title, citation, entity, assembly, crystallization 정보.
6. `final_review_v1/cifs/`: 83개 deposited mmCIF.gz. 원래 qFit PDB가 아님.
7. 필요할 때 `selection_v3/RESULTS.md`: 이전 strict 기준의 상세 절차.

## 지금까지 완료한 작업

### 원천과 상동성 감사

- 원 qFit apo–holo 데이터: 1,947 pairs, 1,266 source sequence groups.
- source pair selection은 그룹별 최고 해상도 apo와 최대 3개 holo를 택한다. 모든 가능한 pair를 열거한 것은 아니다.
- 비보조인자 및 양쪽 resolution≤2.5 Å: 3,527 ligand-chain cases / 1,118 source groups.
- train/validation 동일 그룹·직접 entry overlap 차단 후: 1,209 cases / 385 groups.
- full-SEQRES 동등성 확인, qFit chain A와 deposited auth chain AA 등의 alias rescue 완료.
- MMseqs reciprocal query/reference 검사: identity≥0.30, max(qcov,tcov)≥0.50, E≤0.001, sensitivity7.5. literal full-sequence duplicate도 별도 차단.
- 최종 homology-clean: 443 cases / 197 PDB pairs / 136 source groups / **78 candidate homology components**.
- 이전 중간 기록의72 components는 isolated query6개를 빠뜨린 오류였다. strict 최종2 components에는 영향 없음.
- frozen split과 FASTA provenance hash 검증 완료. 근거는 기존 보고서/JSON 참조.

### strict 필터와 완화 감사

기존 strict 기준에서 12 cases/7 source groups/2 homology components만 남았다. 이것은 준비 편의상 환경을 과도하게 제외한 결과였다.

443건에 대한 기존 순차 탈락: 매핑23 → 크기/occupancy142 → CCD분류/원자69 → 포켓매핑21 → backbone/junction/box41 → 환경123 → 최종좌표/공유결합QC12 → 12잔존.

환경 플래그는 다음과 같다:
- 다른 polymer chain이 ligand에서4.5 Å 이내.
- crystal image가 ligand에서4.5 Å 이내.
- metal이6 Å 이내.
- water 외 ≥6 heavy-atom HET residue가6 Å 이내.

이들은 생물학적으로 부적격이라는 뜻이 아니므로 추가 준비 검토로 바꾸었다. 기존 조기 탈락 때문에 생략했던 downstream geometry checks를 **443건에 재실행**했다. Slurm438440, CPU2개,2분38초,예외0,로컬 RDKit2025.09.3 후처리. strict12 case_id 재현 확인.

조건별 결과:
- strict:12cases/2components.
- 환경 제외 해제:67/22.
- +occupancy≥0.5:75/24.
- +heavy atoms6–50:83/25.
- +pocket backbone RMSD≤1.5 Å:87/25.
- **이번 고정 대상: 금속 근접 제외 유지, 나머지 환경은 검토로 전환, occupancy≥0.5, heavy atoms6–50, RMSD≤0.75 Å → 79cases/46pairs/37source groups/23components.**

79건은 `review_records.json`의 `review_scenarios`에 `no_metal_occ_size`가 있는 항목이며 `candidates79.json`으로 고정했다. 추가 상동성 기준은 완화하지 않았다.

계속 적용한 QC: contact4.5 Å에서 매핑/BB원자완전성, shell6.5 Å coverage≥95%, pocket N/CA/C/O fit, peptide junction, ligand CCD heavy-atom-name 완전성, deposited ligand entity, qFit–deposited ligand raw atom-name RMSD≤0.25 Å, explicit covalent/metal linkage와 <1.8 Å protein contact 검사, RDKit graph/allowed element/rotatable≤12. 예외 및 주의는 코드 참조.

### 이번 79건의 검토에서 실제로 한 일

- 79건 및 23개 component별 pair/ligand inventory 생성.
- 해당 apo/holo **83개 deposited mmCIF** 확보: holo46개는 클러스터 기존 cache에서 가져왔고, 빠진 apo37개는 RCSB 공식 파일에서 내려받았다.
- 모든83개에 대한 deposited citation/entity/assembly/crystallization 메타데이터를 `entry_metadata.json`으로 추출.
- 일부 핵심 원논문/RCSB 자료를 검색하고 아래 잠정 발견을 얻었다.
- **biological assembly를 실제 좌표로 재구성하고 crystal contact를 분해하는 계산은 아직 실행하지 않았다. 그 분석 스크립트도 아직 만들지 않았다.**
- 최종 ligand-role 판정, 최종 선정/제외표, 최종 대표 선정 모두 미완료.
- 도킹, PIPPack inference, 새 학습은 전혀 실행하지 않았다.

## 중요한 잠정 발견 — 확정 판정과 혼동하지 말 것

1. **8TAV/ZKR는 instance별 확인이 필수다.** 원논문은 covalent inhibitor N34를 다룬다. deposited entity details에는 boron이 Ser O/His N과 반응했다고 나온다. 하지만 이번 후보는 `ZKR:A:302`, `ZKR:A:303`이고 같은 구조에 ZKR3개가 있다. 반응한 active-site ligand와 표면의 별도 unreacted ligand가 다를 수 있다. 논문 제목만 보고 이 두 instance가 모두 covalent라고 단정하지 않는다. 원좌표·struct_conn·active-site 위치·원논문으로 각각 확인한다. 303은 minor-state 신호 양성이다.
   - https://www.rcsb.org/structure/8TAV
   - https://doi.org/10.1038/s41467-025-60367-3
   - PubMed40447595.
2. **BTB/B3P는 buffer 가능성이 높다.** B3P는4MNK, BTB는6QO9와6AYO에 있다. 4MNK의 연구 제목은 chitinase–(GlcNAc)3이고 후보 B3P가 그 목적 리간드는 아니다. 6QO9는 ribonucleotide reductase metal/redox 연구다. 6AYO의 실제 inhibitor C1Y와 BTB를 구분해야 한다. compound identity만으로 모든 주변 결합의 의미를 확정하지 않는다.
3. **5F6V/5VL는 의미 있는 allosteric fragment라는 근거가 있다.** Ubc9 crystallographic fragment screen의2,2′-biphenol. 단독 soaking 검증 및 같은 위치 유도체 결합이 보고되었다. crystal contact 의존 여부는 미검사.
   - https://www.osti.gov/servlets/purl/1249263 (저자 원고)
   - https://doi.org/10.1002/anie.201511351, PubMed27038327.
4. **4G27/PHU는 CaM–channel interface의 ligand**이고 추가 단백질 사슬이 실제 부위 구성요소일 수 있다. 결정 대칭체 접촉을 단순히 비생물학적이라고 처리하면 안 된다. apo4J9Y와 partner composition/sequence도 확인해야 한다.
   - https://www.rcsb.org/structure/4G27
   - https://doi.org/10.1038/ncomms2017, PubMed22929778.
5. **2WYW/TCL**은 NAD와 함께 결합한다. apo2WYU에서 필요한 NAD 환경이 확보되는지 확인해야 한다. holo cofactor를 무조건 apo로 이식해 통과시키지 않는다.
   - https://doi.org/10.1107/S1744309112033982, PubMed23027736.
6. **3GAM/MXX** quinone reductase2 inhibitor, 주변FAD와 dimer interface를 확인. **6F66/CUT,6F68/4EU** glutathione transferase polyphenol recognition이고 GSH가 있다. 이들은 리간드 의미와 receptor cofactor 유지가 별개 판정이다.
   - https://doi.org/10.1021/jm801335z, PubMed19265439.
   - https://doi.org/10.1038/s41598-018-26601-3, PubMed29855494.
7. **1OI6/TMP**는 EvaD의 nucleotide-binding substrate fragment/analogue 근거를 확인할 후보다. CIF에 citation이2개 있어 단일 `find_value`로는 제목이 안 나왔지만 `entry_metadata.json`에 모두 있다. primary PubMed15159413, DOI10.1074/JBC.M404091200.

이외 구조 제목과 primary citation은 전부 `entry_metadata.json`에 있으므로 기존 다운로드/검색을 처음부터 반복하지 말 것. 미출판 structure도 있으며 문헌 근거 강도를 분명하게 표시한다.

## 23개 component inventory

```
Q06d8f0a1f5ae8603  1  4mnj>4mnk:B3P
Q071e2319538cf827  5  5nxz>4izs:BMD | 5nxz>5ny7:NCA | 5nxz>5nyb:9EE
Q0d511715f969df6f  9  3hre>3huo:PNK/PNN | 3hre>3hvf:PNK | 4pm6>4pm5:CE3 | 4pm8>4pm7:CE3 | 6c7a>6c79:CE3 | 8oe1>8oe5:FLC
Q0f7156008b35f535  1  5f6e>5f6v:5VL
Q18854647e014fdb5  3  4qod>3gam:MXX
Q244ca05b311c4e65  9  5rg6>5rgc:6NT | 5rg8>5rge:6NT | 5rg9>5rgf:6NT | 7k4t>7k4z:6NT | 8fme>8fos:6NT | 8usg>8ush:6NT
Q2cf9871b35d95d27  4  2wyu>2wyw:TCL
Q2e938254667ffc06  1  4j9y>4g27:PHU
Q31b796f7552d4625  2  3g1y>3g22:U5P
Q403d6c2f54b32f0f  1  5fsv>5fsy:AR6
Q4357fad1fd6385d8  8  3dps>2pga:ANU | 3uav>2ac7:ADN | 3uav>3uaw:ADN | 3uav>3uax:NOS
Q4819d46ca698b01c  1  6eeu>6eev:MEV
Q5908f5eb64851452  2  6qo5>6qo9:BTB
Q664a4fbcd0ed5242  1  2iyt>2iys:SKM
Q806f11f1df0d798c  2  8ftp>8tav:ZKR
Q836548e7a7df94aa  5  6aym>6ayo:BTB/C1Y | 6aym>6ayq:TDI
Q92b4b8ac3490e160  4  5hih>5hol:APR | 6wey>7jme:CMP | 7tx4>7tx5:APR | 8sh6>8sh8:AR6
Q94771a75c08826ac  2  5xfz>5xh2:NPO | 7ykp>7yko:9JE
Qb658135f14f2b435  8  8tn1>8tn6:RPB | 8tn1>8tnb:IQM | 8tn1>8tnc:3JD
Qc50f788b8b4c3247  3  6f43>6f66:CUT | 6f43>6f68:4EU
Qc73cfcfab42ac422  1  5a9r>5a9t:YLG
Qee1463f68f4452fc  2  3g6l>3g6m:CFF
Qf44483b2d65dea51  4  1ofn>1oi6:TMP
```

## apo minor→holo major 신호: 아직 최종 선정 전 4건/4 components

| case_id | apo chain | holo chain | ligand | 원자료 selected_residues |
|---|---|---|---|---|
| fe0893e9d4e1f863 | 3hre A | 3huo A | PNN:A:302 | 161 |
| 475fd70298075f1c | 6qo5 A | 6qo9 B | BTB:B:401 | 74,76,82 |
| 8c7d7558c0c71f80 | 8ftp A | 8tav A | ZKR:A:303 | 65 |
| fe02a159def13cce | 5hih A | 5hol A | APR:A:201 | 127,130,154 |

**위 숫자를 곧바로 PDB author residue number로 출력하면 안 된다.** source `i`/cache index의0/1-based convention을 원 코드에서 확인하고 `residue_map`, `holo_key`, `apo_key`와 대응해 실제 residue name/number/insertion code를 복원해야 한다. apo/holo 번호도 각각 표시한다.

정의는 기존 source의 holo-defined4.5 Å contact에서 χ1 rotamer well의 apo minor occupancy≥0.1 및 holo major well 일치이다. 전체χ1–χ4, 완전한 sidechain pose 일치, 동역학적 전이 또는 인과적 conformational selection을 입증하는 것이 아니다. 가능하면 원 occupancy/well을 함께 재검증하고, 잔기별 source row도 보존한다.

최종 집계는 (a) 적격 전체cases, (b) 중복제거된 일반 pilot representatives, (c) 별도 mechanism subset을 구분한다. 일반 representative 선택 때문에 신호양성case가 다른같은계열대표에 밀릴 수 있다. 이를 숨기거나 대표를 사후 바꾸지 않는다.

## 이어서 해야 할 구체적인 검토

1. 79건 각각 ligand instance의 생물학적 역할과 결합 부위를 문헌/CIF로 판정한다. substrate/product/TS analogue/inhibitor/functional fragment는 질문에 맞게 허용하고 additive/surface opportunistic binding은 주 cohort에서 제외 또는 근거보류한다. inhibitor와 buffer가 같은entry에 있을 수 있다.
2. deposited biological assembly operator와 label/auth chain을 이용해 좌표를 재구성한다. assembly1을 무조건 신뢰하지 말고 author/software provenance, alternateassemblies, 생물학적 설명을 기록한다.
3. holo ligand 주변 protein contact를 실제assembly 내부 접촉과 외부crystalmate 접촉으로 분해한다. apo에도 옮긴 ligand reference 주변 contact를 계산하되 apo/holo/qFit/deposition좌표계 정합을 먼저 검증한다. 기존 geometryflag만으로 최종판정하지 않는다.
4. 필요한 partnerchain/cofactor가 apo에 존재하는지, 보존할환경과 제거할성분이 무엇인지 case별 receptor recipe를 명시한다. B/C가 같은 환경을 받도록 한다. crystalpacking이필수인경우는 주cohort제외/별도crystal-context분석으로 표시한다.
5. near-contact 하나만으로 packing의 인과적 의존을 단정하지 말고 contact수/거리/ligand원자분율/주요interactions와 문헌 근거를 함께 기록한다. threshold는 결과에 맞춰 조절하지 않는다.
6. 79개 case_id 모두에 accept/exclude/defer와 이유/출처를 남긴다. 중복구조는 부적격으로 처리하지 말고 eligible reserve인지 구분한다.
7. 적격 homology component마다 사전 명시된 품질 순서(예: referencepose명확성/환경준비적합성, resolution, stableID tie-break)로 대표를 고른다. minor-state 신호나 docking결과는 선정점수에 쓰지 않는다. 목표24개를채우려고추가기준을낮추지 않는다.
8. 최종표와 minor-state 잔기표, 재현 코드/명령·입력hash·카운트·남은한계를 새 실험폴더에 저장한다. 사용자에게는 실제 확정 수와 주요제외이유, minor-state양성잔기를 요약한다.

## 도구/원자료 위치

클러스터:
- Python `/gpfs/jaeohshin/work/envs/multistate-analysis/bin/python` (gemmi0.7.5,numpy,pandas,Bio; RDKit없음)
- MMseqs `/gpfs/deepfold/others/mmseqs/bin/mmseqs`
- qFit archive `/gpfs/jaeohshin/work/datasets/qfit_at_scale/qfit_PDBs.zip` (~9GB; 전체를 다시 다운로드하지 말 것)
- archive member `qfit_PDBs/{pdb_lower}_qFit_010.pdb`
- source rows `results/apo_holo/residues.csv.gz` (~1.9MB; columns group,holo,apo,chain,i,resname,pocket,ligand_class,ligands,cls,apo_mc,holo_mc,apo_apo_cls,res_holo,res_apo)
- replicate cache `reference/frozen/data/qfit_at_scale/replicate_cache/*.jsonl` (~284MB)
- split `runs/2026-10-19_scratch-v2_01/split/`
- source full sequence `data/pdb_seqres.txt.gz`
- previous homology/geometry artifacts `experiments/docking_pilot_20261028/run_v3/`
- CCD/deposition cache `experiments/docking_pilot_20261028/run_v1/download_cache/`와 `work/binding_site/cif/`

데스크탑 python3: gemmi0.6.5,numpy,RDKit2025.09.3; Bio없음. Gemmi 버전 차이에 주의.

Gemmi에서 확인한 API:
- `gemmi.make_assembly(assembly, model, gemmi.HowToNameCopiedChain...)`
- `st.assemblies`; `assembly.generators`; generator `chains`, `subchains`, `operators`; operator `transform` (mat/vec), `name`.
- `UnitCell.find_nearest_pbc_position(ref,pos,image_idx,inverse=False)`
- `NeighborSearch.find_atoms(pos,alt='\0',min_dist=0,radius=...)`
- `Transform.inverse()`
- https://gemmi.readthedocs.io/en/latest/analysis.html
- https://gemmi.readthedocs.io/en/latest/mol.html

아직 이 API로 full contact audit를 구현/검증하지 않았다. 특히 image_idx=0이어도 unit-cell translation으로 다른copy가 가까워질 수 있으므로 단순image_idx만으로crystalcontact를분류하지 않는다. 원chain이같아도 다른symmetrycopy일수있고, 다른chain이생물학적assembly일수있다.

## 해석 주의

- 79건은 최종실험적격수도, 79개독립단백질도 아니다.
- 4개minor-state양성case는 아직최종selection아니다.
- 8TAV를 논문제목만으로 모든ZKRinstance공유결합이라고단정하지 않는다.
- `|dres|≤0.3`은 기존source에서 crystallographic resolution 차이이며 backboneRMSD가 아니다.
- B/C에서같은apo backbone을쓰기때문에 apo–holo RMSD cutoff는 비교통제의필수조건이아니다. 다만이번입력79건에서는≤0.75 Å로고정했다.
- source_group과homology_component는다르다. homology graph의모든query를포함해isolates를빠뜨리지 않는다.
- 표준화된atommapping, protonation, cofactor환경, missingatoms, receptor예산은향후docking단계에서맞춰야한다. 오늘선정만으로엔진준비가검증된것은아니다.
- 기존PIPPack예정seed20260928, Bepoch88/Cepoch111, 같은architecture; 추론설정은실제실행전원기록과대조한다.
- 기존도킹설계5조건(nativeapo/nativeholo/B/C/Btemperaturecontrol)×8 receptor/search slots×3 repeats=case당120calls. 아직실행승인으로취급하지 않는다.

최종 산출물 권장: `case_review_all79.tsv`, `eligible_cases.tsv`, `pilot_manifest.tsv`, `minor_state_residues.tsv`, `assembly_contact_audit.json`, `RESULTS.md`, 재현스크립트 및 inputSHA256. 선정결과를확정하지못한항목은근거없이통과시키지말고명시적으로defer하되, 나머지완료가능한선정은끝까지수행한다.
