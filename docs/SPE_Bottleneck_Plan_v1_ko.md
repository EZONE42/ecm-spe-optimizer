# ECM 저항 분해로 고분자 전해질 병목을 다시 그린다

**프로젝트 계획서 — KG14기 이차전지 아카데미 ECM 프로젝트 「고분자 전해질 셀의 저항 병목 지도」** (작성 기준일 2026-10-28)

**결론부터 말하면, 출처가 있는 전도도 값으로 다시 계산해도 고분자 전해질 막의 벌크 저항(L/σ)은 셀 전체 면적비저항의 약 4\~15%에 그치며, 다섯 개 시나리오 어디에서도 1순위가 아니었다.** 1순위는 건식 PEO에서는 복합 양극 내부의 이온 경로였고, 인시츄(in-situ) 겔에서는 양극 쪽과 리튬 음극 쪽 계면이 비슷한 크기로 나란히 섰다. 다만 계산 과정에서 계획을 바꿀 만한 사실이 하나 드러났다. 리튬 이온 수송수 t+가 낮은 건식 PEO는 정상상태 농도분극 때문에 막 관련 직류 저항이 R\_b/t+까지 커진다. PEO의 t+는 0.06\~0.18로 보고되어 있으므로([Pesko et al. 2017](https://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf)) 25 °C에서는 막이 다시 1순위가 될 수 있다. 따라서 이 프로젝트의 정확한 명제는 \*\*"막의 이온전도도 σ는 병목이 아니다. 진짜 병목은 양극 이온 경로·계면, 그리고 막의 t+이다"\*\*이다. 이를 검증하려고 팀은 공개 원시 EIS 데이터(Khudyshkina 등, Zenodo, CC BY 4.0)를 MATLAB `lsqnonlin`으로 피팅해 R·C(CPE) 표를 만든다. 이어서 DRT(이완시간 분포)로 회로 요소 수를 정하고, 저항마다 Arrhenius/VTF 온도 의존성을 맞춰 메커니즘을 판별한다. 그 결과를 Ω·cm²로 환산해 셀 수준 병목 순위와 최소 필요 전도도 σ\_min을 역산하고, Simulink 2RC 셀 모델에서 개선 레버를 0.5C/1C 방전 곡선과 발열로 비교한다. 마지막으로 인시츄 경화 카솔라이트(양극용 전해질)의 침투율·전환율·잔류 모노머 trade-off를 설계안으로 제시한다. 실험실 측정 없이 2주 안에 끝낼 수 있고, 평가표의 "MATLAB fitting으로 R·C 값 출력"(10점)과 "과업 결과 도출"(20점)에 직접 대응한다. 결과는 중국 GB/T 43568-2026의 액체 함량 분류 같은 산업 쟁점과도 연결된다.

> **표기 규칙** — 이 문서의 수치에는 아래 꼬리표를 붙인다. 발표 자료에서도 같은 규칙을 지킨다.  
> 

| 꼬리표 | 의미 |
| :---- | :---- |
| **\[출처\]** | 링크된 1차 문헌·데이터에서 직접 확인한 값 |
| **\[계산\]** | \[출처\] 값으로 팀이 직접 계산한 값 |
| **\[가정\]** | 출처가 없어 팀이 정한 설계·추정값. 민감도 분석 대상 |
| **\[미검증\]** | 교과서·통념 수준이라 발표 전에 원문 확인이 필요한 값 |

---

## 0\. 한 장 요약 (두괄식)

| 항목 | 내용 |
| :---- | :---- |
| **주제** | 고분자 전해질 셀의 저항 병목 지도: ECM 저항 분해로 "병목 \= 막 이온전도도"라는 통념을 검증 |
| **핵심 주장(Thesis)** | 셀 수준 1순위 병목은 막의 σ가 아니다. 양극 내부 이온 경로·계면, 그리고 (건식 PEO에서는) 막의 낮은 t+로 인한 농도분극이다. 저항별 메커니즘에 근거해 조성·공정·셀 디자인 개선 우선순위를 낸다 |
| **가설** | H1 셀 수준 1순위 병목은 막(벌크 σ)이 아니다 · H2 병목은 온도·두께·침투율에 따라 이동한다 · H3 저항별 온도 의존성(Arrhenius/VTF)과 정전용량으로 메커니즘을 판별할 수 있다 |
| **H1이 기각되어도** | "막이 병목이 되는 조건(온도·두께·t+ 경계값)"을 정량화한 지도 자체가 결론이 되므로 프로젝트는 성립한다 |
| **사전 계산 결과** | 막 벌크 ASR 비중 3.6\~15% \[계산\]. 1순위는 양극 이온 경로(건식 PEO) 또는 양극·음극 계면(겔) \[계산·가정 혼합\]. t+가 0.19(25 °C) 또는 0.075(60 °C)보다 낮으면 막 관련 DC 저항이 1순위로 올라선다 \[계산\] |
| **데이터** | 공개 원시 EIS (Zenodo 7428858, CC BY 4.0) · 전도도 DB (Bradford/CheMixHub, 16,009행) · MD DB (polymernet, HTP-MD) · 문헌 추출값 · 합성 EIS |
| **도구** | MATLAB (`lsqnonlin`, `lsqcurvefit`, `ode45`), Simulink (2RC 셀), Python (DRT, Lin-KK, 선택 ML) |
| **모듈** | M1 이론·데이터 → M2 EIS 피팅 → M3 DRT → M4 Arrhenius/VTF → M5 Ω·cm² 병목 지도·σ\_min → M6 Simulink 레버 비교 → M7 카솔라이트 설계 (+M8 ML 조성 예측) |
| **산출물** | R·C(CPE) 표(신뢰구간 포함), 온도별 누적 ASR 막대 \= "병목 지도", 병목 교차온도, σ\_min·t+\_min 역산표, 0.5C/1C 방전·발열 비교, 토네이도 차트, 카솔라이트 설계 창(window) |
| **일정** | 2026-10-29 착수 → 11-11(수) 15:00 1차 제출 → 11-12(목) 15:00 최종 제출 (PPT/PDF) |

---

## 1\. 배경·문제의식·필요성 — 고분자는 다시 상용 경쟁에 들어왔고, 통념은 검증된 적이 없다

**용어 정의.** \*고분자 전해질(SPE, solid polymer electrolyte)\*은 리튬염을 고분자(대표적으로 PEO, 폴리에틸렌옥사이드)에 녹인 고체 전해질이다. \*겔 고분자 전해질(GPE)\*은 고분자 그물에 액체(가소제·용매)를 머금게 한 전해질이다. *인시츄(in-situ) 중합/고체화*는 액체 전구체(모노머·올리고머·개시제·염)를 셀 안에 주입한 뒤 열·빛·화학 개시로 굳히는 공정이다. \*카솔라이트(catholyte)\*는 복합 양극 내부 기공을 채우는 양극 전용 전해질을 말한다.

고분자·준고체 계열은 실험실 주제를 넘어 상용 경쟁 단계에 들어왔다. Blue Solutions는 2025-11-17 보도자료에서 "advanced polymer electrolytes"를 쓴 4세대(Gen4) 리튬금속 박막 전고체 셀을 **상온에서** 구동한다고 밝혔다. 집행관 입회 전동킥보드 시험에서 주행거리는 \*\*62.9 → 106 km(+68.8%)\*\*로 늘었고, 승용 모빌리티 적용 목표 시점은 \*\*"as of 2029"\*\*로 제시되었다([Blue Solutions PR](https://www.blue-solutions.com/wp-content/uploads/2025/11/pr-2025-07-blue-solutions-demonstrates-70-more-autonomy-with-its-breakthrough-gen4-lithium-metal-batteries.pdf)). Factorial의 FEST는 출처상 "액체/겔을 머금은 고체 매트릭스" 구조의 **준고체(quasi-solid)** 셀이다. 375 Wh/kg, 77 Ah, 600 사이클 이상, 상온 15→90% 충전 18분을 내세운다([Electronic Design](https://electronicdesign.com/markets/automotive/article/55291872/electronic-design-solid-state-battery-cells-for-fast-charging-high-energy-density-coming-in-2026)). 중국에서는 WeLion이 연 **11 GWh**를 가동하고 18 GWh를 증설 중이며([MIT Technology Review, 2026-10-06](https://technologyreview.com/2026/10/06/1144987/2026-climate-tech-companies-to-watch-welion-semi-solid-state-batteries)), 고체–액체 하이브리드 셀 시장은 2025년 약 15.45 GWh에서 2030년 346 GWh로 커질 것으로 추정된다([SMM](https://news.metal.com/newscontent/104120486-solid-state-battery-capacity-part-1-solid–liquid-batteries-a-transition-route-with-strong-supply-and-demand)). 국내 3사는 2024년 기사 기준 Samsung SDI 2027, LGES 2030, SK On 2029의 전고체 일정을 밝혔다. 그러나 **LGES의 고분자 전용 로드맵은 확인된 출처가 없다**([아시아경제](https://view.asiae.co.kr/en/article/2024030811265347893)).

규제 쪽 변화도 있다. 중국 **GB/T 43568-2026**은 2026-07-01 시행된 **권고(추천) 국가표준**이다. 전해액 질량 비율이 5% 미만이면서 120 °C 진공 6시간 건조 후 질량 손실이 0.5% 이하인 셀만 "전고체"로 부르고, 5\~20%는 "고체–액체 하이브리드", 20% 초과는 "액체 전지"로 분류한다([batterydesign.net](https://www.batterydesign.net/legislation-rules-and-regulations/gb-t-43568-2026/)). 인시츄 고체화 셀에서 미반응 모노머와 가소제는 성능을 올리는 수단이다. 동시에 이 기준에서는 "액체"로 계산될 수 있는 부담이기도 하다. M7의 설계 변수가 바로 이 지점에 걸린다.

**문제의식.** 업계와 교과서는 고분자 전해질의 약점을 흔히 "상온 이온전도도가 낮다"로 요약한다. 실제로 반결정 PEO의 상온 σ는 10⁻⁵~~10⁻⁸ S/cm로 인용된다([Chen et al. 2021](https://www.frontiersin.org/articles/10.3389/fchem.2020.592604/pdf)). 하지만 셀의 전압 손실은 막 하나가 아니라 직렬로 이어진 여러 저항의 합이다. 막은 20~~30 µm로 얇게 만들 수 있지만, 그 안의 전해질이 양극 기공의 굴곡 경로를 60 µm 넘게 지나야 할 수도 있다. 막이 실제로 셀 저항의 몇 %를 차지하는지, 온도가 바뀌면 순위가 어떻게 뒤집히는지를 **ECM 저항 분해 → 메커니즘 귀속 → 셀 수준 환산**의 한 흐름으로 보여 준 학생 프로젝트는 찾기 어렵다. 이것이 과업 창의성(1.1)과 목적성(1.2)의 근거다.

**필요성.** 병목을 잘못 짚으면 개발 자원이 엉뚱한 곳에 쓰인다. 막 σ를 10배 올리는 소재 개발은 수년이 걸린다. 그런데 막이 셀 저항의 5%라면 그 성과는 셀 출력에 거의 드러나지 않는다. 반대로 양극 이온 경로가 병목이라면 카솔라이트 침투율이나 양극 공극률 같은 **공정 레버**가 훨씬 싸고 빠른 해법이 된다. 이 계획서는 그 판단을 정량 근거로 내리는 방법을 다룬다.

---

## 2\. 목표·가설·성공 기준 — 판정 규칙을 미리 정해 두어야 결과가 설득력을 갖는다

### 2.1 목표

| 구분 | 목표 |
| :---- | :---- |
| 최종 목표 | 온도(25/60/80 °C)·막 두께(10~~50 µm)·카솔라이트 침투율(0.4~~1.0)에 따른 **셀 수준 저항 병목 지도**와 **개선 레버 우선순위표** 제시 |
| 세부 목표 1 | 공개 원시 EIS를 MATLAB으로 피팅해 저항·CPE·유효정전용량 표를 95% 신뢰구간과 함께 출력 (평가표 2번) |
| 세부 목표 2 | 저항별 활성화 에너지(Ea) 또는 VTF 파라미터(B, T0)로 메커니즘 판별 |
| 세부 목표 3 | Ω·cm² 환산 → 1순위 병목, 병목 교차온도, σ\_min·t+\_min 역산 |
| 세부 목표 4 | Simulink 2RC 셀 모델로 레버별 0.5C/1C 방전 용량·평균전압·발열 비교 |
| 세부 목표 5 | 인시츄 경화 카솔라이트 설계 창(침투율·전환율·잔류 모노머 trade-off) 제시 |

### 2.2 가설과 판정 규칙

| 가설 | 내용 | 채택 기준(정량) | 기각 시 결론 |
| :---- | :---- | :---- | :---- |
| **H1** | 셀 수준 1순위 병목은 막이 아니다 | 몬테카를로(가정값 범위 내 무작위 표본 ≥ 1,000회)에서 막 벌크 ASR이 1순위일 확률 \< 20%, 그리고 5개 기준 시나리오 중 4개 이상에서 막 비중 \< 25%. **고주파(EIS R\_b) 기준과 DC(R\_b/t+) 기준을 분리해 각각 판정** | "막이 1순위가 되는 경계(T, L, t+)" 지도를 제시 |
| **H2** | 병목은 온도·두께·침투율에 따라 이동한다 | 탐색 범위 안에서 1순위 항목이 최소 1회 바뀌고, 그 교차 조건(교차온도 T\*, 교차 두께 L\*, 교차 침투율 φ\*)을 식으로 계산할 수 있음 | 병목이 조건과 무관하다는 것 자체가 강한 결론(단일 레버 집중) |
| **H3** | 온도 의존성으로 메커니즘을 판별할 수 있다 | 저항별 Ea의 95% CI가 서로 겹치지 않거나(예: 벌크 VTF 곡률 유의 vs 계면 Arrhenius), Ea와 유효정전용량 두 지표가 같은 귀속을 지지 | 판별 불가 조건(온도 점 수 부족, 반원 중첩)을 명시하고 필요한 추가 측정을 제안 |

### 2.3 모듈별 성공 기준

| 모듈 | 정량 기준 |
| :---- | :---- |
| M2 피팅 | Lin-KK 잔차 \< 1~~2% (관행값, \[미검증\]) · χ²/ν가 잡음 분산 수준 · 각 파라미터 95% CI \< ±20~~30% · 유효정전용량이 귀속과 물리적으로 맞음 |
| M3 DRT | 피크가 λ 10배 범위와 모든 온도에서 유지될 때만 "실재 피크"로 인정 |
| M4 온도 | Arrhenius 직선성 R² \> 0.98, 또는 VTF 대비 AIC 차이 \> 2로 모델 선택 |
| M5 지도 | 모든 시나리오에 대해 저항 항목별 ASR, 비중, 1순위 확률(몬테카를로) 제시 |
| M6 셀 모델 | 합성 데이터 회수 시험에서 R0 초기 IR 강하 \= I·R0 (오차 \< 1%) · 레버별 1C 방전 용량 차이 정량 |
| M7 설계 | 침투율–전환율–잔류 모노머 trade-off를 식으로 표현하고 GB/T 43568 분류선과 겹쳐 표시 |

---

## 3\. 이론 — 하나의 직렬 저항 사슬로 묶는 여섯 가지 식

이 장의 논리 흐름은 다음과 같다. ① 셀 전압 손실은 직렬 저항의 합이며, 면적으로 정규화한 ASR로 비교한다. ② 각 저항은 고유한 물리 기구(벌크 이온 이동, 계면 전하전달, 다공 전극 이온 경로, 농도분극)를 가진다. ③ 각 기구는 온도 의존성과 특성 주파수·정전용량이 다르다. ④ 그래서 EIS로 분해하고 온도 의존성으로 귀속할 수 있다. ⑤ 귀속된 저항을 셀 조건으로 환산하면 병목 지도가 된다.

### 3.1 면적비저항(ASR)과 전압 손실 — 모든 비교의 공통 단위

ASR  \= R · A                          \[Ω·cm²\]

ASR\_bulk \= L / σ                      \[Ω·cm²\]

ΔV  \= J · ASR\_total                   \[V\]

J   \= (C-rate) · Q\_areal              \[A/cm²\]

R        : 측정 저항 \[Ω\]

A        : 전극(전해질) 면적 \[cm²\]

L        : 전해질 막 두께 \[cm\]   (1 µm \= 1×10⁻⁴ cm)

σ        : 이온전도도 \[S/cm\]

J        : 전류밀도 \[A/cm²\]

C-rate   : 방전 속도 \[h⁻¹\] (1C \= 1시간에 전부 방전)

Q\_areal  : 면적당 용량 \[Ah/cm²\]

ASR은 면적에 무관하므로 크기가 다른 셀끼리 비교할 수 있다. 이 계획의 기준 셀은 Q\_areal \= 3 mAh/cm²이다 \[가정\]. 따라서 0.5C에서 J \= 1.5 mA/cm², 1C에서 3 mA/cm²이다 \[계산\]. ASR 100 Ω·cm²인 셀을 1.5 mA/cm²로 쓰면 0.15 V가 손실된다.

**대칭셀 주의.** Li\|전해질\|Li 대칭셀에서 측정한 계면저항에는 **계면 두 개**가 들어 있다. 계면 하나의 값은 `ASR_int = R_int · A / 2`이다.

### 3.2 벌크 이온 이동: Arrhenius vs VTF — 사슬 운동에 묶였는가

PEO에서 Li⁺는 에테르 산소 4\~6개에 배위된 채, 비정질 영역의 고분자 사슬 분절 운동(segmental motion)을 따라 이동한다 \[미검증, 교과서 수준: Ratner & Shriver, [Chem. Rev. 1988](https://doi.org/10.1021/cr00083a006)\]. 사슬이 움직여야 이온이 움직이므로, 전도도는 단순한 열활성 도약(Arrhenius)이 아니라 자유부피 이론(VTF)을 따른다.

\[Arrhenius\]   σ(T) \= σ0 · exp( −Ea / (R·T) )

              선형화: ln σ \= ln σ0 − (Ea/R)·(1/T)      → 기울기 \= −Ea/R

              변형형: σ·T \= σ0' · exp(−Ea/(R·T))       (어느 형을 썼는지 반드시 명시)

\[VTF\]         σ(T) \= A · T^(−1/2) · exp( −B / (T − T0) )

              선형화(T0 고정): ln(σ·T^½) \= ln A − B · \[1/(T − T0)\]

              B \= Ea,VTF / R   (의사활성화에너지, 진짜 Ea 아님)

              T0 ≈ Tg − 50 K   (경험칙, \[미검증\])

\[겉보기 Ea\]   Ea,app(T) \= −R · d ln σ / d(1/T) \= R · \[ B·T²/(T − T0)² − T/2 \]

σ0, A     : 전지수 인자 \[S/cm\], \[S·cm⁻¹·K^½\] — 캐리어 수와 시도 빈도

Ea        : 활성화 에너지 \[J/mol\] (1 eV \= 96.485 kJ/mol)

R         : 기체상수 8.314 J/(mol·K)   (kB \= 8.617×10⁻⁵ eV/K 를 쓰면 Ea 단위는 eV)

T         : 절대온도 \[K\]

B         : VTF 상수 \[K\]

T0        : 이상 유리전이온도(분절 운동이 멈추는 온도) \[K\]

Tg        : 유리전이온도 \[K\] (DSC로 측정)

**판별법.** ln σ–1/T 그래프에서 국소 기울기로 계산한 Ea,app가 온도와 무관하게 일정하면(±10%) Arrhenius다. 저온으로 갈수록 커지면 VTF(사슬 결합 이동)다. Tm(PEO 약 65 °C) 근처에서 계단처럼 꺾이면 결정 용융 "knee"다. (PEO)16LiTf의 DSC 용융 피크는 약 65 °C, 결정화도는 50.8\~54.9%였다. σ는 30 °C에서 7.6×10⁻⁷(가열)과 2.2×10⁻⁶ S/cm(냉각)로 열이력에 따라 약 3배 달랐다 **\[출처\]**([Chen et al. 2021](https://www.frontiersin.org/articles/10.3389/fchem.2020.592604/pdf)). 그러므로 Tm 위·아래를 하나의 직선으로 강제 피팅하면 안 되고, 가열·냉각 이력을 반드시 기록한다.

**WLF와의 관계.** WLF 식 `log aT = −C1(T − Tref)/(C2 + T − Tref)`는 VTF와 수학적으로 같다. B \= 2.303·C1·C2, T0 \= Tref − C2의 관계가 있고, Tref \= Tg에서 "보편상수" C2 ≈ 51.6 K이므로 T0 ≈ Tg − 50 K라는 경험칙이 나온다 \[미검증: [Williams, Landel, Ferry 1955](https://doi.org/10.1021/ja01619a008)\].

### 3.3 수송수 t+와 농도분극 — "보이지 않는 막 저항"

\*수송수 t+\*는 전체 전류 중 Li⁺가 나르는 비율이다. PEO-LiTFSI에서는 전류 대부분을 음이온(TFSI⁻)이 나른다. 그러면 방전 중 막 양쪽에 염 농도 구배가 생기고, 이것이 추가 전압 손실(농도분극)이 된다.

\[정의\]          t+ \= σ+ / (σ+ \+ σ−)

\[NMR 근사\]      t+,NMR \= D\_Li / (D\_Li \+ D\_F)

\[Bruce–Vincent\] t+,SS \= I\_ss · (ΔV − I\_0·R\_i,0) / \[ I\_0 · (ΔV − I\_ss·R\_i,ss) \]

                I\_0 \= ΔV / (R\_i,0 \+ R\_b,0)        (측정 초기전류 대신 계산값 권장)

                이상극한(계면저항 변화 없음): t+,SS \= I\_ss / I\_0

ΔV          : 대칭셀에 거는 일정 분극 전압 \[V\] (보통 \~10 mV; ΔV ≪ RT/F ≈ 25.7 mV @25 °C)

I\_0, I\_ss   : 초기 전류, 정상상태 전류 \[A\]

R\_i,0, R\_i,ss : 분극 전·후 EIS 계면저항 \[Ω\]

R\_b,0       : 분극 전 벌크저항 \[Ω\]

D\_Li, D\_F   : ⁷Li, ¹⁹F PFG-NMR 자기확산계수 \[cm²/s\]

식은 [Pesko et al. 2017, Eq. 1·7·8](https://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf)에서 가져왔다. 같은 논문의 90 °C PEO(5 kg/mol)-LiTFSI 실측값은 다음과 같다 **\[출처\]**. t+,SS는 r \= \[Li\]/\[EO\] 0.01→0.16에서 **0.18→0.06**으로 줄었고, t+,NMR은 0.17~~0.25였다. 염 확산계수 D는 5.8×10⁻⁸~~1.3×10⁻⁷ cm²/s, 전도도 최대값은 r \= 0.08(EO:Li \= 12.5)에서 **2×10⁻³ S/cm**였다. 농축용액 이론으로 구한 진짜 t+는 고농도에서 **음수**까지 내려갔다. 즉 Bruce–Vincent 값은 엄밀히는 "전류 분율"이며 진짜 수송수와 같지 않다([arXiv 2203.01659](https://arxiv.org/pdf/2203.01659)).

**핵심 유도: 농도분극 저항 R\_diff.** 이상극한 식 t+ \= I\_ss/I\_0에 I\_0 \= ΔV/R\_b, I\_ss \= ΔV/(R\_b \+ R\_diff)를 넣으면 다음이 나온다.

t+ \= R\_b / (R\_b \+ R\_diff)   ⇒   R\_diff \= R\_b · (1 − t+) / t+

R\_DC,membrane \= R\_b \+ R\_diff \= R\_b / t+

R\_b    : 막 벌크저항(고주파 EIS에서 보이는 값) \[Ω·cm²\]

R\_diff : 정상상태 농도분극 저항(저주파·DC에서만 보이는 값) \[Ω·cm²\]

t+ \= 0.1이면 DC에서 막은 EIS에 보이는 R\_b의 **10배** 저항으로 작동한다 \[계산\]. 이 저항은 고주파 EIS에서는 보이지 않으므로 "보이지 않는 막 저항"이라 부를 수 있다. 정상상태에 도달하는 시간은 대략 τ ≈ L²/(π²D)이다. L \= 30 µm, D \= 1×10⁻⁷ cm²/s이면 약 9초다 \[계산\]. 0.5C 방전(약 2시간) 동안은 사실상 계속 정상상태이므로 셀 모델에 반드시 넣어야 한다.

**한계전류.** 이상 이원 전해질의 1차원 정상상태 한계전류밀도는 다음과 같다(계수 표기는 \[미검증\], [Chazalviel 1990](https://doi.org/10.1103/PhysRevA.42.7355)).

J\_lim \= 2 · z · c0 · F · D\_salt / \[ (1 − t+) · L \]

c0     : 초기 염 농도 \[mol/cm³\]

F      : 패러데이 상수 96,485 C/mol

D\_salt : 앰비폴라 염 확산계수 \[cm²/s\]

z      : 전하수 (Li⁺ \= 1\)

D \= 1×10⁻⁷ cm²/s, c0 ≈ 1.5×10⁻³ mol/cm³ \[가정\], t+ \= 0.1, L \= 50 µm를 넣으면 J\_lim ≈ 6.4 mA/cm²이다 \[계산\]. 0.5C 기준 J \= 1.5 mA/cm²보다는 높지만, 60 °C 이하에서는 D가 더 작으므로 여유가 줄어든다.

**Li⁺ 전도도 지표.** 비교에는 σ 대신 σ·t+를 쓴다. 90 °C PEO-LiTFSI는 σ \= 2×10⁻³ S/cm, t+ ≈ 0.1로 σ\_Li ≈ 2×10⁻⁴ S/cm이다. 상온 PDOL 겔은 σ \= 3.62×10⁻⁴ S/cm, t+ \= 0.61로 σ\_Li ≈ 2.2×10⁻⁴ S/cm이다 \[계산; 원자료 [Pesko 2017](https://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf), [Xin et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11173723)\]. 상온 겔의 Li⁺ 전도도가 90 °C PEO와 비슷하다는 뜻이다.

### 3.4 Nernst–Einstein과 이온성(ionicity) — 캐리어 수의 문제인가, 이동도의 문제인가

σ\_NE \= (F² / (R·T)) · ( c+·z+²·D+ \+ c−·z−²·D− )

Ionicity \= σ\_EIS / σ\_NE  (≤ 1\)        Haven ratio H\_R \= σ\_NE / σ\_EIS

c±   : 양·음이온 몰농도 \[mol/cm³\]

D±   : 자기확산계수 \[cm²/s\] (PFG-NMR 또는 MD)

MD에서의 D: MSD(t) \= ⟨\|r(t) − r(0)\|²⟩ \= 6·D·t  (log MSD–log t 기울기 1 구간에서만)

Ionicity가 1보다 작으면 이온쌍이나 클러스터가 같이 움직여 전하를 덜 나른다는 뜻이다(표준식, \[미검증\]). 이온 간 상관을 모두 포함한 일반식(Onsager 계수)과 기준틀 문제는 [arXiv 2203.01659](https://arxiv.org/pdf/2203.01659)가 정리했다. 이 프로젝트에서는 Khudyshkina 데이터의 PFG-NMR(29\~73 °C)과 EIS σ를 비교해 이 지표를 계산한다(M4 보조 분석).

### 3.5 계면 저항 — 국소 장벽이라 Ea가 크고 Arrhenius를 따른다

계면에서는 Li⁺가 탈용매화(desolvation)하고, 전하를 주고받고(charge transfer), 서로 다른 상의 경계를 건너야 한다. 이런 과정은 사슬 운동보다 국소 에너지 장벽에 지배된다. 그래서 대체로 Arrhenius를 따르고 Ea가 벌크보다 크다. 출처로 확인한 근거는 다음과 같다.

| 근거 | 값 | 출처 |
| :---- | :---- | :---- |
| PEO-LiCF3SO3 / LLT 세라믹 계면 Li⁺ 전달 Ea | **97.6 kJ/mol (≈1.01 eV)**, LLT 벌크 28.3·입계 36.4 kJ/mol 대비 3배 이상 | [Abe, Iriyama, Ogumi 2004](https://www.mrs-j.org/pub/tmrsj/vol29_no1/vol29_no1_009.pdf) |
| 고분자/세라믹 계면저항 | **1.2 kΩ·cm² @30 °C**, 3층 셀 저항의 66%, DMC 첨가 후 사실상 0 | [Chen et al., ACS Energy Lett. 2019 (ORNL)](https://impact.ornl.gov/en/publications/determining-and-minimizing-resistance-for-ion-transport-at-the-po/) |
| 복합 전극 이온저항 Ea | **≈ 20 kJ/mol**, 전해질 이온저항과 비슷해 전하전달(Ea 더 큼)과 구분 가능 | [Fukutsuka et al., JES 2017](https://doi.org/10.1149%2F2.1381714jes) |
| 인시츄 PDOL 겔 Li\|Li 계면 | 33.4 Ω(조립 직후) → 53.6 Ω에서 안정. 액체 전해질은 61.1 → 212.3 Ω (7일) | [Xin et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11173723) |

PDOL 논문은 셀 면적을 밝히지 않았다. CR2032의 리튬 전극(지름 약 16 mm, A ≈ 2.0 cm²)으로 \[가정\]하고 대칭셀 계면 두 개로 나누면 계면 하나당 **약 33\~54 Ω·cm²**이다 \[계산·가정\]. PEO 산화로 생기는 저항성 CEI(양극 계면막)는 고전압 NMC에서 문제가 된다고 알려져 있다. 다만 산화 개시 전위 3.8\~4.0 V vs Li/Li⁺는 \[미검증\]이다.

### 3.6 복합 양극의 이온 경로 — Bruggeman과 전송선 모델(TLM)

복합 양극 안의 전해질은 활물질과 도전재 사이를 굽이굽이 돌아간다. 이 때문에 유효 전도도가 크게 줄어든다.

\[MacMullin\]    N\_M \= σ0 / σ\_eff \= τ / ε

\[Bruggeman\]    τ \= ε^(1 − α)   ⇒   σ\_eff \= σ0 · ε^α      (구형 입자 이상값 α \= 1.5)

\[침투율 반영\]  ε\_eff \= ε · φ\_inf                          (\[가정\] 모델: 기공 중 전해질이 채운 분율)

\[이온저항\]     R\_ion·A \= d / σ\_eff \= d · τ / (σ0 · ε)   \[Ω·cm²\]

ε       : 양극 공극률(전해질이 차지할 수 있는 부피분율) \[-\]

φ\_inf   : 카솔라이트 침투율(기공 중 실제로 채워진 분율) \[-\]

τ       : 굴곡도 인자 \[-\]

α       : Bruggeman 지수 \[-\] (1.5보다 크면 경로가 더 나쁨)

d       : 양극 두께 \[cm\]

σ0      : 양극 내부 전해질 자체 전도도 \[S/cm\] (막의 σ와 다를 수 있음)

전송선 모델(TLM)은 다공 전극을 "이온 저항 사다리 \+ 계면 임피던스 가로대"로 본다. 차단(blocking) 조건 임피던스는 다음과 같다([BioLogic AN-BAT-011](https://icpms.cz/labrulez-bucket-strapi-h3hsga3/AN_BAT_011_3da1cad8af/AN-BAT-011.pdf), 방법 원전 [Landesfeind et al. 2016](https://doi.org/10.1149/2.1141607jes)).

Z\_TLM(ω) \= sqrt( R\_ion / (Q·(jω)^α) ) · coth( sqrt( R\_ion · Q·(jω)^α ) )

같은 구조에서 가로대를 전하전달 저항 R\_ct(양극 계면 전체를 대표하는 면적 ASR)로 바꾸면 직류 극한이 나온다. 이 식은 팀이 유도한 것이다.

R\_cathode,DC \= sqrt(R\_ion · R\_ct) · coth( sqrt(R\_ion / R\_ct) )

극한 1 (R\_ct ≫ R\_ion, 반응 균일):   coth x ≈ 1/x \+ x/3  ⇒  R\_cathode ≈ R\_ct \+ R\_ion/3

극한 2 (R\_ct ≪ R\_ion, 반응 국재):   R\_cathode ≈ sqrt(R\_ion · R\_ct)

반응 침투 깊이: λ / d ≈ sqrt(R\_ct / R\_ion)   (λ/d ≪ 1이면 양극 대부분이 반응에 참여 못 함)

R\_ion/3은 차단 조건 Nyquist의 저주파 실수축 오프셋으로 알려져 있으며, 위 coth 전개에서 수학적으로 나온다(이번 조사에서 이를 명시한 원문은 찾지 못함). **사전 계산에서 확인한 함정이 있다.** 이전 계산처럼 항상 R\_ion/3을 쓰면 극한 2 조건(건식 PEO 25 °C)에서 양극 저항을 약 3배 과대평가한다. 이 계획서에서는 coth 일반식을 쓴다.

### 3.7 EIS 요소와 정전용량 — "어느 반원이 무엇인가"를 정하는 규칙

\*EIS(전기화학 임피던스 분광법)\*는 작은 교류 신호를 주파수별로 걸어 임피던스 Z(ω)를 재는 방법이다. \*ECM(등가회로모델)\*은 이 스펙트럼을 R, C, CPE, Warburg 같은 요소의 조합으로 표현한다.

저항            Z \= R

축전기          Z \= 1 / (jωC)

인덕터          Z \= jωL                              (\>100 kHz 케이블 잡음)

CPE             Z \= 1 / ( Q·(jω)^n ),  0 \< n ≤ 1     (n=1 이상 축전기, n=0.5 Warburg 유사)

R‖CPE           Z \= R / (1 \+ R·Q·(jω)^n),   f\_peak \= 1 / (2π·(R·Q)^(1/n))

반무한 Warburg   Z \= A\_W / sqrt(jω)                   (45° 직선)

유한 Warburg(투과형, Ws)  Z \= R\_W · tanh( sqrt(jωτ\_W) ) / sqrt(jωτ\_W)   (대칭셀 농도분극)

유한 Warburg(반사형, Wo)  Z \= R\_W · coth( sqrt(jωτ\_W) ) / sqrt(jωτ\_W)

\[CPE → 유효정전용량\]

 Hsu–Mansfeld(막 두께 방향 분포): C\_eff \= R^((1−n)/n) · Q^(1/n)

 Brug(표면 분포, 이중층):          C\_eff \= ( R\_s·R\_p/(R\_s+R\_p) )^((1−n)/n) · Q^(1/n)

Q       : CPE 계수 \[F·s^(n−1)\] — 정전용량이 아님\!

n       : CPE 지수 \[-\]

R\_s, R\_p: 직렬저항, 병렬저항 \[Ω\]

τ\_W     : 확산 시간상수 \[s\]

CPE 변환식의 출처는 [Zahner IM7 문서](https://doc.zahner.de/im7/en/applications/constant_phase_element/index.html)이며, 원전은 Brug 1984, Hsu & Mansfeld 2001, Hirschorn 2010이다. 측정된 C\_eff 크기로 반원의 정체를 정한다. Irvine–Sinclair–West(1990) 분류는 원문 표를 확인하지 못해 \[미검증\]이지만, 대략 다음과 같다.

| C\_eff (면적·두께 정규화 전) | 귀속 |
| :---- | :---- |
| \~10⁻¹² \~ 10⁻¹¹ F | 벌크 (기하 정전용량) |
| \~10⁻¹¹ \~ 10⁻⁸ F | 입계·표면층·고분자/세라믹 계면 |
| \~10⁻⁹ \~ 10⁻⁷ F | SEI/CEI 같은 표면 피막 |
| \~10⁻⁷ \~ 10⁻⁵ F | 전극–전해질 이중층·전하전달 |
| ≳10⁻⁴ F, 45° 꼬리 | 확산·농도분극 |

**벌크 반원 위치를 계산으로 예측하기(팀 유도).** 막의 기하 정전용량은 C\_g \= ε0·εr·A/L이다. 벌크 반원의 꼭짓점 주파수는 형상과 무관하게 다음으로 정해진다.

f\_bulk \= 1 / (2π · R\_b · C\_g) \= σ / (2π · ε0 · εr)

ε0 \= 8.854×10⁻¹⁴ F/cm,  εr(PEO) ≈ 5\~10 \[미검증\]

예) A \= 1 cm², L \= 100 µm, εr \= 5  →  C\_g ≈ 44 pF

    σ \= 1×10⁻⁵ S/cm → f\_bulk ≈ 3.6 MHz;  σ \= 1×10⁻³ S/cm → f\_bulk ≈ 360 MHz

> 정정: 연구 노트의 "C\_g ≈ 4~~9 pF"는 10배 작게 계산된 값이다. 올바른 값은 ≈ 44~~89 pF이다 \[계산\].

σ가 10⁻⁵ S/cm만 넘어도 벌크 반원은 일반 측정창(≤1 MHz) 밖으로 나간다. 따라서 고온·겔 스펙트럼에서는 **R\_b를 고주파 실수축 절편으로 읽어야** 한다. Abe 등도 353 K에서 PEO 벌크 반원이 측정창 밖에 있어 절편만 썼다([Abe et al. 2004](https://www.mrs-j.org/pub/tmrsj/vol29_no1/vol29_no1_009.pdf)).

### 3.8 DRT — 반원이 몇 개인지 데이터가 말하게 한다

\*DRT(distribution of relaxation times, 이완시간 분포)\*는 임피던스를 무수한 RC 요소의 연속 분포로 풀어, 시간상수 τ마다 얼마의 저항이 있는지(γ)를 보여 준다.

Z(ω) \= jωL \+ R\_∞ \+ ∫ γ(ln τ) / (1 \+ jωτ) d(ln τ)

이산화·정규화:  min\_γ≥0  \|\| W(A·γ − Z) \|\|² \+ λ · \|\| L\_d·γ \|\|²

피크 k:  f\_k \= 1/(2π τ\_k),   R\_k \= ∫\_peak γ d ln τ,   C\_k \= τ\_k / R\_k

γ     : 분포 함수 \[Ω\]

λ     : Tikhonov 정규화 강도 \[-\]  (GCV/mGCV로 선택)

L\_d   : 1차 또는 2차 미분 행렬

W     : 가중 행렬 (1/\|Z\|)

Ciucci 그룹의 비교 연구에서는 λ 선택법 가운데 **GCV와 modified GCV가 가장 정확**했다([ChemRxiv 2022](https://chemrxiv.org/engage/chemrxiv/article-details/636fbfc5fbfd38d985dbf1bb)). 도구로는 pyDRTtools(Python GUI)와 DRTtools(MATLAB GUI)가 있다([github.com/ciuccislab](https://github.com/ciuccislab)). Li\|SPE\|Li에서 예상되는 순서는 τ가 짧은 쪽부터 벌크(고온에서는 R\_∞), SEI/계면(10⁴~~10² Hz), 전하전달(10²~~1 Hz), 확산(\<1 Hz, 넓은 꼬리)이다. 확산 꼬리와 τ창 가장자리 피크는 RC 요소로 세지 않는다.

### 3.9 2RC 셀 모델과 발열 — 저항 지도를 방전 곡선으로 번역

V\_t   \= OCV(SOC) − I·R0 − V1 − V2

dV1/dt \= −V1/(R1·C1) \+ I/C1

dV2/dt \= −V2/(R2·C2) \+ I/C2

dSOC/dt \= −I / (3600 · Q\_Ah)

Q\_gen \= I²·R0 \+ V1²/R1 \+ V2²/R2       (+ 선택: 가역열 I·T·dOCV/dT)

\[열\]  m·c\_p·dT/dt \= Q\_gen − h·A\_s·(T − T\_amb)

\[온도 연동\] R\_i(T) \= R\_i,ref · exp\[ (Ea\_i/R)·(1/T − 1/T\_ref) \]

\[셀 환산\]   R\_i,cell \= ASR\_i / A\_cell,   C\_i,cell \= C\_i,areal · A\_cell   (τ\_i \= R\_i·C\_i 불변)

V\_t    : 단자전압 \[V\]          I : 방전 전류 \[A\] (양수)

R0     : 순간 저항(막 R\_b \+ 접촉 \+ 양극 이온경로의 빠른 성분) \[Ω\]

R1,C1  : 계면(ms\~s 시간상수) \[Ω\], \[F\]

R2,C2  : 확산·농도분극(10\~100 s 시간상수) \[Ω\], \[F\]

Q\_Ah   : 셀 용량 \[Ah\]          A\_cell : 총 전극 면적 \[cm²\]

m, c\_p : 셀 질량 \[kg\], 비열 \[J/(kg·K)\]   h·A\_s : 열전달 계수×면적 \[W/K\]

식은 표준 Thevenin 2RC다([MathWorks Simscape Battery `ecm`](https://www.mathworks.com/help/simscape-battery/ref/simscape.battery.parameters.ecm.html) 참고. 2RC 식 자체를 쓴 MathWorks 페이지는 확인하지 못함). R\_i(T) 연동식이 M4에서 얻은 Ea를 M6 셀 모델로 넘기는 다리 역할을 한다.

### 3.10 병목 교차온도 — H2를 식으로 만든다 (팀 유도)

두 저항 i, j의 Ea가 다르면, 어느 온도 T\*에서 크기가 같아진다.

1/T\* \= 1/T\_ref − R · ln( R\_i,ref / R\_j,ref ) / ( Ea\_i − Ea\_j )

R\_i,ref, R\_j,ref : 기준온도 T\_ref에서의 ASR \[Ω·cm²\]

Ea\_i, Ea\_j       : 각 저항의 활성화 에너지 \[J/mol\]

Ea가 큰 저항(계면, 반결정 벌크)은 저온에서 더 빨리 커진다. 그래서 T\*보다 낮은 온도에서는 Ea가 큰 쪽이 병목이 된다. VTF 저항에는 국소 Ea,app(T)를 넣어 반복 계산한다.

---

## 4\. 저항 → 메커니즘 → 진단 지문 → 개선 레버 매핑 — 각 저항이 남기는 "지문"

이 표가 프로젝트의 지도 범례이며, 팀 시트 `8_저항×메커니즘` 탭과 일대일로 대응시킨다.

| 저항 항목 | 지배 메커니즘 | EIS 지문 (주파수·C\_eff) | 온도 지문 | 기타 지문 | 개선 레버 (근거) |
| :---- | :---- | :---- | :---- | :---- | :---- |
| **R\_b 막 벌크** | 분절 운동 결합 Li⁺ 이동, 결정화, 이온쌍 | 고주파 반원(pF) 또는 실수축 절편 | VTF 곡률. PEO는 Tm 근처 knee | DSC Tg, Xc, 가열·냉각 이력 | 박막화(8.6 µm PI/PEO 막, 2.3×10⁻⁴ S/cm @30 °C, [Wan 2019](https://www.nature.com/articles/s41565-019-0465-3)), 가소제·저Tg 망(VEC 추가로 0.16→1.37 mS/cm, Ea 0.43→0.26 eV, [Tang 2023](https://www.nature.com/articles/s41467-023-37997-6)) |
| **R\_diff 막 농도분극** | 낮은 t+ (음이온이 전류를 나름) | 저주파 Ws(유한 Warburg) | 염 확산계수 D(T)를 따름 | Bruce–Vincent t+, 전류 의존 | 단일이온 전도체(t+ \> 0.85, 1.3×10⁻⁵ S/cm @60 °C, \[미검증\] [Bouchet 2013](https://doi.org/10.1038/nmat3602)), 고 t+ 인시츄 겔(VC/EGDA+HFBMA t+ 0.82, [RSC 2025](https://pubs.rsc.org/en/content/articlehtml/2025/qm/d5qm00116a)) |
| **R\_int,a 음극 계면(Li\|SPE)** | SEI, 전하전달, 접촉 손실 | 중주파 반원(nF~~µF), n ≈ 0.7~~0.9 \[미검증\] | Arrhenius, 벌크보다 큰 Ea (예상) | 휴지시간에 따라 증가, 압력 의존 | 인시츄 겔 젖음(PDOL 계면 7일 안정, [Xin 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11173723)), 이중층 구조(음극 쪽 환원 안정층) |
| **R\_int,c 양극 계면·CEI** | 탈용매화·전하전달, PEO 산화 피막 | 중\~저주파 반원(µF) | Arrhenius, Ea 큼 | 상한전압 4.0 V vs 4.3 V 비교 시 증가 | 불소화 가교제(ESW 4.55→5.08 V, CEI 2→4 nm, [Tang 2023](https://www.nature.com/articles/s41467-023-37997-6)), 양극 쪽 산화안정층 |
| **R\_ion,c 양극 내부 이온 경로** | 낮은 σ0 × 굴곡도 × 미충전 기공 | 고주파 45° 직선(TLM), 오프셋 R\_ion/3 | 전해질 σ와 같은 Ea(≈20 kJ/mol 예, [Fukutsuka 2017](https://doi.org/10.1149%2F2.1381714jes)) | 양극 두께·공극률·침투율에 비례 | **인시츄 경화 카솔라이트**(침투율↑), 전구체 침투 후 경화(KRICT E-ICP: 200사이클 유지율 75% vs 22%, [Newswise 2026](https://www.newswise.com/articles/view/850751/?sc=rssn)) |
| **R\_pc 고분자/세라믹 계면** (복합 전해질일 때) | 이종상 Li⁺ 전달 장벽 | 별도 반원(10⁻¹¹\~10⁻⁸ F) | Ea ≈ 1 eV ([Abe 2004](https://www.mrs-j.org/pub/tmrsj/vol29_no1/vol29_no1_009.pdf)) | 필러 함량에 비례 | 연속 3D·정렬 세라믹 망, 계면 개질(DMC로 1.2 kΩ·cm² → ≈0, [Chen 2019](https://impact.ornl.gov/en/publications/determining-and-minimizing-resistance-for-ion-transport-at-the-po/)) |

**Ea 기반 판별 작업 규칙** (가설 H3의 판정 규칙. 수치 경계는 \[출처\] 예시 \+ \[미검증\] 통념으로 만든 작업 가설)

| 관측 | 귀속 |
| :---- | :---- |
| Ea \< 약 0.3 eV, Arrhenius 직선 | 겔·액체 유사 경로 (PDOL 0.16 eV \[출처\]) |
| ln σ–1/T 곡률 유의, VTF 적합 우수 | 분절 운동 결합 벌크 |
| Ea ≳ 0.5\~1 eV, Arrhenius, C\_eff ≥ nF | 계면 전달 (세라믹 계면 1.0 eV \[출처\]) |
| Tm 아래 Ea 큼 \+ 열이력 | 반결정 PEO 벌크. 계면과 혼동 주의 → C\_eff로 재확인 |

---

## 5\. 데이터 계획 — 무엇을, 어디서, 어디에 쓰는가

### 5.1 데이터 계획 표

| ID | 데이터명 | 종류 | 형식·규모 | 라이선스 | 적용 대상 (모듈 / 저항) | 적용 방식 | 한계 | 대체안 |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| **D1** | Khudyshkina et al. 2023 원시 EIS (BPE-LiTFSI/NaTFSI/KTFSI, 20:1·15:1·10:1) | 실측 EIS | txt 9개(각 ≈332 kB, 다수 스펙트럼 추정) \+ EIS\_README | **CC BY 4.0** ([Zenodo 7428858](https://zenodo.org/records/7428858)) | M2·M3·M4 / R\_b, (대칭셀이면) R\_int | 고주파 절편·R‖CPE 피팅 → σ \= L/(R\_b·A) → σ(T) VTF/Arrhenius → 염 농도 효과 | 온도 스윕 여부·셀 형상·열 형식 \[미검증\] (README 확인 필요). 블록공중합체이므로 순수 PEO와 다름 | D8 디지타이징, D9 합성 EIS |
| **D2** | 같은 기록의 Bruce–Vincent 세트 (T+\_EIS 분극 전후 \+ T+\_CA) | 실측 EIS \+ 정전위 전류 | txt (EIS 약 1.4 kB, CA 약 351 kB) | CC BY 4.0 | M2·M4·M5 / R\_diff, R\_int | Pesko Eq.8로 t+ 계산 → R\_diff \= R\_b(1−t+)/t+ | 15:1 한 조성뿐 | 문헌 t+ (Pesko 0.06\~0.18) |
| **D3** | 같은 기록의 PFG-NMR (⁷Li, ¹⁹F; 29\~73 °C) | 확산계수 | txt | CC BY 4.0 | M4 / R\_b 메커니즘 | Nernst–Einstein σ\_NE, ionicity, t+,NMR 교차검증 | 15:1만 | 문헌값 |
| **D4** | 같은 기록의 DSC | 열분석 | txt | CC BY 4.0 | M4 / R\_b | Tg → VTF의 T0 \= Tg − 50 K 고정 | — | 문헌 Tg |
| **D5** | Bradford et al. 2023 고분자 전해질 전도도 DB (CheMixHub 경유) | 문헌 집계 σ(T) | CSV, 원자료 16,009행·218 DOI / 정제 11,350행 | 데이터 CC BY 4.0, 코드 MIT ([GitHub chemixhub](https://github.com/chemcognition-lab/chemixhub), [arXiv 2506.12231](https://arxiv.org/pdf/2506.12231)) | M1·M5·M8 / R\_b | 고Mw PEO-LiTFSI 기준값(앵커) 추출, M8 학습 데이터 | 문헌 라벨 잡음(100배 단위 오류 사례, [Wang 2026](https://arxiv.org/pdf/2609.01621)). 계면·양극 정보 없음 | 원 논문 그림 대조 |
| **D6** | polymernet (Xie et al. 2022\) MD 전도도 | 시뮬레이션 σ | CSV, 5 ns 7,123행 / 50 ns 117행 | **LICENSE 파일 없음** → 연구 인용용으로만 ([GitHub](https://github.com/txie-93/polymernet)) | M8 / R\_b | 고분자 골격 상대 순위, 실험 앵커로 보정 | 단위 S/cm 추정 \[미검증\]. LiTFSI만, 고정 조건 | D7 |
| **D7** | HTP-MD (TRI) | MD σ, D\_Li, D\_TFSI, t+ | 웹 CSV 내보내기, ≈6,057 고분자, 353 K, 1.5 mol/kg | 라이선스 명시 없음 ([arXiv 2208.01692](https://arxiv.org/pdf/2208.01692), [htpmd.matr.io](https://www.htpmd.matr.io)) | M8 / R\_b, R\_diff | σ·t+ 순위로 후보 선별 | 80 °C 한 조건. 절대값 편향 | D6 |
| **D8** | 문헌 추출값 (Pesko 2017, Xin 2024, Tang 2023, Abe 2004, Chen 2019, Fukutsuka 2017, RSC 2025 리뷰) | 수치 | 팀 시트 `15_문헌수집양식` | 각 논문 인용 | M4·M5·M7 / 전 항목 | 시나리오 파라미터·Ea 기준값·전환율 | 셀 면적 미기재(PDOL) 등 맥락 손실 | 범위로 처리 → 몬테카를로 |
| **D9** | SI Nyquist 그래프 디지타이징 | 그림 → 수치 | CSV (Z′, Z″) | 인용(공정 이용) | M2 / 전체 셀·인시츄 셀 | 절편·반원 지름만 피팅 | 주파수 정보 손실 → CPE n 신뢰 불가 | D10 |
| **D10** | 합성 EIS (팀 생성) | 모델 생성 | MATLAB/CSV | 팀 소유 | M2·M3 검증, M5 민감도 | 알려진 파라미터 \+ 1% 잡음 → 피팅 회수 시험 | 실측 아님 | — |
| **D11** | 기준 셀 설계값 | 설계 파라미터 | 팀 시트 `9_ASR계산기` | — | M5·M6 | 양극 60 µm, ε 0.3, 3 mAh/cm², α 1.5, A\_cell | 전부 \[가정\] | 민감도·몬테카를로 |
| **D12** | OCV–SOC 표 | 셀 특성 곡선 | 표 | 교재 `ECM_공개데이터셋` 목록에서 선택 (팀 확인 필요) | M6 | 2RC 모델의 OCV(SOC) 조회표 | 고분자 셀 고유 OCV 아님 | NMC/LFP 대표 곡선 |
| (후보) | KIT RDR 고분자 대칭셀 임피던스 / ChemRxiv 자동화 PEO 플랫폼 보충자료 | 실측 EIS | 미확인 | 미확인 ([TIB](https://service.tib.eu/ldmservice/dataset/rdr-doi-10-35097-rhsdekozhadetmgk), [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/675cccadf9980725cf041455)) | M2 | 접근되면 D1 보강 | 봇 차단으로 내용 확인 못 함 | — |

### 5.2 출처로 확인한 앵커값 — 이전 가정값을 대체한다

D5에서 고Mw(≥100 kg/mol) PEO, LiTFSI 단독, 필러·제2고분자 없음으로 거르면 147행·13 DOI가 남는다. EO:Li는 1000/(44.05·몰랄농도)로 환산했다 **\[출처: CheMixHub CSV 로컬 분석\]**([GitHub](https://github.com/chemcognition-lab/chemixhub)).

| 온도 | EO:Li \= 20 중앙값 (S/cm) | 범위·n | 이전 계획 가정값 | 판정 |
| :---- | :---- | :---- | :---- | :---- |
| 25 °C | 1.16×10⁻⁵, 8.44×10⁻⁵ (n=2) | 열이력 의존 큼 | 1×10⁻⁵ | 하한 근처 → **유지(보수적)** |
| 30 °C | 9.4×10⁻⁶ \~ 1.3×10⁻⁴ (n=3) | — | — | 범위로 보고 |
| 60 °C | **4.8×10⁻⁴** | 4.4\~7.1×10⁻⁴, n=4 | 3×10⁻⁴ | **교체** |
| 80 °C | **1.03×10⁻³** | 0.89\~1.28×10⁻³, n=4 | — | 신규 |
| 겔 25 °C | PDOL 3.62×10⁻⁴ ([Xin 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11173723)) · 불소화 UV망 1.37×10⁻³ ([Tang 2023](https://www.nature.com/articles/s41467-023-37997-6)) · LiPF6/DOL 3.8×10⁻³ ([RSC 2025](https://pubs.rsc.org/en/content/articlehtml/2025/qm/d5qm00116a)) | — | 1×10⁻³ | 범위 3.6×10⁻⁴ \~ 3.8×10⁻³ |

60~~80 °C에서 전도도 최대는 EO:Li 12~~20에 있었다(80 °C: EO:Li 12에서 1.48×10⁻³ S/cm). 이는 Pesko의 90 °C 최적 r \= 0.08(EO:Li 12.5)과 맞아 앵커가 내부적으로 일관된다 \[계산\].

### 5.3 데이터 취급 규칙

모든 수치는 `15_문헌수집양식` 탭에 기록한다. 기록 항목은 DOI, 그림 번호, 단위, 온도, 열이력, 두께, 전극 종류, 꼬리표(\[출처\]/\[계산\]/\[가정\]/\[미검증\])다. Wang 등은 문헌 데이터의 오류를 네 유형(본문–그림 불일치, 주석 모호성, 단위 불일치, 측정 맥락 누락)으로 분류했다. 겔 전해질 논문 83편 중 10.8%에서 본문–그림 불일치가 있었고, 100배 오류 하나가 ML 데이터베이스까지 전파된 사례도 있었다([arXiv 2609.01621](https://arxiv.org/pdf/2609.01621)). 그러므로 앵커값은 원 논문 그림과 한 번 더 대조한다. 디지타이징은 WebPlotDigitizer 같은 도구로 점마다 2\~3회 반복하고 평균·표준편차를 남긴다(도구 버전·라이선스는 \[미검증\]).

---

## 6\. 단계별 방법 M1–M8 — 입력에서 검증까지

모듈 사이의 데이터 흐름은 다음과 같다.

M1(이론·데이터) ──► M2(EIS 피팅: R, Q, n, C\_eff ± CI) ◄── M3(DRT: 요소 수)

                         │

                         ▼

                 M4(R\_k(T) → Arrhenius/VTF → Ea\_k, 메커니즘)

                         │

                         ▼

                 M5(Ω·cm² 환산 → 병목 지도, T\*, σ\_min, t+\_min, 몬테카를로)

                         │

                         ▼

                 M6(Simulink 2RC: 레버별 방전·발열, 토네이도)

                         │

                         ▼

                 M7(카솔라이트 설계 창) ◄── M8(선택: ML 새 조성 σ 예측)

### M1. 이론·데이터 수집 (10-29 \~ 11-01)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | 이 계획서 3장, D1\~D12 |
| 절차 | ① Zenodo 7428858 "Download all" → EIS\_README.txt, T+\_EIS\_CA\_README.txt, 논문 실험절에서 두께·면적·온도·측정기 형식 확보 ② CheMixHub CSV 필터링(5.2절 조건) → 앵커표 ③ 문헌 추출값 7편을 `15_문헌수집양식`에 입력 ④ 교재 `ECM_공개데이터셋` 목록에서 OCV–SOC 표 확보 |
| 출력 | 데이터 인벤토리 표, 앵커값 표, 파일 형식 메모 |
| 검증 기준 | 모든 수치에 꼬리표와 DOI가 있음. 앵커값이 원문 그림과 일치 |
| 리스크 | Zenodo README에 온도 정보 없음 → 논문 본문에서 확보. 그래도 없으면 D9·D10으로 온도 시리즈 구성 |

\# M1: CheMixHub에서 PEO-LiTFSI 앵커 추출 (열 이름은 실제 CSV 헤더로 확인 후 수정)

import pandas as pd

df \= pd.read\_csv('datasets/polymer-electrolyte/raw\_data/PolymerElectrolyteData.csv')

peo \= df\[(df\['S1 SMILES'\]=='\[Cu\]OCC\[Au\]') & (df\['S1 Mw'\].fillna(df\['S1 Mn'\])\>=1e5)

         & (df\['Salt1'\]=='LiTFSI') & df\['Salt2'\].isna() & df\['S2 SMILES'\].isna()

         & df\['Inorganic filler'\].isna()\]

peo\['EO\_Li'\] \= 1000/(44.05\*peo\['Salt1 Molality'\])

anchor \= (peo.assign(EOr=peo\['EO\_Li'\].round())

             .groupby(\['Temperature','EOr'\])\['Conductivity'\]

             .agg(\['median','min','max','count'\]))

### M2. MATLAB EIS 등가회로 피팅 (11-01 \~ 11-05) — 평가표 2번의 핵심

| 항목 | 내용 |
| :---- | :---- |
| 입력 | D1·D2 스펙트럼(f, Z′, Z″), D9, D10 |
| 절차 | ① 데이터 정리: 고→저 주파수 정렬, 유도성(Z″\>0 고주파) 점 제거 ② Lin-KK로 잔차 \> 1\~2% 점 제거 ([KIT Lin-KK](https://www.iam.kit.edu/et/english/Lin-KK.php) 또는 impedance.py `linKK`) ③ M3 DRT로 요소 수 결정 ④ Nyquist에서 초기값 추정: R\_s \= 고주파 절편, R\_1 \= 반원 지름, Q\_1 ≈ 1/(R\_1·(2πf\_peak)^n), n0 \= 0.85 ⑤ `lsqnonlin`으로 복소 비선형 최소제곱(실수·허수 잔차 쌓기, 모듈러스 가중 1/\|Z\|, 양수 파라미터 log10 변환) ⑥ Jacobian으로 표준오차·95% CI 계산 ⑦ C\_eff, f\_peak 계산 → 3.7절 표로 귀속 ⑧ ASR \= R·A (대칭셀 계면은 /2) |
| 회로 | 차단셀 SS\|SPE\|SS: `L–R_b–CPE_dl` (고온) 또는 `R_s–(R_b‖CPE_g)–CPE_dl` (저온) · 대칭셀 Li\|SPE\|Li: `R_b–(R_int‖CPE_int)–Ws` · 전체셀: `R_0–(R_int,a‖CPE)–(R_ct,c‖CPE_dl)–W` (+양극 TLM 선택) |
| 출력 | **R·Q·n·C\_eff·f\_peak·ASR ± 95% CI 표**, Nyquist·Bode 피팅 그림, 잔차 그림(주파수별 %) |
| 검증 기준 | 잔차 무작위·\<2%, χ²/ν ≈ 잡음분산, CI \< ±20\~30%, D10 합성 데이터 회수 시험에서 참값이 95% CI 안 |
| 리스크 | 반원 중첩으로 상관이 큰 파라미터 → DRT로 초기값, 요소 축소, CI 보고. 측정창 밖 벌크 반원 → 절편법 |

% M2: 대칭셀 Li\|SPE\|Li 피팅  Rs \+ (R1\|\|CPE1) \+ (R2\|\|CPE2) \+ Ws

% 데이터: f \[Hz\], Z \= Zre \+ 1j\*Zim (열벡터, 고→저 주파수)

w   \= 2\*pi\*f;

cpe \= @(Q,n,w) 1./(Q.\*(1j\*w).^n);

par \= @(Za,Zb) Za.\*Zb./(Za+Zb);

Zws \= @(Rw,tw,w) Rw.\*tanh(sqrt(1j\*w\*tw))./sqrt(1j\*w\*tw);

model \= @(p,w) p(1) \+ par(p(2),cpe(p(3),p(4),w)) ...

                    \+ par(p(5),cpe(p(6),p(7),w)) \+ Zws(p(8),p(9),w);

islog \= logical(\[1 1 1 0 1 1 0 1 1\]);           % n1, n2만 선형

toP \= @(x) (islog.\*10.^x \+ \~islog.\*x);

wt  \= 1./abs(Z);                                 % 모듈러스 가중

res \= @(x) \[real(model(toP(x),w)-Z).\*wt; imag(model(toP(x),w)-Z).\*wt\];

p0  \= \[5 200 1e-9 0.9 800 1e-6 0.8 300 10\];      % Nyquist에서 읽은 초기값

x0  \= p0; x0(islog) \= log10(p0(islog));

lb  \= x0-3; ub \= x0+3; lb(\~islog)=\[0.5 0.5\]; ub(\~islog)=\[1 1\];

opt \= optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',2e4, ...

                   'FunctionTolerance',1e-12,'StepTolerance',1e-12);

\[x,resnorm,r,\~,\~,\~,J\] \= lsqnonlin(res,x0,lb,ub,opt);

p   \= toP(x);

nu  \= numel(r)-numel(x);  s2 \= resnorm/nu;               % chi2/nu

Cx  \= s2\*inv(full(J'\*J));  sx \= sqrt(diag(Cx))';         % 피팅좌표 표준오차

sp  \= sx; sp(islog) \= p(islog).\*log(10).\*sx(islog);      % 선형 파라미터로 전파

CI95 \= 1.96\*sp;

C1 \= (p(2)^(1-p(4))\*p(3))^(1/p(4));                      % Hsu–Mansfeld C\_eff

C2 \= (p(5)^(1-p(7))\*p(6))^(1/p(7));

fpk1 \= 1/(2\*pi\*(p(2)\*p(3))^(1/p(4)));

A \= 2.0;                                                 % \[cm²\] 실제 셀 면적으로 교체

ASR\_int1 \= p(2)\*A/2;  ASR\_int2 \= p(5)\*A/2;               % 대칭셀: 계면 1개당

T \= table({'Rs';'R1';'Q1';'n1';'R2';'Q2';'n2';'Rw';'tw'}, p', CI95', ...

          'VariableNames',{'Param','Value','CI95'});

disp(T); fprintf('C1=%.2e F, C2=%.2e F, chi2/nu=%.2e\\n',C1,C2,s2);

% Optimization Toolbox가 없으면 fminsearch(@(x) sum(res(x).^2), x0) 사용

% M2 검증: 합성 데이터 회수 시험 (D10)

ptrue \= \[5 200 1e-9 0.9 800 1e-6 0.8 300 10\];

f \= logspace(6,-2,70)'; w \= 2\*pi\*f;

Zs \= model(ptrue,w); Z \= Zs.\*(1 \+ 0.01\*(randn(size(Zs))+1j\*randn(size(Zs))));

% 위 피팅 실행 후 all(abs(p-ptrue) \< CI95) 이면 통과

### M3. DRT로 회로 요소 수 결정 (11-02 \~ 11-04)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | Lin-KK 통과 스펙트럼 |
| 절차 | ① λ를 10⁻⁵\~10⁻¹로 스윕하며 DRT 계산 ② pyDRTtools GUI에서 GCV/mGCV로 λ 선택(설치법·입력 형식은 \[미검증\], README 확인) 또는 아래 자체 코드 ③ λ 10배 범위·모든 온도에서 유지되는 피크만 "실재 피크"로 인정 ④ 피크 수 \= R‖CPE 수, 저주파 꼬리 있으면 Warburg 추가 ⑤ 피크 f\_k(T) 추적 → ln f\_k vs 1/T 기울기가 다르면 Ea가 다른 별개 과정 |
| 출력 | 온도별 DRT 그림, 피크 표(τ\_k, f\_k, R\_k, C\_k \= τ\_k/R\_k), 채택 회로 |
| 검증 기준 | DRT 피크 수와 ECM 요소 수 일치, 피크 C\_k와 ECM C\_eff 귀속 일치. R 값은 DRT 면적이 아니라 ECM 피팅값으로 보고 |
| 리스크 | λ가 작으면 가짜 진동 피크, 크면 피크 병합 → λ 스윕 결과 첨부. 확산 꼬리·τ창 가장자리 피크는 세지 않음 |

\# M3: 최소 Tikhonov-DRT (비음수 최소제곱, 2차 미분 정규화). pyDRTtools 결과와 교차 확인용

import numpy as np

from scipy.optimize import nnls

def drt(f, Z, lam=1e-3, M=80, tmin=None, tmax=None):

    w \= 2\*np.pi\*f

    tmin \= tmin or 0.1/w.max(); tmax \= tmax or 10/w.min()

    tau \= np.logspace(np.log10(tmin), np.log10(tmax), M)

    dl  \= np.log(tau\[1\]/tau\[0\])

    WT  \= np.outer(w, tau)

    Are \= dl/(1+WT\*\*2);  Aim \= \-WT\*dl/(1+WT\*\*2)

    K   \= np.block(\[\[np.ones((len(f),1)), Are\],\[np.zeros((len(f),1)), Aim\]\])  \# \[R\_inf, γ\]

    z   \= np.concatenate(\[Z.real, Z.imag\])

    W   \= np.concatenate(\[1/np.abs(Z)\]\*2)                                    \# 모듈러스 가중

    L2  \= np.diff(np.eye(M), 2, axis=0); L2 \= np.hstack(\[np.zeros((M-2,1)), L2\])

    Kaug \= np.vstack(\[K\*W\[:,None\], np.sqrt(lam)\*L2\])

    zaug \= np.concatenate(\[z\*W, np.zeros(M-2)\])

    x, \_ \= nnls(Kaug, zaug)

    return tau, x\[1:\], x\[0\]          \# tau \[s\], gamma \[Ω\], R\_inf \[Ω\]

for lam in \[1e-5, 1e-4, 1e-3, 1e-2, 1e-1\]:

    tau, g, Rinf \= drt(f, Z, lam)    \# 피크 위치가 λ에 따라 유지되는지 비교

### M4. VTF/Arrhenius 메커니즘 귀속 (11-04 \~ 11-06)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | 온도별 R\_b, R\_int, R\_ct (M2), Tg (D4), D\_Li·D\_F (D3) |
| 절차 | ① σ \= L/(R\_b·A) ② ln σ vs 1000/T 그림, 인접점 Ea,app(T) 계산 ③ Tm(≈60\~65 °C) 위·아래 구간 분리 ④ Arrhenius(ln(σT) vs 1/T 직선)와 VTF(3변수 비선형, 점 6개 미만이면 T0 \= Tg − 50 K 고정)를 모두 피팅 → AIC 비교 ⑤ 계면 저항은 ln(1/R\_k) 또는 ln(T/R\_k) vs 1/T로 Ea\_k ⑥ 4장 판별 규칙 적용(Ea \+ C\_eff 이중 확인) ⑦ 보조: σ\_NE와 ionicity 계산 |
| 출력 | 저항별 Ea ± CI 또는 (A, B, T0) 표, Arrhenius 그림(저항별 색), 메커니즘 귀속표 |
| 검증 기준 | Ea CI 비중첩 여부, AIC 차이 \> 2, 귀속이 C\_eff 판정과 일치 |
| 리스크 | 온도 점 부족 → T0 고정 / Arrhenius만 보고. 반결정 벌크 Ea(≈1 eV, \[미검증\])와 계면 Ea 혼동 → 구간 분리 \+ C\_eff |

% M4: Arrhenius / VTF 피팅 (K 단위 온도, S/cm 단위 전도도)

T   \= \[303 313 323 333 343 353\]';  sig \= L./(Rb\*A);

Rg  \= 8.314;

pA  \= polyfit(1./T, log(sig.\*T), 1);  Ea\_kJ \= \-pA(1)\*Rg/1000;     % σT 형 Arrhenius

vtf \= @(b,T) b(1).\*T.^-0.5.\*exp(-b(2)./(T-b(3)));

\[b,rn,res,\~,\~,\~,Jv\] \= lsqcurvefit(vtf,\[10 1000 170\],T,sig,\[0 100 100\],\[1e6 5000 250\]);

EaVTF\_kJ \= b(2)\*Rg/1000;

% T0 고정판 (Tg는 D4 DSC; PEO Tg ≈ 213 K → T0 ≈ 163 K, \[미검증\])

T0 \= Tg \- 50;

pV \= polyfit(1./(T-T0), log(sig.\*sqrt(T)), 1);  B \= \-pV(1);

% 국소 겉보기 Ea

Ea\_app \= \-Rg\*diff(log(sig))./diff(1./T)/1000;                     % kJ/mol

% 모델 선택 (AIC, 잔차 정규분포 가정)

n \= numel(T);

aic \= @(rss,k) n\*log(rss/n) \+ 2\*k;

rssA \= sum((log(sig) \- (polyval(pA,1./T) \- log(T))).^2);

rssV \= sum((log(sig) \- log(vtf(b,T))).^2);

fprintf('AIC Arr=%.2f  VTF=%.2f\\n', aic(rssA,2), aic(rssV,3));

% 계면 저항 R\_k에도 동일: pk \= polyfit(1./T, log(1./Rk), 1); Ea\_k \= \-pk(1)\*Rg/1000;

### M5. Ω·cm² 셀 수준 환산 · 1순위 병목 · σ\_min 역산 (11-05 \~ 11-07)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | M2·M4 결과, D8 문헌값, D11 기준 셀 |
| 절차 | ① 저항 항목별 ASR(T) 생성: R\_b \= L/σ(T), R\_diff \= R\_b(1−t+)/t+, R\_cathode \= TLM coth 일반식, R\_int,a ② 시나리오별 누적 막대 \= **병목 지도** ③ **고주파 기준(R\_b만)과 DC 기준(R\_b/t+) 두 장을 따로** 작성 ④ 교차온도 T\*(3.10절), 반응 침투 깊이 λ/d ⑤ σ\_min·t+\_min 역산 ⑥ \[가정\] 파라미터를 범위 안에서 log-균등 표본 추출하는 몬테카를로 → 1순위 확률 |
| 출력 | 병목 지도 그림, 1순위 확률표, σ\_min/t+\_min 역산표, 교차 조건 표 (`9_ASR계산기` 탭 갱신) |
| 검증 기준 | ASR 합 × J \= 셀 모델 초기 과전압과 일치. 몬테카를로 결과가 표본 수 2배에서 ±2%p 이내로 수렴 |
| 리스크 | 계면값 출처 부족 → 범위로 처리하고 결과를 "조건부"로 표시 |

**사전 계산 결과 (출처 σ 반영, TLM 일반식 사용)** — 양극 60 µm, ε 0.3, α 1.5, 3 mAh/cm², 0.5C(1.5 mA/cm²) \[가정\]

| 시나리오 | σ (S/cm) | L (µm) | φ\_inf | R\_b | R\_cathode (TLM) | R\_int,a | 합계 | 막 비중 | 1순위 | λ/d | 0.5 V 손실 허용 J |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| B 건식 PEO 25 °C | 1.16×10⁻⁵ \[출처\] | 30 | 0.4 \[가정\] | 258.6 | 1,366 (R\_int,c 150 \[가정\]) | 60 \[가정\] | 1,685 | **15.4%** | 양극(81%) | 0.11 | 0.30 mA/cm² (≈C/10) |
| C 건식 PEO 60 °C | 4.8×10⁻⁴ \[출처\] | 30 | 0.6 | 6.2 | 83.8 (40) | 15 | 105 | **5.9%** | 양극(80%) | 0.49 | 4.8 |
| D 건식 PEO 80 °C | 1.03×10⁻³ \[출처\] | 30 | 0.6 | 2.9 | 62.7 (40) | 15 | 81 | **3.6%** | 양극(78%) | 0.72 | 6.2 |
| E PDOL 겔 25 °C | 3.62×10⁻⁴ \[출처\] | 20 | 0.9 | 5.5 | 55.8 (25) | 54 \[출처·가정 면적\] | 115 | **4.8%** | 양극 ≈ 음극 계면(48% / 47%) | 0.46 | 4.3 |
| F 고전도 겔 25 °C | 1×10⁻³ \[가정 원안\] | 20 | 0.9 | 2.0 | 37.9 (25) | 10 | 50 | **4.0%** | 양극(76%) | 0.76 | 10.0 |

단위는 모두 Ω·cm²(J는 mA/cm²)이며 \[계산\]이다. 이전 원안(R\_ion/3 고정, σ 가정값)의 결과는 막 비중 3.9\~6.6%, 25 °C 건식 PEO 양극 4,811 Ω·cm²였다. 출처 σ로 바꿔도 결론의 방향은 같다. 하지만 TLM 일반식을 쓰면 25 °C 건식 PEO의 양극 저항이 약 1/3로 줄고 막 비중이 15%로 올라간다. 대신 반응 침투 깊이가 양극 두께의 11%(약 7 µm)에 그친다. 즉 **저온 건식 PEO의 병목은 전압 손실보다 "양극 대부분이 쓰이지 못하는 용량 손실"로 나타난다**.

**DC 기준 — t+가 판을 바꾼다** \[계산; t+ 범위는 [Pesko 2017](https://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf)의 90 °C 값 0.06\~0.18과 PDOL 0.61\]

| 시나리오 | t+ 경계값 t+\* \= R\_b / R\_cathode | 막 DC(R\_b/t+) 비중 @ t+ \= 0.06 | @ 0.15 | @ 0.61 |
| :---- | :---- | :---- | :---- | :---- |
| B 건식 PEO 25 °C | **0.19** | 75% (1순위) | 55% (1순위) | 23% |
| C 건식 PEO 60 °C | **0.075** | 51% (1순위) | 30% | 9% |
| D 건식 PEO 80 °C | 0.046 | 38% | 20% | 6% |
| E PDOL 겔 25 °C | 0.099 | 46% | 25% | 8% (실제 t+ \= 0.61) |
| F 고전도 겔 25 °C | 0.053 | 41% | 22% | 6% |

t+가 t+\*보다 낮으면 막 관련 DC 저항이 양극을 넘어선다. 건식 PEO는 25 °C에서 t+\*가 0.19라서 보고된 t+ 범위 대부분이 경계값 아래에 있다. 따라서 **막 σ 기준으로는 H1이 채택되지만, 막 t+ 기준으로는 25 °C 건식 PEO에서 H1이 기각될 수 있다.** 이 이분법이 발표의 핵심 메시지가 된다. 겔(t+ 0.61)은 어느 기준으로도 막이 병목이 아니다. 단, 양극 내부에도 농도분극이 있는데 이 모델에는 들어 있지 않다. R\_diff 식도 이상 이원 전해질 가정에 기반한다. 그러므로 이 표는 "막이 병목이 될 수 있는 상한"으로 해석한다.

**σ\_min 역산표** \[계산\]

막:      σ\_min \= L / ASR\_target

양극:    σ0,min \= d / ( 3 · ASR\_target · (ε·φ\_inf)^1.5 )     (R\_ct ≫ R\_ion 극한)

t+:      t+\_min \= R\_b / R\_cathode   (막 DC 저항 ≤ 양극 저항이 되는 조건)

| 막 두께 L | ASR\_target 5 | 10 | 20 Ω·cm² |
| :---- | :---- | :---- | :---- |
| 20 µm | 4×10⁻⁴ | 2×10⁻⁴ | 1×10⁻⁴ S/cm |
| 30 µm | 6×10⁻⁴ | 3×10⁻⁴ | 1.5×10⁻⁴ |
| 50 µm | 1×10⁻³ | 5×10⁻⁴ | 2.5×10⁻⁴ |

| 침투율 φ\_inf (ε \= 0.3, d \= 60 µm) | 양극 R\_ion/3 ≤ 10 Ω·cm² | ≤ 20 Ω·cm² |
| :---- | :---- | :---- |
| 0.4 | 4.8×10⁻³ S/cm | 2.4×10⁻³ |
| 0.6 | 2.6×10⁻³ | 1.3×10⁻³ |
| 0.9 | 1.4×10⁻³ | 7.1×10⁻⁴ |
| 1.0 | 1.2×10⁻³ | 6.1×10⁻⁴ |

두 표를 나란히 놓으면 결론이 선명해진다. 같은 셀에서 **양극 안 전해질에 요구되는 σ는 막에 요구되는 σ보다 5\~20배 높다**. 실용 목표 ASR(예: "막 ≤ 10 Ω·cm²")은 출처를 확보하지 못했으므로 \[미검증\]이며, 발표에서는 5/10/20 Ω·cm² 세 기준을 모두 보여 준다.

% M5: 병목 지도 \+ 몬테카를로 (가정 파라미터는 log-균등)

d=60e-4; eps=0.3; J=1.5e-3; N=5000;

lu \= @(a,b,n) 10.^(log10(a)+(log10(b)-log10(a))\*rand(n,1));

sig \= 4.8e-4\*ones(N,1); L \= 30e-4;                 % 시나리오 C (σ는 \[출처\])

phi \= 0.4+0.5\*rand(N,1);  Rct \= lu(15,150,N);  Ra \= lu(5,100,N);  tp \= 0.06+0.12\*rand(N,1);

Rb  \= L./sig;  Rion \= d./(sig.\*(eps\*phi).^1.5);

Rcat \= sqrt(Rion.\*Rct)./tanh(sqrt(Rion./Rct));

M\_hf \= \[Rb Rcat Ra\];          M\_dc \= \[Rb./tp Rcat Ra\];

\[\~,i1\] \= max(M\_hf,\[\],2);  \[\~,i2\] \= max(M\_dc,\[\],2);

P\_hf \= histcounts(i1,0.5:1:3.5)/N;   P\_dc \= histcounts(i2,0.5:1:3.5)/N;

fprintf('1순위 확률 HF \[막 양극 음극\] \= %s\\n', mat2str(P\_hf,2));

fprintf('1순위 확률 DC \[막 양극 음극\] \= %s\\n', mat2str(P\_dc,2));

% 누적 막대: bar(categorical(scen), \[Rb Rcat Ra\], 'stacked')

### M6. Simulink 2RC 셀 모델로 레버 비교 (11-06 \~ 11-09)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | M5 ASR 항목, M4 Ea, D12 OCV 표, A\_cell·Q\_Ah \[가정\] |
| 매핑 | R0 \= R\_b \+ 접촉 \+ 양극 TLM의 빠른 성분 · R1·C1 \= 계면(τ ms\~s) · R2·C2 \= 농도분극(R2 ≈ R\_diff, τ2 ≈ τ\_W 또는 L²/(π²D)). **이 매핑은 모델링 선택이므로 보고서에 명시** |
| 레버 시나리오 | 기준 / 막 σ ×10 / 막 두께 ½ / t+ 0.15→0.6 / 침투율 0.4→0.9 / 양극 계면 R\_int,c ½ / 음극 계면 ½ / 운전온도 25→60 °C |
| 절차 | ① MATLAB `ode45`로 먼저 구현·검증 ② 같은 식을 Simulink 블록으로 구현 ③ 0.5C·1C CC 방전, 차단전압까지 ④ `Simulink.SimulationInput` 배열 \+ `sim`/`parsim`으로 시나리오 일괄 실행 ⑤ ±50% 단일 변수 민감도 → 토네이도 |
| 출력 | 0.5C/1C 방전 곡선 겹침 그림, 레버별 방전 용량·평균전압·IR 강하·총 발열\[Wh\]·최대 발열\[W\] 표, 토네이도 차트 |
| 검증 기준 | t \= 0⁺에서 전압 강하 \= I·R0 (오차 \< 1%). ode45와 Simulink 결과 일치(전압 RMS 차 \< 1 mV). 저전류 극한에서 용량 → Q\_Ah |
| 리스크 | Simscape Battery `ecm`/`fitECM`은 R2025a 이상 필요 → 아카데미 라이선스 확인([MathWorks](https://www.mathworks.com/help/simscape-battery/ref/simscape.battery.parameters.ecm.html)). 없으면 기본 Simulink 블록만 사용 |

% M6: 2RC ODE (Simulink 구현 전 기준 해)

P.Q=2.0; P.A=100;  Vcut=3.0;                   % Ah, cm² (총 전극면적), V  \[가정\]

ASR=\[R0a R1a R2a\]; Ca=\[C1a C2a\];               % Ω·cm², F/cm²  (M5에서)

R=ASR/P.A; C=Ca\*P.A;

ocv=@(s) interp1(soc\_tab,ocv\_tab,s,'linear','extrap');

for Crate=\[0.5 1\]

  I=Crate\*P.Q;

  f=@(t,x)\[ \-I/(3600\*P.Q); \-x(2)/(R(2)\*C(1))+I/C(1); \-x(3)/(R(3)\*C(2))+I/C(2) \];

  ev=@(t,x) deal(ocv(x(1))-I\*R(1)-x(2)-x(3)-Vcut,1,-1);

  \[t,x\]=ode45(f,\[0 3600/Crate\*1.2\],\[1 0 0\],odeset('Events',ev,'MaxStep',5));

  V=ocv(x(:,1))-I\*R(1)-x(:,2)-x(:,3);

  Qgen=I^2\*R(1)+x(:,2).^2/R(2)+x(:,3).^2/R(3);  % W

  cap=I\*t(end)/3600; E=trapz(t,V\*I)/3600;       % Ah, Wh

  Eloss=trapz(t,Qgen)/3600;                     % Wh 발열

end

**Simulink 조립 순서 (Simscape 불필요)**

1. `Constant`(I) → `Gain`(−1/(3600·Q)) → `Integrator`(초기값 1\) \= SOC  
2. SOC → `1-D Lookup Table`(OCV)  
3. RC 가지 2개: `Transfer Fcn` R/(R·C·s \+ 1)에 입력 I → V1, V2  
4. `Sum`: V\_t \= OCV − I·R0 − V1 − V2  
5. `Compare To Constant`(V\_t \< Vcut) → `Stop Simulation`  
6. `MATLAB Function`: Q\_gen \= I²R0 \+ V1²/R1 \+ V2²/R2 → `Integrator` \= 누적 발열  
7. (선택) 열 블록: m·c\_p·dT/dt를 `Integrator`로 구현하고, `MATLAB Function`으로 R\_i(T)를 계산해 되먹임  
8. 스크립트로 시나리오 실행: `in(k)=Simulink.SimulationInput('cell2RC'); in(k)=in(k).setVariable('R0',...)` → `out=sim(in)`

% M6: 토네이도 (1C 방전 용량 변화 기준)

names \= {'R\_b','R\_diff','R\_int,a','R\_int,c','R\_ion,c'};

base \= runCell(p0);  rows \= zeros(numel(names),2);

for k=1:numel(names)

  lo=p0; lo.(matlab.lang.makeValidName(names{k}))=0.5\*p0.(matlab.lang.makeValidName(names{k}));

  hi=p0; hi.(matlab.lang.makeValidName(names{k}))=1.5\*p0.(matlab.lang.makeValidName(names{k}));

  rows(k,:)=\[runCell(lo)-base, runCell(hi)-base\];

end

\[\~,ix\]=sort(abs(rows(:,2)-rows(:,1))); 

barh(rows(ix,:),'stacked'); yticklabels(names(ix)); xline(0,'k');

xlabel('\\Delta 1C 방전용량 \[Ah\]'); title('레버 민감도 (±50%)');

% runCell: 위 ODE를 감싼 함수(파라미터 구조체 → 1C 용량 반환)

### M7. 해결안 설계 — 인시츄 경화 카솔라이트 (11-07 \~ 11-09)

M5에서 1순위가 양극 이온 경로·계면으로 나오면 레버는 **카솔라이트**다. 공정 개념은 다음과 같다. 양극 슬러리(또는 건조된 양극)에 저점도 모노머/올리고머 \+ 리튬염 \+ 개시제 전구체를 섞거나 함침한다. 이어서 열 또는 양이온 개시로 기공 안에서 굳힌다. 어둡고 두꺼운 양극에서는 UV가 깊이 들어가지 못하므로, 양극 내부 경화에는 열개시(AIBN 약 40\~60 °C)나 상온 양이온 개환중합이 현실적이다(공학적 추론). UV 경화는 자립막에 적합하다.

**trade-off의 구조** (식은 팀 설계 모델 \[가정\], 수치 근거는 \[출처\])

\[침투\]    φ\_inf \= f(점도 η, 기공 크기, 함침 시간·압력)       η가 낮을수록 φ\_inf ↑

          → ε\_eff \= ε·φ\_inf ↑ → R\_ion,c ∝ 1/(ε·φ\_inf)^1.5 ↓

          → R\_int,c ≈ R\_int,c0 / φ\_inf  (\[가정\]: 젖은 활물질 면적에 반비례)

\[전환\]    w\_res \= w\_mon · (1 − X)                            잔류 모노머 질량분율

\[액체분\]  f\_liq,cell \= ( m\_plast \+ m\_mon·(1 − X) ) / m\_cell   (GB/T 43568 판정 대상)

\[수축\]    ΔV/V ∝ X · (모노머 종류별 수축계수)                 \[미검증\] 개환 \< 아크릴

φ\_inf : 침투율 \[-\]      η : 전구체 점도 \[Pa·s\]      X : 모노머 전환율 \[-\]

w\_mon : 전구체 중 모노머 질량분율 \[-\]               m\_plast : 가소제 질량 \[g\]

m\_cell: 셀 질량 \[g\]

**레버별 근거 수치** **\[출처: [RSC Mater. Chem. Front. 2025 리뷰](https://pubs.rsc.org/en/content/articlehtml/2025/qm/d5qm00116a)\]**

| 계 | 개시 방식 | 전환율 X | σ (RT, mS/cm) | t+ | ESW (V) |
| :---- | :---- | :---- | :---- | :---- | :---- |
| LiPF6 개시 DOL | 양이온 개환 | 85.2% | 3.8 | – | – |
| DOL \+ YSZ 필러 | 양이온 개환 | 98.5% | – | – | \>4.9 |
| EGDA/AIBN | 라디칼 | ≈70% | – | – | – |
| EGDA/AIBN in EC/DMC | 라디칼 | ≈45\~55% | – | – | – |
| EGDA/DME | 라디칼 | – | 3.6 | 0.55 | 4.7 |
| VEC/PETEA in PVDF-HFP | 라디칼 | – | 6.9 | 0.51 | 5.3 |
| VC/EGDA \+ HFBMA | 라디칼 | – | 0.63 | 0.82 | 5.1 |

리뷰가 정리한 패턴은 분명하다. **전도도가 높은 계는 모두 액체나 가소제를 남긴다.** 양이온 개환은 전환율(약 85~~98%)이 라디칼(약 45~~70%)보다 높다. 대신 DOL계는 고전압에 약해서 NCM622 셀이 100사이클에 그쳤다. 불소화 UV 망(Tang 2023)은 광개시제 0.6 mol%로 1.37 mS/cm를 냈고, 4.5 V NCM523에서 200사이클 후 약 90%를 유지했다([Nature Commun.](https://www.nature.com/articles/s41467-023-37997-6)). 다만 100\~200 µm 두께 막에서 얻은 결과다. 경화 전후 계면저항(Ω·cm²), 경화 수축률(%), 잔류 모노머 정량값은 이번 조사에서 확보하지 못했다. 그래서 M7의 결과는 \*\*설계 창(조건 영역)\*\*으로 제시하고 정량 예측은 하지 않는다.

| 항목 | 내용 |
| :---- | :---- |
| 입력 | M5 σ0,min(φ\_inf) 표, M6 토네이도, 위 표 |
| 절차 | ① 축 φ\_inf(0.4~~1.0)와 X(0.5~~1.0) 위에 등고선 3종: R\_cathode (M5 식), f\_liq,cell (위 식, w\_mon·m\_cell은 \[가정\] 범위), "σ0 ≥ σ0,min" 경계 ② GB/T 43568 경계선(액체 5%, 20%)을 겹침 ③ 화학 후보 2\~3종(양이온 DOL계 / 불소화 아크릴계 / 이중층: 양극 쪽 산화안정 \+ 음극 쪽 에테르) 배치 ④ 공정 조건표(개시 방식·온도·함침 압력·시간) 제시 |
| 출력 | **카솔라이트 설계 창 그림**, 후보 화학 비교표, 우선순위 3개 |
| 검증 기준 | 설계 창 안의 조건을 M6에 넣었을 때 1C 용량이 기준 대비 개선되고, 막 σ ×10 레버보다 큰지 비교 |
| 리스크 | 수축·잔류 모노머 수치 없음 → 범위로 처리하고 "실측 필요 변수"로 명시(고찰 슬라이드) |

### M8 (선택). 실측 없는 새 조성의 σ 예측 (11-07 \~ 11-09, 여유 시)

| 항목 | 내용 |
| :---- | :---- |
| 입력 | D5 (실험 σ), D6·D7 (MD σ, t+) |
| 절차 | ① D5에서 고분자 SMILES·Mw·염·몰랄농도·Tg·T → log σ 회귀(랜덤포레스트 또는 VTF형 특징). 무작위 분할 대신 DOI 단위로 묶어 교차검증해 정보 누설 차단 ② D7 MD σ에 앵커 보정 계수 k \= σ\_exp(PEO, 80 °C, EO:Li ≈ 15\) / σ\_MD(PEO, 353 K) 곱하기(HTP-MD 1.5 mol/kg \= PEO 기준 EO:Li ≈ 15 \[계산\], 80 °C 실험 앵커 ≈ 1.05×10⁻³ S/cm @EO:Li 16\) ③ 후보 고분자의 σ·t+ 순위 → 상위 3개를 M5 지도에 넣어 "새 조성이면 병목이 어떻게 바뀌나" 시연 |
| 출력 | 예측 σ ± 예측구간, 후보별 병목 지도 |
| 검증 기준 | DOI 묶음 교차검증 log σ RMSE 보고, PEO 앵커 재현 |
| 리스크 | MD는 80 °C·1.5 m 한 조건뿐이고 절대값 편향이 있음 → 상대 순위로만 사용. polymernet은 라이선스 미표기 → 인용만. 문헌 라벨 잡음 |

---

## 7\. 시험(분석) 조건 설계 — 실측은 없지만 "조건 설계"는 있다

평가표 1.4는 "시험 조건 설계·이해"를 본다. 이 프로젝트의 "시험"은 공개 데이터 분석과 수치 실험이다. 따라서 아래 표를 "시험 설명 → 시험 장비·SW → 시험 과정"의 근거로 쓴다.

| 시험 | 목적 | 독립 변수 (수준) | 고정 조건 | 측정·출력 | 반복·통계 | 판정 |
| :---- | :---- | :---- | :---- | :---- | :---- | :---- |
| T1 KK 검증 | 데이터 품질 | 스펙트럼별 | — | 잔차(%) vs f | 전 스펙트럼 | 잔차 \< 1\~2% |
| T2 ECM 피팅 | R·C 추출 | 회로 후보 2\~3종 | 가중 1/\|Z\| | R, Q, n, C\_eff ± CI, χ²/ν | 초기값 5종 다중 시작 | CI \< ±30%, AIC |
| T3 합성 회수 | 피팅법 검증 | 잡음 0.5/1/2% | 참값 고정 | 회수 오차 | 100회 | 참값이 95% CI 안 ≥ 90% |
| T4 DRT λ 스윕 | 요소 수 | λ \= 10⁻⁵\~10⁻¹ (5수준) | — | 피크 수·위치 | 온도별 | 10배 범위 안정 |
| T5 온도 의존 | 메커니즘 | T (데이터 가용 전부, 목표 ≥ 6점) | 조성별 | Ea, (A, B, T0) | 구간(Tm 상·하) 분리 | 4장 규칙 |
| T6 조성 효과 | 염 농도 | EO:Li 20/15/10 (D1) | 같은 T | σ, Ea | — | 최적점 위치 |
| T7 병목 지도 | H1·H2 | T 25/60/80 °C · L 10/20/30/50 µm · φ\_inf 0.4/0.6/0.9/1.0 · t+ 0.06/0.15/0.61 | 양극 60 µm, ε 0.3, 3 mAh/cm² | 항목별 ASR, 1순위 | 몬테카를로 5,000 | 2.2절 |
| T8 셀 방전 | 레버 효과 | 레버 8종 × 0.5C/1C | OCV 표, Vcut | 용량, 평균전압, 발열 | 결정론적 | 기준 대비 Δ |
| T9 민감도 | 레버 순위 | 각 파라미터 ±50% | 나머지 기준값 | Δ1C 용량 | — | 토네이도 최장 막대 |
| T10 설계 창 | 해결안 | φ\_inf × X 격자 (13 × 11\) | 화학 후보별 | R\_cathode, f\_liq | — | GB/T 경계 대비 |

**시험 장비·SW.** MATLAB(Optimization Toolbox: `lsqnonlin`, `lsqcurvefit`), Simulink(+선택 Simscape Battery R2025a 이상), Python 3(numpy, scipy, pandas, impedance.py, pyDRTtools), KIT Lin-KK(무료), WebPlotDigitizer(디지타이징). 원 데이터의 측정 장비(전위가변기, 셀 형상)는 Khudyshkina 논문 실험절에서 인용해 "원 시험 조건"으로 소개한다.

---

## 8\. 일정 (2026-10-29 \~ 11-12)과 역할 분담

### 8.1 일정

| 날짜 | 요일 | 작업 | 산출물 | 비고 |
| :---- | :---- | :---- | :---- | :---- |
| 10-29 | 목 | 착수 회의, 역할 확정, Zenodo·CheMixHub 다운로드, README 확인 | 데이터 인벤토리 | — |
| 10-30 | 금 | M1 이론 정리(3장 슬라이드 초안), 앵커표, 문헌 7편 입력 | `15_문헌수집양식` | ADsP 응시자는 이론 검토만 |
| 10-31 | 토 | M2 피팅 코드 합성 데이터로 검증 (T3) | 회수 시험 결과 | **ADsP 시험일 — 해당 팀원 작업 없음** |
| 11-01 | 일 | M2 실데이터 피팅 시작, M3 DRT 환경 구축 | 첫 R·C 표 | — |
| 11-02 | 월 | M2·M3 병행: 회로 확정 | DRT 그림, 채택 회로 | — |
| 11-03 | 화 | **외부 실습(에너웨이브)** — 실습에서 EIS 측정 조건·장비 질문 목록 준비 | 질문지·메모 | 프로젝트 작업 최소화 |
| 11-04 | 수 | M2 전 온도·조성 피팅 완료, M4 시작 | R·C 전체 표 (CI 포함) | 중간 점검 1 |
| 11-05 | 목 | M4 Arrhenius/VTF, 메커니즘 귀속, M5 시작 | Ea 표, 귀속표 | — |
| 11-06 | 금 | M5 병목 지도·몬테카를로·σ\_min 역산, M6 ODE 구현 | 병목 지도 그림 | 중간 점검 2 (H1·H2 판정) |
| 11-07 | 토 | M6 Simulink 구현·검증, M7 설계 창, (M8) | 방전 곡선, 토네이도 | — |
| 11-08 | 일 | M6 시나리오 일괄 실행, M7 완료 | 레버 비교표, 설계 창 | — |
| 11-09 | 월 | 결과 동결(freeze), 슬라이드 본작업 | 슬라이드 v1 | 이후 분석 변경 금지 |
| 11-10 | 화 | 슬라이드 v2, 내부 리허설, 평가표 자가 점검 | 슬라이드 v2 | — |
| **11-11** | **수** | **15:00 1차 제출 (PPT/PDF)** | 1차본 | 오전 최종 교정 |
| **11-12** | **목** | 피드백 반영, **15:00 최종 제출** | 최종본 | — |

### 8.2 역할 분담 템플릿 (이름은 팀이 채움)

| 역할 | 담당 | 주 모듈 | 산출물 책임 | 보조 |
| :---- | :---- | :---- | :---- | :---- |
| PM·이론 총괄 | (이름) | M1, M4 | 이론 슬라이드, 메커니즘 귀속표, 평가표 대응표 | M7 |
| EIS 피팅 리드 | (이름) | M2, M3 | R·C 표(CI), 잔차·DRT 그림, MATLAB 코드 | M4 |
| 셀 모델 리드 | (이름) | M5, M6 | 병목 지도, 몬테카를로, Simulink 모델, 토네이도 | M2 |
| 해결안·데이터 리드 | (이름) | M7, (M8), 데이터 관리 | 설계 창, 후보 비교, 데이터 인벤토리·출처 | M1 |

모든 팀원이 공통으로 맡는 일도 있다. 각자 자기 모듈 슬라이드 2\~3장을 쓰고, 다른 사람 슬라이드의 꼬리표를 교차 점검한다. 성실도(평가표 3번, 8점)의 증거로 회의록과 Git 또는 공유 드라이브 커밋 로그를 남긴다.

---

## 9\. 평가표 대응 — 배점별로 무엇을 보여 줄 것인가

| 평가 항목 | 배점 | 이 프로젝트의 대응 | 증거 슬라이드 |
| :---- | :---- | :---- | :---- |
| 1.1 과업 창의성 | 10 | 통념("막 σ가 병목") 정면 검증. "σ가 아니라 t+"라는 이분법 결과. 병목 교차온도 T\*·t+\* 같은 새 지표 | 표지 메시지, 결론 |
| 1.2 과업 목적성 | 15 | 두괄식 결론, 산업 배경(Blue Solutions·Factorial·WeLion·GB/T 43568), 정량 목표·가설·판정 규칙 | 선정배경, 목표·가설 |
| 1.3 과목 연관성 | 5 | Arrhenius/VTF, Bruce–Vincent, Nernst–Einstein, Bruggeman/TLM, CPE·DRT (훈련 내용인 ECM·진단과 직결), 논문 14편 이상 인용 | 이론 3장 |
| 1.4 과업 진행 과정 | 15 | 가이드 순서(주제→선정배경→이론→시험설명→장비·SW→시험과정→결과→결론→고찰→팀원소개) 준수, 7장 시험 조건 설계표 | 전체 구성 |
| 1.5 시나리오 수립 당위성 | 10 | 시나리오 B\~F의 σ는 \[출처\] 앵커, 가정값은 범위·몬테카를로로 처리, 레버 8종은 4장 메커니즘에서 도출 | 시나리오 설계표 |
| 1.6 결과보고서 완성도 | 7 | 가이드 목차, 그림 중심, 꼬리표 일관성, 출처 목록 | 전체 |
| 1.7 과업 결과 도출 | 20 | 병목 지도·1순위 확률·σ\_min/t+\_min 역산·방전/발열 비교·토네이도·카솔라이트 설계 창 | 결과 4\~6장 |
| 2 MATLAB 등가회로 | 10 | **`lsqnonlin` 피팅 R·C(CPE)·C\_eff 표(95% CI)** \+ Simulink 2RC로 목적(레버 우선순위) 도출 | M2 표, M6 그림 |
| 3 성실도 | 8 | 일정 준수, 회의록, 코드 저장소, 중간 점검 2회 | 수행 경과 |
| 4 가점 | 10 | (선택) M8 ML 예측, 합성 데이터 회수 시험, 몬테카를로 불확실성, 오픈데이터·코드 공개 | 차별화 슬라이드 |

---

## 10\. 보고서 슬라이드 구성안 — 그림 한 장이 문장 하나를 말한다

제출은 PPT/PDF만 받으며, 글 벽(text-wall)은 피한다. 모든 슬라이드는 상단 한 줄 결론과 핵심 그림 1개로 구성한다.

| \# | 가이드 목차 | 슬라이드 제목(메시지) | 핵심 그림 |
| :---- | :---- | :---- | :---- |
| 1 | 표지 | 고분자 전해질의 병목은 막의 σ가 아니다 | 누적 ASR 막대 한 장(시나리오 C) |
| 2 | 개요: 주제·선정배경 | 고분자가 상용 경쟁에 돌아왔다 | Blue Solutions·Factorial·WeLion 타임라인 \+ GB/T 43568 분류 막대 |
| 3 | 개요: 기획의도 | 통념 vs 검증 질문 3개(H1\~H3) | 가설 카드 3장 |
| 4 | 개요: 특화포인트·차별화 | EIS 분해 → 메커니즘 → 셀 환산을 한 흐름으로 | 모듈 흐름도 M1→M7 |
| 5 | 개요: 훈련내용 연관성 | 수업의 ECM·진단을 고분자 셀로 확장 | 수업 항목 ↔ 모듈 매핑표 |
| 6 | 이론 ① | 셀 손실 \= 직렬 저항의 합(ASR) | 셀 단면 \+ 저항 사슬 그림 |
| 7 | 이론 ② | 사슬에 묶인 이온은 VTF로 휜다 | Arrhenius vs VTF 그래프 \+ Ea,app |
| 8 | 이론 ③ | t+가 낮으면 막이 DC에서 10배 커진다 | R\_b vs R\_b/t+ 막대, Bruce–Vincent 개념도 |
| 9 | 이론 ④ | 양극 기공 속 이온은 굴곡 경로를 간다 | Bruggeman·TLM 사다리 그림, λ/d |
| 10 | 시험 설명 | 반원의 정체는 정전용량이 알려 준다 | Nyquist \+ C\_eff 귀속 표 |
| 11 | 장비·재료·개발환경 | 공개 데이터·MATLAB·Simulink·Python | 데이터 계획 요약표 |
| 12 | 시험 과정 ① | KK 검증 → DRT → CNLS 피팅 | 파이프라인 \+ 잔차 그림 |
| 13 | 시험 과정 ② | 시험 조건 설계표 | 7장 표 요약 |
| 14 | 결과 ① | MATLAB 피팅 R·C 표 (95% CI) | 표 \+ Nyquist 피팅 |
| 15 | 결과 ② | DRT 피크 수 \= 회로 요소 수 | 온도별 DRT |
| 16 | 결과 ③ | 저항마다 Ea가 다르다 → 메커니즘 | 저항별 Arrhenius 그래프 |
| 17 | 결과 ④ | 병목 지도: 온도·두께·침투율별 1순위 | 누적 ASR 막대 매트릭스 \+ 1순위 확률 |
| 18 | 결과 ⑤ | σ 기준 vs t+ 기준: H1의 이분법 | HF/DC 두 지도 나란히 \+ t+\* 곡선 |
| 19 | 결과 ⑥ | 레버 비교: 0.5C/1C 방전·발열 | 방전 곡선 겹침 \+ 토네이도 |
| 20 | 결과 ⑦ | 카솔라이트 설계 창 | φ\_inf × X 등고선 \+ GB/T 경계 |
| 21 | 결론 | 개선 우선순위 Top 3 | 우선순위표 (레버·효과·근거) |
| 22 | 고찰 | 한계와 실측이 필요한 변수 | 꼬리표 통계 \+ 추가 실험 제안 |
| 23 | 활용방안·기대효과 | 소재 R\&D 의사결정 도구 | 활용 시나리오 3칸 |
| 24 | 팀 구성·역할 / 수행 경과 / 자체 평가 | 일정 준수와 자가 평가 | 간트 차트 \+ 평가표 자가점수 |
| 25 | 출처 | — | 참고문헌 |

---

## 11\. 리스크와 대응

| 리스크 | 가능성 | 영향 | 대응 |
| :---- | :---- | :---- | :---- |
| Zenodo EIS에 온도 시리즈·셀 형상 정보가 없음 | 중 | M4 약화 | 논문 실험절 확인 → 없으면 D9 디지타이징 \+ D10 합성으로 온도 시리즈 보강하고 "방법 시연"으로 위치 설정 |
| 원시 데이터가 차단셀뿐이라 계면 분해 불가 | 중 | R\_int 출처 부족 | D2 Bruce–Vincent 분극 전후 대칭셀 EIS 사용, PDOL·ORNL 문헌값과 범위 처리 |
| 반원 중첩으로 파라미터 상관 큼 | 높음 | CI 커짐 | DRT 초기값, 요소 축소, CI를 그대로 보고(정직성 자체가 평가 요소) |
| 계면·양극 파라미터 대부분이 \[가정\] | 높음 | 결론 신뢰도 | 몬테카를로 1순위 확률로 보고, "조건부 결론" 명시, 경계값(t+\*, σ\_min) 중심으로 결론 서술 |
| R\_diff 식(이상 용액)과 양극 농도분극 누락 | 중 | DC 지도 왜곡 | "막 병목 상한"으로 해석, 고찰에 명시 |
| Simscape Battery 미설치 | 중 | M6 지연 | 기본 Simulink 블록 구현을 기본안으로 |
| 일정 충돌(ADsP 10-31, 실습 11-03) | 확정 | 작업 손실 2인·일 | 일정표에 반영 완료, 주말 작업 분산 |
| 문헌 수치 단위·라벨 오류 | 중 | 앵커 오차 | 원 그림 대조, 중앙값·범위 사용 |
| 범위 과다(M8) | 중 | 본 모듈 품질 저하 | M8은 11-07까지 M1\~M6 완료 시에만 착수 |

---

## 12\. 기대효과와 가치

**산업.** 병목 지도는 소재 개발팀에 "막 σ를 올리는 데 투자할지, 카솔라이트 침투율을 올리는 데 투자할지"를 묻는 질문에 정량 답을 준다. σ\_min·t+\_min 역산표는 소재 사양서의 목표값으로 바로 쓸 수 있다. 카솔라이트 설계 창은 GB/T 43568 분류(액체 5%·20% 경계)와 성능을 한 그림에 겹쳐 놓으므로, "전고체" 표기를 노리는 인시츄 셀 개발에 유용하다. 같은 방법은 KRICT가 황화물 전고체에 적용한 전구체 함침–경화(E-ICP) 같은 다른 계로도 옮길 수 있다([Newswise 2026](https://www.newswise.com/articles/view/850751/?sc=rssn)).

**학술.** 공개 원시 EIS를 KK 검증 → DRT → CNLS → 온도 귀속 → TLM 셀 환산까지 재현 가능한 코드로 엮은 사례는 드물다. "σ 기준 vs t+ 기준"이라는 이분법은 고분자 전해질 평가 지표를 σ·t+로 옮겨야 한다는 주장을 숫자로 뒷받침한다.

**교육.** 화학 전공 팀이 전기화학 이론(수송, 계면, 다공 전극)과 데이터 분석(비선형 회귀, 불확실성, 정규화 역문제)과 시스템 모델(Simulink)을 한 과제에서 연결해 본다.

**커리어.** 결과물은 배터리 진단·셀 설계·소재 R\&D 직무 포트폴리오에 바로 넣을 수 있다. EIS 피팅 코드, DRT, 셀 모델, 설계 의사결정이 한 프로젝트에 모두 들어 있고, 공개 데이터와 코드라 면접에서 시연할 수도 있다.

---

## 13\. 용어집

| 용어 | 정의 |
| :---- | :---- |
| SPE (고분자 고체 전해질) | 리튬염을 고분자(예: PEO)에 녹인 고체 전해질 |
| GPE (겔 고분자 전해질) | 고분자 망에 액체를 머금은 전해질 |
| PEO | 폴리에틸렌옥사이드, –CH₂CH₂O– 반복 단위(EO) |
| EO:Li | 에테르 산소 수 대 Li 수의 비. r \= \[Li\]/\[EO\]의 역수 |
| LiTFSI | 리튬 비스(트리플루오로메탄설포닐)이미드, 해리가 잘 되는 리튬염 |
| 인시츄 중합/고체화 | 셀 안에 액체 전구체를 넣고 그 자리에서 굳히는 공정 |
| 카솔라이트 | 복합 양극 내부 기공을 채우는 양극용 전해질 |
| 침투율 φ\_inf | 양극 기공 중 전해질이 실제로 채운 부피 분율 |
| 전환율 X | 모노머 중 고분자로 바뀐 분율 |
| EIS | 전기화학 임피던스 분광법, 주파수별 교류 임피던스 측정 |
| Nyquist 그림 | 가로 Z′(실수부), 세로 −Z″(허수부)로 그린 임피던스 그림 |
| ECM | 등가회로모델, R·C·CPE·Warburg 조합으로 임피던스를 표현 |
| CPE | 정위상 요소, 이상적이지 않은 축전기. Z \= 1/(Q(jω)^n) |
| C\_eff | CPE를 실제 정전용량으로 환산한 값 \[F\] |
| Warburg | 확산 임피던스 요소(45° 직선 또는 유한 길이형) |
| DRT | 이완시간 분포, 시간상수별 저항 분포로 임피던스를 분해 |
| Tikhonov 정규화 / λ | 해를 매끄럽게 만드는 벌점 항과 그 강도 |
| GCV | 일반화 교차검증, λ 자동 선택법 |
| Lin-KK | Kramers–Kronig 관계로 데이터의 선형성·인과성·시불변성을 검사 |
| CNLS | 복소 비선형 최소제곱 피팅 |
| ASR | 면적비저항 R·A \[Ω·cm²\] |
| Arrhenius | σ \= σ0 exp(−Ea/RT), 열활성 도약 |
| VTF | σ \= A T^(−½) exp(−B/(T−T0)), 사슬 운동 결합 이동 |
| WLF | VTF와 등가인 고분자 이동 인자 식 |
| Ea | 활성화 에너지 \[kJ/mol 또는 eV\] |
| Tg / T0 / Tm | 유리전이온도 / 이상 유리전이온도 / 용융온도 |
| 분절 운동 | 고분자 사슬의 국소 부분 운동 |
| t+ (수송수) | 전류 중 Li⁺가 나르는 비율 |
| Bruce–Vincent 법 | 대칭셀 정전위 분극 \+ EIS로 t+를 구하는 법 |
| 농도분극 / R\_diff | 염 농도 구배에 의한 추가 과전압 / 그 등가 저항 |
| 한계전류 J\_lim | 전극 근처 염 농도가 0이 되는 전류밀도 |
| Nernst–Einstein | 확산계수로 전도도를 예측하는 식 |
| Ionicity / Haven ratio | 실측 σ와 NE 예측 σ의 비 / 그 역수 |
| PFG-NMR | 펄스장 경사 NMR, 자기확산계수 측정 |
| MSD | 평균제곱변위, MD에서 확산계수를 구하는 데 사용 |
| 굴곡도 τ / MacMullin 수 N\_M | 이온 경로의 우회 정도 / σ0/σ\_eff \= τ/ε |
| Bruggeman 관계 | σ\_eff \= σ0 ε^α |
| TLM (전송선 모델) | 다공 전극을 저항 사다리로 본 임피던스 모델 |
| 반응 침투 깊이 λ | 양극 두께 중 실제로 반응하는 깊이 |
| SEI / CEI | 음극 / 양극 표면 피막 |
| ESW | 전기화학 안정 창(분해 없이 쓸 수 있는 전위 범위) |
| 2RC 모델 | R0 \+ RC 2개 \+ OCV로 구성한 셀 등가회로 |
| OCV | 개방회로전압 |
| 토네이도 차트 | 변수별 ±변화의 결과 영향을 막대로 정렬한 민감도 그림 |
| 몬테카를로 | 불확실한 입력을 무작위로 뽑아 결과 분포를 얻는 방법 |
| AIC | 모델 선택 지표(적합도 \+ 파라미터 수 벌점) |
| GB/T 43568-2026 | 중국 EV 고체전지 용어·분류 권고 국가표준 |

---

## 14\. 출처

**이론·메커니즘**

- Chen et al., Front. Chem. 8, 592604 (2021) — (PEO)16LiTf DSC·σ·열이력: [https\://www\.frontiersin.org/articles/10.3389/fchem.2020.592604/pdf](https://www.frontiersin.org/articles/10.3389/fchem.2020.592604/pdf)  
- Pesko et al., J. Electrochem. Soc. 164, E3569 (2017) — PEO-LiTFSI t+, σ 최적점, Bruce–Vincent 식: [https\://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf](https://escholarship.org/content/qt2630j6bs/qt2630j6bs.pdf)  
- "Transference number in polymer electrolytes: mind the reference-frame gap", arXiv 2203.01659: [https\://arxiv.org/pdf/2203.01659](https://arxiv.org/pdf/2203.01659)  
- "Bruce–Vincent transference numbers from molecular dynamics simulations", arXiv 2302.07017: [https\://arxiv.org/pdf/2302.07017](https://arxiv.org/pdf/2302.07017)  
- Abe, Iriyama, Ogumi, Trans. MRS-J 29, 9 (2004) — 고분자/세라믹 계면 Ea: [https\://www\.mrs-j.org/pub/tmrsj/vol29\_no1/vol29\_no1\_009.pdf](https://www.mrs-j.org/pub/tmrsj/vol29_no1/vol29_no1_009.pdf)  
- Chen et al., ACS Energy Lett. 4, 1080 (2019) — 고분자/세라믹 계면저항 1.2 kΩ·cm²: [https\://impact.ornl.gov/en/publications/determining-and-minimizing-resistance-for-ion-transport-at-the-po/](https://impact.ornl.gov/en/publications/determining-and-minimizing-resistance-for-ion-transport-at-the-po/)  
- Xin et al., Molecules 29, 2454 (2024) — 인시츄 PDOL 겔: [https\://pmc.ncbi.nlm.nih.gov/articles/PMC11173723](https://pmc.ncbi.nlm.nih.gov/articles/PMC11173723)  
- Fukutsuka, Miyazaki et al., J. Electrochem. Soc. (2017) — 복합 전극 이온저항 Ea: [https\://doi.org/10.1149%2F2.1381714jes](https://doi.org/10.1149%2F2.1381714jes)  
- Landesfeind et al., J. Electrochem. Soc. 163, A1373 (2016) — EIS 굴곡도: [https\://doi.org/10.1149/2.1141607jes](https://doi.org/10.1149/2.1141607jes)  
- BioLogic AN-BAT-011 — TLM 임피던스식·굴곡도 예: [https\://icpms.cz/labrulez-bucket-strapi-h3hsga3/AN\_BAT\_011\_3da1cad8af/AN-BAT-011.pdf](https://icpms.cz/labrulez-bucket-strapi-h3hsga3/AN_BAT_011_3da1cad8af/AN-BAT-011.pdf)  
- ETH, J. Electrochem. Soc. 165 (2018) — EIS 굴곡도 vs 3D 단층촬영: [https\://www\.research-collection.ethz.ch/items/5b68f741-3987-4cbb-8485-dbeb1ea75d63](https://www.research-collection.ethz.ch/items/5b68f741-3987-4cbb-8485-dbeb1ea75d63)  
- \[미검증\] Ratner & Shriver, Chem. Rev. 88, 109 (1988): [https\://doi.org/10.1021/cr00083a006](https://doi.org/10.1021/cr00083a006) · Williams, Landel, Ferry, JACS 77, 3701 (1955): [https\://doi.org/10.1021/ja01619a008](https://doi.org/10.1021/ja01619a008) · Chazalviel, Phys. Rev. A 42, 7355 (1990): [https\://doi.org/10.1103/PhysRevA.42.7355](https://doi.org/10.1103/PhysRevA.42.7355) · Irvine, Sinclair, West, Adv. Mater. 2, 132 (1990): [https\://doi.org/10.1002/adma.19900020304](https://doi.org/10.1002/adma.19900020304)

**EIS·DRT·셀 모델 방법**

- Zahner IM7 — CPE 유효정전용량: [https\://doc.zahner.de/im7/en/applications/constant\_phase\_element/index.html](https://doc.zahner.de/im7/en/applications/constant_phase_element/index.html)  
- Hirschorn PhD thesis (UF 2010): [https\://www\.che.ufl.edu/PDF/orazem/Hirschorn-PdD-UF-2010.pdf](https://www.che.ufl.edu/PDF/orazem/Hirschorn-PdD-UF-2010.pdf)  
- KIT Lin-KK: [https\://www\.iam.kit.edu/et/english/Lin-KK.php](https://www.iam.kit.edu/et/english/Lin-KK.php)  
- impedance.py validation: [https\://impedancepy.readthedocs.io/en/latest/validation.html](https://impedancepy.readthedocs.io/en/latest/validation.html)  
- Ciucci lab (pyDRTtools, DRTtools): [https\://github.com/ciuccislab](https://github.com/ciuccislab)  
- DRT 정규화 파라미터 선택 (ChemRxiv 2022): [https\://chemrxiv.org/engage/chemrxiv/article-details/636fbfc5fbfd38d985dbf1bb](https://chemrxiv.org/engage/chemrxiv/article-details/636fbfc5fbfd38d985dbf1bb)  
- MathWorks Simscape Battery `ecm`: [https\://www\.mathworks.com/help/simscape-battery/ref/simscape.battery.parameters.ecm.html](https://www.mathworks.com/help/simscape-battery/ref/simscape.battery.parameters.ecm.html)  
- MathWorks Equivalent Circuit Battery 블록: [https\://www\.mathworks.com/help/autoblks/ref/equivalentcircuitbattery.html](https://www.mathworks.com/help/autoblks/ref/equivalentcircuitbattery.html)

**해결안·레버**

- In-situ 중합 전해질 리뷰, Mater. Chem. Front. 2025 (DOI 10.1039/D5QM00116A): [https\://pubs.rsc.org/en/content/articlehtml/2025/qm/d5qm00116a](https://pubs.rsc.org/en/content/articlehtml/2025/qm/d5qm00116a)  
- Tang et al., Nat. Commun. 14, 2301 (2023): [https\://www\.nature.com/articles/s41467-023-37997-6](https://www.nature.com/articles/s41467-023-37997-6)  
- Wan et al., Nat. Nanotechnol. (2019): [https\://www\.nature.com/articles/s41565-019-0465-3](https://www.nature.com/articles/s41565-019-0465-3)  
- KRICT E-ICP (Energy Storage Materials 2026, 보도): [https\://www\.newswise.com/articles/view/850751/?sc=rssn](https://www.newswise.com/articles/view/850751/?sc=rssn)  
- \[미검증\] Bouchet et al., Nat. Mater. 12, 452 (2013): [https\://doi.org/10.1038/nmat3602](https://doi.org/10.1038/nmat3602)

**데이터**

- Khudyshkina et al. 데이터셋, Zenodo 7428858: [https\://zenodo.org/records/7428858](https://zenodo.org/records/7428858)  
- CheMixHub (Bradford et al. 2023 DB 포함): [https\://github.com/chemcognition-lab/chemixhub](https://github.com/chemcognition-lab/chemixhub) · [https\://arxiv.org/pdf/2506.12231](https://arxiv.org/pdf/2506.12231)  
- Bradford et al., ACS Cent. Sci. (2023): [https\://pmc.ncbi.nlm.nih.gov/articles/PMC9951296](https://pmc.ncbi.nlm.nih.gov/articles/PMC9951296)  
- polymernet (Xie et al. 2022): [https\://github.com/txie-93/polymernet](https://github.com/txie-93/polymernet)  
- HTP-MD: [https\://arxiv.org/pdf/2208.01692](https://arxiv.org/pdf/2208.01692) · [https\://www\.htpmd.matr.io](https://www.htpmd.matr.io)  
- Wang et al., "When Literature Data Mislead AI in Materials Discovery", arXiv 2609.01621: [https\://arxiv.org/pdf/2609.01621](https://arxiv.org/pdf/2609.01621)  
- KIT RDR 후보: [https\://service.tib.eu/ldmservice/dataset/rdr-doi-10-35097-rhsdekozhadetmgk](https://service.tib.eu/ldmservice/dataset/rdr-doi-10-35097-rhsdekozhadetmgk) · ChemRxiv 자동화 플랫폼: [https\://chemrxiv.org/engage/chemrxiv/article-details/675cccadf9980725cf041455](https://chemrxiv.org/engage/chemrxiv/article-details/675cccadf9980725cf041455)

**산업·규제**

- Blue Solutions Gen4 보도자료 (2025-11-17): [https\://www\.blue-solutions.com/wp-content/uploads/2025/11/pr-2025-07-blue-solutions-demonstrates-70-more-autonomy-with-its-breakthrough-gen4-lithium-metal-batteries.pdf](https://www.blue-solutions.com/wp-content/uploads/2025/11/pr-2025-07-blue-solutions-demonstrates-70-more-autonomy-with-its-breakthrough-gen4-lithium-metal-batteries.pdf)  
- Factorial FEST (Electronic Design): [https\://electronicdesign.com/markets/automotive/article/55291872/electronic-design-solid-state-battery-cells-for-fast-charging-high-energy-density-coming-in-2026](https://electronicdesign.com/markets/automotive/article/55291872/electronic-design-solid-state-battery-cells-for-fast-charging-high-energy-density-coming-in-2026)  
- WeLion (MIT Technology Review, 2026-10-06): [https\://technologyreview.com/2026/10/06/1144987/2026-climate-tech-companies-to-watch-welion-semi-solid-state-batteries](https://technologyreview.com/2026/10/06/1144987/2026-climate-tech-companies-to-watch-welion-semi-solid-state-batteries)  
- SMM 고체–액체 전지 생산능력: [https\://news.metal.com/newscontent/104120486-solid-state-battery-capacity-part-1-solid–liquid-batteries-a-transition-route-with-strong-supply-and-demand](https://news.metal.com/newscontent/104120486-solid-state-battery-capacity-part-1-solid–liquid-batteries-a-transition-route-with-strong-supply-and-demand)  
- GB/T 43568-2026: [https\://www\.batterydesign.net/legislation-rules-and-regulations/gb-t-43568-2026/](https://www.batterydesign.net/legislation-rules-and-regulations/gb-t-43568-2026/)  
- 국내 3사 전고체 일정 (아시아경제, 2024): [https\://view.asiae.co.kr/en/article/2024030811265347893](https://view.asiae.co.kr/en/article/2024030811265347893)

---

## 결론

사전 계산을 마치면 프로젝트의 질문 자체가 바뀐다. 처음 질문은 "막이 병목인가?"였다. 출처 있는 σ로 계산한 결과, 막의 이온전도도는 어떤 시나리오에서도 셀 저항의 15%를 넘지 않았다. 그러니 소재 개발 우선순위를 막 σ에 두는 것은 셀 수준에서는 근거가 약하다. 대신 두 질문이 새로 떠올랐다. 하나는 "양극 안 전해질이 기공을 얼마나 채우고 얼마나 잘 통하는가"이고, 다른 하나는 "막이 Li⁺를 얼마나 선택적으로 나르는가(t+)"이다. 양극 안 전해질에는 막보다 5\~20배 높은 σ가 필요하다. t+가 경계값(25 °C 건식 PEO에서 0.19)보다 낮으면 막은 EIS에 보이지 않는 농도분극으로 다시 병목이 된다. 따라서 발표의 결론은 "막을 무시하라"가 아니라 \*\*"막의 σ에서 막의 t+와 양극 카솔라이트로 지표를 옮겨라"\*\*가 된다.

이 결론의 약한 고리도 분명하다. 계면 저항과 양극 구조 파라미터의 상당 부분은 아직 \[가정\]이다. 그래서 팀이 2주 동안 해야 할 일의 핵심은 두 가지다. 첫째, M2\~M4에서 공개 원시 EIS로 계면 저항과 Ea를 출처 있는 값으로 바꾼다. 둘째, 남는 불확실성은 몬테카를로 1순위 확률과 경계값(t+*, σ\_min, φ\_inf*)으로 정직하게 보여 준다. 그렇게 하면 H1이 채택되든 기각되든 "어떤 조건에서 무엇이 병목인가"를 답하는 지도가 남는다. 다음 단계에서 실측이 가장 필요한 변수는 카솔라이트의 경화 수축률과 잔류 모노머이며, 이것을 고찰 슬라이드의 제안으로 삼는다.