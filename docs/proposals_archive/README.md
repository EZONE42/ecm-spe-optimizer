# ECM 프로젝트 주제 탐색 및 후보 기획안 아카이브 (`ECM_주제_퇴고`)

이 디렉토리는 최종 팀 프로젝트 주제(`Universal_ECM_Tool` 및 `BESS 디지털트윈`)가 확정되기 전, **2RC 등가회로 모델(ECM)을 산업 현장의 다양한 문제 해결에 어떻게 적용할 수 있는가**를 치밀하게 설계하고 검증했던 3가지 도메인 기획안을 보존한 아카이브입니다.

각 기획안은 단순 아이디어를 넘어 **배경 규제·산업 표준, 수학적 지배방정식, 검증 프로토콜(데이터 분할 규칙), 베이스라인 비교 지표(KPI)** 까지 완비하고 있어, 향후 별도의 캡스톤 프로젝트, 산학 과제 제안서, 또는 직무 면접(포트폴리오)의 소재로 즉시 재활용할 수 있습니다.

---

## 수록된 기획안 목록

### 1. Fault-Aware 배터리 디지털 트윈 및 AX 진단 게이트 (ESS 운영 · 제조 검사 AX)
- **파일:**
  - [`01_Fault_Aware_Digital_Twin_Proposal.pdf`](01_Fault_Aware_Digital_Twin_Proposal.pdf) (연구 제안서, 3p)
  - [`02_Fault_Aware_Digital_Twin_Execution_Plan.pdf`](02_Fault_Aware_Digital_Twin_Execution_Plan.pdf) (10일 실행계획서 및 MVP 사양서, 4p)
- **핵심 문제의식:**
  - 디지털 트윈(2RC ECM)을 현장 실측 데이터로 지속 보정(Adaptive Update)할 때, **전압 센서 오프셋·드리프트·노이즈를 실제 배터리 열화(내부저항 $R_0$ 증가)로 오인하여 모델 파라미터 자체를 오염시키는 문제(Model Contamination)** 가 발생합니다.
- **핵심 솔루션:**
  - 무부하/펄스/휴지 구간의 전압 잔차(Residual) 시간적 패턴을 AI(Random Forest / XGBoost + Conformal 불확실성)로 분류하여, **[정상 갱신 / 열화 재확인 후 갱신 / 센서 점검 권고 및 갱신 보류]** 를 결정하는 **진단-갱신 연동형 신뢰도 게이트(Reliability Gate)** 를 설계했습니다.
- **평가 설계 (3대 베이스라인 비교):**
  - `B0 (무갱신)` vs `B1 (무조건 갱신)` vs `B2 (단순 잔차 임계값)` vs `제안 (AI 신뢰도 게이트)`
  - 핵심 KPI: 오갱신률(False Update Rate), 파라미터 드리프트 오차($|\hat{R}_0 - R_{0,\text{true}}|/R_{0,\text{true}}$), 실제 열화 추종 지연 시간(Tracking Delay).

---

### 2. AI 데이터센터 BBU(Battery Backup Unit) 저온 출력 병목 진단 및 셀 어레이 설계
- **파일:** [`03_BBU_Arrhenius_ECM_Plan.tsv`](03_BBU_Arrhenius_ECM_Plan.tsv) (파라메트릭 자동 계산 설계 시트)
- **핵심 문제의식:**
  - AI 데이터센터의 GPU 클러스터 전력 진동 및 정전 시 즉각 투입되는 BBU 배터리가 저온 환경에서 출력이 급감할 때, **옴 저항($R_0$), 전하이동 저항($R_1$), 확산 저항($R_2$) 중 어느 성분이 출력을 가로막는 진짜 범인인지** 규명하고 최소 병렬 셀 수를 산출합니다.
- **핵심 솔루션:**
  - 온도별 2RC 파라미터의 아레니우스 기울기(활성화에너지 $E_a$)를 분리 추출하여 저온 출력 병목 소자를 특정하고, 미학습 온도($10\,^\circ\text{C}$) 홀드아웃 검증에서 **테이블 보간 vs 아레니우스 물리식 vs 순수 AI vs 물리+AI 하이브리드(PINN 형태)** 4종 모델을 비교합니다.
  - 허용 C-rate($\le 6\text{C}$) 제약 하에 GPU 부하(셀당 $17.0\text{ A}$)를 견디기 위한 최소 병렬 셀 수($69\text{개}$)를 자동 역산합니다.

---

### 3. UAM(도심항공모빌리티, eVTOL) 고출력 배터리팩 열전파 디지털 트윈 및 CSFL 안전 착륙 시간 기여도 분석
- **파일:** [`04_UAM_Thermal_Runaway_ECM_Proposal.tsv`](04_UAM_Thermal_Runaway_ECM_Proposal.tsv) (연구 제안 및 8주 로드맵)
- **핵심 문제의식:**
  - 전기차(GB 38031-2020, 5분 탈출 대기)와 달리 비행 중인 eVTOL은 멈춰서 내릴 수 없으므로, 유럽 항공안전청(EASA)의 **CSFL(Continued Safe Flight and Landing, 열폭주 발생 시 착륙까지 버티는 시간)** 인증 기준을 만족해야 합니다. 그러나 항공기 특성상 방화벽 무게를 무작정 늘릴 수 없습니다.
- **핵심 솔루션:**
  - 카네기멜런대(CMU) 공개 eVTOL 데이터(18650 셀 22개, 21,392 사이클)의 착륙 직후 휴지 구간에 2RC 피팅을 적용하고, **Bernardi 발열식과 결합한 전기-열 연성(Electro-Thermal Coupled) 2RC 모델** 및 5-셀 열전파 시뮬레이션을 구축합니다.
  - 수천 건의 가상 실험을 학습한 AI 대리 모델(Surrogate Model)과 전역 민감도 분석(SALib)을 통해 **[분리막 내열성 vs 셀 간 차단재 두께 vs 냉각 세기] 중 단위 무게($\text{kg}$)당 안전 착륙 시간을 가장 많이 늘려주는 설계 변수 순위**와 설계 지도(Design Map)를 도출합니다.
