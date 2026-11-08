# ECM-SPE: Solid Polymer Electrolyte Composition Optimization via Equivalent Circuit Modeling

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2026b-orange)](https://www.mathworks.com/)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-green)](https://creativecommons.org/licenses/by/4.0/)

> **핵심 한 줄 요약:** 막 이온전도도($\sigma$) 단일 지표가 아닌, **셀 전체의 4대 저항(막 저항 · 양극 다공성 전송선 저항 · 전이수 기반 농도분극 저항 · 리튬 계면 저항)을 통합한 등가회로 모델(ECM) 전압 손실($\Delta V$)을 채점기(Scorer)로 사용**하여 고분자 전해질(PEO-LiTFSI / BPE-LiTFSI)의 최적 염 비율($\text{EO:Li} = 15:1 \sim 21:1$)과 1순위 병목(농도분극 $44\%$)을 규명하고, 분자 구조 기반 머신러닝(ML)의 외삽 한계(유사 염 교체 $\approx 4\times$ 오차 vs 이종 고분자 교체 $\approx 22\times$ 오차)를 물리화학 메커니즘으로 검증한 연구 패키지입니다.

---

## 1. 연구 배경 및 문제의식 (Problem Statement)

고분자 전해질(Solid Polymer Electrolyte, SPE)은 누액과 발화 위험을 줄일 수 있는 차세대 전고체/준고체 배터리 소재이지만, 상온 이온전도도($\sigma$)와 리튬 양이온 전이수($t^+$)가 낮다는 한계가 있습니다.

기존 소재 탐색 연구는 주로 **단독 막의 이온전도도($\sigma$)가 가장 높은 조성**을 찾는 데 집중합니다. 그러나 실제 배터리 셀 내부에서는 다음 세 가지 빈틈이 발생합니다:
1. **소재 물성과 셀 성능 간의 번역 부재:** 막의 $\sigma$가 높아도 양극 복합전극 내부의 구불구불한 이온 경로(Tortuosity)와 계면 저항이 지배적일 수 있습니다.
2. **보이지 않는 저항, 농도분극($R_{\text{conc}}$):** 고주파 교류 임피던스(EIS)에서는 보이지 않지만, 직류 충·방전 시에는 낮은 전이수($t^+ \approx 0.06\sim 0.18$)로 인해 음이온이 한쪽에 쌓이는 농도분극이 막 저항을 사실상 $(1/t^+ - 1)$배 증폭시킵니다.
3. **머신러닝(ML) 예측의 신뢰 경계 미규명:** 공개 데이터베이스로 학습한 ML 모델이 학습하지 않은 새로운 염이나 고분자 사슬을 만났을 때 **어디까지 맞고 왜 틀리는지**에 대한 정량적 경계가 부재합니다.

본 프로젝트는 **조성 $\rightarrow$ 물성($\sigma, t^+$) $\rightarrow$ 4-저항 ECM 파라미터 $\rightarrow$ 셀 전압 손실($\Delta V$)** 로 이어지는 End-to-End 파이프라인을 구축하고, KIT 공개 원자료(Zenodo 7428858)와 MIT CheMixHub DB(16,009행, 217편 논문)를 결합해 이를 검증했습니다.

---

## 2. 아키텍처 및 리포지토리 구조 (Repository Structure)

```text
ecm-spe-electrolyte-optimization/
├── src/                                  # Python 분석·머신러닝·최적화 파이프라인
│   ├── load_db.py                        # [M1] CheMixHub DB SMILES 정규화 및 계열별 집계
│   ├── features.py                       # [M2] RDKit 기반 고분자 반복단위(7종)·염 음이온(5종) 분자 표현자 추출
│   ├── transfer_test.py                  # [M5] 물질 제외 교차 검증 (A:B 학습 -> A:C 염 교체 / D:B 고분자 교체 6사례)
│   ├── cause_check.py                    # [M5] 예측 실패 원인 분석 (PEC/PTMC 농도 거동, Tg 영향, PEO 결정화 구간)
│   ├── optimize_demo.py                  # [M4/M6] PEO:LiTFSI 비율 탐색 + 전송선(TLM) ECM 채점기
│   ├── eis_anchor.py                     # [M3] KIT BPE 전해질 원시 EIS 486 스펙트럼 일괄 피팅 및 VTF/Arrhenius 해석
│   ├── bv_analysis.py                    # [M3] Bruce-Vincent 전이수(t+) 및 대칭셀 계면 저항(Kasa 원 피팅) 산출
│   ├── external_test.py                  # [M5] CheMixHub 학습 ML 모델의 KIT BPE 실험 데이터 외부 검증
│   └── format_eis_drt.py                 # 원시 EIS 텍스트 데이터를 DRT/회로 피팅용 CSV로 변환
├── matlab/                               # MATLAB 등가회로 피팅 스크립트 (Python 결과와 1:1 교차 검증)
│   ├── run_eis_anchor.m                  # 차단전극 EIS R_b 절편/lsqnonlin 피팅 및 가열·냉각 아레니우스 플롯
│   └── run_transference.m                # 대칭셀 EIS + CA 기반 Bruce-Vincent t+ 및 계면 ASR 산출
├── results/                              # 사전 실행된 정량 결과 테이블 및 플롯 (재현 기준점)
│   ├── eis_anchor_arrhenius.png          # 가열/냉각 사이클별 BPE 전해질 아레니우스 플롯
│   ├── eis_anchor_results.csv            # 486개 EIS 스펙트럼 전체 피팅 결과
│   ├── eis_anchor_sigma_by_T.csv         # 시료·온도·가열/냉각 방향별 평균 전도도 및 표준편차
│   ├── eis_anchor_temperature_fits.csv   # 시료별 활성화에너지(Ea) 및 VTF 파라미터(A, B, T0, R²)
│   ├── li_tplus.csv                      # BPE-LiTFSI (15:1, 55°C) Bruce-Vincent 전이수 및 계면 저항
│   ├── bv_eis_intercepts.csv             # 대칭셀 분극 전후 EIS 저항 분리표
│   ├── bv_results_NaK.csv                # Na/K 비교군 전이수 민감도 분석 결과
│   ├── transfer_results.csv              # 6대 교차 검증 케이스 오차 배수 및 순위상관 요약표
│   ├── external_test_results.csv         # KIT BPE 외부 데이터에 대한 ML 예측 vs 실측 비교표
│   └── optimize_demo.csv                 # 25°C / 60°C EO:Li 비율(4~40)에 따른 셀 저항 예산 및 전압 손실
├── webapp/                               # 인터랙티브 조성 최적화 인포그래픽 & 웹 뷰어 (React + Design Component)
│   ├── Main.dc.html                      # 브라우저에서 바로 여는 인터랙티브 대시보드
│   ├── support.js
│   └── vendor/
└── docs/                                 # 상세 기획서, 결과보고서 발표 원고 및 후보 주제 아카이브
    ├── SPE_Optimization_Master_Plan_ko.md       # 프로젝트 마스터 계획서 (수식·이론·검증·직무 연결 완비)
    ├── SPE_Bottleneck_Plan_v1_ko.md             # 고분자 전해질 병목 분석 초기 기획서 (v1)
    ├── SPE_Final_Report_Slides_Text_ko.txt      # 최종 결과보고서 발표 슬라이드 전체 텍스트 (30p 분량)
    └── proposals_archive/                       # [아카이브] 주제 선정 과정에서 설계한 타 도메인 ECM 기획안 3종
        ├── README.md                            # 후보 주제 아카이브 안내 및 비교
        ├── 01_Fault_Aware_Digital_Twin_Proposal.pdf
        ├── 02_Fault_Aware_Digital_Twin_Execution_Plan.pdf
        ├── 03_BBU_Arrhenius_ECM_Plan.tsv
        └── 04_UAM_Thermal_Runaway_ECM_Proposal.tsv
```

---

## 3. 핵심 수식 및 4-저항 ECM 채점 모델 (Methodology)

본 프로젝트는 단순 막 저항만 보는 대신, 셀 단위의 면적비저항($\text{ASR}, \Omega\cdot\text{cm}^2$)을 4가지 물리적 성분으로 분해하여 합산합니다:

$$\Delta V(x) = i \cdot \left[ R_{\text{m}}(x) + R_{\text{cath}}(x) + R_{\text{conc}}(x) + R_{\text{int}}(x) \right]$$

1. **막 벌크 저항 ($R_{\text{m}}$):**
   $$R_{\text{m}} = \frac{L_{\text{sep}}}{\sigma(x, T)}$$
2. **양극 내부 다공성 이온 경로 저항 ($R_{\text{cath}}$, de Levie 전송선 모델 · Bruggeman 보정):**
   $$\sigma_{\text{eff}} = \sigma \cdot \varepsilon^{1.5}, \quad R_{\text{ion}} = \frac{L_{\text{cath}}}{\sigma_{\text{eff}}}, \quad R_{\text{cath}} = \sqrt{R_{\text{ion}} R_{\text{ct}}} \coth\left(\sqrt{\frac{R_{\text{ion}}}{R_{\text{ct}}}}\right)$$
3. **농도분극 저항 ($R_{\text{conc}}$, Bruce–Vincent 전이수 연계):**
   $$R_{\text{conc}} \approx R_{\text{m}} \left( \frac{1}{t^+(x)} - 1 \right), \quad t^+ = \frac{I_{\text{ss}}(\Delta V - I_0 R_{i,0})}{I_0(\Delta V - I_{\text{ss}} R_{i,\text{ss}})}$$
4. **전극-전해질 계면 저항 ($R_{\text{int}}$):**
   $$\text{Li}|\text{SPE}|\text{Li 대칭셀 EIS 반원 지름 } R_i \text{로부터 단면 기준 } R_{\text{int}} = \frac{R_i \cdot A}{2}$$

---

## 4. 주요 검증 결과 요약 (Key Findings)

### 4.1 KIT 블록공중합체 전해질(BPE) 실측 앵커 및 최적 염 비율 (`results/eis_anchor_*.csv`, `li_tplus.csv`)
- **온도별 전도도($\sigma$):** $5\sim 85\,^\circ\text{C}$ 전 구간에서 $\text{EO:Li} = 15:1$ 조성이 $20:1$ 및 $10:1$ 대비 가장 높은 이온전도도를 기록했습니다 ($55\,^\circ\text{C}$ 기준 $15:1 = 6.1\times 10^{-5}\text{ S/cm}$ vs $20:1 = 2.7\times 10^{-5}$, $10:1 = 3.2\times 10^{-5}$).
- **DSC 열분석을 통한 물리적 원인 규명:**
  - $20:1$ (염 부족): 냉각 시 결정화 발생 $\rightarrow$ 저온 전도도 급락
  - $10:1$ (염 과다): 이온-고분자 가교로 유리전이온도($T_g$) 상승($-53\,^\circ\text{C} \rightarrow -41\,^\circ\text{C}$) $\rightarrow$ 사슬 분절운동 둔화
  - $15:1$ (최적 균형): 비정질 상태 유지 + 낮은 $T_g$
- **셀 전압 손실($\Delta V$) 및 1순위 병목:** $55\,^\circ\text{C}, 0.5\text{C}$ 기준 $\text{EO:Li} = 15:1$의 전압 손실은 $0.99\text{ V}$로 타 조성($1.71\sim 1.84\text{ V}$)의 절반 수준이나, **전체 저항의 $44\%$가 농도분극($t^+ = 0.147$)에서 발생**하여 향후 소재 개발의 1순위 타깃이 $\sigma$가 아닌 $t^+$ 개선임을 정량적으로 입증했습니다.

### 4.2 머신러닝(HistGradientBoosting + RDKit) 물질 교차 검증 (`results/transfer_results.csv`)
특정 물질 조합을 학습 데이터에서 완전히 배제한 뒤(Zero-shot Group-out) 예측력을 검증하여 **모델의 사용 가능 범위(Applicability Domain)** 를 획정했습니다:

| 검증 구분 | 제외 및 예측 대상 계열 | 데이터 수 (행 / 논문) | 평균 오차 배수 (Fold Error) | 농도별 순위상관 ($\rho$) | 최적 비율 일치율 | 판정 및 물리화학적 원인 |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **A:C 염 교체** | `PEO : LiFSI` (이미드 유사체) | 45 / 4 | **3.8×** | 0.09 | - | 절대 전도도 오차 ~4배 수준으로 탐색 가능 |
| **A:C 염 교체** | `PEO : LiBETI` (장쇄 불소 이미드) | 83 / 2 | **4.6×** | **0.95** | **7 / 8** | 비율 경향과 최적 농도 위치를 정확히 추종 |
| **A:C 염 교체** | `PEO : LiTfO` (설폰산염) | 763 / 28 | **3.7×** | **0.80** | 34 / 67 | 다른 음이온 계열에서도 농도 경향성 유지 |
| **D:B 고분자 교체** | `P(EO-MO) : LiTFSI` (유사도 0.47) | 72 / 1 | **2.1×** | 0.26 | 0 / 9 | 에테르 산소 골격이 유사해 절대 오차는 작음 |
| **D:B 고분자 교체** | `PEC : LiTFSI` (카보네이트, 0.38) | 65 / 4 | **21.9×** | **-0.20** | 3 / 10 | **실패:** PEO와 달리 고농도($13.9\text{ m}$)까지 $\sigma$ 단조 증가 ($\rho=1.0$) |
| **D:B 고분자 교체** | `PTMC : LiTFSI` (카보네이트, 0.29) | 98 / 3 | **20.7×** | 0.23 | 4 / 12 | **실패:** 사슬 운동-이온 이동 결합 메커니즘 상이 |

> **결론:** 기존 공개 DB 기반 ML 모델은 **동일 고분자 골격 내 염·농도·온도 최적화(오차 $\approx 4\times$)** 에만 신뢰성 있게 사용할 수 있으며, 카보네이트계 등 **이온 전도 메커니즘이 다른 신규 고분자 탐색(오차 $> 20\times$)** 시에는 반드시 소량의 앵커 실험 또는 분자동역학(MD) 보정이 선행되어야 합니다.

---

## 5. 실행 방법 (Quick Start)

### 5.1 Python 환경 설정 및 파이프라인 실행
```bash
pip install pandas numpy scipy scikit-learn rdkit

# 1) KIT 원시 EIS 데이터 분석 (데이터 폴더 경로를 인자로 전달)
python src/eis_anchor.py /path/to/Zenodo7428858/EIS
python src/bv_analysis.py /path/to/Zenodo7428858/T+_EIS_CA

# 2) CheMixHub DB 기반 ML 학습·교차검증·조성 최적화
# (환경변수 CHEMIXHUB_CSV 또는 data/PolymerElectrolyteData.csv 배치 후 실행)
PYTHONPATH=src python src/load_db.py
PYTHONPATH=src python src/transfer_test.py
PYTHONPATH=src python src/cause_check.py
PYTHONPATH=src python src/optimize_demo.py
PYTHONPATH=src python src/external_test.py
```

### 5.2 MATLAB 스크립트 실행
MATLAB에서 `matlab/` 폴더로 이동한 뒤 원시 데이터 파일과 함께 아래 스크립트를 실행합니다:
```matlab
run('matlab/run_eis_anchor.m')       % 486개 스펙트럼 R_b 피팅 및 Arrhenius/VTF 분석
run('matlab/run_transference.m')     % Bruce-Vincent 전이수 및 계면 ASR 계산
```

### 5.3 인터랙티브 인포그래픽 뷰어 실행
```bash
cd webapp
python -m http.server 8080
# 브라우저에서 http://localhost:8080/Main.dc.html 접속
```

---

## 6. 데이터 출처 및 라이선스 (References & Credits)

- **KIT BPE Electrolyte Dataset:** Khudyshkina et al., *Electrochim. Acta* 454, 142421 (2023) — [Zenodo 7428858](https://doi.org/10.5281/zenodo.7428858) (CC BY 4.0)
- **MIT CheMixHub Polymer Electrolyte DB:** Bradford et al., *ACS Cent. Sci.* 9(2), 206 (2023) — [github.com/chemcognition-lab/chemixhub](https://github.com/chemcognition-lab/chemixhub) (CC BY 4.0)
- **PEO-LiTFSI Transport Properties:** Pesko et al., *J. Electrochem. Soc.* 164(11), E3569 (2017)
