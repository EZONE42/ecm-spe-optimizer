# ECM-SPE Optimizer: Solid Polymer Electrolyte Composition Optimization via Equivalent Circuit Modeling

[![CI](https://github.com/EZONE42/ecm-spe-optimizer/actions/workflows/ci.yml/badge.svg)](https://github.com/EZONE42/ecm-spe-optimizer/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/tests-22%20passed-brightgreen)](https://github.com/EZONE42/ecm-spe-optimizer/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

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
ecm-spe-optimizer/
├── ecm_spe/                              # [Core Package] 표준 파이썬 연구 라이브러리
│   ├── __init__.py                       # 패키지 버전(1.0.0) 및 주요 API export
│   ├── ecm.py                            # 4-저항 ECM, de Levie 전송선 모델(TLM), 농도분극 해석
│   ├── transference.py                   # Kasa 원 피팅(circle_fit), Bruce-Vincent t⁺, 계면 ASR
│   ├── features.py                       # RDKit 분자 기술자 추출기 및 오프라인 폴백 딕셔너리
│   ├── ml.py                             # HistGradientBoosting 전도도 대리모델 및 외삽 평가기
│   ├── cli.py                            # 명령줄 인터페이스 (optimize, transference, report, all)
│   └── __main__.py                       # python -m ecm_spe 엔트리포인트
├── notebooks/                            # 재현 가능한 튜토리얼 및 인터랙티브 분석
│   ├── README.md                         # 튜토리얼 실행 가이드
│   └── 01_quickstart_ecm_spe_optimization.py # 4대 저항 분해 및 최적화 재현 스크립트
├── tests/                                # pytest 단위/회귀 테스트 스위트 (22개 테스트)
│   ├── conftest.py                       # 합성 나이퀴스트 반원 및 셀 파라미터 픽스처
│   ├── test_ecm.py                       # TLM 물리 극한, 저항 예산, 농도분극 44% 검증
│   ├── test_transference.py              # Kasa 원 복원, Bruce-Vincent t⁺, 대칭셀 계면 ASR
│   ├── test_features.py                  # RDKit 고분자/염 기술자 및 Tanimoto 유사도
│   ├── test_regression.py                # results/*.csv 실측 결과 일치 회귀 검증
│   └── test_cli.py                       # CLI 서브커맨드 스모크 테스트
├── .github/workflows/
│   └── ci.yml                            # GitHub Actions 자동 CI 매트릭스 (Python 3.10/3.11)
├── pyproject.toml                        # PEP 517/621 빌드 설정 및 CLI 콘솔 스크립트 (ecm-spe)
├── Makefile                              # 개발자 자동화 명령어 (install, test, run, clean)
├── CITATION.cff                          # 표준 학술 인용 메타데이터 (단일 연구자: Jiwon Lee)
├── LICENSE                               # MIT 라이선스 (EZONE42 / Jiwon Lee)
├── pytest.ini                            # pytest 경로 및 설정
├── data/                                 # 원시 데이터셋 메타데이터 및 스키마
│   └── README.md                         # KIT Zenodo, MIT CheMixHub, Pesko 2017 출처 및 단위
├── src/                                  # 레거시 파이프라인 스크립트 (load_db, transfer_test 등)
├── matlab/                               # MATLAB 등가회로 피팅 스크립트 (Python 결과와 1:1 교차 검증)
├── results/                              # 정량 결과 테이블 및 플롯 (재현 기준점)
├── webapp/                               # 인터랙티브 인포그래픽 & 웹 뷰어 (HTML/JS)
└── docs/                                 # 상세 기술 계획서 및 리포트 자료
```
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

### 5.1 패키지 설치 및 CLI 실행 (추천)
```bash
# ecm_spe 패키지 개발 모드 설치
pip install -e .

# 1) 셀 4-저항 ECM 및 EO:Li 비율 최적화 스위프 실행
ecm-spe optimize --temp 55 --crate 0.5

# 2) KIT BPE-LiTFSI 대칭셀 Bruce-Vincent 전이수 및 계면 저항 분석
ecm-spe transference

# 3) 농도분극 병목(44%) 및 ML 외삽 한계 종합 분석 리포트 출력
ecm-spe report

# 4) 전체 테스트 스위트 실행 (22개 테스트 통과)
pytest -v
```

### 5.2 재현 튜토리얼 스크립트 실행
```bash
python notebooks/01_quickstart_ecm_spe_optimization.py
```

### 5.3 레거시 파이프라인 및 MATLAB 스크립트
```bash
# Python 레거시 분석 스크립트
PYTHONPATH=src python src/optimize_demo.py
PYTHONPATH=src python src/transfer_test.py

# MATLAB 등가회로 피팅 및 전이수 계산
# run('matlab/run_eis_anchor.m')
# run('matlab/run_transference.m')
```

### 5.4 인터랙티브 인포그래픽 뷰어
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

---

## 7. 연구자 및 인용 (Author & Citation)

**Author:** 이지원 (Jiwon Lee / [@EZONE42](https://github.com/EZONE42))

본 프로젝트의 코드 또는 분석 결과를 인용하실 경우 아래의 BibTeX 형식을 사용해 주시기 바랍니다:

```bibtex
@software{Lee_ECM_SPE_2026,
  author = {Lee, Jiwon},
  title = {{ECM-SPE Optimizer: Solid Polymer Electrolyte Composition Optimization via Equivalent Circuit Modeling}},
  year = {2026},
  url = {https://github.com/EZONE42/ecm-spe-optimizer},
  version = {1.0.0}
}
```
