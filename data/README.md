# Data Provenance, Schemas, and Units

This directory houses or references the experimental and database files utilized by **ECM-SPE**.
Because raw academic datasets are large (>20 MB), raw CSV/ZIP files are excluded from git tracking via `.gitignore`.
Precomputed benchmark summaries are version-controlled in `results/`, and offline test fixtures are provided in `tests/fixtures/`.

---

## 1. Primary Datasets

### A. MIT CheMixHub Polymer Electrolyte Database
* **Reference**: Bradford et al., *ACS Central Science* 9(2), 206–218 (2023).
* **Repository**: [https://github.com/chemcognition-lab/chemixhub](https://github.com/chemcognition-lab/chemixhub)
* **License**: Creative Commons Attribution 4.0 International (CC BY 4.0)
* **Expected Location**: `data/PolymerElectrolyteData.csv` (or set `CHEMIXHUB_CSV` environment variable)
* **Scope**: 16,009 rows compiled from 217 peer-reviewed journal papers.
* **Key Columns & Units**:
  | Column Name | Physical Meaning | Unit | Notes |
  | :--- | :--- | :---: | :--- |
  | `S1 SMILES` | Polymer repeat unit SMILES | string | Cu/Au capping tags denote repeat boundary |
  | `Salt1 SMILES` | Salt anion/cation SMILES | string | Dot-separated ionic pairs |
  | `Temperature (oC)` | Measurement temperature | °C | Converted to Kelvin ($T_K$) |
  | `log Conductivity (S/cm)`| Ionic conductivity logarithm | $\log_{10}(\text{S/cm})$ | Target $y$ for ML regression |
  | `Salt1 Molality (mol salt/kg polymer)` | Salt concentration | $\text{mol salt / kg polymer}$ | $m$ (molality) |
  | `S1 Mn or Mw` | Polymer molecular weight | $\text{g/mol}$ | Weight/number average MW |
  | `DOI` | Literature source DOI | string | Grouping key for GroupKFold validation |

### B. KIT Block Polymer Electrolyte (BPE) Experimental Dataset
* **Reference**: Khudyshkina et al., *Electrochimica Acta* 454, 142421 (2023).
* **Zenodo Repository**: [https://doi.org/10.5281/zenodo.7428858](https://doi.org/10.5281/zenodo.7428858)
* **License**: Creative Commons Attribution 4.0 International (CC BY 4.0)
* **Scope**: 486 raw EIS spectra across 9 cell samples, temperature sweeps ($5\sim 85\,^\circ\text{C}$), heating/cooling cycles, and symmetric Li\|BPE\|Li potentiostatic chronoamperometry (CA).
* **Key Parameters**:
  * Electrolyte thickness ($L$): $70 \sim 145\,\mu\text{m}$ (sample-specific)
  * Electrode geometric surface area ($A$): $0.785\,\text{cm}^2$ (blocking SS) or $0.5026 / 2.0106\,\text{cm}^2$ (symmetric)
  * Potentiostatic polarization ($\Delta V$): $50\,\text{mV} = 0.05\,\text{V}$

### C. PEO-LiTFSI Transport Properties Reference
* **Reference**: Pesko et al., *Journal of The Electrochemical Society* 164(11), E3569 (2017).
* **DOI**: [10.1149/2.0581711jes](https://doi.org/10.1149/2.0581711jes)
* **Scope**: Measured transference numbers ($t^+ \approx 0.06\sim 0.18$) and thermodynamic factors across salt concentrations in high-molecular-weight PEO.

---

## 2. Converting Units

To convert between molality ($m$, $\text{mol/kg}$) and stoichiometric ratio ($\text{EO:Li}$):
$$\text{EO:Li} = \frac{1}{m \times M_{\text{EO}}}$$
where $M_{\text{EO}} = 44.05 \times 10^{-3}\,\text{kg/mol}$ is the molar mass of the ethylene oxide repeat unit ($-\text{CH}_2\text{CH}_2\text{O}-$).

Example:
* $\text{EO:Li} = 15:1 \implies m = \frac{1}{15 \times 0.04405} \approx 1.513\,\text{mol/kg}$
* $\text{EO:Li} = 20:1 \implies m = \frac{1}{20 \times 0.04405} \approx 1.135\,\text{mol/kg}$
