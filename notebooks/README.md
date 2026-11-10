# Interactive Walkthrough Notebooks

This directory contains executable reproduction scripts and interactive walkthroughs for researchers, evaluators, and hiring managers.

## Available Tutorials

1. **`01_quickstart_ecm_spe_optimization.py`**:
   - Single-script demonstration of 4-resistance cell breakdown ($R_{\text{sep}}$, $R_{\text{cath}}$, $R_{\text{conc}}$, $R_{\text{int}}$).
   - Reproduces the **$44\%$ concentration polarization bottleneck**.
   - Sweeps $\text{EO:Li}$ salt ratio grid from $4:1$ to $40:1$ at $25\,^\circ\text{C}$ and $60\,^\circ\text{C}$ to find the minimum overpotential ($\Delta V$) design.

## How to Run

```bash
# Run directly from repository root:
python notebooks/01_quickstart_ecm_spe_optimization.py
```
