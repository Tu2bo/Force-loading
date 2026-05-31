# Force Loading on Molecular Clutches Governs the Stability of Cell Lamellipodia

This repository contains the experimental data and computational procedures associated with the PNAS article:

Ruihao Xue, Lezi Kang, Yonggang Chen, Haoxiang Yang, Hongyuan Jiang, and Ze Gong.
**Force loading on molecular clutches governs the stability of cell lamellipodia.**
*Proceedings of the National Academy of Sciences* (2026).
DOI: https://doi.org/10.1073/pnas.2604349123

The archived dataset and computational procedures are also available on Zenodo:
https://doi.org/10.5281/zenodo.18495893

## Files

### `MainTime.m`

Main simulation script. This script runs stochastic Monte Carlo simulations and mean-field ODE simulations of the BR-integrated motor–clutch model, and analyzes the cell spreading radius, loading rate, and force magnitude.

### `MC_BR_motor_clutch_model.m`

Stochastic Monte Carlo implementation of the BR-integrated motor–clutch model.

### `MF_BR_motor_clutch_model.m`

Mean-field ODE approximation of the BR-integrated motor–clutch model. The state variables include the total adhesion force, the number of bound clutches, and the cell spreading radius.

### `Combined_Analysis.m`

Post-processing and statistical analysis of oscillation data, including oscillation amplitude, period, and power spectral density (PSD) analysis.

### `All_Data_Combined.xlsx`

Experimental data file used by `Combined_Analysis.m`. This file should be placed in the same directory as the analysis script.

## Usage

Run the main simulation script:

```matlab
MainTime
```

Run the oscillation data analysis script:

```matlab
Combined_Analysis
```

## Citation

If you use this code or data, please cite:

Xue R, Kang L, Chen Y, Yang H, Jiang H, Gong Z.
**Force loading on molecular clutches governs the stability of cell lamellipodia.**
*Proceedings of the National Academy of Sciences* (2026).
https://doi.org/10.1073/pnas.2604349123
