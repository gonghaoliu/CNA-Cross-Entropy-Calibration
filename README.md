# Cross-Entropy Calibration for Nonprobability Health Surveys

This repository contains the R code for the manuscript **“Implicit Doubly Robust Cross-Entropy Calibration and Inference for Nonprobability Health Surveys with Aggregated Benchmarks.”** The code includes the CNA data analysis, the primary simulation, and supplementary simulations.

## Files

| File | Purpose |
| --- | --- |
| `cna_result.R` | CNA data analysis and population calibration results, details can be found in section 5 and supplementary section 4.|
| `main_simulation.R` | Primary simulation comparing exponential tilting (ET) with alternative calibration specifications, details can be found in section 6 and supplementary section 5. |
| `simulation_categorical_covariates_only.R` | Supplementary simulation comparing poststratification, grouped inverse probability weighting (IPW), and ET with categorical covariates using joint population cell counts, details can be found in supplementary section 6.3. |
| `simulation_kappa1.R` | Supplementary simulation using different selection-strength setting, \(\kappa=1\), details can be found in supplementary section 6.1. |
| `simulation_reduced_basis.R` | Supplementary analysis of a reduced calibration basis under nonlinear selection, details can be found in supplementary section 6.2. |

All the analyses are performed with R 4.3.3.
