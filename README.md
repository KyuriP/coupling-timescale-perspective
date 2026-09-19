<h1 align="center">Rethinking Group Differences in Psychopathology Networks</h1>
<p align="center"><em>The Role of Slowly Varying Contextual Conditions</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/language-R-276DC3?logo=r&logoColor=white" alt="Language: R">
  <img src="https://img.shields.io/badge/license-MIT-blue" alt="License: MIT">
  <img src="https://img.shields.io/badge/status-preprint%20coming%20soon-yellow" alt="Status: preprint coming soon">
</p>

<p align="center">
  <img src="figs/revised/figure_01_slow_fast_concept.png" width="600" alt="Slower contextual conditions may shift symptom activation and/or symptom-symptom interactions in the fast symptom layer">
</p>

This repository contains the analysis and figure-generation code for the manuscript:

> **Rethinking Group Differences in Psychopathology Networks: The Role of Slowly Varying Contextual Conditions**

The project develops a **slow–fast perspective** on group differences in psychopathology networks. Depressive symptoms are treated as a relatively fast-changing system, while social, economic, psychosocial, health-related, and lifestyle conditions form a more slowly varying contextual layer. The central question is whether group differences in depressive symptom expression are reflected in symptom–symptom interactions, symptom-activation parameters, or both.

## Study overview

The empirical illustration uses baseline data from the HELIUS cohort in Amsterdam and nine binary PHQ-9 symptom indicators. Contextual conditions are summarized in a composite measure called *Slow Risk Load*.

The analysis includes:

- descriptive comparisons of depressive symptom levels and Slow Risk Load;
- a 10,000-permutation Network Comparison Test;
- exact multigroup Ising models in the \(\{0,1\}\) parameterization;
- separate tests of interaction and activation invariance;
- Wald confidence intervals for group differences in symptom activation;
- symptom-specific logistic regressions using continuous Slow Risk Load;
- appendix diagnostics and sensitivity analyses.

Because there are only nine binary symptoms, the Ising partition function can be evaluated exactly over all \(2^9=512\) possible symptom configurations.

## Repository structure

```text
coupling-timescale-perspective/
├── R/                         # Data preparation and statistical analyses
├── fig_scripts/
│   └── revised/               # Current figure-generation scripts (01-05, A1)
├── figs/
│   └── revised/               # Current manuscript figures — everything else in figs/ is an earlier iteration, kept for provenance
├── results/
│   ├── network_invariance_exact/
│   │   └── slow_risk/
│   │       ├── exact_multigroup_ising/
│   │       └── nct_10000/
│   ├── appendix/              # Appendix tables and figures
│   └── figure_data/           # Data used to generate figures
├── archive/                   # Superseded and exploratory code
└── coupling-timescale-perspective.Rproj
```

The current definitive multigroup analysis is stored under:

```text
results/network_invariance_exact/
```

Older exploratory or package-specific implementations are retained in `archive/` for provenance but are not part of the manuscript pipeline.

## Analysis workflow

Run all scripts from the repository root. The main analysis order is:

```r
source("R/00_setup.R")
source("R/01_prepare_analysis_data.R")
source("R/03_prepare_binary_data.R")
source("R/04_network_invariance_exact.R")
source("R/05_nct_final_10000.R")
source("R/06_activation_difference_ci.R")
source("R/07_build_appendix_outputs.R")
```

`R/02_model_helpers.R` contains helper functions and is sourced by the analysis scripts.

The main outputs are written to:

```text
results/network_invariance_exact/slow_risk/exact_multigroup_ising/
results/network_invariance_exact/slow_risk/nct_10000/
results/appendix/
```

## Data availability

The individual-level HELIUS data are *not* included in this repository because of participant privacy and study-governance requirements. Researchers may apply for access through the HELIUS study's established data-access procedures.

The repository contains code and non-sensitive derived outputs needed to document the analyses and reproduce manuscript figures once authorized data are available in the expected local format.


## Software

The analyses are written in R. Core packages used across the pipeline include:

```text
dplyr, tidyr, stringr, ggplot2, forcats, patchwork, ggrepel,
purrr, broom, haven, skimr, bootnet, IsingSampler, qgraph,
ggpubr, NetworkComparisonTest, knitr
```

Install any missing packages before running the pipeline. Package versions will be documented at preprint release.

## Reproducibility notes

- Run scripts from the repository root.
- The exact Ising analysis uses the \(\{0,1\}\) parameterization without a separate inverse-temperature parameter.
- The final NCT uses 10,000 permutations, \(\gamma=0.25\), the AND rule, and Holm adjustment across 36 edge-specific tests.
- Random seeds are fixed within the analysis scripts.
- Generated outputs should not be edited manually unless explicitly noted in the manuscript workflow.

## Citation

A formal citation and preprint link will be added when the manuscript is posted.

For now, please cite this repository as:

```text
Park, K. et al. Rethinking Group Differences in Psychopathology Networks:
The Role of Slowly Varying Contextual Conditions.
GitHub repository: https://github.com/KyuriP/coupling-timescale-perspective
```
