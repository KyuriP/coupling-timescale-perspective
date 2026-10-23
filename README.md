<h1 align="center">Rethinking Group Differences in Psychopathology Networks</h1>
<p align="center"><em>The Role of Slowly Varying Contextual Conditions</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/language-R-276DC3?logo=r&logoColor=white" alt="Language: R">
  <img src="https://img.shields.io/badge/language-Python-3776AB?logo=python&logoColor=white" alt="Language: Python">
  <img src="https://img.shields.io/badge/license-MIT-blue" alt="License: MIT">
  <a href="https://osf.io/preprints/psyarxiv/u8wf5_v1"><img src="https://img.shields.io/badge/preprint-PsyArXiv-1B75BC" alt="Preprint: available on PsyArXiv"></a>
</p>

<p align="center">
  <img src="figs/slowfast.png" width="460" alt="Coupled slow–fast perspective linking contextual conditions and depressive symptoms">
</p>

This repository contains the analysis and figure-generation code for the manuscript:

> **Rethinking Group Differences in Psychopathology Networks: The Role of Slowly Varying Contextual Conditions**

The project develops a **slow–fast perspective** on group differences in psychopathology networks. Depressive symptoms are treated as a relatively fast-changing system, while social, economic, psychosocial, health-related, and lifestyle conditions form a more slowly varying contextual layer. The central question is whether group differences in depressive symptom expression are reflected in symptom–symptom interactions, symptom-activation parameters, or both.

## Study overview

The empirical illustration uses baseline data from the HELIUS cohort in Amsterdam and nine binary PHQ-9 symptom indicators. Contextual conditions are summarized in a composite measure called *Slow Risk Load*.

The manuscript uses continuous Slow Risk Load (SRL) in an exact moderated Ising model. Age, gender, and Dutch/non-Dutch background enter the symptom-activation equations; Dutch/non-Dutch background also enters the interaction equations. The primary analytic sample contains 23,270 participants. The likelihood is evaluated over all 512 configurations of the nine binary symptoms.

The code also covers the reduced-overlap SRL, the original SRL coding, separate and joint domain scores, a higher symptom-presence cutoff, ethnicity sensitivity analyses, exact median-split comparisons, a 10,000-permutation Network Comparison Test (NCT), and a four-category mixed graphical model (MGM) cross-check.

## Current analysis workflow

Run commands from the repository root, with authorized HELIUS data at `data/helius.sav`. The file must contain the original HELIUS variable names listed in `revision/revision_analyses.py`.

Install Python dependencies in a virtual environment:

```sh
python3 -m venv revision/.venv
source revision/.venv/bin/activate
python -m pip install -r revision/requirements.txt
python revision/revision_analyses.py --selftest
```

The self-test requires no HELIUS data. It checks derivatives and the equivalence between a binary group moderator and separate group models.

Run the manuscript analyses and create the local inputs for R:

```sh
python revision/revision_analyses.py --sections A,C,C0,C0b,Corig,B,Dses,Dstr,Dhlt,Djoint,E,Esplit,F,FDutch,FnonDutch,mgtables,export,rexport
```

Then run, in R:

```r
source("revision/nct_srl_split.R")
source("revision/mgm_ordinal_check.R")
source("revision/make_figures_current_srl.R")
```

The full NCT and MGM runs are slow. Setting `QUICK_TEST <- TRUE` in the MGM script runs a 3,000-person subsample and writes to `quick_test/`.

`--sections` accepts comma-separated section names listed at the top of the Python script. `rexport` needs the section C estimates: when both are requested, C runs first. `--sections all` also runs `G` (the two-mechanism example in Supplementary Table S11) and `figdata` (aggregate descriptive summaries).

## Files and outputs

| Location | Purpose |
|---|---|
| `revision/revision_analyses.py` | Current SRL construction and exact moderated/multigroup Ising analyses |
| `revision/nct_srl_split.R` | NCT for the SRL median split |
| `revision/mgm_ordinal_check.R` | MGM cross-checks (0-3 items and binary items) |
| `revision/make_figures_current_srl.R` | Manuscript figures (runs the plotting scripts on the current results) |
| `results/revision_analyses/` | Current aggregate estimates, diagnostics and summaries |
| `results/revision_analyses/nct/` | Current NCT summaries and network matrices |
| `results/figure_data/` | Aggregate data used by figure scripts |
| `figs/revised/*current*` | Current descriptive and supplementary figures |
| `figs/revised/figure_04_srl_activation_effects.*` | SRL effects on activation and interactions |
| `R/` | Original SRL workflow and shared R setup/helpers |
| `fig_scripts/revised/` | Plotting templates reused by the current wrapper |
| `archive/` | Earlier and exploratory code/aggregate outputs; not the current analysis entry point |

The conceptual figures (`figs/three_panel_network_analysis_infographic.png`, `figs/slowfast.png`) are not generated by code.

The scripts `R/01_prepare_analysis_data.R` to `R/07_build_appendix_outputs.R` use the alternative SRL coding and write to `results/network_invariance_exact/` and `results/appendix/`. The setup and plotting code in `R/` is also used by the figure scripts.

## Interpretation and reproducibility

- The exact Ising model uses the {0,1} parameterization without a separate inverse-temperature parameter. An activation coefficient is a conditional log-odds effect when all other symptoms are absent; model-implied marginal probabilities combine activation and interaction terms.
- Primary SRL construction uses 15 indicators across three domains. A domain needs at least half its indicators observed, and SRL needs at least two domains. Work/labour indicators use binary adversity coding; activity and pack-years are winsorized at the 1st/99th percentiles. See `prepare()` for the exact recodes.
- The low/high SRL cutoff is the median among participants with computable SRL and PHQ total, before the covariate and item-completeness exclusions. Ties enter the low group. This yields 11,714 low and 11,556 high participants in the analytic sample.
- NCT uses 10,000 permutations, gamma = 0.25, the AND rule, and Holm adjustment over 36 edge tests. Unadjusted edge p-values use `(b + 1)/(B + 1)`, as in the package.
- In the MGM, PHQ-9 items are categorical nodes with four levels (the order is not used). `weightsAgg` values are absolute weights without sign.
- Random seeds are set in the scripts.

## Data availability

Individual-level HELIUS data are subject to participant privacy and study-governance requirements. Researchers may apply through the HELIUS study's established data-access procedures. Only code and aggregate results are intended for distribution.

Raw data and participant-level exports belong in the ignored `data/` directory. R model objects (`.rds`, `.RData`, `.rda`) can contain the data and are not tracked.

## Software

Python dependencies are listed in `revision/requirements.txt`. The current R workflow uses `mgm`, `NetworkComparisonTest`, and the plotting/setup packages:

```text
dplyr, tidyr, stringr, ggplot2, forcats, patchwork, ggrepel,
purrr, broom, haven, skimr, bootnet, IsingSampler, qgraph,
ggpubr, scales, ggtext
```

The appendix scripts also use `knitr`. Analyses were run with NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.2 and pyreadstat, and in R with `mgm` 1.2.15 and `NetworkComparisonTest` 2.2.2.

## Citation

The preprint is available on PsyArXiv:

> Park, K. et al. Rethinking Group Differences in Psychopathology Networks: A Slow-Fast Perspective on Context and Symptom Activation. *PsyArXiv*.
> https://osf.io/preprints/psyarxiv/u8wf5_v1

To cite this repository specifically:

```text
Park, K. et al. Rethinking Group Differences in Psychopathology Networks:
The Role of Slowly Varying Contextual Conditions.
GitHub repository: https://github.com/KyuriP/coupling-timescale-perspective
```
