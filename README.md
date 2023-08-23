# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging

### Bonnie E. Shook-Sa, Paul N. Zivich, Samuel P. Rosin, Jessie K. Edwards, Adaora A. Adimora, Michael G. Hudgens, and Stephen R. Cole

**Citation**: Shook-Sa BE, Zivich PN, Rosin SP, Edwards JK, Adimora AA, Hudgens MG, and Cole SR. "Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging." *arXiv* arXiv:2305.00845.
[![arXiv](https://img.shields.io/badge/arXiv-2305.00845-b31b1b.svg)](https://arxiv.org/abs/2305.00845)
--------------------------------

## Abstract

While randomized controlled trials (RCTs) are critical for establishing the efficacy of new therapies, there are limitations regarding what comparisons can be made directly from trial data. RCTs are limited to a small number of comparator arms and often compare a new therapeutic to a standard of care which has already proven efficacious. It is sometimes of interest to estimate the efficacy of the new therapy relative to a treatment that was not evaluated in the same trial, such as a placebo or an alternative therapy that was evaluated in a different trial. Such dual-study comparisons are challenging because of potential differences between trial populations that can affect the outcome. In this paper, two bridging estimators are considered that allow for comparisons of treatments evaluated in different trials, accounting for measured differences in trial populations. A "multi-span" estimator leverages a shared arm between two trials, while a "single-span" estimator does not require a shared arm. A diagnostic statistic that compares the outcome in the standardized shared arms is provided. The two estimators are compared in simulations, where both estimators demonstrate minimal empirical bias and nominal confidence interval coverage when the identification assumptions are met. The estimators are applied to data from the AIDS Clinical Trials Group 320 and 388 to compare the efficacy of two-drug versus four-drug antiretroviral therapy on CD4 cell counts among persons with advanced HIV. The single-span approach requires fewer identification assumptions and was more efficient in simulations and the application.

--------------------------------

## File Manifesto

### Python
The `Python/` path Python code to replicate the simulation study presented in the paper.

### R
The `R/` path contains R programs to replicate the simulatoin study presented in the paper.
