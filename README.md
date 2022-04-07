# stbench

Functions to fit various spatio-temporal GLMMs in Template Model Builder (TMB), and produce fully Bayesian benchmarked estimates. Additional functionality includes various benchmarking approaches applied to posterior samples.

## Overview

**stbench** is an open-source R package for producing fully Bayesian benchmarked estimates for hierarchical models in TMB. The package currently supports spatial and spatio-temporal modeling of binomial count data of under-5 mortality (U5MR), spatial models for generic binomial count data, and spatial models for normally distributed outcomes. 

Current model implementations:

* U5MR
  * Single-survey 
    * Space-only
      * Binomial benched/unbenched
    * Time-only
      * Binomial benched/unbenched
    * Space-time
      * Binomial benched/unbenched
      * BetaBinomial benched/unbenched
  * Multi-survey
    * Space-time
      * Binomial benched/unbenched

* Generic binary outcomes
  * Single survey
    * Space-only
      * Binomial benched/unbenched

* Normally distributed outcomes
  * Single survey
    * Space-only
      * benched/unbenched
      
Details can be found in the functions `fit_u5mr`, `fit_binary`, and `fit_normal`, respectively.

Current functions available that apply benchmarking methods to samples from a distribution include:

* `benchmark_sampler`: Fully Bayesian benchmarking via a rejection sampler, described in [Okonek and Wakefield, 2022](https://arxiv.org/abs/2203.12195)

* `benchmark_bayesest`: Constrained Bayes estimate approach to benchmarking, described in [Datta et al., 2011](https://doi.org/10.1007/s11749-010-0218-y)

## Installation

The current development version can be installed using `devtools::install_github()`:

```R
devtools::install_github(repo="taylorokonek/stbench")
```


