# `inference`: A Python Package for Statistical Inference

## Overview

The `inference` package is a collection of Python modules implementing a variety of methods targeting the statistical inference problems—and the statistical modeling *style*—of the physical sciences. Our own discipline is astronomy, and our choice of problems and methods most directly targets the needs of astronomers, but many tools here may be of use to other physical scientists.

We adopt a narrow view of statistical inference for this package: we strive to provide tools that not only provide a "best" estimate or choice in a problem, but that also _**quantify uncertainties**_. For many scientific purposes, merely knowing what's "best" is not enough; we also need to know the range of possibilities that are *nearly as good* (and thus plausible). That's the kind of data science the `inference` package aims to address.

## A Library and a Framework

The modules comprising the `inference` package fall into two categories:

- **Library modules:** These modules implement tools that are largely self-contained. They fall into three categories:

  - *Statistical tools:* These modules implement statistical methods targeting specific, focused problems, such as estimating rates from Poisson-distributed event counts, or fitting a line to data with errors in both abscissa and ordinate.
  - *Signal models:* These modules implement functions and classes for calculating predictions from selected widely-used models in astronomy and physics, such as Keplerian models for orbital motion, or cosmological models for the distribution of luminous sources in space-time.
  - *Utilities:* These modules implement algorithms that are useful for statistical computation, but which may be useful in other settings; examples include multi-dimensional numerical integration and Monte Carlo algorithms.

- **Parametric Inference Engine (PIE):** These modules comprise a framework facilitating exploring the parameter spaces of statistical models for data, for three different general parametric inference paradigms: *chi-squared* (more accurately, weighted least squares), *maximum likelihood*, and *Bayesian*.

For more information about `inference`, see our (preliminary and evolving) [documentation site](http://inference.astro.cornell.edu/inference/index.html).

