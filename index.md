# because

### Bayesian Estimation of Causal Effects

**because** provides a unified framework for specifying and fitting
Bayesian structural equation models in `R` using
[JAGS](http://mcmc-jags.sourceforge.net).

> **Please beware that this package and its functions are still in
> beta.** Functionality may change, and bugs may occur — feedback is
> very welcome!

## Features

**because** simplifies the process of running complex Bayesian
Structural Equation Models by automatically generating JAGS code from
standard R formulas. Key features include:

- **Automatic JAGS Code Generation**: Builds models directly from a list
  of structural equations.
- **Generalized Covariance Structures**: Supports standard SEMs,
  Phylogenetic SEMS (PhyBaSE) as well as Spatial SEMs and Animal Models
  via custom covariance matrices or pedigrees.
- **Phylogenetic Uncertainty**: Incorporates uncertainty by sampling
  across a set of trees.
- **Missing Data Support**: Handles missing values in both response and
  predictor variables using phylogenetic imputation.
- **Measurement Error**: Accounts for within-species variability or
  measurement error.
- **Multiple Response Distributions**: Supports Gaussian, Binomial,
  Multinomial, Ordinal, Poisson, and Negative Binomial distributions.
- **Categorical Predictors**: Automatic handling of factor variables
  with dummy variable expansion.
- **Latent Variables**: Support for modeling induced correlations from
  latent common causes.
- **Model Validation**: Tools for d-separation testing to validate
  causal hypotheses.
- **Parallel Computing**: Run MCMC chains in parallel on multi-core
  systems for faster computation.

## Installation

To install the **stable release** (`v0.9.6`), run:

``` r
remotes::install_github("because-pkg/because@v0.9.6", build_vignettes = TRUE)
```

To install the **latest development version** (unstable), run:

``` r
remotes::install_github("because-pkg/because", build_vignettes = TRUE)
```

## Prerequisites

Before using phybaseR, you need to have **JAGS** (Just Another Gibbs
Sampler) installed on your machine.

- **macOS**: `brew install jags` or download from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Windows**: Download installer from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Linux**: `sudo apt-get install jags`.
