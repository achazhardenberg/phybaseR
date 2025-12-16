# because

### Bayesian Estimation of Causal Effects ![](reference/figures/because_hex.svg)

**because** provides a unified framework for specifying and fitting
Bayesian structural equation models using
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

To install the **stable release** (`v0.9.5`), run:

``` r
remotes::install_github("achazhardenberg/because@v0.9.5", build_vignettes = TRUE)
```

To install the **latest development version** (unstable), run:

``` r
remotes::install_github("achazhardenberg/because", build_vignettes = TRUE)
```

## Prerequisites

Before using phybaseR, you need to have **JAGS** (Just Another Gibbs
Sampler) installed on your machine.

- **macOS**: `brew install jags` or download from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Windows**: Download installer from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Linux**: `sudo apt-get install jags`.

## Documentation & Examples

For detailed examples and step-by-step tutorials on how to use all
features of the package, please refer to the package vignettes:

``` r
vignette("01_getting_started", package = "because")
```

You can also download the full **[Reference Manual
(PDF)](https://achazhardenberg.github.io/because/inst/pdf/because_reference_manual.pdf)**
which describes all functions in the package.

**PDF Vignettes:** - [Getting
Started](https://achazhardenberg.github.io/because/inst/pdf/01_getting_started.pdf) -
[D-Separation
Tests](https://achazhardenberg.github.io/because/inst/pdf/02_dseparation.pdf) -
[Phylogenetic
Models](https://achazhardenberg.github.io/because/inst/pdf/03_phylogenetic_models.pdf)

## Citation

To cite the **because** package in publications, please use:

> von Hardenberg, A. (2025). because: Bayesian Estimation of Causal
> Effects. R package version 0.9.4.
> <https://github.com/achazhardenberg/because>

If you use **because** to run Phylogenetic Bayesian Structural Equation
Models (PhyBaSE) in your research, please *also* cite:

> von Hardenberg, A. and Gonzalez-Voyer, A. (2025). PhyBaSE: A Bayesian
> approach to Phylogenetic Structural Equation Models. *Methods in
> Ecology and Evolution*. <https://doi.org/10.1111/2041-210X.70044>
