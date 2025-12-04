# phybaseR

**phybaseR** provides tools for specifying and fitting Phylogenetic Bayesian structural equation models (PhyBaSE) using [JAGS](http://mcmc-jags.sourceforge.net).

> ⚠️ **Please beware that this package and its functions are still in beta.**
> Functionality may change, and bugs may occur — feedback is very welcome!

## ✨ Features

**phybaseR** simplifies the process of running complex phylogenetic path analyses by automatically generating JAGS code from standard R formulas. Key features include:

-   **Automatic JAGS Code Generation**: Builds models directly from a list of structural equations.
-   **Phylogenetic Uncertainty**: Incorporates uncertainty by sampling across a set of trees.
-   **Missing Data Support**: Handles missing values in both response and predictor variables using phylogenetic imputation.
-   **Measurement Error**: Accounts for within-species variability or measurement error.
-   **Multiple Response Distributions**: Supports Gaussian, Binomial, Multinomial, Ordinal, Poisson, and Negative Binomial distributions.
-   **Categorical Predictors**: Automatic handling of factor variables with dummy variable expansion.
-   **Latent Variables**: Support for modeling induced correlations from latent common causes.
-   **Model Validation**: Tools for d-separation testing to validate causal hypotheses.
-   **Optimized Performance**: Fast random effects formulation (4.6× speedup) for large datasets.
-   **Parallel Computing**: Run MCMC chains in parallel on multi-core systems for faster computation.


## 📦 Installation

To install the package from GitHub, run the following command:

```r
# Install from GitHub (if using devtools or remotes)
remotes::install_github("achazhardenberg/phybaseR", build_vignettes = TRUE)
```

## 🛠 Prerequisites

Before using phybaseR, you need to have **JAGS** (Just Another Gibbs Sampler) installed on your machine.

-   **macOS**: `brew install jags` or download from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Windows**: Download installer from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Linux**: `sudo apt-get install jags`.

## � Documentation & Examples

For detailed examples and a step-by-step tutorial on how to use all features of the package, please refer to the package vignette:

```r
vignette("phybaseR_tutorial", package = "phybaseR")
```

The vignette covers:
-   Basic model specification and running
-   Handling missing data
-   Incorporating measurement error
-   Multiple response distributions (Gaussian, Binomial, Multinomial, Ordinal, Poisson, Negative Binomial)
-   Categorical predictors with automatic dummy coding
-   Using latent variables
-   Causal model validation with d-separation
-   Optimized performance mode for large datasets
-   Parallel computing for faster MCMC


## 📖 Citation

If you use **phybaseR** in your research, please cite both the package and the methodological paper:

**Package:**
> von Hardenberg, A. and Gonzalez-Voyer, A. (2025). phybaseR: An R package to easily create and run PhyBaSE models. R package version 1.0. https://github.com/achazhardenberg/phybaseR

**Methodological Paper:**
> von Hardenberg, A. and Gonzalez-Voyer, A. (2025). PhyBaSE: A Bayesian approach to Phylogenetic Structural Equation Models. *Methods in Ecology and Evolution*. https://doi.org/10.1111/2041-210X.70044

