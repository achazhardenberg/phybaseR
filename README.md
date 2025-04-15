# phybaseR

**phybaseR** provides tools for specifying and fitting Phylogenetic Bayesian structural equation models (PhyBaSE) using [JAGS](http://mcmc-jags.sourceforge.net) via the `R2jags` package.

## ✨ Functions

- **`phybase_model()`** - Automatically builds a PhyBaSE model in JAGS from a list of structural equations.
- **`phybase_run()`** - Runs the model using `R2jags`. The user only needs to provide the phylogenetic tree, the data (as a list), and a list with the structural equations.
- **`phybase_dsep()`** - Starting from the list of structural equations, calculates the implied conditional independencies and provides a new list of equations to test them. This list can then be fed to `phybase_run()` for testing.

## 📦 Installation

To install the package from GitHub, run the following command:

```r
# Install from GitHub (if using devtools or remotes)
remotes::install_github("achazhardenberg/phybaser")

