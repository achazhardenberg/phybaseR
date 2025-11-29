# phybaseR

**phybaseR** provides tools for specifying and fitting Phylogenetic Bayesian structural equation models (PhyBaSE) using [JAGS](http://mcmc-jags.sourceforge.net) via the `R2jags` package. 

> ⚠️ **Please beware that this package and its functions are still in beta.**  
> Functionality may change, and bugs may occur — feedback is very welcome!

## ✨ Functions

- **`phybase_model()`** - Automatically builds a PhyBaSE model in JAGS from a list of structural equations.
- **`phybase_run()`** - Runs the model using `R2jags`. The user only needs to provide the phylogenetic tree, the data (as a list), and a list with the structural equations.
- **`phybase_dsep()`** - Starting from the list of structural equations, calculates the implied conditional independencies and provides a new list of equations to test them. This list can then be fed to `phybase_run()` for testing.

## 📦 Installation

To install the package from GitHub, run the following command:

```r
# Install from GitHub (if using devtools or remotes)
remotes::install_github("achazhardenberg/phybaser")
```
## 🛠 Prerequisites

Before using phybaseR, you need to have JAGS (Just Another Gibbs Sampler) installed on your machine, as the package relies on it for Bayesian model fitting.

## 📥 Install JAGS

### macOS

You can install JAGS via Homebrew:
```bash
brew install jags
```
Or download the installer from the official website:
[http://mcmc-jags.sourceforge.net](http://mcmc-jags.sourceforge.net)

### Windows

Download and run the installer from:
[http://mcmc-jags.sourceforge.net](http://mcmc-jags.sourceforge.net)

### Linux (Debian/Ubuntu)

```bash
sudo apt-get install jags
```

Make sure to restart R after installing.

## 🧑‍💻 Example

Here’s an example workflow for using phybaseR:

```r
# Load required libraries
library(phybaseR)

# Load example data
data("rhino.tree")  # Example phylogenetic tree
data("rhino.dat")   # Example life-history data on Rhinograds

# Define the structural equations for the causal model (this is for model 8 in Gonzalez-Voyer & von Hardenberg (2013)
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

# Create a PhyBaSE model in the JAGS language
mod8.jg <- phybase_model(equations)  # Create a JAGS model for inspection or modification
cat(mod8.jg)  # Print the JAGS model for inspection

# Prepare data 
mod8.dat <- list(
  BM = rhino.dat$BM,
  LS = rhino.dat$LS,
  NL = rhino.dat$NL,
  DD = rhino.dat$DD,
  RS = rhino.dat$RS
)

# Run the PhyBaSE model
mod8.mcmc <- phybase_run(
  data = mod8.dat, 
  tree = rhino.tree, 
  equations = equations, 
  n.iter = 1000, 
  n.burnin = 500, 
  n.thin = 10, 
  n.chains = 3
)

# Print the summary of the MCMC object
mod8.mcmc

# Test for conditional independence
ind8.eq <- phybase_dsep(equations)  # Extract the minimum basis set of independence equations

# Create a JAGS model for testing conditional independencies
mod8.ind <- phybase_model(ind8.eq)
cat(mod8.ind)  # Print the JAGS model for inspection

# Run the PhyBaSE model with dsep=TRUE so it will provide only the betas needed to test for conditional independencies
mod8.ind.mcmc <- phybase_run(
  data = mod8.dat, 
  tree = rhino.tree, 
  equations = ind8.eq, 
  n.iter = 1000, 
  n.burnin = 500, 
  n.thin = 10, 
  n.chains = 3,
  dsep = TRUE
)

# Print the summary of the MCMC object for testing conditional independencies
mod8.ind.mcmc
```
## 📖 Citation

The implemented models are described in:

von Hardenberg, A. and Gonzalez-Voyer, A. (2025).
PhyBaSE: A Bayesian approach to Phylogenetic Structural Equation Models.
*Methods in Ecology and Evolution*. https://doi.org/10.1111/2041-210X.70044

Code for all models described in the paper is available at: [https://github.com/achazhardenberg/phybase](https://github.com/achazhardenberg/phybase)
