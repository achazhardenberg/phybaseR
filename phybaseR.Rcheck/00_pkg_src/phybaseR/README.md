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
## � Measurement Error / Variability

You can account for measurement error or within-species variability by providing standard errors for your traits.
To do this:
1. Ensure your data list includes the standard error for the variable (suffixed with `_se`, e.g., `BM_se`).
2. Specify the variable name in the `variability` argument of `phybase_run()`.

```r
# Example with measurement error on Body Mass (BM)
# Add standard errors to the data
mod8.dat$BM_se <- rep(0.1, 100) # Example SE

# Run the model specifying variability
mod8.me <- phybase_run(
  data = mod8.dat,
  tree = rhino.tree,
  equations = equations,
  variability = "BM" # Treat BM as having measurement error
)
```

### Repeated Measures

If you have repeated measures for a trait (e.g., multiple individuals per species), you can provide the raw data as a matrix (rows = species, columns = replicates). The model will estimate the measurement error variance from the data.

1. Provide the data as a matrix in your data list.
2. Specify the variable name in `variability`.

```r
# Example with repeated measures for Body Mass (BM)
# mod8.dat$BM is a matrix of size N_species x N_replicates
mod8.reps <- phybase_run(
  data = mod8.dat, # Contains BM as a matrix
  tree = rhino.tree,
  equations = equations,
  variability = "BM" # Inferred as 'reps' type because BM is a matrix
)
```

You can also explicitly specify the type if needed:
```r
phybase_run(..., variability = c(BM = "reps", LS = "se"))
```

> **Note**: This feature is fully compatible with phylogenetic uncertainty. If you provide a list of trees (e.g., `multiPhylo`) to `phybase_run()`, the model will account for both measurement error and phylogenetic uncertainty simultaneously.

## ✅ Missing Data Support

PhyBaSE **fully supports missing data** (NA values) in both response and predictor variables.

### How it works
When missing data is detected:
1.  **Response Variables**: PhyBaSE automatically switches to a **Latent Variable (GLMM)** formulation for that variable. This preserves phylogenetic signal during imputation by modeling the error term as:
    $$ Y = \mu + \epsilon_{phylo} + \epsilon_{residual} $$
    where $\epsilon_{phylo}$ tracks the phylogeny and $\epsilon_{residual}$ is independent error.

2.  **Predictor-Only Variables**: If a root predictor has missing values, PhyBaSE automatically assigns it a prior (by adding `Variable ~ 1` to the model) and uses the same GLMM formulation to impute missing values while accounting for phylogenetic signal.

### Example
```r
# Data with NAs in both predictors and responses
data_list <- list(
  X = c(2.1, NA, 2.3, 1.9, 2.0),   # Missing in predictor
  Y = c(1.2, 0.9, NA, 1.4, 1.5)    # Missing in response
)

# Just run it!
fit <- phybase_run(
  data = data_list, 
  tree = tree, 
  equations = list(Y ~ X)
)

# PhyBaSE will:
# 1. Detect missing X -> Add "X ~ 1" automatically
# 2. Use GLMM imputation for both X and Y
# 3. Estimate parameters using the full dataset
```

> **Note**: For variables with missing data, the model estimates `lambda` (phylogenetic signal) derived from the variance components of the latent variable model.


## 🎲 Binomial Variables

PhyBaSE supports binary response variables (0/1) using a phylogenetic logistic regression approach. This is useful for modeling traits like presence/absence or behavioral states.

The model uses:
- A **logit link function** to relate predictors to the probability of success
- A **Bernoulli likelihood** for the observed binary data  
- **Phylogenetic correlation** modeled on the latent scale via a random effect

**Important**: Binomial variables should be child nodes only in your DAG (i.e., they should not be used as predictors for other variables).

```r
# Example: Gregariousness (0/1) depends on Body Mass
data_list <- list(
  BM = log_body_mass,
  Gregarious = binary_gregarious  # 0 or 1
)

equations <- list(Gregarious ~ BM)

# Specify the distribution
fit <- phybase_run(
  data = data_list,
  tree = tree,
  equations = equations,
  distribution = c(Gregarious = "binomial")
)
```

## 📖 Citation

The implemented models are described in:

von Hardenberg, A. and Gonzalez-Voyer, A. (2025).
PhyBaSE: A Bayesian approach to Phylogenetic Structural Equation Models.
*Methods in Ecology and Evolution*. https://doi.org/10.1111/2041-210X.70044

Code for all models described in the paper is available at: [https://github.com/achazhardenberg/phybase](https://github.com/achazhardenberg/phybase)
