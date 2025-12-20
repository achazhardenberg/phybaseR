# 06. Automated Mediation Analysis

## Introduction

Causal inference often requires investigating *how* an effect occurs.
Does $X$ affect $Y$ directly, or does it work through a mediator $M$?

The `because` package provides a fully Bayesian automated mediation
analysis tool,
[`because_mediation()`](https://because-pkg.github.io/because/reference/because_mediation.md),
which decomposes the Total Effect of an exposure on an outcome into: 1.
**Direct Effect**: The effect of $\left. X\rightarrow Y \right.$ not
mediated by other variables in the graph. 2. **Indirect Effect(s)**: The
effect propagated through intermediate variables
($\left. X\rightarrow M\rightarrow Y \right.$).

This is calculated by multiplying the posterior distributions of
coefficients along each path, preserving full uncertainty
quantification.

## Example: Rhino Model (`sem8`)

Once again we will use the **Rhinogradentia** dataset to investigate the
causal drivers of **Dispersal Distance (DD)**.

### 1. Load Package and Data

``` r
library(because)

# Load package data
# Note: In a vignette context before installation, data() might be tricky,
# but assuming standard execution, these should be available.
data(rhino.dat)
data(rhino.tree)
```

### 2. Fit the Structural Equation Model

We define the `sem8` structure from Gonzalez-Voyer & von Hardenberg
(2014): **Body Mass (BM)** increases **Nose Length (NL)** which is also
influenced by Range Size (RS). \* **Nose Length (NL)** increases **Range
Size (RS)**. \* **Body Mass (BM)** *also* directly affects **Litter Size
(LS)**.

We want to investigate the effect of **Body Mass (BM)** on **Dispersal
Distance (DD)**. There is no direct link `BM -> DD` in the model
definition, so any effect MUST be mediated by **Nose Length (NL)**.

``` r
# Define the structural equations (sem8)
sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

# Fit the model
fit <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = rhino.tree,
    id_col = "SP",# Map species names
  n.iter = 1000  #We run a short run just for                        demonstration run it longer for                    real analyses!
)
#> Converted data.frame to list with 5 variables: LS, BM, NL, RS, DD
#> Standardizing tree (max_bt: 3.31 -> 1.0)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 300
#>    Unobserved stochastic nodes: 16
#>    Total graph size: 22546
#> 
#> Initializing model
```

### 3. Perform Mediation Analysis

We analyze the effect of **Body Mass (BM)** on **Dispersal Distance
(DD)**.

``` r
# Run Automated Mediation Analysis for BM -> DD
med_results <- because_mediation(fit, exposure = "BM", outcome = "DD")
```

#### Inspect the Summary

The summary table provides the posterior mean, standard deviation, and
credible intervals for the Total, Direct, and Total Indirect effects.

``` r
med_results$summary
#>                        Type     Mean         SD     Lower     Upper
#> 2.5%           Total Effect 0.252227 0.04780124 0.1730653 0.3454253
#> 2.5%1         Direct Effect 0.000000 0.00000000 0.0000000 0.0000000
#> 2.5%2 Total Indirect Effect 0.252227 0.04780124 0.1730653 0.3454253
```

**Interpretation:** \* **Total Effect**: The overall causal impact of BM
on DD. \* **Direct Effect**: The path `BM -> DD`. \* **Indirect
Effect**: The path `BM -> NL -> DD` (plus any others if they existed).

If the **Indirect Effect** is substantial and its credible interval
excludes zero, we have evidence of mediation.

#### Inspect Individual Paths

If there were multiple indirect paths (e.g., `BM -> LS -> NL -> RS`),
`because_mediation` would list each one separately.

``` r
med_results$paths
#>                Path     Type     Mean         SD     Lower     Upper
#> 2.5% BM -> NL -> DD Indirect 0.252227 0.04780124 0.1730653 0.3454253
```

This table shows exactly which causal pathways contribute to the
relationship.

### Technical Note

The function calculates the indirect effect as the product of
coefficients along the path:
$$\text{Indirect} = \beta_{BM\rightarrow NL} \times \beta_{NL\rightarrow DD}$$
This assumes linear relationships. For non-linear or generalized linear
models (e.g., Binomial), this product rule is an approximation of the
average causal effect on the linear predictor scale.
