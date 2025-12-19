# Custom Priors and Mechanistic Modeling

## Introduction

One of the strengths of Bayesian modeling is the ability to incorporate
prior knowledge into your analysis. In the `because` package, we
generally use “uninformative” or weakly informative priors by default to
let the data speak for itself. However, there are many cases where you
might want to inject specific information:

1.  **Mechanistic Knowledge**: You know a parameter (like a growth rate)
    must be positive or lies within a biologically feasible range.
2.  **Previous Studies**: You have estimates from a previous
    meta-analysis or experiment that you want to use as a starting
    point.
3.  **Measurement Error**: You know the measurement error of an
    instrument and want to fix or constrain the residual variance.

This vignette demonstrates how to use the `priors` argument in
[`because()`](https://because-pkg.github.io/because/reference/because.md)
to override default priors and inject your own mechanistic constraints.

## Example 1: Comparing Default and Custom Priors

For this example, let’s simulate **Absolute Growth Rate** in **grams per
day (g/day)** for hypothetical juvenile birds responding to ambient
temperature.

``` r
library(because)

set.seed(123)
N <- 100
# Ambient temperature in Celsius (Mean 20C)
Temp_Raw <- rnorm(N, mean = 20, sd = 5)

# Center it so '0' represents the mean (20C).
# This ensures the intercept represents growth at average conditions,
# rather than at 0°C (where the bird would be frozen solid!).
Temp_Centered <- Temp_Raw - 20

# True relationship:
# Baseline growth at mean temp is 10 g/day.
# Every degree of warming above mean adds 0.5 g/day.
Growth_g_day <- 0.5 * Temp_Centered + 10 + rnorm(N, sd = 2)

df <- data.frame(Temp_Centered, Growth_g_day)
```

### Default Model

Let’s fit a standard model first. By default, `because` assigns wide
Gaussian priors (`dnorm(0, 1.0E-6)`) to intercepts (`alpha`) and slopes
(`beta`), and Gamma priors (`dgamma(1, 1)`) to precisions (`tau`).

``` r
fit_default <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df
)
#> No tree provided. Running standard (non-phylogenetic) SEM.
#> Converted data.frame to list with 2 variables: Growth_g_day, Temp_Centered
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 10509
#> 
#> Initializing model

summary(fit_default)
#>                                  Mean    SD Naive SE Time-series SE  2.5%   50%
#> alphaGrowth_g_day               9.793 0.196    0.004          0.004 9.401 9.793
#> beta_Growth_g_day_Temp_Centered 0.479 0.043    0.001          0.001 0.393 0.479
#> sigmaGrowth_g_day               1.940 0.138    0.003          0.003 1.692 1.933
#>                                  97.5% Rhat n.eff
#> alphaGrowth_g_day               10.178    1  2912
#> beta_Growth_g_day_Temp_Centered  0.568    1  3187
#> sigmaGrowth_g_day                2.219    1  2902
#> 
#> DIC:
#> Mean deviance:  417.4 
#> penalty 3.02 
#> Penalized deviance: 420.5
```

### Custom Prior Model

Now, suppose we have strong prior knowledge from scaling theory or
previous literature: 1. **Intercept (Alpha)**: At **mean temperatures**
(20°C), metabolic constraints suggest growth should be around **10
g/day**. 2. **Slope (Beta)**: We expect a positive effect of
temperature, likely around **0.5 g/day/°C**.

The parameter names usually follow this convention:

**Intercepts**: `alphaResponseVariable` (e.g., `alphaGrowth_g_day`)

**Slopes**: `beta_Response_Predictor` (e.g.,
`beta_Growth_g_day_Temp_Centered`)

**Residual precision:** `tau_e_Response`

``` r
# Define our custom priors
my_priors <- list(
    # Strong prior on intercept: Mean 10, Precision 100 (SD = 0.1)
    alphaGrowth_g_day = "dnorm(10, 100)",

    # Informative prior on slope: Mean 0.5, Precision 400 (SD = 0.05)
    beta_Growth_g_day_Temp_Centered = "dnorm(0.5, 400)"
)

fit_custom <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    priors = my_priors
)
#> No tree provided. Running standard (non-phylogenetic) SEM.
#> Converted data.frame to list with 2 variables: Growth_g_day, Temp_Centered
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 10511
#> 
#> Initializing model

summary(fit_custom)
#>                                  Mean    SD Naive SE Time-series SE  2.5%   50%
#> alphaGrowth_g_day               9.957 0.090    0.002          0.002 9.780 9.956
#> beta_Growth_g_day_Temp_Centered 0.485 0.033    0.001          0.001 0.420 0.485
#> sigmaGrowth_g_day               1.934 0.133    0.002          0.002 1.697 1.924
#>                                  97.5% Rhat n.eff
#> alphaGrowth_g_day               10.135    1  2818
#> beta_Growth_g_day_Temp_Centered  0.549    1  3000
#> sigmaGrowth_g_day                2.225    1  3150
#> 
#> DIC:
#> Mean deviance:  417 
#> penalty 1.804 
#> Penalized deviance: 418.8
```

Notice how the credible intervals for the custom model will be tighter
and centered closer to our priors, especially if the data were sparse
(small N).

### Visualizing the Impact

We can plot the posterior estimates side-by-side to see the “shrinkage”
effect of our informative priors.

``` r
# We can use the helper function plot_posterior() to compare models.
# By passing a list of models, they are overlaid on the same plot.
# The 'parameter' argument uses partial matching (regex), so "Growth_g_day"
# will match both the intercept (alphaGrowth_g_day) and the slope (beta_Growth_g_day...).

plot_posterior(
    model = list("Default" = fit_default, "Custom" = fit_custom),
    parameter = "Growth_g_day"
)
```

![](05_custom_priors_files/figure-html/plot_compare-1.png)

## Example 2: Mechanistic Constraints

In many biological contexts, a negative effect is not just unlikely, it
is **physically impossible**.

For example, consider the relationship between **Body Mass** and
**Metabolic Rate** (Kleiber’s Law). A larger animal simply cannot
consume *less* energy than a smaller one (a negative slope), all else
equal. The physics of life requires costs to scale positively with mass.

However, if our sample size is small and measurement error is huge, we
might accidentally estimate a negative slope. We can prevent this by
enforcing a positive prior.

``` r
# Simulate data following Kleiber's Law: MR = a * Mass^0.75
# Taking logs: log(MR) = log(a) + 0.75 * log(Mass)
set.seed(42)
Mass <- runif(30, 10, 100)
# True scaling exponent is 0.75
# We add enough noise that the estimated slope might be negative by chance
Log_Mass <- log(Mass)
Log_MR <- 1.5 + 0.75 * Log_Mass + rnorm(30, sd = 2.5)

df_kleiber <- data.frame(Log_Mass, Log_MR)

# Prior: The scaling exponent (slope) must be positive
# Physics dictates that metabolic cost increases with mass.
positive_prior <- list(
    beta_Log_MR_Log_Mass = "dnorm(0, 1) T(0, )"
)

# 1. Fit Default Model (Unconstrained) - Log-Log regression
fit_default_kleiber <- because(
    equations = list(Log_MR ~ Log_Mass),
    data = df_kleiber
)
#> No tree provided. Running standard (non-phylogenetic) SEM.
#> Converted data.frame to list with 2 variables: Log_MR, Log_Mass
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 1059
#> 
#> Initializing model

# 2. Fit Mechanistic Model (Constrained)
fit_mech_kleiber <- because(
    equations = list(Log_MR ~ Log_Mass),
    data = df_kleiber,
    priors = positive_prior
)
#> No tree provided. Running standard (non-phylogenetic) SEM.
#> Converted data.frame to list with 2 variables: Log_MR, Log_Mass
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 30
#>    Unobserved stochastic nodes: 3
#>    Total graph size: 1059
#> 
#> Initializing model

# Compare Estimates
summary(fit_default_kleiber)
#>                       Mean    SD Naive SE Time-series SE   2.5%   50%  97.5%
#> alphaLog_MR          3.471 4.246    0.078          0.269 -5.066 3.533 11.793
#> beta_Log_MR_Log_Mass 0.051 1.038    0.019          0.065 -2.002 0.033  2.074
#> sigmaLog_MR          2.932 0.395    0.007          0.008  2.279 2.885  3.814
#>                       Rhat n.eff
#> alphaLog_MR          1.003   250
#> beta_Log_MR_Log_Mass 1.003   253
#> sigmaLog_MR          1.000  2638
#> 
#> DIC:
#> Mean deviance:  151.1 
#> penalty 3.128 
#> Penalized deviance: 154.3
summary(fit_mech_kleiber)
#>                       Mean    SD Naive SE Time-series SE   2.5%   50% 97.5%
#> alphaLog_MR          1.151 1.917    0.035          0.059 -3.237 1.476 3.925
#> beta_Log_MR_Log_Mass 0.623 0.456    0.008          0.014  0.030 0.533 1.679
#> sigmaLog_MR          2.912 0.382    0.007          0.007  2.271 2.879 3.800
#>                      Rhat n.eff
#> alphaLog_MR             1  1045
#> beta_Log_MR_Log_Mass    1  1011
#> sigmaLog_MR             1  2918
#> 
#> DIC:
#> Mean deviance:  150.5 
#> penalty 2.23 
#> Penalized deviance: 152.7

# Visualize: Unconstrained vs. Truncated
plot_posterior(
    list(Unconstrained = fit_default_kleiber, Mechanistic = fit_mech_kleiber),
    parameter = "beta_Log_MR_Log_Mass",
    density_args = list(Mechanistic = list(from = 0))
)
```

![](05_custom_priors_files/figure-html/mechanistic_constraint-1.png)

## Example 3: Informing Residual Variance components

We usually estimate the residual variance (`sigma` or its inverse `tau`)
from the data. However, you might have prior knowledge about the
expected noise level.

> \[!NOTE\] **Distinction from `variability` argument**: The
> `variability` argument in
> [`because()`](https://because-pkg.github.io/because/reference/because.md)
> is designed for **Measurement Error** in predictors
> (Error-in-Variables) or when you have repeated measures per
> individual.
>
> Custom priors on `tau_e`, shown below, are for providing information
> about the **Residual Variance** of the response variable (which
> includes both Process Error and unrecognized Observation Error).

For example, if you know the precision of your scale is roughly 2 g,
this gives you a strong expectation for the residual standard deviation
(sigma).

``` r
# We believe sigma is around 2 (grams/day).
# tau = 1/2^2 = 0.25.
# A Gamma(shape, rate) prior can be tuned.
# Mean of Gamma = shape/rate. Var = shape/rate^2.
# Let's try to center it around 0.25. Gamma(10, 40) -> Mean = 0.25.

variance_prior <- list(
    # Prior for residual precision of Growth_g_day
    tau_e_Growth_g_day = "dgamma(10, 40)"
)

fit_variance <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    priors = variance_prior,
    n.iter = 1000,
    quiet = TRUE
)

# Inspect the estimated sigma (derived from tau)
summary(fit_variance)
#>                                  Mean    SD Naive SE Time-series SE  2.5%   50%
#> alphaGrowth_g_day               9.776 0.204    0.013          0.013 9.399 9.791
#> beta_Growth_g_day_Temp_Centered 0.480 0.046    0.003          0.003 0.396 0.482
#> sigmaGrowth_g_day               1.961 0.135    0.009          0.009 1.702 1.958
#>                                  97.5%  Rhat n.eff
#> alphaGrowth_g_day               10.134 0.995   240
#> beta_Growth_g_day_Temp_Centered  0.577 0.999   282
#> sigmaGrowth_g_day                2.245 0.999   240
#> 
#> DIC:
#> Mean deviance:  417.4 
#> penalty 2.892 
#> Penalized deviance: 420.2

# Visualize the posterior for the residual standard deviation
plot_posterior(fit_variance, parameter = "sigma")
```

![](05_custom_priors_files/figure-html/variance_prior-1.png)

## Default Priors Reference

It is important to know what you are overriding. `because` aims to use
**weakly informative** defaults that provide minimal regularization
while allowing the data to dominate the posterior.

| Parameter Type        | Parameter Name | Default Prior        | Description                                                       |
|:----------------------|:---------------|:---------------------|:------------------------------------------------------------------|
| **Intercepts**        | `alpha_*`      | `dnorm(0, 1.0E-6)`   | Wide Normal (Precision 1e-6 = Variance 1,000,000)                 |
| **Coefficients**      | `beta_*`       | `dnorm(0, 1.0E-6)`   | Wide Normal                                                       |
| **Precision**         | `tau_*`        | `dgamma(1, 1)`       | Gamma(1,1). Weakly informative for precision/variance components. |
| **Phylo Signal**      | `lambda_*`     | `dunif(0, 1)`        | Uniform on \[0,1\]                                                |
| **Zero-Inflation**    | `psi_*`        | `dbeta(1, 1)`        | Beta(1,1) (Uniform on \[0,1\])                                    |
| **Ordinal Cutpoints** | `cutpoint_*`   | `dnorm(0, 1e-6)`     | First cutpoint fixed, others relative or specifically ordered     |
| **NegBinomial Size**  | `r_*`          | `dgamma(0.01, 0.01)` | Wide Gamma for dispersion parameter                               |

> **Precision vs Variance**: JAGS uses precision $\tau = 1/\sigma^{2}$.
> A prior of `dnorm(0, 1.0E-6)` means a normal distribution with mean 0
> and precision $0.000001$, which corresponds to a variance of
> $1,000,000$ (SD = 1000). This is very flat.

## Parameter Names Reference

To use custom priors, you need to know the exact internal name of the
parameter in the JAGS code.

Common patterns:

- `alpha_{Response}`: Intercept
- `beta_{Response}_{Predictor}`: Regression coefficient
- `tau_e_{Response}`: Residual precision
- `sigma_{Response}`: Residual standard deviation (derived, usually not
  set directly as prior, set `tau` instead)
- `lambda_{Response}`: Phylogenetic signal (0-1)
- `psi_{Response}`: Zero-inflation probability

If you are unsure, run a quick model with `n.adapt=0, n.iter=0` (just
compilation) and check the generated model code:

``` r
fit_check <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    n.iter = 0, quiet = TRUE
)

# Print the first few lines of the JAGS model
cat(substr(fit_check$model_code, 1, 1000))
#> model {
#>   # Dummy usage of ID to prevent warnings for unused data
#>   dummy_ID <- ID[1,1]
#>   # Structural equations
#>   for (i in 1:N) {
#> 
#>     muGrowth_g_day[i] <- alphaGrowth_g_day + beta_Growth_g_day_Temp_Centered*Temp_Centered[i]
#>   }
#>   # Multivariate normal likelihoods
#>   for (i in 1:N) {
#>     Growth_g_day[i] ~ dnorm(muGrowth_g_day[i], tau_e_Growth_g_day)
#>     log_lik_Growth_g_day[i] <- logdensity.norm(Growth_g_day[i], muGrowth_g_day[i], tau_e_Growth_g_day)
#>   }
#>   # Priors for structural parameters
#>   alphaGrowth_g_day ~ dnorm(0, 1.0E-6)
#>   tau_e_Growth_g_day ~ dgamma(1, 1)
#>   sigmaGrowth_g_day <- 1/sqrt(tau_e_Growth_g_day)
#>   beta_Growth_g_day_Temp_Centered ~ dnorm(0, 1.0E-6)
#> }
```

\`\`\`
