# Calculate WAIC with Standard Errors for a Because Model

Calculates the Widely Applicable Information Criterion (WAIC) with
standard errors for a fitted Because model using pointwise
log-likelihoods.

## Usage

``` r
because_waic(model)
```

## Arguments

- model:

  A fitted model object of class `"because"` returned by
  [`because`](https://achazhardenberg.github.io/because/reference/because.md)
  with `WAIC = TRUE`.

## Value

A data frame with columns `Estimate` and `SE` containing:

- elpd_waic:

  Expected log pointwise predictive density (higher is better)

- p_waic:

  Effective number of parameters

- waic:

  The WAIC value (lower is better for model comparison)

The returned object also has a `pointwise` attribute containing
individual observation contributions for model comparison.

## Details

This function is automatically called by
[`because`](https://achazhardenberg.github.io/because/reference/because.md)
when `WAIC = TRUE`. It can also be called manually. If the model was not
originally fitted with `WAIC = TRUE` (so pointwise log-likelihoods are
missing), this function will automatically refit the model (using a
short MCMC run) to compute them.

## WAIC Definition

The Widely Applicable Information Criterion (WAIC) is calculated as:
\$\$WAIC = -2 \times (lppd - p\_{waic})\$\$ where:

- \\lppd = \sum\_{i=1}^N \log(\frac{1}{S} \sum\_{s=1}^S
  \exp(log\\lik\_{is}))\\ is the log pointwise predictive density

- \\p\_{waic} = \sum\_{i=1}^N \text{var}(log\\lik\_{is})\\ is the
  effective number of parameters

**WAIC Algorithm**:

1.  **lpd** (log pointwise predictive density): For each observation
    \\i\\, compute \\\log(\text{mean}(\exp(\text{log\\lik}\_i)))\\
    across MCMC samples

2.  **p_waic**: For each observation \\i\\, compute
    \\\text{var}(\text{log\\lik}\_i)\\ across MCMC samples

3.  **elpd_waic**: \\\text{lpd}\_i - \text{p\\waic}\_i\\ for each
    observation

4.  **waic**: \\-2 \times \sum \text{elpd\\waic}\_i\\

## Standard Errors

Standard errors for WAIC are calculated using the pointwise
contributions: \$\$SE(WAIC) = \sqrt{N \times \text{var}(waic_i)}\$\$
where \\waic_i = -2 \times (lppd_i - p\_{waic,i})\\.

## References

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413-1432.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Fit model with WAIC monitoring
  fit <- because(data, tree, equations, WAIC = TRUE)

  # View WAIC with standard errors
  fit$WAIC
  #             Estimate   SE
  # elpd_waic   -617.3   12.4
  # p_waic        12.3    3.1
  # waic        1234.5   24.8

  # Compare two models
  fit1$WAIC
  fit2$WAIC
  # Model with lower WAIC is preferred
  # Difference is significant if |WAIC1 - WAIC2| > 2 * sqrt(SE1^2 + SE2^2)
} # }
```
