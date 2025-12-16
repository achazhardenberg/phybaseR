# Calculate Pagel's Lambda from Variance Components

This function calculates Pagel's lambda for each response variable in a
fitted `because` model by extracting the posterior samples of the
phylogenetic and residual variance components.

## Usage

``` r
because_lambda(model, prob = 0.95)
```

## Arguments

- model:

  A fitted model object of class `"because"`.

- prob:

  A numeric value specifying the probability mass for the credibility
  intervals (default 0.95).

## Value

A data frame containing the summary statistics for the derived lambda
parameter(s):

- Mean:

  Posterior mean

- SD:

  Posterior standard deviation

- Median:

  Posterior median

- LowerCI:

  Lower bound of the credibility interval

- UpperCI:

  Upper bound of the credibility interval

## Details

In the optimized formulation of `because`, the phylogenetic signal
(\\\lambda\\) is not always explicitly monitored for complex models with
multiple structures. However, it can be derived post-hoc from the
estimated standard deviations of the phylogenetic (\\\sigma\_{phylo}\\)
and residual (\\\sigma\_{res}\\) components: \$\$\lambda =
\frac{\sigma\_{phylo}^2}{\sigma\_{phylo}^2 + \sigma\_{res}^2}\$\$

This function performs this calculation on the full posterior chains to
provide valid Bayesian estimates and credibility intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- because(...)
  because_lambda(fit)
} # }
```
