# Compare Models Using LOO-CV

Wrapper for
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
to compare multiple Because models.

## Usage

``` r
because_loo_compare(...)
```

## Arguments

- ...:

  Two or more `loo` objects from
  [`because_loo()`](https://achazhardenberg.github.io/because/reference/because_loo.md).

## Value

A comparison table ranking models by expected out-of-sample predictive
accuracy.

## Examples

``` r
if (FALSE) { # \dontrun{
  loo1 <- because_loo(fit1)
  loo2 <- because_loo(fit2)
  because_loo_compare(loo1, loo2)
} # }
```
