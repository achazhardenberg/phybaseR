# Ensure Log-Likelihoods are Available

Internal helper to check if a model object has pointwise
log-likelihoods. If not, it automatically refits the model (with reduced
settings) to calculate them.

## Usage

``` r
ensure_log_lik(model, quiet = FALSE)
```

## Arguments

- model:

  A fitted `because` model object.

- quiet:

  Logical; if `TRUE`, suppresses messages.

## Value

A matrix of pointwise log-likelihoods (samples x observations).
