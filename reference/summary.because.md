# Summary for Because Model

Summarizes the output of a Because model run.

## Usage

``` r
# S3 method for class 'because'
summary(
  object,
  show_internal = FALSE,
  show_nodes = FALSE,
  show_random = FALSE,
  ...
)
```

## Arguments

- show_internal:

  Logical. If `TRUE`, shows internal parameters created for
  deterministically defined nodes (e.g., `beta_..._det_...`). Defaults
  to `FALSE`.

- show_nodes:

  Logical. If `TRUE`, shows latent node values (e.g., `Age[1]`).
  Defaults to `FALSE` to prevent clutter when `monitor="all"`.

- show_random:

  Logical. If `TRUE`, shows random effect estimates (e.g., `u_...`).
  Defaults to `FALSE`.

- ...:

  Additional arguments passed to
  [`summary.mcmc`](https://rdrr.io/pkg/coda/man/summary.mcmc.html).

## Value

A summary object containing statistics for the monitored parameters. If
`dsep = TRUE` was used in `because`, the summary focuses on the
conditional independence tests.
