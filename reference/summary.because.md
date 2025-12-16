# Summary for Because Model

Summarizes the output of a Because model run.

## Usage

``` r
# S3 method for class 'because'
summary(object, ...)
```

## Arguments

- object:

  A fitted model object of class `"because"`.

- ...:

  Additional arguments passed to
  [`summary.mcmc`](https://rdrr.io/pkg/coda/man/summary.mcmc.html).

## Value

A summary object containing statistics for the monitored parameters. If
`dsep = TRUE` was used in `because`, the summary focuses on the
conditional independence tests.
