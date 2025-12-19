# Plot Posterior Distributions for Model Comparison

Visualizes posterior density distributions for selected parameters,
allowing easy comparison between multiple models (e.g., Default vs.
Custom priors) or inspection of a single model.

## Usage

``` r
plot_posterior(
  models,
  parameter,
  col = NULL,
  lwd = 2,
  legend_pos = "topleft",
  ...
)
```

## Arguments

- models:

  A `because` object or a named list of `because` objects. If a list is
  provided, the names are used for the legend. Example:
  `list("Default" = fit1, "Custom" = fit2)`.

- parameter:

  Character string. A regular expression or exact name of the
  parameter(s) to plot. Example: `"beta"`, `"^alpha"`, `"beta_Y_X"`.
  Plots all matching parameters.

- col:

  Optional vector of colors. Defaults to a standard palette (black,
  blue, red, etc.).

- lwd:

  Line width (default = 2).

- legend_pos:

  Legend position (default = "topleft"). Set to `NULL` to suppress.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible NULL. Produces a plot.
