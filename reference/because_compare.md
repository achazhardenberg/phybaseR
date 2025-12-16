# Compare Because Models

A unified function to either (1) compare previously fitted models, or
(2) run multiple model specifications in parallel and then compare them.

## Usage

``` r
because_compare(
  ...,
  model_specs = NULL,
  data = NULL,
  tree = NULL,
  n.cores = 1,
  cl = NULL,
  sort = TRUE
)
```

## Arguments

- ...:

  For comparing fitted models: individual fitted model objects of class
  `"because"`. For running models: additional arguments passed to
  [`because`](https://because-pkg.github.io/because/reference/because.md)
  (e.g., `n.iter`).

- model_specs:

  A named list of model specifications to run (Mode 2). Each element
  should be a list containing arguments for `because`. Alternatively,
  this argument can accept the first fitted model object (Mode 1).

- data:

  The dataset (required for Mode 2). Alternatively, the second fitted
  model object (Mode 1).

- tree:

  The phylogenetic tree (optional for Mode 2). Alternatively, the third
  fitted model object (Mode 1).

- n.cores:

  Number of cores for parallel execution (Mode 2). Default is 1.

- cl:

  Optional cluster object (Mode 2).

- sort:

  Logical. If `TRUE` (default), sort comparison table by WAIC.

## Value

If comparing fitted models: A class `"because_comparison"` object (data
frame) with WAIC rankings.

If running models: A list containing:

- results:

  List of fitted model objects.

- comparison:

  The comparison data frame.

## Details

**Mode 1: Compare Fitted Models** Call `because_compare(fit1, fit2)` or
`because_compare(models = list(fit1, fit2))`. Extracts WAIC (with SE)
from each model and ranks them.

**Mode 2: Run and Compare** Call
`because_compare(model_specs = list(m1=..., m2=...), data=data, tree=tree)`.
This runs the models in parallel and returns the comparison.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Mode 1: Compare existing fits
  because_compare(fit1, fit2)

  # Mode 2: Run and compare
  specs <- list(m1 = list(equations = list(Y ~ X)), m2 = list(equations = list(Y ~ X + Z)))
  res <- because_compare(specs, data = df, tree = tr, n.cores = 2)
  print(res$comparison)
} # }
```
