# Identify Reusable D-Separation Tests

Scans a list of previously fitted 'because' models to find d-separation
tests that have already been run and can be reused. Reusability requires
strict data identity and identical test formula + required support
equations.

## Usage

``` r
find_reusable_tests(
  current_tests,
  current_equations,
  reuse_models,
  current_data,
  family = NULL,
  quiet = FALSE
)
```

## Arguments

- current_tests:

  List of formula objects representing the d-separation tests to run.

- current_equations:

  List of model formulas for the current model (support equations).

- reuse_models:

  List of 'because' model objects to scan for reusable results.

- current_data:

  The data object (list or environment) being used for the current run.

- family:

  Named vector of families (needed for dependency resolution).

- quiet:

  Logical; suppress messages (default = FALSE).

## Value

A list with two components:

- `found`: A list of d-separation result objects found in previous
  models.

- `missing_indices`: A numeric vector of indices in `current_tests` that
  need to be run.
