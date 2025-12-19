# Extract Deterministic Terms from Formulas

Scans a list of formulas for terms that require deterministic nodes in
JAGS (interactions, I() calls, etc.)

## Usage

``` r
extract_deterministic_terms(equations)
```

## Arguments

- equations:

  List of formulas

## Value

List of deterministic term definitions
