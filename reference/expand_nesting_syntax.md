# Expand Nesting Syntax in Random Effects

Converts lme4-style nesting (1\|A/B) to (1\|A) + (1\|A:B)

## Usage

``` r
expand_nesting_syntax(rhs)
```

## Arguments

- rhs:

  Character string, RHS of random effects formula

## Value

Expanded character string
