# Sanitize Term Name for JAGS

Converts complex R terms into valid JAGS variable names

## Usage

``` r
sanitize_term_name(term)
```

## Arguments

- term:

  Character string (e.g., "A:B", "I(A^2)")

## Value

Sanitized string (e.g., "A_x_B", "A_pow2")
