# Extract Imputed Values for Missing Data

Extracts the posterior distributions of missing values that were imputed
by the model.

## Usage

``` r
extract_imputed(object, id_col = NULL)
```

## Arguments

- object:

  A `because` model object.

- id_col:

  Optional character vector of IDs (e.g., species names) corresponding
  to the rows of the data. If the data in the model object already has
  names (e.g., from a phylogenetic model), these are used automatically.
  If provided, this vector must have the same length as the number of
  observations in the model (N).

## Value

A data frame containing:

- `Variable`: Name of the variable with missing data.

- `ID`: The identifier (e.g., Species) for the observation (if
  available).

- `RowIndex`: The original row index in the data.

- `Mean`: Posterior mean of the imputed value.

- `SD`: Posterior standard deviation.

- `Q2.5`: 2.5% quantile (lower credible interval).

- `Q50`: Median.

- `Q97.5`: 97.5% quantile (upper credible interval).

## Details

This function identifies which values in the original data were missing
(`NA`) and looks up the corresponding imputed nodes in the posterior
samples.

Note: The imputed values are only available if they were monitored
during the run. You must use `monitor = "all"` (or include the variable
names in `monitor`) when running
[`because()`](https://because-pkg.github.io/because/reference/because.md)
to ensure these values are saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'fit' is a because model run with missing data and monitor="all"
imputed_values <- extract_imputed(fit)
head(imputed_values)
} # }
```
