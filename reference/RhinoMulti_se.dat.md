# Rhinograd Means and Standard Errors

A simulated dataset containing the mean and standard error of varying
life-history traits for 100 Rhinograd species. This dataset is in "wide"
format (one row per species), suitable for measurement error models
using `variability = "se"`.

## Usage

``` r
RhinoMulti_se.dat
```

## Format

A data frame with 100 rows and 11 columns:

- SP:

  Character. Species names matching the tip labels in `rhino.tree`.

- BM:

  Numeric. Mean Body mass

- BM_se:

  Numeric. Standard error of Body mass (SD / sqrt(N))

- LS:

  Numeric. Mean Litter size

- LS_se:

  Numeric. Standard error of Litter size

- NL:

  Numeric. Mean Nose length

- NL_se:

  Numeric. Standard error of Nose length

- DD:

  Numeric. Mean Dispersal distance

- DD_se:

  Numeric. Standard error of Dispersal distance

- RS:

  Numeric. Mean Range size

- RS_se:

  Numeric. Standard error of Range size

## Source

Calculated from `RhinoMulti.dat`.

## Details

This dataset is derived from `RhinoMulti.dat` by calculating the mean
and standard error for each trait per species (N=10 replicates). It is
designed for use with the `variability = "se"` option in
[`because()`](https://because-pkg.github.io/because/reference/because.md),
where precision is modeled as `1/SE^2`.

## See also

[`rhino.dat`](https://because-pkg.github.io/because/reference/rhino.dat.md),
[`RhinoMulti.dat`](https://because-pkg.github.io/because/reference/RhinoMulti.dat.md)
