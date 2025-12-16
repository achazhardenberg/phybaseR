# Rhinograd Simulated Repeated Measures Data

A simulated dataset containing repeated measurements of life-history
traits for 100 Rhinograd species. This dataset is in "long" format, with
10 observations per species, suitable for testing measurement error
models using `variability = "reps"`.

## Usage

``` r
RhinoMulti.dat
```

## Format

A data frame with 1000 rows (10 replicates x 100 species) and 6 columns:

- SP:

  Character. Species names matching the tip labels in `rhino.tree`. Each
  species appears 10 times.

- BM:

  Numeric. Body mass (repeated measurements)

- LS:

  Numeric. Litter size (repeated measurements)

- NL:

  Numeric. Nose length (repeated measurements)

- DD:

  Numeric. Dispersal distance (repeated measurements)

- RS:

  Numeric. Range size (repeated measurements)

## Source

Simulated based on Gonzalez-Voyer & von Hardenberg (2014) data.

## Details

This dataset was simulated by adding white noise to the original
`rhino.dat` species values. It represents a scenario where each species
has 10 independent observations for each trait. This allows using the
`variability = "reps"` option in
[`because()`](https://achazhardenberg.github.io/because/reference/because.md)
to explicitly model measurement error.

The underlying causal model is the same as for `rhino.dat`.

## See also

[`rhino.dat`](https://achazhardenberg.github.io/because/reference/rhino.dat.md),
[`RhinoMulti_se.dat`](https://achazhardenberg.github.io/because/reference/RhinoMulti_se.dat.md)
