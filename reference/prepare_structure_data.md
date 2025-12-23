# Prepare Structure Data

Modules implement this to process the structure object into data for
JAGS.

## Usage

``` r
prepare_structure_data(structure, data, ...)
```

## Arguments

- structure:

  The structure object.

- data:

  The model data.

- ...:

  Additional arguments.

## Value

A named list of data to be passed to JAGS (e.g., list(VCV = ...)).
