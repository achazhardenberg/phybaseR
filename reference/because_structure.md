# Create a Custom Covariance Structure for use with because

This helper function allows users to easily create custom covariance
structures by providing just a function that computes the precision (or
covariance) matrix. The function handles all the S3 method registration
automatically.

## Usage

``` r
because_structure(name, precision_fn, description = NULL)
```

## Arguments

- name:

  Character string; the name of your structure (e.g., "spatial_knn").
  This will be used to create the S3 class and methods.

- precision_fn:

  A function that takes your structure data and returns a precision
  matrix (N x N). The first return value should be the precision matrix.
  Additional list elements can be returned and will be passed to JAGS as
  data.

- description:

  Optional description of the structure for documentation.

## Value

A constructor function that creates structure objects of your custom
class.

## Details

This function creates:

- A constructor function (returned) for creating structure objects

- S3 methods for `jags_structure_definition` (generates JAGS code)

- S3 methods for `prepare_structure_data` (prepares data for JAGS)

The precision function should return either:

- A single matrix (the precision matrix)

- A list with `Prec` (the precision matrix) and optionally other data

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a spatial distance-decay structure
spatial_decay <- because_structure(
  name = "spatial_decay",
  precision_fn = function(coords, decay_rate = 0.1) {
    dist_mat <- as.matrix(dist(coords))
    W <- exp(-decay_rate * dist_mat)
    diag(W) <- 0
    D <- diag(rowSums(W))
    D - 0.99 * W  # Precision matrix
  }
)

# Use it
my_struct <- spatial_decay(coords = my_coords, decay_rate = 0.2)
fit <- because(Y ~ X, data = data, structure = my_struct)
} # }
```
