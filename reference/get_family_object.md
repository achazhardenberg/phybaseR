# Get Family Object for S3 Dispatch

Converts a family name string into a family class object for S3
dispatch. This is the gatekeeper that ensures module packages are
installed before allowing specialized families to be used.

## Usage

``` r
get_family_object(family_name)
```

## Arguments

- family_name:

  The name of the family (e.g., "occupancy", "gaussian").

## Value

An object of class `because_family_<name>` and `because_family`.
