# Format Data for Because Analysis

Converts data from long format (one row per observation) to the list
format required by
[`because`](https://achazhardenberg.github.io/because/reference/because.md).

## Usage

``` r
because_format_data(data, species_col = "SP", tree)
```

## Arguments

- data:

  A data.frame in long format with one row per observation.

- species_col:

  Name of the column containing species identifiers (default: "SP").

- tree:

  A phylogenetic tree (class `phylo`). Required to determine species
  order.

## Value

A named list where each element is either:

- A numeric vector (if all species have exactly 1 observation)

- A numeric matrix with species in rows and replicates in columns

Species are ordered to match `tree$tip.label`.

## Details

This function handles:

- Different numbers of replicates per species (creates rectangular
  matrix with NA padding)

- Missing values (NA)

- Automatic alignment with phylogenetic tree tip labels

When species have different numbers of replicates, the function creates
a matrix with dimensions (number of species) x (maximum number of
replicates). Species with fewer replicates are padded with NA values.

Species in the tree but not in the data will have all NA values. Species
in the data but not in the tree will be excluded with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example data in long format
data_long <- data.frame(
  SP = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3"),
  BM = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8),
  NL = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)
)

tree <- ape::read.tree(text = "(sp1:1,sp2:1,sp3:1);")
data_list <- because_format_data(data_long, species_col = "SP", tree = tree)

# Use with because
fit <- because(data = data_list, tree = tree, equations = list(NL ~ BM))
} # }
```
