# Rhinograd phylogenetic tree

A phylogenetic tree (phylo object) for 100 simulated Rhinograd species.
This tree is used to demonstrate phylogenetic Bayesian structural
equation modeling with the becauseR package.

## Usage

``` r
rhino.tree
```

## Format

A phylo object (from the ape package) with 100 tips and 99 internal
nodes. The tree has been scaled so that the root age is 1.0.

## Source

Simulated data for Gonzalez-Voyer and von Hardenberg (2014)

## Details

Rhinograds are hypothetical mammals used as examples in phylogenetic
comparative methods. This simulated tree represents the evolutionary
relationships among 100 Rhinograd species.

## References

Gonzalez-Voyer, A., & von Hardenberg, A. (2014). An introduction to
phylogenetic path analysis. In: Modern phylogenetic comparative methods
and their application in evolutionary biology (pp. 201–229). Springer.

## Examples

``` r
data(rhino.tree)
plot(rhino.tree)

```
