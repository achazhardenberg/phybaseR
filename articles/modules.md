# Extension Packages

`because` is designed as a modular framework. While the core package
handles standard Bayesian SEMs, specialized functionality is provided
through extension packages.

## Available packages

### Phylogenetic Structural Equation Models (`because.phybase`)

The **because.phybase** package adds support for phylogenetic covariance
structures, enabling “Phylogenetic Bayesian Structural Equation
Modeling” (PhyBaSE).

- **Features**: Single tree (`phylo`) and multiple trees (`multiPhylo`)
  for uncertainty in phylogeny.
- **Repository**:
  [github.com/because-pkg/because.phybase](https://github.com/because-pkg/because.phybase)
- **Documentation**: [Reference
  Site](https://because-pkg.github.io/because.phybase/) (Coming Soon)

To install:

``` r
remotes::install_github("because-pkg/because.phybase")
```

## Planned packages

### Occupancy (`because.occupancy`)

*Status: In Development* Adds support for hierarchical occupancy models
(single-season, multi-season) directly within the causal framework.

### Spatial (`because.spatial`)

*Status: Planned* Will add support for spatially structured residuals
(CAR/SAR models).
