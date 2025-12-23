# Extension Packages

`because` is designed as a modular framework. While the core package
handles standard Bayesian SEMs, specialized functionality is provided
through extension packages.

## Available packages

### Phylogenetic Structural Equation Models (`because.phybase`)

The [because.phybase](https://because-pkg.github.io/because.phybase/)
package adds support for phylogenetic covariance structures, enabling
the fitting of Phylogenetic Bayesian Structural Equation Models
(PhyBaSE; von Hardenberg & Gonzalez-Voyer, 2025).

*Features*

- **Phylogenetic Covariance**: Incorporates phylogenetic relatedness
  into `because` models.
- **Phylogenetic Uncertainty**: Supports `multiPhylo` objects to account
  for uncertainty in tree topology or branch lengths.
- **Measurement error**: Accounting for measurement error in traits both
  providing repeated measures as well as specifying known measurement
  error variances.
- **Phylogenetic missing data imputation**: informs imputation of
  missing data using phylogenetic relationships.
- **Seamless Integration**: Designed to work transparently with
  `because` via S3 method dispatch.

## Planned packages

### Occupancy (`because.occupancy`)

*Status: In Development* Adds support for hierarchical occupancy models
(single-season, multi-season) directly within the causal framework.

### Spatial (`because.spatial`)

*Status: Planned* Will add support for spatially structured residuals
(CAR/SAR models).
