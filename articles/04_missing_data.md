# Missing Data Imputation in because

The `because` package handles missing values (NAs) automatically using
Bayesian imputation. Unlike frequentist methods that often require
deleting rows with missing values (list-wise deletion), `because` treats
missing values as unknown parameters to be estimated by the model.

This vignette demonstrates how `because` handles missing values in
both: 1. **Response variables**: Imputed based on predictors and
phylogenetic correlations. 2. **Predictor variables**: Automatically
treated as response variables to enable imputation.

## Setup

First, let’s load the package and the standard `rhino` dataset.

``` r
library(because)
data(rhino.dat)
data(rhino.tree)
```

## Simulating Missing Data

We will artificially introduce missing values (`NA`) into the dataset to
demonstrate how `because` handles them. We’ll add missingness to: \*
`LS` (Litter Size): A response variable in our model. \* `BM` (Body
Mass): A predictor variable.

``` r
# Create a copy of the data
rhino_na <- rhino.dat

# Introduce 10% missing values in Litter Size (Response)
set.seed(123)
na_indices_ls <- sample(seq_len(nrow(rhino_na)), size = 10)
rhino_na$LS[na_indices_ls] <- NA

# Introduce 5% missing values in Body Mass (Predictor)
na_indices_bm <- sample(seq_len(nrow(rhino_na)), size = 5)
rhino_na$BM[na_indices_bm] <- NA

# Verify missingness
sum(is.na(rhino_na$LS))
sum(is.na(rhino_na$BM))
```

## The Model

We define the same structural equation model (sem8) as in previous
vignettes: \* Body Mass (`BM`) affects Litter Size (`LS`) and Nose
Length (`NL`). \* Range Size (`RS`) affects Nose Length (`NL`). \* Nose
Length (`NL`) affects Dispersal Distance (`DD`).

``` r
sem8.eq <- list(
  LS ~ BM,
  NL ~ BM + RS,
  DD ~ NL
)
```

## Running the Model

When we run
[`because()`](https://because-pkg.github.io/because/reference/because.md)
with missing data, it automatically detects the `NA` values and adjusts
the JAGS model structure.

``` r
sem8_na.bc <- because(
  equations = sem8.eq,
  data = rhino_na,
  structure = rhino.tree
)

summary(sem8_na.bc)

#                  Mean    SD Naive SE Time-series SE   2.5%    50% 97.5%
# alphaBM         1.283 0.545    0.010          0.029  0.170  1.294 2.348
# alphaDD         0.574 0.280    0.005          0.007  0.030  0.566 1.133
# alphaLS         0.346 0.407    0.007          0.015 -0.496  0.362 1.113
# alphaNL        -0.014 0.330    0.006          0.013 -0.659 -0.007 0.655
# beta_DD_NL      0.538 0.076    0.001          0.001  0.394  0.539 0.687
# beta_LS_BM      0.478 0.092    0.002          0.002  0.296  0.479 0.656
# beta_NL_BM      0.457 0.070    0.001          0.002  0.316  0.458 0.595
# beta_NL_RS      0.601 0.067    0.001          0.001  0.468  0.602 0.732
# lambdaBM        0.809 0.088    0.002          0.003  0.592  0.827 0.931
# lambdaDD        0.455 0.130    0.002          0.003  0.213  0.453 0.705
# lambdaLS        0.656 0.116    0.002          0.003  0.409  0.668 0.848
# lambdaNL        0.722 0.089    0.002          0.002  0.525  0.731 0.868
# sigma_BM_phylo  1.750 0.260    0.005          0.008  1.244  1.745 2.255
# sigma_BM_res    0.812 0.135    0.002          0.004  0.581  0.799 1.099
# sigma_DD_phylo  0.793 0.171    0.003          0.004  0.497  0.781 1.158
# sigma_DD_res    0.858 0.088    0.002          0.002  0.690  0.854 1.033
# sigma_LS_phylo  1.188 0.208    0.004          0.005  0.811  1.175 1.614
# sigma_LS_res    0.837 0.113    0.002          0.002  0.629  0.829 1.067
# sigma_NL_phylo  1.037 0.150    0.003          0.004  0.772  1.029 1.366
# sigma_NL_res    0.628 0.079    0.001          0.002  0.488  0.622 0.798
#                 Rhat n.eff
# alphaBM        1.004   356
# alphaDD        1.000  1437
# alphaLS        1.002   734
# alphaNL        1.000   623
# beta_DD_NL     1.000  2910
# beta_LS_BM     1.000  2366
# beta_NL_BM     1.000  2198
# beta_NL_RS     1.002  2819
# lambdaBM       1.000  1152
# lambdaDD       1.000  2211
# lambdaLS       1.000  1630
# lambdaNL       1.003  1421
# sigma_BM_phylo 1.000  1183
# sigma_BM_res   1.000  1493
# sigma_DD_phylo 1.000  2275
# sigma_DD_res   1.000  2853
# sigma_LS_phylo 1.000  1742
# sigma_LS_res   1.000  2156
# sigma_NL_phylo 1.002  1307
# sigma_NL_res   1.002  2183
# 
# DIC:
# Mean deviance:  884.1 
# penalty 189.7 
# Penalized deviance: 1074 
```

### What happened automatically?

You will notice two important messages if you run the above model:

1.  **“Detected missing data in predictor-only variables: BM”** JAGS
    cannot impute missing values in a “fixed” covariate data node. To
    solve this, `because` automatically adds an intercept-only equation
    for the missing predictor (e.g., `BM ~ 1`). This converts `BM` into
    a stochastic node (a response variable) so that JAGS can sample its
    missing values during MCMC.

2.  **“Using Latent Variable (GLMM) approach…”** When variables have
    missing data, `because` handles them carefully to preserve
    phylogenetic signal. Instead of simple imputation, it uses an
    element-wise likelihood that integrates the phylogenetic error
    structure, ensuring that the relationships are estimated using the
    full phylogenetic information available.

## Inspecting Imputed Values

Since the missing values are now parameters in the model, we can
actually inspect their estimated values (posterior distributions). To do
this, you must run
[`because()`](https://because-pkg.github.io/because/reference/because.md)
with the argument `monitor = "all"`. This ensures that the imputed
stochastic nodes are saved in the output.

You can use the helper function
[`extract_imputed()`](https://because-pkg.github.io/because/reference/extract_imputed.md)
to easily retrieve these values:

``` r
# Run strictly monitoring 'all' parameters to get imputed values
sem8_na_monitored.bc <- because(
  equations = sem8.eq,
  data = rhino_na,
  structure = rhino.tree,
  monitor = "all"
)

# Extract posterior summaries for imputed values
imputed_values <- extract_imputed(sem8_na_monitored.bc)

# View the first few imputed values
head(imputed_values)

#       Variable  ID RowIndex         Mean        SD       Q2.5        Q50
# Mean        LS s14       14 -0.004710686 0.9342152 -1.8430269 0.02246534
# Mean1       LS s25       25  1.678043988 1.0829263 -0.4576528 1.73566859
# Mean2       LS s31       31  0.561576970 1.1507743 -1.6971686 0.57087689
# Mean3       LS s42       42  1.162445522 1.2150872 -1.3452741 1.13780737
# Mean4       LS s43       43  2.781955491 1.0852266  0.7474752 2.79039031
# Mean5       LS s50       50  2.762965171 1.1253402  0.5564907 2.71909773
#          Q97.5
# Mean  1.740653
# Mean1 3.646009
# Mean2 2.760267
# Mean3 3.589646
# Mean4 4.957594
# Mean5 4.982451
```

This returns a tidy data frame with the `Mean`, `SD`, and Credible
Intervals (`Q2.5`, `Q97.5`) for each missing data point, labeled by the
Species ID.

The variable `LS` had missing values. In the JAGS model, the missing
values for `LS` are treated as latent variables. The variable `BM`
(missing predictor) was converted to a response variable (`BM ~ 1`) and
also imputed.
