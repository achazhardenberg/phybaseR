# Rhinograd life-history data

A dataset containing simulated life-history trait data for 100 Rhinograd
species. This data is used to demonstrate phylogenetic Bayesian
structural equation modeling with the becauseR package.

## Usage

``` r
rhino.dat
```

## Format

A data frame with 100 rows (one per species) and 6 columns:

- SP:

  Character. Species names matching the tip labels in rhino.tree

- BM:

  Numeric. Body mass (standardised)

- LS:

  Numeric. Litter size (standardised)

- NL:

  Numeric. Nose length (standardised)

- DD:

  Numeric. Dispersal distance (standardised)

- RS:

  Numeric. Range size (standardised)

## Source

Simulated data for Gonzalez-Voyer & von Hardenberg (2014)

## Details

This simulated dataset represents life-history traits for hypothetical
Rhinograd species. The data can be used to test causal hypotheses about
life-history evolution using phylogenetic structural equation models.

Example causal model (Model 8 from Gonzalez-Voyer and von Hardenberg
(2014)):

- Body mass (BM) affects litter size (LS)

- Body mass (BM) and range size (RS) affect nose length (NL)

- Nose length (NL) affects dispersal distance (DD)

## References

Gonzalez-Voyer, A., & von Hardenberg, A. (2014). An introduction to
phylogenetic path analysis. In: Modern phylogenetic comparative methods
and their application in evolutionary biology (pp. 201–229). Springer.

## Examples

``` r
data(rhino.dat)
head(rhino.dat)
#>   SP         BM        NL          LS         DD         RS
#> 1 s1 -0.7669046 -2.017574 -0.09476713 -0.7922619 -3.3926965
#> 2 s2 -1.0097044 -1.742466 -2.14633494  1.7643208 -1.0467189
#> 3 s3 -1.2252812 -2.468553 -0.63206843 -1.8061550 -2.8596606
#> 4 s4  1.3305352  0.340559 -0.13426489  0.9488741 -0.2786894
#> 5 s5  2.4923437  1.224005  0.60823620  1.6241513  0.2208142
#> 6 s6  1.4251486  2.969227  3.08553399  2.2524221  2.4268047
summary(rhino.dat)
#>       SP                  BM                NL                 LS          
#>  Length:100         Min.   :-2.3373   Min.   :-2.66773   Min.   :-2.14633  
#>  Class :character   1st Qu.: 0.4526   1st Qu.: 0.01542   1st Qu.:-0.09736  
#>  Mode  :character   Median : 1.5854   Median : 1.10201   Median : 1.20535  
#>                     Mean   : 1.6824   Mean   : 1.05258   Mean   : 1.23445  
#>                     3rd Qu.: 2.7206   3rd Qu.: 2.14295   3rd Qu.: 2.26928  
#>                     Max.   : 5.4368   Max.   : 5.02897   Max.   : 5.69822  
#>        DD                RS           
#>  Min.   :-2.9747   Min.   :-3.392697  
#>  1st Qu.: 0.1442   1st Qu.:-0.001407  
#>  Median : 1.0129   Median : 0.879031  
#>  Mean   : 1.0841   Mean   : 0.965484  
#>  3rd Qu.: 2.1425   3rd Qu.: 2.057613  
#>  Max.   : 4.4773   Max.   : 4.942733  
```
