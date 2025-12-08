## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Individuals nested within sites
# fit <- phybase_run(
#     equations = list(
#         weight ~ age + temperature
#     ),
#     data = data,
#     random = ~ (1 | site), # Random intercept for each site
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Crossed random effects
# fit <- phybase_run(
#     equations = list(
#         response ~ treatment
#     ),
#     data = data,
#     random = ~ (1 | individual) + (1 | year), # Both independent
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Plots nested within sites
# fit <- phybase_run(
#     equations = list(
#         biomass ~ rainfall
#     ),
#     data = data,
#     random = ~ (1 | site / plot), # Expands to (1|site) + (1|site:plot)
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Different random structures for different equations
# equations <- list(
#     weight ~ age + (1 | individual), # Weight varies by individual
#     age ~ cohort # Age is fixed
# )
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Age has a quadratic effect on weight
# fit <- phybase_run(
#     equations = list(
#         weight ~ age + I(age^2)
#     ),
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Cubic polynomial
# fit <- phybase_run(
#     equations = list(
#         y ~ x + I(x^2) + I(x^3)
#     ),
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         response ~ x + I(x^2) + z + I(z^2)
#     ),
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Effect of age depends on sex
# fit <- phybase_run(
#     equations = list(
#         weight ~ age * sex # Expands to: age + sex + age:sex
#     ),
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Interactive effects
# fit <- phybase_run(
#     equations = list(
#         growth ~ temperature * rainfall
#     ),
#     data = data,
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Sex as a predictor (automatically dummy-coded)
# fit <- phybase_run(
#     equations = list(
#         weight ~ age + sex
#     ),
#     data = data, # sex should be a factor
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Habitat type with 3 levels: forest, grassland, wetland
# fit <- phybase_run(
#     equations = list(
#         abundance ~ temperature + habitat
#     ),
#     data = data, # habitat is a factor with 3 levels
#     n.chains = 2,
#     n.iter = 2000
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         # Polynomial + interaction + random effects
#         weight ~ age + I(age^2) + age * sex + (1 | site),
# 
#         # Multiple predictors with random effects
#         survival ~ weight + temperature + (1 | year)
#     ),
#     data = data,
#     random = ~ (1 | individual), # Global random effect
#     n.chains = 3,
#     n.iter = 5000
# )

