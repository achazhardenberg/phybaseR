## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Your model
# equations <- list(
#     z ~ x,
#     y ~ z
# )
# 
# # Implied independence: x ⊥ y | z
# # "x and y are independent given z"
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     dsep = TRUE # Test d-separation!
# )

## ----eval=FALSE---------------------------------------------------------------
# # Hypothesized model
# equations <- list(
#     range_size ~ body_size,
#     population ~ range_size
# )
# 
# fit <- phybase_run(equations, data, dsep = TRUE)
# 
# # D-sep result:
# # population _||_ body_size | range_size
# # P(~0) = 0.02  ← FAIL!
# 
# # Interpretation: body_size affects population even controlling for range_size
# # Fix: Add direct path
# equations_fixed <- list(
#     range_size ~ body_size,
#     population ~ range_size + body_size # Added!
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         body_mass ~ age + sex
#     ),
#     data = data,
#     distribution = list(body_mass = "gaussian") # Default
# )

## ----eval=FALSE---------------------------------------------------------------
# # Binary: survival (0 = dead, 1 = alive)
# fit <- phybase_run(
#     equations = list(
#         survival ~ age + condition
#     ),
#     data = data,
#     distribution = list(survival = "binomial")
# )
# 
# # Proportions: hatching success
# data$N_trials <- 10 # Clutch size
# data$N_success <- ... # Number hatched
# 
# fit <- phybase_run(
#     equations = list(
#         hatching_success ~ temperature
#     ),
#     data = data,
#     distribution = list(hatching_success = "binomial")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Number of offspring
# fit <- phybase_run(
#     equations = list(
#         offspring_count ~ age + body_size
#     ),
#     data = data,
#     distribution = list(offspring_count = "poisson")
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         parasite_load ~ host_size
#     ),
#     data = data,
#     distribution = list(parasite_load = "negbinomial")
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         lifespan ~ body_mass
#     ),
#     data = data,
#     distribution = list(lifespan = "gamma")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Habitat choice: Forest, Grassland, Wetland (no inherent order)
# data$habitat <- factor(data$habitat,
#     levels = c("Forest", "Grassland", "Wetland")
# )
# 
# fit <- phybase_run(
#     equations = list(
#         habitat ~ temperature + rainfall
#     ),
#     data = data,
#     distribution = list(habitat = "multinomial")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Body condition: Poor < Fair < Good < Excellent
# data$condition <- factor(data$condition,
#     levels = c("Poor", "Fair", "Good", "Excellent"),
#     ordered = TRUE
# )
# 
# fit <- phybase_run(
#     equations = list(
#         condition ~ age + food_availability
#     ),
#     data = data,
#     distribution = list(condition = "ordinal")
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         survival ~ age, # Binomial
#         offspring ~ survival, # Poisson
#         body_mass ~ age # Gaussian
#     ),
#     data = data,
#     distribution = list(
#         survival = "binomial",
#         offspring = "poisson",
#         body_mass = "gaussian"
#     )
# )

## ----eval=FALSE---------------------------------------------------------------
# summary(fit)
# 
# # Check Rhat column:
# #           Estimate  Rhat
# # beta_x      0.45   1.002  ✓ Good
# # beta_y      0.32   1.15   ✗ Poor - needs more iterations

## ----eval=FALSE---------------------------------------------------------------
# # Low n.eff:
# #           n.eff
# # beta_x      45   ✗ Too low
# 
# # Solution: Thin more or run longer
# fit <- phybase_run(
#     equations, data,
#     n.iter = 20000,
#     n.thin = 5 # Save every 5th iteration
# )

## ----eval=FALSE---------------------------------------------------------------
# # Extract samples
# samples <- fit$samples
# 
# # Plot traces
# plot(samples[, "beta_x"], type = "l")
# 
# # Good: "fuzzy caterpillar" - well-mixed, stationary
# # Bad: trends, getting stuck, different chains separated

## ----eval=FALSE---------------------------------------------------------------
# fit1 <- phybase_run(equations1, data, DIC = TRUE)
# fit2 <- phybase_run(equations2, data, DIC = TRUE)
# 
# fit1$DIC # 245.3
# fit2$DIC # 238.1  ← Better
# 
# # Difference > 5-10: Strong evidence for better model

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations, data,
#     WAIC = TRUE,
#     n.chains = 2 # Requires 2+ chains
# )
# 
# fit$WAIC$WAIC # Overall WAIC
# fit$WAIC$lppd # Log pointwise predictive density

## ----eval=FALSE---------------------------------------------------------------
# # Standardize predictors beforehand
# data$age_std <- scale(data$age)
# data$mass_std <- scale(data$mass)
# 
# fit <- phybase_run(
#     equations = list(
#         survival ~ age_std + mass_std
#     ),
#     data = data
# )
# 
# # Now beta estimates are comparable!

