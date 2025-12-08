## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# data <- data.frame(
#     x = rnorm(100),
#     y = rnorm(100),
#     z = rnorm(100)
# )
# 
# # "quality" is latent (not in data)
# fit <- phybase_run(
#     equations = list(
#         quality ~ x, # x affects quality
#         y ~ quality, # quality affects y
#         z ~ quality # quality affects z
#     ),
#     data = data,
#     latent = c("quality"), # Explicit declaration (optional)
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # "body_condition" measured by three indicators
# fit <- phybase_run(
#     equations = list(
#         body_condition ~ mass, # Cause
#         fat_reserves ~ body_condition, # Indicator 1
#         muscle_mass ~ body_condition, # Indicator 2
#         immune_function ~ body_condition # Indicator 3
#     ),
#     data = data,
#     latent = c("body_condition"),
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Intelligence causes test scores
# equations <- list(
#     intelligence ~ education,
#     math_score ~ intelligence,
#     verbal_score ~ intelligence
# )

## ----eval=FALSE---------------------------------------------------------------
# # SES composed of income, education, occupation
# equations <- list(
#     SES ~ income + education + occupation,
#     health ~ SES
# )

## ----eval=FALSE---------------------------------------------------------------
# # Model with latent confounder
# equations <- list(
#     latent ~ x1,
#     x2 ~ latent,
#     x3 ~ latent,
#     y ~ x2 + x3
# )
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data, # "latent" not in data
#     dsep = TRUE
# )
# 
# # M-separation tests check:
# # - Are x2 and x3 correlated? (induced by latent)
# # - Conditional independencies given latent structure

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     latent = c("quality"),
#     standardize_latent = TRUE, # Scale latent to mean=0, SD=1
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # "environment" affects both traits
# # Not in data (unmeasured!)
# equations <- list(
#     environment ~ rainfall,
#     trait_a ~ environment,
#     trait_b ~ environment,
#     fitness ~ trait_a + trait_b
# )
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data, # environment is latent
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# equations <- list(
#     # Two latents: resources and predation
#     resources ~ rainfall,
#     predation ~ forest_cover,
# 
#     # Observable consequences
#     food_availability ~ resources,
#     shelter ~ resources,
#     mortality ~ predation,
#     vigilance ~ predation,
# 
#     # Fitness outcome
#     reproduction ~ food_availability + mortality
# )
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     latent = c("resources", "predation"),
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Nested latent structure
# equations <- list(
#     quality ~ environment, # quality caused by environment
#     environment ~ temperature, # Both latent!
#     survival ~ quality,
#     reproduction ~ quality
# )
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     latent = c("quality", "environment"),
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Good: Well-justified theoretical construct
# latent <- c("habitat_quality")
# 
# # Bad: Adding latents just to improve fit
# latent <- c("latent1", "latent2", "latent3") # No theory!

## ----eval=FALSE---------------------------------------------------------------
# # Good: 3+ indicators
# equations <- list(
#     quality ~ x,
#     indicator1 ~ quality,
#     indicator2 ~ quality,
#     indicator3 ~ quality
# )
# 
# # Weak: Only 1 indicator
# equations <- list(
#     quality ~ x,
#     indicator1 ~ quality # Underidentified!
# )

## ----eval=FALSE---------------------------------------------------------------
# # Use more iterations
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     latent = latent_vars,
#     n.chains = 4, # More chains
#     n.iter = 10000, # More iterations
#     n.burnin = 5000 # Longer burn-in
# )

## ----eval=FALSE---------------------------------------------------------------
# # Model WITH latent
# fit_latent <- phybase_run(equations_with_latent, data, latent = "quality")
# 
# # Model WITHOUT latent (direct paths instead)
# fit_direct <- phybase_run(equations_direct, data)
# 
# # Compare DIC/WAIC
# fit_latent$DIC
# fit_direct$DIC

