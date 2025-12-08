## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Survival: 0 = died, 1 = survived
# data <- data.frame(
#     individual_id = 1:100,
#     survival = rbinom(100, 1, 0.7), # Binary: 0 or 1
#     age = rnorm(100, 5, 2),
#     body_mass = rnorm(100, 50, 10)
# )
# 
# fit <- phybase_run(
#     equations = list(
#         survival ~ age + body_mass
#     ),
#     data = data,
#     distribution = list(survival = "binomial"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Hatching success: eggs hatched / eggs laid
# data <- data.frame(
#     nest_id = 1:50,
#     eggs_laid = round(runif(50, 5, 12)),
#     temperature = rnorm(50, 25, 3)
# )
# data$eggs_hatched <- rbinom(50, data$eggs_laid, plogis(-2 + 0.3 * data$temperature))
# data$prop_hatched <- data$eggs_hatched / data$eggs_laid
# 
# # phybaseR needs: N (trials) and k (successes)
# data$N_trials <- data$eggs_laid
# data$N_success <- data$eggs_hatched
# 
# fit <- phybase_run(
#     equations = list(
#         prop_hatched ~ temperature
#     ),
#     data = data,
#     distribution = list(prop_hatched = "binomial"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Number of offspring
# data <- data.frame(
#     female_id = 1:80,
#     offspring_count = rpois(80, 3),
#     age = rnorm(80, 4, 1.5),
#     territory_quality = rnorm(80, 0, 1)
# )
# 
# fit <- phybase_run(
#     equations = list(
#         offspring_count ~ age + territory_quality
#     ),
#     data = data,
#     distribution = list(offspring_count = "poisson"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Parasite load (often overdispersed)
# data <- data.frame(
#     host_id = 1:100,
#     parasite_count = rnbinom(100, mu = 10, size = 2), # Overdispersed!
#     host_size = rnorm(100, 100, 20),
#     immune_function = rnorm(100, 0, 1)
# )
# 
# fit <- phybase_run(
#     equations = list(
#         parasite_count ~ host_size + immune_function
#     ),
#     data = data,
#     distribution = list(parasite_count = "negbinomial"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Lifespan (always positive, often right-skewed)
# data <- data.frame(
#     individual_id = 1:120,
#     lifespan = rgamma(120, shape = 2, rate = 0.5),
#     body_mass = rnorm(120, 50, 10),
#     environmental_stress = rnorm(120, 0, 1)
# )
# 
# fit <- phybase_run(
#     equations = list(
#         lifespan ~ body_mass + environmental_stress
#     ),
#     data = data,
#     distribution = list(lifespan = "gamma"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Habitat selection: Forest, Grassland, Wetland
# data <- data.frame(
#     territory_id = 1:100,
#     habitat = sample(c("Forest", "Grassland", "Wetland"), 100, replace = TRUE),
#     temperature = rnorm(100, 20, 5),
#     precipitation = rnorm(100, 800, 200),
#     elevation = rnorm(100, 500, 100)
# )
# 
# # IMPORTANT: Must be a factor!
# data$habitat <- factor(data$habitat,
#     levels = c("Forest", "Grassland", "Wetland")
# )
# 
# fit <- phybase_run(
#     equations = list(
#         habitat ~ temperature + precipitation + elevation
#     ),
#     data = data,
#     distribution = list(habitat = "multinomial"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Body condition: Poor < Fair < Good < Excellent
# data <- data.frame(
#     individual_id = 1:150,
#     condition = sample(c("Poor", "Fair", "Good", "Excellent"), 150, replace = TRUE),
#     food_availability = rnorm(150, 0, 1),
#     age = rnorm(150, 5, 2),
#     parasite_load = rpois(150, 3)
# )
# 
# # IMPORTANT: Must be ordered factor!
# data$condition <- factor(data$condition,
#     levels = c("Poor", "Fair", "Good", "Excellent"),
#     ordered = TRUE
# )
# 
# fit <- phybase_run(
#     equations = list(
#         condition ~ food_availability + age + parasite_load
#     ),
#     data = data,
#     distribution = list(condition = "ordinal"),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = list(
#         survival ~ age, # Binomial
#         offspring ~ survival, # Poisson
#         body_mass ~ age, # Gaussian
#         lifespan ~ body_mass, # Gamma
#         habitat ~ temperature # Multinomial
#     ),
#     data = data,
#     distribution = list(
#         survival = "binomial",
#         offspring = "poisson",
#         body_mass = "gaussian", # Default, can omit
#         lifespan = "gamma",
#         habitat = "multinomial"
#     ),
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Poisson/negative binomial: Check overdispersion
# # Compare variance to mean in residuals
# 
# # Binomial: Plot observed vs predicted probabilities
# 
# # Ordinal: Check proportional odds assumption
# # (effects should be similar across cutpoints)

## ----eval=FALSE---------------------------------------------------------------
# # Try both Poisson and negative binomial
# fit_pois <- phybase_run(..., distribution = list(count = "poisson"), DIC = TRUE)
# fit_nb <- phybase_run(..., distribution = list(count = "negbinomial"), DIC = TRUE)
# 
# fit_pois$DIC # 543.2
# fit_nb$DIC # 512.8  ← Better (lower)
# 
# # Negative binomial fits better (ΔIC > 10)

