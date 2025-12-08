## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Individual-level data (finest grain)
# individual_data <- data.frame(
#     individual_id = 1:100,
#     site = rep(c("A", "B"), each = 50),
#     year = rep(c(2020, 2021), 50),
#     body_mass = rnorm(100, 50, 10),
#     age = rnorm(100, 5, 2)
# )
# 
# # Site-year level data (only 4 unique combinations!)
# site_year_data <- data.frame(
#     site = c("A", "A", "B", "B"),
#     year = c(2020, 2021, 2020, 2021),
#     temperature = c(15, 16, 17, 18),
#     rainfall = c(800, 850, 900, 950)
# )
# 
# # Specify hierarchical structure
# fit <- phybase_run(
#     data = list(
#         individual = individual_data,
#         site_year = site_year_data
#     ),
#     levels = list(
#         individual = c("body_mass", "age"),
#         site_year = c("temperature", "rainfall")
#     ),
#     hierarchy = "site_year > individual",
#     link_vars = c("site", "year"),
#     equations = list(
#         body_mass ~ age + temperature,
#         age ~ rainfall
#     ),
#     dsep = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Plot-level: soil, vegetation
# # Individual-level: plant traits
# 
# fit <- phybase_run(
#     data = list(
#         plot = plot_data,
#         individual = plant_data
#     ),
#     levels = list(
#         plot = c("soil_pH", "canopy_cover"),
#         individual = c("height", "leaf_area")
#     ),
#     hierarchy = "plot > individual",
#     link_vars = "plot_id"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Year-level: climate
# # Visit-level: measurements
# 
# fit <- phybase_run(
#     data = list(
#         year = yearly_climate,
#         visit = individual_visits
#     ),
#     levels = list(
#         year = c("temperature", "rainfall"),
#         visit = c("weight", "condition")
#     ),
#     hierarchy = "year > visit",
#     link_vars = "year"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Region > Site > Quadrat
# 
# fit <- phybase_run(
#     data = list(
#         region = region_data,
#         site = site_data,
#         quadrat = quadrat_data
#     ),
#     levels = list(
#         region = c("climate_zone"),
#         site = c("elevation", "soil_type"),
#         quadrat = c("species_count", "biomass")
#     ),
#     hierarchy = "region > site > quadrat",
#     link_vars = c("region_id", "site_id")
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     data = list(
#         site = site_data,
#         individual = individual_data
#     ),
#     levels = list(
#         site = c("habitat_quality"),
#         individual = c("body_mass", "age")
#     ),
#     hierarchy = "site > individual",
#     link_vars = "site_id",
#     random = ~ (1 | site_id), # Random site effects
#     equations = list(
#         body_mass ~ age + habitat_quality
#     )
# )

