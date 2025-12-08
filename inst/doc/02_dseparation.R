## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # Hypothesized model: X → Z → Y
# equations <- list(
#     z ~ x,
#     y ~ z
# )
# 
# # DAG implies: X ⊥ Y | Z
# # "X and Y are independent given Z"

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     dsep = TRUE # Test d-separation!
# )

## ----eval=FALSE---------------------------------------------------------------
# # Based on theory, not data!
# equations <- list(
#     range_size ~ body_size,
#     population ~ range_size + body_size
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     dsep = TRUE, # Critical!
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# summary(fit)

## ----eval=FALSE---------------------------------------------------------------
# # Updated model
# equations_v2 <- list(
#     range_size ~ body_size,
#     population ~ range_size + body_size,
#     dispersal ~ body_size # Added!
# )
# 
# fit_v2 <- phybase_run(equations_v2, data, dsep = TRUE)

