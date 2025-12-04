# Debug script to print JAGS model for negbinomial

source("R/phybase_model.R")

# Simple test
model_out <- phybase_model(
    equations = list(Y ~ X),
    distribution = c(Y = "negbinomial"),
    optimize = TRUE,
    multi.tree = FALSE
)

cat("JAGS Model:\n")
cat(model_out$model)
cat("\n")
