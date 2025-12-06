library(phybaseR)
library(ape)

# Reproduce MAG bug and print model
set.seed(123)
N <- 10 # Small for easier debugging
tree <- rtree(N)

df <- data.frame(
    SP = tree$tip.label,
    Y1 = rnorm(N, mean = 10, sd = 2),
    Y2 = rnorm(N, mean = 8, sd = 1.5)
)

cat("Generating MAG model...\n")

# Just generate the model without running
model_output <- phybaseR:::phybase_model(
    eq_list = list(
        list(response = "Y1", predictors = "L1"),
        list(response = "Y2", predictors = "L1")
    ),
    all_vars = c("Y1", "Y2", "L1"),
    response_vars = c("Y1", "Y2"),
    latent = "L1",
    latent_method = "correlations",
    structure = tree,
    structure_names = "phylo",
    independent = FALSE,
    optimise = TRUE,
    multi.tree = FALSE
)

cat("\n=== JAGS Model Code ===\n")
cat(model_output$model_string)
cat("\n\n=== Checking for TAU_phylo definitions ===\n")
if (grepl("TAU_phylo_Y1", model_output$model_string)) {
    cat("✓ TAU_phylo_Y1 defined\n")
} else {
    cat("✗ TAU_phylo_Y1 NOT defined\n")
}

if (grepl("TAU_phylo_Y2", model_output$model_string)) {
    cat("✓ TAU_phylo_Y2 defined\n")
} else {
    cat("✗ TAU_phylo_Y2 NOT defined\n")
}
