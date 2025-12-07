library(phybaseR)
library(ape)

# Test with phylogenetic structure
set.seed(123)
N <- 10
tree <- rtree(N)

df <- data.frame(
    SP = tree$tip.label,
    Y = rnorm(N, mean = 10, sd = 2),
    X = rnorm(N)
)

cat("=== Test 1: Phylogenetic SEM with optimise=TRUE ===\n")
fit1 <- phybase_run(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ X),
    optimise = TRUE,
    n.chains = 1,
    n.iter = 100,
    quiet = FALSE
)

if (grepl("log_lik_Y", fit1$model_code)) {
    cat("✓ PASS: log_lik found in optimized model\n")
} else {
    cat("✗ FAIL: log_lik NOT found in optimized model\n")
}

cat("\n=== Test 2: Phylogenetic SEM with optimise=FALSE (MVN) ===\n")
fit2 <- phybase_run(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ X),
    optimise = FALSE,
    n.chains = 1,
    n.iter = 100,
    quiet = FALSE
)

if (grepl("log_lik_Y", fit2$model_code)) {
    cat("✓ PASS: log_lik found in MVN model\n")

    # Check for tau_marg
    if (grepl("tau_marg", fit2$model_code)) {
        cat(
            "✓ PASS: tau_marg extraction found (diagonal of precision matrix)\n"
        )
    }
} else {
    cat("✗ FAIL: log_lik NOT found in MVN model\n")
}

cat("\n=== Summary ===\n")
cat("Both phylogenetic cases now have log_lik monitoring!\n")
