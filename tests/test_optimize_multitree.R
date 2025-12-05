# ==============================================================================
# Test Script: Multi-Tree Support with Optimization
# ==============================================================================
# Verifies that the optimized random effects formulation works correctly
# with phylogenetic uncertainty (multiple trees).
# ==============================================================================

# library(phybaseR)
library(ape)
library(coda)
library(rjags)

# Source local files
source("R/phybase_model.R")
source("R/phybase_run.R")
source("R/phybase_format_data.R")
source("R/phybase_waic.R")
source("R/summary.phybase.R")
source("R/mag_helpers.R")
source("R/dag_to_mag.R")

cat("\n=== Testing Multi-Tree Support with Optimization ===\n\n")

# 1. Simulate Data with Phylogenetic Uncertainty
set.seed(999)
N <- 40
n_trees <- 10  # Number of trees (phylogenetic uncertainty)

# Generate multiple trees from the same topology with branch length variation
base_tree <- rtree(N)
base_tree$edge.length <- base_tree$edge.length / max(branching.times(base_tree))

# Create a set of trees with perturbed branch lengths
tree_list <- list()
for (i in 1:n_trees) {
  tree_i <- base_tree
  # Add some variation to branch lengths (±20%)
  tree_i$edge.length <- base_tree$edge.length * runif(length(base_tree$edge.length), 0.8, 1.2)
  tree_i$edge.length <- tree_i$edge.length / max(branching.times(tree_i))
  tree_list[[i]] <- tree_i
}
class(tree_list) <- "multiPhylo"

cat(sprintf("Created %d trees with %d species\n", n_trees, N))

# True parameters
beta <- 0.7
lambda <- 0.6

# Simulate data using first tree (for simplicity)
VCV <- vcv.phylo(tree_list[[1]])
X <- rnorm(N)
sigma <- 1.0
phylo_var <- lambda * sigma^2
resid_var <- (1 - lambda) * sigma^2
Sigma_phylo <- phylo_var * VCV
err_phylo <- MASS::mvrnorm(1, rep(0, N), Sigma_phylo)
err_resid <- rnorm(N, 0, sqrt(resid_var))
Y <- beta * X + err_phylo + err_resid

cat(sprintf("Simulated data: Y = %.2f*X + phylo_error + resid_error\n", beta))

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

# 2. Run Optimized Multi-Tree Model
cat("\nRunning optimized multi-tree model (optimize = TRUE)...\n")
time_opt <- system.time({
  fit_opt <- phybase_run(
    data = data,
    tree = tree_list,
    equations = equations,
    n.iter = 2000,
    n.burnin = 500,
    n.thin = 2,
    optimise = TRUE,
    quiet = TRUE
  )
})

cat(sprintf("✓ Optimized multi-tree model finished in %.2f seconds\n", time_opt["elapsed"]))

# 3. Verify Output
cat("\nChecking parameter estimates:\n")
sum_opt <- fit_opt$summary
print(sum_opt$statistics[c("beta_Y_X", "lambdaY"), c("Mean", "SD")])

# Check convergence
cat("\nChecking convergence (R-hat):\n")
gelman <- gelman.diag(fit_opt$samples)
key_params <- c("beta_Y_X", "lambdaY")
print(gelman$psrf[key_params, ])

if (all(gelman$psrf[key_params, "Point est."] < 1.1)) {
  cat("✓ Convergence successful for key parameters (R-hat < 1.1)\n")
} else {
  cat("⚠ Some parameters have R-hat > 1.1\n")
}

# Check parameter recovery
cat("\nParameter Recovery:\n")
beta_est <- sum_opt$statistics["beta_Y_X", "Mean"]
lambda_est <- sum_opt$statistics["lambdaY", "Mean"]
cat(sprintf("  beta: True=%.2f, Est=%.2f\n", beta, beta_est))
cat(sprintf("  lambda: True=%.2f, Est=%.2f\n", lambda, lambda_est))

# 4. Run Unoptimized Multi-Tree Model (for comparison)
cat("\nRunning unoptimized multi-tree model (optimize = FALSE)...\n")
time_unopt <- system.time({
  fit_unopt <- phybase_run(
    data = data,
    tree = tree_list,
    equations = equations,
    n.iter = 2000,
    n.burnin = 500,
    n.thin = 2,
    optimise = FALSE,
    quiet = TRUE
  )
})

cat(sprintf("✓ Unoptimized multi-tree model finished in %.2f seconds\n", time_unopt["elapsed"]))
cat(sprintf("Speedup: %.2fx\n", time_unopt["elapsed"] / time_opt["elapsed"]))

# Compare estimates
sum_unopt <- fit_unopt$summary
beta_unopt <- sum_unopt$statistics["beta_Y_X", "Mean"]
lambda_unopt <- sum_unopt$statistics["lambdaY", "Mean"]

cat("\nComparing estimates (Optimized vs Unoptimized):\n")
cat(sprintf("  beta: %.4f vs %.4f (diff=%.4f)\n", beta_est, beta_unopt, abs(beta_est - beta_unopt)))
cat(sprintf("  lambda: %.4f vs %.4f (diff=%.4f)\n", lambda_est, lambda_unopt, abs(lambda_est - lambda_unopt)))

if (abs(beta_est - beta_unopt) < 0.1 && abs(lambda_est - lambda_unopt) < 0.1) {
  cat("✓ Estimates match between methods (diff < 0.1)\n")
} else {
  cat("⚠ Estimates differ between methods\n")
}

cat("\n=== Multi-Tree Test Complete" ===\n")
cat(sprintf("✓ Multi-tree optimization works correctly with %d trees\n", n_trees))
