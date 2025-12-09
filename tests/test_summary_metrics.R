library(becauseR)
library(ape)

# Simulate data
set.seed(123)
N <- 50
tree <- rtree(N)
data <- data.frame(
    X = rnorm(N),
    Y = rnorm(N)
)
data$Y <- 0.5 * data$X + rnorm(N)


# Generate model code directly to inspect it
model_out <- becauseR:::because_model(
    equations = list(Y ~ X),
    multi.tree = FALSE,
    distribution = NULL,
    variability = NULL,
    vars_with_na = NULL,
    induced_correlations = NULL
)

cat("--- Generated JAGS Model ---\n")
cat(model_out$model)
cat("\n--------------------------\n")

# Run model with 2 chains to enable Rhat calculation
fit <- because(
    data = data,
    tree = tree,
    equations = list(Y ~ X),
    n.chains = 2,
    n.iter = 1000,
    n.burnin = 500,
    quiet = TRUE
)

# Capture summary output
cat("--- Summary Output ---\n")
summ <- summary(fit)

# Check if Rhat and n.eff are present in the output object (if invisible) or printed
# Since summary prints, we can also inspect the returned object if I modified it to return the combined table
# In my modification:
# For standard summary, I return `invisible(combined)` where combined has Rhat and n.eff columns.

if ("Rhat" %in% colnames(summ)) {
    cat("\n[PASS] Rhat column found in summary object.\n")
    print(head(summ[, c("Rhat", "n.eff")]))
} else {
    cat("\n[FAIL] Rhat column NOT found in summary object.\n")
}

if ("n.eff" %in% colnames(summ)) {
    cat("\n[PASS] n.eff column found in summary object.\n")
} else {
    cat("\n[FAIL] n.eff column NOT found in summary object.\n")
}
