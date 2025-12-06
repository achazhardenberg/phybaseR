# ==============================================================================
# Parameter Comparison: JAGS vs NIMBLE
# ==============================================================================
# Run this after completing the experimental_nimble_benchmark.R script
# to compare parameter estimates
# ==============================================================================

# Assuming samples_jags and nimble_samples are already in your environment
# If not, re-run the benchmark script first

library(coda)

cat("\n=== PARAMETER COMPARISON: JAGS vs NIMBLE ===\n\n")

# Extract summaries
jags_summary <- summary(samples_jags)
nimble_summary <- summary(nimble_samples)

# Parameters to compare
params <- c(
    "betaBM",
    "betaBM2",
    "betaRS",
    "betaNL",
    "lambdaLS",
    "lambdaNL",
    "lambdaDD"
)

# Create comparison table
cat(sprintf(
    "%-12s %12s %12s %12s %10s\n",
    "Parameter",
    "JAGS Mean",
    "NIMBLE Mean",
    "Difference",
    "% Diff"
))
cat(rep("-", 70), "\n", sep = "")

results <- data.frame(
    Parameter = params,
    JAGS_Mean = numeric(length(params)),
    NIMBLE_Mean = numeric(length(params)),
    Diff = numeric(length(params)),
    Pct_Diff = numeric(length(params))
)

for (i in seq_along(params)) {
    p <- params[i]

    # JAGS mean
    jags_mean <- jags_summary$statistics[p, "Mean"]

    # NIMBLE mean (combine chains if multiple)
    if (is.list(nimble_samples)) {
        # Multiple chains - combine
        all_samples <- do.call(rbind, nimble_samples)
        nimble_mean <- mean(all_samples[, p])
    } else {
        # Single chain
        nimble_mean <- mean(nimble_samples[, p])
    }

    diff <- abs(jags_mean - nimble_mean)
    pct_diff <- (diff / abs(jags_mean)) * 100

    results$JAGS_Mean[i] <- jags_mean
    results$NIMBLE_Mean[i] <- nimble_mean
    results$Diff[i] <- diff
    results$Pct_Diff[i] <- pct_diff

    cat(sprintf(
        "%-12s %12.4f %12.4f %12.4f %9.2f%%\n",
        p,
        jags_mean,
        nimble_mean,
        diff,
        pct_diff
    ))
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 70), "\n\n", sep = "")

cat(sprintf("Maximum absolute difference: %.4f\n", max(results$Diff)))
cat(sprintf("Mean absolute difference:    %.4f\n", mean(results$Diff)))
cat(sprintf("Maximum percent difference:  %.2f%%\n", max(results$Pct_Diff)))

if (max(results$Diff) < 0.05) {
    cat("\n✓ Excellent agreement (max diff < 0.05)\n")
    cat("  NIMBLE and JAGS are producing statistically equivalent results.\n")
} else if (max(results$Diff) < 0.10) {
    cat("\n✓ Good agreement (max diff < 0.10)\n")
    cat("  Differences are within expected MCMC variation.\n")
} else {
    cat("\n⚠ Notable differences detected (max diff > 0.10)\n")
    cat("  Consider running longer chains or checking convergence.\n")
}

# Optional: print full comparison
cat("\n\nFull comparison table:\n")
print(results, row.names = FALSE)


cat("\n=== CONVERGENCE DIAGNOSTICS ===\n\n")

cat("JAGS R-hat:\n")
print(gelman.diag(samples_jags))

cat("\nNIMBLE R-hat:\n")
nimble_coda <- as.mcmc.list(lapply(nimble_samples, as.mcmc))
print(gelman.diag(nimble_coda))
