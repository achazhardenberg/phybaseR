library(because)
set.seed(42)

# Recreate Example 2 Data
Rainfall <- runif(50, 0, 10)
# Biologically, rainfall increases biomass.
# But we simulate weak/noisy data where slope could be estimated as negative without prior.
Biomass <- 0.1 * Rainfall + rnorm(50, sd = 2)
df_bio <- data.frame(Rainfall, Biomass)

# Enforce positive slope using truncation
positive_prior <- list(
    beta_Biomass_Rainfall = "dnorm(0, 1) T(0, )"
)

message("Fitting truncated model...")
fit <- because(
    equations = list(Biomass ~ Rainfall),
    data = df_bio,
    priors = positive_prior,
    n.iter = 1000,
    quiet = TRUE
)

# Extract samples for the slope
samples <- as.matrix(fit$samples)[, "beta_Biomass_Rainfall"]

message("Summary of Slope Samples:")
print(summary(samples))

min_val <- min(samples)
message(paste("Minimum sample value:", min_val))

if (min_val >= 0) {
    message(
        "SUCCESS: All samples are non-negative. The plot overlap is a density smoothing artifact."
    )
} else {
    message("FAILURE: Some samples are negative! Truncation failed.")
}

# Check density default behavior
d <- density(samples)
message(paste("Density plot x-min:", min(d$x)))
