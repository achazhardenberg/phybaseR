library(because)
library(testthat)

message("Testing Custom Priors Implementation...")

# Mock data
set.seed(123)
N <- 50
df <- data.frame(
    x = rnorm(N),
    y = rnorm(N)
)

# Test 1: Standard Gaussian priors (Alpha/Beta/Tau)
equations <- list(y ~ x)
priors <- list(
    alphay = "dnorm(100, 1000)", # Strong intercept prior
    beta_y_x = "dnorm(5, 500)", # Strong slope prior
    tau_e_y = "dgamma(10, 10)" # Strong precision prior
)

message("Generating model with custom priors...")

# Run because() minimal run
fit <- because(
    equations = equations,
    data = df,
    priors = priors,
    n.chains = 1,
    n.iter = 50,
    n.adapt = 10,
    quiet = TRUE
)

# Extract model string
# fix$model is the JAGS object, fit$model_code is the string
if (!is.null(fit$model_code)) {
    model_str <- fit$model_code
} else {
    model_str <- fit$model
}

message(paste("Model string class:", class(model_str)))
message(paste("Model string length:", length(model_str)))

# Verification Checks
# We look for the exact strings in the generated model code
# distinct lines or collapsed block, any(grepl) handles both
found_alpha <- any(grepl("alphay ~ dnorm(100, 1000)", model_str, fixed = TRUE))
found_beta <- any(grepl("beta_y_x ~ dnorm(5, 500)", model_str, fixed = TRUE))
found_tau <- any(grepl("tau_e_y ~ dgamma(10, 10)", model_str, fixed = TRUE))

if (found_alpha) {
    message("SUCCESS: Custom Alpha prior found.")
} else {
    message("FAILURE: Custom Alpha prior NOT found.")
}
if (found_beta) {
    message("SUCCESS: Custom Beta prior found.")
} else {
    message("FAILURE: Custom Beta prior NOT found.")
}
if (found_tau) {
    message("SUCCESS: Custom Tau prior found.")
} else {
    message("FAILURE: Custom Tau prior NOT found.")
}

if (found_alpha && found_beta && found_tau) {
    message("\nALL CHECKS PASSED.")
} else {
    stop("Verification Failed.")
}
