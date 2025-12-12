library(rjags)
# Load simulation function
source("tests/poc_location_scale/simulate_data.R")

# 1. Generate Data
cat("Simulating data...\n")
sim <- simulate_plsm_data(N = 100, seed = 42)
data_list <- sim$data

# 2. Prepare for JAGS
# Invert C covariance matrix for precision
invC <- solve(data_list$C)
zeros <- rep(0, data_list$N)

jags_data <- list(
    y = data_list$y,
    x = data_list$x,
    N = data_list$N,
    invC = invC,
    zeros = zeros
)

# 3. Initialization
init_values <- function() {
    list(
        alpha_m = rnorm(1, 0, 0.1),
        beta_m = rnorm(1, 0, 0.1),
        alpha_v = rnorm(1, 0, 0.1),
        beta_v = rnorm(1, 0, 0.1),
        sigma_phy_m = runif(1, 0.1, 1),
        sigma_phy_v = runif(1, 0.1, 1)
    )
}

# 4. Run Model
cat("Running JAGS model...\n")
model_file <- "tests/poc_location_scale/model_plsm.bug"

jm <- jags.model(
    file = model_file,
    data = jags_data,
    inits = init_values,
    n.chains = 3,
    n.adapt = 1000
)

update(jm, n.iter = 1000)

samples <- coda.samples(
    jm,
    variable.names = c(
        "alpha_m",
        "beta_m",
        "alpha_v",
        "beta_v",
        "lambda_m_est",
        "lambda_v_est"
    ),
    n.iter = 5000
)

# 5. Results
cat("\n--- Simulation Parameters ---\n")
print(unlist(sim$params))

cat("\n--- Estimated Parameters (Summary) ---\n")
print(summary(samples)$statistics[, c("Mean", "SD", "Time-series SE")])

cat("\n--- Quantiles ---\n")
print(summary(samples)$quantiles)

# Check coverage
means <- summary(samples)$statistics[, "Mean"]
# Simple check logic could go here
