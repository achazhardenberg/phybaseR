library(rjags)
source("tests/poc_location_scale/simulate_data_correlated.R")

cat("Simulating Correlated Data (N=50 for speed)...\n")
# Reduced N because 100 -> 200x200 matrix inversion in JAGS might be slow
sim <- simulate_plsm_data_correlated(N = 50, seed = 999, rho = 0.7)
data_list <- sim$data

# Check if rho matters
# cor(sim$latent$u_m, sim$latent$u_v)

invC <- solve(data_list$C)
zeros_2N <- rep(0, 2 * data_list$N)

jags_data <- list(
    y = data_list$y,
    x = data_list$x,
    N = data_list$N,
    invC = invC,
    zeros_2N = zeros_2N
)

init_values <- function() {
    list(
        alpha_m = rnorm(1),
        beta_m = rnorm(1),
        alpha_v = rnorm(1),
        beta_v = rnorm(1),
        sigma_phy_m = runif(1, 0.5, 1),
        sigma_phy_v = runif(1, 0.5, 1),
        rho = runif(1, -0.5, 0.5)
    )
}

cat("Running Correlated PLSM JAGS Model...\n")
jm <- jags.model(
    file = "tests/poc_location_scale/model_plsm_correlated.bug",
    data = jags_data,
    inits = init_values,
    n.chains = 3,
    n.adapt = 1000
)

update(jm, n.iter = 1000)

samples <- coda.samples(
    jm,
    variable.names = c("beta_m", "beta_v", "sigma_phy_m", "sigma_phy_v", "rho"),
    n.iter = 5000
)

print(summary(samples))

cat("\nTrue Rho:", sim$params$rho, "\n")
cat("Estimated Rho Mean:", summary(samples)$statistics["rho", "Mean"], "\n")
