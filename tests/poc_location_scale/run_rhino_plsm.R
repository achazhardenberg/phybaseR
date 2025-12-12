library(rjags)
library(ape)
library(because) # Load package for data

# 1. Load Data
data(rhino.dat)
data(rhino.tree)

# 2. Select Variables for Test: LS (Lifespan) ~ BM (Body Mass)
# Hypothesis: Does Body Mass affect the mean AND variance of Lifespan?
y_raw <- rhino.dat$LS
x_raw <- rhino.dat$BM
species <- rhino.dat$SP

# Match with tree (just to be safe, though usually pre-cleaned)
common_sp <- intersect(species, rhino.tree$tip.label)
# Sort data by tree tip label order
rhino.tree <- keep.tip(rhino.tree, common_sp)
rhino.dat_sorted <- rhino.dat[match(rhino.tree$tip.label, rhino.dat$SP), ]

y <- scale(rhino.dat_sorted$LS) # Standardize
x <- scale(rhino.dat_sorted$BM) # Standardize
N <- length(y)

# 3. Prepare Phylogeny
C <- vcv(rhino.tree)
C <- C / max(C) # Standardize tree height to 1
invC <- solve(C)
zeros <- rep(0, N)

jags_data <- list(
    y = as.vector(y),
    x = as.vector(x),
    N = N,
    invC = invC,
    zeros = zeros
)

# 4. Initialization
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

# 5. Run Model
cat("Running PLSM on Rhinogradentia Data (LS ~ BM)...\n")
cat("Hypothesis: Body Mass affects the mean and variance of Lifespan.\n\n")

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

# 6. Results
cat("\n--- Estimated Parameters ---\n")
print(summary(samples)$statistics[, c("Mean", "SD", "Time-series SE")])

cat("\n--- Quantiles ---\n")
print(summary(samples)$quantiles)

# Does BM affect variance? Check beta_v
beta_v_samples <- unlist(samples[, "beta_v"])
prob_pos <- mean(beta_v_samples > 0)
cat("\n--- Hypothesis Test (beta_v != 0) ---\n")
cat("P(beta_v > 0) =", prob_pos, "\n")
cat("P(beta_v < 0) =", 1 - prob_pos, "\n")

if (prob_pos > 0.975 || prob_pos < 0.025) {
    cat("Result: SIGNIFICANT effect of Body Mass on Lifespan variance.\n")
} else {
    cat("Result: No significant effect of Body Mass on Lifespan variance.\n")
}
