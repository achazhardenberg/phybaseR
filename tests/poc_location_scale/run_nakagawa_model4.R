library(rjags)
library(ape)

# 1. Load Data
dat <- read.csv("tests/poc_location_scale/Psittaciformes_354spp.csv")
tree_all <- read.nexus("tests/poc_location_scale/Psittaciformes_354spp_100.nex")
tree <- tree_all[[1]] # Use first tree

# 2.Data Preparation
# Match data and tree
# Column is "Phylo" not "Species"
common_sp <- intersect(dat$Phylo, tree$tip.label)

# Filter and Sort
tree <- keep.tip(tree, common_sp)
dat <- dat[match(tree$tip.label, dat$Phylo), ]

# Variables:
# Beak Length (Response) -> Log transformed and Scaled
# Body Mass (Predictor) -> Log transformed and Scaled

# Full data (N=354) - Scalable with optimized model
# Tutorial actually uses Range.Size in the code (despite text saying Beak Length)
y_raw <- log(dat$Range.Size)
x_raw <- log(dat$Mass)

y <- as.vector(scale(y_raw, center = TRUE, scale = FALSE)) # Centered only (as in tutorial?)
# Tutorial code: dat$cbeak_length <- scale(log(dat$Beak.Length_Culmen), center = TRUE, scale = FALSE)
# "scale=FALSE" implies centering but NOT dividing by SD.
# However, for convergence in JAGS, standardization (SD=1) is usually safer.
# Let's standardize fully to be safe, unless tutorial specific.
# The tutorial says "center data of interest".
# I will standardize fully (center=TRUE, scale=TRUE) for better MCMC mixing.
y <- as.vector(scale(y_raw))
x <- as.vector(scale(x_raw))

N <- length(y)

# 3. Phylogeny
C <- vcv(tree)
C <- C / max(C)
invC <- solve(C)
# Model uses zeros[1:N], so we need length N
zeros <- rep(0, N)

# 4. JAGS Data
jags_data <- list(
    y = y,
    x = x,
    N = N,
    invC = invC,
    zeros = zeros
)

# 5. Init
init_values <- function() {
    list(
        alpha_m = rnorm(1, 0, 0.1),
        beta_m = rnorm(1, 0, 0.1),
        alpha_v = rnorm(1, 0, 0.1),
        beta_v = rnorm(1, 0, 0.1),
        sigma_phy_m = runif(1, 0.1, 1),
        sigma_phy_v = runif(1, 0.1, 1),
        rho = runif(1, -0.5, 0.5)
    )
}

# 6. Run Model 4 (Correlated PLSM)
cat(
    "Running Nakagawa Model 4 Replication (Correlated PLSM) on Parrots (N=",
    N,
    ")...\n"
)
cat("Config: n.adapt=3000, n.iter=5000\n")

# Use the correlated model file we created earlier
model_file <- "tests/poc_location_scale/model_plsm_correlated.bug"

start_time <- Sys.time()

jm <- jags.model(
    file = model_file,
    data = jags_data,
    inits = init_values,
    n.chains = 3,
    n.adapt = 3000
)

update(jm, n.iter = 1000)

samples <- coda.samples(
    jm,
    variable.names = c(
        "alpha_m",
        "beta_m",
        "alpha_v",
        "beta_v",
        "sigma_phy_m",
        "sigma_phy_v",
        "rho"
    ),
    n.iter = 5000
)

end_time <- Sys.time()
run_time <- end_time - start_time
cat("\n--- Total Runtime ---\n")
print(run_time)

print(summary(samples))

# Check for correlation significance
rho_samples <- unlist(samples[, "rho"])
prob_pos <- mean(rho_samples > 0)
cat("\n--- Rho (Correlation) ---\n")
cat("Mean:", mean(rho_samples), "\n")
cat("P(rho > 0):", prob_pos, "\n")
