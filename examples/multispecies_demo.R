# Example: Causal Multi-Species Occupancy Model
# ---------------------------------------------
# This script demonstrates the 'Auto-Stacking' and 'd-separation'
# features of the 'because' package for multispecies data.

# Examples from the project root
devtools::load_all(".")

# 1. Simulate Multispecies Data
# -----------------------------
set.seed(123)
n_sites <- 40
n_reps <- 3 # Multiple visits per site

# Environmental Covariates (Independent)
habitat <- rnorm(n_sites) # Driven by biology
wind <- rnorm(n_sites) # Driven by meteorology

# Species 1: Occupancy driven by Habitat, Detection by Wind
psi1 <- plogis(0.5 + 2.0 * habitat)
y1_true <- rbinom(n_sites, 1, psi1)
y1_obs <- matrix(NA, n_sites, n_reps)
for (i in 1:n_sites) {
    y1_obs[i, ] <- rbinom(n_reps, 1, y1_true[i] * plogis(0 - 1.5 * wind[i]))
}

# Species 2: Occupancy also driven by Habitat (Independent of Sp1)
psi2 <- plogis(-1.0 + 1.5 * habitat)
y2_true <- rbinom(n_sites, 1, psi2)
y2_obs <- matrix(NA, n_sites, n_reps)
for (i in 1:n_sites) {
    y2_obs[i, ] <- rbinom(n_reps, 1, y2_true[i] * plogis(0.5))
}

# Data format: List of matrices for Y, vectors for covariates
data_list <- list(
    Y = list(BlueTit = y1_obs, GreatTit = y2_obs),
    Habitat = habitat,
    Wind = wind
)

# 2. Fit Causal Occupancy Model
# -----------------------------
# 'because' will:
# 1. Automatically stack 'BlueTit' and 'GreatTit' into a hierarchical dataset.
# 2. Add a random effect (1|SpeciesID) to account for species differences.
# 3. Use MAG to test conditional independencies (d-separation).
# 4. Inject an 'Observation Layer' (Y_obs) for proper structural testing.

cat("\n--- Running Multi-Species Causal Model ---\n")
results <- because(
    equations = list(
        Y ~ Habitat, # Causal driver of state
        p_Y ~ Wind # Causal driver of detection
    ),
    data = data_list,
    distribution = c(Y = "occupancy"),
    dsep = TRUE, # Enable d-separation tests
    n.iter = 2500, # Sampling iterations
    quiet = FALSE
)

# 3. Review Results
# -----------------
cat("\n--- Model Summary ---\n")
print(summary(results))

# Look closely at the "d-separation Tests" section in the output.
# It should confirm that Habitat and Wind are independent, verifying
# that your environmental predictors aren't confounded.
