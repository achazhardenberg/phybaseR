library(because)
library(rjags)

set.seed(42)

# 1. Simulate Data (Null Case: Prey and Predator Independent, but clustered by Site)
N_sites <- 50
N_reps <- 5 # Visits per site
sigma_site <- 1.5 # Significant site variation

sites <- factor(rep(1:N_sites, each = N_reps))
u_site <- rnorm(N_sites, 0, sigma_site)
p_pred <- plogis(-1 + u_site[sites]) # Predator occupancy
p_prey <- plogis(-1 + u_site[sites]) # Prey occupancy (Same site preference, but independent process)

# Observations
y_pred <- rbinom(length(sites), 1, p_pred)
y_prey <- rbinom(length(sites), 1, p_prey)

data <- data.frame(
    Site = sites,
    Pred = y_pred,
    Prey = y_prey
)

# 2. Run 'because' with Random Effect (GLMM)
# We expect beta_Pred to be ~ 0
message("Fitting GLMM...")
mod <- because(
    equations = list(Prey ~ Pred),
    random = ~ 1 | Site,
    family = c(Prey = "binomial"), # Named vector for family
    data = data,
    n.iter = 1000,
    quiet = FALSE
)

# 3. Check Results
print(summary(mod))

# 4. Check if random effects were actually included in JAGS code
cat("\n--- Checking Random Effects ---\n")
cat("Random argument passed:", deparse(~ 1 | Site), "\n")
cat("Family:", "binomial", "\n")
cat("\n--- JAGS Model ---\n")
cat(mod$model_string)
