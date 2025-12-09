library(becauseR)
library(ape)

set.seed(123)
N <- 30
tree <- rtree(N)

# Simulate binomial data
X <- rnorm(N)
p_true <- plogis(-0.5 + 0.8 * X)
Y_binary <- rbinom(N, 1, p_true)

data_bin <- data.frame(
    X = X,
    Gregarious = Y_binary
)

cat("=== Testing binomial model (eigenvalue issue) ===\n")
cat("Running with 3 chains...\n")

fit_bin <- because(
    data = data_bin,
    tree = tree,
    equations = list(Gregarious ~ X),
    distribution = c(Gregarious = "binomial"),
    n.chains = 3,
    n.iter = 2000,
    n.burnin = 1000,
    quiet = FALSE
)

cat("\n✓ Success! No eigenvalue error\n")
print(summary(fit_bin))
