library(phybaseR)
library(ape)

# Simulate data
set.seed(123)
N <- 50
tree <- rtree(N)
data <- data.frame(
    X = rnorm(N),
    Y = rnorm(N)
)
data$Y <- 0.5 * data$X + rnorm(N)

cat("=== TEST 1: 3 Chains (Should show Rhat and n.eff) ===\n")
fit3 <- phybase_run(
    data = data,
    tree = tree,
    equations = list(Y ~ X),
    n.chains = 3,
    n.iter = 1000,
    n.burnin = 500,
    quiet = TRUE
)
print(summary(fit3))

cat("\n=== TEST 2: 1 Chain (Should show n.eff, NO Rhat) ===\n")
fit1 <- phybase_run(
    data = data,
    tree = tree,
    equations = list(Y ~ X),
    n.chains = 1,
    n.iter = 1000,
    n.burnin = 500,
    quiet = TRUE
)
print(summary(fit1))
