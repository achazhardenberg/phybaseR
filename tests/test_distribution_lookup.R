# Test script to reproduce distribution lookup error

source("R/because_format_data.R")
source("R/because.R")
source("R/because_model.R")

library(ape)
library(MASS)

set.seed(123)
N <- 30
tree <- rtree(N)

# Simulate data
VCV <- ape::vcv.phylo(tree)
X <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = VCV))
Y <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = VCV))

# Ordinal response
eta <- -0.8 * X # Linear predictor
cutpoints <- c(-1, 0, 1, 2)
Q <- matrix(0, nrow = N, ncol = 4)
for (i in 1:N) {
    for (k in 1:4) {
        Q[i, k] <- 1 / (1 + exp(-(cutpoints[k] - eta[i])))
    }
}
P <- matrix(0, nrow = N, ncol = 5)
P[, 1] <- Q[, 1]
for (k in 2:4) {
    P[, k] <- Q[, k] - Q[, k - 1]
}
P[, 5] <- 1 - Q[, 4]
IUCN <- apply(P, 1, function(p) sample(1:5, 1, prob = p))

data_ordinal <- list(
    X = log(abs(X) + 1),
    Y = log(abs(Y) + 1),
    IUCN = IUCN,
    K_IUCN = 5
)

cat("Running model with ordinal + gaussian...\n")
fit_ordinal <- because(
    data = data_ordinal,
    tree = tree,
    n.burnin = 100,
    n.iter = 500,
    n.chains = 1,
    equations = list(IUCN ~ X, Y ~ X),
    family = c(IUCN = "ordinal") # Y distribution not specified!
)

cat("\n✓ SUCCESS: Model ran without error\n")
summary(fit_ordinal)
