# Prototype for Multinomial Phylogenetic Model
# Goal: Verify parameter recovery for a 3-category response using latent variables

library(rjags)
library(ape)
library(MASS)

set.seed(123)

# 1. Simulate Data
N <- 100
K <- 3 # Number of categories
tree <- rtree(N)
C <- vcv(tree)

# Predictor
X <- rnorm(N)

# Parameters for K-1 latent variables (Category 1 is reference)
# Latent 1 (Cat 2 vs Cat 1)
beta1_0 <- 0.5
beta1_1 <- 1.2
lambda1 <- 0.8 # Phylogenetic signal

# Latent 2 (Cat 3 vs Cat 1)
beta2_0 <- -0.5
beta2_1 <- -0.8
lambda2 <- 0.5

# Simulate latent phylogenetic effects
# Independent evolution for now
phy1 <- mvrnorm(1, mu = rep(0, N), Sigma = lambda1 * C + diag(0.1, N))
phy2 <- mvrnorm(1, mu = rep(0, N), Sigma = lambda2 * C + diag(0.1, N))

# Linear predictors
L1 <- beta1_0 + beta1_1 * X + phy1
L2 <- beta2_0 + beta2_1 * X + phy2

# Softmax to probabilities
# L_ref (Cat 1) = 0
P <- matrix(0, N, K)
denom <- 1 + exp(L1) + exp(L2)
P[, 1] <- 1 / denom
P[, 2] <- exp(L1) / denom
P[, 3] <- exp(L2) / denom

# Sample categories
Y <- numeric(N)
for (i in 1:N) {
    Y[i] <- sample(1:K, 1, prob = P[i, ])
}

# 2. JAGS Model
model_string <- "
model {
  # Priors
  beta1_0 ~ dnorm(0, 0.01)
  beta1_1 ~ dnorm(0, 0.01)
  beta2_0 ~ dnorm(0, 0.01)
  beta2_1 ~ dnorm(0, 0.01)
  
  # Precision for phylogenetic effects
  tau1 ~ dgamma(1, 1)
  tau2 ~ dgamma(1, 1)
  
  # Latent phylogenetic effects
  # Using multivariate normal for phylogenetic structure
  # For efficiency, we can use the inverse VCV matrix
  phy1[1:N] ~ dmnorm(mu_zeros, tau1 * invC)
  phy2[1:N] ~ dmnorm(mu_zeros, tau2 * invC)
  
  for (i in 1:N) {
    # Latent linear predictors
    L[i, 1] <- 0 # Reference category
    L[i, 2] <- beta1_0 + beta1_1 * X[i] + phy1[i]
    L[i, 3] <- beta2_0 + beta2_1 * X[i] + phy2[i]
    
    # Softmax
    for (k in 1:K) {
      exp_L[i, k] <- exp(L[i, k])
    }
    sum_exp_L[i] <- sum(exp_L[i, 1:K])
    
    for (k in 1:K) {
      p[i, k] <- exp_L[i, k] / sum_exp_L[i]
    }
    
    # Likelihood
    Y[i] ~ dcat(p[i, 1:K])
  }
}
"

# 3. Fit Model
data_jags <- list(
    Y = Y,
    X = X,
    N = N,
    K = K,
    invC = solve(C),
    mu_zeros = rep(0, N)
)

model <- jags.model(
    textConnection(model_string),
    data = data_jags,
    n.chains = 3
)
update(model, 1000)
samples <- coda.samples(
    model,
    c("beta1_0", "beta1_1", "beta2_0", "beta2_1"),
    n.iter = 3000
)

# 4. Check Recovery
sum_stats <- summary(samples)
print(sum_stats$statistics)

# True values
cat("\nTrue Values:\n")
cat("beta1_0:", beta1_0, "\n")
cat("beta1_1:", beta1_1, "\n")
cat("beta2_0:", beta2_0, "\n")
cat("beta2_1:", beta2_1, "\n")

# Check if true values are within 95% CI
ci <- sum_stats$quantiles[, c(1, 5)]
print(ci)
