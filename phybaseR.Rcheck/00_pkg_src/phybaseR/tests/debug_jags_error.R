library(rjags)
library(ape)

set.seed(123)
N <- 50
tree <- rtree(N)
VCV <- vcv(tree)
ID <- diag(N)
Y <- rnorm(N)
X <- rnorm(N)

data <- list(
  Y = Y,
  X = X,
  N = N,
  VCV = VCV,
  ID = ID
)

model_string <- "
model {
  # Structural equations
  for (i in 1:N) {
    muY[i] <- alphaY + beta_Y_X*X[i]
  }
  # Multivariate normal likelihoods
  Y[1:N] ~ dmnorm(muY[], TAUy)
  # Priors for structural parameters
  alphaY ~ dnorm(0, 1.0E-6)
  lambdaY ~ dunif(0, 1)
  tauY ~ dgamma(1, 1)
  beta_Y_X ~ dnorm(0, 1.0E-6)
  
  # Covariance structure for responses
  MlamY <- lambdaY*VCV + (1-lambdaY)*ID
  TAUy <- tauY*inverse(MlamY)
  
  # Predictor priors for imputation (X)
  for (i in 1:N) {
    muX[i] <- 0
  }
  X[1:N] ~ dmnorm(muX[], TAUx)
  lambdaX ~ dunif(0, 1)
  tauX ~ dgamma(1, 1)
  MlamX <- lambdaX*VCV + (1 - lambdaX)*ID
  TAUx <- tauX*inverse(MlamX)
}
"

cat("Compiling model...\n")
model <- jags.model(textConnection(model_string), data = data, n.chains = 1)
cat("Success!\n")
