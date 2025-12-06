# Script to run test and capture JAGS model
Sys.setenv(NOT_CRAN = 'true')
devtools::load_all()
library(testthat)

# Run the test but catch errors and save model files
set.seed(123)
N <- 100
K <- 3
tree <- ape::rtree(N)
C <- ape::vcv(tree)

X <- rnorm(N)

# Simulate data (same as test)
phy1 <- MASS::mvrnorm(1, mu = rep(0, N), Sigma = 0.8 * C + diag(0.1, N))
L1 <- 0.5 + 1.2 * X + phy1

phy2 <- MASS::mvrnorm(1, mu = rep(0, N), Sigma = 0.5 * C + diag(0.1, N))
L2 <- -0.5 - 0.8 * X + phy2

P <- matrix(0, N, K)
denom <- 1 + exp(L1) + exp(L2)
P[, 1] <- 1 / denom
P[, 2] <- exp(L1) / denom
P[, 3] <- exp(L2) / denom

Y <- numeric(N)
for (i in 1:N) {
    Y[i] <- sample(1:K, 1, prob = P[i, ])
}
Y <- factor(Y)

data_list <- list(Y = Y, X = X)

tryCatch(
    {
        res <- phybase_run(
            equations = list(Y ~ X),
            data = data_list,
            tree = tree,
            distribution = c(Y = "multinomial"),
            n.iter = 1000,
            n.burnin = 500,
            n.chains = 2,
            quiet = TRUE
        )
    },
    error = function(e) {
        cat("\n=== TEST ERROR ===\n")
        cat(e$message, "\n")

        # Find latest model file
        model_files <- list.files(
            tempdir(),
            pattern = "\\.jg$",
            full.names = TRUE
        )
        if (length(model_files) > 0) {
            latest <- model_files[which.max(file.mtime(model_files))]
            cat("\n=== JAGS MODEL FROM TEST ===\n")
            cat(readLines(latest), sep = "\n")
            cat("\n=== END MODEL ===\n")
        }
    }
)
