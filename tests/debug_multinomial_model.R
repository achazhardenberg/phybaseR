# Debug script to see the JAGS model code for multinomial
library(phybaseR)
library(ape)
library(MASS)

set.seed(123)
N <- 100 # Match test conditions
tree <- rtree(N)
C <- vcv(tree)
X <- rnorm(N)

# Simulate Y with K=3 categories
phy1 <- mvrnorm(1, mu = rep(0, N), Sigma = 0.8 * C + diag(0.1, N))
phy2 <- mvrnorm(1, mu = rep(0, N), Sigma = 0.5 * C + diag(0.1, N))
L1 <- 0.5 + 1.2 * X + phy1
L2 <- -0.5 - 0.8 * X + phy2

P <- matrix(0, N, 3)
denom <- 1 + exp(L1) + exp(L2)
P[, 1] <- 1 / denom
P[, 2] <- exp(L1) / denom
P[, 3] <- exp(L2) / denom

Y <- numeric(N)
for (i in 1:N) {
    Y[i] <- sample(1:3, 1, prob = P[i, ]) # K=3 categories
}
Y <- factor(Y)

data_list <- list(Y = Y, X = X)

# Try to run - it will fail, but we can catch the model code
tryCatch(
    {
        res <- phybase_run(
            equations = list(Y ~ X),
            data = data_list,
            tree = tree,
            distribution = c(Y = "multinomial"),
            n.iter = 100,
            n.burnin = 50,
            n.chains = 1,
            quiet = FALSE
        )
    },
    error = function(e) {
        cat("\n=== ERROR ===\n")
        cat(e$message, "\n")

        # Try to read the model file if it was created
        model_files <- list.files(
            tempdir(),
            pattern = "\\.jg$",
            full.names = TRUE
        )
        if (length(model_files) > 0) {
            latest <- model_files[which.max(file.mtime(model_files))]
            cat("\n=== JAGS MODEL CODE (from", basename(latest), ") ===\n")
            cat(readLines(latest), sep = "\n")
            cat("\n=== END MODEL CODE ===\n")
        }
    }
)
