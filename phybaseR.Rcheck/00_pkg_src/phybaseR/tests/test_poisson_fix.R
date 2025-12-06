library(phybaseR)

# Quick Poisson test
set.seed(123)
N <- 35
df <- data.frame(
    SP = paste0("sp", 1:N),
    Y = rpois(N, lambda = 5),
    X = rnorm(N)
)

cat("Testing Poisson model with WAIC...\n")
fit <- phybase_run(
    data = df,
    id_col = "SP",
    equations = list(Y ~ X),
    distribution = list(Y = "poisson"),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 1000,
    quiet = TRUE
)

cat("\n✓ SUCCESS! Poisson model compiles\n\n")
cat("WAIC Results:\n")
print(fit$WAIC)
