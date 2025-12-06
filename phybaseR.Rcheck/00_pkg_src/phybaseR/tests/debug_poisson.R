library(phybaseR)

# Reproduce Poisson bug
set.seed(123)
N <- 35
df <- data.frame(
    SP = paste0("sp", 1:N),
    Y = rpois(N, lambda = 5),
    X = rnorm(N)
)

cat("Testing Poisson model...\n")
fit <- phybase_run(
    data = df,
    id_col = "SP",
    equations = list(Y ~ X),
    distribution = list(Y = "poisson"),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 500,
    quiet = FALSE
)

cat("\nModel code:\n")
cat(fit$model_code)
