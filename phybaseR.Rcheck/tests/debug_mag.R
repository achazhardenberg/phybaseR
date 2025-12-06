library(phybaseR)
library(ape)

# Reproduce MAG bug
set.seed(123)
N <- 35
tree <- rtree(N)

df <- data.frame(
    SP = tree$tip.label,
    Y1 = rnorm(N, mean = 10, sd = 2),
    Y2 = rnorm(N, mean = 8, sd = 1.5)
)

cat("Testing MAG model with latent variable...\n")
fit <- phybase_run(
    data = df,
    structure = tree,
    id_col = "SP",
    latent = "L1",
    equations = list(
        Y1 ~ L1,
        Y2 ~ L1
    ),
    latent_method = "correlations",
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 500,
    quiet = FALSE
)

cat("\nSuccess!\n")
