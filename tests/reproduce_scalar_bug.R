# Scalar Bug Reproduction
if (requireNamespace("devtools")) {
    devtools::load_all(".")
} else {
    library(because)
}

set.seed(42)
N <- 50
Y <- matrix(rbinom(N * 2, 1, 0.5), ncol = 2)

# Case A: Explicit N
message("\n--- CASE A: Explicit N ---")
data_A <- list(Y = Y, N = N)
res_A <- because(
    equations = list(Y ~ 1, p_Y ~ 1),
    family = c(Y = "occupancy"),
    data = data_A,
    quiet = TRUE,
    n.chains = 1,
    n.iter = 50
)
params_A <- rownames(res_A$summary$statistics)
psi_A <- grep("psi", params_A, value = TRUE)
n_psi_A <- length(psi_A)
message("A Params Count: ", n_psi_A)
if (n_psi_A == 50) {
    message("SUCCESS: Case A has 50 psi params.")
} else {
    message("FAILURE: Case A has ", n_psi_A, " psi params.")
}

# Case B: Implicit N (via X)
message("\n--- CASE B: Implicit N (via X) ---")
data_B <- list(Y = Y, X = rnorm(N))
res_B <- because(
    equations = list(Y ~ 1, p_Y ~ 1),
    family = c(Y = "occupancy"),
    data = data_B,
    quiet = TRUE,
    n.chains = 1,
    n.iter = 50
)
params_B <- rownames(res_B$summary$statistics)
n_psi_B <- length(grep("psi", params_B))
message("B Params Count: ", n_psi_B)
