library(because)
set.seed(42)
N <- 400
J <- 3

psi_Pred <- 0.4
z_Pred <- rbinom(N, 1, psi_Pred)
p_Pred <- 0.7
y_Pred <- matrix(NA, N, J)
for (i in 1:N) {
    for (j in 1:J) {
        y_Pred[i, j] <- rbinom(1, 1, z_Pred[i] * p_Pred)
    }
}

psi_Prey <- 0.6
z_Prey <- rbinom(N, 1, psi_Prey)
alpha_p_Prey <- 1.0
beta_p_Prey_Pred <- -2.0
p_Prey <- plogis(alpha_p_Prey + beta_p_Prey_Pred * z_Pred)
y_Prey <- matrix(NA, N, J)
for (i in 1:N) {
    for (j in 1:J) {
        y_Prey[i, j] <- rbinom(1, 1, z_Prey[i] * p_Prey[i])
    }
}

res <- because(
    equations = list(
        Predator ~ 1,
        p_Predator ~ 1,
        Prey ~ 1,
        p_Prey ~ Predator
    ),
    data = list(Predator = y_Pred, Prey = y_Prey, N = N),
    family = c(Predator = "occupancy", Prey = "occupancy"),
    priors = list(
        beta_p_Prey_Predator = "dnorm(0, 0.5)",
        alpha_p_Prey = "dnorm(0, 0.5)"
    ),
    n.iter = 5000,
    quiet = TRUE
)

est <- summary(res)$results["beta_p_Prey_Predator", "Mean"]
cat("Estimated Beta:", est, "\n")
