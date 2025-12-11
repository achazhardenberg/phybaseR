# Replicate MCMCglmm SShorns Example with Because

library(because)
if (require(MCMCglmm)) {
    data(SShorns)
} else {
    stop("MCMCglmm package required for SShorns data")
}

# Data Inspection
head(SShorns)
table(SShorns$horn, SShorns$sex)

# ==============================================================================
# Model 1: Intercept Only
# MCMCglmm: horn ~ trait - 1, family = "categorical"
# Because:  horn ~ 1, distribution = "multinomial"
# ==============================================================================

message("\n--- Fits Model 1: Intercept Only ---")
eq_1 <- list(horn ~ 1)

fit_1 <- because(
    equations = eq_1,
    data = SShorns,
    distribution = "multinomial", # Will auto-fix residual variance to 1
    n.chains = 3,
    n.iter = 5000,
    quiet = FALSE
)

summary(fit_1)

# Comparison Note:
# MCMCglmm estimates are typically Probit or requiring rescaling (c2 constant) to match Logit.
# Because uses Logit link (dbern/dcat via latent logit).
# However, if MCMCglmm uses "categorical" family, it often uses a specific link/latent variable formulation.
# The user prompt mentions "rescale the intercepts as if the residual covariance matrix was zero" using c2 = (16 * sqrt(3)/(15 * pi))^2 for comparison.
# This implies MCMCglmm might be using a Logit approximation or Probit.
# Because implementation uses standard JAGS dmnorm/dcat or similar.
# Wait, because_model.R for multinomial uses:
#   Multinomial: K categories
#   L[i, 1] <- 0 (Baseline)
#   L[i, k] <- alpha[k] + beta*X
#   p[i, 1:K] <- exp(L[i, 1:K]) / sum(exp(L))
#   y[i] ~ dcat(p[i, 1:K])
# This is a standard Multinomial Logit choice.
# MCMCglmm "categorical" uses a latent variable threshold formulation (Probit usually) or Logit.
# The prompt mentions "residual variance is not identified... set to... I + J".
# Because with 'fix_residual_variance=1' sets tau_e=1 (var=1).
# We'll see how the estimates compare.

# ==============================================================================
# Model 2: Sex Effect
# MCMCglmm: horn ~ trait + sex - 1
# Because:  horn ~ sex
# ==============================================================================

message("\n--- Fits Model 2: Sex Effect ---")
eq_2 <- list(horn ~ sex)

fit_2 <- because(
    equations = eq_2,
    data = SShorns,
    distribution = "multinomial",
    n.chains = 3,
    n.iter = 5000,
    quiet = TRUE
)

summary(fit_2)
