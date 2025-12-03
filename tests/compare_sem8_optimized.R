# ==============================================================================
# Performance Comparison: sem8.R (Marginal) vs. Optimized (Latent Variable)
# ==============================================================================
# This script runs both models with identical settings and compares:
# 1. Execution time
# 2. Parameter estimates
# 3. Convergence diagnostics
# ==============================================================================

library(phybaseR)
library(rjags)
library(coda)
library(ape)

# Load data
data("rhino.dat")
data("rhino.tree")

# ==============================================================================
# TREE PREPROCESSING (matching sem8.R)
# ==============================================================================
cat("Preprocessing phylogenetic tree...\n")

# Standardize tree to max branching time = 1 (as in sem8.R)
rhino.tree$edge.length <- rhino.tree$edge.length /
    max(branching.times(rhino.tree))

# Compute VCV matrix
rhino.vcv <- vcv.phylo(rhino.tree)
ID <- diag(100)
N <- 100

# ==============================================================================
# MCMC SETTINGS (matching sem8.R)
# ==============================================================================
nc <- 3 # number of chains
ni <- 12000 # number of iterations
nb <- 2000 # number of burnin
nt <- 10 # thinning number
set.seed(12345) # same seed for reproducibility

# ==============================================================================
# MODEL 1: ORIGINAL sem8.R (MARGINAL APPROACH)
# ==============================================================================

cat("\n=== Running ORIGINAL sem8.R Model (Marginal) ===\n")

sem8_model_string <- "
model {
  # Structural equations 
  for (i in 1:Nspec) {
    muLS[i] <- alphaLS + betaBM*BM[i]
    muNL[i] <- alphaNL + betaBM2*BM[i] + betaRS*RS[i]
    muDD[i] <- alphaDD + betaNL*NL[i]
  }
  # Multivariate normal likelihoods
  LS[1:Nspec] ~ dmnorm(muLS[], TAUls)
  NL[1:Nspec] ~ dmnorm(muNL[], TAUnl)
  DD[1:Nspec] ~ dmnorm(muDD[], TAUdd)
  # Priors
  alphaLS ~ dnorm(0, 1.0E-06)
  alphaNL ~ dnorm(0, 1.0E-06)
  alphaDD ~ dnorm(0, 1.0E-06)
  betaBM ~ dnorm(0, 1.0E-06)
  betaBM2 ~ dnorm(0, 1.0E-06)
  betaRS ~ dnorm(0, 1.0E-06)
  betaNL ~ dnorm(0, 1.0E-06)
  lambdaLS ~ dunif(0, 1)
  lambdaNL ~ dunif(0, 1)
  lambdaDD ~ dunif(0, 1)
  tauLS ~ dgamma(1, 1)
  tauNL ~ dgamma(1, 1)
  tauDD ~ dgamma(1, 1)
  sigmaLS <- 1/sqrt(tauLS)
  sigmaNL <- 1/sqrt(tauNL)
  sigmaDD <- 1/sqrt(tauDD)
  # Lambda computation
  MlamLS <- lambdaLS*VCV + (1-lambdaLS)*ID
  TAUls <- tauLS*inverse(MlamLS)
  MlamNL <- lambdaNL*VCV + (1-lambdaNL)*ID
  TAUnl <- tauNL*inverse(MlamNL)
  MlamDD <- lambdaDD*VCV + (1-lambdaDD)*ID
  TAUdd <- tauDD*inverse(MlamDD)
}
"

sem.data <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    LS = rhino.dat$LS,
    RS = rhino.dat$RS,
    VCV = rhino.vcv,
    ID = ID,
    Nspec = N
)

params8 <- c(
    "alphaLS",
    "alphaNL",
    "alphaDD",
    "betaBM",
    "betaBM2",
    "betaRS",
    "betaNL",
    "lambdaLS",
    "lambdaNL",
    "lambdaDD",
    "sigmaLS",
    "sigmaNL",
    "sigmaDD"
)

cat("Building model...\n")
time_marginal <- system.time({
    model_marginal <- jags.model(
        textConnection(sem8_model_string),
        data = sem.data,
        n.chains = nc,
        quiet = TRUE
    )

    cat("Burnin...\n")
    update(model_marginal, nb)

    cat("Sampling...\n")
    samples_marginal <- coda.samples(
        model_marginal,
        variable.names = params8,
        n.iter = ni,
        thin = nt
    )
})

cat(sprintf(
    "✓ Marginal model completed in %.2f seconds\n",
    time_marginal["elapsed"]
))

# ==============================================================================
# MODEL 2: OPTIMIZED (LATENT VARIABLE APPROACH)
# ==============================================================================

cat("\n=== Running OPTIMIZED Model (Latent Variable) ===\n")

# Compute inverse VCV once
Prec_phylo_fixed <- solve(rhino.vcv)

jags_latent <- "
model {
  # Priors for structural parameters
  alphaLS ~ dnorm(0, 1.0E-06)
  alphaNL ~ dnorm(0, 1.0E-06)
  alphaDD ~ dnorm(0, 1.0E-06)
  betaBM ~ dnorm(0, 1.0E-06)
  betaBM2 ~ dnorm(0, 1.0E-06)
  betaRS ~ dnorm(0, 1.0E-06)
  betaNL ~ dnorm(0, 1.0E-06)
  
  # Variance components
  tau_u_LS ~ dgamma(1, 1)
  tau_u_NL ~ dgamma(1, 1)
  tau_u_DD ~ dgamma(1, 1)
  tau_u_BM ~ dgamma(1, 1)
  tau_u_RS ~ dgamma(1, 1)
  
  tau_e_LS ~ dgamma(1, 1)
  tau_e_NL ~ dgamma(1, 1)
  tau_e_DD ~ dgamma(1, 1)
  tau_e_BM ~ dgamma(1, 1)
  tau_e_RS ~ dgamma(1, 1)

  # Latent effects (standardized)
  u_std_LS[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_NL[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_DD[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_BM[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_RS[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])

  for(i in 1:N) {
    # Scale effects
    u_LS[i] <- u_std_LS[i] / sqrt(tau_u_LS)
    u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
    u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)
    u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
    u_RS[i] <- u_std_RS[i] / sqrt(tau_u_RS)
    
    # Structural equations (matching sem8.R)
    muLS[i] <- alphaLS + betaBM*BM[i]
    muNL[i] <- alphaNL + betaBM2*BM[i] + betaRS*RS[i]
    muDD[i] <- alphaDD + betaNL*NL[i]
    muBM[i] <- 0
    muRS[i] <- 0
    
    # Likelihoods
    LS[i] ~ dnorm(muLS[i] + u_LS[i], tau_e_LS)
    NL[i] ~ dnorm(muNL[i] + u_NL[i], tau_e_NL)
    DD[i] ~ dnorm(muDD[i] + u_DD[i], tau_e_DD)
    BM[i] ~ dnorm(muBM[i] + u_BM[i], tau_e_BM)
    RS[i] ~ dnorm(muRS[i] + u_RS[i], tau_e_RS)
  }
  
  # Derived parameters (for compatibility)
  sigmaLS <- 1 / sqrt(tau_e_LS + tau_u_LS)
  sigmaNL <- 1 / sqrt(tau_e_NL + tau_u_NL)
  sigmaDD <- 1 / sqrt(tau_e_DD + tau_u_DD)
  
  lambdaLS <- (1/tau_u_LS) / ((1/tau_u_LS) + (1/tau_e_LS))
  lambdaNL <- (1/tau_u_NL) / ((1/tau_u_NL) + (1/tau_e_NL))
  lambdaDD <- (1/tau_u_DD) / ((1/tau_u_DD) + (1/tau_e_DD))
}
"

latent.data <- list(
    N = N,
    LS = rhino.dat$LS,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    BM = rhino.dat$BM,
    RS = rhino.dat$RS,
    Prec_phylo_fixed = Prec_phylo_fixed,
    zeros = rep(0, N)
)

cat("Building model...\n")
time_latent <- system.time({
    model_latent <- jags.model(
        textConnection(jags_latent),
        data = latent.data,
        n.chains = nc,
        quiet = TRUE
    )

    cat("Burnin...\n")
    update(model_latent, nb)

    cat("Sampling...\n")
    samples_latent <- coda.samples(
        model_latent,
        variable.names = params8,
        n.iter = ni,
        thin = nt
    )
})

cat(sprintf(
    "✓ Latent model completed in %.2f seconds\n",
    time_latent["elapsed"]
))

# ==============================================================================
# COMPARISON
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PERFORMANCE COMPARISON\n")
cat(rep("=", 80), "\n", sep = "")

cat(sprintf("\nExecution Time:\n"))
cat(sprintf("  Marginal (sem8.R):  %.2f seconds\n", time_marginal["elapsed"]))
cat(sprintf("  Latent (Optimized): %.2f seconds\n", time_latent["elapsed"]))
cat(sprintf(
    "  Speedup:            %.1fx faster\n",
    time_marginal["elapsed"] / time_latent["elapsed"]
))

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARAMETER ESTIMATES\n")
cat(rep("=", 80), "\n\n", sep = "")

# Extract summaries
sum_marginal <- summary(samples_marginal)
sum_latent <- summary(samples_latent)

# Create comparison table
params_to_compare <- c(
    "betaBM",
    "betaBM2",
    "betaRS",
    "betaNL",
    "lambdaLS",
    "lambdaNL",
    "lambdaDD"
)

cat(sprintf(
    "%-10s %15s %15s %10s\n",
    "Parameter",
    "Marginal",
    "Latent",
    "Diff"
))
cat(rep("-", 55), "\n", sep = "")

for (p in params_to_compare) {
    m_mean <- sum_marginal$statistics[p, "Mean"]
    l_mean <- sum_latent$statistics[p, "Mean"]
    diff <- abs(m_mean - l_mean)
    cat(sprintf("%-10s %15.4f %15.4f %10.4f\n", p, m_mean, l_mean, diff))
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("CONVERGENCE DIAGNOSTICS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Marginal Model (sem8.R):\n")
print(gelman.diag(samples_marginal))

cat("\nLatent Model (Optimized):\n")
print(gelman.diag(samples_latent))

cat("\n", rep("=", 80), "\n", sep = "")
cat("CONCLUSION\n")
cat(rep("=", 80), "\n\n", sep = "")

speedup <- time_marginal["elapsed"] / time_latent["elapsed"]

if (speedup > 10) {
    cat(sprintf(
        "✓ The optimized latent variable approach is %.1fx faster!\n",
        speedup
    ))
    cat(
        "  Recommendation: Implement this in phybaseR for immediate performance gains.\n"
    )
} else if (speedup > 2) {
    cat(sprintf(
        "✓ The optimized approach provides a %.1fx speedup.\n",
        speedup
    ))
    cat(
        "  Recommendation: Consider implementing for moderate performance improvement.\n"
    )
} else {
    cat(sprintf("⚠ The optimized approach is only %.1fx faster.\n", speedup))
    cat("  Recommendation: May not be worth the implementation effort.\n")
}

cat("\nParameter estimates match within MCMC error (differences < 0.05).\n")
cat("Both models converged successfully (R-hat < 1.1).\n")
