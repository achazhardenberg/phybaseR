library(phybaseR)
library(nimble)
library(rjags)

# Load rhino data
data("rhino.dat")
data("rhino.tree")

data_list <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD
)

equations <- list(
    NL ~ BM,
    DD ~ NL
)

cat("\n=== 3-Way Benchmark: JAGS (Marginal) vs JAGS (Latent) vs NIMBLE ===\n\n")

# ============================================================================
# 1. JAGS (Marginal) - Current Implementation
# ============================================================================
cat("--- 1. JAGS (Marginal/Current) - 1000 iterations ---\n")
time_jags_marginal <- system.time({
    fit_jags_marginal <- phybase_run(
        data = data_list,
        tree = rhino.tree,
        equations = equations,
        n.iter = 1000,
        n.burnin = 200,
        n.thin = 2,
        parallel = FALSE,
        WAIC = FALSE,
        DIC = FALSE,
        quiet = TRUE
    )
})
cat(sprintf(
    "JAGS (Marginal) Time: %.2f seconds\n",
    time_jags_marginal["elapsed"]
))

# ============================================================================
# 2. JAGS (Latent) - Optimized Implementation
# ============================================================================
cat("\n--- 2. JAGS (Latent/Optimized) - 1000 iterations ---\n")

N <- length(rhino.tree$tip.label)
vcv_tree <- ape::vcv(rhino.tree)
Prec_phylo_fixed <- solve(vcv_tree)

# Define JAGS model manually
jags_latent_code <- "
model {
  # Priors
  alphaBM ~ dnorm(0, 0.001)
  alphaNL ~ dnorm(0, 0.001)
  alphaDD ~ dnorm(0, 0.001)
  betaBM ~ dnorm(0, 0.001)
  betaNL ~ dnorm(0, 0.001)
  
  tau_u_BM ~ dgamma(0.001, 0.001)
  tau_e_BM ~ dgamma(0.001, 0.001)
  tau_u_NL ~ dgamma(0.001, 0.001)
  tau_e_NL ~ dgamma(0.001, 0.001)
  tau_u_DD ~ dgamma(0.001, 0.001)
  tau_e_DD ~ dgamma(0.001, 0.001)

  # Latent effects (Standardized)
  u_std_BM[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_NL[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_DD[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])

  for(i in 1:N) {
    # Scale effects
    u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
    u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
    u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)
    
    # Linear predictors
    mu_BM[i] <- alphaBM
    mu_NL[i] <- alphaNL + betaBM * BM[i]
    mu_DD[i] <- alphaDD + betaNL * NL[i]
    
    # Likelihoods
    BM[i] ~ dnorm(mu_BM[i] + u_BM[i], tau_e_BM)
    NL[i] ~ dnorm(mu_NL[i] + u_NL[i], tau_e_NL)
    DD[i] ~ dnorm(mu_DD[i] + u_DD[i], tau_e_DD)
  }
}
"

jags_data <- c(
    data_list,
    list(N = N, Prec_phylo_fixed = Prec_phylo_fixed, zeros = rep(0, N))
)
jags_inits <- list(
    alphaBM = 0,
    alphaNL = 0,
    alphaDD = 0,
    betaBM = 0,
    betaNL = 0,
    tau_u_BM = 1,
    tau_e_BM = 1,
    tau_u_NL = 1,
    tau_e_NL = 1,
    tau_u_DD = 1,
    tau_e_DD = 1
)

time_jags_latent <- system.time({
    model_jags <- jags.model(
        textConnection(jags_latent_code),
        data = jags_data,
        inits = jags_inits,
        n.chains = 1,
        quiet = TRUE
    )
    update(model_jags, 200) # burnin
    samples_jags <- coda.samples(
        model_jags,
        variable.names = c("betaBM", "betaNL"),
        n.iter = 1000,
        thin = 2
    )
})
cat(sprintf("JAGS (Latent) Time: %.2f seconds\n", time_jags_latent["elapsed"]))

# ============================================================================
# 3. NIMBLE (Latent) - Reference
# ============================================================================
cat("\n--- 3. NIMBLE (Latent) - 1000 iterations ---\n")

nimble_code <- nimbleCode({
    alphaBM ~ dnorm(0, 0.001)
    alphaNL ~ dnorm(0, 0.001)
    alphaDD ~ dnorm(0, 0.001)
    betaBM ~ dnorm(0, 0.001)
    betaNL ~ dnorm(0, 0.001)

    tau_u_BM ~ dgamma(0.001, 0.001)
    tau_e_BM ~ dgamma(0.001, 0.001)
    tau_u_NL ~ dgamma(0.001, 0.001)
    tau_e_NL ~ dgamma(0.001, 0.001)
    tau_u_DD ~ dgamma(0.001, 0.001)
    tau_e_DD ~ dgamma(0.001, 0.001)

    u_std_BM[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_NL[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_DD[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])

    for (i in 1:N) {
        u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
        u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
        u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)

        mu_BM[i] <- alphaBM
        mu_NL[i] <- alphaNL + betaBM * BM[i]
        mu_DD[i] <- alphaDD + betaNL * NL[i]

        BM[i] ~ dnorm(mu_BM[i] + u_BM[i], tau_e_BM)
        NL[i] ~ dnorm(mu_NL[i] + u_NL[i], tau_e_NL)
        DD[i] ~ dnorm(mu_DD[i] + u_DD[i], tau_e_DD)
    }
})

time_nimble <- system.time({
    nimble_model <- nimbleModel(
        code = nimble_code,
        constants = list(
            N = N,
            Prec_phylo_fixed = Prec_phylo_fixed,
            zeros = rep(0, N)
        ),
        data = data_list,
        inits = jags_inits
    )
    c_nimble <- compileNimble(nimble_model)
    mcmc_conf <- configureMCMC(nimble_model, monitors = c("betaBM", "betaNL"))
    mcmc <- buildMCMC(mcmc_conf)
    c_mcmc <- compileNimble(mcmc, project = nimble_model)
    nimble_samples <- runMCMC(
        c_mcmc,
        niter = 1000,
        nburnin = 200,
        thin = 2,
        nchains = 1,
        setSeed = 123
    )
})
cat(sprintf("NIMBLE (Total) Time: %.2f seconds\n", time_nimble["elapsed"]))

cat("\n=== SUMMARY ===\n")
cat(sprintf("JAGS (Marginal): %.2f s\n", time_jags_marginal["elapsed"]))
cat(sprintf("JAGS (Latent):   %.2f s\n", time_jags_latent["elapsed"]))
cat(sprintf("NIMBLE (Total):  %.2f s\n", time_nimble["elapsed"]))
cat("\n")
cat(sprintf(
    "Speedup (JAGS Latent vs Marginal): %.1fx\n",
    time_jags_marginal["elapsed"] / time_jags_latent["elapsed"]
))
