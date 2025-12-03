library(phybaseR)
library(nimble)

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

cat("\n=== JAGS vs NIMBLE: Corrected Model Benchmark ===\n\n")

# ============================================================================
# 1. Run with JAGS (Marginal Model)
# ============================================================================
cat("--- Running JAGS (Marginal Model, 2000 iterations) ---\n")
time_jags <- system.time({
    fit_jags <- phybase_run(
        data = data_list,
        tree = rhino.tree,
        equations = equations,
        n.iter = 2000,
        n.burnin = 500,
        n.thin = 5,
        parallel = FALSE,
        WAIC = FALSE,
        DIC = TRUE,
        quiet = TRUE
    )
})
cat(sprintf("JAGS Time: %.2f seconds\n", time_jags["elapsed"]))

# ============================================================================
# 2. Build NIMBLE Model (Latent Variable Equivalent)
# ============================================================================
cat("\n--- Building NIMBLE Model (Latent Variable Equivalent) ---\n")

N <- length(rhino.tree$tip.label)
vcv_tree <- ape::vcv(rhino.tree)
# Pre-compute inverse VCV once!
Prec_phylo_fixed <- solve(vcv_tree)

nimble_code <- nimbleCode({
    # Priors for intercepts and slopes
    alphaBM ~ dnorm(0, 0.001)
    alphaNL ~ dnorm(0, 0.001)
    alphaDD ~ dnorm(0, 0.001)
    betaBM ~ dnorm(0, 0.001)
    betaNL ~ dnorm(0, 0.001)

    # Priors for variance components (equivalent to tau/lambda)
    # tau_u = precision of phylogenetic effect
    # tau_e = precision of residual error
    tau_u_BM ~ dgamma(0.001, 0.001)
    tau_e_BM ~ dgamma(0.001, 0.001)

    tau_u_NL ~ dgamma(0.001, 0.001)
    tau_e_NL ~ dgamma(0.001, 0.001)

    tau_u_DD ~ dgamma(0.001, 0.001)
    tau_e_DD ~ dgamma(0.001, 0.001)

    # Latent phylogenetic effects
    # u ~ N(0, (1/tau_u) * V)  <==> u ~ dmnorm(0, prec = tau_u * V_inv)
    # We use the scaled precision matrix: tau_u * Prec_phylo_fixed

    # Note: In NIMBLE, we can't easily multiply matrix by scalar in distribution parameter
    # So we define the precision matrix node or use a loop.
    # Efficient way: u[1:N] ~ dmnorm(zeros[1:N], prec = Prec_u_BM[1:N, 1:N])
    # But updating Prec_u_BM every step is O(N^2).
    # BETTER: Standardize u_raw ~ N(0, I) and transform?
    # u = L * u_raw where L is Cholesky of V.
    # Then u ~ N(0, V). Then scale by sigma_u.
    # u_final = (1/sqrt(tau_u)) * L * u_raw.
    # This is O(N^2) for matrix-vector mult.

    # Let's try the direct dmnorm approach first, NIMBLE might optimize it.
    # Actually, if we use the formulation:
    # u_std[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    # u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
    # This gives u_BM ~ N(0, (1/tau_u_BM) * V). Correct!

    u_std_BM[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_NL[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_DD[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])

    for (i in 1:N) {
        # Scale the phylogenetic effects
        u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
        u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
        u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)

        # Linear predictors
        mu_BM[i] <- alphaBM
        mu_NL[i] <- alphaNL + betaBM * BM[i]
        mu_DD[i] <- alphaDD + betaNL * NL[i]

        # Likelihoods (Residual error)
        # y ~ N(mu + u, 1/tau_e)
        BM[i] ~ dnorm(mu_BM[i] + u_BM[i], tau_e_BM)
        NL[i] ~ dnorm(mu_NL[i] + u_NL[i], tau_e_NL)
        DD[i] ~ dnorm(mu_DD[i] + u_DD[i], tau_e_DD)
    }

    # Derived parameters (to match JAGS output)
    # lambda = sigma_u^2 / (sigma_u^2 + sigma_e^2)
    # sigma_u^2 = 1/tau_u
    # sigma_e^2 = 1/tau_e
    lambda_BM <- (1 / tau_u_BM) / ((1 / tau_u_BM) + (1 / tau_e_BM))
    lambda_NL <- (1 / tau_u_NL) / ((1 / tau_u_NL) + (1 / tau_e_NL))
    lambda_DD <- (1 / tau_u_DD) / ((1 / tau_u_DD) + (1 / tau_e_DD))
})

time_compile <- system.time({
    nimble_model <- nimbleModel(
        code = nimble_code,
        constants = list(
            N = N,
            Prec_phylo_fixed = Prec_phylo_fixed,
            zeros = rep(0, N)
        ),
        data = data_list,
        inits = list(
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
            tau_e_DD = 1,
            u_std_BM = rep(0, N),
            u_std_NL = rep(0, N),
            u_std_DD = rep(0, N)
        )
    )
    c_nimble <- compileNimble(nimble_model)
    mcmc_conf <- configureMCMC(
        nimble_model,
        monitors = c("betaBM", "betaNL", "lambda_BM", "lambda_NL", "lambda_DD")
    )
    mcmc <- buildMCMC(mcmc_conf)
    c_mcmc <- compileNimble(mcmc, project = nimble_model)
})
cat(sprintf("Compilation time: %.2f seconds\n", time_compile["elapsed"]))

# ============================================================================
# 3. Run NIMBLE
# ============================================================================
cat("\n--- Running NIMBLE (2000 iterations) ---\n")
time_nimble <- system.time({
    nimble_samples <- runMCMC(
        c_mcmc,
        niter = 2000,
        nburnin = 500,
        thin = 5,
        nchains = 1,
        setSeed = 123
    )
})
cat(sprintf("NIMBLE Time (sampling): %.2f seconds\n", time_nimble["elapsed"]))

# ============================================================================
# 4. Results
# ============================================================================
cat("\n=== RESULTS ===\n")
cat(sprintf("JAGS Total:        %.2f seconds\n", time_jags["elapsed"]))
cat(sprintf("NIMBLE Compile:    %.2f seconds\n", time_compile["elapsed"]))
cat(sprintf("NIMBLE Sample:     %.2f seconds\n", time_nimble["elapsed"]))
cat(sprintf(
    "NIMBLE Total:      %.2f seconds\n",
    time_nimble["elapsed"] + time_compile["elapsed"]
))
cat("\n")
cat(sprintf(
    "Speedup (sampling only): %.1fx\n",
    time_jags["elapsed"] / time_nimble["elapsed"]
))

# Check estimates
cat("\n--- Parameter Estimates (NIMBLE) ---\n")
print(round(apply(nimble_samples, 2, mean), 4))
