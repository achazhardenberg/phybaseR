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

cat("\n=== JAGS vs NIMBLE Performance Benchmark ===\n\n")

# ============================================================================
# 1. Run with JAGS
# ============================================================================
cat("--- Running JAGS (2000 iterations) ---\n")
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
# 2. Build NIMBLE model
# ============================================================================
cat("\n--- Building NIMBLE Model ---\n")

N <- length(rhino.tree$tip.label)
vcv_tree <- ape::vcv(rhino.tree)
Prec_phylo <- solve(vcv_tree)

nimble_code <- nimbleCode({
    for (i in 1:N) {
        mu_phylo_BM[i] <- alphaBM
        mu_phylo_NL[i] <- alphaNL + betaBM * BM[i]
    }

    BM_phylo[1:N] ~ dmnorm(mu_phylo_BM[1:N], prec = Prec_phylo_BM[1:N, 1:N])
    NL_phylo[1:N] ~ dmnorm(mu_phylo_NL[1:N], prec = Prec_phylo_NL[1:N, 1:N])

    for (i in 1:N) {
        BM[i] ~ dnorm(BM_phylo[i], tau_res_BM)
        NL[i] ~ dnorm(NL_phylo[i], tau_res_NL)

        mu_DD[i] <- alphaDD + betaNL * NL[i]
        DD_phylo[i] ~ dnorm(mu_DD[i], tau_phylo_DD)
        DD[i] ~ dnorm(DD_phylo[i], tau_res_DD)
    }

    alphaBM ~ dnorm(0, 0.001)
    alphaNL ~ dnorm(0, 0.001)
    alphaDD ~ dnorm(0, 0.001)
    betaBM ~ dnorm(0, 0.001)
    betaNL ~ dnorm(0, 0.001)

    tau_phylo_BM ~ dgamma(0.001, 0.001)
    tau_phylo_NL ~ dgamma(0.001, 0.001)
    tau_phylo_DD ~ dgamma(0.001, 0.001)
    tau_res_BM ~ dgamma(0.001, 0.001)
    tau_res_NL ~ dgamma(0.001, 0.001)
    tau_res_DD ~ dgamma(0.001, 0.001)
})

time_compile <- system.time({
    nimble_model <- nimbleModel(
        code = nimble_code,
        constants = list(
            N = N,
            Prec_phylo_BM = Prec_phylo,
            Prec_phylo_NL = Prec_phylo
        ),
        data = data_list,
        inits = list(
            alphaBM = 0,
            alphaNL = 0,
            alphaDD = 0,
            betaBM = 0,
            betaNL = 0,
            tau_phylo_BM = 1,
            tau_phylo_NL = 1,
            tau_phylo_DD = 1,
            tau_res_BM = 1,
            tau_res_NL = 1,
            tau_res_DD = 1,
            BM_phylo = rep(0, N),
            NL_phylo = rep(0, N),
            DD_phylo = rep(0, N)
        )
    )
    c_nimble <- compileNimble(nimble_model)
    mcmc_conf <- configureMCMC(nimble_model, print = FALSE)
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
# 4. Summary
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
cat(sprintf(
    "Speedup (total):         %.1fx\n",
    time_jags["elapsed"] / (time_nimble["elapsed"] + time_compile["elapsed"])
))
