library(phybaseR)
library(nimble)
library(rjags)

# Load rhino data
data("rhino.dat")
data("rhino.tree")

# Prepare data
data_list <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD
)

equations <- list(
    NL ~ BM,
    DD ~ NL
)

cat("\n=== BENCHMARK: JAGS vs NIMBLE ===\n\n")

# ============================================================================
# 1. Run with JAGS (current implementation)
# ============================================================================
cat("--- Running JAGS ---\n")
time_jags <- system.time({
    fit_jags <- phybase_run(
        data = data_list,
        tree = rhino.tree,
        equations = equations,
        n.iter = 5000,
        n.burnin = 1000,
        n.thin = 5,
        parallel = FALSE,
        WAIC = FALSE,
        DIC = TRUE,
        quiet = TRUE
    )
})

cat("JAGS Time:", time_jags["elapsed"], "seconds\n")
cat("JAGS DIC:", fit_jags$DIC$deviance + fit_jags$DIC$penalty, "\n\n")

# ============================================================================
# 2. Build equivalent NIMBLE model
# ============================================================================
cat("--- Building NIMBLE Model ---\n")

# Get unique tip labels and create mapping
N <- length(rhino.tree$tip.label)

# Build NIMBLE model code
nimble_code <- nimbleCode({
    # Phylogenetic precision matrices (will be provided as constants)
    # Prec_phylo_BM[1:N, 1:N] and Prec_phylo_NL[1:N, 1:N]

    # Mean vectors for phylogenetic component
    for (i in 1:N) {
        mu_phylo_BM[i] <- alphaBM
        mu_phylo_NL[i] <- alphaNL + betaBM * BM[i]
    }

    # Phylogenetic component (latent variables)
    BM_phylo[1:N] ~ dmnorm(mu_phylo_BM[1:N], prec = Prec_phylo_BM[1:N, 1:N])
    NL_phylo[1:N] ~ dmnorm(mu_phylo_NL[1:N], prec = Prec_phylo_NL[1:N, 1:N])

    # Residual component
    for (i in 1:N) {
        BM[i] ~ dnorm(BM_phylo[i], tau_res_BM)
        NL[i] ~ dnorm(NL_phylo[i], tau_res_NL)
    }

    # Response variable (DD)
    for (i in 1:N) {
        mu_DD[i] <- alphaDD + betaNL * NL[i]
        DD_phylo[i] ~ dnorm(mu_DD[i], tau_phylo_DD)
    }

    for (i in 1:N) {
        DD[i] ~ dnorm(DD_phylo[i], tau_res_DD)
    }

    # Priors
    alphaBM ~ dnorm(0, 0.001)
    alphaNL ~ dnorm(0, 0.001)
    alphaDD ~ dnorm(0, 0.001)
    betaBM ~ dnorm(0, 0.001)
    betaNL ~ dnorm(0, 0.001)

    # Precision priors (convert to tau = 1/sigma^2)
    tau_phylo_BM ~ dgamma(0.001, 0.001)
    tau_phylo_NL ~ dgamma(0.001, 0.001)
    tau_phylo_DD ~ dgamma(0.001, 0.001)
    tau_res_BM ~ dgamma(0.001, 0.001)
    tau_res_NL ~ dgamma(0.001, 0.001)
    tau_res_DD ~ dgamma(0.001, 0.001)

    # Phylogenetic signal (lambda-like parameters, derived)
    lambda_BM <- tau_phylo_BM / (tau_phylo_BM + tau_res_BM)
    lambda_NL <- tau_phylo_NL / (tau_phylo_NL + tau_res_NL)
    lambda_DD <- tau_phylo_DD / (tau_phylo_DD + tau_res_DD)
})

# Prepare phylogenetic precision matrices
library(ape)
vcv_tree <- vcv(rhino.tree)
Prec_phylo <- solve(vcv_tree)

# Prepare NIMBLE data and constants
nimble_constants <- list(
    N = N,
    Prec_phylo_BM = Prec_phylo,
    Prec_phylo_NL = Prec_phylo,
    Prec_phylo_DD = Prec_phylo
)

nimble_data <- list(
    BM = data_list$BM,
    NL = data_list$NL,
    DD = data_list$DD
)

# Initial values
nimble_inits <- list(
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

cat("Compiling NIMBLE model... (this may take 10-30 seconds)\n")
time_compile <- system.time({
    nimble_model <- nimbleModel(
        code = nimble_code,
        constants = nimble_constants,
        data = nimble_data,
        inits = nimble_inits
    )

    # Compile model
    c_nimble_model <- compileNimble(nimble_model)

    # Configure MCMC
    mcmc_conf <- configureMCMC(nimble_model, print = FALSE)

    # Build MCMC
    mcmc <- buildMCMC(mcmc_conf)

    # Compile MCMC
    c_mcmc <- compileNimble(mcmc, project = nimble_model)
})

cat("Compilation time:", time_compile["elapsed"], "seconds\n\n")

# ============================================================================
# 3. Run NIMBLE
# ============================================================================
cat("--- Running NIMBLE ---\n")
time_nimble <- system.time({
    nimble_samples <- runMCMC(
        c_mcmc,
        niter = 5000,
        nburnin = 1000,
        thin = 5,
        nchains = 1,
        setSeed = 123
    )
})

cat("NIMBLE Time (sampling only):", time_nimble["elapsed"], "seconds\n")
cat(
    "NIMBLE Time (total with compilation):",
    time_nimble["elapsed"] + time_compile["elapsed"],
    "seconds\n\n"
)

# ============================================================================
# 4. Compare Results
# ============================================================================
cat("\n=== COMPARISON ===\n\n")

cat("Timing:\n")
cat(sprintf("  JAGS: %.2f seconds\n", time_jags["elapsed"]))
cat(sprintf("  NIMBLE (compile): %.2f seconds\n", time_compile["elapsed"]))
cat(sprintf("  NIMBLE (sample): %.2f seconds\n", time_nimble["elapsed"]))
cat(sprintf(
    "  NIMBLE (total): %.2f seconds\n",
    time_nimble["elapsed"] + time_compile["elapsed"]
))
cat(sprintf(
    "  Speedup (sampling only): %.2fx\n",
    time_jags["elapsed"] / time_nimble["elapsed"]
))

cat("\n\nParameter Estimates (means):\n")
cat("\nJAGS:\n")
print(fit_jags$summary$statistics[
    c("alphaBM", "betaBM", "betaNL", "lambda_BM", "lambda_NL"),
    "Mean"
])

cat("\nNIMBLE:\n")
nimble_summary <- apply(nimble_samples, 2, mean)
print(nimble_summary[c(
    "alphaBM",
    "betaBM",
    "betaNL",
    "lambda_BM",
    "lambda_NL"
)])

cat("\n\n=== CONCLUSION ===\n")
speedup <- time_jags["elapsed"] / time_nimble["elapsed"]
if (speedup > 1.5) {
    cat(
        "✅ NIMBLE is significantly faster (",
        sprintf("%.2fx", speedup),
        " speedup)\n",
        sep = ""
    )
} else if (speedup > 1.1) {
    cat(
        "⚠️  NIMBLE is moderately faster (",
        sprintf("%.2fx", speedup),
        " speedup)\n",
        sep = ""
    )
} else {
    cat("❌ NIMBLE is not faster for this model\n")
}

cat(
    "\nNote: Compilation overhead is ",
    sprintf("%.2f", time_compile["elapsed"]),
    " seconds\n",
    sep = ""
)
cat("      For repeated runs, only sampling time matters.\n")
cat("      Could cache compiled models for major speedup!\n")
