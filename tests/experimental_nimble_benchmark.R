# ==============================================================================
# NIMBLE vs JAGS Benchmark: Full sem8 Model
# ==============================================================================
# Experimental comparison to assess whether NIMBLE provides additional speedup
# beyond the optimized JAGS implementation (random effects formulation).
# ==============================================================================

library(phybaseR)
library(nimble)
library(rjags)
library(coda)
library(ape)

# Load data
data("rhino.dat")
data("rhino.tree")

# Tree preprocessing (matching sem8.R)
rhino.tree$edge.length <- rhino.tree$edge.length /
    max(branching.times(rhino.tree))
rhino.vcv <- vcv.phylo(rhino.tree)
ID <- diag(100)
N <- 100
Prec_phylo_fixed <- solve(rhino.vcv)

# MCMC settings
nc <- 3
ni <- 12000
nb <- 2000
nt <- 10
set.seed(12345)

cat("\n=== NIMBLE vs JAGS: sem8 Model Benchmark ===\n\n")

# ==============================================================================
# 1. JAGS (Random Effects - Optimized)
# ==============================================================================

cat("--- 1. JAGS (Random Effects) ---\n")

jags_re_model <- "
model {
  # Priors
  alphaLS ~ dnorm(0, 1.0E-06)
  alphaNL ~ dnorm(0, 1.0E-06)
  alphaDD ~ dnorm(0, 1.0E-06)
  betaBM ~ dnorm(0, 1.0E-06)
  betaBM2 ~ dnorm(0, 1.0E-06)
  betaRS ~ dnorm(0, 1.0E-06)
  betaNL ~ dnorm(0, 1.0E-06)

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

  u_std_LS[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_NL[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_DD[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_BM[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])
  u_std_RS[1:N] ~ dmnorm(zeros[1:N], Prec_phylo_fixed[1:N, 1:N])

  for(i in 1:N) {
    u_LS[i] <- u_std_LS[i] / sqrt(tau_u_LS)
    u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
    u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)
    u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
    u_RS[i] <- u_std_RS[i] / sqrt(tau_u_RS)

    muLS[i] <- alphaLS + betaBM*BM[i]
    muNL[i] <- alphaNL + betaBM2*BM[i] + betaRS*RS[i]
    muDD[i] <- alphaDD + betaNL*NL[i]
    muBM[i] <- 0
    muRS[i] <- 0

    LS[i] ~ dnorm(muLS[i] + u_LS[i], tau_e_LS)
    NL[i] ~ dnorm(muNL[i] + u_NL[i], tau_e_NL)
    DD[i] ~ dnorm(muDD[i] + u_DD[i], tau_e_DD)
    BM[i] ~ dnorm(muBM[i] + u_BM[i], tau_e_BM)
    RS[i] ~ dnorm(muRS[i] + u_RS[i], tau_e_RS)
  }

  lambdaLS <- (1/tau_u_LS) / ((1/tau_u_LS) + (1/tau_e_LS))
  lambdaNL <- (1/tau_u_NL) / ((1/tau_u_NL) + (1/tau_e_NL))
  lambdaDD <- (1/tau_u_DD) / ((1/tau_u_DD) + (1/tau_e_DD))
}
"

jags_data <- list(
    N = N,
    LS = rhino.dat$LS,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    BM = rhino.dat$BM,
    RS = rhino.dat$RS,
    Prec_phylo_fixed = Prec_phylo_fixed,
    zeros = rep(0, N)
)

cat("Building and sampling...\n")
time_jags_re <- system.time({
    model_jags <- jags.model(
        textConnection(jags_re_model),
        data = jags_data,
        n.chains = nc,
        quiet = TRUE
    )
    update(model_jags, nb)
    samples_jags <- coda.samples(
        model_jags,
        variable.names = c(
            "betaBM",
            "betaBM2",
            "betaRS",
            "betaNL",
            "lambdaLS",
            "lambdaNL",
            "lambdaDD"
        ),
        n.iter = ni,
        thin = nt
    )
})

cat(sprintf("✓ JAGS completed in %.2f seconds\n", time_jags_re["elapsed"]))

# ==============================================================================
# 2. NIMBLE (Random Effects)
# ==============================================================================

cat("\n--- 2. NIMBLE (Random Effects) ---\n")

nimble_code <- nimbleCode({
    # Priors
    alphaLS ~ dnorm(0, sd = 1000)
    alphaNL ~ dnorm(0, sd = 1000)
    alphaDD ~ dnorm(0, sd = 1000)
    betaBM ~ dnorm(0, sd = 1000)
    betaBM2 ~ dnorm(0, sd = 1000)
    betaRS ~ dnorm(0, sd = 1000)
    betaNL ~ dnorm(0, sd = 1000)

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

    u_std_LS[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_NL[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_DD[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_BM[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])
    u_std_RS[1:N] ~ dmnorm(zeros[1:N], prec = Prec_phylo_fixed[1:N, 1:N])

    for (i in 1:N) {
        u_LS[i] <- u_std_LS[i] / sqrt(tau_u_LS)
        u_NL[i] <- u_std_NL[i] / sqrt(tau_u_NL)
        u_DD[i] <- u_std_DD[i] / sqrt(tau_u_DD)
        u_BM[i] <- u_std_BM[i] / sqrt(tau_u_BM)
        u_RS[i] <- u_std_RS[i] / sqrt(tau_u_RS)

        muLS[i] <- alphaLS + betaBM * BM[i]
        muNL[i] <- alphaNL + betaBM2 * BM[i] + betaRS * RS[i]
        muDD[i] <- alphaDD + betaNL * NL[i]
        muBM[i] <- 0
        muRS[i] <- 0

        LS[i] ~ dnorm(muLS[i] + u_LS[i], tau = tau_e_LS)
        NL[i] ~ dnorm(muNL[i] + u_NL[i], tau = tau_e_NL)
        DD[i] ~ dnorm(muDD[i] + u_DD[i], tau = tau_e_DD)
        BM[i] ~ dnorm(muBM[i] + u_BM[i], tau = tau_e_BM)
        RS[i] ~ dnorm(muRS[i] + u_RS[i], tau = tau_e_RS)
    }

    lambdaLS <- (1 / tau_u_LS) / ((1 / tau_u_LS) + (1 / tau_e_LS))
    lambdaNL <- (1 / tau_u_NL) / ((1 / tau_u_NL) + (1 / tau_e_NL))
    lambdaDD <- (1 / tau_u_DD) / ((1 / tau_u_DD) + (1 / tau_e_DD))
})

nimble_data <- list(
    LS = rhino.dat$LS,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    BM = rhino.dat$BM,
    RS = rhino.dat$RS
)

nimble_constants <- list(
    N = N,
    Prec_phylo_fixed = Prec_phylo_fixed,
    zeros = rep(0, N)
)

nimble_inits <- list(
    alphaLS = 0,
    alphaNL = 0,
    alphaDD = 0,
    betaBM = 0,
    betaBM2 = 0,
    betaRS = 0,
    betaNL = 0,
    tau_u_LS = 1,
    tau_e_LS = 1,
    tau_u_NL = 1,
    tau_e_NL = 1,
    tau_u_DD = 1,
    tau_e_DD = 1,
    tau_u_BM = 1,
    tau_e_BM = 1,
    tau_u_RS = 1,
    tau_e_RS = 1,
    u_std_LS = rep(0, N),
    u_std_NL = rep(0, N),
    u_std_DD = rep(0, N),
    u_std_BM = rep(0, N),
    u_std_RS = rep(0, N)
)

cat("Building NIMBLE model (this may take ~30s)...\n")
time_compile <- system.time({
    nimble_model <- nimbleModel(
        code = nimble_code,
        constants = nimble_constants,
        data = nimble_data,
        inits = nimble_inits
    )
    c_nimble <- compileNimble(nimble_model)

    mcmc_conf <- configureMCMC(
        nimble_model,
        monitors = c(
            "betaBM",
            "betaBM2",
            "betaRS",
            "betaNL",
            "lambdaLS",
            "lambdaNL",
            "lambdaDD"
        ),
        print = FALSE
    )
    mcmc <- buildMCMC(mcmc_conf)
    c_mcmc <- compileNimble(mcmc, project = nimble_model)
})

cat(sprintf("✓ Compiled in %.2f seconds\n", time_compile["elapsed"]))

cat("Sampling...\n")
time_nimble_sample <- system.time({
    nimble_samples <- runMCMC(
        c_mcmc,
        niter = 60000,
        nburnin = 10000,
        thin = nt,
        nchains = nc,
        setSeed = c(12345, 12346, 12347) # One seed per chain
    )
})

cat(sprintf("✓ Sampled in %.2f seconds\n", time_nimble_sample["elapsed"]))

# ==============================================================================
# 3. Comparison
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RESULTS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Timing:\n")
cat(sprintf("  JAGS (Random Effects):  %.2f s\n", time_jags_re["elapsed"]))
cat(sprintf("  NIMBLE Compile:         %.2f s\n", time_compile["elapsed"]))
cat(sprintf(
    "  NIMBLE Sample:          %.2f s\n",
    time_nimble_sample["elapsed"]
))
cat(sprintf(
    "  NIMBLE Total:           %.2f s\n",
    time_compile["elapsed"] + time_nimble_sample["elapsed"]
))
cat("\n")

speedup_vs_compile <- time_jags_re["elapsed"] /
    (time_compile["elapsed"] + time_nimble_sample["elapsed"])
speedup_vs_sample <- time_jags_re["elapsed"] / time_nimble_sample["elapsed"]

cat(sprintf("Speedup (vs total):     %.2fx\n", speedup_vs_compile))
cat(sprintf("Speedup (vs sampling):  %.2fx\n", speedup_vs_sample))

cat("\n", rep("=", 70), "\n", sep = "")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep = "")

if (speedup_vs_compile > 2) {
    cat(
        "✓ NIMBLE provides substantial additional speedup beyond JAGS optimization.\n"
    )
    cat(
        "  Migration to NIMBLE would be worthwhile for performance-critical users.\n"
    )
} else if (speedup_vs_compile > 1.5) {
    cat("~ NIMBLE provides moderate additional speedup.\n")
    cat("  Consider offering NIMBLE as an optional backend.\n")
} else {
    cat("✗ NIMBLE does not provide meaningful speedup beyond optimized JAGS.\n")
    cat("  Stick with JAGS optimization (no C++ compiler dependency).\n")
}

cat("\nNote: For repeated runs of the same model, NIMBLE's compilation cost\n")
cat("      is amortized, making the sampling-only speedup more relevant.\n")
