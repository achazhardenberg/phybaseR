library(phybaseR)

# Simulate bovids data for testing
set.seed(123)
n <- 50

bovids.data <- list(
    BR = rnorm(n),
    BM = rnorm(n),
    S = rnorm(n),
    G = rnorm(n),
    L = rnorm(n)
)

# Create a simple tree
bovidsCut.tree <- ape::rtree(n)

# User's equations with latent variable
bovids.equations_lat <- list(BR ~ Lat, BM ~ Lat, S ~ BM, G ~ BR, L ~ BR)

cat("\n=== Testing MAG Latent Variable Fix ===\n")
tryCatch(
    {
        fit_bovids_CB2_lat <- phybase_run(
            data = bovids.data,
            tree = bovidsCut.tree,
            equations = bovids.equations_lat,
            n.iter = 1000,
            n.burnin = 500,
            latent = c("Lat"),
            WAIC = FALSE,
            parallel = FALSE # Start with sequential for debugging
        )

        cat("\n✓ Model fitted successfully!\n")
        print(summary(fit_bovids_CB2_lat))
    },
    error = function(e) {
        cat("\n✗ ERROR:\n")
        print(e)
    }
)
