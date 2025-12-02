library(phybaseR)

data("rhino.dat")
data("rhino.tree")

data_list <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    LS = rhino.dat$LS,
    RS = rhino.dat$RS
)

equations <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

# Reproduce the user's workflow
cat("\n=== Running phybase_run (Parallel) ===\n")
fit_rhino <- phybase_run(
    data = data_list,
    tree = rhino.tree,
    equations = equations,
    WAIC = FALSE,
    parallel = TRUE,
    n.cores = 2,
    ic_recompile = FALSE, # Force recompilation in phybase_waic
    n.iter = 1000,
    n.burnin = 500
)

cat("\n=== Running phybase_waic with n.burnin=200 ===\n")
tryCatch(
    {
        # We expect recompilation to use this burn-in
        waic_res <- phybase_waic(fit_rhino, n.burnin = 200)
        print(waic_res)

        cat("\n=== Checking R-hat ===\n")
        sum_stats <- summary(fit_rhino)

        if (is.matrix(sum_stats) || is.data.frame(sum_stats)) {
            if ("Rhat" %in% colnames(sum_stats)) {
                cat("\nR-hat column found in summary statistics!\n")
                print(head(sum_stats[, "Rhat"]))
            } else {
                cat(
                    "\nWARNING: R-hat column NOT found in summary statistics.\n"
                )
            }
        } else {
            print(head(sum_stats$statistics))
        }
    },
    error = function(e) {
        cat("\nERROR CAUGHT:\n")
        print(e)
    }
)
