library(phybaseR)
library(ape)

data("rhino.dat")
data("rhino.tree")

data_list <- list(
    BM = rhino.dat$BM,
    NL = rhino.dat$NL,
    DD = rhino.dat$DD,
    LS = rhino.dat$LS,
    RS = rhino.dat$RS
)
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

# Run with dsep=TRUE
tryCatch(
    {
        fit_rhino_dsep <- phybase_run(
            data = data_list,
            tree = rhino.tree,
            equations = equations,
            dsep = TRUE,
            n.iter = 100, # Short run for speed
            n.burnin = 50
        )
        print(summary(fit_rhino_dsep))
    },
    error = function(e) {
        print(e)
    }
)
