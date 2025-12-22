# Reproduce Hierarchical Error in D-separation for Occupancy Models

devtools::load_all(".")

# Create multispecies data
set.seed(123)
N_sites <- 20
N_reps <- 3
species <- c("Sp1", "Sp2")

Y_list <- list(
    Sp1 = matrix(rbinom(N_sites * N_reps, 1, 0.4), N_sites, N_reps),
    Sp2 = matrix(rbinom(N_sites * N_reps, 1, 0.6), N_sites, N_reps)
)

data_list <- list(
    Y = Y_list,
    Habitat = rnorm(N_sites),
    Wind = rnorm(N_sites)
)

print("--- Testing multispecies occupancy with dsep=TRUE ---")
tryCatch(
    {
        res <- because(
            equations = list(
                Y ~ Habitat,
                p_Y ~ Wind
            ),
            data = data_list,
            family = c(Y = "occupancy"),
            dsep = TRUE,
            n.iter = 500, # Sufficient for convergence in toy model
            quiet = FALSE,
            verbose = TRUE
        )
        print("SUCCESS")
        print(summary(res))
    },
    error = function(e) {
        print(paste("CAUGHT ERROR:", e$message))
    }
)
