# Reproduce the Hierarchical Error in D-separation (User Case)

devtools::load_all(".")

# Simulated user data
set.seed(123)
N_sites <- 30
N_reps <- 3

# data_list for occupancy: a list containing a matrix Y and vectors for covariates
data_list <- list(
    Y = matrix(rbinom(N_sites * N_reps, 1, 0.5), N_sites, N_reps),
    Habitat = rnorm(N_sites),
    Wind = rnorm(N_sites)
)

print("--- Testing Occupancy + D-sep + List Data ---")
tryCatch(
    {
        results <- because(
            equations = list(
                Y ~ Habitat,
                p_Y ~ Wind
            ),
            data = data_list,
            family = c(Y = "occupancy"),
            dsep = TRUE,
            n.iter = 100, # Fast
            quiet = FALSE,
            verbose = TRUE
        )
        print("--- Basis Set ---")
        print(results$dsep)
        print("SUCCESS")
    },
    error = function(e) {
        print(paste("CAUGHT ERROR:", e$message))
        # traceback()
    }
)
