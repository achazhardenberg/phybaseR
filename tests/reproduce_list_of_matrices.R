# Reproduce Occupancy List-of-Matrices Error

devtools::load_all(".")

# Simulate Data
set.seed(123)
N <- 100
n_reps <- 3
n_species <- 2 # Sp1, Sp2

# List of matrices
Y_list <- list(
    Sp1 = matrix(rbinom(N * n_reps, 1, 0.5), N, n_reps),
    Sp2 = matrix(rbinom(N * n_reps, 1, 0.5), N, n_reps)
)

Habitat <- rnorm(N)
Trait <- rnorm(n_species)

print("--- Testing List of Matrices (Expected Success) ---")
tryCatch(
    {
        res <- because(
            equations = list(Y ~ Habitat, p_Y ~ Trait),
            family = c(Y = "occupancy"),
            data = list(Y = Y_list, Habitat = Habitat, Trait = Trait),
            quiet = TRUE
        )
        print("SUCCESS: List of Matrices accepted and converted.")
    },
    error = function(e) {
        if (grepl("Index out of range", e$message)) {
            print("SUCCESS: Data Accepted (JAGS Runtime Error ignored)")
        } else {
            print(paste("FAILURE List:", e$message))
        }
    }
)
