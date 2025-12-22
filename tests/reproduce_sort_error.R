# Reproduce sort.int error in because()

devtools::load_all(".")

# Simulate Data
set.seed(123)
N <- 100
n_reps <- 3
Predator <- matrix(rbinom(N * n_reps, 1, 0.5), N, n_reps)
Prey <- matrix(rbinom(N * n_reps, 1, 0.5), N, n_reps)

# For occupancy, we technically need N_reps vector usually, or because() handles it?
# because() expects a list if not hierarchical?
# The user passed `data = data_list`.
# Let's assume standard structure:
data_list <- list(
    Predator = Predator,
    Prey = Prey,
    N = N
)

print("--- Data List Prepared ---")
str(data_list)

print("--- Running because() ---")
tryCatch(
    {
        res <- because(
            equations = list(
                Predator ~ 1,
                p_Predator ~ 1,
                Prey ~ 1,
                p_Prey ~ Predator # Should link to z_Predator
            ),
            family = c(Predator = "occupancy", Prey = "occupancy"),
            data = data_list,
            quiet = TRUE
        )
        print("SUCCESS")
    },
    error = function(e) {
        print(paste("FAILURE:", e$message))
        # trace output if possible?
        traceback()
    }
)
