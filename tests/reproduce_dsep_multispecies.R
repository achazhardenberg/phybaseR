# Reproduce multispecies d-separation power

devtools::load_all(".")

# Simulated multispecies data
set.seed(42)
N_sites <- 50
N_reps <- 3

# Two species influenced by Habitat and Wind
data_list <- list(
    # Y is a list of matrices
    Y = list(
        Sp1 = matrix(rbinom(N_sites * N_reps, 1, 0.4), N_sites, N_reps),
        Sp2 = matrix(rbinom(N_sites * N_reps, 1, 0.6), N_sites, N_reps)
    ),
    Habitat = rnorm(N_sites),
    Wind = rnorm(N_sites)
)

print("--- Testing Multispecies D-sep ---")
# Using the new Auto-Stacking feature
results <- because(
    equations = list(
        Y ~ Habitat,
        p_Y ~ Wind
    ),
    data = data_list,
    family = c(Y = "occupancy"),
    dsep = TRUE,
    n.iter = 100,
    quiet = FALSE
)

print("--- D-sep result Summary ---")
print(summary(results)$dsep)
