# Test Stacked DataFrame Input (unmarked-style)

devtools::load_all(".")

# Simulate Data
set.seed(123)
N_sites <- 50
N_reps <- 3
N_species <- 5

# Create Stacked Data
Y_mat <- matrix(
    rbinom(N_sites * N_species * N_reps, 1, 0.5),
    N_sites * N_species,
    N_reps
)
Habitat <- rnorm(N_sites * N_species)
Species_Trait <- rnorm(N_sites * N_species)
SpeciesID <- factor(rep(1:N_species, each = N_sites))

# Dataframe with Matrix Column
df <- data.frame(
    Habitat = Habitat,
    Trait = Species_Trait,
    SpeciesID = SpeciesID
)
df$Y <- Y_mat

print("--- Testing Stacked DataFrame (unmarked-style) ---")
tryCatch(
    {
        res <- because(
            equations = list(
                Y ~ Habitat + Trait,
                p_Y ~ Trait
            ),
            family = c(Y = "occupancy"),
            random = ~ (1 | SpeciesID),
            data = df,
            n.iter = 1000,
            n.burnin = 500,
            quiet = TRUE
        )
        print("SUCCESS: Stacked DataFrame accepted and model run.")
        summ <- summary(res)
        if (any(grepl("SpeciesID", rownames(summ$statistics)))) {
            print(
                "Verified: Random effect parameters for SpeciesID found in output."
            )
        }
    },
    error = function(e) {
        print(paste("FAILURE Stacked DF:", e$message))
    }
)
