# Reproduce DAG Edge Missing Bug
# User reports that psi_Dingo is tested for independence with its own covariates

library(because)

# Create dummy data matching user structure
set.seed(123)
N <- 100
SiteID <- rep(1:20, each = 5)
cover_high <- rnorm(N)
protected <- rbinom(N, 1, 0.5)

# Simulate Dingo and Fox (occupancy)
# Just dummy data, structure matters for d-sep
Dingo <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
Fox <- matrix(rbinom(N * 3, 1, 0.5), N, 3)

wooster_data_full <- list(
    Dingo = Dingo,
    Fox = Fox,
    cover_high = cover_high,
    protected = protected,
    SiteID = as.factor(SiteID)
)

print("Running d-separation test to inspect basis set...")

tryCatch(
    {
        m4a_cover_rand_ds <- because(
            equations = list(
                psi_Dingo ~ protected + cover_high, # Explicit psi_
                p_Dingo ~ 1,
                z_Fox ~ protected + cover_high, # Explicit z_
                p_Fox ~ psi_Dingo + cover_high
            ),
            family = c(Dingo = "occupancy", Fox = "occupancy"),
            data = wooster_data_full,
            random = ~ (1 | SiteID),
            quiet = FALSE,
            dsep = TRUE
        )

        # Check if Test 1 (cover_high ~ protected) used binomial
        # (Testing first element of dsep_results)
        # Actually, we can just look at the output messages if we add a print in because.R
    },
    error = function(e) {
        print(e)
    }
)
