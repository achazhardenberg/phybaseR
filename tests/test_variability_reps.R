library(ape)
library(devtools)
load_all(".")
library(rjags)

# 1. Setup Data (Same as user's structure)
sp <- paste0("s", 1:5)
tree <- rcoal(5, tip.label = sp)

# Repeated measures (3 per species)
set.seed(123)
dat <- data.frame(
    SP = rep(sp, each = 3),
    Y = rnorm(15),
    X = rnorm(15)
)

# 2. Run 'because' with variability="reps"
# This should trigger auto-formatting from long 'dat' to matrices
cat("\n--- Testing variability='reps' Auto-Format ---\n")
tryCatch(
    {
        mod_reps <- because(
            equations = list(Y ~ X),
            data = dat,
            id_col = "SP",
            structure = tree,
            variability = "reps", # Global setting test
            n.adapt = 100,
            n.iter = 100,
            quiet = FALSE
        )
        print("Model with variability='reps' ran successfully!")
        print(summary(mod_reps))

        # Check if model string contains the loop for replicates
        # We expect something like 'Y_obs[i,j] ~ dnorm(Y[i], ...)'
        if (grepl("dnorm\\(.*, .*tau.*\\)", mod_reps$model)) {
            print("Model string contains expected Observation Model structure.")
        }
    },
    error = function(e) {
        print(paste("variability='reps' FAILED:", e$message))
    }
)
