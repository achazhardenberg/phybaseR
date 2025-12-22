# Minimal Reproduction of User's D-Sep Call
if (requireNamespace("devtools", quietly = TRUE)) {
    # devtools::load_all(".")
    # We want to test INSTALLED package to match user environment
    library(because)
}

set.seed(123)
N <- 100
# Simulate simple occupancy data
protected <- rbinom(N, 1, 0.5)
cover_high <- rnorm(N)
Dingo <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
Fox <- matrix(rbinom(N * 3, 1, 0.5), N, 3)

# Data frame format as user likely has
data_sim <- list(
    Dingo = Dingo,
    Fox = Fox,
    protected = protected,
    cover_high = cover_high
)

cat("=== Running Debug Reproduction ===\n")

tryCatch(
    {
        m <- because(
            equations = list(
                Dingo ~ protected + cover_high,
                p_Dingo ~ 1,
                Fox ~ protected + cover_high,
                p_Fox ~ cover_high
            ),
            family = c(Dingo = "occupancy", Fox = "occupancy"),
            data = data_sim,
            dsep = TRUE,
            n.iter = 0, # Don't fit, just generate tests
            quiet = FALSE
        )
    },
    error = function(e) {
        cat("Error during run: ", e$message, "\n")
    }
)
