# Reproduction: N_SiteID Unknown in D-Sep
library(because)

set.seed(123)
N <- 40 # 40 sites
J <- 3 # 3 visits
# Create SiteID
SiteID <- paste0("S", 1:N)
# Replicate for observations?
# User model is occupancy. Dingo_obs is matrix?
# wooster_data_full likely has SiteID as vector of length N (sites)?
# Or N*J?
# If occupancy family, data usually is N rows (sites). Detection history is matrix.
# Cover_high and protected are site covariates (N).
# SiteID is vector of N.

# Data
protected <- rbinom(N, 1, 0.5)
cover_high <- rnorm(N)
Dingo <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
Fox <- matrix(rbinom(N * 3, 1, 0.5), N, 3)
names(Dingo) <- NULL # Matrix

data <- list(
    Dingo = Dingo,
    Dingo_obs = Dingo, # pre-formatted
    Fox = Fox,
    Fox_obs = Fox,
    protected = protected,
    cover_high = cover_high,
    SiteID = as.factor(SiteID)
)

cat("=== Running Reproduction ===\n")

# Model mimicking user's
# random = ~ (1|SiteID)
# dsep test: cover_high ~ protected
# This test treats cover_high (N length) as response.
# Global random applies (1|SiteID).

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
            data = data,
            random = ~ (1 | SiteID),
            dsep = TRUE,
            n.iter = 100,
            quiet = FALSE
        )
    },
    error = function(e) {
        cat("Catch Error: ", e$message, "\n")
    }
)
