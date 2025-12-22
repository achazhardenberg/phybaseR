# Verification: Granular Reuse Safety
devtools::load_all(".")

# Data
set.seed(42)
N <- 100
# Occupancy variables must be matrices (N x J)
Dingo <- matrix(rbinom(N * 3, 1, 0.5), nrow = N)
Fox <- matrix(rbinom(N * 3, 1, 0.5), nrow = N)
protected <- rbinom(N, 1, 0.5)
cover_high <- rnorm(N)
# Use list to allow matrices
data_list <- list(
    Dingo = Dingo,
    Fox = Fox,
    protected = protected,
    cover_high = cover_high
)

# Model 3 (Base)
eq3 <- list(
    Dingo ~ 1,
    p_Dingo ~ psi_Fox + protected,
    Fox ~ 1,
    p_Fox ~ psi_Dingo + cover_high
)

message("--- Running M3 (Base) ---")
m3 <- because(
    eq3,
    data = data_list,
    family = c(Dingo = "occupancy", Fox = "occupancy"),
    dsep = TRUE,
    quiet = FALSE,
    n.iter = 100,
    n.adapt = 100
)

# Model 4 (Refined)
# p_Fox changed
eq4 <- list(
    Dingo ~ 1,
    p_Dingo ~ psi_Fox + protected,
    Fox ~ 1,
    p_Fox ~ psi_Dingo + protected + cover_high # CHANGED
)

message("--- Running M4 (Reuse M3) ---")
m4 <- because(
    eq4,
    data = data_list,
    family = c(Dingo = "occupancy", Fox = "occupancy"),
    dsep = TRUE,
    quiet = FALSE,
    n.iter = 100,
    n.adapt = 100,
    reuse_models = list(m3)
)

test_str <- "cover_high ~ protected"
find_test <- function(m) {
    for (i in seq_along(m$dsep_tests)) {
        if (grepl("cover_high ~ protected", deparse(m$dsep_tests[[i]]))) {
            return(i)
        }
    }
    return(NULL)
}

idx3 <- find_test(m3)
idx4 <- find_test(m4)

if (!is.null(idx3) && !is.null(idx4)) {
    s3 <- m3$dsep_results[[idx3]]$samples
    s4 <- m4$dsep_results[[idx4]]$samples

    if (identical(s3, s4)) {
        message("SUCCESS: Test reused despite global equation change!")
    } else {
        message("FAILURE: Test re-run (Global check prevented reuse).")
    }
} else {
    message("Test not found in one of the models.")
}
