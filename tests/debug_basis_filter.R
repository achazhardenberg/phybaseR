library(becauseR)

# Simple test
equations <- list(
    Y ~ X,
    X ~ L1,
    Z ~ L1
)
latent <- "L1"

# Call with random effects
result <- because_dsep(
    equations,
    latent = latent,
    random_terms = list(
        list(response = "Y", group = "GroupA", type = "1"),
        list(response = "X", group = "GroupA", type = "1")
    ),
    quiet = FALSE
)

cat("\n\n=== RESULT STRUCTURE ===\n")
str(result$tests, max.level = 1)
