library(becauseR)

# Check package version
cat("Package namespace functions:\n")
cat("equations_to_dag signature:\n")
print(args(becauseR:::equations_to_dag))

# Simple test
equations <- list(
    Y ~ X + (1 | GroupA),
    X ~ L1,
    Z ~ L1
)
latent <- "L1"
random_terms <- list(
    list(response = "Y", group = "GroupA", type = "1")
)

cat("\n\n=== Testing equations_to_dag directly ===\n")
grouping_vars <- unique(sapply(random_terms, function(x) x$group))
cat("Grouping vars to exclude:", paste(grouping_vars, collapse = ", "), "\n")

dag <- becauseR:::equations_to_dag(equations, exclude_vars = grouping_vars)
cat("\nDAG variables:\n")
print(rownames(dag))
cat("\nDAG matrix:\n")
print(dag)
