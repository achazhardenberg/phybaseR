# =============================================================================
# Extensibility Demo: Easy Custom Extensions with Helper Functions
# =============================================================================
# This demo shows how EASY it is to extend 'because' using the helper functions.
# You only need to provide the mathematical core - the helpers handle everything else!
# =============================================================================

library(because)

# =============================================================================
# PART 1: Custom Structure using because_structure() helper
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 1: Creating a Custom Spatial Structure\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Create a spatial distance-decay structure with ONE function call!
# Just provide the precision matrix function - that's it!
spatial_decay <- because_structure(
    name = "spatial_decay",
    precision_fn = function(coords, decay_rate = 0.1) {
        # Calculate distance matrix
        dist_mat <- as.matrix(dist(coords))

        # Weight by exponential decay
        W <- exp(-decay_rate * dist_mat)
        diag(W) <- 0

        # Create precision matrix (D - rho*W style)
        D <- diag(rowSums(W))
        D - 0.99 * W # Return precision matrix
    },
    description = "Exponential distance decay spatial structure"
)

# =============================================================================
# PART 2: Simulate Data and Fit Model
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 2: Fitting a Model with the Custom Structure\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Simulate spatial data
set.seed(123)
n <- 40 # Smaller for fast demo

# Random spatial coordinates
coords <- matrix(runif(n * 2), ncol = 2)

# True spatial effect (smooth surface)
true_spatial <- sin(coords[, 1] * 3) + cos(coords[, 2] * 3)

# Simulate response
beta_true <- 0.7
X <- rnorm(n)
Y <- beta_true * X + true_spatial + rnorm(n, sd = 0.3)

# Create data frame
sim_data <- data.frame(Y = Y, X = X)

# Create structure using our helper-generated constructor
my_struct <- spatial_decay(coords, decay_rate = 0.2)

cat("Structure created!\n")
cat("  Class:", paste(class(my_struct), collapse = ", "), "\n")
cat("  N:", my_struct$n, "\n")
cat(
    "  Precision matrix:",
    dim(my_struct$Prec)[1],
    "x",
    dim(my_struct$Prec)[2],
    "\n\n"
)

# Fit the model
cat("Fitting model: Y ~ X with spatial_decay structure\n")
cat("True beta:", beta_true, "\n\n")

fit <- because(
    equations = list(Y ~ X),
    data = sim_data,
    structure = my_struct,
    n.iter = 1000,
    n.chains = 2,
    quiet = FALSE
)

# =============================================================================
# PART 3: Results
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 3: Results\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

print(summary(fit))

cat("\n", paste(rep("-", 50), collapse = ""), "\n")
cat("SUMMARY\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
cat("True beta:", beta_true, "\n")

# Extract from MCMC samples
samples <- as.matrix(fit$samples)
beta_est <- mean(samples[, "beta_Y_X"])
sigma_est <- mean(samples[, "sigma_Y_spatial_decay"])

cat("Estimated beta:", round(beta_est, 3), "\n")
cat("Spatial structure variance:", round(sigma_est, 3), "\n")

# =============================================================================
# CONCLUSION
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("WHAT WE LEARNED\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("With 'because', extending to custom structures is as simple as:\n\n")

cat("1. Define your precision matrix function:\n\n")
cat("   my_structure <- because_structure(\n")
cat("       name = 'my_custom',\n")
cat("       precision_fn = function(data) {\n")
cat("           # Your math here - return NxN precision matrix\n")
cat("       }\n")
cat("   )\n\n")

cat("2. Use it in any because() model:\n\n")
cat("   struct <- my_structure(your_data)\n")
cat("   because(Y ~ X, data = data, structure = struct)\n\n")

cat("That's it! No package creation, no S3 method writing,\n")
cat("no JAGS expertise needed - just provide the math!\n\n")
