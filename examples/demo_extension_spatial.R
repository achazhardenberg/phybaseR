# =============================================================================
# Extensibility Demo: Custom Spatial CAR Module for 'because'
# =============================================================================
# This script demonstrates how EASY it is to extend the 'because' package
# with custom covariance structures. You only need to know:
#   1. R S3 methods (basic)
#   2. JAGS model syntax
#
# We create a "spatial_knn" structure for k-nearest-neighbor spatial models.
# No separate package installation required - just define the methods!
# =============================================================================

library(because)

# =============================================================================
# STEP 1: Define a Constructor for Your Custom Structure
# =============================================================================
# This creates an object with a class that because can dispatch on.

spatial_knn <- function(coords, k = 5) {
    # Calculate k-nearest neighbors using base R
    n <- nrow(coords)

    # Distance matrix
    dist_mat <- as.matrix(dist(coords))

    # For each point, find k nearest neighbors
    nn <- t(apply(dist_mat, 1, function(row) {
        order(row)[2:(k + 1)] # Exclude self (index 1 after sorting)
    }))

    # Build adjacency list for JAGS CAR model
    adj <- as.vector(t(nn))
    num <- rep(k, n)
    weights <- rep(1, length(adj)) # Equal weights

    # Also compute a precision matrix for the standard because workflow
    # Using exponential decay distance weighting
    W <- exp(-dist_mat / mean(dist_mat))
    diag(W) <- 0
    D <- diag(rowSums(W))
    Prec <- D - 0.99 * W # CAR-like precision matrix

    structure(
        list(
            coords = coords,
            k = k,
            adj = adj,
            num = num,
            weights = weights,
            n = n,
            Prec = Prec
        ),
        class = c("spatial_knn", "because_structure")
    )
}

# =============================================================================
# STEP 2: Implement the S3 Methods
# =============================================================================

#' Generate JAGS Code for Spatial KNN Structure
#'
#' @param structure A spatial_knn object
#' @param variable_name Name of the error term
#' @param ... Additional arguments (ignored)
#' @return List with setup_code and error_prior
jags_structure_definition.spatial_knn <- function(
    structure,
    variable_name = "u_spatial",
    ...
) {
    # Standard multivariate normal with precision matrix
    setup_code <- c(
        "    # Spatial KNN Structure (distance-weighted precision)"
    )

    error_prior <- paste0(
        "    ",
        variable_name,
        "[1:N] ~ dmnorm(zeros[1:N], tau_spatial_knn * Prec_spatial_knn[1:N, 1:N])"
    )

    return(list(
        setup_code = setup_code,
        error_prior = error_prior
    ))
}

#' Prepare Data for Spatial KNN Structure
#'
#' @param structure A spatial_knn object
#' @param data The model data
#' @param optimize Whether to use optimized formulation
#' @param ... Additional arguments
#' @return List with structure_object and data_list
prepare_structure_data.spatial_knn <- function(
    structure,
    data,
    optimize = TRUE,
    ...
) {
    # Return the precision matrix - this is what because expects
    data_list <- list(
        Prec_spatial_knn = structure$Prec
    )

    return(list(
        structure_object = structure,
        data_list = data_list
    ))
}

# =============================================================================
# STEP 3: Test It!
# =============================================================================

# Simulate spatial data
set.seed(42)
n <- 50 # Smaller for fast demo

# Random coordinates
coords <- data.frame(
    x = runif(n, 0, 10),
    y = runif(n, 0, 10)
)

# True spatial effect (smooth surface)
true_spatial <- sin(coords$x / 2) + cos(coords$y / 2)

# Simulate response with known effect
beta_true <- 0.5
X <- rnorm(n)
Y <- beta_true * X + true_spatial + rnorm(n, sd = 0.3)

# Create data
sim_data <- data.frame(Y = Y, X = X)

# Create our custom spatial structure
spatial_struct <- spatial_knn(as.matrix(coords), k = 5)

cat("\n=== Custom Spatial Structure Created ===\n")
cat("Class:", class(spatial_struct)[1], "\n")
cat("N sites:", spatial_struct$n, "\n")
cat("K neighbors:", spatial_struct$k, "\n")
cat("Precision matrix dimensions:", dim(spatial_struct$Prec), "\n")

# -----------------------------------------------------------------------------
# Fit the model using because!
# -----------------------------------------------------------------------------
cat("\n=== Fitting Model with Custom Spatial Structure ===\n")
cat("Model: Y ~ X with spatial_knn covariance structure\n")
cat("True beta:", beta_true, "\n\n")

# Fit the model
fit <- because(
    equations = list(Y ~ X),
    data = sim_data,
    structure = spatial_struct, # Our custom structure!
    n.iter = 1000,
    n.chains = 2,
    quiet = FALSE
)

# Show results
cat("\n=== Model Results ===\n")
print(summary(fit))

cat("\n=== Demo Complete! ===\n")
cat("This demonstrates that external developers can extend 'because'\n")
cat("by simply defining two S3 methods. No package modification needed!\n")
cat("\nThe estimated beta coefficient should be close to", beta_true, "\n")
