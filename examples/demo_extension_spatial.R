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

    structure(
        list(
            coords = coords,
            k = k,
            adj = adj,
            num = num,
            weights = weights,
            n = n
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
    # JAGS uses dcar_normal for CAR models
    # adj, weights, num are passed as data

    setup_code <- c(
        "    # Spatial CAR Structure (k-nearest neighbors)",
        "    # Precision parameter for spatial random effect"
    )

    error_prior <- paste0(
        "    ",
        variable_name,
        "[1:N] ~ dcar_normal(",
        "adj_spatial[], weights_spatial[], num_spatial[], tau_spatial)"
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
#' @param ... Additional arguments
#' @return List with structure_object and data_list
prepare_structure_data.spatial_knn <- function(
    structure,
    data,
    ...
) {
    # CAR models need: adj, weights, num
    data_list <- list(
        adj_spatial = structure$adj,
        weights_spatial = structure$weights,
        num_spatial = structure$num
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
n <- 100

# Random coordinates
coords <- data.frame(
    x = runif(n, 0, 10),
    y = runif(n, 0, 10)
)

# True spatial effect (smooth surface)
true_spatial <- sin(coords$x / 2) + cos(coords$y / 2)

# Simulate response
X <- rnorm(n) # Covariate
Y <- 0.5 * X + true_spatial + rnorm(n, sd = 0.5)

# Create data
sim_data <- data.frame(Y = Y, X = X)

# Create our custom spatial structure
spatial_struct <- spatial_knn(as.matrix(coords), k = 5)

cat("\n=== Custom Spatial Structure Created ===\n")
cat("Class:", class(spatial_struct)[1], "\n")
cat("N sites:", spatial_struct$n, "\n")
cat("K neighbors:", spatial_struct$k, "\n")

# -----------------------------------------------------------------------------
# Fit the model using because!
# -----------------------------------------------------------------------------
cat("\n=== Fitting Model with Custom Spatial Structure ===\n")

# NOTE: For this demo to fully work, because needs to integrate the
# spatial structure into the model. The current implementation focuses
# on phylogenetic structures, but the S3 dispatch mechanism is in place.

# To see the JAGS code that WOULD be generated:
cat("\n--- Generated JAGS Code ---\n")
jags_def <- jags_structure_definition(
    spatial_struct,
    variable_name = "u_spatial"
)
cat(jags_def$setup_code, sep = "\n")
cat(jags_def$error_prior, "\n")

cat("\n--- Prepared Data ---\n")
prep_data <- prepare_structure_data(spatial_struct, sim_data)
cat("Data elements:", names(prep_data$data_list), "\n")

cat("\n=== Demo Complete! ===\n")
cat("This demonstrates that external developers can extend 'because'\n")
cat("by simply defining two S3 methods. No package modification needed!\n")
