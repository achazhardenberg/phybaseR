## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# library(geosphere) # For distance calculations
# 
# # Data with coordinates
# data <- data.frame(
#     site_id = 1:50,
#     latitude = runif(50, 40, 50),
#     longitude = runif(50, -10, 10),
#     temperature = rnorm(50, 15, 3),
#     species_richness = rpois(50, 20),
#     habitat_quality = rnorm(50, 0, 1)
# )
# 
# # Calculate geographic distances
# coords <- cbind(data$longitude, data$latitude)
# distances <- distm(coords) / 1000 # Convert to km
# 
# # Create spatial covariance matrix (exponential decay)
# spatial_sigma <- 100 # Spatial scale (km)
# spatial_cov <- exp(-distances / spatial_sigma)
# 
# # Format for phybaseR
# rownames(spatial_cov) <- data$site_id
# colnames(spatial_cov) <- data$site_id
# rownames(data) <- data$site_id
# 
# # Fit model with spatial structure
# fit <- phybase_run(
#     equations = list(
#         species_richness ~ temperature + habitat_quality
#     ),
#     data = data,
#     tree = spatial_cov, # Spatial covariance!
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Sharp initial decline, long tail
# spatial_cov <- exp(-distances / sigma)

## ----eval=FALSE---------------------------------------------------------------
# # Smooth, faster decay at distance
# spatial_cov <- exp(-(distances / sigma)^2)

## ----eval=FALSE---------------------------------------------------------------
# library(geoR)
# 
# # More flexible, controls smoothness
# spatial_cov <- matern(distances, phi = sigma, kappa = 1.5)

## ----eval=FALSE---------------------------------------------------------------
# # Exploratory: variogram
# library(gstat)
# library(sp)
# 
# coordinates(data) <- ~ longitude + latitude
# variogram_emp <- variogram(species_richness ~ 1, data)
# plot(variogram_emp)
# 
# # Fit variogram to estimate range
# vgm_fit <- fit.variogram(variogram_emp, vgm("Exp"))
# sigma_estimate <- vgm_fit$range[2]

## ----eval=FALSE---------------------------------------------------------------
# library(kinship2) # Or nadiv, pedigreemm
# 
# # Pedigree structure
# pedigree <- data.frame(
#     id = 1:100,
#     sire = c(NA, NA, rep(1:10, each = 10), ...), # Father IDs
#     dam = c(NA, NA, rep(11:20, each = 10), ...) # Mother IDs
# )
# 
# # Phenotype data
# data <- data.frame(
#     individual = 1:100,
#     body_size = rnorm(100, 50, 10),
#     litter_size = rpois(100, 5),
#     maternal_care = rnorm(100, 0, 1),
#     food_availability = rnorm(100, 0, 1)
# )
# 
# # Calculate additive genetic relatedness matrix (A-matrix)
# ped_obj <- pedigree(
#     id = pedigree$id,
#     dadid = pedigree$sire,
#     momid = pedigree$dam
# )
# A_matrix <- kinship(ped_obj) * 2 # A = 2 × kinship
# 
# # Ensure matching
# rownames(A_matrix) <- pedigree$id
# colnames(A_matrix) <- pedigree$id
# rownames(data) <- data$individual
# 
# # Quantitative genetics SEM!
# fit <- phybase_run(
#     equations = list(
#         body_size ~ food_availability + maternal_care,
#         litter_size ~ body_size + maternal_care
#     ),
#     data = data,
#     tree = A_matrix, # Genetic relatedness
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Variance components from model
# # h² = V_genetic / (V_genetic + V_residual)

## ----eval=FALSE---------------------------------------------------------------
# # Genetic merit of individuals
# # Used in animal/plant breeding

## ----eval=FALSE---------------------------------------------------------------
# # Are traits genetically correlated?
# fit <- phybase_run(
#     equations = list(
#         trait1 ~ environment,
#         trait2 ~ environment
#     ),
#     data = data,
#     tree = A_matrix
# )
# # Correlation in genetic effects

## ----eval=FALSE---------------------------------------------------------------
# # Additive + dominance
# library(nadiv)
# D_matrix <- makeD(pedigree) # Dominance relatedness
# 
# # Would need custom JAGS to include both A and D

## ----eval=FALSE---------------------------------------------------------------
# library(igraph)
# 
# # Social network adjacency matrix
# # 1 = individuals interact, 0 = no interaction
# social_adj <- matrix(0, 30, 30)
# # Fill with observed interactions
# social_adj[1, 2] <- 1
# social_adj[2, 1] <- 1
# # ... etc
# 
# # Convert to covariance
# # Option 1: Direct adjacency as correlation
# social_cov <- social_adj
# diag(social_cov) <- 1
# 
# # Option 2: Normalized Laplacian
# graph <- graph_from_adjacency_matrix(social_adj, mode = "undirected")
# laplacian <- laplacian_matrix(graph, normalized = TRUE)
# social_cov <- solve(laplacian + diag(0.1, nrow(laplacian))) # Regularize
# 
# # Behavioral data
# data <- data.frame(
#     individual = 1:30,
#     dominance_rank = sample(1:30),
#     mating_success = rpois(30, 3),
#     foraging_efficiency = rnorm(30, 100, 20),
#     aggression = rnorm(30, 0, 1)
# )
# 
# rownames(social_cov) <- data$individual
# colnames(social_cov) <- data$individual
# rownames(data) <- data$individual
# 
# # SEM with social structure
# fit <- phybase_run(
#     equations = list(
#         mating_success ~ dominance_rank + aggression,
#         foraging_efficiency ~ dominance_rank
#     ),
#     data = data,
#     tree = social_cov,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Instead of (or in addition to) covariance structure
# g <- graph_from_adjacency_matrix(social_adj)
# 
# data$degree <- degree(g) # Number of connections
# data$betweenness <- betweenness(g) # Bridging position
# data$eigenvector <- eigen_centrality(g)$vector # Influence
# 
# # Use as predictors
# fit <- phybase_run(
#     equations = list(
#         mating_success ~ degree + betweenness
#     ),
#     data = data,
#     tree = social_cov # Still account for network structure
# )

## ----eval=FALSE---------------------------------------------------------------
# # Data ordered by time
# data <- data.frame(
#     time = 1:100,
#     measurement = arima.sim(list(ar = 0.7), n = 100),
#     temperature = rnorm(100, 15, 3),
#     rainfall = rnorm(100, 50, 20)
# )
# 
# # AR(1) correlation matrix
# rho <- 0.7 # Autocorrelation parameter
# n <- nrow(data)
# temporal_cov <- rho^abs(outer(1:n, 1:n, "-"))
# 
# rownames(temporal_cov) <- data$time
# colnames(temporal_cov) <- data$time
# rownames(data) <- data$time
# 
# fit <- phybase_run(
#     equations = list(
#         measurement ~ temperature + rainfall
#     ),
#     data = data,
#     tree = temporal_cov,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # From data
# acf(data$measurement, plot = FALSE)$acf[2] # Lag-1 autocorrelation
# 
# # Or use nlme/lme4 first
# library(nlme)
# temp_fit <- gls(measurement ~ temperature + rainfall,
#     correlation = corAR1(form = ~time),
#     data = data
# )
# rho_est <- coef(temp_fit$modelStruct$corStruct, unconstrained = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # Species functional traits
# traits <- data.frame(
#     species = paste0("sp", 1:50),
#     bill_length = rnorm(50, 20, 5),
#     wing_length = rnorm(50, 100, 20),
#     body_mass = rnorm(50, 50, 15)
# )
# 
# # Euclidean distance in trait space
# trait_matrix <- as.matrix(traits[, -1])
# trait_dist <- dist(trait_matrix)
# trait_dist_matrix <- as.matrix(trait_dist)
# 
# # Convert distance to similarity (Gaussian kernel)
# sigma_trait <- median(trait_dist_matrix)
# functional_sim <- exp(-(trait_dist_matrix / sigma_trait)^2)
# 
# rownames(functional_sim) <- traits$species
# colnames(functional_sim) <- traits$species
# 
# # Use in SEM
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = functional_sim,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Genomic relationship matrix from SNPs
# library(rrBLUP)
# 
# # SNP data (individuals × markers)
# snps <- matrix(sample(0:2, 100 * 1000, replace = TRUE), 100, 1000)
# rownames(snps) <- paste0("ind", 1:100)
# 
# # Calculate genomic relationship matrix
# G_matrix <- A.mat(snps - 1) # Center & scale
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = G_matrix,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Can't pass two matrices to tree argument
# # Use one as tree, other as random effect
# 
# fit <- phybase_run(
#     equations = list(
#         trait ~ environment + (1 | site)
#     ),
#     data = data,
#     tree = phylo_tree, # Phylogenetic structure
#     random = ~ (1 | site), # Spatial grouping
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Phylogeny + Genetic + Spatial
# fit <- phybase_run(
#     equations = list(
#         trait ~ environment + (1 | population) + (1 | site)
#     ),
#     data = data,
#     tree = phylo_tree, # Species phylogeny
#     random = ~ (1 | population) + (1 | site), # Other structures
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# validate_covariance_matrix <- function(mat, data) {
#     # Symmetry
#     if (!isSymmetric(mat)) {
#         stop("Matrix is not symmetric")
#     }
# 
#     # Positive semi-definite
#     eigs <- eigen(mat, only.values = TRUE)$values
#     if (any(eigs < -1e-8)) { # Allow tiny numerical error
#         warning("Matrix may not be positive semi-definite")
#     }
# 
#     # Dimensions
#     if (nrow(mat) != nrow(data)) {
#         stop("Matrix dimensions don't match data")
#     }
# 
#     # Names
#     if (!all(rownames(mat) %in% rownames(data))) {
#         stop("Matrix rownames don't match data rownames")
#     }
# 
#     message("✓ Matrix passed validation")
#     return(TRUE)
# }
# 
# # Use before fitting
# validate_covariance_matrix(my_cov_matrix, data)

## ----eval=FALSE---------------------------------------------------------------
# # 1. Add small value to diagonal
# mat_fixed <- mat + diag(0.01, nrow(mat))
# 
# # 2. Nearest positive semi-definite matrix
# library(Matrix)
# mat_fixed <- nearPD(mat)$mat
# 
# # 3. Check for errors in calculation

## ----eval=FALSE---------------------------------------------------------------
# # Use sparse matrices
# library(Matrix)
# sparse_cov <- as(cov_matrix, "dgCMatrix")
# 
# # Optimize mode
# fit <- phybase_run(..., optimise = TRUE)
# 
# # Reduce matrix size (subsample)

