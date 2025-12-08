## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# library(phybaseR)
# library(ape)
# 
# # Load or create phylogeny
# tree <- read.tree("my_phylogeny.nwk")
# 
# # Species data (row names = species names matching tree)
# data <- data.frame(
#     species = tree$tip.label,
#     body_size = rnorm(length(tree$tip.label), 100, 20),
#     range_size = rnorm(length(tree$tip.label), 500, 100),
#     population = rnorm(length(tree$tip.label), 1000, 200),
#     row.names = "species"
# )
# 
# # Phylogenetic SEM
# fit <- phybase_run(
#     equations = list(
#         range_size ~ body_size,
#         population ~ range_size + body_size
#     ),
#     data = data,
#     tree = tree, # Include phylogeny!
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Option 1: Row names
# rownames(data) <- data$species
# data$species <- NULL # Remove column
# 
# # Option 2: Specify id_col
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree,
#     id_col = "species" # Column name with species IDs
# )

## ----eval=FALSE---------------------------------------------------------------
# # Prune tree to match data
# library(ape)
# tree_pruned <- keep.tip(tree, rownames(data))
# 
# # Or use geiger for matching
# library(geiger)
# matched <- treedata(tree, data, sort = TRUE)
# tree_clean <- matched$phy
# data_clean <- matched$data

## ----eval=FALSE---------------------------------------------------------------
# is.ultrametric(tree) # Should be TRUE
# is.binary(tree) # Should be TRUE
# all(tree$edge.length > 0) # Should be TRUE

## ----eval=FALSE---------------------------------------------------------------
# # Load multiple trees (e.g., from BEAST posterior)
# trees <- read.nexus("posterior_trees.nex")
# 
# # Sample subset for computational efficiency
# sampled_trees <- sample(trees, size = 100)
# 
# # Fit model
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = sampled_trees, # List of trees!
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree,
#     optimise = FALSE # Full MCMC
# )

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree,
#     optimise = TRUE # Faster!
# )

## ----eval=FALSE---------------------------------------------------------------
# # Model without phylogeny
# fit_nophy <- phybase_run(
#     equations = equations,
#     data = data,
#     # No tree!
#     n.chains = 3,
#     n.iter = 5000
# )
# 
# # Model with phylogeny
# fit_phy <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree, # With tree
#     n.chains = 3,
#     n.iter = 5000
# )
# 
# # Compare estimates
# summary(fit_nophy)
# summary(fit_phy)

## ----eval=FALSE---------------------------------------------------------------
# library(geiger)
# 
# # Reduce phylogenetic signal
# tree_lambda <- rescale(tree, "lambda", 0.5) # λ = 0.5
# 
# # Use transformed tree
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree_lambda
# )

## ----eval=FALSE---------------------------------------------------------------
# library(phytools)
# 
# # Pagel's lambda for each variable
# phylosig(tree, data$body_size, method = "lambda")
# phylosig(tree, data$range_size, method = "lambda")

## ----eval=FALSE---------------------------------------------------------------
# library(sp)
# library(geosphere)
# 
# # Your data with coordinates
# data <- data.frame(
#     site_id = 1:50,
#     latitude = runif(50, 40, 50),
#     longitude = runif(50, -10, 10),
#     temperature = rnorm(50, 15, 3),
#     species_richness = rpois(50, 20)
# )
# 
# # Create spatial distance matrix
# coords <- cbind(data$longitude, data$latitude)
# distances <- distm(coords) / 1000 # km
# 
# # Convert to covariance (exponential decay)
# spatial_sigma <- 100 # Spatial scale parameter (km)
# spatial_cov <- exp(-distances / spatial_sigma)
# 
# # Ensure proper format (row names match data)
# rownames(spatial_cov) <- data$site_id
# colnames(spatial_cov) <- data$site_id
# rownames(data) <- data$site_id
# 
# # Use as "tree" argument!
# fit <- phybase_run(
#     equations = list(
#         species_richness ~ temperature
#     ),
#     data = data,
#     tree = spatial_cov, # Spatial covariance matrix
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# library(kinship2) # or nadiv
# 
# # Pedigree data
# pedigree <- data.frame(
#     id = 1:100,
#     sire = c(NA, NA, 1, 1, 2, 2, 3, ...), # Father IDs
#     dam = c(NA, NA, 2, 2, 1, 1, 4, ...) # Mother IDs
# )
# 
# # Phenotype data
# data <- data.frame(
#     individual = 1:100,
#     body_size = rnorm(100, 50, 10),
#     litter_size = rpois(100, 5),
#     food_availability = rnorm(100, 0, 1)
# )
# 
# # Calculate additive genetic relatedness matrix (A-matrix)
# ped <- pedigree(
#     id = pedigree$id,
#     dadid = pedigree$sire,
#     momid = pedigree$dam
# )
# A_matrix <- kinship(ped) * 2 # A = 2 × kinship
# 
# # Ensure matching
# rownames(A_matrix) <- pedigree$id
# colnames(A_matrix) <- pedigree$id
# rownames(data) <- data$individual
# 
# # Quantitative genetics SEM!
# fit <- phybase_run(
#     equations = list(
#         body_size ~ food_availability,
#         litter_size ~ body_size
#     ),
#     data = data,
#     tree = A_matrix, # Genetic relatedness
#     dsep = TRUE,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# library(igraph)
# 
# # Social network adjacency matrix
# # (who interacts with whom)
# social_network <- matrix(0, 30, 30)
# # ... fill with interaction data ...
# 
# # Convert to correlation matrix (or use directly)
# social_cov <- cov2cor(social_network %*% t(social_network))
# 
# # Behavioral data
# data <- data.frame(
#     individual = 1:30,
#     dominance_rank = 1:30,
#     mating_success = rpois(30, 3),
#     foraging_efficiency = rnorm(30, 100, 20)
# )
# 
# rownames(social_cov) <- data$individual
# colnames(social_cov) <- data$individual
# rownames(data) <- data$individual
# 
# # Account for social structure
# fit <- phybase_run(
#     equations = list(
#         mating_success ~ dominance_rank,
#         foraging_efficiency ~ dominance_rank
#     ),
#     data = data,
#     tree = social_cov,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Functional similarity (trait-based)
# # Morphological distances
# # Ecological niche overlap
# # Genomic similarity (SNPs)
# # Temporal autocorrelation
# 
# # Just ensure:
# # 1. Symmetric matrix
# # 2. Positive semi-definite
# # 3. Row/column names match data
# # 4. Diagonal typically = 1 (correlation) or max (covariance)
# 
# custom_cov <- your_similarity_function(data)
# rownames(custom_cov) <- data$id
# colnames(custom_cov) <- data$id
# 
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = custom_cov,
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Phylogeny for species + random effects for space
# fit <- phybase_run(
#     equations = list(
#         trait ~ environment + (1 | site)
#     ),
#     data = data,
#     tree = phylo_tree, # Phylogenetic correlation
#     random = ~ (1 | site), # Spatial grouping
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Phylogenetic + spatial
# # Genetic + environmental
# # Multiple levels of structure

## ----eval=FALSE---------------------------------------------------------------
# # Is it symmetric?
# isSymmetric(my_matrix)
# 
# # Is it positive semi-definite?
# all(eigen(my_matrix)$values >= -1e-8) # Allow tiny negative for numerical error
# 
# # Are dimensions correct?
# nrow(my_matrix) == nrow(data)
# 
# # Do names match?
# all(rownames(my_matrix) %in% rownames(data))

## ----eval=FALSE---------------------------------------------------------------
# # Account for phylogeny AND repeated measures
# fit <- phybase_run(
#     equations = list(
#         trait ~ environment + (1 | population)
#     ),
#     data = data,
#     tree = tree,
#     random = ~ (1 | species), # Additional random effect
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Species-level phylogeny with individual-level data
# fit <- phybase_run(
#     data = list(
#         species = species_data,
#         individual = individual_data
#     ),
#     levels = list(
#         species = c("range_size", "habitat"),
#         individual = c("body_mass", "age")
#     ),
#     hierarchy = "species > individual",
#     link_vars = "species_id",
#     tree = tree, # Applied at species level
#     equations = equations
# )

## ----eval=FALSE---------------------------------------------------------------
# # What's in tree?
# tree$tip.label
# 
# # What's in data?
# rownames(data) # or data$species
# 
# # Find mismatches
# setdiff(tree$tip.label, rownames(data)) # In tree but not data
# setdiff(rownames(data), tree$tip.label) # In data but not tree

## ----eval=FALSE---------------------------------------------------------------
# library(ape)
# tree_ultra <- chronos(tree) # Date the tree

## ----eval=FALSE---------------------------------------------------------------
# tree_binary <- multi2di(tree) # Randomly resolve

## ----eval=FALSE---------------------------------------------------------------
# fit <- phybase_run(
#     equations = equations,
#     data = data,
#     tree = tree,
#     optimise = TRUE, # Much faster for large trees
#     n.chains = 3,
#     n.iter = 5000
# )

## ----eval=FALSE---------------------------------------------------------------
# plot(tree)

