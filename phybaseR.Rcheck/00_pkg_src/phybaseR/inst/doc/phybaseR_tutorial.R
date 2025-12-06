## ----setup--------------------------------------------------------------------
library(phybaseR)

## ----load_ape-----------------------------------------------------------------
library(ape)

## ----basic_sim----------------------------------------------------------------
set.seed(123)
N <- 50
tree <- rtree(N)

# Simulate traits and create a data.frame
my_data <- data.frame(
    species = tree$tip.label,
    X = rTraitCont(tree, model = "BM", sigma = 1),
    Y = 0.5 + 0.8 * rTraitCont(tree, model = "BM", sigma = 1) +
        rTraitCont(tree, model = "BM", sigma = 0.5)
)

## ----basic_run, eval=FALSE----------------------------------------------------
# # Define equations
# equations <- list(Y ~ X)
# 
# # Run model - pass the data.frame directly!
# fit <- phybase_run(
#     data = my_data,
#     id_col = "species", # Column with species names
#     tree = tree,
#     equations = equations,
#     n.iter = 5000,
#     n.burnin = 1000,
#     n.chains = 3
# )
# 
# # View summary
# summary(fit)

## ----mediation_sim------------------------------------------------------------
# Simulate mediator M
M <- 0.3 + 0.6 * my_data$X + rTraitCont(tree, model = "BM", sigma = 0.5)
# Simulate Y affected by both X and M
Y_med <- 0.5 + 0.4 * my_data$X + 0.5 * M + rTraitCont(tree, model = "BM", sigma = 0.5)

# Create data.frame with all variables
med_data <- data.frame(
    species = tree$tip.label,
    X = my_data$X,
    M = M,
    Y = Y_med
)

## ----mediation_run, eval=FALSE------------------------------------------------
# # Define DAG equations
# equations_med <- list(
#     M ~ X,
#     Y ~ X + M
# )
# 
# fit_med <- phybase_run(
#     data = med_data,
#     id_col = "species",
#     tree = tree,
#     equations = equations_med
# )
# 
# summary(fit_med)

## ----me_se, eval=FALSE--------------------------------------------------------
# # Assume we have SEs for X
# X_se <- rep(0.1, N)
# data_se <- list(X = X, Y = Y, X_se = X_se)
# 
# fit_se <- phybase_run(
#     data = data_se,
#     tree = tree,
#     equations = list(Y ~ X),
#     variability = c(X = "se") # Tell PhyBaSE that X has SEs
# )
# 
# summary(fit_se)

## ----me_reps, eval=FALSE------------------------------------------------------
# # Simulate unequal repeated measures (ragged array)
# # Some species have 2 reps, some have 3
# max_reps <- 3
# X_obs <- matrix(NA, nrow = N, ncol = max_reps)
# 
# for (i in 1:N) {
#     # Randomly decide if species has 2 or 3 reps
#     n_i <- sample(2:3, 1)
# 
#     # Simulate data
#     vals <- rnorm(n_i, mean = X[i], sd = 0.2)
# 
#     # Fill matrix (padding with NA for missing reps)
#     X_obs[i, 1:n_i] <- vals
# }
# 
# fit_reps <- phybase_run(
#     data = list(X = X_obs, Y = Y),
#     tree = tree,
#     equations = list(Y ~ X),
#     variability = c(X = "reps")
# )
# 
# summary(fit_reps)
# 
# # PhyBaSE automatically handles the NAs and counts valid replicates per species.

## ----binomial, eval=FALSE-----------------------------------------------------
# # Simulate binary trait
# prob <- 1 / (1 + exp(-(0.5 + 0.8 * X)))
# BinaryTrait <- rbinom(N, 1, prob)
# 
# data_bin <- list(X = X, Bin = BinaryTrait)
# 
# fit_bin <- phybase_run(
#     data = data_bin,
#     tree = tree,
#     equations = list(Bin ~ X),
#     distribution = c(Bin = "binomial"), # Specify distribution
#     WAIC = TRUE, # Use WAIC for model comparison (recommended for non-Gaussian models)
#     DIC = FALSE # DIC may produce NaN penalty with latent variable models
# )
# 
# summary(fit_bin)

## ----multinomial, eval=FALSE--------------------------------------------------
# # Simulate multinomial trait with phylogenetic signal (3 categories)
# # Example: Diet type (Carnivore, Herbivore, Omnivore)
# library(MASS)
# 
# VCV <- ape::vcv.phylo(tree)
# 
# # Phylogenetically structured errors for categories 2 and 3
# # (category 1 is the reference)
# lambda_2 <- 0.6 # Phylogenetic signal for category 2
# tau_u_2 <- 1 / 0.5
# tau_e_2 <- 1 / 0.3
# u_std_2 <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# err_2 <- u_std_2 / sqrt(tau_u_2) + rnorm(N, 0, sqrt(1 / tau_e_2))
# 
# lambda_3 <- 0.7 # Phylogenetic signal for category 3
# tau_u_3 <- 1 / 0.5
# tau_e_3 <- 1 / 0.3
# u_std_3 <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# err_3 <- u_std_3 / sqrt(tau_u_3) + rnorm(N, 0, sqrt(1 / tau_e_3))
# 
# # Latent utilities for each category (with phylogenetic signal)
# L1 <- rep(0, N) # Reference category (Carnivore)
# L2 <- -1 + 0.5 * X + err_2 # Category 2 (Herbivore)
# L3 <- 1 - 0.5 * X + err_3 # Category 3 (Omnivore)
# 
# # Softmax transformation to get probabilities
# exp_L1 <- exp(L1)
# exp_L2 <- exp(L2)
# exp_L3 <- exp(L3)
# sum_exp <- exp_L1 + exp_L2 + exp_L3
# 
# P1 <- exp_L1 / sum_exp # P(Carnivore)
# P2 <- exp_L2 / sum_exp # P(Herbivore)
# P3 <- exp_L3 / sum_exp # P(Omnivore)
# 
# # Sample multinomial outcome
# MultiTrait <- apply(cbind(P1, P2, P3), 1, function(p) sample(1:3, 1, prob = p))
# 
# data_multi <- list(X = X, Multi = MultiTrait) # K is auto-detected
# 
# fit_multi <- phybase_run(
#     data = data_multi,
#     tree = tree,
#     equations = list(Multi ~ X),
#     distribution = c(Multi = "multinomial"),
#     optimise = TRUE
# )
# 
# summary(fit_multi)
# 
# # Check recovery of phylogenetic signal for each category
# # lambda parameters should be ~ 0.6 and 0.7
# fit_multi$summary$statistics[grep("lambda", rownames(fit_multi$summary$statistics)), ]

## ----ordinal, eval=FALSE------------------------------------------------------
# # Simulate IUCN-like conservation status with phylogenetic signal
# # DAG: BodyMass -> IUCN <- HabitatLoss
# #      BodyMass -> HabitatLoss
# # Implied independence: IUCN ⊥ BodyMass | HabitatLoss
# library(MASS)
# 
# # Phylogenetically structured predictors
# VCV <- ape::vcv.phylo(tree)
# 
# # Body Mass (exogenous)
# lambda_mass <- 0.6
# Sigma_mass <- lambda_mass * VCV + (1 - lambda_mass) * diag(N)
# BodyMass <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = Sigma_mass))
# 
# # Habitat Loss (caused by body mass)
# lambda_habitat <- 0.5
# tau_u_habitat <- 1 / 0.4
# tau_e_habitat <- 1 / 0.3
# u_std_habitat <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# HabitatLoss <- 1.2 + 0.7 * BodyMass + u_std_habitat / sqrt(tau_u_habitat) +
#     rnorm(N, 0, sqrt(1 / tau_e_habitat))
# 
# # IUCN status (ordinal response, caused by both predictors)
# lambda_IUCN <- 0.75
# tau_u_IUCN <- 1 / 0.5
# tau_e_IUCN <- 1 / 0.3
# u_std_IUCN <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# phylo_error <- u_std_IUCN / sqrt(tau_u_IUCN) + rnorm(N, 0, sqrt(1 / tau_e_IUCN))
# 
# # Linear predictor: larger mass + more habitat loss -> higher threat
# eta <- -0.5 * BodyMass + 0.8 * HabitatLoss + phylo_error
# cutpoints <- c(-1, 0, 1, 2) # Thresholds between categories
# 
# # Calculate probabilities using cumulative logit
# Q <- matrix(0, nrow = N, ncol = 4)
# for (i in 1:N) {
#     for (k in 1:4) {
#         Q[i, k] <- 1 / (1 + exp(-(cutpoints[k] - eta[i])))
#     }
# }
# 
# # Category probabilities
# P <- matrix(0, nrow = N, ncol = 5)
# P[, 1] <- Q[, 1]
# for (k in 2:4) {
#     P[, k] <- Q[, k] - Q[, k - 1]
# }
# P[, 5] <- 1 - Q[, 4]
# 
# # Sample ordinal outcome (1=LC, 2=NT, 3=VU, 4=EN, 5=CR)
# IUCN <- apply(P, 1, function(p) sample(1:5, 1, prob = p))
# 
# data_ordinal <- list(
#     BodyMass = BodyMass,
#     HabitatLoss = HabitatLoss,
#     IUCN = IUCN,
#     K_IUCN = 5
# )
# 
# # Fit the DAG
# fit_ordinal <- phybase_run(
#     data = data_ordinal,
#     tree = tree,
#     equations = list(
#         IUCN ~ BodyMass + HabitatLoss,
#         HabitatLoss ~ BodyMass
#     ),
#     distribution = c(IUCN = "ordinal"),
#     optimise = TRUE
# )
# 
# summary(fit_ordinal)
# 
# # Check recovery of phylogenetic signal
# fit_ordinal$summary$statistics["lambdaIUCN", ]
# 
# # Test d-separation: IUCN ⊥ BodyMass | HabitatLoss
# fit_dsep <- phybase_run(
#     data = data_ordinal,
#     tree = tree,
#     equations = list(
#         IUCN ~ BodyMass + HabitatLoss,
#         HabitatLoss ~ BodyMass
#     ),
#     distribution = c(IUCN = "ordinal"),
#     dsep = TRUE,
#     optimise = TRUE
# )
# 
# # Check if betaBodyMass is close to zero in the d-sep test
# summary(fit_dsep$samples)
# # If 95% CI of betaBodyMass includes zero -> independence supported

## ----poisson, eval=FALSE------------------------------------------------------
# # Simulate count data with phylogenetic signal
# # Example: Clutch size (number of eggs) in birds
# # DAG: BodyMass -> ClutchSize <- Latitude
# #      BodyMass -> Latitude
# # Implied independence: ClutchSize ⊥ BodyMass | Latitude
# library(MASS)
# 
# # Phylogenetically structured predictors
# VCV <- ape::vcv.phylo(tree)
# 
# # Body mass (log-scale, exogenous)
# lambda_bm <- 0.8
# Sigma_bm <- lambda_bm * VCV + (1 - lambda_bm) * diag(N)
# BodyMass <- as.vector(mvrnorm(1, mu = rep(4, N), Sigma = Sigma_bm))
# 
# # Latitude (caused by body mass - Bergmann's rule proxy)
# lambda_lat <- 0.6
# tau_u_lat <- 1 / 0.4
# tau_e_lat <- 1 / 0.3
# u_std_lat <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# Latitude <- 30 - 0.4 * BodyMass + u_std_lat / sqrt(tau_u_lat) +
#     rnorm(N, 0, sqrt(1 / tau_e_lat))
# 
# # Phylogenetically structured error term for clutch size
# lambda_clutch <- 0.7 # True phylogenetic signal in clutch size
# tau_u <- 1 / 0.5 # Phylogenetic variance
# tau_e <- 1 / 0.3 # Residual variance
# u_std <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# u <- u_std / sqrt(tau_u)
# epsilon <- rnorm(N, 0, sqrt(1 / tau_e))
# phylo_error <- u + epsilon
# 
# # True model: log(E[ClutchSize]) = 1.5 + 0.03*Latitude + phylo_error
# # Larger clutches at higher latitudes, body mass affects ONLY through latitude
# log_clutch <- 1.5 + 0.03 * Latitude + phylo_error
# ClutchSize <- rpois(N, lambda = exp(log_clutch))
# 
# data_poisson <- list(
#     BodyMass = BodyMass,
#     Latitude = Latitude,
#     ClutchSize = ClutchSize
# )
# 
# # Fit the DAG
# fit_poisson <- phybase_run(
#     data = data_poisson,
#     tree = tree,
#     equations = list(
#         ClutchSize ~ Latitude,
#         Latitude ~ BodyMass
#     ),
#     distribution = c(ClutchSize = "poisson"),
#     optimise = TRUE
# )
# 
# summary(fit_poisson)
# 
# # Check recovery of phylogenetic signal
# # lambdaClutchSize should be ~ 0.7
# fit_poisson$summary$statistics["lambdaClutchSize", ]
# 
# # Test d-separation: ClutchSize ⊥ BodyMass | Latitude
# fit_dsep <- phybase_run(
#     data = data_poisson,
#     tree = tree,
#     equations = list(
#         ClutchSize ~ Latitude,
#         Latitude ~ BodyMass
#     ),
#     distribution = c(ClutchSize = "poisson"),
#     dsep = TRUE,
#     optimise = TRUE
# )
# 
# # Check if betaBodyMass is close to zero in the d-sep test
# summary(fit_dsep$samples)
# # If 95% CI of betaBodyMass includes zero -> independence supported

## ----negbinomial, eval=FALSE--------------------------------------------------
# # Simulate heavily overdispersed count data with phylogenetic signal
# # DAG: HostSize -> Parasites <- Stress
# #      HostSize -> Stress
# # Implied independence: Parasites ⊥ HostSize | Stress
# library(MASS)
# 
# # Phylogenetically structured predictors
# VCV <- ape::vcv.phylo(tree)
# 
# # Host Size (exogenous)
# lambda_size <- 0.4
# Sigma_size <- lambda_size * VCV + (1 - lambda_size) * diag(N)
# HostSize <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = Sigma_size))
# 
# # Stress (caused by host size)
# lambda_stress <- 0.5
# tau_u_stress <- 1 / 0.4
# tau_e_stress <- 1 / 0.3
# u_std_stress <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# Stress <- 1.5 + 0.6 * HostSize + u_std_stress / sqrt(tau_u_stress) +
#     rnorm(N, 0, sqrt(1 / tau_e_stress))
# 
# # Phylogenetically structured error term for parasites
# lambda_parasites <- 0.8 # True phylogenetic signal
# tau_u <- 1 / 0.6
# tau_e <- 1 / 0.4
# u_std <- as.vector(mvrnorm(1, mu = rep(0, N), Sigma = solve(VCV)))
# u <- u_std / sqrt(tau_u)
# epsilon <- rnorm(N, 0, sqrt(1 / tau_e))
# phylo_error <- u + epsilon
# 
# # True model with overdispersion
# # Stress increases parasite load (HostSize affects load ONLY through Stress)
# r_true <- 2 # Size parameter (smaller = more overdispersed)
# log_mu <- 2 + 0.7 * Stress + phylo_error
# mu <- exp(log_mu)
# parasite_count <- rnbinom(N, size = r_true, mu = mu)
# 
# data_nb <- list(
#     HostSize = HostSize,
#     Stress = Stress,
#     Parasites = parasite_count
# )
# 
# # Fit the DAG
# fit_nb <- phybase_run(
#     data = data_nb,
#     tree = tree,
#     equations = list(
#         Parasites ~ Stress,
#         Stress ~ HostSize
#     ),
#     distribution = c(Parasites = "negbinomial"),
#     optimise = TRUE
# )
# 
# summary(fit_nb)
# 
# # Check recovery of phylogenetic signal
# # lambdaParasites should be ~ 0.8
# fit_nb$summary$statistics["lambdaParasites", ]
# 
# # Test d-separation: Parasites ⊥ HostSize | Stress
# fit_dsep <- phybase_run(
#     data = data_nb,
#     tree = tree,
#     equations = list(
#         Parasites ~ Stress,
#         Stress ~ HostSize
#     ),
#     distribution = c(Parasites = "negbinomial"),
#     dsep = TRUE,
#     optimise = TRUE
# )
# 
# # Check if betaHostSize is close to zero in the d-sep test
# summary(fit_dsep$samples)
# # If 95% CI of betaHostSize includes zero -> independence supported

## ----categorical, eval=FALSE--------------------------------------------------
# # Data with categorical predictor
# data <- data.frame(
#     Species = tree$tip.label,
#     BodyMass = rnorm(N),
#     Diet = sample(c("Carnivore", "Herbivore", "Omnivore"), N, replace = TRUE),
#     Habitat = factor(sample(c("Forest", "Grassland"), N, replace = TRUE))
# )
# 
# # phybase_format_data automatically creates dummy variables
# data_list <- phybase_format_data(data, tree = tree)
# # Messages:
# # Categorical variable 'Diet' expanded to 2 dummy variable(s) | Reference: 'Carnivore'
# # Categorical variable 'Habitat' expanded to 1 dummy variable(s) | Reference: 'Forest'
# 
# # Use categorical variables DIRECTLY in equations - auto-expands!
# fit <- phybase_run(
#     data = data_list,
#     tree = tree,
#     equations = list(BodyMass ~ Diet + Habitat) # Auto-expands to dummies!
# )
# # Messages:
# # Expanded 'Diet' to: Diet_Herbivore, Diet_Omnivore
# # Expanded 'Habitat' to: Habitat_Grassland
# 
# summary(fit)

## ----missing, eval=FALSE------------------------------------------------------
# # Introduce missing values
# X_miss <- X
# X_miss[c(1, 5, 10)] <- NA # Missing in predictor
# Y_miss <- Y
# Y_miss[c(2, 6, 12)] <- NA # Missing in response
# 
# data_miss <- list(X = X_miss, Y = Y_miss)
# 
# # Just run it! No extra setup needed.
# fit_miss <- phybase_run(
#     data = data_miss,
#     tree = tree,
#     equations = list(Y ~ X)
# )
# 
# summary(fit_miss)
# 
# # PhyBaSE automatically detects NAs and handles them.

## ----multitree, eval=FALSE----------------------------------------------------
# # For demonstration, simulate multiple trees
# num_trees <- 20
# trees <- lapply(1:num_trees, function(x) rtree(N))
# 
# fit_multi <- phybase_run(
#     data = data_list,
#     tree = trees, # Pass list of trees
#     equations = equations
# )
# 
# summary(fit_multi)

## ----parallel, eval=FALSE-----------------------------------------------------
# fit_parallel <- phybase_run(
#     data = data_list,
#     tree = tree,
#     equations = equations,
#     n.iter = 10000,
#     n.chains = 4,
#     parallel = TRUE,
#     n.cores = 4 # Use 4 cores for 4 chains
# )

## ----parallel_compare, eval=FALSE---------------------------------------------
# # Define multiple models
# models <- list(
#     model1 = list(equations = list(Y ~ X)),
#     model2 = list(equations = list(Y ~ 1)) # Null model
# )
# 
# # Run all models in parallel
# results <- phybase_compare(
#     models = models,
#     data = data_list,
#     tree = tree,
#     n.iter = 5000,
#     n.cores = 2 # Run 2 models at a time
# )
# 
# # Compare WAIC
# print(results$comparison)

## ----waic_example, eval=FALSE-------------------------------------------------
# fit <- phybase_run(
#     data = data,
#     structure = tree,
#     id_col = "SP",
#     equations = list(Y ~ X),
#     WAIC = TRUE, # Enable WAIC calculation
#     n.chains = 2, # At least 2 chains recommended
#     n.iter = 2000
# )
# 
# # View WAIC results with standard errors
# fit$WAIC
# #            Estimate   SE
# # elpd_waic   -617.3  12.4
# # p_waic        12.3   3.1
# # waic        1234.5  24.8

## ----waic_comparison, eval=FALSE----------------------------------------------
# # Fit competing models
# fit_complex <- phybase_run(..., WAIC = TRUE)
# fit_simple <- phybase_run(..., WAIC = TRUE)
# 
# # Compare results
# comp <- phybase_compare(fit_complex, fit_simple)
# print(comp)
# 
# # Output:
# #              WAIC   SE dWAIC  dSE p_waic weight
# # fit_complex 212.8 12.4   0.0  0.0   12.3  0.98
# # fit_simple  220.0 11.8   7.2  3.1   10.5  0.02

## ----batch_comparison, eval=FALSE---------------------------------------------
# # Define model specifications
# specs <- list(
#     Full = list(equations = list(Y ~ X1 + X2)),
#     Reduced = list(equations = list(Y ~ X1)),
#     Null = list(equations = list(Y ~ 1))
# )
# 
# # Run batch in parallel
# results <- phybase_compare(
#     model_specs = specs,
#     data = data,
#     tree = tree,
#     n.cores = 3 # Parallel execution
# )
# 
# # View comparison table
# print(results$comparison)

## ----dsep_example, eval=FALSE-------------------------------------------------
# # Define a path model: A -> B -> C
# # This implies: A is independent of C, conditional on B
# A <- rTraitCont(tree)
# B <- 0.5 + 0.8 * A + rTraitCont(tree)
# C <- 0.5 + 0.8 * B + rTraitCont(tree)
# 
# dsep_data <- data.frame(
#     species = tree$tip.label,
#     A = A, B = B, C = C
# )
# 
# equations_dsep <- list(
#     B ~ A,
#     C ~ B
# )
# 
# # Run d-separation tests by setting dsep = TRUE
# fit_dsep <- phybase_run(
#     data = dsep_data,
#     id_col = "species",
#     tree = tree,
#     equations = equations_dsep,
#     dsep = TRUE # <--- This triggers d-sep testing
# )
# 
# # The summary will show the conditional independence tests
# summary(fit_dsep)

## ----latent_mag, eval=FALSE---------------------------------------------------
# # Suppose we have a latent variable L that affects both X and Y
# # L -> X, L -> Y (L is unmeasured)
# 
# equations <- list(
#     X ~ L,
#     Y ~ L
# )
# 
# fit_mag <- phybase_run(
#     data = data_list,
#     tree = tree,
#     equations = equations,
#     latent = "L",
#     latent_method = "correlations", # Default, can omit
#     n.iter = 10000
# )
# 
# # Check induced correlations
# fit_mag$induced_correlations
# # [[1]]
# # [1] "X" "Y"
# 
# # The correlation parameter rho_X_Y is estimated
# summary(fit_mag)

## ----latent_explicit, eval=FALSE----------------------------------------------
# fit_explicit <- phybase_run(
#     data = data_list,
#     tree = tree,
#     equations = list(X ~ L, Y ~ L),
#     latent = "L",
#     latent_method = "explicit", # Model L as a node
#     n.iter = 10000
# )
# 
# # Estimates betaX_L and betaY_L
# summary(fit_explicit)

## ----independent, eval=FALSE--------------------------------------------------
# fit_indep <- phybase_run(
#     data = data_list,
#     structure = NULL, # Independent model
#     equations = equations
# )

## ----spatial, eval=FALSE------------------------------------------------------
# # Create a spatial precision matrix
# spatial_prec <- solve(exp(-dist_matrix / range))
# 
# fit_spatial <- phybase_run(
#     data = data_list,
#     structure = spatial_prec, # Pass matrix directly
#     equations = equations
# )

## ----animal, eval=FALSE-------------------------------------------------------
# # Using MCMCglmm to get A-inverse
# library(MCMCglmm)
# inv_A <- inverseA(pedigree)$Ainv
# 
# fit_animal <- phybase_run(
#     data = data_list,
#     structure = as.matrix(inv_A),
#     equations = equations
# )

## ----multi_struct, eval=FALSE-------------------------------------------------
# fit_multi <- phybase_run(
#     data = data_list,
#     structure = list(
#         phylo = tree,
#         spatial = spatial_cov_matrix
#     ),
#     equations = equations
# )

