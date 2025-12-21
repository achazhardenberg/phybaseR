# Demo: Hierarchical Data Analysis with because
# ===============================================
# This script demonstrates how to use the new hierarchical data support features.
# - Auto-detection of levels
# - Cross-level interactions
# - Polynomials in hierarchical models

library(because)
set.seed(42)

# 1. Simulate Hierarchical Data
# -----------------------------
# Hierarchy: Species (Level 1) -> Individual (Level 2)
# We want to model Individual Size as a function of:
# - Species Trait (Macro-level predictor)
# - Individual Age (Micro-level predictor, non-linear)

# Level 1: Species (N = 20)
N_Species <- 20
df_species <- data.frame(
    SpeciesID = paste0("sp", 1:N_Species),
    Trait = rnorm(N_Species, mean = 10, sd = 2)
)

# Level 2: Individuals (N = 400)
# Each species has ~20 individuals
N_Ind <- 400
df_individuals <- data.frame(
    IndID = 1:N_Ind,
    SpeciesID = sample(df_species$SpeciesID, N_Ind, replace = TRUE),
    Age = runif(N_Ind, 1, 10)
)

# Generate Outcome (Size)
# True Model: Size ~ 2.5 * Trait + 0.5 * Age^2 + Noise
# Note: Trait is from Species level, Age is from Individual level
species_map <- setNames(df_species$Trait, df_species$SpeciesID)
trait_mapped <- species_map[df_individuals$SpeciesID]

mu <- 2.5 * trait_mapped + 0.5 * df_individuals$Age^2
df_individuals$Size <- rnorm(N_Ind, mean = mu, sd = 2)

# 2. Specify the Model
# --------------------
# We use standard formula syntax.
# because() will figure out which variable belongs to which level.
eq_list <- list(
    Size ~ Trait + Age + I(Age^2)
)

# 3. Fit the Model (Auto-Detection)
# ---------------------------------
# Simply pass the dataframes in a list.
# The function will match 'SpeciesID' to link them.
print("Fitting hierarchical model...")

fit <- because(
    equations = eq_list,
    data = list(
        Species = df_species, # Level 1 Data
        Individuals = df_individuals # Level 2 Data
    ),
    n.iter = 2000
)

# 4. Inspect Results
summary(fit)

# Interpretation:
# - beta_Size_Trait: Should be close to 2.5
# - beta_Size_Age_pow2: Should be close to 0.5
# - N-inflation check: The model used separate loops for Species (N=20) and Individuals (N=400).

# 5. Diagnostic Output (What happened under the hood?)
# ----------------------------------------------------
# The package auto-detected:
# Level 'Species' (N=20): Trait
# Level 'Individuals' (N=400): Size, Age
# Hierarchy: Species > Individuals
# Polynomial 'Age^2' was correctly handled at the Individual level.
