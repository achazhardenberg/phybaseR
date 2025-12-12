# Replication of Mendes et al. (2025) Simulation
# Purpose: Test 'because' package against the simulation scenarios from:
# "The use, misuse and opportunities for structural equation modelling (SEM) in wildlife ecology"

# library(because)
devtools::load_all()
library(tidyverse)

# --- 1. Simulation Function (Adapted directly from Mendes et al.) ---
create_count_sim_data = function(
    nlandscapes = 80, # Number of simulated landscapes
    nsurveys = 240, # Total number of surveys
    zero_inflation = 0 # Zero-inflation parameter
) {
    ### Unchangeable parameters
    effort_min = 50
    effort_max = 5000

    # Boar parameters
    boar_intercept_mean = -2.6
    boar_intercept_sd = 0.96
    boar_crop_beta = 0.5
    boar_NB2_dispersion_parameter = 0.93
    boar_zero_inflation = zero_inflation

    # Herbivore parameters
    herbivore_intercept_mean = -2.8
    herbivore_intercept_sd = 0.59
    herbivore_boar_beta = 0.0005
    herbivore_human_activity_beta = -0.47
    herbivore_NB2_dispersion_parameter = 1.96
    herbivore_zero_inflation = zero_inflation

    # Predator parameters
    predator_intercept_mean = -3.16
    predator_intercept_sd = 0.64
    predator_human_activity_beta = -0.19
    predator_NB2_dispersion_parameter = 1.74
    predator_zero_inflation = zero_inflation

    sim_data = as.data.frame(matrix(NA, nrow = nsurveys, ncol = 12))
    colnames(sim_data) = c(
        "surveyID",
        "landscapeID",
        "effort",
        "forest_cover",
        "crops",
        "human_activity",
        "boars",
        "herbivores",
        "predators",
        "boar_intercept",
        "herbivore_intercept",
        "predator_intercept"
    )

    sim_data$surveyID = 1:nsurveys
    sim_data$landscapeID = as.factor(sort(rep(
        1:nlandscapes,
        ceiling(nsurveys / nlandscapes)
    )[1:nsurveys]))

    sim_data$effort = as.integer(runif(nsurveys, effort_min, effort_max))
    sim_data$log_effort = log(sim_data$effort) # Pre-calculate for model

    ### make the covariates
    sim_data$human_activity = runif(nsurveys)
    sim_data$crops = runif(nsurveys)
    sim_data$forest_cover = runif(nsurveys)

    ### Boars
    intercepts = rnorm(nlandscapes, boar_intercept_mean, boar_intercept_sd)
    sim_data$boar_intercept = rep(NA, nsurveys)
    for (i in 1:nlandscapes) {
        sim_data$boar_intercept[sim_data$landscapeID == i] = intercepts[i]
    }

    for (i in 1:nsurveys) {
        eta = sim_data$boar_intercept[i] +
            boar_crop_beta * sim_data$crops[i] +
            log(sim_data$effort[i])
        sim_data$boars[i] = rnbinom(
            1,
            mu = exp(eta),
            size = boar_NB2_dispersion_parameter
        )
    }

    ### Herbivores
    intercepts = rnorm(
        nlandscapes,
        herbivore_intercept_mean,
        herbivore_intercept_sd
    )
    sim_data$herbivore_intercept = rep(NA, nsurveys)
    for (i in 1:nlandscapes) {
        sim_data$herbivore_intercept[sim_data$landscapeID == i] = intercepts[i]
    }

    for (i in 1:nsurveys) {
        eta = sim_data$herbivore_intercept[i] +
            herbivore_boar_beta * sim_data$boars[i] +
            herbivore_human_activity_beta * sim_data$human_activity[i] +
            log(sim_data$effort[i])
        sim_data$herbivores[i] = rnbinom(
            1,
            mu = exp(eta),
            size = herbivore_NB2_dispersion_parameter
        )
    }

    ### Predators
    intercepts = rnorm(
        nlandscapes,
        predator_intercept_mean,
        predator_intercept_sd
    )
    sim_data$predator_intercept = rep(NA, nsurveys)
    for (i in 1:nlandscapes) {
        sim_data$predator_intercept[sim_data$landscapeID == i] = intercepts[i]
    }

    for (i in 1:nsurveys) {
        eta = sim_data$predator_intercept[i] +
            predator_human_activity_beta * sim_data$human_activity[i] +
            log(sim_data$effort[i])
        sim_data$predators[i] = rnbinom(
            1,
            mu = exp(eta),
            size = predator_NB2_dispersion_parameter
        )
    }

    ### Zero Inflation
    # NOTE: 'because' uses standard Negative Binomial, so high zero-inflation might be mapped to overdispersion
    # but won't be perfectly modeled as a separate process.
    if (zero_inflation > 0) {
        zinf_boar = 1 - as.numeric(rbinom(nsurveys, 1, p = boar_zero_inflation))
        sim_data$boars = sim_data$boars * zinf_boar

        zinf_herb = 1 -
            as.numeric(rbinom(nsurveys, 1, p = herbivore_zero_inflation))
        sim_data$herbivores = sim_data$herbivores * zinf_herb

        zinf_pred = 1 -
            as.numeric(rbinom(nsurveys, 1, p = predator_zero_inflation))
        sim_data$predators = sim_data$predators * zinf_pred
    }

    return(sim_data)
}

# --- 2. Generate Data ---
set.seed(42)
cat("Generating simulated data...\n")
# Using default parameters from paper's "Section 1"
# nlandscapes = 80, nsurveys = 240, zero_inflation = 0 (Start with 0 to test baseline)
data <- create_count_sim_data(
    nlandscapes = 80,
    nsurveys = 240,
    zero_inflation = 0.3
)

# --- 3. Run 'because' Model ---
cat("Fitting 'because' model (ZINB)...\n")
start_time <- Sys.time()

# Model structure from paper:
# Boars ~ Crops (+ random landscape)
# Herbivores ~ Boars + HumanActivity (+ random landscape)
# Predators ~ HumanActivity (+ random landscape)
# All have offset(log(effort)), which we model as a covariate 'log_effort'

mod <- because(
    equations = list(
        boars ~ crops + log_effort,
        herbivores ~ boars + human_activity + log_effort,
        predators ~ human_activity + log_effort
    ),
    data = data,
    distribution = list(
        boars = "zinb",
        herbivores = "zinb",
        predators = "zinb"
    ),
    random = ~ (1 | landscapeID), # Global random effect for landscape
    n.iter = 5000,
    n.burnin = 1000,
    n.chains = 3,
    quiet = TRUE
)
end_time <- Sys.time()
duration <- end_time - start_time
cat(paste(
    "Model runtime:",
    round(as.numeric(duration, units = "mins"), 2),
    "mins\n"
))

# --- 4. Compare Results ---
cat("\n--- Results Comparison ---\n")
print(summary(mod))

# Extract betas
# Extract betas
sum_obj <- summary(mod)
if (
    inherits(sum_obj, "summary.because") ||
        (is.list(sum_obj) && "results" %in% names(sum_obj))
) {
    sum_tab <- sum_obj$results
} else {
    sum_tab <- sum_obj
}

betas <- sum_tab[grep("beta", rownames(sum_tab)), "Mean"]
names(betas) <- rownames(sum_tab)[grep("beta", rownames(sum_tab))]

# True values
true_vals <- c(
    "beta_boars_crops" = 0.5,
    "beta_herbivores_boars" = 0.0005,
    "beta_herbivores_human_activity" = -0.47,
    "beta_predators_human_activity" = -0.19,
    "beta_log_effort" = 1.0 # Theoretically should be 1.0
)

cat("\n--- True vs Estimated Coefficients ---\n")
cat(sprintf(
    "%-30s %-10s %-10s %-10s\n",
    "Parameter",
    "True",
    "Estimated",
    "Bias"
))
cat(sprintf(
    "%-30s %-10s %-10s %-10s\n",
    "---------",
    "----",
    "---------",
    "----"
))

# Mapping 'because' parameter names to Truth
# because names: beta_RESPONSE_PREDICTOR
mapping <- list(
    "beta_boars_crops" = "beta_boars_crops",
    "beta_herbivores_boars" = "beta_herbivores_boars",
    "beta_herbivores_human_activity" = "beta_herbivores_human_activity",
    "beta_predators_human_activity" = "beta_predators_human_activity"
    # log_effort betas are separate for each response
)

for (true_name in names(mapping)) {
    model_name <- mapping[[true_name]]
    if (model_name %in% names(betas)) {
        est <- betas[[model_name]]
        truth <- true_vals[[true_name]]
        bias <- est - truth
        cat(sprintf(
            "%-30s %-10.4f %-10.4f %-10.4f\n",
            true_name,
            truth,
            est,
            bias
        ))
    }
}

# Check explicit log_effort coefficients
cat("\n--- Offset Checks (log_effort ~ 1.0) ---\n")
effort_betas <- betas[grep("log_effort", names(betas))]
for (name in names(effort_betas)) {
    est <- effort_betas[[name]]
    cat(sprintf("%-30s %-10.4f %-10.4f %-10.4f\n", name, 1.0, est, est - 1.0))
}

cat("\nDone.\n")
