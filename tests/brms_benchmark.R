# Benchmark for brms on Mendes Simulation
# Purpose: Measure runtime of brms model for comparison with 'because'

library(brms)
library(tidyverse)
# Using cmdstanr backend as per tutorial preferences if available
# options(brms.backend = "cmdstanr") # defaulting to rstan if cmdstanr not set up

# --- 1. Simulation Function (Copy from tests/mendes_comparison.R) ---
create_count_sim_data = function(
    nlandscapes = 80,
    nsurveys = 240,
    zero_inflation = 0
) {
    effort_min = 50
    effort_max = 5000
    boar_intercept_mean = -2.6
    boar_intercept_sd = 0.96
    boar_crop_beta = 0.5
    boar_NB2_dispersion_parameter = 0.93
    boar_zero_inflation = zero_inflation
    herbivore_intercept_mean = -2.8
    herbivore_intercept_sd = 0.59
    herbivore_boar_beta = 0.0005
    herbivore_human_activity_beta = -0.47
    herbivore_NB2_dispersion_parameter = 1.96
    herbivore_zero_inflation = zero_inflation
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
    sim_data$effort_scaled = scale(sim_data$effort)

    sim_data$human_activity = runif(nsurveys)
    sim_data$crops = runif(nsurveys)
    sim_data$forest_cover = runif(nsurveys)

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
sim_data <- create_count_sim_data(
    nlandscapes = 80,
    nsurveys = 240,
    zero_inflation = 0
)

# --- 3. Define brms Model (Exact syntax from tutorial) ---
cat("Compiling and Running brms model...\n")
cat("(This may take a while due to compilation...)\n")

# Model definitions from tutorial
mod1 <- bf(
    boars ~ human_activity +
        forest_cover +
        crops +
        effort_scaled +
        (1 | landscapeID),
    family = negbinomial()
)
mod2 <- bf(
    herbivores ~ boars +
        human_activity +
        forest_cover +
        crops +
        effort_scaled +
        (1 | landscapeID),
    family = negbinomial()
)
mod3 <- bf(
    predators ~ herbivores +
        boars +
        human_activity +
        forest_cover +
        effort_scaled +
        (1 | landscapeID),
    family = negbinomial()
)

# Timing the run
start_time <- Sys.time()

# Using exact settings from tutorial: iter=4000, warmup=1500, chains=4
# Tutorial uses threads=2, cores=6. We'll match threads=2 but cap cores to 4 to be safe/fair if running locally.
brms_count <- brm(
    mod1 + mod2 + mod3 + set_rescor(FALSE),
    prior = c(
        prior(student_t(4, 0, 5), coef = "crops", resp = "boars"),
        prior(student_t(4, 0, 5), coef = "forest_cover", resp = "boars"),
        prior(student_t(4, 0, 5), coef = "human_activity", resp = "boars"),
        prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "boars"),
        prior(student_t(4, 0, 5), coef = "boars", resp = "herbivores"),
        prior(student_t(4, 0, 5), coef = "crops", resp = "herbivores"),
        prior(student_t(4, 0, 5), coef = "forest_cover", resp = "herbivores"),
        prior(student_t(4, 0, 5), coef = "human_activity", resp = "herbivores"),
        prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "herbivores"),
        prior(student_t(4, 0, 5), coef = "boars", resp = "predators"),
        prior(student_t(4, 0, 5), coef = "herbivores", resp = "predators"),
        prior(student_t(4, 0, 5), coef = "forest_cover", resp = "predators"),
        prior(student_t(4, 0, 5), coef = "human_activity", resp = "predators"),
        prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "predators")
    ),
    data = sim_data,
    iter = 4000,
    warmup = 1500,
    threads = threading(2),
    cores = 4,
    chains = 4,
    init = 0,
    control = list(max_treedepth = 10, adapt_delta = 0.80),
    refresh = 500
)

end_time <- Sys.time()
duration <- end_time - start_time

cat("\n--- brms Benchmark Results ---\n")
print(duration)
cat("\nSummary:\n")
print(summary(brms_count))
