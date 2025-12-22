rm(list = ls())

wd = getwd()
wd

file_storade_wd = "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/file_storade_wd"


##### Prepare the environment #####

library(tidyverse)
library(piecewiseSEM)
library(glmmTMB)
library(lavaan)
library(seminr)
library(brms)
library(kSamples)
library(because)


library(cmdstanr)
#set_cmdstan_path(path = "set path to cmdstan, if needed")

options(
  file_refit = "on_change",
  brms.backend = "cmdstanr",
  mc.cores = parallel::detectCores(),
  threads = 4
)

options(scipen = 999)


### Iterations###

Inter.num = 10 #Set the number of iterations


### Simulation Scenarios

nlandscapes_values = c(80) #c(10, 20, 40, 80, 120)
nsurveys_value = c(600) #c(200, 400, 600, 800, 1000, 1200)
zero_inflation_value = c(0.3) #c(0, 0.3, 0.5)


##### Section 1 - Simulation function and parameters #####

create_count_sim_data = function(
  nlandscapes = 80, # Number of simulated landscapes
  nsurveys = 240, #Total number of surveys
  zero_inflation = 0 #Zero-inflation
) {
  ###Unchangeable parameters

  effort_min = 50 #Minimum effort per survey
  effort_max = 5000 #Maximum effort per survey
  # boar parameters
  boar_intercept_mean = -2.6 #for the random intercepts
  boar_intercept_sd = 0.96 #for of the random intercepts
  boar_crop_beta = 0.5 #Crop effect on boar
  boar_NB2_dispersion_parameter = 0.93 #Negative binomial dispersion parameter
  boar_zero_inflation = zero_inflation #zero inflation (0 to 1; 1 = 100%)
  #herbivore parameters
  herbivore_intercept_mean = -2.8 #for the random intercepts
  herbivore_intercept_sd = 0.59 #for the random intercepts
  herbivore_boar_beta = 0.0005 #boar effect on herbivores
  herbivore_human_activity_beta = -0.47 #Human activity effect on herbivores
  herbivore_NB2_dispersion_parameter = 1.96 #Negative binomial dispersion parameter
  herbivore_zero_inflation = zero_inflation #zero inflation (0 to 1; 1 = 100%)
  #predator parameters
  predator_intercept_mean = -3.16 #for the random intercepts
  predator_intercept_sd = 0.64 #for the random intercepts
  predator_human_activity_beta = -0.19 #Human activity effect on predators
  predator_NB2_dispersion_parameter = 1.74 #Negative binomial dispersion parameter
  predator_zero_inflation = zero_inflation #zero inflation (0 to 1; 1 = 100%)

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

  sim_data$surveyID = 1:nsurveys #make the survey ID
  sim_data$landscapeID = as.factor(sort(rep(
    1:nlandscapes,
    ceiling(nsurveys / nlandscapes)
  )[1:nsurveys])) #Make the landscape IDs.

  sim_data$effort = as.integer(runif(nsurveys, effort_min, effort_max)) #make the effort as an integer within the predefined limits

  ### make the covariates
  sim_data$human_activity = runif(nsurveys)
  sim_data$crops = runif(nsurveys)
  sim_data$forest_cover = runif(nsurveys)

  ### make the boar capture data
  intercepts = rnorm(nlandscapes, boar_intercept_mean, boar_intercept_sd)
  sim_data$boar_intercept = rep(NA, nsurveys)

  for (i in 1:nlandscapes) {
    sim_data$boar_intercept[sim_data$landscapeID == i] = intercepts[i]
  }
  rm(intercepts)

  for (i in 1:nsurveys) {
    #eta = X %*% beta   #eta is log(mu)!
    #y = rnbinom(length(eta), mean = exp(eta), size = theta)

    eta = sim_data$boar_intercept[i] +
      boar_crop_beta * sim_data$crops[i] +
      log(sim_data$effort[i])
    sim_data$boars[i] = rnbinom(
      1,
      mu = exp(eta),
      size = boar_NB2_dispersion_parameter
    )
  }
  rm(i)

  ### make the herbivore capture data
  intercepts = rnorm(
    nlandscapes,
    herbivore_intercept_mean,
    herbivore_intercept_sd
  )
  sim_data$herbivore_intercept = rep(NA, nsurveys)

  for (i in 1:nlandscapes) {
    sim_data$herbivore_intercept[sim_data$landscapeID == i] = intercepts[i]
  }
  rm(intercepts)

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
  rm(i)

  ### make the carnivore capture data
  intercepts = rnorm(
    nlandscapes,
    predator_intercept_mean,
    predator_intercept_sd
  )
  sim_data$predator_intercept = rep(NA, nsurveys)

  for (i in 1:nlandscapes) {
    sim_data$predator_intercept[sim_data$landscapeID == i] = intercepts[i]
  }
  rm(intercepts)

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
  rm(i)

  ### Add zero Inflation (This is done here because it is a detection process)
  zinf = 1 - as.numeric(rbinom(nsurveys, 1, p = boar_zero_inflation))
  sim_data$boars = sim_data$boars * zinf

  zinf = 1 - as.numeric(rbinom(nsurveys, 1, p = herbivore_zero_inflation))
  sim_data$herbivores = sim_data$herbivores * zinf

  zinf = 1 - as.numeric(rbinom(nsurveys, 1, p = predator_zero_inflation))
  sim_data$predators = sim_data$predators * zinf

  ### Make the RAI
  sim_data$boars_rai = (sim_data$boars / sim_data$effort) * 100
  sim_data$herbivores_rai = (sim_data$herbivores / sim_data$effort) * 100
  sim_data$predators_rai = (sim_data$predators / sim_data$effort) * 100

  sim_data$effort_scaled = scale(sim_data$effort)

  ### Make the Log RAI
  sim_data$log_boars_rai = log(sim_data$boars_rai + 1)
  sim_data$log_herbivores_rai = log(sim_data$herbivores_rai + 1)
  sim_data$log_predators_rai = log(sim_data$predators_rai + 1)

  ### Make the Log Count
  sim_data$log_boars = log(sim_data$boars + 1)
  sim_data$log_herbivores = log(sim_data$herbivores + 1)
  sim_data$log_predators = log(sim_data$predators + 1)

  return(sim_data)
}


### Create simulated the data
#sim_data = create_count_sim_data()
#head(sim_data)

##### Section 2 - Making a log of all scenarios and SEM approaches to be tested #####

### Make the log
model_list = c(
  paste("nlandscapes", nlandscapes_values, sep = ""),
  paste("nsurveys", nsurveys_value, sep = ""),
  paste("zero_inflation", zero_inflation_value, sep = "")
)

approaches = c(
  "brms_zi_count",
  "brms_zi_med_count",
  #"brms_rai",
  #"brms_zi_med_rai",
  #"piecewiseSEM_count",
  #"piecewiseSEM_zi_count",
  #"piecewiseSEM_zi_med_count",
  #"piecewiseSEM_rai",
  #"piecewiseSEM_zi_med_rai",
  "lavaan_count",
  "lavaan_zi_med_count",
  #"lavaan_rai",
  #"lavaan_zi_med_rai",
  "SEMinR_count",
  #"SEMinR_rai",
  "SEMinR_zi_med_count",
  #"SEMinR_zi_med_rai",
  "because_zinb_count",
  "because_zi_med_count"
)


model_list_log = data.frame(
  nID = seq(1, length(approaches) * length(model_list)),
  "target_parameter" = sort(rep(model_list, length(approaches))), #Change to "test_parameter"?
  "approach" = rep(approaches, length(model_list))
) %>%
  mutate(
    model_ID = paste(target_parameter, approach, sep = "-"),
    nlandscapes = 80, # default number of simulated landscapes
    nsurveys = 240, #default number of surveys
    zero_inflation = 0 #default zero inflation
  )


### Add the values which differ from the default
for (i in nlandscapes_values) {
  model_list_log$nlandscapes[
    model_list_log$target_parameter == paste("nlandscapes", i, sep = "")
  ] = i
}

for (i in nsurveys_value) {
  model_list_log$nsurveys[
    model_list_log$target_parameter == paste("nsurveys", i, sep = "")
  ] = i
}

for (i in zero_inflation_value) {
  model_list_log$zero_inflation[
    model_list_log$target_parameter == paste("zero_inflation", i, sep = "")
  ] = i
}

rm(model_list, approaches)


model_output_list = list()


##### Section 3 - Main loop (Simulate, test, save, repeat) #####

for (i in 1:nrow(model_list_log)) {
  m = model_list_log$model_ID[i]

  if (model_list_log$approach[i] == "brms_count") {
    ##### BRMS COUNT #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    # Create the data for the first model
    while (length(output) < 1) {
      tryCatch(
        {
          sim_data = create_count_sim_data(
            nlandscapes = model_list_log$nlandscapes[i],
            nsurveys = model_list_log$nsurveys[i],
            zero_inflation = model_list_log$zero_inflation[i]
          )

          #Create the first model
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

          brms_count <- brm(
            mod1 + mod2 + mod3 + set_rescor(FALSE),

            prior = c(
              prior(student_t(4, 0, 5), coef = "crops", resp = "boars"),
              prior(student_t(4, 0, 5), coef = "forest_cover", resp = "boars"),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "boars"
              ),
              prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "boars"),

              prior(student_t(4, 0, 5), coef = "boars", resp = "herbivores"),
              prior(student_t(4, 0, 5), coef = "crops", resp = "herbivores"),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "herbivores"
              ),

              prior(student_t(4, 0, 5), coef = "boars", resp = "predators"),
              prior(
                student_t(4, 0, 5),
                coef = "herbivores",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "predators"
              )
            ),

            data = sim_data,
            iter = 4000,
            warmup = 1500,
            threads = 2,
            cores = 6,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.80)
          )

          #save the output
          output[[length(output) + 1]] = summary(brms_count)
        },
        error = function(e) {}
      )
    }

    #REPEAT IT USING UPDATE
    while (length(output) < Inter.num) {
      print(paste(m, "- iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          brms_count = update(brms_count, newdata = sim_data)

          output[[length(output) + 1]] = summary(brms_count)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "brms_zi_count") {
    ##### BRMS ZI COUNT #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    # Create the data for the first model
    while (length(output) < 1) {
      tryCatch(
        {
          sim_data = create_count_sim_data(
            nlandscapes = model_list_log$nlandscapes[i],
            nsurveys = model_list_log$nsurveys[i],
            zero_inflation = model_list_log$zero_inflation[i]
          )

          #Create the first model
          mod1 <- bf(
            boars ~ human_activity +
              forest_cover +
              crops +
              effort_scaled +
              (1 | landscapeID),
            family = zero_inflated_negbinomial()
          )
          mod2 <- bf(
            herbivores ~ boars +
              human_activity +
              forest_cover +
              crops +
              effort_scaled +
              (1 | landscapeID),
            family = zero_inflated_negbinomial()
          )
          mod3 <- bf(
            predators ~ herbivores +
              boars +
              human_activity +
              forest_cover +
              effort_scaled +
              (1 | landscapeID),
            family = zero_inflated_negbinomial()
          )

          brms_count_zi.dist <- brm(
            mod1 + mod2 + mod3 + set_rescor(FALSE),

            prior = c(
              prior(student_t(4, 0, 5), coef = "crops", resp = "boars"),
              prior(student_t(4, 0, 5), coef = "forest_cover", resp = "boars"),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "boars"
              ),
              prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "boars"),

              prior(student_t(4, 0, 5), coef = "boars", resp = "herbivores"),
              prior(student_t(4, 0, 5), coef = "crops", resp = "herbivores"),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "herbivores"
              ),

              prior(student_t(4, 0, 5), coef = "boars", resp = "predators"),
              prior(
                student_t(4, 0, 5),
                coef = "herbivores",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "predators"
              )
            ),

            data = sim_data,
            iter = 4000,
            warmup = 1500,
            threads = 2,
            cores = 6,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.80)
          )

          #save the output
          output[[length(output) + 1]] = summary(brms_count_zi.dist)
        },
        error = function(e) {}
      )
    }

    #REPEAT IT USING UPDATE
    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          brms_count_zi.dist = update(brms_count_zi.dist, newdata = sim_data)

          output[[length(output) + 1]] = summary(brms_count_zi.dist)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "brms_zi_med_count") {
    ##### BRMS ZI COUNT (Jiang´s approach) #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    # Create the data for the first model
    while (length(output) < 1) {
      tryCatch(
        {
          sim_data = create_count_sim_data(
            nlandscapes = model_list_log$nlandscapes[i],
            nsurveys = model_list_log$nsurveys[i],
            zero_inflation = model_list_log$zero_inflation[i]
          )
          sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
          sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

          #Create the first model
          mod <-
            bf(
              boars ~ human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = zero_inflated_negbinomial()
            ) +
            bf(
              herbivores ~ boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = zero_inflated_negbinomial()
            ) +
            bf(
              predators ~ herbivores +
                herbivores_nonzero +
                boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                effort_scaled +
                (1 | landscapeID),
              family = zero_inflated_negbinomial()
            ) +
            set_rescor(FALSE)

          brms_count_zi.jiang <- brm(
            mod,

            prior = c(
              prior(student_t(4, 0, 5), coef = "crops", resp = "boars"),
              prior(student_t(4, 0, 5), coef = "forest_cover", resp = "boars"),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "boars"
              ),
              prior(student_t(4, 0, 5), coef = "effort_scaled", resp = "boars"),

              prior(student_t(4, 0, 5), coef = "boars", resp = "herbivores"),
              prior(
                student_t(4, 0, 5),
                coef = "boars_nonzero",
                resp = "herbivores"
              ),
              prior(student_t(4, 0, 5), coef = "crops", resp = "herbivores"),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "herbivores"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "herbivores"
              ),

              prior(student_t(4, 0, 5), coef = "boars", resp = "predators"),
              prior(
                student_t(4, 0, 5),
                coef = "boars_nonzero",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "herbivores",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "herbivores_nonzero",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "predators"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "effort_scaled",
                resp = "predators"
              )
            ),

            data = sim_data,
            iter = 4000,
            warmup = 1500,
            threads = 2,
            cores = 6,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.80)
          )

          #save the output
          output[[length(output) + 1]] = summary(brms_count_zi.jiang)
        },
        error = function(e) {}
      )
    }

    #REPEAT IT USING UPDATE
    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          brms_count_zi.jiang = update(brms_count_zi.jiang, newdata = sim_data)

          output[[length(output) + 1]] = summary(brms_count_zi.jiang)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "brms_rai") {
    ##### BRMS RAI #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    # Create the data for the first model
    while (length(output) < 1) {
      tryCatch(
        {
          sim_data = create_count_sim_data(
            nlandscapes = model_list_log$nlandscapes[i],
            nsurveys = model_list_log$nsurveys[i],
            zero_inflation = model_list_log$zero_inflation[i]
          )

          #Create the first model
          mod1 <- bf(
            log_boars_rai ~ human_activity +
              forest_cover +
              crops +
              (1 | landscapeID),
            family = student()
          )
          mod2 <- bf(
            log_herbivores_rai ~ log_boars_rai +
              human_activity +
              forest_cover +
              crops +
              (1 | landscapeID),
            family = student()
          )
          mod3 <- bf(
            log_predators_rai ~ log_herbivores_rai +
              log_boars_rai +
              human_activity +
              forest_cover +
              (1 | landscapeID),
            family = student()
          )

          brms_RAI <- brm(
            mod1 + mod2 + mod3 + set_rescor(FALSE),
            prior = c(
              prior(student_t(4, 0, 5), coef = "crops", resp = "logboarsrai"),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "logboarsrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "logboarsrai"
              ),

              prior(
                student_t(4, 0, 5),
                coef = "log_boars_rai",
                resp = "logherbivoresrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "crops",
                resp = "logherbivoresrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "logherbivoresrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "logherbivoresrai"
              ),

              prior(
                student_t(4, 0, 5),
                coef = "log_boars_rai",
                resp = "logpredatorsrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "log_herbivores_rai",
                resp = "logpredatorsrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "forest_cover",
                resp = "logpredatorsrai"
              ),
              prior(
                student_t(4, 0, 5),
                coef = "human_activity",
                resp = "logpredatorsrai"
              )
            ),

            data = sim_data,
            iter = 4000,
            warmup = 1500,
            threads = 2,
            cores = 6,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.90)
          )

          #save the output
          output[[length(output) + 1]] = summary(brms_RAI)
        },
        error = function(e) {}
      )
    }

    #REPEAT IT USING UPDATE
    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          brms_RAI = update(brms_RAI, newdata = sim_data)

          output[[length(output) + 1]] = summary(brms_RAI)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "brms_zi_rai_jiang") {
    ##### BRMS ZI RAI (Jiang´s approach) #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    # Create the data for the first model

    sim_data = create_count_sim_data(
      nlandscapes = model_list_log$nlandscapes[i],
      nsurveys = model_list_log$nsurveys[i],
      zero_inflation = model_list_log$zero_inflation[i]
    )
    sim_data$boars_nonzero = as.numeric(sim_data$log_boars_rai > 0)
    sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

    #Create the first model
    mod <-
      bf(
        log_boars_rai ~ human_activity +
          forest_cover +
          crops +
          (1 | landscapeID),
        family = student()
      ) +
      bf(
        log_herbivores_rai ~ log_boars_rai +
          boars_nonzero +
          human_activity +
          forest_cover +
          crops +
          (1 | landscapeID),
        family = student()
      ) +
      bf(
        log_predators_rai ~ log_herbivores_rai +
          herbivores_nonzero +
          log_boars_rai +
          boars_nonzero +
          human_activity +
          forest_cover +
          (1 | landscapeID),
        family = student()
      ) +
      set_rescor(FALSE)

    brms_ZI_RAI <- brm(
      mod,

      prior = c(
        prior(student_t(4, 0, 5), coef = "crops", resp = "logboarsrai"),
        prior(student_t(4, 0, 5), coef = "forest_cover", resp = "logboarsrai"),
        prior(
          student_t(4, 0, 5),
          coef = "human_activity",
          resp = "logboarsrai"
        ),

        prior(
          student_t(4, 0, 5),
          coef = "log_boars_rai",
          resp = "logherbivoresrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "boars_nonzero",
          resp = "logherbivoresrai"
        ),
        prior(student_t(4, 0, 5), coef = "crops", resp = "logherbivoresrai"),
        prior(
          student_t(4, 0, 5),
          coef = "forest_cover",
          resp = "logherbivoresrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "human_activity",
          resp = "logherbivoresrai"
        ),

        prior(
          student_t(4, 0, 5),
          coef = "log_boars_rai",
          resp = "logpredatorsrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "boars_nonzero",
          resp = "logpredatorsrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "log_herbivores_rai",
          resp = "logpredatorsrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "herbivores_nonzero",
          resp = "logpredatorsrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "forest_cover",
          resp = "logpredatorsrai"
        ),
        prior(
          student_t(4, 0, 5),
          coef = "human_activity",
          resp = "logpredatorsrai"
        )
      ),

      data = sim_data,
      iter = 4000,
      warmup = 1500,
      threads = 2,
      cores = 6,
      chains = 4,
      control = list(max_treedepth = 10, adapt_delta = 0.80)
    )

    #save the output
    output[[length(output) + 1]] = summary(brms_ZI_RAI)

    #REPEAT IT USING UPDATE
    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$log_boars_rai > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          brms_ZI_RAI = update(brms_ZI_RAI, newdata = sim_data)

          output[[length(output) + 1]] = summary(brms_ZI_RAI)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "piecewiseSEM_count") {
    ##### PieacewiseSEM COUNT  #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          psem_mod_count <- psem(
            glmmTMB(
              boars ~ human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              data = sim_data
            ),
            glmmTMB(
              herbivores ~ boars +
                human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              data = sim_data
            ),
            glmmTMB(
              predators ~ herbivores +
                boars +
                human_activity +
                forest_cover +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              data = sim_data
            ),
            sim_data
          )

          output[[length(output) + 1]] = summary(
            psem_mod_count,
            standardize.type = "Mendard.OE",
            conserve = T
          )
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "piecewiseSEM_zi_count") {
    ##### PieacewiseSEM ZI COUNT  #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          psem_mod_zi_count <- psem(
            glmmTMB(
              boars ~ human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            glmmTMB(
              herbivores ~ boars +
                human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            glmmTMB(
              predators ~ herbivores +
                boars +
                human_activity +
                forest_cover +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            sim_data
          )

          output[[length(output) + 1]] = summary(
            psem_mod_zi_count,
            standardize.type = "Mendard.OE",
            conserve = T
          )
        },
        error = function(e) {
          print(e)
        }
      )
    }
  } else if (model_list_log$approach[i] == "piecewiseSEM_zi_med_count") {
    ##### PieacewiseSEM ZI COUNT (Jiang´s approach) #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          psem_mod_count_zi.jiang <- psem(
            glmmTMB(
              boars ~ human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            glmmTMB(
              herbivores ~ boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                crops +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            glmmTMB(
              predators ~ herbivores +
                herbivores_nonzero +
                boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                effort_scaled +
                (1 | landscapeID),
              family = nbinom2,
              zi = ~1,
              data = sim_data
            ),
            sim_data
          )

          output[[length(output) + 1]] = summary(
            psem_mod_count_zi.jiang,
            standardize.type = "Mendard.OE",
            conserve = T
          )
        },
        error = function(e) {
          print(e)
        }
      )
    }
  } else if (model_list_log$approach[i] == "piecewiseSEM_rai") {
    ##### PieacewiseSEM RAI #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          psem_mod_rai <- psem(
            glmmTMB(
              log_boars_rai ~ human_activity +
                forest_cover +
                crops +
                (1 | landscapeID),
              data = sim_data
            ),
            glmmTMB(
              log_herbivores_rai ~ log_boars_rai +
                human_activity +
                forest_cover +
                crops +
                (1 | landscapeID),
              data = sim_data
            ),
            glmmTMB(
              log_predators_rai ~ log_herbivores_rai +
                log_boars_rai +
                human_activity +
                forest_cover +
                (1 | landscapeID),
              data = sim_data
            ),
            sim_data
          )

          output[[length(output) + 1]] = summary(
            psem_mod_rai,
            standardize.type = "Mendard.OE",
            conserve = T
          )
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "piecewiseSEM_zi_med_rai") {
    ##### PieacewiseSEM ZI RAI (Jiang´s approach) #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          psem_mod_rai.jiang <- psem(
            glmmTMB(
              log_boars_rai ~ human_activity +
                forest_cover +
                crops +
                (1 | landscapeID),
              data = sim_data
            ),
            glmmTMB(
              log_herbivores_rai ~ log_boars_rai +
                boars_nonzero +
                human_activity +
                forest_cover +
                crops +
                (1 | landscapeID),
              data = sim_data
            ),
            glmmTMB(
              log_predators_rai ~ log_herbivores_rai +
                herbivores_nonzero +
                log_boars_rai +
                boars_nonzero +
                human_activity +
                forest_cover +
                (1 | landscapeID),
              data = sim_data
            ),
            sim_data
          )

          output[[length(output) + 1]] = summary(
            psem_mod_rai.jiang,
            standardize.type = "Mendard.OE",
            conserve = T
          )
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "lavaan_count") {
    ##### Lavaan Count #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          Lavaan_counts_mod <- "
          log_boars ~ human_activity + forest_cover + crops + effort
          log_herbivores ~ log_boars + human_activity + forest_cover +  crops + effort
          log_predators ~ log_herbivores + log_boars + human_activity + forest_cover + effort
        "

          Lavaan_counts_fit = sem(
            Lavaan_counts_mod,
            meanstructure = TRUE,
            cluster = "landscapeID",
            estimator = "MLM",
            data = sim_data
          )

          output[[length(output) + 1]] = summary(Lavaan_counts_fit)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "lavaan_zi_med_count") {
    ##### Lavaan Count #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          Lavaan_counts_mod_zi.jiang <- "
          log_boars ~ human_activity + forest_cover + crops + effort
          log_herbivores ~ log_boars + boars_nonzero + human_activity + forest_cover +  crops + effort
          log_predators ~ log_herbivores + herbivores_nonzero + log_boars + boars_nonzero + human_activity + forest_cover + effort
        "

          Lavaan_counts_fit_zi.jiang = sem(
            Lavaan_counts_mod_zi.jiang,
            meanstructure = TRUE,
            cluster = "landscapeID",
            estimator = "MLM",
            data = sim_data
          )

          output[[length(output) + 1]] = summary(Lavaan_counts_fit_zi.jiang)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "lavaan_rai") {
    ##### Lavaan RAI #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          Lavaan_rai_mod <- "
          log_boars_rai ~ human_activity + forest_cover + crops
          log_herbivores_rai ~ log_boars_rai + human_activity + forest_cover +  crops
          log_predators_rai ~ log_herbivores_rai + log_boars_rai + human_activity + forest_cover
        "

          Lavaan_rai_fit = sem(
            Lavaan_rai_mod,
            meanstructure = TRUE,
            cluster = "landscapeID",
            estimator = "MLM",
            data = sim_data
          )

          output[[length(output) + 1]] = summary(Lavaan_rai_fit)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "lavaan_zi_med_rai") {
    ##### Lavaan RAI (Jiang´s approach) #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          Lavaan_rai_zi.jiang_mod <- "
          log_boars_rai ~ human_activity + forest_cover + crops
          log_herbivores_rai ~ log_boars_rai + boars_nonzero + human_activity + forest_cover +  crops
          log_predators_rai ~ log_herbivores_rai + herbivores_nonzero + log_boars_rai + boars_nonzero + human_activity + forest_cover
        "

          Lavaan_rai_zi.jiang_fit = sem(
            Lavaan_rai_zi.jiang_mod,
            meanstructure = TRUE,
            cluster = "landscapeID",
            estimator = "MLM",
            data = sim_data
          )

          output[[length(output) + 1]] = summary(Lavaan_rai_zi.jiang_fit)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "SEMinR_count") {
    ##### SEMinR COUNT #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          SEMinR_count_measurements <- constructs(
            composite("boars", single_item("boars")),
            composite("herbivores", single_item("herbivores")),
            composite("predators", single_item("predators")),
            composite("human_activity", single_item("human_activity")),
            composite("forest_cover", single_item("forest_cover")),
            composite("crops", single_item("crops")),
            composite("effort", single_item("effort"))
          )

          SEMinR_count_struc_model = relationships(
            paths(from = "human_activity", to = "boars"),
            paths(from = "forest_cover", to = "boars"),
            paths(from = "crops", to = "boars"),
            paths(from = "effort", to = "boars"),

            paths(from = "boars", to = "herbivores"),
            paths(from = "human_activity", to = "herbivores"),
            paths(from = "forest_cover", to = "herbivores"),
            paths(from = "crops", to = "herbivores"),
            paths(from = "effort", to = "herbivores"),

            paths(from = "boars", to = "predators"),
            paths(from = "herbivores", to = "predators"),
            paths(from = "human_activity", to = "predators"),
            paths(from = "forest_cover", to = "predators"),
            paths(from = "effort", to = "predators")
          )

          SEMinR_count_fit = estimate_pls(
            data = sim_data,
            measurement_model = SEMinR_count_measurements,
            structural_model = SEMinR_count_struc_model,
            inner_weights = path_weighting
          )

          SEMinR_count_fit_bootstrap = bootstrap_model(
            SEMinR_count_fit,
            nboot = 10000
          )

          output[[length(output) + 1]] = summary(SEMinR_count_fit_bootstrap)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "SEMinR_rai") {
    ##### SEMinR RAI #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )

      tryCatch(
        {
          SEMinR_RAI_measurements <- constructs(
            composite("log_boars", single_item("log_boars_rai")),
            composite("log_herbivores", single_item("log_herbivores_rai")),
            composite("log_predators", single_item("log_predators_rai")),
            composite("human_activity", single_item("human_activity")),
            composite("forest_cover", single_item("forest_cover")),
            composite("crops", single_item("crops"))
          )

          SEMinR_RAI_struc_model = relationships(
            paths(from = "human_activity", to = "log_boars"),
            paths(from = "forest_cover", to = "log_boars"),
            paths(from = "crops", to = "log_boars"),

            paths(from = "log_boars", to = "log_herbivores"),
            paths(from = "human_activity", to = "log_herbivores"),
            paths(from = "forest_cover", to = "log_herbivores"),
            paths(from = "crops", to = "log_herbivores"),

            paths(from = "log_boars", to = "log_predators"),
            paths(from = "log_herbivores", to = "log_predators"),
            paths(from = "human_activity", to = "log_predators"),
            paths(from = "forest_cover", to = "log_predators")
          )

          SEMinR_RAI_fit = estimate_pls(
            data = sim_data,
            measurement_model = SEMinR_RAI_measurements,
            structural_model = SEMinR_RAI_struc_model,
            inner_weights = path_weighting
          )

          SEMinR_RAI_fit_bootstrap = bootstrap_model(
            SEMinR_RAI_fit,
            nboot = 10000
          )

          output[[length(output) + 1]] = summary(SEMinR_RAI_fit_bootstrap)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "because_zinb_count") {
    ##### Because ZINB Count #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$log_effort = log(sim_data$effort)

      tryCatch(
        {
          mod <- because(
            equations = list(
              boars ~ crops + human_activity + forest_cover + log_effort,
              herbivores ~ boars +
                human_activity +
                forest_cover +
                crops +
                log_effort,
              predators ~ herbivores +
                boars +
                human_activity +
                forest_cover +
                log_effort
            ),
            data = sim_data,
            distribution = list(
              boars = "zinb",
              herbivores = "zinb",
              predators = "zinb"
            ),
            random = ~ (1 | landscapeID), # Global random effect for landscape
            n.iter = 5000,
            n.burnin = 1000,
            n.chains = 3,
            parallel = TRUE,
            n.cores = 3,
            quiet = TRUE
          )

          output[[length(output) + 1]] = summary(mod)
        },
        error = function(e) {}
      )
    }
  } else if (model_list_log$approach[i] == "because_zi_med_count") {
    ##### Because ZI Mediation Count #####

    output = list()
    print(paste(m, "iteration 1", sep = " "))

    while (length(output) < Inter.num) {
      print(paste(m, "iteration", length(output) + 1, sep = " "))

      ### Create the data
      sim_data = create_count_sim_data(
        nlandscapes = model_list_log$nlandscapes[i],
        nsurveys = model_list_log$nsurveys[i],
        zero_inflation = model_list_log$zero_inflation[i]
      )
      sim_data$log_effort = log(sim_data$effort)
      sim_data$boars_nonzero = as.numeric(sim_data$boars > 0)
      sim_data$herbivores_nonzero = as.numeric(sim_data$herbivores > 0)

      tryCatch(
        {
          mod <- because(
            equations = list(
              boars ~ crops + human_activity + forest_cover + log_effort,
              herbivores ~ boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                crops +
                log_effort,
              predators ~ herbivores +
                herbivores_nonzero +
                boars +
                boars_nonzero +
                human_activity +
                forest_cover +
                log_effort
            ),
            data = sim_data,
            distribution = list(
              boars = "zinb",
              herbivores = "zinb",
              predators = "zinb"
            ),
            random = ~ (1 | landscapeID), # Global random effect for landscape
            n.iter = 5000,
            n.burnin = 1000,
            n.chains = 3,
            parallel = TRUE,
            n.cores = 3,
            quiet = TRUE
          )

          output[[length(output) + 1]] = summary(mod)
        },
        error = function(e) {}
      )
    }
  }

  ##### SAVE the output list

  saveRDS(output, file = paste(file_storade_wd, "/", m, ".rds", sep = ""))
  #readRDS(file = "my_data.rds")
}


##### Section 4 - Combine the model results #####

### List model files

setwd(file_storade_wd)

model_output = list.files()
model_output = model_output[str_detect(model_output, ".rds")]
model_output


### Make template table

results = as.data.frame(matrix(NA, ncol = 7, nrow = 0))
colnames(results) = c(
  "model_ID",
  "target_parameter",
  "approach",
  "path",
  "estimate",
  "SE",
  "significance"
)


model_significance = as.data.frame(matrix(NA, ncol = 2, nrow = 0))
colnames(model_significance) = c("model_ID", "significance")

unique(model_list_log$approach)


#Make a custom function
path_flip = function(x) {
  return(paste(
    str_split(x, "  ->  ")[[1]][2],
    str_split(x, "  ->  ")[[1]][1],
    sep = "_"
  ))
}


### Extract the results

for (it in 1:length(model_output)) {
  i = model_output[it]

  #Load the file
  dt = read_rds(i)

  model_ID = str_replace(i, ".rds", "")
  target_parameter = model_list_log$target_parameter[
    model_list_log$model_ID == model_ID
  ]
  approach = model_list_log$approach[model_list_log$model_ID == model_ID]

  results_temp = results[0, ]

  ###brms_count
  if (str_detect(i, "brms_count")) {
    #select the model type

    dt = dt[map_int(dt, length) == 18] #the summary of this model should contain 18 items. if it fails during estimating, it returns 8 items
    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$fixed #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = rownames(t),
          "SE" = Est.Error,
          "significance" = case_when(
            `l-95% CI` > 0 & `u-95% CI` > 0 ~ TRUE,
            `l-95% CI` < 0 & `u-95% CI` < 0 ~ TRUE,
            is.na(`l-95% CI`) | is.na(`u-95% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"brms_zi_count"
  } else if (str_detect(i, "brms_zi_count")) {
    #select the model type

    dt = dt[map_int(dt, length) == 18] #the summary of this model should contain 18 items. if it fails during estimating, it returns 8 items
    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$fixed #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = rownames(t),
          "SE" = Est.Error,
          "significance" = case_when(
            `l-95% CI` > 0 & `u-95% CI` > 0 ~ TRUE,
            `l-95% CI` < 0 & `u-95% CI` < 0 ~ TRUE,
            is.na(`l-95% CI`) | is.na(`u-95% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"brms_zi_med_count"
  } else if (str_detect(i, "brms_zi_med_count")) {
    #select the model type

    dt = dt[map_int(dt, length) == 18] #the summary of this model should contain 18 items. if it fails during estimating, it returns 8 items
    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$fixed #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = rownames(t),
          "SE" = Est.Error,
          "significance" = case_when(
            `l-95% CI` > 0 & `u-95% CI` > 0 ~ TRUE,
            `l-95% CI` < 0 & `u-95% CI` < 0 ~ TRUE,
            is.na(`l-95% CI`) | is.na(`u-95% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"brms_rai"
  } else if (str_detect(i, "brms_rai")) {
    #select the model type

    dt = dt[map_int(dt, length) == 18] #the summary of this model should contain 18 items. if it fails during estimating, it returns 8 items
    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$fixed #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = rownames(t),
          "SE" = Est.Error,
          "significance" = case_when(
            `l-95% CI` > 0 & `u-95% CI` > 0 ~ TRUE,
            `l-95% CI` < 0 & `u-95% CI` < 0 ~ TRUE,
            is.na(`l-95% CI`) | is.na(`u-95% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results))

      t$path = str_replace(t$path, "logboarsrai_crops", "boars_crops")
      t$path = str_replace(
        t$path,
        "logherbivoresrai_log_boars_rai",
        "herbivores_boars"
      )
      t$path = str_replace(
        t$path,
        "logherbivoresrai_human_activity",
        "herbivores_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logpredatorsrai_human_activity",
        "predators_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logboarsrai_human_activity",
        "boars_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logpredatorsrai_log_herbivores_rai",
        "predators_herbivores"
      )

      t = t %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"brms_zi_med_rai"
  } else if (str_detect(i, "brms_zi_med_rai")) {
    #select the model type

    dt = dt[map_int(dt, length) == 18] #the summary of this model should contain 18 items. if it fails during estimating, it returns 8 items
    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$fixed #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = rownames(t),
          "SE" = Est.Error,
          "significance" = case_when(
            `l-95% CI` > 0 & `u-95% CI` > 0 ~ TRUE,
            `l-95% CI` < 0 & `u-95% CI` < 0 ~ TRUE,
            is.na(`l-95% CI`) | is.na(`u-95% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results))

      t$path = str_replace(t$path, "logboarsrai_crops", "boars_crops")
      t$path = str_replace(
        t$path,
        "logherbivoresrai_log_boars_rai",
        "herbivores_boars"
      )
      t$path = str_replace(
        t$path,
        "logherbivoresrai_human_activity",
        "herbivores_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logpredatorsrai_human_activity",
        "predators_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logboarsrai_human_activity",
        "boars_human_activity"
      )
      t$path = str_replace(
        t$path,
        "logpredatorsrai_log_herbivores_rai",
        "predators_herbivores"
      )

      t = t %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"piecewiseSEM_count"
  } else if (str_detect(i, "piecewiseSEM_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$coefficients #extract the estimate table

      t = t[, !colnames(t) == ""]

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = paste(t$Response, t$Predictor, sep = "_"),
          "SE" = Std.Error,
          "significance" = case_when(
            P.Value <= 0.05 ~ TRUE,
            is.na(P.Value) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      t2 = data.frame(
        "model_ID" = model_ID,
        "significance" = dt[b][[1]]$Cstat$P.Value
      )
      model_significance = rbind(model_significance, t2)
    }

    ###"piecewiseSEM_zi_count"
  } else if (str_detect(i, "piecewiseSEM_zi_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$coefficients #extract the estimate table

      t = t[, !colnames(t) == ""]

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = paste(t$Response, t$Predictor, sep = "_"),
          "SE" = Std.Error,
          "significance" = case_when(
            P.Value <= 0.05 ~ TRUE,
            is.na(P.Value) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      t2 = data.frame(
        "model_ID" = model_ID,
        "significance" = dt[b][[1]]$Cstat$P.Value
      )
      model_significance = rbind(model_significance, t2)
    }

    ###"piecewiseSEM_zi_med_count"
  } else if (str_detect(i, "piecewiseSEM_zi_med_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$coefficients #extract the estimate table

      t = t[, !colnames(t) == ""]

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          "path" = paste(t$Response, t$Predictor, sep = "_"),
          "SE" = Std.Error,
          "significance" = case_when(
            P.Value <= 0.05 ~ TRUE,
            is.na(P.Value) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      t2 = data.frame(
        "model_ID" = model_ID,
        "significance" = dt[b][[1]]$Cstat$P.Value
      )
      model_significance = rbind(model_significance, t2)
    }

    ###"piecewiseSEM_rai"
  } else if (str_detect(i, "piecewiseSEM_rai")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$coefficients #extract the estimate table

      t = t[, !colnames(t) == ""]

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          Response = str_replace(Response, "log_", ""),
          Response = str_replace(Response, "_rai", ""),
          Predictor = str_replace(Predictor, "log_", ""),
          Predictor = str_replace(Predictor, "_rai", ""),
          "path" = paste(Response, Predictor, sep = "_"),
          "SE" = Std.Error,
          "significance" = case_when(
            P.Value <= 0.05 ~ TRUE,
            is.na(P.Value) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      t2 = data.frame(
        "model_ID" = model_ID,
        "significance" = dt[b][[1]]$Cstat$P.Value
      )
      model_significance = rbind(model_significance, t2)
    }

    ###"piecewiseSEM_zi_med_rai"
  } else if (str_detect(i, "piecewiseSEM_zi_med_rai")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$coefficients #extract the estimate table

      t = t[, !colnames(t) == ""]

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Estimate,
          Response = str_replace(Response, "log_", ""),
          Response = str_replace(Response, "_rai", ""),
          Predictor = str_replace(Predictor, "log_", ""),
          Predictor = str_replace(Predictor, "_rai", ""),
          "path" = paste(Response, Predictor, sep = "_"),
          "SE" = Std.Error,
          "significance" = case_when(
            P.Value <= 0.05 ~ TRUE,
            is.na(P.Value) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      t2 = data.frame(
        "model_ID" = model_ID,
        "significance" = dt[b][[1]]$Cstat$P.Value
      )
      model_significance = rbind(model_significance, t2)
    }

    ###"lavaan_count"
  } else if (str_detect(i, "lavaan_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$pe #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = est,
          rhs = str_replace(rhs, "log_", ""),
          lhs = str_replace(lhs, "log_", ""),
          "path" = paste(lhs, rhs, sep = "_"),
          "SE" = se,
          "significance" = case_when(
            pvalue <= 0.05 ~ TRUE,
            is.na(pvalue) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      pvalue = dt[b][[1]]$test$satorra.bentler$pvalue
      if (!is.null(pvalue)) {
        t2 = data.frame(
          "model_ID" = model_ID,
          "significance" = dt[b][[1]]$test$satorra.bentler$pvalue
        )
        model_significance = rbind(model_significance, t2)
      }
    }

    ###"lavaan_zi_med_count"
  } else if (str_detect(i, "lavaan_zi_med_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$pe #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = est,
          rhs = str_replace(rhs, "log_", ""),
          lhs = str_replace(lhs, "log_", ""),
          "path" = paste(lhs, rhs, sep = "_"),
          "SE" = se,
          "significance" = case_when(
            pvalue <= 0.05 ~ TRUE,
            is.na(pvalue) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      pvalue = dt[b][[1]]$test$satorra.bentler$pvalue
      if (!is.null(pvalue)) {
        t2 = data.frame(
          "model_ID" = model_ID,
          "significance" = dt[b][[1]]$test$satorra.bentler$pvalue
        )
        model_significance = rbind(model_significance, t2)
      }
    }

    ###"lavaan_rai"
  } else if (str_detect(i, "lavaan_rai")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$pe #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = est,
          rhs = str_replace(rhs, "log_", ""),
          lhs = str_replace(lhs, "log_", ""),
          rhs = str_replace(rhs, "_rai", ""),
          lhs = str_replace(lhs, "_rai", ""),
          "path" = paste(lhs, rhs, sep = "_"),
          "SE" = se,
          "significance" = case_when(
            pvalue <= 0.05 ~ TRUE,
            is.na(pvalue) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      pvalue = dt[b][[1]]$test$satorra.bentler$pvalue
      if (!is.null(pvalue)) {
        t2 = data.frame(
          "model_ID" = model_ID,
          "significance" = dt[b][[1]]$test$satorra.bentler$pvalue
        )
        model_significance = rbind(model_significance, t2)
      }
    }

    ###"lavaan_zi_med_rai"
  } else if (str_detect(i, "lavaan_zi_med_rai")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = dt[b][[1]]$pe #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = est,
          rhs = str_replace(rhs, "log_", ""),
          lhs = str_replace(lhs, "log_", ""),
          rhs = str_replace(rhs, "_rai", ""),
          lhs = str_replace(lhs, "_rai", ""),
          "path" = paste(lhs, rhs, sep = "_"),
          "SE" = se,
          "significance" = case_when(
            pvalue <= 0.05 ~ TRUE,
            is.na(pvalue) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }
      results_temp = rbind(results_temp, t)

      pvalue = dt[b][[1]]$test$satorra.bentler$pvalue
      if (!is.null(pvalue)) {
        t2 = data.frame(
          "model_ID" = model_ID,
          "significance" = dt[b][[1]]$test$satorra.bentler$pvalue
        )
        model_significance = rbind(model_significance, t2)
      }
    }

    ###"SEMinR_count"
  } else if (str_detect(i, "SEMinR_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = as.data.frame(dt[b][[1]]$bootstrapped_paths) #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = `Bootstrap Mean`,
          "path" = rownames(t),
          "path" = map_chr(path, path_flip),
          "SE" = `Bootstrap SD`,
          "significance" = case_when(
            `2.5% CI` > 0 & `97.5% CI` > 0 ~ TRUE,
            `2.5% CI` < 0 & `97.5% CI` < 0 ~ TRUE,
            is.na(`2.5% CI`) | is.na(`97.5% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"SEMinR_count"
  } else if (str_detect(i, "SEMinR_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = as.data.frame(dt[b][[1]]$bootstrapped_paths) #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = `Bootstrap Mean`,
          "path" = rownames(t),
          "path" = map_chr(path, path_flip),
          "SE" = `Bootstrap SD`,
          "significance" = case_when(
            `2.5% CI` > 0 & `97.5% CI` > 0 ~ TRUE,
            `2.5% CI` < 0 & `97.5% CI` < 0 ~ TRUE,
            is.na(`2.5% CI`) | is.na(`97.5% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }

    ###"SEMinR_rai"
  } else if (str_detect(i, "SEMinR_rai")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = as.data.frame(dt[b][[1]]$bootstrapped_paths) #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = `Bootstrap Mean`,
          "path" = rownames(t),
          "path" = map_chr(path, path_flip),
          path = str_replace_all(path, "log_", ""),
          "SE" = `Bootstrap SD`,
          "significance" = case_when(
            `2.5% CI` > 0 & `97.5% CI` > 0 ~ TRUE,
            `2.5% CI` < 0 & `97.5% CI` < 0 ~ TRUE,
            is.na(`2.5% CI`) | is.na(`97.5% CI`) ~ NA,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }
  } else if (str_detect(i, "because_zinb_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = as.data.frame(dt[[b]]$results) #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Mean,
          "path" = rownames(t),
          "path" = str_replace(path, "beta_", ""), # remove beta_ prefix
          "SE" = SD,
          "significance" = case_when(
            `2.5%` > 0 & `97.5%` > 0 ~ TRUE,
            `2.5%` < 0 & `97.5%` < 0 ~ TRUE,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }
  } else if (str_detect(i, "because_zi_med_count")) {
    #select the model type

    for (b in 1:length(dt)) {
      # Loop per individual model summaries

      t = as.data.frame(dt[[b]]$results) #extract the estimate table

      t = t %>% #Turn it into our desirable format
        mutate(
          "model_ID" = model_ID,
          "target_parameter" = target_parameter,
          "approach" = approach,
          "estimate" = Mean,
          "path" = rownames(t),
          "path" = str_replace(path, "beta_", ""), # remove beta_ prefix
          "SE" = SD,
          "significance" = case_when(
            `2.5%` > 0 & `97.5%` > 0 ~ TRUE,
            `2.5%` < 0 & `97.5%` < 0 ~ TRUE,
            .default = FALSE
          )
        ) %>%
        select(colnames(results)) %>%
        filter(
          path %in%
            c(
              "boars_crops",
              "herbivores_boars",
              "herbivores_human_activity",
              "predators_human_activity",
              "boars_human_activity",
              "predators_herbivores"
            )
        ) #Select paths of interest

      if (!nrow(t) == 6) {
        t = t[0, ]
      }

      rownames(t) = NULL
      results_temp = rbind(results_temp, t)
    }
  }

  results = rbind(results, results_temp)

  print(paste(it, " of ", length(model_output), sep = ""))
}


### Add the true value of simulated parameters
results$true_value = 0
results$true_value[results$path == "boars_crops"] = 0.5
results$true_value[results$path == "herbivores_boars"] = 0.0005
results$true_value[results$path == "herbivores_human_activity"] = -0.47
results$true_value[results$path == "predators_human_activity"] = -0.19


### calculate the true error

results = results %>%
  mutate(
    estimate = as.numeric(estimate),
    SE = as.numeric(SE),
    "true_error" = estimate - true_value,
    "approach_path" = paste(approach, path, sep = " - ")
  )


setwd(wd)


#-----------------------------------------------------------------------------------------------------------------------------------------#
##### Make the Plots #####
#-----------------------------------------------------------------------------------------------------------------------------------------#

#Rename the models
plot_data = results %>%
  mutate("approach_path_plot_names" = path)

plot_data$approach_path_plot_names = plot_data$approach_path_plot_names %>%
  str_replace("boars_crops", "Boars ~ crops") %>%
  str_replace("predators_herbivores", "Predators ~ herbivores") %>%
  str_replace("herbivores_human_activity", "Herbivores ~ humans") %>%
  str_replace("boars_human_activity", "Boars ~ humans") %>%
  str_replace("herbivores_boars", "Herbivores ~ boars") %>%
  str_replace("predators_human_activity", "Predators ~ humans")


plot_data$approach_path = factor(
  plot_data$approach_path,
  levels = sort(unique(plot_data$approach_path), decreasing = T)
)


labels = rep(
  c(
    "Predators ~ humans",
    "Predators ~ herbivores",
    "Herbivores ~ humans",
    "Herbivores ~ boars",
    "Boars ~ humans",
    "Boars ~ crops"
  ),
  16
)


#### Select the simulated paths
plot_data_paths = plot_data[
  !plot_data$approach_path_plot_names %in%
    c("Boars ~ humans", "Predators ~ herbivores"),
]
label_paths = rep(
  c(
    "Predators ~ humans",
    "Herbivores ~ humans",
    "Herbivores ~ boars",
    "Boars ~ crops"
  ),
  16
)

#### Select paths which were not directly simulated (coefficient should be 0)
plot_data_no_paths = plot_data[
  plot_data$approach_path_plot_names %in%
    c("Boars ~ humans", "Predators ~ herbivores"),
]
label_no_paths = rep(
  c(
    "Predators ~ humans",
    "Predators ~ herbivores",
    "Herbivores ~ humans",
    "Herbivores ~ boars",
    "Boars ~ humans",
    "Boars ~ crops"
  ),
  16
)


#ZI0
plot_data_paths_zi0 = plot_data_paths[
  plot_data_paths$target_parameter == "zero_inflation0",
]
ggplot(
  plot_data_paths_zi0,
  aes(x = approach_path, y = true_error, fill = approach)
) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(-2.5, 2.5) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width") +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Estimates_zi0.jpeg", width = 20, height = 30, units = "cm", scale = 1)

#ZI0 Errors
plot_data_paths_zi0 = plot_data_paths[
  plot_data_paths$target_parameter == "zero_inflation0",
]
ggplot(plot_data_paths_zi0, aes(x = approach_path, y = SE, fill = approach)) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(0, 3.2) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width") +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Errors_zi0.jpeg", width = 20, height = 30, units = "cm", scale = 1)

#ZI0.3
plot_data_paths_zi0.3 = plot_data_paths[
  plot_data_paths$target_parameter == "zero_inflation0.3",
]
ggplot(
  plot_data_paths_zi0.3,
  aes(x = approach_path, y = true_error, fill = approach)
) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(-2.5, 2.5) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width") +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Estimates_zi0.3.jpeg", width = 20, height = 30, units = "cm", scale = 1)

#ZI0.3 Errors
plot_data_paths_zi0.3 = plot_data_paths[
  plot_data_paths$target_parameter == "zero_inflation0.3",
]
ggplot(plot_data_paths_zi0.3, aes(x = approach_path, y = SE, fill = approach)) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(0, 3.2) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width") +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Errors_zi0.3.jpeg", width = 20, height = 30, units = "cm", scale = 1)

#ZI0.5
plot_data_paths_zi0.5 = plot_data_paths[
  plot_data_paths$target_parameter == "zero_inflation0.5",
]
ggplot(
  plot_data_paths_zi0.5,
  aes(x = approach_path, y = true_error, fill = approach)
) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(-2.5, 2.5) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width", ) +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Estimates_zi0.5.jpeg", width = 20, height = 30, units = "cm", scale = 1)

ggplot(plot_data_paths_zi0.5, aes(x = approach_path, y = SE, fill = approach)) +
  geom_hline(yintercept = 0, color = "grey20") +
  ylim(0, 3.2) +
  scale_fill_manual(
    values = c(
      "grey40",
      "white",
      "dodgerblue4",
      "firebrick3",
      "forestgreen",
      "pink",
      "yellow",
      "lightskyblue",
      "darkorange4",
      "darkviolet",
      "springgreen2",
      "navajowhite",
      "cyan2",
      "darkorange",
      "violet",
      "grey80"
    )
  ) +

  geom_violin(scale = "width") +
  scale_x_discrete(labels = label_paths) +
  stat_summary(fun = mean, geom = "point", shape = 3, size = 2) +
  coord_flip() +
  theme_bw()

#ggsave("Errors_zi0.5.jpeg", width = 20, height = 30, units = "cm", scale = 1)

##### MAKE A TABLE WITH THE MAIN RESULTS

make_result_table = function(x, y, zi) {
  # plot_data_zi0, model_significance_zi0, "0%"

  result_table = as.data.frame(matrix(
    NA,
    ncol = 14,
    nrow = length(unique(x$approach))
  ))
  colnames(result_table) = c(
    "approach",
    "zero_inflation",
    "deviation_mean",
    "deviation_median",
    "deviation_sd",
    "error_mean",
    "error_median",
    "error_sd",
    "model_performace",
    "CropPigs",
    "Pigsherb",
    "HumanHerb",
    "HumanCarn",
    "false_paths"
  )
  result_table$approach = unique(x$approach)

  result_table$zero_inflation = zi

  x_no_paths = x[x$true_value == 0, ] #not simulated paths (i.e., should have coefficient zero)
  x = x[x$true_value != 0, ] #Simulated paths

  colnames(x)
  ###FILL the table

  for (i in as.character(unique(x$approach))) {
    result_table$deviation_mean[
      result_table$approach == i
    ] = mean(abs(x$true_error[x$approach == i]))
    result_table$deviation_median[
      result_table$approach == i
    ] = median(abs(x$true_error[x$approach == i]))
    result_table$deviation_sd[result_table$approach == i] = sd(x$true_error[
      x$approach == i
    ])

    result_table$error_mean[result_table$approach == i] = mean(
      x$SE[x$approach == i],
      na.rm = T
    )
    result_table$error_median[result_table$approach == i] = median(
      x$SE[x$approach == i],
      na.rm = T
    )
    result_table$error_sd[result_table$approach == i] = sd(
      x$SE[x$approach == i],
      na.rm = T
    )

    result_table$CropPigs[result_table$approach == i] = 100 *
      sum(
        x$significance[x$path == "boars_crops" & x$approach == i],
        na.rm = T
      ) /
      length(na.omit(x$significance[x$path == "boars_crops" & x$approach == i]))
    result_table$Pigsherb[result_table$approach == i] = 100 *
      sum(
        x$significance[x$path == "herbivores_boars" & x$approach == i],
        na.rm = T
      ) /
      length(na.omit(x$significance[
        x$path == "herbivores_boars" & x$approach == i
      ]))
    result_table$HumanHerb[result_table$approach == i] = 100 *
      sum(
        x$significance[x$path == "herbivores_human_activity" & x$approach == i],
        na.rm = T
      ) /
      length(na.omit(x$significance[
        x$path == "herbivores_human_activity" & x$approach == i
      ]))
    result_table$HumanCarn[result_table$approach == i] = 100 *
      sum(
        x$significance[x$path == "predators_human_activity" & x$approach == i],
        na.rm = T
      ) /
      length(na.omit(x$significance[
        x$path == "predators_human_activity" & x$approach == i
      ]))

    result_table$false_paths[result_table$approach == i] = 100 *
      sum(
        x_no_paths$significance[
          x_no_paths$path %in%
            c("boars_human_activity", "predators_herbivores") &
            x_no_paths$approach == i
        ],
        na.rm = T
      ) /
      length(na.omit(x_no_paths$significance[
        x_no_paths$path %in%
          c("boars_human_activity", "predators_herbivores") &
          x_no_paths$approach == i
      ]))

    if (str_detect(i, "lavaan")) {
      result_table$model_performace[result_table$approach == i] = 100 *
        sum(y[str_detect(y$model_ID, i), ]$significance <= 0.05, na.rm = T) /
        length(na.omit(y[str_detect(y$model_ID, i), ]$significance <= 0.05))
    }
    if (str_detect(i, "piecewiseSEM")) {
      result_table$model_performace[result_table$approach == i] = 100 *
        sum(y[str_detect(y$model_ID, i), ]$significance > 0.05, na.rm = T) /
        length(na.omit(y[str_detect(y$model_ID, i), ]$significance <= 0.05))
    }
  }

  return(result_table)
}


plot_data_zi0 = plot_data[plot_data$target_parameter == "zero_inflation0", ]
model_significance_zi0 = model_significance[
  str_detect(model_significance$model_ID, "zero_inflation0-"),
]

result_table_zi0 = make_result_table(
  plot_data_zi0,
  model_significance_zi0,
  "0%"
)


plot_data_zi0.3 = plot_data[plot_data$target_parameter == "zero_inflation0.3", ]
model_significance_zi0.3 = model_significance[
  str_detect(model_significance$model_ID, "zero_inflation0.3"),
]

result_table_zi0.3 = make_result_table(
  plot_data_zi0.3,
  model_significance_zi0.3,
  "30%"
)


plot_data_zi0.5 = plot_data[plot_data$target_parameter == "zero_inflation0.5", ]
model_significance_zi0.5 = model_significance[
  str_detect(model_significance$model_ID, "zero_inflation0.5"),
]

result_table_zi0.5 = make_result_table(
  plot_data_zi0.5,
  model_significance_zi0.5,
  "50%"
)


###Merge

result_table = rbind(result_table_zi0, result_table_zi0.3, result_table_zi0.5)


result_table = result_table = result_table[
  order(result_table$approach, result_table$zero_inflation),
]


### Make some rounding and select some columns

result_table = result_table %>%
  mutate(
    deviation_mean = round(deviation_mean, 3),
    error_mean = round(error_mean, 3),
    CropPigs = round(CropPigs, 1),
    Pigsherb = round(Pigsherb, 1),
    HumanHerb = round(HumanHerb, 1),
    HumanCarn = round(HumanCarn, 1),
    false_paths = round(false_paths, 1),
    upCI = round(deviation_mean + (1.96 * error_mean), 3),
    loCI = round(deviation_mean - (1.96 * error_mean), 3)
  ) %>%
  select(c(
    approach,
    zero_inflation,
    deviation_mean,
    error_mean,
    upCI,
    loCI,
    model_performace,
    CropPigs,
    Pigsherb,
    HumanHerb,
    HumanCarn,
    false_paths
  ))


#write.csv(result_table, "Result table.csv")

##### MAKE THE LINE PLOTS

###Summarize the data

dt = plot_data_paths %>%
  group_by(model_ID) %>%
  nest() %>%
  mutate(
    "dev.from.true.coef" = map_dbl(
      .x = data,
      .f = ~ mean(abs(.x$true_error), na.rm = T)
    ),
    "stat.power" = map_dbl(.x = data, .f = ~ mean(.x$significance, na.rm = T)),
    "estimated.error" = map_dbl(.x = data, .f = ~ mean(.x$SE, na.rm = T)),
    "target_parameter" = map_chr(
      .x = data,
      .f = ~ unique(.x$target_parameter, na.rm = T)
    ),
    "approach" = map_chr(.x = data, .f = ~ unique(.x$approach, na.rm = T)),
    "package" = approach,
    "approach2" = approach
  ) %>%
  select(-data)


#Set the lines
dt$package[str_detect(dt$package, "brms")] = "brms"
dt$package[str_detect(dt$package, "lavaan")] = "lavaan"
dt$package[str_detect(dt$package, "piecewiseSEM")] = "piecewiseSEM"
dt$package[str_detect(dt$package, "SEMinR")] = "SEMinR"


dt$approach2[str_detect(dt$approach, "_count")] = "counts"
dt$approach2[str_detect(dt$approach, "_rai")] = "RAI"
dt$approach2[str_detect(dt$approach, "_zi_count")] = "ZI_counts"
dt$approach2[str_detect(dt$approach, "_zi_J_count")] = "ZI.med_counts"
dt$approach2[str_detect(dt$approach, "_zi_count_jiang")] = "ZI.med_counts"
dt$approach2[str_detect(dt$approach, "_zi_rai_jiang")] = "ZI.med_RAI"
dt$approach2[str_detect(dt$approach, "_zi_J_rai")] = "ZI.med_RAI"


#Zero-Inflation

dt_zi = dt %>%
  filter(str_detect(model_ID, "zero_inflation")) %>%
  mutate("zi" = 0)

dt_zi$zi[str_detect(dt_zi$model_ID, "0.3")] = 30
dt_zi$zi[str_detect(dt_zi$model_ID, "0.5")] = 50


ggplot(
  dt_zi,
  aes(
    x = zi,
    y = dev.from.true.coef,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Deviation from the True Coefficient") +
  xlab("Zero Inflation (%)") +
  scale_x_continuous(breaks = c(0, 30, 50)) +
  theme_light()

#ggsave("Deviation from True Coefficient and Zero-Inflation.jpeg", width = 14, height = 10, units = "cm", scale = 1)

ggplot(
  dt_zi,
  aes(
    x = zi,
    y = stat.power,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Statistical Power") +
  xlab("Zero Inflation (%)") +
  scale_x_continuous(breaks = c(0, 30, 50)) +
  theme_light()

#ggsave("Statistical Power and Zero-Inflation.jpeg", width = 14, height = 10, units = "cm", scale = 1)

# Sample Size

dt_ns = dt %>%
  filter(str_detect(model_ID, "nsurveys"))

dt_ns$ns = as.numeric(str_replace(dt_ns$target_parameter, "nsurveys", ""))


ggplot(
  dt_ns,
  aes(
    x = ns,
    y = dev.from.true.coef,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Deviation from the True Coefficient") +
  xlab("Sample size") +
  scale_x_continuous(breaks = sort(unique(dt_ns$ns))) +
  theme_light()

#ggsave("Deviation from True Coefficient and Sample Size.jpeg", width = 14, height = 10, units = "cm", scale = 1)

ggplot(
  dt_ns,
  aes(
    x = ns,
    y = stat.power,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Statistical Power") +
  xlab("Sample size") +
  scale_x_continuous(breaks = sort(unique(dt_ns$ns))) +
  theme_light()

#ggsave("Statistical Power and Sample Size.jpeg", width = 14, height = 10, units = "cm", scale = 1)

# Number of landscapes (Random variables)

dt_nl = dt %>%
  filter(str_detect(model_ID, "nlandscapes"))

dt_nl$nl = as.numeric(str_replace(dt_nl$target_parameter, "nlandscapes", ""))


ggplot(
  dt_nl,
  aes(
    x = nl,
    y = dev.from.true.coef,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Deviation from the True Coefficient") +
  xlab("Number of Landscapes (random variable levels)") +
  scale_x_continuous(breaks = sort(unique(dt_nl$nl))) +
  theme_light()

#ggsave("Deviation from True Coefficient and Number of Random Variables.jpeg", width = 14, height = 10, units = "cm", scale = 1)

ggplot(
  dt_nl,
  aes(
    x = nl,
    y = stat.power,
    group = approach,
    col = package,
    linetype = approach2,
    shape = approach2
  )
) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    values = c("tomato", "olivedrab3", "cornflowerblue", "goldenrod")
  ) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15)) +
  scale_linetype_manual(
    values = c("solid", "dashed", "solid", "dashed", "solid")
  ) +
  ylab("Statistical Power") +
  xlab("Number of Landscapes (random variable levels)") +
  scale_x_continuous(breaks = sort(unique(dt_nl$nl))) +
  theme_light()

#ggsave("Statistical Power and Number of Random Variables.jpeg", width = 14, height = 10, units = "cm", scale = 1)

### Kolmogorov-Smirnov Tests (not exhaustive)
unique(plot_data_paths_zi0$approach)


ks.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_count"
  ]
)

ks.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_count"
  ]
)

ks.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"]
)


### Anderson-Darling k-Sample Test (It is just a multi-variable alternative for the Kolmogorov-Smirnov Test)

ad.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_count"
  ]
)


ks.test(
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "brms_zi_med_rai"
  ]
)


ad.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_count"
  ],
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "brms_zi_med_rai"
  ]
)


ks.test(
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_rai"
  ]
)

ks.test(
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "brms_zi_med_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_rai"
  ]
)


ad.test(
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "brms_zi_med_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_rai"
  ]
)


ad.test(
  plot_data_paths_zi0$true_error[plot_data_paths_zi0$approach == "brms_rai"],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "brms_zi_med_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_rai"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_med_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_zi_count"
  ],
  plot_data_paths_zi0$true_error[
    plot_data_paths_zi0$approach == "piecewiseSEM_count"
  ]
)


ks.test(
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "lavaan_zi_med_rai"
  ],
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "brms_zi_med_rai"
  ]
)

ks.test(
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "lavaan_zi_med_rai"
  ],
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "piecewiseSEM_zi_med_rai"
  ]
)

ad.test(
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "lavaan_zi_med_rai"
  ],
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "brms_zi_med_rai"
  ],
  plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "piecewiseSEM_zi_med_rai"
  ]
)


plot_data_paths_ns1200 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys1200",
]
plot_data_paths_ns1000 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys1000",
]
plot_data_paths_ns800 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys800",
]
plot_data_paths_ns600 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys600",
]
plot_data_paths_ns400 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys400",
]
plot_data_paths_ns200 = plot_data_paths[
  plot_data_paths$target_parameter == "nsurveys200",
]


ks.test(
  abs(plot_data_paths_ns1200$true_error[
    plot_data_paths_ns1200$approach == "piecewiseSEM_zi_med_count"
  ]),
  abs(plot_data_paths_ns200$true_error[
    plot_data_paths_ns200$approach == "piecewiseSEM_zi_med_count"
  ])
)


ad.test(
  abs(plot_data_paths_ns1200$true_error[
    plot_data_paths_ns1200$approach == "SEMinR_rai"
  ]),
  abs(plot_data_paths_ns1000$true_error[
    plot_data_paths_ns1000$approach == "SEMinR_rai"
  ]),
  abs(plot_data_paths_ns800$true_error[
    plot_data_paths_ns800$approach == "SEMinR_rai"
  ]),
  abs(plot_data_paths_ns600$true_error[
    plot_data_paths_ns600$approach == "SEMinR_rai"
  ]),
  abs(plot_data_paths_ns400$true_error[
    plot_data_paths_ns400$approach == "SEMinR_rai"
  ]),
  abs(plot_data_paths_ns200$true_error[
    plot_data_paths_ns200$approach == "SEMinR_rai"
  ])
)


ad.test(
  abs(plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "lavaan_zi_med_rai"
  ]),
  abs(plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "brms_zi_med_rai"
  ]),
  abs(plot_data_paths_zi0.5$true_error[
    plot_data_paths_zi0.5$approach == "piecewiseSEM_zi_med_rai"
  ])
)


plot_data_paths_nl10 = plot_data_paths[
  plot_data_paths$target_parameter == "nlandscapes10",
]
plot_data_paths_nl20 = plot_data_paths[
  plot_data_paths$target_parameter == "nlandscapes20",
]
plot_data_paths_nl40 = plot_data_paths[
  plot_data_paths$target_parameter == "nlandscapes40",
]
plot_data_paths_nl80 = plot_data_paths[
  plot_data_paths$target_parameter == "nlandscapes80",
]
plot_data_paths_nl120 = plot_data_paths[
  plot_data_paths$target_parameter == "nlandscapes120",
]


ks.test(
  abs(plot_data_paths_nl80$true_error[
    plot_data_paths_nl80$approach == "lavaan_zi_med_rai"
  ]),
  abs(plot_data_paths_nl10$true_error[
    plot_data_paths_nl10$approach == "lavaan_zi_med_rai"
  ])
)

ks.test(
  abs(plot_data_paths_nl120$true_error[
    plot_data_paths_nl120$approach == "lavaan_zi_med_count"
  ]),
  abs(plot_data_paths_nl10$true_error[
    plot_data_paths_nl10$approach == "lavaan_zi_med_count"
  ])
)

