phybase_model <- function(equations) {

  # Track beta names to ensure uniqueness
  beta_counter <- list()
  response_counter <- list()

  # Parse equations
  eq_list <- lapply(equations, function(eq) {
    response <- as.character(formula(eq))[2]
    predictors <- attr(terms(formula(eq)), "term.labels")
    list(response = response, predictors = predictors)
  })

  # Start model
  model_lines <- c("model {", "  # Structural equations", "  for (i in 1:N) {")

  for (j in seq_along(eq_list)) {
    eq <- eq_list[[j]]
    response <- eq$response
    predictors <- eq$predictors

    # Count and assign unique suffix for the response variable
    response_count <- response_counter[[response]] %||% 0
    response_count <- response_count + 1
    response_counter[[response]] <- response_count
    suffix <- if (response_count == 1) "" else as.character(response_count)

    mu <- paste0("mu", response, suffix)
    alpha <- paste0("alpha", response, suffix)
    response_var <- paste0(response, suffix)

    linpred <- alpha
    for (pred in predictors) {
      key <- paste(response, pred, suffix, sep = "_")
      if (!key %in% names(beta_counter)) {
        base <- paste0("beta", pred)
        count <- sum(grepl(paste0("^", base), unlist(beta_counter)))
        beta_name <- if (count == 0) base else paste0(base, count + 1)
        beta_counter[[key]] <- beta_name
      }
      beta_name <- beta_counter[[key]]
      linpred <- paste0(linpred, " + ", beta_name, "*", pred, "[i]")
    }
    model_lines <- c(model_lines, paste0("    ", mu, "[i] <- ", linpred))
  }

  model_lines <- c(model_lines, "  }", "  # Multivariate normal likelihoods")

  # Add likelihoods
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      response_var <- paste0(response, suffix)
      mu <- paste0("mu", response, suffix)
      tau <- paste0("TAU", tolower(response), suffix)
      model_lines <- c(model_lines, paste0("  ", response_var, "[1:N] ~ dmnorm(", mu, "[], ", tau, ")"))
    }
  }

  model_lines <- c(model_lines, "  # Priors")

  # Priors for alpha and lambda, tau, sigma
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      model_lines <- c(
        model_lines,
        paste0("  alpha", response, suffix, " ~ dnorm(0, 1.0E-6)"),
        paste0("  lambda", response, suffix, " ~ dunif(0, 1)"),
        paste0("  tau", response, suffix, " ~ dgamma(1, 1)"),
        paste0("  sigma", response, suffix, " <- 1/sqrt(tau", response, suffix, ")")
      )
    }
  }

  # Priors for unique betas
  unique_betas <- unique(unlist(beta_counter))
  for (beta in unique_betas) {
    model_lines <- c(model_lines, paste0("  ", beta, " ~ dnorm(0, 1.0E-6)"))
  }

  # Covariance structure
  model_lines <- c(model_lines, "  # Covariance structure")
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      model_lines <- c(
        model_lines,
        paste0("  Mlam", response, suffix, " <- lambda", response, suffix, "*VCV + (1-lambda", response, suffix, ")*ID"),
        paste0("  TAU", tolower(response), suffix, " <- tau", response, suffix, "*inverse(Mlam", response, suffix, ")")
      )
    }
  }

  model_lines <- c(model_lines, "}")
  jags_model_string <- paste(model_lines, collapse = "\n")
  return(jags_model_string)
}
