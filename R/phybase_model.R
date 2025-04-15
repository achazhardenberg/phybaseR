#' Generate a JAGS model string for Phylogenetic Bayesian SEM (PhyBaSE)
#'
#' This function builds the model code to be passed to JAGS based on a set of structural equations.
#' It supports both single and multiple phylogenetic trees (to account for phylogenetic uncertainty).
#'
#' @param equations A list of model formulas (one per structural equation), e.g., \code{list(Y ~ X1 + X2, Z ~ Y)}.
#' @param multi.tree Logical. If \code{TRUE}, generates a model that samples from a set of phylogenetic variance-covariance matrices (\code{multiVCV}) to account for phylogenetic uncertainty. Defaults to \code{FALSE}.
#'
#' @return A character string containing the JAGS model code.
#'
#' @details
#' The generated model includes:
#' \itemize{
#'   \item Linear predictors and multivariate normal likelihoods for each response variable.
#'   \item Priors for intercepts (\code{alpha}), slopes (\code{beta}), lambda parameters (\code{lambda}), and residual precisions (\code{tau}).
#'   \item Phylogenetic covariance modeled via a single \code{VCV} matrix (when \code{multi.tree = FALSE}) or a 3D array \code{multiVCV[,,K]} with categorical sampling across trees (when \code{multi.tree = TRUE}).
#' }
#'
#' @examples
#' eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
#' cat(phybase_model(eqs, multi.tree = TRUE))
#'
#' @export
phybase_model <- function(equations, multi.tree = FALSE) {

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
  if (multi.tree) {
    model_lines <- c(model_lines,
                     "",
                     "  for (k in 1:Ntree) {",
                     "    p[k] <- 1 / Ntree",
                     "  }",
                     "  K ~ dcat(p[])")
  }

  # Covariance structure
  model_lines <- c(model_lines, "  # Covariance structure")

  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      if (multi.tree) {
        model_lines <- c(
          model_lines,
          paste0("  Mlam", response, suffix, " <- lambda", response, suffix, "*multiVCV[,,K] + (1-lambda", response, suffix, ")*ID"),
          paste0("  TAU", tolower(response), suffix, " <- tau", response, suffix, "*inverse(Mlam", response, suffix, ")")
        )
      }
      else {
        model_lines <- c(
          model_lines,
          paste0("  Mlam", response, suffix, " <- lambda", response, suffix, "*VCV + (1-lambda", response, suffix, ")*ID"),
          paste0("  TAU", tolower(response), suffix, " <- tau", response, suffix, "*inverse(Mlam", response, suffix, ")")
        )
      }
    }
  }
  model_lines <- c(model_lines, "}")
  jags_model_string <- paste(model_lines, collapse = "\n")
  return(jags_model_string)
}
