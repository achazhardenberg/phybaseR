#' Generate a JAGS model string for Phylogenetic Bayesian SEM (PhyBaSE)
#'
#' This function builds the model code to be passed to JAGS based on a set of structural equations.
#' It supports both single and multiple phylogenetic trees (to account for phylogenetic uncertainty).
#' Missing values are handled both in the response and predictor variables treating all of them as stochastic nodes.
#'
#' @param equations A list of model formulas.
#' @param multi.tree Logical; if \code{TRUE}, incorporates phylogenetic uncertainty by sampling across a set of trees.
#' @param variability Optional character vector or named character vector of variable names that have measurement error or within-species variability.
#'   If named, the names should be the variable names and the values should be the type of variability: "se" (for mean and standard error) or "reps" (for repeated measures).
#'   If unnamed, it defaults to "se" for all specified variables.
#'   \itemize{
#'     \item "se": Expects \code{Var_mean} and \code{Var_se} in the data. The model fixes observation error: \code{Var_mean ~ dnorm(Var, 1/Var_se^2)}.
#'     \item "reps": Expects \code{Var_obs} (matrix) and \code{N_reps_Var} (vector) in the data. The model estimates observation error: \code{Var_obs[i,j] ~ dnorm(Var[i], Var_tau)}.
#'   }
#' @param distribution Optional named character vector specifying the distribution for response variables.
#'   Default is "gaussian" for all variables. Supported values: "gaussian", "binomial".
#'   For "binomial" variables, the model uses a logit link and a Bernoulli likelihood, with phylogenetic correlation modeled on the latent scale.
#' @param vars_with_na Optional character vector of response variable names that have missing data.
#'   These variables will use element-wise likelihoods instead of multivariate normal.
#'
#' @return A character string containing the JAGS model code.
#'
#' @details
#' The generated model includes:
#' \itemize{
#'   \item Linear predictors and multivariate normal likelihoods for each response variable.
#'   \item Priors for intercepts (\code{alpha}), slopes (\code{beta}), lambda parameters (\code{lambda}), and residual precisions (\code{tau}).
#'   \item Phylogenetic covariance modeled via a single \code{VCV} matrix (when \code{multi.tree = FALSE}) or a 3D array \code{multiVCV[,,K]} with categorical sampling across trees (when \code{multi.tree = TRUE}).
#'   \item (Optional) Observation models for variables with measurement error:
#'     \itemize{
#'       \item Type "se": \code{Var_mean ~ dnorm(Var, 1/Var_se^2)}
#'       \item Type "reps": \code{Var_obs[i,j] ~ dnorm(Var[i], Var_tau)}
#'     }
#'   \item (Optional) Generalized linear mixed models for non-Gaussian responses (e.g., binomial).
#'   \item (Optional) Element-wise likelihoods for response variables with missing data.
#' }
#'
#' @examples
#' eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
#' cat(phybase_model(eqs, multi.tree = TRUE))
#'
#' @export
#' @importFrom stats formula terms setNames
#'
phybase_model <- function(
  equations,
  multi.tree = FALSE,
  variability = NULL,
  distribution = NULL,
  vars_with_na = NULL
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  beta_counter <- list()
  response_counter <- list()

  # Parse equations
  eq_list <- lapply(equations, function(eq) {
    response <- as.character(formula(eq))[2]
    predictors <- attr(terms(formula(eq)), "term.labels")
    list(response = response, predictors = predictors)
  })

  # Track all variables (for imputation)
  all_vars <- unique(unlist(lapply(eq_list, function(eq) {
    c(eq$response, eq$predictors)
  })))

  # Handle variability argument
  variability_list <- list()
  if (!is.null(variability)) {
    if (is.null(names(variability))) {
      # Default to "se" if unnamed
      variability_list <- setNames(rep("se", length(variability)), variability)
    } else {
      variability_list <- variability
    }
  }

  # Handle distribution argument
  dist_list <- list()
  if (!is.null(distribution)) {
    dist_list <- distribution
  }

  # Start model
  model_lines <- c("model {", "  # Structural equations", "  for (i in 1:N) {")

  for (j in seq_along(eq_list)) {
    eq <- eq_list[[j]]
    response <- eq$response
    predictors <- eq$predictors
    dist <- dist_list[[response]] %||% "gaussian"

    # Count and assign unique suffix for the response variable
    response_count <- response_counter[[response]] %||% 0
    response_count <- response_count + 1
    response_counter[[response]] <- response_count
    suffix <- if (response_count == 1) "" else as.character(response_count)

    alpha <- paste0("alpha", response, suffix)

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

    if (dist == "gaussian") {
      mu <- paste0("mu", response, suffix)
      model_lines <- c(model_lines, paste0("    ", mu, "[i] <- ", linpred))
    } else if (dist == "binomial") {
      # Binomial: logit(p) = linpred + error
      # error ~ dmnorm(0, TAU)
      # We define the error mean here as 0
      mu_err <- paste0("mu_err_", response, suffix)
      err <- paste0("err_", response, suffix)
      p <- paste0("p_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    ", mu_err, "[i] <- 0"),
        paste0("    logit(", p, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", response, "[i] ~ dbern(", p, "[i])")
      )
    } else {
      stop(paste("Unknown distribution:", dist))
    }
  }

  model_lines <- c(model_lines, "  }", "  # Multivariate normal likelihoods")

  # Likelihoods for responses
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"
      response_var <- paste0(response, suffix)

      tau <- paste0("TAU", tolower(response), suffix)

      if (dist == "gaussian") {
        mu <- paste0("mu", response, suffix)
        tau_scalar <- paste0("tau", response, suffix)

        # Check if this variable has missing data
        if (!is.null(vars_with_na) && response %in% vars_with_na) {
          # Use Latent Variable (GLMM) approach for missing data
          # Y[i] ~ dnorm(mu[i] + err[i], tau_res)
          # err[1:N] ~ dmnorm(0, tau_phylo * inv(VCV))

          err <- paste0("err_", response, suffix)
          mu_err <- paste0("mu_err_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0(
              "  # GLMM likelihood for missing data (preserves phylo signal)"
            ),
            paste0("  for (i in 1:N) {"),
            paste0("    ", mu_err, "[i] <- 0"),
            paste0(
              "    ",
              response_var,
              "[i] ~ dnorm(",
              mu,
              "[i] + ",
              err,
              "[i], ",
              tau_res,
              ")"
            ),
            paste0("  }")
          )
        } else {
          # Standard MVN for complete data
          model_lines <- c(
            model_lines,
            paste0("  ", response_var, "[1:N] ~ dmnorm(", mu, "[], ", tau, ")")
          )
        }
      } else if (dist == "binomial") {
        # For binomial, the error term has the phylogenetic structure
        err <- paste0("err_", response, suffix)
        mu_err <- paste0("mu_err_", response, suffix)
        model_lines <- c(
          model_lines,
          paste0("  ", err, "[1:N] ~ dmnorm(", mu_err, "[], ", tau, ")")
        )
      }
    }
  }

  # Measurement error / Variability
  if (length(variability_list) > 0) {
    model_lines <- c(model_lines, "  # Measurement error / Variability")

    for (var in names(variability_list)) {
      type <- variability_list[[var]]

      # Ensure the variable is in the model
      if (!var %in% all_vars) {
        warning(paste(
          "Variable",
          var,
          "specified in 'variability' but not found in equations."
        ))
        next
      }

      if (type == "se") {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "_mean[i] ~ dnorm(",
            var,
            "[i], ",
            var,
            "_tau_obs[i])"
          ),
          paste0(
            "    ",
            var,
            "_tau_obs[i] <- 1/(",
            var,
            "_se[i] * ",
            var,
            "_se[i])"
          ),
          paste0("  }")
        )
      } else if (type == "reps") {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0("    for (j in 1:N_reps_", var, "[i]) {"),
          paste0(
            "      ",
            var,
            "_obs[i, j] ~ dnorm(",
            var,
            "[i], ",
            var,
            "_tau)"
          ),
          paste0("    }"),
          paste0("  }"),
          paste0("  ", var, "_tau ~ dgamma(1, 1)"),
          paste0("  ", var, "_sigma <- 1/sqrt(", var, "_tau)")
        )
      } else {
        warning(paste("Unknown variability type:", type, "for variable", var))
      }
    }
  }

  model_lines <- c(model_lines, "  # Priors for structural parameters")

  # Priors for alpha, lambda, tau, sigma
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      model_lines <- c(
        model_lines,
        paste0("  alpha", response, suffix, " ~ dnorm(0, 1.0E-6)")
      )

      # Only generate lambda/tau priors if NOT using GLMM (missing data)
      # For GLMM, we generate specific priors in the covariance section
      if (is.null(vars_with_na) || !response %in% vars_with_na) {
        model_lines <- c(
          model_lines,
          paste0("  lambda", response, suffix, " ~ dunif(0, 1)"),
          paste0("  tau", response, suffix, " ~ dgamma(1, 1)"),
          paste0(
            "  sigma",
            response,
            suffix,
            " <- 1/sqrt(tau",
            response,
            suffix,
            ")"
          )
        )
      }
    }
  }

  # Priors for regression coefficients
  unique_betas <- unique(unlist(beta_counter))
  for (beta in unique_betas) {
    model_lines <- c(model_lines, paste0("  ", beta, " ~ dnorm(0, 1.0E-6)"))
  }

  # Phylogenetic tree selection (if multi.tree)
  if (multi.tree) {
    model_lines <- c(
      model_lines,
      "",
      "  for (k in 1:Ntree) {",
      "    p[k] <- 1 / Ntree",
      "  }",
      "  K ~ dcat(p[])"
    )
  }

  # Covariance structures for responses
  model_lines <- c(model_lines, "  # Covariance structure for responses")
  for (response in names(response_counter)) {
    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"

      # Skip TAU matrix for variables with missing data (using element-wise) or binomial error terms
      use_glmm <- (!is.null(vars_with_na) && response %in% vars_with_na)

      if (dist == "gaussian" && !use_glmm) {
        if (multi.tree) {
          model_lines <- c(
            model_lines,
            paste0(
              "  Mlam",
              response,
              suffix,
              " <- lambda",
              response,
              suffix,
              "*multiVCV[,,K] + (1-lambda",
              response,
              suffix,
              ")*ID"
            ),
            paste0(
              "  TAU",
              tolower(response),
              suffix,
              " <- tau",
              response,
              suffix,
              "*inverse(Mlam",
              response,
              suffix,
              ")"
            )
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0(
              "  Mlam",
              response,
              suffix,
              " <- lambda",
              response,
              suffix,
              "*VCV + (1-lambda",
              response,
              suffix,
              ")*ID"
            ),
            paste0(
              "  TAU",
              tolower(response),
              suffix,
              " <- tau",
              response,
              suffix,
              "*inverse(Mlam",
              response,
              suffix,
              ")"
            )
          )
        }
      } else if (dist == "gaussian" && use_glmm) {
        # GLMM covariance for latent error term
        # err ~ dmnorm(0, tau_phylo * inv(VCV))

        err <- paste0("err_", response, suffix)
        mu_err <- paste0("mu_err_", response, suffix)
        tau_phylo <- paste0("tau_phylo_", response, suffix)
        tau_res <- paste0("tau_res_", response, suffix)

        # Priors
        model_lines <- c(
          model_lines,
          paste0("  ", tau_res, " ~ dgamma(1, 1)"),
          paste0("  ", tau_phylo, " ~ dgamma(1, 1)"),
          # Calculate lambda for reporting
          paste0(
            "  lambda",
            response,
            suffix,
            " <- (1/",
            tau_phylo,
            ") / ((1/",
            tau_phylo,
            ") + (1/",
            tau_res,
            "))"
          )
        )

        if (multi.tree) {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_phylo_",
              response,
              suffix,
              " <- ",
              tau_phylo,
              " * inverse(multiVCV[,,K])"
            ),
            paste0(
              "  ",
              err,
              "[1:N] ~ dmnorm(",
              mu_err,
              "[], TAU_phylo_",
              response,
              suffix,
              ")"
            )
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_phylo_",
              response,
              suffix,
              " <- ",
              tau_phylo,
              " * inverse(VCV)"
            ),
            paste0(
              "  ",
              err,
              "[1:N] ~ dmnorm(",
              mu_err,
              "[], TAU_phylo_",
              response,
              suffix,
              ")"
            )
          )
        }
      } else if (dist == "binomial") {
        # Binomial uses TAU for error term
        if (multi.tree) {
          model_lines <- c(
            model_lines,
            paste0(
              "  Mlam",
              response,
              suffix,
              " <- lambda",
              response,
              suffix,
              "*multiVCV[,,K] + (1-lambda",
              response,
              suffix,
              ")*ID"
            ),
            paste0(
              "  TAU",
              tolower(response),
              suffix,
              " <- inverse(tau",
              response,
              suffix,
              "*Mlam",
              response,
              suffix,
              "[,])"
            )
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0(
              "  Mlam",
              response,
              suffix,
              " <- lambda",
              response,
              suffix,
              "*VCV + (1-lambda",
              response,
              suffix,
              ")*ID"
            ),
            paste0(
              "  TAU",
              tolower(response),
              suffix,
              " <- inverse(tau",
              response,
              suffix,
              "*Mlam",
              response,
              suffix,
              "[,])"
            )
          )
        }
      }
    }
  }

  # Imputation priors for predictors (those not modeled as responses)
  model_lines <- c(model_lines, "  # Predictor priors for imputation")
  non_response_vars <- setdiff(all_vars, names(response_counter))
  for (var in non_response_vars) {
    model_lines <- c(
      model_lines,
      paste0("  for (i in 1:N) {"),
      paste0("    mu", var, "[i] <- 0"),
      paste0("  }"),
      paste0("  ", var, "[1:N] ~ dmnorm(mu", var, "[], TAU", tolower(var), ")"),
      paste0("  lambda", var, " ~ dunif(0, 1)"),
      paste0("  tau", var, " ~ dgamma(1, 1)"),
      paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
    )

    if (multi.tree) {
      model_lines <- c(
        model_lines,
        paste0(
          "  Mlam",
          var,
          " <- lambda",
          var,
          "*multiVCV[,,K] + (1 - lambda",
          var,
          ")*ID"
        ),
        paste0("  TAU", tolower(var), " <- tau", var, "*inverse(Mlam", var, ")")
      )
    } else {
      model_lines <- c(
        model_lines,
        paste0(
          "  Mlam",
          var,
          " <- lambda",
          var,
          "*VCV + (1 - lambda",
          var,
          ")*ID"
        ),
        paste0("  TAU", tolower(var), " <- tau", var, "*inverse(Mlam", var, ")")
      )
    }
  }

  model_lines <- c(model_lines, "}")
  jags_model_string <- paste(model_lines, collapse = "\n")
  return(jags_model_string)
}
