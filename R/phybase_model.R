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
#'   Default is "gaussian" for all variables. Supported values: "gaussian", "binomial", "multinomial".
#'   For "binomial" variables, the model uses a logit link and a Bernoulli likelihood, with phylogenetic correlation modeled on the latent scale.
#' @param vars_with_na Optional character vector of response variable names that have missing data.
#'   These variables will use element-wise likelihoods instead of multivariate normal.
#' @param induced_correlations Optional list of variable pairs with induced correlations
#'   from latent variables. Each element should be a character vector of length 2 specifying
#'   the pair of variables that share a latent common cause.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{model}: A character string containing the JAGS model code.
#'   \item \code{parameter_map}: A data frame mapping response variables to their predictors and parameter names.
#' }
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
#' cat(phybase_model(eqs, multi.tree = TRUE)$model)
#'
#' @param optimise Logical. If TRUE (default), use random effects formulation for 4.6× speedup.
#'   If FALSE, use original marginal covariance formulation.
#' @export
#' @importFrom stats formula terms setNames sd
#' @importFrom utils combn
#'
phybase_model <- function(
  equations,
  multi.tree = FALSE,
  variability = NULL,
  distribution = NULL,
  vars_with_na = NULL,
  induced_correlations = NULL,
  latent = NULL,
  standardize_latent = TRUE,
  optimise = TRUE,
  structure_names = NULL
) {
  # Helper: returns b if a is NULL or if a is a list element that doesn't exist
  `%||%` <- function(a, b) {
    tryCatch(if (!is.null(a)) a else b, error = function(e) b)
  }

  independent <- is.null(structure_names) || length(structure_names) == 0

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

  # Identify variables involved in induced correlations
  correlated_vars <- if (!is.null(induced_correlations)) {
    unique(unlist(induced_correlations))
  } else {
    character(0)
  }

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

  param_map <- list()

  # Start model
  model_lines <- c(
    "model {",
    "  # Dummy usage of ID to prevent warnings in GLMM-only models",
    "  dummy_ID <- ID[1,1]",
    "  # Structural equations",
    "  for (i in 1:N) {"
  )

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
        # Use consistent naming: beta_Response_Predictor
        beta_name <- paste0("beta_", response, suffix, "_", pred)
        beta_counter[[key]] <- beta_name
      }
      beta_name <- beta_counter[[key]]
      linpred <- paste0(linpred, " + ", beta_name, "*", pred, "[i]")

      # Store in parameter map
      param_map[[length(param_map) + 1]] <- list(
        response = response,
        predictor = pred,
        parameter = beta_name,
        equation_index = j
      )
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
    } else if (dist == "multinomial") {
      # Multinomial: K categories
      # L[i, 1] <- 0
      # L[i, k] <- alpha[k] + beta*X + err[i, k]

      K_var <- paste0("K_", response)
      err <- paste0("err_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Multinomial linear predictor for ", response),
        paste0("    L_", response, "[i, 1] <- 0"),
        paste0("    for (k in 2:", K_var, ") {")
      )

      # Linear predictor for k-th dimension
      # Note: We use array indexing [k] for parameters
      linpred_k <- paste0("alpha_", response, "[k]")

      for (pred in predictors) {
        # We force a specific name for multinomial betas to ensure array usage
        beta_name <- paste0("beta_", response, "_", pred)
        linpred_k <- paste0(linpred_k, " + ", beta_name, "[k] * ", pred, "[i]")

        # Map (only once)
        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(beta_counter)) {
          beta_counter[[key]] <- beta_name # Mark as used
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = paste0(beta_name, "[]"),
            equation_index = j
          )
        }
      }

      linpred_k <- paste0(linpred_k, " + ", err, "[i, k]")

      model_lines <- c(
        model_lines,
        paste0("      L_", response, "[i, k] <- ", linpred_k),
        "    }"
      )

      # Softmax
      model_lines <- c(
        model_lines,
        paste0("    # Softmax for ", response),
        paste0("    for (k in 1:", K_var, ") {"),
        paste0(
          "      exp_L_",
          response,
          "[i, k] <- exp(L_",
          response,
          "[i, k])"
        ),
        "    }",
        paste0(
          "    sum_exp_L_",
          response,
          "[i] <- sum(exp_L_",
          response,
          "[i, 1:",
          K_var,
          "])"
        ),
        paste0("    for (k in 1:", K_var, ") {"),
        paste0(
          "      p_",
          response,
          "[i, k] <- exp_L_",
          response,
          "[i, k] / sum_exp_L_",
          response,
          "[i]"
        ),
        "    }",
        paste0(
          "    ",
          response,
          "[i] ~ dcat(p_",
          response,
          "[i, 1:",
          K_var,
          "])"
        )
      )
    } else if (dist == "ordinal") {
      # Ordinal: Cumulative Logit (Proportional Odds)
      # P(Y <= k) = logit^(-1)(cutpoint[k] - eta)
      # eta = linpred + error

      K_var <- paste0("K_", response)
      err <- paste0("err_", response, suffix)
      eta <- paste0("eta_", response, suffix)

      # Linear predictor (eta)
      # Note: No intercept in eta (intercept is absorbed into cutpoints)
      linpred_no_int <- "0"
      for (pred in predictors) {
        beta_name <- paste0("beta_", response, "_", pred)
        linpred_no_int <- paste0(
          linpred_no_int,
          " + ",
          beta_name,
          " * ",
          pred,
          "[i]"
        )

        # Map
        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(beta_counter)) {
          beta_counter[[key]] <- beta_name
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = beta_name,
            equation_index = j
          )
        }
      }

      model_lines <- c(
        model_lines,
        paste0("    # Ordinal linear predictor for ", response),
        paste0("    ", eta, "[i] <- ", linpred_no_int, " + ", err, "[i]"),

        # Cumulative probabilities
        paste0("    for (k in 1:(", K_var, "-1)) {"),
        paste0(
          "      logit(Q_",
          response,
          "[i, k]) <- cutpoint_",
          response,
          "[k] - ",
          eta,
          "[i]"
        ),
        "    }",

        # Category probabilities
        paste0("    p_", response, "[i, 1] <- Q_", response, "[i, 1]"),
        paste0("    for (k in 2:(", K_var, "-1)) {"),
        paste0(
          "      p_",
          response,
          "[i, k] <- Q_",
          response,
          "[i, k] - Q_",
          response,
          "[i, k-1]"
        ),
        "    }",
        paste0(
          "    p_",
          response,
          "[i, ",
          K_var,
          "] <- 1 - Q_",
          response,
          "[i, ",
          K_var,
          "-1]"
        ),

        # Likelihood
        paste0(
          "    ",
          response,
          "[i] ~ dcat(p_",
          response,
          "[i, 1:",
          K_var,
          "])"
        )
      )
    } else if (dist == "poisson") {
      # Poisson: log(μ) = linpred + error
      # Naturally handles overdispersion via epsilon

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Poisson log link for ", response),
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", response, "[i] ~ dpois(", mu, "[i])")
      )
    } else if (dist == "negbinomial") {
      # Negative Binomial: log(μ) = linpred + error
      # Y ~ NegBin(p, r) where p = r/(r+μ) and r = size parameter

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      p <- paste0("p_", response, suffix)
      r <- paste0("r_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Negative Binomial log link for ", response),
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])"),
        paste0("    ", response, "[i] ~ dnegbin(", p, "[i], ", r, ")")
      )
    } else {
      stop(paste("Unknown distribution:", dist))
    }
  }

  model_lines <- c(model_lines, "  }", "  # Multivariate normal likelihoods")

  # Likelihoods for responses
  for (response in names(response_counter)) {
    # Skip if involved in induced correlations (handled separately)
    if (response %in% correlated_vars) {
      next
    }

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"
      response_var <- paste0(response, suffix)

      tau <- paste0("TAU", tolower(response), suffix)

      if (dist == "gaussian") {
        mu <- paste0("mu", response, suffix)
        tau_scalar <- paste0("tau", response, suffix)

        # Check if this variable has missing data
        if (!is.null(vars_with_na) && response %in% vars_with_na && !optimise) {
          # Use Latent Variable (GLMM) approach for missing data (Only if optimisation disabled)
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
          if (independent) {
            # Independent Model (No random effects)
            # Y[i] ~ dnorm(mu[i], tau)
            # Note: We use tau_e for consistency with optimized model naming if desired,
            # but tau is standard for simple normal. Let's use tau_e to distinguish from matrix TAU
            tau_e <- paste0("tau_e_", response, suffix)

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:N) {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i], ",
                tau_e,
                ")"
              ),
              paste0("  }")
            )
          } else if (optimise) {
            # Optimized Random Effects Formulation (Additive)
            additive_terms <- ""

            for (s_name in structure_names) {
              s_suffix <- if (length(structure_names) > 1) {
                paste0("_", s_name)
              } else {
                ""
              }

              u_std <- paste0("u_std_", response, suffix, s_suffix)
              u <- paste0("u_", response, suffix, s_suffix)
              tau_u <- paste0("tau_u_", response, suffix, s_suffix)

              prec_name <- paste0("Prec_", s_name)
              prec_idx <- paste0(prec_name, "[1:N, 1:N]")
              if (multi.tree && s_name == "phylo") {
                prec_idx <- paste0(prec_name, "[1:N, 1:N, K]")
              }

              model_lines <- c(
                model_lines,
                paste0(
                  "  ",
                  u_std,
                  "[1:N] ~ dmnorm(zeros[1:N], ",
                  prec_idx,
                  ")"
                ),
                paste0(
                  "  for (i in 1:N) { ",
                  u,
                  "[i] <- ",
                  u_std,
                  "[i] / sqrt(",
                  tau_u,
                  ") }"
                )
              )
              additive_terms <- paste0(additive_terms, " + ", u, "[i]")
            }

            tau_e <- paste0("tau_e_", response, suffix)

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:N) {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i]",
                additive_terms,
                ", ",
                tau_e,
                ")"
              ),
              paste0("  }")
            )
          } else {
            # Standard MVN for complete data (Marginal)
            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                response_var,
                "[1:N] ~ dmnorm(",
                mu,
                "[1:N], ",
                tau,
                ")"
              )
            )
          }
        }
      } else if (dist == "binomial") {
        # For binomial, the error term has the phylogenetic structure
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Binomial: Just error term (epsilon)
          # err[i] ~ dnorm(0, tau_e)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0("  # Independent residual error for binomial: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]"),
            paste0("  }")
          )
        } else if (optimise) {
          # Random Effects Formulation
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree: use Prec_phylo_fixed[,,K] instead of Prec_phylo_fixed[,]
          prec_index <- if (multi.tree) {
            "Prec_phylo_fixed[1:N, 1:N, K]"
          } else {
            "Prec_phylo_fixed[1:N, 1:N]"
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for binomial: ", response),
            paste0("  ", u_std, "[1:N] ~ dmnorm(zeros[1:N], ", prec_index, ")"),
            paste0("  for (i in 1:N) {"),
            paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", u, "[i] + ", epsilon, "[i]"),
            paste0("  }")
          )
        } else {
          # Marginal Formulation (original)
          mu_err <- paste0("mu_err_", response, suffix)
          model_lines <- c(
            model_lines,
            paste0("  ", err, "[1:N] ~ dmnorm(", mu_err, "[], ", tau, ")")
          )
        }
      } else if (dist == "multinomial") {
        # Multinomial error terms: err[1:N, k]
        # Independent phylogenetic effects for each k (2..K)
        err <- paste0("err_", response, suffix)
        K_var <- paste0("K_", response)

        if (independent) {
          # Define variables for independent mode
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0(
              "  # Independent residual error for multinomial: ",
              response
            ),
            paste0("  for (k in 2:", K_var, ") {"),
            paste0("    for (i in 1:N) {"),
            paste0("      ", epsilon, "[i, k] ~ dnorm(0, ", tau_e, "[k])"),
            paste0("      ", err, "[i, k] <- ", epsilon, "[i, k]"),
            paste0("    }"),
            paste0("  }")
          )
        } else if (optimise) {
          # Optimised Random Effects Formulation for Multinomial
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree: use Prec_phylo_fixed[,,K] instead of Prec_phylo_fixed[,]
          prec_index <- if (multi.tree) {
            "Prec_phylo_fixed[1:N, 1:N, K]"
          } else {
            "Prec_phylo_fixed[1:N, 1:N]"
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for multinomial: ", response),
            paste0("  for (k in 2:", K_var, ") {"),
            paste0(
              "    ",
              u_std,
              "[1:N, k] ~ dmnorm(zeros[1:N], ",
              prec_index,
              ")"
            ),
            paste0("    for (i in 1:N) {"),
            paste0(
              "      ",
              u,
              "[i, k] <- ",
              u_std,
              "[i, k] / sqrt(",
              tau_u,
              "[k])"
            ),
            paste0("      ", epsilon, "[i, k] ~ dnorm(0, ", tau_e, "[k])"),
            paste0(
              "      ",
              err,
              "[i, k] <- ",
              u,
              "[i, k] + ",
              epsilon,
              "[i, k]"
            ),
            paste0("    }"),
            paste0("  }")
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0("  # Multinomial phylogenetic errors for ", response),
            paste0("  for (k in 2:", K_var, ") {"),
            paste0(
              "    ",
              err,
              "[1:N, k] ~ dmnorm(zero_vec[], TAU_",
              tolower(response),
              "_",
              suffix,
              "[,,k])"
            ),
            "  }"
          )
        }
      } else if (dist == "ordinal") {
        # Ordinal error term: err[1:N]
        # Single phylogenetic effect (unlike multinomial with K-1 effects)
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Ordinal: Just error term (epsilon)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0("  # Independent residual error for ordinal: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]"),
            paste0("  }")
          )
        } else {
          # Random Effects Formulation (Default/Optimised)
          # Note: Ordinal currently only supports optimized formulation in this codebase structure
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree
          prec_index <- if (multi.tree) {
            "Prec_phylo_fixed[1:N, 1:N, K]"
          } else {
            "Prec_phylo_fixed[1:N, 1:N]"
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for ordinal: ", response),
            paste0("  ", u_std, "[1:N] ~ dmnorm(zeros[1:N], ", prec_index, ")"),
            paste0("  for (i in 1:N) {"),
            paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", u, "[i] + ", epsilon, "[i]"),
            paste0("  }")
          )
        }
      } else if (dist == "poisson") {
        # Poisson error term: err[1:N]
        # Single phylogenetic effect (like ordinal)
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Poisson: Just error term (epsilon/overdispersion)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0("  # Independent residual error for Poisson: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]"),
            paste0("  }")
          )
        } else if (optimise) {
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree
          prec_index <- if (multi.tree) {
            "Prec_phylo_fixed[1:N, 1:N, K]"
          } else {
            "Prec_phylo_fixed[1:N, 1:N]"
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Poisson: ", response),
            paste0("  ", u_std, "[1:N] ~ dmnorm(zeros[1:N], ", prec_index, ")"),
            paste0("  for (i in 1:N) {"),
            paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", u, "[i] + ", epsilon, "[i]"),
            paste0("  }")
          )
        }
      } else if (dist == "negbinomial") {
        # Negative Binomial error term: err[1:N]
        # Single phylogenetic effect (like Poisson/ordinal)
        err <- paste0("err_", response, suffix)

        if (optimise) {
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree
          prec_index <- if (multi.tree) {
            "Prec_phylo_fixed[1:N, 1:N, K]"
          } else {
            "Prec_phylo_fixed[1:N, 1:N]"
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Negative Binomial: ", response),
            paste0("  ", u_std, "[1:N] ~ dmnorm(zeros[1:N], ", prec_index, ")"),
            paste0("  for (i in 1:N) {"),
            paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", u, "[i] + ", epsilon, "[i]"),
            paste0("  }")
          )
        }
      }
    }
  }

  # Handle Induced Correlations (Latent Variables)
  if (!is.null(induced_correlations)) {
    model_lines <- c(model_lines, "  # Induced Correlations (Latent Variables)")

    # Track variables and their error terms
    # vars_error_terms[[var]] <- c("term1", "term2", ...)
    vars_error_terms <- list()
    processed_params <- c() # Track tau_res/sigma_res/tau_phylo/tau_obs definitions

    # 1. Process pairs to generate correlated error terms
    for (pair in induced_correlations) {
      var1 <- pair[1]
      var2 <- pair[2]

      # Initialize lists if needed
      if (is.null(vars_error_terms[[var1]])) {
        vars_error_terms[[var1]] <- c()
      }
      if (is.null(vars_error_terms[[var2]])) {
        vars_error_terms[[var2]] <- c()
      }

      # Define pair-specific names
      res_err <- paste0("err_res_", var1, "_", var2)
      tau_res_matrix <- paste0("TAU_res_", var1, "_", var2)
      cov_matrix <- paste0("cov_", var1, "_", var2)

      # Define Residual Precision/Sigma Parameters (Deduplicated)
      # sigma_res_var is shared across all pairs involving var
      for (v in c(var1, var2)) {
        if (!(v %in% processed_params)) {
          model_lines <- c(
            model_lines,
            paste0("  tau_res_", v, " ~ dgamma(1, 1)"),
            paste0("  sigma_res_", v, " <- 1/sqrt(tau_res_", v, ")")
          )
          processed_params <- c(processed_params, v)
        }
      }

      # Define Correlation and Covariance Matrix
      model_lines <- c(
        model_lines,
        paste0("  # Correlated residuals between ", var1, " and ", var2),
        paste0(
          "  z_",
          var1,
          "_",
          var2,
          " ~ dnorm(0, 1)  # Fisher Z transformation"
        ),
        paste0("  rho_", var1, "_", var2, " <- tanh(z_", var1, "_", var2, ")"),
        # Restore Diagonal Elements
        paste0("  ", cov_matrix, "[1, 1] <- 1/tau_res_", var1),
        paste0("  ", cov_matrix, "[2, 2] <- 1/tau_res_", var2),
        paste0(
          "  ",
          cov_matrix,
          "[1, 2] <- rho_",
          var1,
          "_",
          var2,
          " * sigma_res_",
          var1,
          " * sigma_res_",
          var2
        ),
        paste0("  ", cov_matrix, "[2, 1] <- ", cov_matrix, "[1, 2]"),

        # Calculate Precision Matrix
        paste0(
          "  ",
          tau_res_matrix,
          "[1:2, 1:2] <- inverse(",
          cov_matrix,
          "[1:2, 1:2])"
        ),

        # Generate correlated error terms from MVN
        paste0("  for (i in 1:N) {"),
        paste0(
          "    ",
          res_err,
          "[i, 1:2] ~ dmnorm(zero_vec[1:2], ",
          tau_res_matrix,
          "[1:2, 1:2])"
        ),
        paste0("  }")
      )

      # Record error terms for each variable
      # Usage: err_res_var1_var2[i, 1] for var1, [i, 2] for var2
      vars_error_terms[[var1]] <- c(
        vars_error_terms[[var1]],
        paste0(res_err, "[i, 1]")
      )
      vars_error_terms[[var2]] <- c(
        vars_error_terms[[var2]],
        paste0(res_err, "[i, 2]")
      )
    }

    # 2. Generate Likelihood for each variable (summing error terms)
    for (var in names(vars_error_terms)) {
      err_terms <- vars_error_terms[[var]]
      suffix <- if ((response_counter[[var]] %||% 0) > 1) "1" else ""

      # Define Phylogenetic Error (if structure exists)
      phylo_term <- ""
      if (length(structure_names) > 0) {
        err_phylo <- paste0("err_phylo_", var)
        model_lines <- c(
          model_lines,
          paste0("  tau_phylo_", var, " ~ dgamma(1, 1)"),
          paste0(
            "  ",
            err_phylo,
            "[1:N] ~ dmnorm(zero_vec[], TAU_phylo_",
            var,
            ")"
          )
        )
        phylo_term <- paste0(" + ", err_phylo, "[i]")
      }

      # Define Observation Precision (Constant High Precision)
      model_lines <- c(
        model_lines,
        paste0("  tau_obs_", var, " <- 10000")
      )

      # Build sum string
      sum_res_errs <- paste(err_terms, collapse = " + ")

      model_lines <- c(
        model_lines,
        paste0("  for (i in 1:N) {"),
        paste0(
          "    ",
          var,
          "[i] ~ dnorm(mu",
          var,
          suffix,
          "[i]",
          phylo_term,
          " + ",
          sum_res_errs,
          ", tau_obs_",
          var,
          ")"
        ),
        paste0("  }")
      )
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
    # Skip correlated vars (handled above)
    if (response %in% correlated_vars) {
      next
    }

    # Skip multinomial, ordinal, poisson, and negbinomial (handled separately)
    dist <- dist_list[[response]] %||% "gaussian"
    if (
      dist == "multinomial" ||
        dist == "ordinal" ||
        dist == "poisson" ||
        dist == "negbinomial"
    ) {
      next
    }

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      model_lines <- c(
        model_lines,
        paste0("  alpha", response, suffix, " ~ dnorm(0, 1.0E-6)")
      )

      # Only generate lambda/tau priors if NOT using GLMM (missing data)
      # For GLMM, we generate specific priors in the covariance section
      # Only generate lambda/tau priors if NOT using GLMM (missing data),
      # UNLESS we are optimizing (which uses standard priors even for missing data)
      if (is.null(vars_with_na) || !response %in% vars_with_na || optimise) {
        if (independent) {
          # Independent Priors (only tau_e)
          # We can also generate sigma for monitoring convenience
          model_lines <- c(
            model_lines,
            paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
            paste0(
              "  sigma",
              response,
              suffix,
              " <- 1/sqrt(tau_e_",
              response,
              suffix,
              ")"
            )
          )
        } else if (optimise) {
          if (length(structure_names) > 1) {
            # Multiple Structures: Estimate independent variance components
            model_lines <- c(
              model_lines,
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)")
            )

            for (s_name in structure_names) {
              s_suffix <- if (length(structure_names) > 1) {
                paste0("_", s_name)
              } else {
                ""
              }
              tau_u <- paste0("tau_u_", response, suffix, s_suffix)

              model_lines <- c(
                model_lines,
                paste0("  ", tau_u, " ~ dgamma(1, 1)"),
                paste0(
                  "  sigma_",
                  response,
                  suffix,
                  s_suffix,
                  " <- 1/sqrt(",
                  tau_u,
                  ")"
                )
              )
            }
            model_lines <- c(
              model_lines,
              paste0(
                "  sigma_",
                response,
                "_res <- 1/sqrt(tau_e_",
                response,
                suffix,
                ")"
              )
            )
          } else {
            # Random Effects Priors (tau_u, tau_e) - Single Legacy
            model_lines <- c(
              model_lines,
              paste0("  tau_u_", response, suffix, " ~ dgamma(1, 1)"),
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
              # Derived parameters for backward compatibility
              paste0(
                "  lambda",
                response,
                suffix,
                " <- (1/tau_u_",
                response,
                suffix,
                ") / ((1/tau_u_",
                response,
                suffix,
                ") + (1/tau_e_",
                response,
                suffix,
                "))"
              ),
              paste0(
                "  sigma",
                response,
                suffix,
                " <- sqrt(1/tau_u_",
                response,
                suffix,
                " + 1/tau_e_",
                response,
                suffix,
                ")"
              )
            )
          }
        } else {
          # Marginal Priors (lambda, tau)
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
  }

  # Priors for multinomial parameters (arrays)
  for (response in names(response_counter)) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial") {
      K_var <- paste0("K_", response)
      if (independent) {
        model_lines <- c(
          model_lines,
          paste0("  # Independent Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),
          paste0("    tau_e_", response, "[k] ~ dgamma(1, 1)"),
          "  }"
        )
      } else if (optimise) {
        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),
          paste0("    tau_u_", response, "[k] ~ dgamma(1, 1)"),
          paste0("    tau_e_", response, "[k] ~ dgamma(1, 1)"),
          # Derived lambda
          paste0(
            "    lambda_",
            response,
            "[k] <- (1/tau_u_",
            response,
            "[k]) / ((1/tau_u_",
            response,
            "[k]) + (1/tau_e_",
            response,
            "[k]))"
          ),
          "  }"
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),
          paste0("    lambda_", response, "[k] ~ dunif(0, 1)"),
          paste0("    tau_", response, "[k] ~ dgamma(1, 1)"),
          "  }"
        )
      }

      # Betas (arrays)
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  for (k in 2:", K_var, ") {"),
              paste0("    ", beta_name, "[k] ~ dnorm(0, 1.0E-6)"),
              "  }"
            )
          }
        }
      }
    }
  }

  # Priors for ordinal parameters (cutpoints + variance components)
  for (response in names(response_counter)) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "ordinal") {
      K_var <- paste0("K_", response)

      # Loop over response instances (if there are repeats)
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, suffix, " (Ordinal)"),
          # Ordered cutpoints using delta transformation
          paste0("  cutpoint_raw_", response, suffix, "[1] ~ dnorm(0, 0.1)"),
          paste0(
            "  cutpoint_",
            response,
            suffix,
            "[1] <- cutpoint_raw_",
            response,
            suffix,
            "[1]"
          ),
          paste0("  for (k in 2:(", K_var, "-1)) {"),
          paste0("    cutpoint_raw_", response, suffix, "[k] ~ dnorm(0, 0.1)"),
          paste0(
            "    cutpoint_",
            response,
            suffix,
            "[k] <- cutpoint_",
            response,
            suffix,
            "[k-1] + exp(cutpoint_raw_",
            response,
            suffix,
            "[k])"
          ),
          "  }",
          # Variance components
          if (independent) {
            # Independent: Only tau_e
            c(
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
              paste0("  # No tau_u or lambda for independent model")
            )
          } else {
            # Structured: tau_u and tau_e and lambda
            c(
              paste0("  tau_u_", response, suffix, " ~ dgamma(1, 1)"),
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
              # Derived lambda
              paste0(
                "  lambda_",
                response,
                suffix,
                " <- (1/tau_u_",
                response,
                suffix,
                ") / ((1/tau_u_",
                response,
                suffix,
                ") + (1/tau_e_",
                response,
                suffix,
                "))"
              )
            )
          }
        )
      }

      # Betas for ordinal predictors
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  ", beta_name, " ~ dnorm(0, 1.0E-6)")
            )
          }
        }
      }
    }
  }

  # Priors for Poisson parameters (variance components)
  for (response in names(response_counter)) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "poisson") {
      # Loop over response instances (if there are repeats)
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, suffix, " (Poisson)"),
          # Variance components
          if (independent) {
            c(
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)")
            )
          } else {
            c(
              paste0("  tau_u_", response, suffix, " ~ dgamma(1, 1)"),
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
              # Derived lambda
              paste0(
                "  lambda_",
                response,
                suffix,
                " <- (1/tau_u_",
                response,
                suffix,
                ") / ((1/tau_u_",
                response,
                suffix,
                ") + (1/tau_e_",
                response,
                suffix,
                "))"
              )
            )
          }
        )
      }

      # Betas and intercepts for Poisson predictors
      model_lines <- c(
        model_lines,
        paste0("  alpha", response, " ~ dnorm(0, 1.0E-6)")
      )

      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  ", beta_name, " ~ dnorm(0, 1.0E-6)")
            )
          }
        }
      }
    }
  }

  # Priors for Negative Binomial parameters (variance components + size)
  for (response in names(response_counter)) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "negbinomial") {
      # Loop over response instances (if there are repeats)
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, suffix, " (Negative Binomial)"),
          # Variance components
          if (independent) {
            c(
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)")
            )
          } else {
            c(
              paste0("  tau_u_", response, suffix, " ~ dgamma(1, 1)"),
              paste0("  tau_e_", response, suffix, " ~ dgamma(1, 1)"),
              # Derived lambda
              paste0(
                "  lambda_",
                response,
                suffix,
                " <- (1/tau_u_",
                response,
                suffix,
                ") / ((1/tau_u_",
                response,
                suffix,
                ") + (1/tau_e_",
                response,
                suffix,
                "))"
              )
            )
          },
          # Size parameter (controls overdispersion)
          paste0(
            "  r_",
            response,
            suffix,
            " ~ dgamma(0.01, 0.01)  # Vague prior for size"
          )
        )
      }

      # Betas and intercepts for Negative Binomial predictors
      model_lines <- c(
        model_lines,
        paste0("  alpha", response, " ~ dnorm(0, 1.0E-6)")
      )

      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  ", beta_name, " ~ dnorm(0, 1.0E-6)")
            )
          }
        }
      }
    }
  }

  # Priors for correlated vars alphas (intercepts)
  for (var in correlated_vars) {
    model_lines <- c(model_lines, paste0("  alpha", var, " ~ dnorm(0, 1.0E-6)"))
  }

  # Priors for regression coefficients
  unique_betas <- unique(unlist(beta_counter))

  # Exclude multinomial betas (arrays)
  multinomial_betas <- c()
  for (response in names(response_counter)) {
    if ((dist_list[[response]] %||% "gaussian") == "multinomial") {
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            multinomial_betas <- c(
              multinomial_betas,
              paste0("beta_", response, "_", pred)
            )
          }
        }
      }
    }
  }
  unique_betas <- setdiff(unique_betas, multinomial_betas)

  for (beta in unique_betas) {
    model_lines <- c(model_lines, paste0("  ", beta, " ~ dnorm(0, 1.0E-6)"))
  }

  # Zero vectors for multivariate normals
  # Check if we need zero_vec (for induced_correlations OR multinomial)
  need_zero_vec <- !is.null(induced_correlations)
  for (response in names(response_counter)) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial") {
      need_zero_vec <- TRUE
      break
    }
  }

  if (need_zero_vec) {
    model_lines <- c(
      model_lines,
      "  for(k in 1:N) { zero_vec[k] <- 0 }"
    )
  }

  if (!is.null(induced_correlations)) {
    model_lines <- c(
      model_lines,
      "  zero_vec_2[1] <- 0",
      "  zero_vec_2[2] <- 0"
    )
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
    # Skip correlated vars (handled above)
    if (response %in% correlated_vars) {
      next
    }

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"

      # Skip TAU matrix for variables with missing data (using element-wise) or binomial error terms
      # Only use legacy GLMM blocking if optimisation implies marginal approach (i.e. optimise=FALSE)
      use_glmm <- (!is.null(vars_with_na) &&
        response %in% vars_with_na &&
        !optimise)

      if (dist == "gaussian" && !use_glmm) {
        if (!optimise) {
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
        # GLMM covariance for binomial error term
        if (!optimise) {
          # Marginal Formulation
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
        }
        # When optimise=TRUE, skip - using random effects instead
      } else if (dist == "multinomial") {
        # Multinomial covariance
        # We need TAU[,,k] for each k
        K_var <- paste0("K_", response)

        if (!optimise) {
          # k=1 is reference category (fixed to identity)
          # k>=2 have estimated phylogenetic signal
          if (multi.tree) {
            model_lines <- c(
              model_lines,
              paste0("  # Covariance matrices for multinomial"),
              "  # Reference category k=1",
              "  for (i in 1:N) {",
              "    for (j in 1:N) {",
              paste0("      Mlam_", response, "[i,j,1] <- ID[i,j]"),
              "    }",
              "  }",
              paste0(
                "  TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N,1:N,1] <- ID[1:N,1:N]"
              ),
              "  # Estimated categories k>=2",
              paste0("  for (k in 2:", K_var, ") {"),
              "    for (i in 1:N) {",
              "      for (j in 1:N) {",
              paste0(
                "        Mlam_",
                response,
                "[i,j,k] <- lambda_",
                response,
                "[k]*multiVCV[i,j,K] + (1-lambda_",
                response,
                "[k])*ID[i,j]"
              ),
              "      }",
              "    }",
              paste0(
                "    TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N,1:N,k] <- inverse(Mlam_",
                response,
                "[,,k])"
              ),
              "  }"
            )
          } else {
            model_lines <- c(
              model_lines,
              paste0("  # Covariance matrices for multinomial"),
              "  # Reference category k=1",
              paste0(
                "  TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N,1:N,1] <- ID[1:N,1:N]"
              ),
              "  # Estimated categories k>=2",
              paste0("  for (k in 2:", K_var, ") {"),
              paste0(
                "    TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N, 1:N, k] <- inverse(lambda_",
                response,
                "[k] * VCV + (1 - lambda_",
                response,
                "[k]) * ID)"
              ),
              "  }"
            )
          }
        }
      }
    }
  }

  # Covariance for correlated vars (phylogenetic part)
  # Only use VCV approach when optimise=FALSE; optimised models use eigendecomposition
  if (!is.null(induced_correlations) && !optimise) {
    for (var in correlated_vars) {
      if (multi.tree) {
        model_lines <- c(
          model_lines,
          paste0(
            "  TAU_phylo_",
            var,
            " <- tau_phylo_",
            var,
            " * inverse(multiVCV[,,K])"
          )
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  TAU_phylo_", var, " <- tau_phylo_", var, " * inverse(VCV)")
        )
      }
    }
  }

  # Imputation priors for predictors (those not modeled as responses)
  model_lines <- c(model_lines, "  # Predictor priors for imputation")
  non_response_vars <- setdiff(all_vars, names(response_counter))
  for (var in non_response_vars) {
    # Check if this is a latent variable
    is_latent <- !is.null(latent) && var %in% latent

    model_lines <- c(
      model_lines,
      paste0("  for (i in 1:N) {"),
      paste0("    mu", var, "[i] <- 0"),
      paste0("  }")
    )

    if (independent) {
      # Independent imputation (i.i.d normal)
      if (is_latent && standardize_latent) {
        # Use N(0,1) prior for standardized latent variable
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], tau_e_", var, ")"),
          paste0("  }")
        )
      }
    } else if (optimise) {
      # Optimized Random Effects Formulation for Predictors (Additive)

      # If latent variable with standardize_latent = TRUE, use simple N(0,1) prior
      # ignoring structure (assumes latent is standardized white noise)
      if (is_latent && standardize_latent) {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        # Standard random effects formulation
        additive_terms <- ""

        for (s_name in structure_names) {
          s_suffix <- if (length(structure_names) > 1) {
            paste0("_", s_name)
          } else {
            ""
          }

          u_std <- paste0("u_std_", var, s_suffix)
          u <- paste0("u_", var, s_suffix)
          tau_u <- paste0("tau_u_", var, s_suffix)

          prec_name <- paste0("Prec_", s_name)
          prec_idx <- paste0(prec_name, "[1:N, 1:N]")
          if (multi.tree && s_name == "phylo") {
            prec_idx <- paste0(prec_name, "[1:N, 1:N, K]")
          }

          model_lines <- c(
            model_lines,
            paste0("  ", u_std, "[1:N] ~ dmnorm(zeros[1:N], ", prec_idx, ")"),
            paste0(
              "  for (i in 1:N) { ",
              u,
              "[i] <- ",
              u_std,
              "[i] / sqrt(",
              tau_u,
              ") }"
            )
          )
          additive_terms <- paste0(additive_terms, " + ", u, "[i]")
        }

        tau_e <- paste0("tau_e_", var)

        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(mu",
            var,
            "[i]",
            additive_terms,
            ", ",
            tau_e,
            ")"
          ),
          paste0("  }")
        )
      }
    } else {
      # Standard MVN (Marginal)
      # If latent variable with standardize_latent = TRUE, use simple N(0,1) prior
      if (is_latent && standardize_latent) {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0(
            "  ",
            var,
            "[1:N] ~ dmnorm(mu",
            var,
            "[1:N], TAU",
            tolower(var),
            ")"
          )
        )
      }
    }

    # For latent variables, fix tau = 1 (standardize)
    # For observed predictors, estimate tau
    # Skip if standardize_latent = TRUE (variance already set in N(0,1) prior)
    if (is_latent && !standardize_latent) {
      if (optimise) {
        model_lines <- c(
          model_lines,
          paste0("  # Latent variable: standardized (var = 1)"),
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          # In random effects: Var = (1/tau_u)*V + (1/tau_e)*I
          # We want Var = lambda*V + (1-lambda)*I
          # So 1/tau_u = lambda => tau_u = 1/lambda
          # And 1/tau_e = 1-lambda => tau_e = 1/(1-lambda)
          paste0("  tau_u_", var, " <- 1/lambda", var),
          paste0("  tau_e_", var, " <- 1/(1-lambda", var, ")"),
          paste0("  sigma", var, " <- 1") # Fixed to 1
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  # Latent variable: standardized (var = 1)"),
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          paste0("  tau", var, " <- 1  # Fixed for identification"),
          paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
        )
      }
    } else if (is_latent && standardize_latent) {
      # No variance parameters needed - already specified in N(0,1) prior
      model_lines <- c(
        model_lines,
        paste0(
          "  # Latent variable: fully standardized with N(0,1) prior (no variance parameters)"
        )
      )
    } else {
      if (independent) {
        # Independent Predictor Prior
        model_lines <- c(
          model_lines,
          paste0("  tau_e_", var, " ~ dgamma(1, 1)"),
          paste0("  sigma", var, " <- 1/sqrt(tau_e_", var, ")")
        )
      } else if (optimise) {
        if (length(structure_names) > 1) {
          # Multiple Structures: Estimate independent variance components
          model_lines <- c(
            model_lines,
            paste0("  tau_e_", var, " ~ dgamma(1, 1)")
          )

          for (s_name in structure_names) {
            s_suffix <- paste0("_", s_name)
            tau_u <- paste0("tau_u_", var, s_suffix)
            model_lines <- c(
              model_lines,
              paste0("  ", tau_u, " ~ dgamma(1, 1)"),
              paste0("  sigma_", var, s_suffix, " <- 1/sqrt(", tau_u, ")")
            )
          }
          model_lines <- c(
            model_lines,
            paste0("  sigma_", var, "_res <- 1/sqrt(tau_e_", var, ")")
          )
        } else {
          # Single Structure (Legacy behavior with lambda partitioning)
          # This generates tau_u_var and tau_e_var matching the single-structure Gaussian block
          model_lines <- c(
            model_lines,
            paste0("  lambda", var, " ~ dunif(0, 1)"),
            paste0("  tau", var, " ~ dgamma(1, 1)"),
            paste0("  tau_u_", var, " <- tau", var, "/lambda", var),
            paste0("  tau_e_", var, " <- tau", var, "/(1-lambda", var, ")"),
            paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
          )
        }
      } else {
        model_lines <- c(
          model_lines,
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          paste0("  tau", var, " ~ dgamma(1, 1)"),
          paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
        )
      }
    }

    # Only generate TAU matrix with VCV when there is a structure defined
    if (!optimise && length(structure_names) > 0) {
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
          paste0(
            "  TAU",
            tolower(var),
            " <- tau",
            var,
            "*inverse(Mlam",
            var,
            ")"
          )
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
          paste0(
            "  TAU",
            tolower(var),
            " <- tau",
            var,
            "*inverse(Mlam",
            var,
            ")"
          )
        )
      }
    }
  }

  model_lines <- c(model_lines, "}")
  model_string <- paste(model_lines, collapse = "\n")

  # Convert param_map to data frame
  param_map_df <- do.call(
    rbind,
    lapply(param_map, as.data.frame, stringsAsFactors = FALSE)
  )

  return(list(model = model_string, parameter_map = param_map_df))
}
