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
#' @param optimize Logical. If TRUE (default), use random effects formulation for 4.6× speedup.
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
  optimize = TRUE
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
        base <- paste0("beta", pred)
        count <- sum(grepl(paste0("^", base), unlist(beta_counter)))
        beta_name <- if (count == 0) base else paste0(base, count + 1)
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
          if (optimize) {
            # Optimized Random Effects Formulation (4.6x faster)
            # u_std ~ dmnorm(0, Prec_phylo_fixed)
            # u = u_std / sqrt(tau_u)
            # Y ~ dnorm(mu + u, tau_e)

            u_std <- paste0("u_std_", response, suffix)
            u <- paste0("u_", response, suffix)
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
              paste0(
                "  ",
                u_std,
                "[1:N] ~ dmnorm(zeros[1:N], ",
                prec_index,
                ")"
              ),
              paste0("  for (i in 1:N) {"),
              paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i] + ",
                u,
                "[i], ",
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

        if (optimize) {
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

        if (optimize) {
          # Optimized Random Effects Formulation for Multinomial
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

        # Random Effects Formulation (optimize-only)
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
    }
  }

  # Handle Induced Correlations (Latent Variables)
  if (!is.null(induced_correlations)) {
    model_lines <- c(model_lines, "  # Induced Correlations (Latent Variables)")

    # Process each pair
    for (pair in induced_correlations) {
      var1 <- pair[1]
      var2 <- pair[2]

      # We assume only one instance of each variable for now (no repeated measures logic for latents yet)
      suffix1 <- if ((response_counter[[var1]] %||% 0) > 1) "1" else ""
      suffix2 <- if ((response_counter[[var2]] %||% 0) > 1) "1" else ""

      # We model this using the GLMM approach:
      # var = mu + phylo_err + res_err
      # phylo_err is independent (standard phylo model)
      # res_err is correlated between var1 and var2

      # 1. Phylogenetic errors (independent)
      err_phylo1 <- paste0("err_phylo_", var1)
      err_phylo2 <- paste0("err_phylo_", var2)

      model_lines <- c(
        model_lines,
        paste0(
          "  ",
          err_phylo1,
          "[1:N] ~ dmnorm(zero_vec[], TAU_phylo_",
          var1,
          ")"
        ),
        paste0(
          "  ",
          err_phylo2,
          "[1:N] ~ dmnorm(zero_vec[], TAU_phylo_",
          var2,
          ")"
        )
      )

      # 2. Correlated residual errors
      # We need a loop for this as it's i.i.d across species but correlated within species
      res_err <- paste0("err_res_", var1, "_", var2)
      tau_res_matrix <- paste0("TAU_res_", var1, "_", var2)

      model_lines <- c(
        model_lines,
        paste0("  for (i in 1:N) {"),
        paste0(
          "    ",
          res_err,
          "[i, 1:2] ~ dmnorm(zero_vec_2[], ",
          tau_res_matrix,
          "[1:2, 1:2])"
        ),

        # 3. Observation model
        # var[i] ~ dnorm(mu[i] + phylo[i] + res[i], high_precision)
        # Note: If variables have measurement error, we should use that precision.
        # For now, assuming "perfect" observation of the sum components (standard GLMM)
        # But JAGS needs stochastic node for data.
        # Standard trick: var[i] ~ dnorm(mean, tau_obs)
        # Here mean = mu + phylo + res

        # For var1
        paste0(
          "    ",
          var1,
          "[i] ~ dnorm(mu",
          var1,
          suffix1,
          "[i] + ",
          err_phylo1,
          "[i] + ",
          res_err,
          "[i, 1], tau_obs_",
          var1,
          ")"
        ),

        # For var2
        paste0(
          "    ",
          var2,
          "[i] ~ dnorm(mu",
          var2,
          suffix2,
          "[i] + ",
          err_phylo2,
          "[i] + ",
          res_err,
          "[i, 2], tau_obs_",
          var2,
          ")"
        ),
        paste0("  }")
      )

      # Priors for the correlated residuals
      # Construct 2x2 precision matrix from variances and correlation
      model_lines <- c(
        model_lines,
        paste0(
          "  # Priors for correlated residuals between ",
          var1,
          " and ",
          var2
        ),
        paste0("  rho_", var1, "_", var2, " ~ dunif(-1, 1)"),
        paste0("  tau_res_", var1, " ~ dgamma(1, 1)"),
        paste0("  tau_res_", var2, " ~ dgamma(1, 1)"),
        paste0("  sigma_res_", var1, " <- 1/sqrt(tau_res_", var1, ")"),
        paste0("  sigma_res_", var2, " <- 1/sqrt(tau_res_", var2, ")"),

        # Covariance matrix construction
        paste0("  cov_", var1, "_", var2, "[1, 1] <- 1/tau_res_", var1),
        paste0("  cov_", var1, "_", var2, "[2, 2] <- 1/tau_res_", var2),
        paste0(
          "  cov_",
          var1,
          "_",
          var2,
          "[1, 2] <- rho_",
          var1,
          "_",
          var2,
          " * sigma_res_",
          var1,
          " * sigma_res_",
          var2
        ),
        paste0(
          "  cov_",
          var1,
          "_",
          var2,
          "[2, 1] <- cov_",
          var1,
          "_",
          var2,
          "[1, 2]"
        ),

        paste0(
          "  ",
          tau_res_matrix,
          "[1:2, 1:2] <- inverse(cov_",
          var1,
          "_",
          var2,
          "[1:2, 1:2])"
        )
      )

      # Priors for phylogenetic errors (standard)
      model_lines <- c(
        model_lines,
        paste0("  tau_phylo_", var1, " ~ dgamma(1, 1)"),
        paste0("  tau_phylo_", var2, " ~ dgamma(1, 1)")
      )

      # Observation precision (high if no measurement error, or estimated)
      # For now, we fix it high to treat the decomposition as exact-ish
      # Or better: treat it as the residual error if we didn't have the decomposition
      # Actually, in GLMMs, usually: Y = Fixed + Random + Error
      # Here: Y = Mu + Phylo + (Correlated_Res) + (Measurement_Error)
      # If no measurement error, Correlated_Res IS the residual.
      # So we can set tau_obs to be very high (effectively 0 variance)
      # BUT JAGS fails with infinite precision for data.
      # So we use a large value, e.g., 10000 * observed precision?
      # OR: We can just say the residual IS the correlated part.
      # But we defined res[i] ~ dmnorm.
      # So Y[i] IS deterministic given res[i]? No.
      # Let's use a fixed high precision for the "observation" link
      model_lines <- c(
        model_lines,
        paste0("  tau_obs_", var1, " <- 10000"),
        paste0("  tau_obs_", var2, " <- 10000")
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

    # Skip multinomial and ordinal (handled separately)
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial" || dist == "ordinal") {
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
      if (is.null(vars_with_na) || !response %in% vars_with_na) {
        if (optimize) {
          # Random Effects Priors (tau_u, tau_e)
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
      if (optimize) {
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

  # Priors for correlated vars alphas (intercepts)
  for (var in correlated_vars) {
    model_lines <- c(model_lines, paste0("  alpha", var, " ~ dnorm(0, 1.0E-6)"))
  }

  # Priors for regression coefficients
  unique_betas <- unique(unlist(beta_counter))
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
      use_glmm <- (!is.null(vars_with_na) && response %in% vars_with_na)

      if (dist == "gaussian" && !use_glmm) {
        if (!optimize) {
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
        if (!optimize) {
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
        # When optimize=TRUE, skip - using random effects instead
      } else if (dist == "multinomial") {
        # Multinomial covariance
        # We need TAU[,,k] for each k
        K_var <- paste0("K_", response)

        if (!optimize) {
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
  if (!is.null(induced_correlations)) {
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

    if (optimize) {
      # Optimized Random Effects Formulation for Predictors
      u_std <- paste0("u_std_", var)
      u <- paste0("u_", var)
      tau_u <- paste0("tau_u_", var)
      tau_e <- paste0("tau_e_", var)

      # Handle multi-tree: use Prec_phylo_fixed[,,K] instead of Prec_phylo_fixed[,]
      prec_index <- if (multi.tree) {
        "Prec_phylo_fixed[1:N, 1:N, K]"
      } else {
        "Prec_phylo_fixed[1:N, 1:N]"
      }

      model_lines <- c(
        model_lines,
        paste0(
          "  ",
          u_std,
          "[1:N] ~ dmnorm(zeros[1:N], ",
          prec_index,
          ")"
        ),
        paste0("  for (i in 1:N) {"),
        paste0("    ", u, "[i] <- ", u_std, "[i] / sqrt(", tau_u, ")"),
        paste0(
          "    ",
          var,
          "[i] ~ dnorm(mu",
          var,
          "[i] + ",
          u,
          "[i], ",
          tau_e,
          ")"
        ),
        paste0("  }")
      )
    } else {
      # Standard MVN (Marginal)
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

    # For latent variables, fix tau = 1 (standardize)
    # For observed predictors, estimate tau
    # For latent variables, fix tau = 1 (standardize)
    # For observed predictors, estimate tau
    if (is_latent) {
      if (optimize) {
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
    } else {
      if (optimize) {
        model_lines <- c(
          model_lines,
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          paste0("  tau", var, " ~ dgamma(1, 1)"),
          # Transform to random effects parameterization
          # tau_u = tau/lambda
          # tau_e = tau/(1-lambda)
          paste0("  tau_u_", var, " <- tau", var, "/lambda", var),
          paste0("  tau_e_", var, " <- tau", var, "/(1-lambda", var, ")"),
          paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          paste0("  tau", var, " ~ dgamma(1, 1)"),
          paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
        )
      }
    }

    if (!optimize) {
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
