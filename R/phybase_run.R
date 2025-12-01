#' Run a Phylogenetic Bayesian Structural Equation model (PhyBaSE)
#'
#' @param data A named list of data for JAGS.
#' @param tree A single phylogenetic tree of class \code{"phylo"} or a list of trees of class \code{"multiPhylo"} (i.e., a list of \code{"phylo"} objects). If multiple trees are provided, phylogenetic uncertainty is incorporated by sampling among them during model fitting.
#' @param equations A list of model formulas describing the structural equation model.
#' @param monitor Optional character vector of parameter names to monitor. If \code{NULL}, parameters will be selected automatically based on model structure.
#' @param n.chains Number of MCMC chains (default = 3).
#' @param n.iter Total number of MCMC iterations (default = 2000).
#' @param n.burnin Number of burn-in iterations (default = n.iter / 2).
#' @param n.thin Thinning rate (default = max(1, floor((n.iter - n.burnin) / 1000))).
#' @param DIC Logical; whether to compute DIC using \code{dic.samples()} (default = TRUE).
#' @param WAIC Logical; whether to sample values for WAIC and deviance (default = FALSE).
#' @param n.adapt Number of adaptation iterations (default = 100).
#' @param quiet Logical; suppress JAGS output (default = FALSE).
#' @param dsep Logical; if \code{TRUE}, monitor only the first beta in each structural equation (used for d-separation testing).
#' @param variability Optional character vector or named character vector of variable names that have measurement error or within-species variability.
#'   If named, the names should be the variable names and the values should be the type of variability: "se" (for mean and standard error) or "reps" (for repeated measures).
#'   If unnamed, the function attempts to infer the type from the data:
#'   \itemize{
#'     \item If \code{Var_se} exists in data, type is "se".
#'     \item If \code{Var} is a matrix or \code{Var_obs} exists, type is "reps".
#'   }
#' @param distribution Optional named character vector specifying the distribution for response variables.
#'   Default is "gaussian" for all variables. Supported values: "gaussian", "binomial".
#'   Example: \code{distribution = c(Gregarious = "binomial")}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If specified, the model will account for induced correlations among observed
#'   variables that share these latent common causes.
#' @param parallel Logical; if \code{TRUE}, run MCMC chains in parallel (default = FALSE).
#'   Note: Requires \code{n.cores > 1} to take effect.
#' @param n.cores Integer; number of CPU cores to use for parallel chains (default = 1).
#'   Only used when \code{parallel = TRUE}.
#' @param cl Optional; a cluster object created by \code{parallel::makeCluster()}.
#'   If \code{NULL}, a cluster will be created and destroyed automatically.
#' @param ic_recompile Logical; if \code{TRUE} and \code{parallel = TRUE}, recompile the model
#'   after parallel chains to compute DIC/WAIC (default = TRUE).
#'   This adds a small sequential overhead but enables information criteria calculation.
#'
#' @return A list of class \code{"phybase"} with model output and diagnostics.
#' @export
#' @importFrom ape vcv.phylo branching.times
#' @importFrom rjags jags.model coda.samples dic.samples jags.samples
#' @importFrom stats na.omit update formula terms setNames
#' @importFrom coda gelman.diag effectiveSize
#' @import coda
phybase_run <- function(
  data,
  tree,
  equations,
  monitor = NULL,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = floor(n.iter / 2),
  n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
  DIC = TRUE,
  WAIC = FALSE,
  n.adapt = 100,
  quiet = FALSE,
  dsep = FALSE,
  variability = NULL,
  distribution = NULL,
  latent = NULL,
  parallel = FALSE,
  n.cores = 1,
  cl = NULL,
  ic_recompile = TRUE
) {
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(tree)) {
    stop("Argument 'tree' must be provided.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
  }

  # Ensure data is a list (crucial for adding matrices like VCV)
  data <- as.list(data)

  # Handle tree(s)
  is_multiple <- inherits(tree, "multiPhylo") ||
    (is.list(tree) && all(sapply(tree, inherits, "phylo")))

  if (is_multiple) {
    for (i in seq_along(tree)) {
      tree[[i]]$edge.length <- tree[[i]]$edge.length /
        max(ape::branching.times(tree[[i]]))
    }
    multiVCV <- sapply(tree, ape::vcv.phylo, simplify = "array")
    ID <- diag(nrow(multiVCV[,, 1]))
    VCV <- NULL
    N <- dim(multiVCV)[1]
    data$multiVCV <- multiVCV
    data$Ntree <- dim(multiVCV)[3]
  } else {
    tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
    VCV <- ape::vcv.phylo(tree)
    ID <- diag(nrow(VCV))
    N <- nrow(VCV)
    data$VCV <- VCV
  }

  data$ID <- ID
  data$N <- N

  # Handle multinomial data
  if (!is.null(distribution)) {
    for (var in names(distribution)) {
      if (distribution[[var]] == "multinomial") {
        if (!var %in% names(data)) {
          stop(paste("Multinomial variable", var, "not found in data."))
        }

        val <- data[[var]]
        if (is.factor(val)) {
          K <- nlevels(val)
          data[[var]] <- as.integer(val)
        } else {
          # Assume it's already integer or character
          val <- as.factor(val)
          K <- nlevels(val)
          data[[var]] <- as.integer(val)
        }

        if (K < 3) {
          warning(paste(
            "Multinomial variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        # Pass K to JAGS
        data[[paste0("K_", var)]] <- K
      }
    }
  }

  # Check for missing data
  all_vars <- unique(unlist(lapply(equations, function(eq) {
    c(all.vars(eq[[3]]), all.vars(eq[[2]]))
  })))

  response_vars <- unique(sapply(equations, function(eq) all.vars(eq[[2]])))
  predictor_only_vars <- setdiff(all_vars, response_vars)

  # Detect variables with missing data
  response_vars_with_na <- character(0)
  predictor_vars_with_na <- character(0)

  for (var in all_vars) {
    if (var %in% names(data)) {
      var_data <- data[[var]]
      if (
        !is.matrix(var_data) && any(is.na(var_data)) && !all(is.na(var_data))
      ) {
        if (var %in% response_vars) {
          response_vars_with_na <- c(response_vars_with_na, var)
        } else {
          predictor_vars_with_na <- c(predictor_vars_with_na, var)
        }
      }
    }
  }

  # Handle predictor-only variables with missing data
  if (length(predictor_vars_with_na) > 0) {
    if (!quiet) {
      message(
        "Note: Detected missing data in predictor-only variables: ",
        paste(predictor_vars_with_na, collapse = ", "),
        "\nAutomatically adding intercept-only equations (e.g., ",
        predictor_vars_with_na[1],
        " ~ 1) to enable imputation."
      )
    }

    # Add intercept-only equations
    for (var in predictor_vars_with_na) {
      new_eq <- stats::as.formula(paste(var, "~ 1"))
      equations <- c(equations, list(new_eq))
    }

    # Treat them as responses now
    response_vars_with_na <- c(response_vars_with_na, predictor_vars_with_na)
  }

  # For variables with missing data, we'll use the GLMM (Latent Variable) approach
  if (length(response_vars_with_na) > 0 && !quiet) {
    message(
      "Note: Using Latent Variable (GLMM) approach for variables with missing data to preserve phylogenetic signal."
    )
  }

  # Handle variability data
  variability_list <- list()
  if (!is.null(variability)) {
    # If unnamed, infer types
    if (is.null(names(variability))) {
      vars <- variability
      types <- rep(NA, length(vars))
      names(types) <- vars
    } else {
      vars <- names(variability)
      types <- variability
    }

    for (var in vars) {
      type <- types[[var]]

      # Infer type if missing
      if (is.na(type)) {
        if (paste0(var, "_se") %in% names(data)) {
          type <- "se"
        } else if (
          is.matrix(data[[var]]) || paste0(var, "_obs") %in% names(data)
        ) {
          type <- "reps"
        } else {
          stop(paste(
            "Could not infer variability type for",
            var,
            "- provide '_se' for summary stats or matrix for repeated measures."
          ))
        }
        types[[var]] <- type
      }

      variability_list[[var]] <- type

      if (type == "se") {
        # Check if SE exists
        if (!paste0(var, "_se") %in% names(data)) {
          stop(paste(
            "Variable",
            var,
            "specified as 'se' type but",
            paste0(var, "_se"),
            "not found in data."
          ))
        }

        # Handle mean
        mean_name <- paste0(var, "_mean")
        if (!mean_name %in% names(data)) {
          if (var %in% names(data)) {
            # Rename var to var_mean
            data[[mean_name]] <- data[[var]]
            data[[var]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var,
              "specified as 'se' type but neither",
              var,
              "nor",
              mean_name,
              "found in data."
            ))
          }
        } else {
          # If both exist, ensure var is removed (it's a latent parameter now)
          if (var %in% names(data)) data[[var]] <- NULL
        }
      } else if (type == "reps") {
        # Handle repeated measures
        obs_name <- paste0(var, "_obs")
        nrep_name <- paste0("N_reps_", var)

        if (!obs_name %in% names(data)) {
          if (var %in% names(data) && is.matrix(data[[var]])) {
            # Rename var to var_obs
            data[[obs_name]] <- data[[var]]
            data[[var]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var,
              "specified as 'reps' type but neither",
              var,
              "(as matrix) nor",
              obs_name,
              "found in data."
            ))
          }
        }

        # Calculate N_reps if not provided
        if (!nrep_name %in% names(data)) {
          mat <- data[[obs_name]]
          # Count non-NA values per row
          n_reps <- apply(mat, 1, function(x) sum(!is.na(x)))
          data[[nrep_name]] <- n_reps

          # Compact matrix (move non-NAs to left) to ensure 1:N_reps indexing works
          compact_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
          for (i in seq_len(nrow(mat))) {
            vals <- mat[i, !is.na(mat[i, ])]
            if (length(vals) > 0) {
              compact_mat[i, seq_along(vals)] <- vals
            }
          }
          data[[obs_name]] <- compact_mat
        }

        # Ensure var is removed (latent)
        if (var %in% names(data)) data[[var]] <- NULL
      }
    }
  }

  # Handle d-sep logic
  dsep_tests <- NULL
  induced_cors <- NULL

  if (dsep) {
    if (!quiet) {
      if (!is.null(latent)) {
        message("Generating m-separation tests (MAG with latent variables)...")
      } else {
        message("Generating d-separation tests...")
      }
    }

    dsep_result <- phybase_dsep(equations, latent = latent)

    # Extract tests and correlations
    if (!is.null(latent)) {
      dsep_tests <- dsep_result$tests
      induced_cors <- dsep_result$correlations

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Found ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      dsep_tests <- dsep_result
    }

    if (length(dsep_tests) == 0) {
      warning(
        "No d-separation tests implied by the model (model is saturated)."
      )
    } else {
      # Use d-sep tests as equations
      equations <- dsep_tests
    }
  }

  # JAGS model code
  model_output <- phybase_model(
    equations,
    multi.tree = is_multiple,
    variability = variability_list,
    distribution = distribution,
    vars_with_na = response_vars_with_na,
    induced_correlations = induced_cors
  )

  model_string <- model_output$model
  parameter_map <- model_output$parameter_map

  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # Monitor parameters
  if (is.null(monitor)) {
    lines <- unlist(strsplit(model_string, "\n"))

    extract_names <- function(pattern) {
      out <- grep(pattern, lines, value = TRUE)
      out <- grep("(<-|~)", out, value = TRUE)
      # Extract parameter name, including array indices like [k]
      # This regex captures word characters followed by optional array notation
      gsub("^\\s*(\\w+(?:\\[\\w+\\])?)\\s*(<-|~).*", "\\1", out)
    }
    monitor <- unique(c(
      extract_names("^\\s*beta"),
      extract_names("^\\s*alpha"),
      extract_names("^\\s*lambda"),
      extract_names("^\\s*tau"),
      extract_names("^\\s*rho")
    ))
  }

  # Add response variables
  matches <- regmatches(
    model_string,
    gregexpr("\\b(\\w+)\\[1:N\\]\\s*~", model_string)
  )[[1]]
  response_vars <- unique(gsub("\\[1:N\\]\\s*~", "", matches))
  for (v in response_vars) {
    # Skip if variable is in variability list (it's latent, not data)
    if (!is.null(variability) && v %in% names(variability_list)) {
      next
    }

    if (!v %in% names(data)) {
      base <- sub("[0-9]+$", "", v)
      if (base %in% names(data)) data[[v]] <- data[[base]]
    }
  }

  # Run MCMC chains (parallel or sequential)
  if (parallel && n.cores > 1 && n.chains > 1) {
    # Parallel execution
    message(sprintf(
      "Running %d chains in parallel on %d cores...",
      n.chains,
      n.cores
    ))

    # Setup cluster if not provided
    if (is.null(cl)) {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    # Helper function to run a single chain
    run_single_chain <- function(
      chain_id,
      model_file,
      data,
      monitor,
      n.burnin,
      n.iter,
      n.thin,
      n.adapt,
      quiet
    ) {
      # Load rjags in each worker
      if (!requireNamespace("rjags", quietly = TRUE)) {
        stop("Package 'rjags' is required for parallel execution.")
      }
      # Explicitly load rjags to ensure modules are available
      loadNamespace("rjags")
      # Compile model for this chain
      model <- jags.model(
        model_file,
        data = data,
        n.chains = 1,
        n.adapt = n.adapt,
        quiet = quiet
      )

      # Burn-in
      update(model, n.iter = n.burnin)

      # Sample
      samples <- coda.samples(
        model,
        variable.names = monitor,
        n.iter = n.iter - n.burnin,
        thin = n.thin
      )

      return(list(samples = samples, model = model))
    }

    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("run_single_chain"), envir = environment())

    # Run chains in parallel
    chain_results <- parallel::parLapply(cl, seq_len(n.chains), function(i) {
      run_single_chain(
        i,
        model_file,
        data,
        monitor,
        n.burnin,
        n.iter,
        n.thin,
        n.adapt,
        quiet
      )
    })

    # Combine samples from all chains
    samples <- coda::mcmc.list(lapply(chain_results, function(x) {
      x$samples[[1]]
    }))

    # Use the first chain's model for DIC/WAIC (they all have the same structure)
    model <- chain_results[[1]]$model
  } else {
    # Sequential execution (default)
    # Compile model
    model <- rjags::jags.model(
      model_file,
      data = data,
      n.chains = n.chains,
      n.adapt = n.adapt,
      quiet = quiet
    )
    update(model, n.iter = n.burnin)

    # Sample posterior
    samples <- rjags::coda.samples(
      model,
      variable.names = monitor,
      n.iter = n.iter - n.burnin,
      thin = n.thin
    )
  }

  # Summarize posterior
  sum_stats <- summary(samples)

  # Initialize result object
  result <- list(
    model = model,
    model_code = model_output$model,
    samples = samples,
    summary = sum_stats,
    monitor = monitor,
    modfile = model_file,
    dsep = dsep,
    dsep_tests = dsep_tests,
    parameter_map = parameter_map,
    induced_correlations = induced_cors
  )

  # Add DIC and WAIC
  # For parallel runs, recompile the model if ic_recompile=TRUE
  if (
    (DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1 && ic_recompile
  ) {
    message("Recompiling model for DIC/WAIC calculation...")

    # Recompile model with 2 chains for IC calculation (DIC requires >=2)
    ic_model <- rjags::jags.model(
      model_file,
      data = data,
      n.chains = 2,
      n.adapt = n.adapt,
      quiet = quiet
    )

    # Short burn-in (use a fraction of original)
    update(ic_model, n.iter = min(n.burnin, 500))

    # Compute DIC
    if (DIC) {
      result$DIC <- rjags::dic.samples(
        ic_model,
        n.iter = min(n.iter - n.burnin, 1000)
      )
    }

    # Compute WAIC
    if (WAIC) {
      waic_samples <- rjags::jags.samples(
        ic_model,
        c("WAIC", "deviance"),
        type = "mean",
        n.iter = min(n.iter - n.burnin, 1000),
        thin = n.thin
      )
      waic_samples$p_waic <- waic_samples$WAIC
      waic_samples$waic <- waic_samples$deviance + waic_samples$p_waic
      tmp <- sapply(waic_samples, sum)
      result$WAIC <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 1)
    }
  } else if ((DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1) {
    # Parallel without recompilation - warn user
    if (DIC) {
      warning(
        "DIC calculation disabled for parallel chains. Set ic_recompile=TRUE to compute DIC."
      )
      result$DIC <- NULL
    }
    if (WAIC) {
      warning(
        "WAIC calculation disabled for parallel chains. Set ic_recompile=TRUE to compute WAIC."
      )
      result$WAIC <- NULL
    }
  } else {
    # Sequential execution - use standard approach
    if (DIC) {
      result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
    }

    if (WAIC) {
      waic_samples <- rjags::jags.samples(
        model,
        c("WAIC", "deviance"),
        type = "mean",
        n.iter = n.iter - n.burnin,
        thin = n.thin
      )
      waic_samples$p_waic <- waic_samples$WAIC
      waic_samples$waic <- waic_samples$deviance + waic_samples$p_waic
      tmp <- sapply(waic_samples, sum)
      result$WAIC <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 1)
    }
  }

  class(result) <- "phybase"
  return(result)
}
