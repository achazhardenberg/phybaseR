#' Run a Phylogenetic Bayesian Structural Equation model (PhyBaSE)
#'
#' @param data Data for the model. Accepts:
#'   \itemize{
#'     \item \code{data.frame}: A data frame with variables as columns. Variables needed
#'           for the model are automatically extracted from the equations. Extra columns
#'           are ignored.
#'     \item \code{list}: A named list where each element is a vector of values
#'           (traditional format for backward compatibility).
#'   }
#' @param equations A list of model formulas describing the structural equation model.
#' @param id_col Character string specifying the column name in a data.frame containing
#'   unit identifiers (species, individuals, sites, etc.). This is used to:
#'   \itemize{
#'     \item Match data rows to tree tip labels (for phylogenetic models)
#'     \item Handle repeated measures (multiple rows per unit share the same covariance)
#'   }
#'   If \code{NULL} (default): uses meaningful row names if available.
#'   Ignored when \code{data} is already a list.
#' @param structure The covariance structure for the model. Accepts:
#'   \itemize{
#'     \item \code{"phylo"} object: Phylogenetic tree (Standard PGLS/PhyloSEM).
#'     \item \code{"multiPhylo"} object: List of trees (incorporates phylogenetic uncertainty).
#'     \item \code{NULL}: Independent model (Standard SEM, no covariance structure).
#'     \item \code{matrix}: Custom covariance or precision matrix (e.g., spatial connectivity, kinship).
#'   }
#' @param tree (Deprecated alias for \code{structure}). A single phylogenetic tree of class
#'   \code{"phylo"} or a list of trees. Use \code{structure} instead for new code.
#' @param monitor Parameter monitoring mode. Options:
#'   \itemize{
#'     \item \code{"interpretable"} (default): Monitor only scientifically meaningful parameters:
#'           intercepts (alpha), regression coefficients (beta), phylogenetic signals (lambda) for
#'          responses, and WAIC terms. Excludes variance components (tau) and auxiliary predictor parameters.
#'     \item \code{"all"}: Monitor all model parameters including variance components and implicit equation parameters.
#'     \item Custom vector: Provide a character vector of specific parameter names to monitor.
#'     \item \code{NULL}: Auto-detect based on model structure (equivalent to "interpretable").
#'   }
#' @param n.chains Number of MCMC chains (default = 3).
#' @param n.iter Total number of MCMC iterations (default = 12500).
#' @param n.burnin Number of burn-in iterations (default = n.iter / 5).
#' @param n.thin Thinning rate (default = 10).
#' @param DIC Logical; whether to compute DIC using \code{dic.samples()} (default = TRUE).
#'   **Note**: DIC penalty will be inflated for models with measurement error or repeated measures
#'   because latent variables are counted as parameters (penalty ~ structural parameters + N).
#'   For model comparison, use WAIC or compare mean deviance across models with similar structure.
#' @param WAIC Logical; whether to sample values for WAIC and deviance (default = FALSE).
#'   WAIC is generally more appropriate than DIC for hierarchical models with latent variables.
#' @param n.adapt Number of adaptation iterations (default = n.iter / 5).
#' @param quiet Logical; suppress JAGS output (default = FALSE).
#' @param dsep Logical; if \code{TRUE}, monitor only the first beta in each structural equation (used for d-separation testing).
#' @param variability Optional specification for variables with measurement error or within-species variability.
#'   **AUTO-DETECTION**: If a variable \code{X} has a corresponding column \code{X_se} (standard errors),
#'   \code{X_obs} (repeated measures), or is provided as a matrix, variability is automatically detected.
#'
#'   **Manual specification** (for non-standard column names):
#'   \itemize{
#'     \item Simple: \code{list(X = "se", Y = "reps")} - uses standard \code{X_se}/\code{Y_obs} naming
#'     \item Custom columns: \code{list(X = list(type = "se", se_col = "X_SD"))} - specify custom column names
#'     \item For SE: \code{se_col} (SE column), \code{mean_col} (mean column, optional)
#'     \item For reps: \code{obs_col} (observations matrix column)
#'   }
#' @param distribution Optional named character vector specifying the distribution for response variables.
#'   Default is "gaussian" for all variables. Supported values:
#'   \itemize{
#'     \item "gaussian" (default)
#'     \item "binomial" (binary data)
#'     \item "multinomial" (unordered categorical > 2 levels)
#'     \item "ordinal" (ordered categorical > 2 levels)
#'     \item "poisson" (count data)
#'     \item "negbinomial" (overdispersed count data)
#'   }
#'   Example: \code{distribution = c(Gregarious = "binomial")}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If specified, the model will account for induced correlations among observed
#'   variables that share these latent common causes.
#' @param latent_method Method for handling latent variables (default = "correlations").
#'   \itemize{
#'     \item \code{"correlations"}: MAG approach - marginalize latent variables and estimate
#'           induced correlations (\code{rho}) between observed variables that share latent parents.
#'     \item \code{"explicit"}: Model latent variables as JAGS nodes and estimate structural
#'           paths from latents to observed variables.
#'   }
#' @param standardize_latent Logical; if \code{TRUE} and \code{latent_method = "explicit"},
#'   adds standardized priors (\code{N(0,1)}) to latent variables to identify scale and location.
#'   This improves convergence and makes regression coefficients interpretable as standardized effects.
#'   Only applicable when using explicit latent variable modeling (default = TRUE).
#' @param parallel Logical; if \code{TRUE}, run MCMC chains in parallel (default = FALSE).
#'   Note: Requires \code{n.cores > 1} to take effect.
#' @param n.cores Integer; number of CPU cores to use for parallel chains (default = 1).
#'   Only used when \code{parallel = TRUE}.
#' @param cl Optional; a cluster object created by \code{parallel::makeCluster()}.
#'   If \code{NULL}, a cluster will be created and destroyed automatically.
#' @param ic_recompile Logical; if \code{TRUE} and \code{parallel = TRUE}, recompile the model
#'   after parallel chains to compute DIC/WAIC (default = TRUE).
#'   This adds a small sequential overhead but enables information criteria calculation.
#' @param optimise Logical; if \code{TRUE} (default), use the optimized random effects formulation
#'   for phylogenetic models. This is significantly faster (5-10x) and more numerically stable.
#'   If \code{FALSE}, use the traditional marginal formulation (slower, but provided for comparison).
#'
#' @return A list of class \code{"phybase"} with model output and diagnostics.
#' @export
#' @importFrom ape vcv.phylo branching.times
#' @importFrom rjags jags.model coda.samples dic.samples jags.samples
#' @importFrom stats na.omit update formula terms setNames start var
#' @importFrom utils capture.output
#' @importFrom coda gelman.diag effectiveSize
#' @import coda
phybase_run <- function(
  data,
  equations,
  id_col = NULL,
  structure = NULL,
  tree = NULL,
  monitor = NULL,
  n.chains = 3,
  n.iter = 12500,
  n.burnin = n.iter / 5,
  n.thin = 10,
  DIC = TRUE,
  WAIC = FALSE,
  n.adapt = n.iter / 5,
  quiet = FALSE,
  dsep = FALSE,
  variability = NULL,
  distribution = NULL,
  latent = NULL,
  latent_method = c("correlations", "explicit"),
  standardize_latent = TRUE,
  parallel = FALSE,
  n.cores = 1,
  cl = NULL,
  ic_recompile = TRUE,
  optimise = TRUE
) {
  # Handle 'structure' alias
  if (is.null(tree) && !is.null(structure)) {
    tree <- structure
  }

  # Check if both are missing (tree is NULL, structure is NULL)
  # This implies Independent Model (tree = NULL is valid for that)
  # But we want to be explicit?
  # If both are NULL, we run as independent model.

  if (is.null(tree)) {
    if (!quiet) {
      message(
        "No structure/tree provided. Running standard (non-phylogenetic) SEM."
      )
    }
  }

  # Validate inputs
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(tree) && !quiet) {
    message("No tree provided. Running standard (non-phylogenetic) SEM.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
  }

  # --- Data Frame Preprocessing ---
  # If data is a data.frame, convert to list format expected by the model
  original_data <- data
  if (is.data.frame(data)) {
    # Extract all variable names from equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Check which variables are in the data frame
    available_vars <- intersect(eq_vars, colnames(data))
    missing_vars <- setdiff(eq_vars, colnames(data))

    # Some "missing" vars might be latent - that's OK
    if (!is.null(latent)) {
      missing_vars <- setdiff(missing_vars, latent)
    }

    if (length(missing_vars) > 0 && length(available_vars) == 0) {
      stop(
        "None of the variables in equations found in data frame. ",
        "Missing: ",
        paste(missing_vars, collapse = ", ")
      )
    }

    if (length(missing_vars) > 0 && !quiet) {
      message(
        "Note: Variables not in data (may be latent/derived): ",
        paste(missing_vars, collapse = ", ")
      )
    }

    # Handle id_col for matching to tree/structure
    row_ids <- NULL
    if (!is.null(id_col)) {
      if (!id_col %in% colnames(data)) {
        stop("id_col '", id_col, "' not found in data frame columns.")
      }
      row_ids <- data[[id_col]]
      # Remove id_col from variables to include (it's metadata, not a model variable)
      available_vars <- setdiff(available_vars, id_col)
    } else {
      # Try to use row names if they're meaningful (not just 1, 2, 3...)
      rn <- rownames(data)
      if (!is.null(rn) && !all(rn == as.character(seq_len(nrow(data))))) {
        row_ids <- rn
      }
    }

    # Also check for variability-related columns (X_se, X_obs patterns)
    se_cols <- grep("_se$", colnames(data), value = TRUE)
    obs_cols <- grep("_obs$", colnames(data), value = TRUE)
    extra_cols <- c(se_cols, obs_cols)

    # Convert to list format
    data_list <- list()
    for (var in c(available_vars, extra_cols)) {
      if (var %in% colnames(original_data)) {
        data_list[[var]] <- original_data[[var]]
      }
    }

    # Set names on vectors if we have row_ids
    if (!is.null(row_ids)) {
      for (var in names(data_list)) {
        if (
          is.vector(data_list[[var]]) &&
            length(data_list[[var]]) == length(row_ids)
        ) {
          names(data_list[[var]]) <- row_ids
        }
      }
    }

    # Preserve any existing attributes
    data_attrs <- attributes(original_data)
    data <- data_list

    if (!quiet) {
      message(
        "Converted data.frame to list with ",
        length(data),
        " variables: ",
        paste(names(data), collapse = ", ")
      )
    }
  } else {
    # Ensure data is a list (crucial for adding matrices like VCV)
    # Preserve attributes (like categorical_vars) which are lost during as.list()
    data_attrs <- attributes(data)
    data <- as.list(data)
  }

  # Restore categorical_vars if present
  if ("categorical_vars" %in% names(data_attrs)) {
    attr(data, "categorical_vars") <- data_attrs$categorical_vars
  }

  # --- Structure Processing ---
  structures <- list()
  is_multiple <- FALSE
  N <- NULL

  # 1. Normalize Input to List
  if (is.null(tree)) {
    # Independent
  } else if (
    inherits(tree, "multiPhylo") ||
      (is.list(tree) && all(sapply(tree, inherits, "phylo")))
  ) {
    structures[["phylo"]] <- tree
    is_multiple <- TRUE
  } else if (inherits(tree, "phylo")) {
    structures[["phylo"]] <- tree
  } else if (is.matrix(tree)) {
    structures[["custom"]] <- tree
  } else if (is.list(tree)) {
    structures <- tree
    # Safety Check for list
    if (any(sapply(structures, inherits, "multiPhylo"))) is_multiple <- TRUE
  } else {
    stop("Unknown structure format. Must be tree, matrix, or list of them.")
  }

  structure_names <- names(structures)
  if (is.null(structure_names) && length(structures) > 0) {
    structure_names <- paste0("Struct", seq_along(structures))
    names(structures) <- structure_names
  }

  # 2. Process Structures
  if (length(structures) == 0) {
    # Independent Logic: Determine N from data
    potential_vectors <- Filter(function(x) is.vector(x) || is.factor(x), data)
    if (length(potential_vectors) > 0) {
      N <- length(potential_vectors[[1]])
    } else {
      stop(
        "Could not determine N from data. Please provide structure or vector data."
      )
    }
  } else {
    # Structured Logic

    # Helper to standardize tree
    standardize_tree <- function(phylo_tree) {
      max_bt <- max(ape::branching.times(phylo_tree))
      if (abs(max_bt - 1.0) > 0.01) {
        if (!quiet) {
          message(sprintf("Standardizing tree (max_bt: %.2f -> 1.0)", max_bt))
        }
        phylo_tree$edge.length <- phylo_tree$edge.length / max_bt
      }
      return(phylo_tree)
    }

    for (s_name in structure_names) {
      obj <- structures[[s_name]]

      if (is_multiple && s_name == "phylo") {
        # Multi-Tree (Legacy)
        obj <- lapply(obj, standardize_tree)
        K_tree <- length(obj)
        N_tree <- length(obj[[1]]$tip.label)

        if (is.null(N)) {
          N <- N_tree
        } else if (N != N_tree) {
          stop("Dimension mismatch in multi-tree")
        }

        # Compute Array of Precisions
        Prec_array <- array(0, dim = c(N, N, K_tree))
        for (k in 1:K_tree) {
          Prec_array[,, k] <- solve(ape::vcv.phylo(obj[[k]]))
        }

        data[[paste0("Prec_", s_name)]] <- Prec_array
        data$Ntree <- K_tree

        # Legacy VCV/multiVCV required for some logic?
        # If 'optimise=FALSE', phybase_model might use VCV.
        # But we are moving towards requiring optimise=TRUE for multiple structures.
        # For legacy single multi-tree, we might want to keep 'multiVCV' for safe Measure.
        if (length(structures) == 1) {
          multiVCV <- array(0, dim = c(N, N, K_tree))
          for (k in 1:K_tree) {
            multiVCV[,, k] <- ape::vcv.phylo(obj[[k]])
          }
          data$multiVCV <- multiVCV
        }
      } else if (inherits(obj, "phylo")) {
        obj <- standardize_tree(obj)
        V <- ape::vcv.phylo(obj)
        P <- solve(V)
        data[[paste0("Prec_", s_name)]] <- P

        if (is.null(N)) {
          N <- nrow(V)
        } else if (nrow(V) != N) {
          stop(paste("Dimension mismatch in", s_name))
        }

        if (length(structures) == 1 && optimise == FALSE) {
          data$VCV <- V
        }
      } else if (is.matrix(obj)) {
        V <- obj
        P <- solve(V)
        data[[paste0("Prec_", s_name)]] <- P

        if (is.null(N)) {
          N <- nrow(V)
        } else if (nrow(V) != N) {
          stop(paste("Dimension mismatch in", s_name))
        }
      }
    }

    data$zeros <- rep(0, N)
  }

  data$ID <- diag(N)
  data$ID2 <- diag(2)
  data$N <- N

  # Handle multinomial and ordinal data
  if (!is.null(distribution)) {
    for (var in names(distribution)) {
      if (distribution[[var]] %in% c("multinomial", "ordinal")) {
        if (!var %in% names(data)) {
          stop(paste(
            distribution[[var]],
            "variable",
            var,
            "not found in data."
          ))
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

        if (distribution[[var]] == "multinomial" && K < 3) {
          warning(paste(
            "Multinomial variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        if (distribution[[var]] == "ordinal" && K < 3) {
          warning(paste(
            "Ordinal variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        # Pass K to JAGS (auto-detected from data)
        K_name <- paste0("K_", var)
        if (!K_name %in% names(data)) {
          data[[K_name]] <- K
          if (!quiet) {
            message(sprintf(
              "Auto-detected K_%s = %d from %s variable '%s'",
              var,
              K,
              distribution[[var]],
              var
            ))
          }
        }
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

  # Auto-detect variability from data column names (user-friendly)
  # Look for patterns: X_se, X_obs or matrix columns
  auto_variability <- list()

  for (var in all_vars) {
    # Skip if already in manual variability specification
    if (!is.null(variability) && var %in% c(names(variability), variability)) {
      next
    }

    # Check for SE pattern (X_se)
    se_name <- paste0(var, "_se")
    if (se_name %in% names(data)) {
      auto_variability[[var]] <- "se"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has standard errors in '%s'",
          var,
          se_name
        ))
      }
    }

    # Check for repeated measures pattern (X_obs or matrix)
    obs_name <- paste0(var, "_obs")
    if (var %in% names(data) && is.matrix(data[[var]])) {
      auto_variability[[var]] <- "reps"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has repeated measures (matrix format)",
          var
        ))
      }
    } else if (obs_name %in% names(data)) {
      auto_variability[[var]] <- "reps"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has repeated measures in '%s'",
          var,
          obs_name
        ))
      }
    }
  }

  # Merge auto-detected with manual specification (manual takes precedence)
  if (length(auto_variability) > 0) {
    if (is.null(variability)) {
      variability <- auto_variability
    } else {
      # Convert variability to named list if needed
      if (is.null(names(variability))) {
        variability <- setNames(rep(NA, length(variability)), variability)
      }
      # Merge: manual overrides auto
      for (var in names(auto_variability)) {
        if (!var %in% names(variability)) {
          variability[[var]] <- auto_variability[[var]]
        }
      }
    }
  }

  # Handle variability data
  variability_list <- list()
  if (!is.null(variability)) {
    for (var_name in names(variability)) {
      var_spec <- variability[[var_name]]

      # Parse specification: can be "se"/"reps" or list(type="se", se_col="X_SD")
      if (is.list(var_spec)) {
        # Extended format with custom column names
        type <- var_spec$type
        custom_se_col <- var_spec$se_col
        custom_obs_col <- var_spec$obs_col
        custom_mean_col <- var_spec$mean_col
      } else {
        # Simple format: just the type
        type <- as.character(var_spec)
        custom_se_col <- NULL
        custom_obs_col <- NULL
        custom_mean_col <- NULL
      }

      # Validate type
      if (!type %in% c("se", "reps")) {
        stop(paste(
          "Invalid variability type for",
          var_name,
          "- must be 'se' or 'reps', got:",
          type
        ))
      }

      variability_list[[var_name]] <- type

      if (type == "se") {
        # Determine SE column name (custom or standard)
        se_col <- custom_se_col %||% paste0(var_name, "_se")
        mean_col <- custom_mean_col %||% paste0(var_name, "_mean")

        # Check if SE column exists
        if (!se_col %in% names(data)) {
          stop(paste(
            "Variable",
            var_name,
            "specified as 'se' type but column",
            se_col,
            "not found in data."
          ))
        }

        # Handle mean column
        if (!mean_col %in% names(data)) {
          if (var_name %in% names(data)) {
            # Rename var to var_mean
            data[[mean_col]] <- data[[var_name]]
            data[[var_name]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'se' type but neither",
              var_name,
              "nor",
              mean_col,
              "found in data."
            ))
          }
        } else {
          # If both exist, ensure var is removed (it's a latent parameter now)
          if (var_name %in% names(data)) data[[var_name]] <- NULL
        }

        # Rename custom column to standard name if needed
        if (se_col != paste0(var_name, "_se")) {
          data[[paste0(var_name, "_se")]] <- data[[se_col]]
        }
      } else if (type == "reps") {
        # Determine obs column name (custom or standard)
        obs_col <- custom_obs_col %||% paste0(var_name, "_obs")
        nrep_name <- paste0("N_reps_", var_name)

        # Check if obs column exists
        if (!obs_col %in% names(data)) {
          if (var_name %in% names(data) && is.matrix(data[[var_name]])) {
            # Rename var to var_obs
            data[[obs_col]] <- data[[var_name]]
            data[[var_name]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'reps' type but column",
              obs_col,
              "(as matrix) not found in data."
            ))
          }
        }

        # Rename custom column to standard name if needed
        if (obs_col != paste0(var_name, "_obs")) {
          data[[paste0(var_name, "_obs")]] <- data[[obs_col]]
        }

        # Calculate N_reps if not provided
        if (!nrep_name %in% names(data)) {
          mat <- data[[paste0(var_name, "_obs")]]
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
          data[[paste0(var_name, "_obs")]] <- compact_mat
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
    # Force WAIC and DIC off for d-separation testing (not needed for conditional independence tests)
    if (WAIC || DIC) {
      if (!quiet) {
        message(
          "Note: WAIC and DIC are not computed for d-separation tests (not needed for conditional independence testing)."
        )
      }
      WAIC <- FALSE
      DIC <- FALSE
    }

    # Auto-detect latent variables: variables in equations but not in data
    if (is.null(latent)) {
      vars_in_equations <- unique(unlist(lapply(equations, all.vars)))
      vars_in_data <- names(data)

      # Find variables that appear in equations but not in data
      potential_latents <- setdiff(vars_in_equations, vars_in_data)

      if (length(potential_latents) > 0) {
        # Auto-detect latent variables
        latent <- potential_latents

        if (!quiet) {
          message(
            "Auto-detected latent variable(s): ",
            paste(latent, collapse = ", "),
            "\n(Variables in equations but not in data will be treated as latent.)\n",
            "Generating m-separation tests for MAG..."
          )
        }
      }
    }

    if (!quiet) {
      if (!is.null(latent)) {
        message("Generating m-separation tests (MAG with latent variables)...")
      } else {
        message("Generating d-separation tests...")
      }
    }

    dsep_result <- phybase_dsep(equations, latent = latent, quiet = !dsep)

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
      stop(
        "No d-separation tests implied by the model (model is saturated). Stopping run."
      )
    }

    # Run tests sequentially to avoid cyclic dependencies in JAGS
    # Decide whether to run tests in parallel
    use_parallel <- parallel && n.cores > 1 && length(dsep_tests) > 1

    if (!quiet) {
      if (use_parallel) {
        message(sprintf(
          "Running %d d-separation tests in parallel on %d cores...",
          length(dsep_tests),
          n.cores
        ))
      } else {
        message(sprintf(
          "Running %d d-separation tests sequentially...",
          length(dsep_tests)
        ))
      }
    }

    combined_samples <- NULL
    combined_map <- NULL

    # Define function to run a single d-sep test
    run_single_dsep_test <- function(i, test_eq, monitor_params) {
      # Run model for this single test
      # We pass dsep=FALSE to treat it as a standard model run
      # We pass parallel=FALSE to avoid nested parallelism
      fit <- phybase_run(
        data = data,
        tree = tree,
        equations = list(test_eq), # Pass as list of 1 equation
        monitor = monitor_params,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        DIC = FALSE, # DIC not needed for d-sep tests
        WAIC = FALSE,
        n.adapt = n.adapt,
        quiet = TRUE, # Suppress output for individual runs
        dsep = FALSE,
        variability = variability,
        distribution = distribution,
        latent = NULL, # D-sep tests are on observed variables only
        latent_method = latent_method,
        parallel = FALSE, # Disable nested parallelism
        n.cores = 1,
        cl = NULL,
        ic_recompile = ic_recompile
      )

      # Extract samples, map, and model
      samples <- fit$samples
      param_map <- fit$parameter_map
      model_string <- fit$model

      # Update equation index in parameter map to match the d-sep test index
      param_map$equation_index <- i

      list(
        samples = samples,
        param_map = param_map,
        model = model_string,
        test_index = i
      )
    }

    # Run tests (parallel or sequential)
    if (use_parallel) {
      # Setup cluster if not provided
      if (is.null(cl)) {
        cl <- parallel::makeCluster(n.cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }

      # Export necessary objects to cluster
      parallel::clusterExport(
        cl,
        c(
          "data",
          "tree",
          "monitor",
          "n.chains",
          "n.iter",
          "n.burnin",
          "n.thin",
          "n.adapt",
          "variability",
          "distribution",
          "latent",
          "latent_method",
          "ic_recompile"
        ),
        envir = environment()
      )

      # Run tests in parallel
      results <- parallel::parLapply(
        cl,
        seq_along(dsep_tests),
        function(i) {
          test_eq <- dsep_tests[[i]]
          # Monitor betas for ALL predictors in the d-sep equation
          # This ensures we capture the relevant test statistic regardless of variable order
          test_vars <- all.vars(test_eq)
          response <- as.character(test_eq)[2]
          predictors <- test_vars[test_vars != response]

          params_to_monitor <- character(0)

          for (predictor in predictors) {
            # Check if predictor is categorical
            is_categorical <- FALSE
            if (!is.null(attr(data, "categorical_vars"))) {
              cat_vars <- attr(data, "categorical_vars")
              if (predictor %in% names(cat_vars)) {
                is_categorical <- TRUE
                dummies <- cat_vars[[predictor]]$dummies
                # Monitor all dummy betas
                # For Gaussian: betaPredictor format
                params_to_monitor <- c(
                  params_to_monitor,
                  paste0("beta", dummies)
                )
              }
            }

            if (!is_categorical) {
              # Standard continuous predictor
              params_to_monitor <- c(
                params_to_monitor,
                paste0("beta_", response, "_", predictor)
              )
            }
          }

          current_monitor <- c(monitor, params_to_monitor) # Combine with global monitor

          # Run model for this d-sep test
          if (!quiet) {
            message(sprintf("  Testing: %s", deparse(test_eq)))
            message(sprintf(
              "    Monitoring: %s",
              paste(params_to_monitor, collapse = ", ")
            ))
          }
          run_single_dsep_test(i, test_eq, current_monitor) # Pass current_monitor
        }
      )
    } else {
      # Sequential execution
      results <- lapply(seq_along(dsep_tests), function(i) {
        test_eq <- dsep_tests[[i]]
        # Monitor betas for ALL predictors in the d-sep equation
        # This ensures we capture the relevant test statistic regardless of variable order
        test_vars <- all.vars(test_eq)
        response <- as.character(test_eq)[2]
        predictors <- test_vars[test_vars != response]

        params_to_monitor <- character(0)

        for (predictor in predictors) {
          # Check if predictor is categorical
          is_categorical <- FALSE
          if (!is.null(attr(data, "categorical_vars"))) {
            cat_vars <- attr(data, "categorical_vars")
            if (predictor %in% names(cat_vars)) {
              is_categorical <- TRUE
              dummies <- cat_vars[[predictor]]$dummies
              # Monitor all dummy betas
              # For Gaussian: betaPredictor format
              params_to_monitor <- c(params_to_monitor, paste0("beta", dummies))
            }
          }

          if (!is_categorical) {
            # Standard continuous predictor
            params_to_monitor <- c(
              params_to_monitor,
              paste0("beta_", response, "_", predictor)
            )
          }
        }

        current_monitor <- c(monitor, params_to_monitor) # Combine with global monitor

        if (!quiet) {
          message(sprintf(
            "  Test %d/%d: %s",
            i,
            length(dsep_tests),
            deparse(test_eq)
          ))
          message(sprintf(
            "    Monitoring: %s",
            paste(params_to_monitor, collapse = ", ")
          ))
        }
        run_single_dsep_test(i, test_eq, current_monitor) # Pass current_monitor
      })
    }

    # Combine results
    combined_models <- list()

    for (result in results) {
      samples <- result$samples
      param_map <- result$param_map
      model_string <- result$model
      i <- result$test_index

      # Store model for this test
      combined_models[[i]] <- model_string

      # Rename parameters to include equation index to avoid collisions
      # e.g., betaRS becomes betaRS_1 for equation 1, betaRS_2 for equation 2
      for (ch in seq_along(samples)) {
        chain <- samples[[ch]]
        colnames_orig <- colnames(chain)

        # Add suffix _i to all beta, alpha, lambda, tau, rho parameters
        new_colnames <- sapply(colnames_orig, function(name) {
          if (grepl("^(beta|alpha|lambda|tau|rho|sigma)", name)) {
            paste0(name, "_", i)
          } else {
            name
          }
        })

        colnames(chain) <- new_colnames
        samples[[ch]] <- chain
      }

      # Update parameter_map to reflect new names
      if (!is.null(param_map) && nrow(param_map) > 0) {
        param_map$parameter <- paste0(param_map$parameter, "_", i)
      }

      # Combine samples (cbind chains)
      if (is.null(combined_samples)) {
        combined_samples <- samples
      } else {
        # Check if dimensions match
        if (niter(combined_samples) != niter(samples)) {
          stop("MCMC iteration mismatch between d-sep tests")
        }
        # Combine chains: for each chain, cbind the variables
        new_samples <- coda::mcmc.list()
        for (ch in 1:nchain(combined_samples)) {
          # Combine matrices
          mat1 <- combined_samples[[ch]]
          mat2 <- samples[[ch]]
          # All columns from mat2 should be new (due to renaming)
          new_mat <- cbind(mat1, mat2)
          new_samples[[ch]] <- coda::mcmc(
            new_mat,
            start = start(mat1),
            thin = thin(mat1)
          )
        }
        combined_samples <- new_samples
      }

      # Combine parameter maps
      if (is.null(combined_map)) {
        combined_map <- param_map
      } else {
        combined_map <- rbind(combined_map, param_map)
      }
    }

    # Return combined result
    result <- list(
      samples = combined_samples,
      parameter_map = combined_map,
      models = combined_models,
      dsep = TRUE,
      dsep_tests = dsep_tests,
      induced_correlations = induced_cors
    )
    class(result) <- "phybase"
    return(result)
  }

  # Handle latent variable method
  if (!is.null(latent)) {
    latent_method <- match.arg(latent_method)

    # Force MAG approach when doing d-separation testing
    if (dsep && latent_method == "explicit") {
      if (!quiet) {
        message(
          "Note: d-separation testing with latent variables requires MAG approach. ",
          "Using latent_method = 'correlations'."
        )
      }
      latent_method <- "correlations"
    }

    if (latent_method == "correlations") {
      # MAG approach: marginalize latents, use induced correlations
      # If not already computed by dsep, compute now
      if (is.null(induced_cors)) {
        dsep_result <- phybase_dsep(equations, latent = latent, quiet = !dsep)
        induced_cors <- dsep_result$correlations
      }

      # Filter out equations involving latent variables
      # Handle equations involving latent variables
      original_eq_count <- length(equations)
      equations <- lapply(equations, function(eq) {
        vars <- all.vars(eq)
        if (any(vars %in% latent)) {
          # Remove latent variables from the formula
          # We construct a new formula string excluding latent vars
          rhs <- labels(terms(eq))
          keep_terms <- rhs[!rhs %in% latent]

          if (length(keep_terms) == 0) {
            # Becomes intercept-only model
            new_eq <- as.formula(paste(as.character(eq)[2], "~ 1"))
          } else {
            # Keep observed predictors
            new_eq <- as.formula(paste(
              as.character(eq)[2],
              "~",
              paste(keep_terms, collapse = " + ")
            ))
          }
          return(new_eq)
        }
        return(eq)
      })

      # We no longer filter them out, so we don't report "removed equations"
      if (!quiet) {
        message(
          "Using MAG approach: marginalized latent variables from structural equations."
        )
      }

      # Check if variables with induced correlations need intercept-only models
      # This happens when a variable is involved in an induced correlation
      # but is not a response variable in any remaining equation
      if (length(induced_cors) > 0) {
        # Get all variables from induced correlations
        vars_with_correlations <- unique(unlist(induced_cors))

        # Get response variables from remaining equations
        response_vars <- sapply(equations, function(eq) {
          as.character(eq)[2]
        })

        # Find variables that need intercept models
        vars_needing_intercept <- setdiff(vars_with_correlations, response_vars)

        if (length(vars_needing_intercept) > 0) {
          # Create intercept-only models: X ~ 1
          intercept_equations <- lapply(vars_needing_intercept, function(v) {
            as.formula(paste(v, "~ 1"))
          })

          # Add to equations list
          equations <- c(equations, intercept_equations)

          if (!quiet) {
            message(
              "Created intercept-only models for ",
              length(vars_needing_intercept),
              " variable(s) with induced correlations: ",
              paste(vars_needing_intercept, collapse = ", ")
            )
          }
        }
      }

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Estimating ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      # Explicit approach: keep all equations, don't use induced correlations
      induced_cors <- NULL

      if (!quiet) {
        message("Using explicit latent variable modeling")
      }
    }
  }

  # Auto-expand categorical variables in equations
  if (!is.null(attr(data, "categorical_vars"))) {
    categorical_vars <- attr(data, "categorical_vars")

    equations <- lapply(equations, function(eq) {
      # Parse formula to get all variables
      vars <- all.vars(eq)

      # Check if any predictors are categorical
      for (var in vars) {
        if (var %in% names(categorical_vars)) {
          # Get dummy variable names
          dummies <- categorical_vars[[var]]$dummies

          # Convert formula to character for manipulation
          eq_str <- deparse(eq)

          # Replace categorical variable with its dummies
          # Match whole word only (avoid partial matches)
          pattern <- paste0("\\b", var, "\\b")
          replacement <- paste(dummies, collapse = " + ")
          eq_str <- gsub(pattern, replacement, eq_str)

          # Convert back to formula
          eq <- as.formula(eq_str)

          if (!quiet) {
            message(sprintf(
              "Expanded '%s' to: %s",
              var,
              paste(dummies, collapse = ", ")
            ))
          }
        }
      }
      return(eq)
    })
  }

  # JAGS model code
  model_output <- phybase_model(
    equations = equations,
    multi.tree = is_multiple,
    variability = variability_list,
    distribution = distribution,
    vars_with_na = response_vars_with_na,
    induced_correlations = induced_cors,
    latent = latent,
    standardize_latent = standardize_latent,
    optimise = optimise,
    structure_names = structure_names
  )

  model_string <- model_output$model
  parameter_map <- model_output$parameter_map

  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # Monitor parameters
  # Handle monitor mode
  monitor_mode <- NULL
  if (
    is.character(monitor) &&
      length(monitor) == 1 &&
      monitor %in% c("interpretable", "all")
  ) {
    monitor_mode <- monitor
    monitor <- NULL # Will be auto-detected based on mode
  }

  # Default mode is "interpretable"
  if (is.null(monitor_mode) && is.null(monitor)) {
    monitor_mode <- "interpretable"
  }

  if (is.null(monitor)) {
    lines <- unlist(strsplit(model_string, "\n"))

    extract_names <- function(pattern) {
      out <- grep(pattern, lines, value = TRUE)
      out <- grep("(<-|~)", out, value = TRUE)
      # Extract parameter name, ignoring array indices like [k]
      # This regex captures only the base variable name
      gsub("^\\s*(\\w+)(?:\\[.*\\])?\\s*(<-|~).*", "\\1", out)
    }

    # Extract all parameters
    all_params <- unique(c(
      extract_names("^\\s*beta"),
      extract_names("^\\s*alpha"),
      extract_names("^\\s*lambda"),
      extract_names("^\\s*tau"),
      extract_names("^\\s*rho"),
      extract_names("^\\s*sigma")
    ))

    # Remove tau_obs_* (deterministic constants, not stochastic parameters)
    all_params <- all_params[!grepl("^tau_obs", all_params)]

    if (!is.null(monitor_mode) && monitor_mode == "interpretable") {
      # Filter to interpretable parameters only
      # Get response variables to distinguish them from predictors
      response_vars <- unique(sapply(equations, function(eq) {
        all.vars(formula(eq)[[2]])
      }))

      # Get variables that have induced correlations (MAG exogenous variables)
      # These are like auxiliary predictors - their intercepts are not interpretable
      mag_exogenous_vars <- character(0)
      if (!is.null(induced_cors) && length(induced_cors) > 0) {
        mag_exogenous_vars <- unique(unlist(induced_cors))
      }

      # Include:
      # - Alphas for RESPONSE variables only, EXCLUDING MAG exogenous variables
      # - All betas (regression coefficients)
      # - Lambdas for RESPONSE variables only (exclude auxiliary predictor lambdas)
      # - Rhos (induced correlations)
      # Exclude:
      # - All taus (variance components)
      # - Lambdas for predictor-only variables
      # - Alphas for predictor-only variables (implicit X ~ 1 equations)
      # - Alphas for MAG exogenous variables (e.g., BR, BM in BR <-> BM)

      monitor <- all_params[
        (grepl("^alpha", all_params) &
          gsub("^alpha", "", all_params) %in% response_vars &
          !gsub("^alpha", "", all_params) %in% mag_exogenous_vars) | # Response intercepts, excluding MAG exogenous
          grepl("^beta", all_params) | # All regression coefficients
          grepl("^rho", all_params) | # Induced correlations
          grepl("^sigma", all_params) | # Variance components
          grepl("^rho", all_params) | # Induced correlations
          (grepl("^lambda", all_params) &
            gsub("^lambda", "", all_params) %in% response_vars &
            !gsub("^lambda", "", all_params) %in% mag_exogenous_vars) # Response lambdas, excluding MAG exogenous
      ]
    } else {
      # monitor_mode == "all" or NULL: include everything
      monitor <- all_params
    }
  }

  # Add pointwise log-likelihood monitoring if WAIC requested
  # (Future: LOO will also use this)
  if (WAIC) {
    # Extract log_lik parameters from model
    log_lik_params <- unique(c(
      extract_names("^\\s*log_lik")
    ))

    if (length(log_lik_params) > 0) {
      monitor <- unique(c(monitor, log_lik_params))
      if (!quiet) {
        message(
          "Monitoring ",
          length(log_lik_params),
          " pointwise log-likelihood parameter(s) for WAIC"
        )
      }
    }
  }

  # Add response variables
  # Use perl=TRUE for robust regex matching of variable names
  matches <- regmatches(
    model_string,
    gregexpr("\\b([a-zA-Z0-9_]+)\\s*\\[1:N\\]\\s*~", model_string, perl = TRUE)
  )[[1]]

  response_vars <- unique(gsub("\\s*\\[1:N\\]\\s*~", "", matches))
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
      # Explicitly set RNG seed to ensure chains are different
      inits_list <- list(
        .RNG.name = "base::Wichmann-Hill",
        .RNG.seed = 12345 + chain_id
      )

      model <- jags.model(
        model_file,
        data = data,
        inits = inits_list,
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
    if (!quiet) {
      message(sprintf("Sampling %d chains in parallel...", n.chains))
    }

    chain_results <- parallel::parLapply(cl, seq_len(n.chains), function(i) {
      res <- run_single_chain(
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
      return(res)
    })

    if (!quiet) {
      message("All chains completed.")
    }

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

    # Disable DIC/WAIC if only 1 chain (rjags requirement)
    if (n.chains < 2 && (DIC || WAIC)) {
      warning("DIC and WAIC require at least 2 chains. Disabling calculation.")
      DIC <- FALSE
      WAIC <- FALSE
    }

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

  # Explicitly calculate R-hat if multiple chains
  if (n.chains > 1) {
    tryCatch(
      {
        # Manual R-hat calculation to avoid coda::gelman.diag issues
        # with parallel chains and R scoping problems
        n_chains <- length(samples)

        # Use base R colnames to get parameter names
        first_chain <- as.matrix(samples[[1]])
        pnames <- colnames(first_chain)
        n_params <- length(pnames)

        psrf <- matrix(NA, nrow = n_params, ncol = 2)
        rownames(psrf) <- pnames
        colnames(psrf) <- c("Point est.", "Upper C.I.")

        # Convert chains to matrices ONCE outside the loop
        chain_matrices <- lapply(samples, as.matrix)

        for (idx in seq_len(n_params)) {
          # Extract column idx from each chain matrix
          vals <- do.call(cbind, lapply(chain_matrices, function(m) m[, idx]))

          # Check for constant chains (variance 0)
          # Use explicit variance calculation to avoid R scoping issues
          chain_vars <- numeric(ncol(vals))
          for (col_idx in seq_len(ncol(vals))) {
            chain_vars[col_idx] <- stats::var(vals[, col_idx])
          }

          if (any(chain_vars < 1e-10)) {
            psrf[idx, 1] <- 1.0 # If constant, Rhat is 1
            next
          }

          # Calculate B/W (Gelman-Rubin statistic)
          n_samples <- nrow(vals)
          chain_means <- colMeans(vals)

          # Between-chain variance
          B <- n_samples * stats::var(chain_means)

          # Within-chain variance
          W <- mean(chain_vars)

          # Estimated variance
          var_plus <- (n_samples - 1) / n_samples * W + B / n_samples

          # R-hat
          rhat <- sqrt(var_plus / W)
          psrf[idx, 1] <- rhat
        }
        # Add R-hat to summary statistics
        # summary(samples) returns a list with 'statistics' and 'quantiles'
        # We want to add R-hat to the statistics matrix

        # Match parameter names
        common_params <- intersect(
          rownames(sum_stats$statistics),
          rownames(psrf)
        )

        if (length(common_params) > 0) {
          sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = NA)
          rhat_col_idx <- which(colnames(sum_stats$statistics) == "Rhat")

          # Use numeric indexing to avoid any strange symbol evaluation
          for (j in seq_along(common_params)) {
            p <- common_params[j]
            row_idx <- which(rownames(sum_stats$statistics) == p)
            psrf_row_idx <- which(rownames(psrf) == p)
            if (length(row_idx) == 1 && length(psrf_row_idx) == 1) {
              sum_stats$statistics[row_idx, rhat_col_idx] <- psrf[
                psrf_row_idx,
                1
              ]
            }
          }
        }
      },
      error = function(e) {
        warning("Could not calculate R-hat: ", e$message)
      }
    )
  }

  # Initialize result object
  result <- list(
    model = model,
    model_code = model_output$model,
    data = data, # Store data for recompilation if needed
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

    # Compute WAIC using pointwise log-likelihoods
    if (WAIC) {
      result$WAIC <- phybase_waic(result)
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
      # WAIC can be computed from pointwise log-likelihoods even with parallel chains
      result$WAIC <- phybase_waic(result)
    }
  } else {
    # Sequential execution - use standard approach
    if (DIC) {
      result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
    }
    # WAIC will be computed after class assignment
  }

  # Assign class before WAIC computation (phybase_waic needs this)
  class(result) <- "phybase"

  # Compute WAIC if requested (must be after class assignment)
  if (WAIC) {
    result$WAIC <- phybase_waic(result)
  }

  return(result)
}
