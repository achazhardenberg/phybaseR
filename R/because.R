#' Run a Bayesian Structural Equation Model
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
#'     \item Link data to external spatial or custom covariance matrices.
#'   }
#'   **Note**: For standard random effects models (e.g. \code{random = ~(1|species)}) where no external structure
#'   (like a tree) is provided, this argument is **not required**. The grouping column is read directly from the data.
#'
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
#'   \strong{Global Setting}:
#'   \itemize{
#'     \item \code{"reps"}: Applies repeat-measures modeling to \strong{all} continuous variables in the equations (except grouping variables). Expects \code{X_obs} matrix or long-format data.
#'     \item \code{"se"}: Applies measurement error modeling to \strong{all} continuous variables. Expects \code{X_se} columns.
#'   }
#'
#'   \strong{Manual Specification} (Named Vector/List):
#'   \itemize{
#'     \item Simple: \code{c(X = "se", Y = "reps")} - mixed types
#'     \item Custom columns: \code{list(X = list(type = "se", se_col = "X_sd"))}
#'     \item For SE: \code{se_col} (SE column), \code{mean_col} (mean column, optional)
#'     \item For reps: \code{obs_col} (observations matrix column)
#'   }
#'
#'   \strong{Auto-Detection}:
#'   If not specified, the package attempts to detect variability based on column names:
#'   \itemize{
#'     \item \code{X_se} -> type="se"
#'     \item \code{X_obs} or matrix column -> type="reps"
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
#'     \item "zip" (zero-inflated poisson): Models excess zeros with probability \code{psi} and counts with mean \code{lambda}.
#'     \item "zinb" (zero-inflated negative binomial): Models excess zeros with probability \code{psi} and overdispersed counts with mean \code{mu} and size \code{r}.
#'   }
#'   The model will estimate a zero-inflation probability parameter \code{psi_Response} for these distributions.
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
#' @param random Optional formula or list of formulas specifying global random effects
#'   applied to all equations (e.g. \code{~(1|species)}).
#' @param levels (Hierarchical Data) A named list mapping variables to their hierarchy levels.
#'   Required if \code{data} is a list of data frames (hierarchical format).
#'   Example: \code{list(individual = c("y", "x"), site = c("z"))}.
#' @param hierarchy (Hierarchical Data) Character string describing the topological ordering of levels
#'   (e.g., \code{"site > individual"}). Required for hierarchical data if not fully inferred from random effects.
#' @param link_vars (Hierarchical Data) Optional named character vector specifying variables used to link
#'   data levels (e.g. \code{c(site = "site_id")}).
#' @param fix_residual_variance Optional named vector for fixing residual variances.
#'   Useful for handling non-identified models or specific theoretical constraints.
#'   Example: \code{c(response_var = 1)}.
#'
#' @return A list of class \code{"because"} with model output and diagnostics.
#' @export
#' @importFrom ape vcv.phylo branching.times
#' @importFrom rjags jags.model coda.samples dic.samples jags.samples
#' @importFrom stats na.omit update formula terms setNames start var
#' @importFrom utils capture.output
#' @importFrom coda gelman.diag effectiveSize
#' @import coda
because <- function(
  equations,
  data,
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
  optimise = TRUE,
  random = NULL, # Global random effects applied to all equations
  levels = NULL, # Hierarchical data: variable-to-level mapping
  hierarchy = NULL, # Hierarchical data: level ordering (e.g., "site > individual")
  link_vars = NULL, # Hierarchical data: variables linking levels
  fix_residual_variance = NULL # Optional: fix residual variance (tau_e) for specific variables
) {
  # --- Input Validation & Setup ---

  # WAIC Validity Check: Conditional WAIC (optimise=FALSE) is misleading
  if (!optimise && WAIC) {
    warning(
      "WAIC is calculated using conditional likelihoods when optimise=FALSE, which produces much higher values than marginal pseudo-likelihoods and is NOT comparable to optimise=TRUE results. Disabling WAIC to avoid confusion. Please use DIC for model comparison in this mode."
    )
    WAIC <- FALSE
  }

  latent_method <- match.arg(latent_method)

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Handle 'structure' alias
  if (is.null(tree) && !is.null(structure)) {
    tree <- structure
  }

  # Check if both are missing (tree is NULL, structure is NULL)
  # This implies Independent Model (tree = NULL is valid for that)

  # If both are NULL, we run as independent model.

  # Tree check moved to line 165

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

  # --- Automatic Data Cleaning (Handle Character/Factor Columns) ---
  data <- preprocess_categorical_vars(data, quiet = quiet)

  # --- Hierarchical Data Detection & Validation ---
  # Data is hierarchical if it's a list (not dataframe) AND levels are provided
  is_hierarchical <- is.list(data) && !is.data.frame(data) && !is.null(levels)
  hierarchical_info <- NULL

  if (is_hierarchical) {
    # Validate hierarchical data structure
    validate_hierarchical_data(data, levels, hierarchy, link_vars)

    # Try to infer hierarchy from random effects if not provided
    if (is.null(hierarchy)) {
      hierarchy <- parse_hierarchy_from_random(random, data)
      if (is.null(hierarchy)) {
        stop(
          "Hierarchical data detected but 'hierarchy' not specified. ",
          "Provide either:\n",
          "  1. 'hierarchy' argument (e.g., \"site_year > individual\"), or\n",
          "  2. Nested random effects (e.g., ~(1|site/individual))"
        )
      }
    }

    # Store hierarchical info for later use
    hierarchical_info <- list(
      data = data,
      levels = levels,
      hierarchy = hierarchy,
      link_vars = link_vars
    )

    if (!quiet) {
      message("Hierarchical data structure detected: ", hierarchy)
    }
  } else {
    # Single-level data - ensure levels/hierarchy not mistakenly provided
    if (!is.null(levels) || !is.null(hierarchy)) {
      warning(
        "'levels' or 'hierarchy' provided but 'data' is not hierarchical. ",
        "Ignoring these arguments."
      )
    }
  }

  # --- Random Effects Parsing ---
  # Extract (1|Group) and update equations to be fixed-effects only
  parsed_random <- extract_random_effects(equations)
  equations <- parsed_random$fixed_equations
  # Start with equation-specific random terms
  random_terms <- parsed_random$random_terms

  # Parse global random argument if provided
  if (!is.null(random)) {
    global_random_terms <- parse_global_random(random, equations)
    # Combine with equation-specific terms
    random_terms <- c(random_terms, global_random_terms)

    # Deduplicate terms (avoid adding same (1|Group) twice for same response)
    # Create unique keys
    if (length(random_terms) > 0) {
      keys <- sapply(random_terms, function(x) {
        paste(x$response, x$group, sep = "|")
      })
      random_terms <- random_terms[!duplicated(keys)]
    }
  }

  # --- Polynomial Term Extraction ---
  # Extract I(var^power) terms and expand formulas
  all_poly_terms <- get_all_polynomial_terms(equations)

  if (!is.null(all_poly_terms)) {
    # Expand formulas to replace I(x^2) with x_pow2
    equations <- lapply(equations, function(eq) {
      poly_terms <- extract_polynomial_terms(eq)
      expand_polynomial_formula(eq, poly_terms)
    })

    if (!quiet) {
      message(
        "Detected ",
        length(all_poly_terms),
        " polynomial term(s): ",
        paste(sapply(all_poly_terms, function(x) x$original), collapse = ", ")
      )
    }

    # Auto-assign polynomial variables to the same level as their base variable
    if (is_hierarchical && !is.null(levels)) {
      for (poly in all_poly_terms) {
        base_var <- poly$base_var
        new_var <- poly$internal_name

        # Find level of base_var
        for (lvl_name in names(levels)) {
          if (base_var %in% levels[[lvl_name]]) {
            # Add the new poly var to this level
            levels[[lvl_name]] <- c(levels[[lvl_name]], new_var)
            break
          }
        }
      }

      # Update the info stored for later data retrieval
      hierarchical_info$levels <- levels
    }
  }
  # Initialize random structures (will be populated later)
  random_structures <- list()
  random_data_updates <- list()

  # Initialize result variables
  dsep_tests <- NULL
  dsep_results <- NULL
  dsep_correlations <- NULL

  # Handle global variability setting (e.g. variability = "reps")
  # If user provides a single string, apply it to all variables in equations
  if (
    !is.null(variability) &&
      is.character(variability) &&
      length(variability) == 1 &&
      is.null(names(variability))
  ) {
    global_type <- variability
    if (global_type %in% c("se", "reps")) {
      message(sprintf(
        "Global variability setting detected: applying '%s' to all variables.",
        global_type
      ))

      # Extract all variables from equations
      all_eq_vars <- unique(unlist(lapply(equations, all.vars)))

      # Exclude grouping variables from random effects
      grouping_vars <- character(0)
      if (!is.null(random)) {
        if (inherits(random, "formula")) {
          random_list <- list(random)
        } else {
          random_list <- random
        }

        for (r in random_list) {
          vars_in_random <- all.vars(r)
          grouping_vars <- c(grouping_vars, vars_in_random)
        }
      }

      # Also exclude Id col if provided
      if (!is.null(id_col)) {
        grouping_vars <- c(grouping_vars, id_col)
      }

      # Variables to apply variability to
      target_vars <- setdiff(all_eq_vars, grouping_vars)

      # Create named vector
      variability <- setNames(
        rep(global_type, length(target_vars)),
        target_vars
      )
    }
  }

  # --- Data Frame Preprocessing ---
  # If data is a data.frame, convert to list format expected by the model
  original_data <- data

  # --- Hierarchical Data Assembly ---
  # If hierarchical data provided, assemble full dataset for main model run
  if (is_hierarchical) {
    # Get all variables needed across all equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add random effect grouping variables
    if (length(random_terms) > 0) {
      random_vars <- unique(sapply(random_terms, function(x) x$group))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

    # Remove latent variables (not in data)
    if (!is.null(latent)) {
      eq_vars <- setdiff(eq_vars, latent)
    }

    # Get assembled dataset with all needed variables
    data <- get_data_for_variables(
      eq_vars,
      hierarchical_info$data,
      hierarchical_info$levels,
      hierarchical_info$hierarchy,
      hierarchical_info$link_vars
    )

    # Update original_data for later use
    original_data <- data

    if (!quiet) {
      message(
        "Assembled hierarchical data: ",
        nrow(data),
        " observations, ",
        ncol(data),
        " variables"
      )
    }
  }

  # --- Random Effects Data Prep (Post-Assembly) ---
  # Create structures for JAGS using the assembled data
  if (length(random_terms) > 0) {
    if (is.null(optimise) & is.null(tree)) {
      optimise <- TRUE
    }

    if (is_hierarchical && !is.data.frame(data)) {
      # Should not happen if assembly worked, but safety check
      stop("Failed to assemble hierarchical data frame for random effects.")
    }

    rand_structs <- create_group_structures(data, random_terms)
    random_structures <- rand_structs$structures
    random_data_updates <- rand_structs$data_updates
  }

  if (is.data.frame(data)) {
    # Check for long format data requiring matrix conversion
    # If any variability specified as 'reps', we attempt to auto-format using because_format_data
    has_reps <- any(grepl("reps", variability))

    if (has_reps) {
      if (missing(tree) || is.null(tree)) {
        # Cannot auto-format without tree to determine species order

        warning(
          "Variability 'reps' specified but no tree provided. Automatic formatting requires a tree to order species rows. Assuming data is already aggregated or user handles index mapping."
        )
      } else if (!is.null(id_col) && id_col %in% names(data)) {
        if (!quiet) {
          message(
            "Detected 'reps' variability and long-format data. Auto-formatting matrices using because_format_data()..."
          )
        }

        # Determine the tree to use (first one if list)
        use_tree <- if (inherits(tree, "multiPhylo")) tree[[1]] else tree

        formatted_list <- because_format_data(
          data,
          species_col = id_col,
          tree = use_tree
        )

        # Now rename the variables that are 'reps' to include '_obs' suffix
        # And keep others as is
        reps_vars <- names(variability)[variability == "reps"]

        final_data_list <- list()
        for (nm in names(formatted_list)) {
          if (nm %in% reps_vars) {
            # This is a matrix of replicates, rename to _obs
            final_data_list[[paste0(nm, "_obs")]] <- formatted_list[[nm]]
          } else {
            # This is a regular variable (vector or matrix depending on because_format_data logic)
            final_data_list[[nm]] <- formatted_list[[nm]]
          }
        }

        # Update data to be this list
        data <- final_data_list

        if (!quiet) {
          message(
            "  Formatted ",
            length(names(formatted_list)),
            " variables as replicate matrices."
          )
        }

        # Check for missing variables that were expected to be formatted
        missing_reps <- setdiff(reps_vars, names(formatted_list))
        if (length(missing_reps) > 0) {
          stop(paste(
            "The following variables were identified for 'reps' processing (from equations) but were NOT found in the data:",
            paste(missing_reps, collapse = ", "),
            "\nPlease check your column names."
          ))
        }
      } else {
        warning(
          "Variability 'reps' specified but 'id_col' missing or not in data. Cannot auto-format long data."
        )
      }
    }
  }

  if (is.data.frame(data)) {
    # Extract all variable names from fixed equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add variables from random terms (grouping factors)
    if (length(random_terms) > 0) {
      random_vars <- unique(sapply(random_terms, function(x) x$group))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

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

    # Add generated dummy variables for categorical predictors
    dummy_vars <- character(0)
    if (!is.null(attr(data, "categorical_vars"))) {
      cat_vars <- attr(data, "categorical_vars")
      dummy_vars <- unlist(lapply(cat_vars, function(x) x$dummies))
    }

    extra_cols <- c(se_cols, obs_cols, dummy_vars)

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

    # Combine data_list with random effect data
    data <- data_list
    if (length(random_data_updates) > 0) {
      for (nm in names(random_data_updates)) {
        if (!is.null(nm) && nm != "") {
          data[[nm]] <- random_data_updates[[nm]]
        }
      }
    }

    # Remove raw grouping variables from data passed to JAGS to avoid "Unused variable" warnings
    if (length(random_terms) > 0) {
      # Use setdiff to avoid error if variable not present (though they should be)
      vars_to_remove <- intersect(names(data), random_vars)
      if (length(vars_to_remove) > 0) {
        data[vars_to_remove] <- NULL
      }
    }

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

  # Compute polynomial values
  # We add them to original_data (for d-sep/residuals) but NOT to 'data' passed to JAGS
  # because JAGS creates deterministic nodes for them (var_pow2 <- var^2)
  if (!is.null(all_poly_terms)) {
    for (poly_term in all_poly_terms) {
      base_var <- poly_term$base_var
      power <- poly_term$power
      internal_name <- poly_term$internal_name

      # Check if base variable exists in data
      if (base_var %in% names(data)) {
        # Compute polynomial: x^2, x^3, etc.
        poly_vals <- data[[base_var]]^power

        # Add to original_data if possible
        if (is.list(original_data) || is.data.frame(original_data)) {
          original_data[[internal_name]] <- poly_vals
        }

        # If hierarchical, we MUST also add it to the source dataframes in hierarchical_info
        # Otherwise d-separation tests (which re-fetch data) won't find the new variable
        if (is_hierarchical && !is.null(hierarchical_info)) {
          # Find which level/dataframe holds the base variable
          for (lvl_name in names(hierarchical_info$data)) {
            if (base_var %in% names(hierarchical_info$data[[lvl_name]])) {
              # Add computed column to this level's dataframe
              source_df <- hierarchical_info$data[[lvl_name]]
              hierarchical_info$data[[lvl_name]][[internal_name]] <- source_df[[
                base_var
              ]]^power
              break
            }
          }
        }
      }
    }
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

        # Legacy VCV/multiVCV required for some logic
        # If 'optimise=FALSE', because_model might use VCV.

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
        if (optimise) {
          data[[paste0("Prec_", s_name)]] <- P
        }

        if (is.null(N)) {
          N <- nrow(V)
        } else if (nrow(V) != N) {
          stop(paste("Dimension mismatch in", s_name))
        }

        if (length(structures) == 1 && !optimise) {
          data$VCV <- V
        }
      } else if (is.matrix(obj)) {
        V <- obj
        if (optimise) {
          data[[paste0("Prec_", s_name)]] <- solve(V)
        }

        if (is.null(N)) {
          N <- nrow(V)
        } else if (nrow(V) != N) {
          stop(paste("Dimension mismatch in", s_name))
        }
      }
    }

    if (optimise) {
      data$zeros <- rep(0, N)
    }
  }

  data$ID <- diag(N)
  if (
    !is.null(latent_method) &&
      latent_method == "correlations" &&
      !is.null(latent) &&
      length(latent) > 0
  ) {
    data$ID2 <- diag(2)
  }
  data$N <- N

  # Handle multinomial and ordinal data
  if (!is.null(distribution)) {
    # If distribution is provided but unnamed, try to auto-assign if there is only one response
    if (is.null(names(distribution))) {
      response_vars <- unique(sapply(equations, function(eq) all.vars(eq[[2]])))
      if (length(distribution) == 1 && length(response_vars) == 1) {
        names(distribution) <- response_vars
        message(sprintf(
          "Auto-assigned distribution '%s' to response variable '%s'",
          distribution,
          response_vars
        ))
      } else if (length(distribution) == length(response_vars)) {
        # Riskier, but if lengths match, assume order
        # Better to warn and ask for names.
        warning(
          "Argument 'distribution' is unnamed. Please provide a named vector like c(Response = 'binomial'). Assuming defaults (Gaussian) for safety."
        )
      } else {
        warning(
          "Argument 'distribution' is unnamed and length does not match response variables. Ignoring."
        )
      }
    }

    # Auto-fix residual variance for non-Gaussian distributions if not specified
    for (var in names(distribution)) {
      dist <- distribution[[var]]
      if (dist %in% c("binomial", "multinomial", "ordinal")) {
        should_fix <- FALSE
        if (is.null(fix_residual_variance)) {
          should_fix <- TRUE
          fix_residual_variance <- c()
        } else if (
          is.numeric(fix_residual_variance) &&
            length(fix_residual_variance) == 1 &&
            is.null(names(fix_residual_variance))
        ) {
          # It's a global fix (e.g. fix=1), so it applies to this var too. No action needed.
          should_fix <- FALSE
        } else if (!var %in% names(fix_residual_variance)) {
          should_fix <- TRUE
        }

        if (should_fix) {
          # Append to fixed variance vector
          new_fix <- setNames(1, var)
          fix_residual_variance <- c(fix_residual_variance, new_fix)

          if (!quiet) {
            message(sprintf(
              "Note: Fixing residual variance of '%s' (%s) to 1 for identifiability.",
              var,
              dist
            ))
          }
        }
      }
    }

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

  if (!is.null(distribution)) {
    for (var_name in names(distribution)) {
      dist_type <- distribution[[var_name]]
      if (dist_type %in% c("zip", "zinb")) {
        zeros_name <- paste0("zeros_", var_name)
        if (is.null(data[[zeros_name]])) {
          data[[zeros_name]] <- rep(0, N)
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
    sd_name <- paste0(var, "_sd")

    if (se_name %in% names(data)) {
      auto_variability[[var]] <- "se"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has standard errors in '%s'",
          var,
          se_name
        ))
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
  }

  # Merge auto-detected with manual specification (manual takes precedence)
  if (length(auto_variability) > 0) {
    if (is.null(variability)) {
      variability <- auto_variability
    } else {
      # Convert variability to named list if needed
      if (is.null(names(variability))) {
        # This case (unnamed vector but length > 1) is ambiguous or unsupported by global logic
        # We'll treat as "names missing" warning or error?
        # For backward compatibility / safety, just set names to values?
        # Actually, existing code: variability <- setNames(rep(NA, ...)) seems wrong if we passed values like c("reps", "se")
        # But previous logic (line 958 in original) did: variability <- setNames(rep(NA, length(variability)), variability)
        # This implied variability was a vector of NAMES.
        # But wait, the doc says variability is "c(Var = 'type')".
        # If unnamed, `variability` contains TYPES? Or NAMES?
        # Original Doc line 11 (approx): "If unnamed, it defaults to 'se' for all specified variables."
        # Meaning: variability = c("Var1", "Var2") -> Var1="se", Var2="se".
        # My new global logic handles length==1 separately.
        # If length > 1 and unnamed, we assume standard behavior (list of vars, default type "se")
        variability <- setNames(rep("se", length(variability)), variability)
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
      }

      # Ensure var is removed (latent)
      if (var_name %in% names(data)) data[[var_name]] <- NULL
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

      # Exclude polynomial internal variables (they're deterministic, not latent)
      if (!is.null(all_poly_terms)) {
        poly_internal_names <- sapply(all_poly_terms, function(x) {
          x$internal_name
        })
        potential_latents <- setdiff(potential_latents, poly_internal_names)
      }

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
        message(
          "Generating m-separation tests (MAG with latent variables)..."
        )
      } else {
        message("Generating d-separation tests...")
      }
    }

    dsep_result <- because_dsep(
      equations,
      latent = latent,
      random_terms = random_terms,
      hierarchical_info = hierarchical_info,
      poly_terms = all_poly_terms,
      categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
        attr(data, "categorical_vars")
      } else {
        NULL
      },
      quiet = !dsep
    )

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
      if (!quiet) {
        message(paste("D-sep test eq:", deparse(test_eq)))
      }

      # Select appropriate dataset for this test
      test_data <- original_data

      if (!is.null(hierarchical_info)) {
        # Extract variables from this test equation
        test_vars <- all.vars(test_eq)

        # Get appropriate dataset for these variables
        test_data <- get_data_for_variables(
          test_vars,
          hierarchical_info$data,
          hierarchical_info$levels,
          hierarchical_info$hierarchy,
          hierarchical_info$link_vars
        )

        if (!quiet) {
          message(
            "  Hierarchical data: using ",
            nrow(test_data),
            " observations for this test"
          )
        }
      }

      # Run model for this single test
      # We pass dsep=FALSE to treat it as a standard model run
      # We pass parallel=FALSE to avoid nested parallelism
      fit <- because(
        data = test_data, # Use selected dataset
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
        quiet = FALSE, # Suppress output for individual runs (temporarily ENABLED)
        dsep = FALSE,
        variability = variability,
        distribution = distribution,
        fix_residual_variance = fix_residual_variance,
        latent = NULL, # D-sep tests are on observed variables only
        latent_method = latent_method,
        parallel = FALSE, # Disable nested parallelism
        n.cores = 1,
        cl = NULL,
        ic_recompile = ic_recompile,
        random = random # Pass global random effects to d-sep tests
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
          "original_data",
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
          # We must exclude random effects grouping variables (which are not fixed predictors)
          parsed_test_eq <- extract_random_effects(list(test_eq))
          fixed_test_eq <- parsed_test_eq$fixed_equations[[1]]

          test_vars <- all.vars(fixed_test_eq)
          response <- as.character(fixed_test_eq)[2]
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
                # For Gaussian: beta_Response_Predictor format
                # Note: 'dummies' are like 'sex_m'
                params_to_monitor <- c(
                  params_to_monitor,
                  paste0("beta_", response, "_", dummies)
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

          # Remove "all" keyword if present, as it causes warnings in d-sep (rjags doesn't support it here)
          if ("all" %in% current_monitor) {
            current_monitor <- setdiff(current_monitor, "all")
          }

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
        parsed_test_eq <- extract_random_effects(list(test_eq))
        fixed_test_eq <- parsed_test_eq$fixed_equations[[1]]

        test_vars <- all.vars(fixed_test_eq)
        response <- as.character(fixed_test_eq)[2]
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
              # For Gaussian: beta_Response_Predictor format
              # Note: 'dummies' are like 'sex_m'
              params_to_monitor <- c(
                params_to_monitor,
                paste0("beta_", response, "_", dummies)
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

        # Remove "all" keyword if present (sequential mode)
        if ("all" %in% current_monitor) {
          current_monitor <- setdiff(current_monitor, "all")
        }

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
        if (coda::niter(combined_samples) != coda::niter(samples)) {
          stop("MCMC iteration mismatch between d-sep tests")
        }
        # Combine chains: for each chain, cbind the variables
        new_samples <- coda::mcmc.list()
        for (ch in 1:coda::nchain(combined_samples)) {
          # Combine matrices
          mat1 <- combined_samples[[ch]]
          mat2 <- samples[[ch]]
          # All columns from mat2 should be new (due to renaming)
          new_mat <- cbind(mat1, mat2)
          new_samples[[ch]] <- coda::mcmc(
            new_mat,
            start = stats::start(mat1),
            thin = coda::thin(mat1)
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
      dsep_results = results, # Store individual test results for summary
      induced_correlations = induced_cors
    )
    class(result) <- "because"
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
        dsep_result <- because_dsep(
          equations,
          latent = latent,
          random_terms = random_terms,
          hierarchical_info = hierarchical_info,
          quiet = !dsep
        )
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
        vars_needing_intercept <- setdiff(
          vars_with_correlations,
          response_vars
        )

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
    new_equations <- list()

    # Use for loop instead of lapply to safely modify data and collect equations
    for (idx in seq_along(equations)) {
      eq <- equations[[idx]]
      # Parse formula to get all variables
      vars <- all.vars(eq)

      # Check if any predictors are categorical
      for (var in vars) {
        if (var %in% names(categorical_vars)) {
          # Check if var is the response (LHS)
          lhs_var <- all.vars(eq[[2]])
          if (var %in% lhs_var) {
            # Skip expansion if it's the response
            next
          }

          # Get dummy variable names
          levels <- categorical_vars[[var]]$levels
          dummies <- categorical_vars[[var]]$dummies

          # If the parent variable is being imputed (is in response_vars_with_na),
          # we MUST define the dummies deterministically in JAGS to link them.
          if (var %in% response_vars_with_na) {
            for (k in 2:length(levels)) {
              dummy_name <- paste0(var, "_", levels[k])

              # Create deterministic equation: dummy ~ I(var == k)
              det_eq_str <- sprintf("%s ~ I(%s == %d)", dummy_name, var, k)
              det_eq <- stats::as.formula(det_eq_str)

              # Only add if not already present
              eq_exists <- any(sapply(c(equations, new_equations), function(e) {
                deparse(e) == deparse(det_eq)
              }))

              if (!eq_exists) {
                new_equations <- c(new_equations, list(det_eq))
                # Also remove the dummy from the data list so JAGS uses the definition
                if (dummy_name %in% names(data)) {
                  data[[dummy_name]] <- NULL
                }
              }
            }
          }

          # Convert formula to character for manipulation
          eq_str <- paste(deparse(eq), collapse = " ")

          # Replace categorical variable with its dummies
          pattern <- paste0("\\b", var, "\\b")
          replacement <- paste(dummies, collapse = " + ")
          eq_str <- gsub(pattern, replacement, eq_str)

          # Convert back to formula
          eq <- as.formula(eq_str)
          equations[[idx]] <- eq

          if (!quiet) {
            message(sprintf(
              "Expanded '%s' to: %s",
              var,
              paste(dummies, collapse = ", ")
            ))
          }
        }
      }
    }
    # Add the new deterministic equations
    equations <- c(equations, new_equations)
  }

  # JAGS model code
  model_output <- because_model(
    equations = equations,
    multi.tree = is_multiple,
    variability = variability_list,
    distribution = distribution,
    vars_with_na = response_vars_with_na,
    induced_correlations = induced_cors,
    latent = latent,
    standardize_latent = standardize_latent,
    optimise = optimise,
    structure_names = structure_names,
    random_structure_names = names(random_structures),
    random_terms = random_terms,
    poly_terms = all_poly_terms,
    categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
      attr(data, "categorical_vars")
    } else {
      NULL
    }
  )

  model_string <- model_output$model
  parameter_map <- model_output$parameter_map

  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # If latent variables are present and this is a standard run (dsep=FALSE),
  # print the MAG structure and basis set for user verification, as requested.
  # Display MAG structure for latent variable models (non-dsep runs)
  if (!dsep && !is.null(latent) && length(latent) > 0 && !quiet) {
    message("--- Latent Variable Structure (MAG) ---")
    # We call because_dsep just for its side effect (printing MAG info).
    # We wrap it in tryCatch to ensure it doesn't block the main run if it fails.
    tryCatch(
      {
        if (length(equations) > 0) {
          invisible(because_dsep(
            equations,
            latent = latent,
            random_terms = random_terms,
            quiet = FALSE
          ))
        }
      },
      error = function(e) {
        # Silently skip if MAG display fails - not critical for main run
        if (!quiet) {
          message("(MAG structure display skipped)")
        }
      }
    )
    message("---------------------------------------")
  }

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
      extract_names("^\\s*sigma"),
      extract_names("^\\s*psi"),
      extract_names("^\\s*r_"),
      extract_names("^\\s*cutpoint")
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
          gsub("^alpha_?", "", all_params) %in% response_vars &
          !gsub("^alpha_?", "", all_params) %in% mag_exogenous_vars) | # Response intercepts, excluding MAG exogenous
          grepl("^beta", all_params) | # All regression coefficients
          grepl("^rho", all_params) | # Induced correlations
          grepl("^sigma", all_params) | # Variance components
          grepl("^rho", all_params) | # Induced correlations
          grepl("^psi", all_params) | # Zero-inflation probability
          grepl("^r_", all_params) | # Negative Binomial size
          grepl("^cutpoint", all_params) | # Ordinal cutpoints
          (grepl("^lambda", all_params) &
            gsub("^lambda_?", "", all_params) %in% response_vars &
            !gsub("^lambda_?", "", all_params) %in% mag_exogenous_vars) # Response lambdas, excluding MAG exogenous
      ]
    } else {
      # monitor_mode == "all" or NULL: include everything
      monitor <- all_params

      # Also include response variables (for imputation inspection)
      # We extract them directly from the equations (which include auto-added intercept models)
      response_vars_all <- unique(sapply(equations, function(eq) {
        all.vars(formula(eq)[[2]])
      }))

      # Only add them if they are in the model (obviously)
      if (length(response_vars_all) > 0) {
        monitor <- unique(c(monitor, response_vars_all))
      }
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
    gregexpr(
      "\\b([a-zA-Z0-9_]+)\\s*\\[1:N\\]\\s*~",
      model_string,
      perl = TRUE
    )
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

      model <- rjags::jags.model(
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
      samples <- rjags::coda.samples(
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
    if (!quiet) {
      # message("--- JAGS Model Code (Pre-Compile) ---")
      # message(paste(model_output$model, collapse = "\n"))
      # message("-------------------------------------")
    }
    # Compile model

    model <- tryCatch(
      {
        rjags::jags.model(
          model_file,
          data = data,
          n.chains = n.chains,
          n.adapt = n.adapt,
          quiet = quiet
        )
      },
      error = function(e) {
        if (!quiet) {
          message("\nCRITICAL JAGS ERROR during compilation:")
          message(e$message)
          message("Check your model code syntax or data dimensions.\n")
        }
        stop(e)
      }
    )
    update(model, n.iter = n.burnin)

    # Disable DIC/WAIC if only 1 chain (rjags requirement)
    if (n.chains < 2 && (DIC || WAIC)) {
      warning(
        "DIC and WAIC require at least 2 chains. Disabling calculation."
      )
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

  # Filter internal parameters (log_lik) from summary parameters
  # We keep them in samples for WAIC calculation but hide them from the summary output
  if (!is.null(sum_stats)) {
    if (is.matrix(sum_stats$statistics)) {
      rows_to_keep <- !grepl("^log_lik", rownames(sum_stats$statistics))
      sum_stats$statistics <- sum_stats$statistics[
        rows_to_keep,
        ,
        drop = FALSE
      ]
      sum_stats$quantiles <- sum_stats$quantiles[rows_to_keep, , drop = FALSE]
    } else {
      # Single parameter case (statistics is a vector)
      # Check if the single parameter is log_lik
      param_name <- colnames(samples[[1]])
      if (length(param_name) == 1 && grepl("^log_lik", param_name)) {
        # If the only parameter is log_lik, return empty stats
        # Or handle appropriately. For now, empty seems safest or just nullify.
        sum_stats <- NULL
      }
    }
  }

  # Initialize result object
  result <- list(
    model = model,
    model_code = model_output$model,
    data = data, # Store data for recompilation if needed
    input = list(
      equations = equations,
      random = random,
      structure = structure,
      data = original_data, # Store original data too for safety
      latent = latent
    ),
    samples = samples,
    summary = sum_stats,
    monitor = monitor,
    modfile = model_file,
    dsep = dsep,
    dsep_tests = dsep_tests,
    dsep_results = dsep_results,
    parameter_map = parameter_map,
    induced_correlations = induced_cors
  )

  # Check if we can store species/unit identifiers for easy reference
  if (!is.null(tree)) {
    # If tree used, data is sorted by tip labels
    if (inherits(tree, "multiPhylo")) {
      result$species_order <- tree[[1]]$tip.label
    } else {
      result$species_order <- tree$tip.label
    }
  } else if (!is.null(id_col) && is.data.frame(original_data)) {
    # If no tree but ID col provided
    result$species_order <- as.character(original_data[[id_col]])
  }

  # Assign class immediately (needed for print/summary/waic methods)
  class(result) <- "because"

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
    # Note: WAIC calculation generally uses the posterior samples already collected.
    # We defer WAIC calculation to the common block at the end to ensure consistency.
    # if (WAIC) {
    #   result$WAIC <- because_waic(result)
    # }
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
      # Defer to common block
      # result$WAIC <- because_waic(result)
    }
  } else {
    # Sequential execution - use standard approach
    if (DIC) {
      result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
    }
    # WAIC will be computed after class assignment
  }

  # Assign class before WAIC computation (because_waic needs this)
  # Already assigned earlier
  # class(result) <- "because"

  # Compute WAIC if requested (must be after class assignment)
  if (WAIC) {
    result$WAIC <- because_waic(result)
  }

  return(result)
}

#' Preprocess categorical variables (character/factor) to integer codes and dummies
#'
#' @param data A data.frame or list of data.frames
#' @param quiet Logical; whether to suppress informational messages
#' @return The modified data object with categorical_vars attribute
#' @keywords internal
preprocess_categorical_vars <- function(data, quiet = FALSE) {
  if (is.null(data)) {
    return(NULL)
  }

  # Recursively process lists of data frames (hierarchical data)
  if (is.list(data) && !is.data.frame(data)) {
    all_cat_vars <- list()

    # Process each element (usually levels in hierarchy)
    for (i in seq_along(data)) {
      if (is.data.frame(data[[i]])) {
        processed <- preprocess_categorical_vars(data[[i]], quiet = quiet)
        data[[i]] <- processed
        # Collect categorical vars metadata
        level_cat_vars <- attr(processed, "categorical_vars")
        if (!is.null(level_cat_vars)) {
          # Use utils::modifyList if available, or manual merge
          for (name in names(level_cat_vars)) {
            all_cat_vars[[name]] <- level_cat_vars[[name]]
          }
        }
      }
    }

    # Attach merged metadata to the top-level list
    attr(data, "categorical_vars") <- all_cat_vars
    return(data)
  }

  # Process single data frame
  if (!is.data.frame(data)) {
    return(data)
  }

  char_cols <- sapply(data, function(x) is.character(x) || is.factor(x))
  if (any(char_cols)) {
    categorical_vars <- list()
    if (!is.null(attr(data, "categorical_vars"))) {
      categorical_vars <- attr(data, "categorical_vars")
    }

    col_names <- names(data)[char_cols]
    for (col in col_names) {
      # Convert to factor first to get levels
      f_vals <- factor(data[[col]])
      levels <- levels(f_vals)

      if (length(levels) < 2) {
        if (!quiet) {
          warning(sprintf(
            "Variable '%s' has < 2 levels. Converting to numeric constant.",
            col
          ))
        }
        data[[col]] <- as.numeric(f_vals)
      } else {
        # Store metadata for model expansion
        categorical_vars[[col]] <- list(
          levels = levels,
          reference = levels[1],
          dummies = paste0(col, "_", levels[-1])
        )

        # Convert to integer codes for JAGS
        data[[col]] <- as.integer(f_vals)

        # Generate Dummy Variables explicitly
        # This is required so JAGS can find 'sex_m' etc.
        for (k in 2:length(levels)) {
          lev_name <- levels[k]
          dummy_col_name <- paste0(col, "_", lev_name)
          # Create binary column
          data[[dummy_col_name]] <- as.numeric(data[[col]] == k)
        }

        if (!quiet) {
          message(sprintf(
            "Converted categorical '%s' to integers (1..%d). Reference: '%s'",
            col,
            length(levels),
            levels[1]
          ))
        }
      }
    }
    attr(data, "categorical_vars") <- categorical_vars
  }

  return(data)
}
