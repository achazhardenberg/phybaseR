#' Extract d-separation statements from a structural equation model
#'
#' This function takes a set of structural equations defining a causal model
#' and returns the conditional independence statements (d-separation or m-separation tests)
#' implied by the model structure. If latent variables are specified, the function
#' uses the MAG (Maximal Ancestral Graph) approach by Shipley and Douma (2021)
#' to account for unmeasured latent variables.
#'
#' @param equations A list of model formulas (one per structural equation),
#'   e.g., \code{list(Y ~ X1 + X2, Z ~ Y)}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If provided, the function converts the DAG to a MAG and returns m-separation tests.
#'
#' @return If \code{latent} is NULL, returns a list of formulas representing
#'   conditional independence tests. If \code{latent} is specified, returns a list with:
#'   \itemize{
#'     \item \code{tests}: List of m-separation test formulas
#'     \item \code{correlations}: List of variable pairs with induced correlations
#'   }
#'
#' @details
#' The function implements the basis set approach to d-separation testing
#' (Shipley 2000, 2009, 2016). For standard DAGs without latent variables, it identifies
#' pairs of non-adjacent variables and creates conditional independence tests.
#'
#' When latent variables are specified, the function uses the DAG-to-MAG conversion
#' (Shipley & Douma 2021) to identify m-separation statements and induced correlations
#' among observed variables that arise from shared latent common causes.
#'
#' @references
#' Shipley, B. (2000). A new inferential test for path models based on
#' directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.
#'
#' Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#' context. Ecology, 90(2), 363-368.
#'
#' Shipley, B. (2016). Cause and Correlation in Biology (2nd ed.).
#' Cambridge University Press.
#'
#' Shipley, B., & Douma, J. C. (2021). Testing Piecewise Structural Equations
#' Models in the Presence of Latent Variables and Including Correlated Errors.
#' Structural Equation Modeling: A Multidisciplinary Journal, 28(4), 582–589.
#' https://doi.org/10.1080/10705511.2020.1871355

#'
#' @examples
#' # Standard DAG
#' equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#' ind_tests <- because_dsep(equations)
#'
#' # With latent variable
#' equations_latent <- list(X ~ Quality, Y ~ Quality)
#' result <- because_dsep(equations_latent, latent = "Quality")
#' # result$tests: m-separation tests
#' # result$correlations: induced correlation between X and Y
#' @param quiet Logical; if FALSE (default), print the basis set and MAG structure.
#'   If TRUE, suppress informational output.
#' @param random_terms Optional list of random effects (group, type) parsed from equations.
#' @param hierarchical_info Internal argument used to pass data hierarchy information
#'   (levels, grouping variables) for future implementation of multilevel d-separation
#'   tests (following Shipley 2009). Currently unused by the d-separation logic.
#' @export
#' @importFrom stats formula terms as.formula
because_dsep <- function(
  equations,
  latent = NULL,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  quiet = FALSE
) {
  # If no latents, use standard DAG d-separation
  if (is.null(latent)) {
    return(dsep_standard(
      equations,
      random_terms = random_terms,
      hierarchical_info = hierarchical_info,
      poly_terms = poly_terms,
      quiet = quiet
    ))
  }

  # With latents: use MAG m-separation
  return(dsep_with_latents(
    equations,
    latent,
    random_terms = random_terms,
    hierarchical_info = hierarchical_info,
    poly_terms = poly_terms,
    quiet = quiet
  ))
}

# Standard d-separation for DAGs (using ggm::basiSet)
dsep_standard <- function(
  equations,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  quiet = FALSE
) {
  # Extract grouping variables from random terms to exclude from DAG
  grouping_vars <- NULL
  if (length(random_terms) > 0) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))
  }

  # Extract polynomial internal variables to exclude
  poly_internal_vars <- NULL
  # Use passed poly_terms if available, otherwise try to detect (fallback)
  if (!is.null(poly_terms)) {
    all_poly_terms <- poly_terms
  } else {
    all_poly_terms <- get_all_polynomial_terms(equations)
  }

  if (!is.null(all_poly_terms)) {
    poly_internal_vars <- sapply(all_poly_terms, function(x) x$internal_name)
  }

  # Combine exclusions
  exclude_vars <- c(grouping_vars, poly_internal_vars)

  # Convert equations to adjacency matrix (DAG)
  dag <- equations_to_dag(equations, exclude_vars = exclude_vars)

  # Use ggm::basiSet to get the correct d-separation basis set
  if (!requireNamespace("ggm", quietly = TRUE)) {
    stop("Package 'ggm' is required for d-separation tests.")
  }

  basis <- ggm::basiSet(dag)

  # Check for polynomial term injection in conditioning sets
  if (!is.null(all_poly_terms) && length(basis) > 0) {
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        new_cond_vars <- cond_vars

        # Check each conditioning variable
        for (cv in cond_vars) {
          # Find any polynomial terms derived from this variable
          for (pt in all_poly_terms) {
            if (pt$base_var == cv) {
              new_cond_vars <- c(new_cond_vars, pt$internal_name)
            }
          }
        }
        new_cond_vars <- unique(new_cond_vars)

        # Reconstruct test vector
        return(c(test[1:2], new_cond_vars))
      } else {
        return(test)
      }
    })
  }

  # Convert basis set to formula list
  # We reuse mag_basis_to_formulas as the format is identical (list of vectors)
  tests <- mag_basis_to_formulas(basis)

  # Append random terms if relevant (same logic as in with_latents)
  if (length(random_terms) > 0 && length(tests) > 0) {
    new_tests <- list()
    for (t_idx in seq_along(tests)) {
      t_eq <- tests[[t_idx]]
      resp <- as.character(t_eq)[2]

      # Find random terms for response
      vocab_rand <- Filter(function(x) x$response == resp, random_terms)

      if (length(vocab_rand) > 0) {
        rand_str <- paste(
          sapply(vocab_rand, function(rt) {
            paste0("(1 | ", rt$group, ")")
          }),
          collapse = " + "
        )

        f_str <- paste(deparse(t_eq), collapse = " ")
        f_str <- paste0(f_str, " + ", rand_str)
        new_eq <- as.formula(f_str)
        # Preserve attribute
        # Re-attach test_var from original
        attr(new_eq, "test_var") <- attr(t_eq, "test_var")
        new_tests[[t_idx]] <- new_eq
      } else {
        new_tests[[t_idx]] <- t_eq
      }
    }
    tests <- new_tests
  }

  # Print basis set if not quiet
  if (!quiet) {
    cat("Basis Set for DAG:", "\n")
    cat(
      "I(X,Y|Z) means X is d-separated from Y given the set Z in the DAG",
      "\n"
    )
    if (length(tests) == 0) {
      cat("No elements in the basis set", "\n")
    } else {
      for (test in tests) {
        cat(format_dsep_test(test), "\n")
      }
    }
  }

  return(tests)
}

# M-separation for MAGs (with latent variables)
dsep_with_latents <- function(
  equations,
  latent,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  quiet = FALSE
) {
  # Extract grouping variables from random terms to exclude from DAG
  grouping_vars <- NULL
  if (length(random_terms) > 0) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))
  }

  # Extract polynomial internal variables to exclude from DAG
  # They're deterministic transformations, not causal nodes
  poly_internal_vars <- NULL
  # Use passed poly_terms if available, otherwise try to detect (fallback)
  if (!is.null(poly_terms)) {
    all_poly_terms <- poly_terms
  } else {
    all_poly_terms <- get_all_polynomial_terms(equations)
  }

  if (!is.null(all_poly_terms)) {
    poly_internal_vars <- sapply(all_poly_terms, function(x) x$internal_name)
  }

  # Combine exclusions
  exclude_vars <- c(grouping_vars, poly_internal_vars)

  # Convert equations to ggm DAG format (excluding grouping & polynomial variables)
  dag <- equations_to_dag(equations, exclude_vars = exclude_vars)

  # Always suppress DAG.to.MAG output - we'll print our own filtered version
  invisible(capture.output(
    {
      mag <- suppressMessages(DAG.to.MAG(dag, latents = latent))
    },
    type = "output"
  ))

  # Extract basis set from MAG
  # Always capture output because basiSet.mag prints to stdout,
  # and we want to print our own modified version (with random effects) later.
  invisible(capture.output(
    {
      basis <- suppressMessages(basiSet.mag(mag))
    },
    type = "output"
  ))

  # Check for polynomial term injection in conditioning sets
  if (!is.null(all_poly_terms) && length(basis) > 0) {
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        new_cond_vars <- cond_vars

        # Check each conditioning variable
        for (cv in cond_vars) {
          # Find any polynomial terms derived from this variable
          for (pt in all_poly_terms) {
            if (pt$base_var == cv) {
              new_cond_vars <- c(new_cond_vars, pt$internal_name)
            }
          }
        }
        new_cond_vars <- unique(new_cond_vars)

        # Reconstruct test vector
        return(c(test[1:2], new_cond_vars))
      } else {
        return(test)
      }
    })
  }

  # Filter out random effect grouping variables from basis set conditioning sets
  # These should not be treated as causal/fixed predictors
  if (length(random_terms) > 0 && !is.null(basis)) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))

    # basis is a list where each element is c(var1, var2, cond_var1, cond_var2, ...)
    # Remove grouping variables from conditioning sets (positions 3+)
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        # Keep var1 and var2, filter conditioning variables
        cond_vars <- test[3:length(test)]
        filtered_cond <- cond_vars[!cond_vars %in% grouping_vars]
        c(test[1:2], filtered_cond)
      } else {
        test
      }
    })
  }

  # Identify variables that are direct children of latent variables
  # These should be predictors (not responses) in independence tests
  latent_children <- character(0)
  if (!is.null(latent) && length(latent) > 0) {
    all_vars <- rownames(dag)
    for (lat in latent) {
      if (lat %in% all_vars) {
        # Find children of this latent (dag[latent, child] == 1)
        children <- all_vars[dag[lat, ] == 1]
        latent_children <- unique(c(latent_children, children))
      }
    }
  }

  # Convert to formula format, with variable ordering preference
  tests <- mag_basis_to_formulas(basis, latent_children = latent_children)

  # Save tests without random effects for clean display
  tests_for_display <- tests

  # Append random terms to MAG tests
  if (length(random_terms) > 0 && length(tests) > 0) {
    new_tests <- list()
    for (t_idx in seq_along(tests)) {
      t_eq <- tests[[t_idx]]
      resp <- as.character(t_eq)[2]

      # Find random terms for response
      vocab_rand <- Filter(function(x) x$response == resp, random_terms)

      if (length(vocab_rand) > 0) {
        rand_str <- paste(
          sapply(vocab_rand, function(rt) {
            paste0("(1 | ", rt$group, ")")
          }),
          collapse = " + "
        )

        # Rebuild formula
        # deparse might wrap lines, be careful
        # Use reliable string construction
        rhs <- labels(terms(t_eq))
        # d-sep output usually has predictors on RHS
        # Reconstruct: Resp ~ Preds + Random
        # Paste to the formula string representation

        # Safe approach: deparse
        f_str <- paste(deparse(t_eq), collapse = " ")
        f_str <- paste0(f_str, " + ", rand_str)

        new_eq <- as.formula(f_str)

        new_tests[[t_idx]] <- new_eq
      } else {
        new_tests[[t_idx]] <- t_eq
      }
    }
    tests <- new_tests
  }

  # Extract bidirected edges (induced correlations)
  correlations <- extract_bidirected_edges(mag)

  # Print basis set if not quiet (our own formatted version)
  # Use tests_for_display (without random effects) to avoid showing grouping vars
  if (!quiet) {
    cat("Basis Set for MAG:", "\n")
    cat(
      "I(X,Y|Z) means X is m-separated from Y given the set Z in the MAG",
      "\n"
    )
    if (length(tests_for_display) == 0) {
      cat("No elements in the basis set", "\n")
    } else {
      for (test in tests_for_display) {
        cat(format_dsep_test(test), "\n")
      }
    }
  }

  return(list(
    tests = tests,
    correlations = correlations,
    mag = mag # Include MAG for reference
  ))
}

# Helper to format a d-sep test for printing
format_dsep_test <- function(test) {
  # Extract variables from formula
  vars <- all.vars(test)
  response <- as.character(test)[2]
  predictors <- setdiff(vars, response)

  if (length(predictors) == 1) {
    # No conditioning set
    return(paste0("I( ", response, " , ", predictors[1], " | ", " )"))
  } else {
    # First predictor is the test variable, rest are conditioning
    test_var <- predictors[1]
    cond_set <- paste(predictors[-1], collapse = ", ")
    return(paste0("I( ", response, " , ", test_var, " | ", cond_set, " )"))
  }
}
