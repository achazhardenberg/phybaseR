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
#' ind_tests <- phybase_dsep(equations)
#'
#' # With latent variable
#' equations_latent <- list(X ~ Quality, Y ~ Quality)
#' result <- phybase_dsep(equations_latent, latent = "Quality")
#' # result$tests: m-separation tests
#' # result$correlations: induced correlation between X and Y
#' @param quiet Logical; if FALSE (default), print the basis set and MAG structure.
#'   If TRUE, suppress informational output.
#' @export
#' @importFrom stats formula terms as.formula
phybase_dsep <- function(equations, latent = NULL, quiet = FALSE) {
  # If no latents, use standard DAG d-separation
  if (is.null(latent)) {
    return(dsep_standard(equations, quiet = quiet))
  }

  # With latents: use MAG m-separation
  return(dsep_with_latents(equations, latent, quiet = quiet))
}

# Standard d-separation for DAGs (original logic)
dsep_standard <- function(equations, quiet = FALSE) {
  # Parse equations to extract parent-child relationships
  parents <- list()
  children <- list()

  for (eq in equations) {
    lhs <- as.character(formula(eq))[2] # Left-hand side (child)
    rhs <- attr(terms(eq), "term.labels") # Right-hand side (parents)

    parents[[lhs]] <- rhs
    for (var in rhs) {
      children[[var]] <- c(children[[var]], lhs)
    }
  }

  # Get the list of all variables (children and parents)
  all_vars <- unique(c(names(parents), unlist(parents)))

  # Initialize a list to store the conditional independence regressions
  cond_indep_regressions <- list()

  # Loop through all pairs of non-adjacent variables
  for (i in 1:(length(all_vars) - 1)) {
    for (j in (i + 1):length(all_vars)) {
      var1 <- all_vars[i]
      var2 <- all_vars[j]

      # Check if var1 and var2 are non-adjacent (no direct edge between them)
      if (!var2 %in% children[[var1]] && !var1 %in% children[[var2]]) {
        # Get the parents of both variables
        parents_var1 <- parents[[var1]]
        parents_var2 <- parents[[var2]]

        # Determine the conditioning set
        if (!is.null(parents_var1) && !is.null(parents_var2)) {
          conditioning_vars <- intersect(parents_var1, parents_var2)
        } else {
          conditioning_vars <- NULL
        }

        # Create the regression formula
        if (length(parents_var1) == 0 && length(parents_var2) == 0) {
          # If both variables have no parents, put var2 on the right-hand side
          reg_formula <- paste(var1, "~", var2)
        } else if (length(parents_var1) == 0) {
          # If var1 has no parents, put var1 on the right-hand side
          reg_formula <- paste(
            var2,
            "~",
            paste(c(var1, conditioning_vars), collapse = " + ")
          )
        } else if (length(parents_var2) == 0) {
          # If var2 has no parents, put var2 on the right-hand side
          reg_formula <- paste(
            var1,
            "~",
            paste(c(var2, conditioning_vars), collapse = " + ")
          )
        } else {
          # Apply the previous logic where the child/grandchild goes on the left-hand side
          if (var1 %in% children[[var2]]) {
            reg_formula <- paste(
              var2,
              "~",
              paste(c(var1, conditioning_vars), collapse = " + ")
            )
          } else {
            reg_formula <- paste(
              var1,
              "~",
              paste(c(var2, conditioning_vars), collapse = " + ")
            )
          }
        }

        # Create the formula
        f <- as.formula(reg_formula)

        # Identify the test variable
        if (length(parents_var1) == 0 && length(parents_var2) == 0) {
          test_var <- var2
        } else if (length(parents_var1) == 0) {
          test_var <- var1
        } else if (length(parents_var2) == 0) {
          test_var <- var2
        } else {
          if (var1 %in% children[[var2]]) {
            test_var <- var1
          } else {
            test_var <- var2
          }
        }

        attr(f, "test_var") <- test_var

        # Add the formula to the list
        cond_indep_regressions[[
          length(cond_indep_regressions) + 1
        ]] <- f
      }
    }
  }

  # Print basis set if not quiet
  if (!quiet) {
    cat("Basis Set for DAG:", "\n")
    cat(
      "I(X,Y|Z) means X is d-separated from Y given the set Z in the DAG",
      "\n"
    )
    if (length(cond_indep_regressions) == 0) {
      cat("No elements in the basis set", "\n")
    } else {
      for (test in cond_indep_regressions) {
        cat(format_dsep_test(test), "\n")
      }
    }
  }

  return(cond_indep_regressions)
}

# M-separation for MAGs (with latent variables)
dsep_with_latents <- function(equations, latent, quiet = FALSE) {
  # Convert equations to ggm DAG format
  dag <- equations_to_dag(equations)

  # Call local DAG.to.MAG (suppress cat output when quiet=TRUE)
  if (quiet) {
    # Capture printed output and assign the result to 'mag'
    invisible(capture.output(
      {
        mag <- suppressMessages(DAG.to.MAG(dag, latents = latent))
      },
      type = "output"
    ))
  } else {
    mag <- suppressMessages(DAG.to.MAG(dag, latents = latent))
  }

  # Extract basis set from MAG (suppress cat output when quiet=TRUE)
  if (quiet) {
    invisible(capture.output(
      {
        basis <- suppressMessages(basiSet.mag(mag))
      },
      type = "output"
    ))
  } else {
    basis <- suppressMessages(basiSet.mag(mag))
  }

  # Convert to formula format
  tests <- mag_basis_to_formulas(basis)

  # Extract bidirected edges (induced correlations)
  correlations <- extract_bidirected_edges(mag)

  # Print basis set if not quiet (our own formatted version)
  if (!quiet) {
    cat("Basis Set for MAG:", "\n")
    cat(
      "I(X,Y|Z) means X is m-separated from Y given the set Z in the MAG",
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
