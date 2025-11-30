#' Extract d-separation statements from a structural equation model
#'
#' This function takes a set of structural equations defining a causal model
#' and returns the conditional independence statements (d-separation or m-separation tests)
#' implied by the model structure. If latent variables are specified, the function
#' uses Shipley's MAG (Maximal Ancestral Graph) approach to account for unmeasured confounders.
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
#' (Shipley 2000, 2009). For standard DAGs without latent variables, it identifies
#' pairs of non-adjacent variables and creates conditional independence tests.
#'
#' When latent variables are specified, the function uses the DAG-to-MAG conversion
#' (Shipley 2016) to identify m-separation statements and induced correlations
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
#'
#' @export
#' @importFrom stats formula terms as.formula
phybase_dsep <- function(equations, latent = NULL) {
  # If no latents, use standard DAG d-separation
  if (is.null(latent)) {
    return(dsep_standard(equations))
  }

  # With latents: use MAG m-separation
  return(dsep_with_latents(equations, latent))
}

# Standard d-separation for DAGs (original logic)
dsep_standard <- function(equations) {
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

        # Combine the parents of both variables
        conditioning_vars <- unique(c(parents_var1, parents_var2))

        # Ensure that if a variable has no parents, it always goes on the right-hand side
        if (length(parents_var1) == 0 && length(parents_var2) == 0) {
          # If both have no parents, put both on the right-hand side
          reg_formula <- paste(
            var1,
            "~",
            paste(c(var2, conditioning_vars), collapse = " + ")
          )
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

  return(cond_indep_regressions)
}

# M-separation for MAGs (with latent variables)
dsep_with_latents <- function(equations, latent) {
  # Convert equations to ggm DAG format
  dag <- equations_to_dag(equations)

  # Call Shipley's DAG.to.MAG (suppressing verbose output)
  mag <- suppressMessages(DAG.to.MAG(dag, latents = latent))

  # Extract basis set from MAG
  basis <- suppressMessages(basiSet.mag(mag))

  # Convert to formula format
  tests <- mag_basis_to_formulas(basis)

  # Extract bidirected edges (induced correlations)
  correlations <- extract_bidirected_edges(mag)

  return(list(
    tests = tests,
    correlations = correlations,
    mag = mag # Include MAG for reference
  ))
}
