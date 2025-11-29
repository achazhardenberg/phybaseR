#' Extract d-separation statements from a structural equation model
#'
#' This function takes a set of structural equations defining a causal model
#' and returns the conditional independence statements (d-separation tests)
#' implied by the model structure. These can be used to test the model's
#' fit using the basis set method.
#'
#' @param equations A list of model formulas (one per structural equation),
#'   e.g., \code{list(Y ~ X1 + X2, Z ~ Y)}.
#'
#' @return A list of formulas representing conditional independence tests.
#'   Each formula represents a regression that should yield a non-significant
#'   coefficient if the model is correctly specified.
#'
#' @details
#' The function implements the basis set approach to d-separation testing
#' (Shipley 2000, 2009). For each pair of non-adjacent variables in the
#' causal graph, it creates a conditional independence test by regressing
#' one variable on the other while conditioning on their common parents.
#'
#' Variables with no parents are preferentially placed on the right-hand
#' side of test equations. When both variables have parents, descendants
#' are placed on the left-hand side.
#'
#' @references
#' Shipley, B. (2000). A new inferential test for path models based on
#' directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.
#'
#' Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#' context. Ecology, 90(2), 363-368.
#'
#' @examples
#' # Define a simple path model
#' equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#'
#' # Get conditional independence tests
#' ind_tests <- phybase_dsep(equations)
#'
#' # Each test is a formula
#' print(ind_tests)
#'
#' @export
phybase_dsep <- function(equations) {
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

        # Add the formula to the list
        cond_indep_regressions[[
          length(cond_indep_regressions) + 1
        ]] <- as.formula(reg_formula)
      }
    }
  }

  return(cond_indep_regressions)
}
