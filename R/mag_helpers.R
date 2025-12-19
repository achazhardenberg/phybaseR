#' Convert becauseR equations to ggm DAG adjacency matrix
#'
#' @param equations List of formulas
#' @param exclude_vars Character vector of variable names to exclude (e.g., grouping variables)
#' @return Named adjacency matrix in ggm format
#' @keywords internal
equations_to_dag <- function(equations, exclude_vars = NULL) {
    # Helper function to check if a term is a random effect term
    is_random_term <- function(term) {
        # Matches patterns like "(1|var)", "(1 | var)", etc.
        grepl("\\(.*\\|.*\\)", term) || grepl("^\\d+\\s*\\|\\s*\\w+$", term)
    }

    # Parse all variables from both sides of equations
    all_vars <- unique(c(
        sapply(equations, function(eq) as.character(eq[[2]])),
        unlist(lapply(equations, function(eq) {
            # Use delete.response to get predictors safely, handling I() and interactions
            tf <- stats::terms(eq)
            rhs <- stats::delete.response(tf)

            # Extract variables
            # Note: this includes grouping variables from random terms (e.g. 'ID' in 1|ID)
            # We filter those out if needed, but since excluding random terms is generally handled
            # by 'exclude_vars' or downstream, getting all PROPER vars is safer than getting "terms".

            vars <- all.vars(rhs)
            return(vars)
        }))
    ))

    # Exclude specified variables (e.g., grouping variables)
    if (!is.null(exclude_vars)) {
        all_vars <- setdiff(all_vars, exclude_vars)
    }

    n <- length(all_vars)
    dag <- matrix(0, n, n, dimnames = list(all_vars, all_vars))

    # Fill adjacency matrix: parent -> child is coded as dag[parent, child] = 1
    for (eq in equations) {
        child <- as.character(eq[[2]])
        # Skip if child was excluded
        if (!child %in% all_vars) {
            next
        }

        # Use all.vars(delete.response(terms(eq))) to extract all predictors reliably
        # This automatically handles I(...), interactions, and functions
        rhs_formula <- stats::delete.response(stats::terms(eq))

        # We need to handle random effects.
        # terms(eq) includes random terms like '1 | group'.
        # all.vars on that returns 'group'.
        # We want to exclude purely random grouping variables from 'parents'.

        # Strategy: Get all vars from RHS. Then verify if they are ONLY in random terms.

        all_rhs_vars <- all.vars(rhs_formula)

        # Identifty vars that are part of random terms
        term_labels <- attr(stats::terms(eq), "term.labels")

        # Use vapply to ensure logical output even if term_labels is empty
        if (length(term_labels) > 0) {
            random_indices <- which(vapply(
                term_labels,
                is_random_term,
                FUN.VALUE = logical(1)
            ))
        } else {
            random_indices <- integer(0)
        }

        if (length(random_indices) > 0) {
            # Extract vars from FIXED terms only
            fixed_terms <- term_labels[-random_indices]

            if (length(fixed_terms) == 0) {
                parents <- character(0)
            } else {
                # Create dummy formula to extract vars from fixed terms
                # Paste terms together: "A + B + I(C^2)"
                dummy_rhs <- paste("~", paste(fixed_terms, collapse = " + "))
                parents <- all.vars(stats::as.formula(dummy_rhs))
            }
        } else {
            parents <- all_rhs_vars
        }

        # Filter out excluded variables from parents
        if (!is.null(exclude_vars)) {
            parents <- setdiff(parents, exclude_vars)
        }

        for (parent in parents) {
            # Only add edge if parent is in all_vars (wasn't excluded)
            if (parent %in% all_vars) {
                dag[parent, child] <- 1
            }
        }
    }

    return(dag)
}

#' Extract bidirected edges (induced correlations) from MAG
#'
#' @param mag MAG adjacency matrix from DAG.to.MAG
#' @return List of variable pairs with induced correlations
#' @keywords internal
extract_bidirected_edges <- function(mag) {
    # Bidirected edges in MAG are coded as 100
    n <- nrow(mag)
    var_names <- rownames(mag)

    correlations <- list()

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (mag[i, j] == 100 && mag[j, i] == 100) {
                # Bidirected edge found
                correlations[[length(correlations) + 1]] <- c(
                    var_names[i],
                    var_names[j]
                )
            }
        }
    }

    return(correlations)
}

#' Convert MAG basis set to becauseR formula format
#'
#' @param basis_set Basis set from basiSet.mag()
#' @param latent_children Optional character vector of variables that are direct children of latents
#' @return List of formulas with test_var attribute
#' @keywords internal
mag_basis_to_formulas <- function(
    basis_set,
    latent_children = NULL,
    categorical_vars = NULL
) {
    if (is.null(basis_set) || length(basis_set) == 0) {
        return(list())
    }

    formulas <- list()

    for (test in basis_set) {
        # test is a vector: c(var1, var2, conditioning_vars...)
        var1 <- test[1]
        var2 <- test[2]
        cond_vars <- if (length(test) > 2) test[3:length(test)] else NULL

        # Apply ordering rule: if var1 is a latent child and var2 is not,
        # swap them so the latent child becomes the predictor (test variable)
        if (!is.null(latent_children)) {
            var1_is_latent_child <- var1 %in% latent_children
            var2_is_latent_child <- var2 %in% latent_children

            # Swap if var1 is latent child but var2 is not
            if (var1_is_latent_child && !var2_is_latent_child) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # New Rule: Favor categorical variables as PREDICTORS (RHS/test_var)
        # If var1 (potential Response) is categorical, and var2 (Predictor) is NOT,
        # SWAP them so var1 becomes the predictor.
        # This allows us to use proper dummy variables in the test instead of
        # modeling the categorical variable as a Gaussian response.
        if (!is.null(categorical_vars)) {
            var1_is_cat <- var1 %in% names(categorical_vars)
            var2_is_cat <- var2 %in% names(categorical_vars)

            if (var1_is_cat && !var2_is_cat) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # Build formula string
        if (is.null(cond_vars) || length(cond_vars) == 0) {
            formula_str <- paste(var1, "~", var2)
        } else {
            formula_str <- paste(
                var1,
                "~",
                var2,
                "+",
                paste(cond_vars, collapse = " + ")
            )
        }

        f <- stats::as.formula(formula_str)
        attr(f, "test_var") <- var2 # The variable being tested for independence

        formulas[[length(formulas) + 1]] <- f
    }

    return(formulas)
}
