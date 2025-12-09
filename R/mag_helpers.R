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
            terms <- attr(stats::terms(eq), "term.labels")
            # Filter out random effect terms
            terms[!sapply(terms, is_random_term)]
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

        parents <- attr(stats::terms(eq), "term.labels")
        # Filter out random effect terms
        parents <- parents[!sapply(parents, is_random_term)]
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
mag_basis_to_formulas <- function(basis_set, latent_children = NULL) {
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
