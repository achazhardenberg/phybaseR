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
    categorical_vars = NULL,
    family = NULL
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

        # FINAL RULE (User Override): Enforce directionality for species parameters.
        # Species parameters (p_, psi_, z_) should always be the RESPONSE in d-sep tests
        # with covariates, as covariates cannot cause parameters (biologically).
        # Regex to identify parameters: starts with p_, psi_, or z_
        is_param <- function(v) grepl("^(p_|psi_|z_)", v)

        if (is_param(var2) && !is_param(var1)) {
            # Predictor is parameter, Response is not. SWAP to make Parameter the Response.
            temp <- var1
            var1 <- var2
            var2 <- temp
        } else if (is_param(var2) && is_param(var1)) {
            # Both are parameters. Tie-breaker rule (User Request):
            # p_Species should be RESPONSE if paired with psi_Species/z_Species
            # Hierarchy: p_ > psi_/z_ > covariate

            is_p1 <- grepl("^p_", var1)
            is_p2 <- grepl("^p_", var2)

            # If Predictor is p_ and Response is NOT p_ (i.e. psi/z), SWAP.
            if (is_p2 && !is_p1) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # Build formula string

        # PREDICTOR RENAMING (User Request): Use psi_Species instead of z_Species/Species as predictor
        # This improves convergence for d-separation tests.
        if (!is.null(family)) {
            # Helper: rename if occupancy
            rename_if_occ <- function(v) {
                # Check if this variable is an occupancy variable
                # (either the base name or maybe already prefixed?)
                # Usually basis set uses base variable names (e.g. "Dingo")

                # Check if 'v' itself is in family as occupancy
                if (!is.na(family[v]) && family[v] == "occupancy") {
                    return(paste0("psi_", v))
                }

                # If it's already p_ or psi_ or z_, keep as is
                if (grepl("^(p_|psi_|z_)", v)) {
                    return(v)
                }

                return(v)
            }

            var2 <- rename_if_occ(var2)
            if (!is.null(cond_vars) && length(cond_vars) > 0) {
                cond_vars <- sapply(cond_vars, rename_if_occ)
            }
        }

        if (is.null(cond_vars) || length(cond_vars) == 0) {
            formula_str <- paste(var1, "~", var2)
        } else {
            # Sort conditioning variables for a canonical representation and stable deduplication
            sorted_cond <- sort(cond_vars)
            formula_str <- paste(
                var1,
                "~",
                var2,
                "+",
                paste(sorted_cond, collapse = " + ")
            )
        }
        f <- stats::as.formula(formula_str)
        attr(f, "test_var") <- var2 # The variable being tested for independence

        formulas[[length(formulas) + 1]] <- f
    }

    # EXCLUSION RULE (User Request): Remove tests between p_Species and psi_Species (or Species itself)
    # for the same species. These are structurally coupled and testing independence is confusing/invalid.

    if (length(formulas) > 0) {
        keep_indices <- rep(TRUE, length(formulas))

        for (i in seq_along(formulas)) {
            f <- formulas[[i]]
            # Handle potential multiple terms on RHS, get the first one (test var)
            # Actually, `f` constructed above is `var1 ~ var2 (+ conds)`
            # But `as.character(f)` yields c("~", "var1", "var2 + conds")
            # We stored `test_var` as attribute! Use that.

            test_var <- attr(f, "test_var")
            resp_var <- as.character(f)[2] # Response is reliable from formula structure

            # Helper to extract species name from p_Species or psi_Species or Species
            get_species <- function(x) {
                x <- sub("^p_", "", x)
                x <- sub("^psi_", "", x)
                x <- sub("^z_", "", x)
                return(x)
            }

            s1 <- get_species(resp_var)
            s2 <- get_species(test_var)

            # Check if one is p_ and the other is a state variable (psi, z, or raw species name)
            is_p1 <- grepl("^p_", resp_var)
            is_p2 <- grepl("^p_", test_var)

            is_state1 <- !is_p1 # Simplified: if not p, it's state (psi, z, or observed)
            is_state2 <- !is_p2

            # Condition: Same species AND one is p, one is state
            if (s1 == s2 && ((is_p1 && is_state2) || (is_state1 && is_p2))) {
                keep_indices[i] <- FALSE
            }
        }

        formulas <- formulas[keep_indices]
    }

    return(formulas)
}
