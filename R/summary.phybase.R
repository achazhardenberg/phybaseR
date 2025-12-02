#' Summary for PhyBaSE Model
#'
#' Summarizes the output of a PhyBaSE model run.
#'
#' @param object A fitted model object of class \code{"phybase"}.
#' @param ... Additional arguments passed to \code{\link[coda]{summary.mcmc}}.
#'
#' @return A summary object containing statistics for the monitored parameters.
#'   If \code{dsep = TRUE} was used in \code{phybase_run}, the summary focuses on
#'   the conditional independence tests.
#'
#' @export
summary.phybase <- function(object, ...) {
    # Use stored summary if available, otherwise calculate it
    if (!is.null(object$summary)) {
        summ <- object$summary
    } else {
        summ <- summary(object$samples, ...)
    }

    # Calculate convergence diagnostics
    n_chains <- coda::nchain(object$samples)

    # Effective sample size
    eff_size <- tryCatch(
        coda::effectiveSize(object$samples),
        error = function(e) return(NULL)
    )

    # Gelman-Rubin diagnostic (Rhat)
    # Use stored Rhat if available (calculated in phybase_run)
    rhat <- NULL
    if ("Rhat" %in% colnames(summ$statistics)) {
        rhat <- summ$statistics[, "Rhat"]
        names(rhat) <- rownames(summ$statistics)
    } else if (n_chains > 1) {
        # Try to calculate if not stored
        rhat <- tryCatch(
            coda::gelman.diag(object$samples, multivariate = FALSE)$psrf[, 1],
            error = function(e) return(NULL)
        )
    }

    # If this was a d-sep run, we want to format the output specifically
    if (!is.null(object$dsep) && object$dsep) {
        cat("PhyBaSE d-separation Tests\n")
        cat("==========================\n\n")

        tests <- object$dsep_tests
        map <- object$parameter_map

        # Create a results table
        results <- data.frame(
            Test = character(),
            Parameter = character(),
            Estimate = numeric(),
            LowerCI = numeric(),
            UpperCI = numeric(),
            Indep = character(),
            P_approx_0 = numeric(),
            Rhat = numeric(),
            n.eff = numeric(),
            stringsAsFactors = FALSE
        )

        # Extract MCMC samples as matrix
        samples_matrix <- as.matrix(object$samples)

        # Track which samples pass all tests (for joint probability)
        n_samples <- nrow(samples_matrix)

        for (i in seq_along(tests)) {
            test_formula <- tests[[i]]
            test_var <- attr(test_formula, "test_var")

            # The response variable in the test formula
            response <- as.character(test_formula)[2]

            # Find the parameter name for this path (response ~ test_var)
            # We filter by equation_index to handle cases where the same response/predictor pair
            # appears in multiple equations (e.g. different d-sep tests)
            param_row <- map[
                map$response == response &
                    map$predictor == test_var &
                    map$equation_index == i,
            ]

            if (nrow(param_row) == 0) {
                warning(paste(
                    "Could not find parameter for test:",
                    format(test_formula)
                ))
                next
            }

            param_name <- param_row$parameter

            # Extract from summary
            if (param_name %in% rownames(summ$quantiles)) {
                est <- summ$statistics[param_name, "Mean"]
                lower <- summ$quantiles[param_name, "2.5%"]
                upper <- summ$quantiles[param_name, "97.5%"]

                # Get diagnostics for this parameter
                p_rhat <- if (!is.null(rhat) && param_name %in% names(rhat)) {
                    rhat[param_name]
                } else {
                    NA
                }
                p_neff <- if (
                    !is.null(eff_size) && param_name %in% names(eff_size)
                ) {
                    eff_size[param_name]
                } else {
                    NA
                }

                # Independence check (if CI includes 0, they are independent)
                indep <- if (lower > 0 || upper < 0) "No" else "Yes"

                # Calculate Bayesian p-value: P(~0)
                # This is the two-tailed probability of crossing zero
                param_samples <- samples_matrix[, param_name]
                n_above <- sum(param_samples > 0)
                n_below <- sum(param_samples < 0)
                p_approx_0 <- 2 * min(n_above, n_below) / n_samples

                rhs <- attr(stats::terms(test_formula), "term.labels")
                cond_vars <- setdiff(rhs, test_var)

                test_str <- paste0(
                    response,
                    " _||_ ",
                    test_var,
                    if (length(cond_vars) > 0) {
                        paste0(" | {", paste(cond_vars, collapse = ","), "}")
                    } else {
                        ""
                    }
                )

                results[nrow(results) + 1, ] <- list(
                    Test = test_str,
                    Parameter = param_name,
                    Estimate = round(est, 3),
                    LowerCI = round(lower, 3),
                    UpperCI = round(upper, 3),
                    Indep = indep,
                    P_approx_0 = round(p_approx_0, 3),
                    Rhat = round(p_rhat, 3),
                    n.eff = round(p_neff, 0)
                )
            }
        }

        print(results, row.names = FALSE)

        cat("\nLegend:\n")
        cat(
            "  Indep: 'Yes' = Conditionally Independent, 'No' = Dependent (based on 95% CI)\n"
        )
        cat(
            "  P(~0): Bayesian probability that effect crosses zero (0-1 scale)\n"
        )
        cat(
            "\nNote: For d-separation, we expect high P(~0) values (close to 1).\n"
        )

        invisible(results)
    } else {
        # Standard summary
        # Combine statistics with diagnostics
        stats_table <- summ$statistics[,
            c("Mean", "SD", "Naive SE", "Time-series SE"),
            drop = FALSE
        ]
        quant_table <- summ$quantiles[, c("2.5%", "50%", "97.5%"), drop = FALSE]

        # Create combined table
        combined <- cbind(stats_table, quant_table)

        # Add Rhat if available
        if (!is.null(rhat)) {
            # Ensure alignment
            rhat_aligned <- rhat[rownames(combined)]
            combined <- cbind(combined, Rhat = rhat_aligned)
        }

        # Add n.eff if available
        if (!is.null(eff_size)) {
            eff_aligned <- eff_size[rownames(combined)]
            combined <- cbind(combined, n.eff = eff_aligned)
        }

        print(combined)

        if (!is.null(object$DIC)) {
            cat("\nDIC:\n")
            print(object$DIC)
        }

        if (!is.null(object$WAIC)) {
            cat("\nWAIC:\n")
            print(object$WAIC)
        }

        invisible(combined)
    }
}
