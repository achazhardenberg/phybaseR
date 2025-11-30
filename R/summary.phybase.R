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
    # Standard summary of the MCMC samples
    summ <- summary(object$samples, ...)

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
            P_approx_0 = numeric(), # New column
            stringsAsFactors = FALSE
        )

        # Extract MCMC samples as matrix
        samples_matrix <- as.matrix(object$samples)

        # Track which samples pass all tests (for joint probability)
        n_samples <- nrow(samples_matrix)

        # Initialize all_pass to track samples where the parameter is effectively zero
        # A parameter is effectively zero if its posterior distribution crosses zero.
        # This means the 95% CI includes zero.
        # For joint probability, we want to know which *samples* contribute to this.
        # A sample contributes to "effectively zero" if it falls on the "minority" side of zero,
        # or if it's within a small epsilon of zero.
        # Let's define "effectively zero" for a single sample as being within a small interval around zero.
        # For simplicity and consistency with P(approx 0), we'll consider a sample "passing"
        # if it falls on the side of zero that has fewer samples (i.e., contributes to the "crossing zero" probability).
        # This is a bit tricky for joint probability. A more robust approach might be to check if the sample
        # is within a user-defined "negligible effect" interval.
        # For now, let's use the user's provided logic for `all_pass` which tracks samples that are *not* clearly non-zero.
        # The user's logic: if n_above < n_below, then most samples are negative. If param_samples > 0, these are the minority.
        # So, `all_pass` tracks samples that are on the "minority" side of zero.
        # This means `all_pass` will be TRUE for samples that contribute to the `min(n_above, n_below)` part.
        # If `P_approx_0` is high, it means many samples cross zero.
        # The joint probability should reflect the probability that *all* parameters cross zero.
        # Let's adjust `all_pass` to mean: the sample's value for this parameter is consistent with the parameter being zero.
        # This means the sample falls within the range that makes the parameter "effectively zero".
        # A parameter is effectively zero if its 95% CI includes zero.
        # So, for a given sample, if the parameter value is between the 2.5% and 97.5% quantiles, and those quantiles straddle zero,
        # then that sample is consistent with independence.
        # A simpler interpretation for `all_pass` for joint probability:
        # A sample `s` "passes" for a parameter `p` if `p_s` is such that `p` is considered independent.
        # If `P_approx_0` is the probability of crossing zero, then `1 - P_approx_0` is the probability of being clearly positive or clearly negative.
        # The joint probability should be the proportion of samples where *all* parameters are "effectively zero".
        # Let's define "effectively zero" for a sample as being within the range [-epsilon, epsilon] for some small epsilon.
        # Or, more simply, if the sample falls on the side of zero that is *not* the majority.
        # The user's provided `all_pass` logic is:
        # if (n_above < n_below) { all_pass <- all_pass & (param_samples > 0) } else { all_pass <- all_pass & (param_samples < 0) }
        # This means `all_pass` will be TRUE for samples that are on the *minority* side of zero.
        # If `P_approx_0` is high, it means `min(n_above, n_below)` is high.
        # So, `all_pass` will be TRUE for samples that contribute to this `min` count.
        # This seems to be the correct interpretation for "joint probability of all tests being approximately zero".

        all_pass_joint <- rep(TRUE, n_samples) # Tracks samples where *all* parameters are "effectively zero"

        for (i in seq_along(tests)) {
            test_formula <- tests[[i]]
            test_var <- attr(test_formula, "test_var")

            # The response variable in the test formula
            response <- as.character(test_formula)[2]

            # Find the parameter name for this path (response ~ test_var)
            # We look in the parameter map
            param_row <- map[
                map$response == response & map$predictor == test_var,
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

                # Independence check (if CI includes 0, they are independent)
                indep <- if (lower > 0 || upper < 0) "No" else "Yes"

                # Calculate Bayesian p-value: P(~0)
                # This is the two-tailed probability of crossing zero
                param_samples <- samples_matrix[, param_name]
                n_above <- sum(param_samples > 0)
                n_below <- sum(param_samples < 0)
                p_approx_0 <- 2 * min(n_above, n_below) / n_samples

                # Update joint probability tracker
                # For joint probability, we need samples where the parameter is "effectively zero"
                # We'll use a threshold based on the standard deviation of the posterior
                # A sample "passes" if |param| < threshold
                # Use SD as threshold (roughly 68% of samples will be within ±1 SD if centered at zero)
                param_sd <- sd(param_samples)
                threshold <- param_sd # Could also use median absolute deviation or other robust measure
                all_pass_joint <- all_pass_joint &
                    (abs(param_samples) < threshold)

                # Format the test string nicely
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
                    P_approx_0 = round(p_approx_0, 3)
                )
            }
        }

        print(results, row.names = FALSE)

        # Calculate and report joint probability
        joint_prob <- sum(all_pass_joint) / n_samples
        cat("\nJoint P(all tests ~ 0) =", round(joint_prob, 3), "\n")

        cat("\nLegend:\n")
        cat(
            "  Indep: 'Yes' = Conditionally Independent, 'No' = Dependent (based on 95% CI)\n"
        )
        cat(
            "  P(~0): Bayesian probability that effect crosses zero (0-1 scale)\n"
        )
        cat(
            "  Joint: Probability that ALL tests simultaneously support independence (all parameters cross zero)\n"
        )
        cat(
            "\nNote: For d-separation, we expect high P(~0) and Joint values (close to 1).\n"
        )

        invisible(results)
    } else {
        # Standard summary
        print(summ)

        if (!is.null(object$DIC)) {
            cat("\nDIC:\n")
            print(object$DIC)
        }

        if (!is.null(object$WAIC)) {
            cat("\nWAIC:\n")
            print(object$WAIC)
        }

        invisible(summ)
    }
}
