#' Calculate Pagel's Lambda from Variance Components
#'
#' This function calculates Pagel's lambda for each response variable in a fitted
#' \code{because} model by extracting the posterior samples of the phylogenetic
#' and residual variance components.
#'
#' @param model A fitted model object of class \code{"because"}.
#' @param prob A numeric value specifying the probability mass for the credibility intervals (default 0.95).
#'
#' @return A data frame containing the summary statistics for the derived lambda parameter(s):
#' \item{Mean}{Posterior mean}
#' \item{SD}{Posterior standard deviation}
#' \item{Median}{Posterior median}
#' \item{LowerCI}{Lower bound of the credibility interval}
#' \item{UpperCI}{Upper bound of the credibility interval}
#'
#' @details
#' In the optimized formulation of \code{because}, the phylogenetic signal (\eqn{\lambda}) is not always
#' explicitly monitored for complex models with multiple structures. However, it can be derived post-hoc
#' from the estimated standard deviations of the phylogenetic (\eqn{\sigma_{phylo}}) and residual (\eqn{\sigma_{res}}) components:
#' \deqn{\lambda = \frac{\sigma_{phylo}^2}{\sigma_{phylo}^2 + \sigma_{res}^2}}
#'
#' This function performs this calculation on the full posterior chains to provide valid Bayesian estimates and credibility intervals.
#'
#' @export
#' @examples
#' \dontrun{
#'   fit <- because(...)
#'   because_lambda(fit)
#' }
because_lambda <- function(model, prob = 0.95) {
    if (!inherits(model, "because")) {
        stop("Input must be a 'because' model object.")
    }

    if (is.null(model$samples)) {
        stop("Model object has no MCMC samples.")
    }

    # Get all parameter names
    params <- coda::varnames(model$samples)

    # Find residual sigma parameters: sigma_[RESP]_res
    res_sigmas <- grep("^sigma_.*_res$", params, value = TRUE)

    if (length(res_sigmas) == 0) {
        message(
            "No residual variance components found. Assuming model was not fit with optimized variance components."
        )
        # Check if explicit lambda exists
        lambdas <- grep("^lambda", params, value = TRUE)
        if (length(lambdas) > 0) {
            message(
                "Found explicit lambda parameters. Returning summary of those instead."
            )
            return(summary(model$samples[, lambdas, drop = FALSE])$statistics)
        }
        return(NULL)
    }

    results <- list()

    # Iterate through each residual sigma to find matching phylogenetic sigma
    for (res_param in res_sigmas) {
        # Parse response name: sigma_[RESP]_res
        # Regex: ^sigma_(.*)_res$
        response <- sub("^sigma_(.*)_res$", "\\1", res_param)

        # Look for matching phylogenetic sigma: sigma_[RESP]_phylo or similar
        # We look for sigma_[RESP]_[STRUCT] where STRUCT is likely 'phylo' or user-defined structure name
        pattern <- paste0("^sigma_", response, "_.*")
        candidates <- grep(pattern, params, value = TRUE)

        # Exclude the residual parameter itself
        candidates <- setdiff(candidates, res_param)

        # Also exclude sigma (total) if it exists
        candidates <- setdiff(candidates, paste0("sigma", response))

        if (length(candidates) == 0) {
            next # No matching structure for this response
        }

        # For each structure component found
        for (phy_param in candidates) {
            # Extract chains
            sigma_res_chain <- unlist(model$samples[, res_param])
            sigma_phy_chain <- unlist(model$samples[, phy_param])

            # Calculate lambda: Var_phy / (Var_phy + Var_res)
            # Remember parameters are SDs (sigma), so we square them
            var_res <- sigma_res_chain^2
            var_phy <- sigma_phy_chain^2

            lambda_chain <- var_phy / (var_phy + var_res)

            # Determine label: lambda_[RESP]
            # If multiple structures, maybe lambda_[RESP]_[STRUCT]
            # phy_param is like sigma_[RESP]_[STRUCT]
            struct_suffix <- sub(
                paste0("^sigma_", response, "_"),
                "",
                phy_param
            )
            label <- paste0("lambda_", response)
            if (struct_suffix != "phylo") {
                label <- paste0(label, "_", struct_suffix)
            }

            # Compute statistics
            est_mean <- mean(lambda_chain)
            est_sd <- sd(lambda_chain)
            est_median <- median(lambda_chain)

            alpha <- 1 - prob
            ci <- quantile(lambda_chain, probs = c(alpha / 2, 1 - alpha / 2))

            results[[label]] <- c(
                Mean = est_mean,
                SD = est_sd,
                Median = est_median,
                LowerCI = ci[1],
                UpperCI = ci[2]
            )
        }
    }

    if (length(results) == 0) {
        message("Could not derive any lambda parameters.")
        return(NULL)
    }

    # Convert to data frame
    df <- do.call(rbind, results)
    colnames(df) <- c(
        "Mean",
        "SD",
        "Median",
        paste0("Lower", prob * 100, "%"),
        paste0("Upper", prob * 100, "%")
    )

    return(as.data.frame(df))
}
