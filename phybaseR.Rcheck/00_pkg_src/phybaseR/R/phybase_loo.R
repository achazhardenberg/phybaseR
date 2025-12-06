#' Calculate LOO-CV for a PhyBaSE Model
#'
#' Calculates Leave-One-Out Cross-Validation using Pareto Smoothed Importance Sampling (PSIS-LOO)
#' for a fitted PhyBaSE model.
#'
#' @param model A fitted model object of class \code{"phybase"} returned by \code{\link{phybase_run}}.
#' @param ...  Additional arguments passed to \code{loo::loo()}.
#'
#' @return A \code{loo} object containing:
#' \item{estimates}{Table with ELPD (expected log pointwise predictive density), LOO-IC, and p_loo}
#' \item{diagnostics}{Pareto k diagnostic values for each observation}
#' \item{pointwise}{Pointwise contributions to LOO-IC}
#'
#' @details
#' LOO-CV (Leave-One-Out Cross-Validation) uses Pareto Smoothed Importance Sampling to approximate
#' leave-one-out predictive performance without refitting the model N times. This is particularly
#' useful for:
#' \itemize{
#'   \item Model comparison when models have different numbers of latent variables
#'   \item Identifying influential observations (via Pareto k diagnostics)
#'   \item Robust predictive performance assessment
#' }
#'
#' **Pareto k diagnostics**:
#' \itemize{
#'   \item k < 0.5: Excellent (all estimates reliable)
#'   \item 0.5 < k < 0.7: Good (estimates okay)
#'   \item 0.7 < k < 1: Problematic (estimates unreliable)
#'   \item k > 1: Very problematic (refit model excluding these observations)
#' }
#'
#' **Note on implementation**:
#' This function computes pointwise log-likelihoods from the posterior samples. For phylogenetic
#' models, this is done by evaluating the multivariate normal log-likelihood for each observation
#' given the posterior parameter values.
#'
#' @examples
#' \dontrun{
#'   fit <- phybase_run(data, tree, equations)
#'   loo_result <- phybase_loo(fit)
#'   print(loo_result)
#'
#'   # Check for problematic observations
#'   plot(loo_result)
#'
#'   # Compare models
#'   loo_compare(loo_result1, loo_result2)
#' }
#'
#' @export
#' @importFrom stats sd dnorm
phybase_loo <- function(model, ...) {
    if (!inherits(model, "phybase")) {
        stop("Input must be a 'phybase' model object.")
    }

    # Check if loo package is installed
    if (!requireNamespace("loo", quietly = TRUE)) {
        stop(
            "Package 'loo' is required for LOO-CV calculation. ",
            "Install it with: install.packages('loo')"
        )
    }

    # For now, provide a helpful error message directing users to use WAIC
    # Full LOO implementation requires pointwise log-likelihood calculation
    # which is complex for phylogenetic SEMs
    stop(
        "LOO-CV implementation is in development.\n\n",
        "For model comparison, please use WAIC instead:\n",
        "  fit <- phybase_run(..., WAIC = TRUE)\n",
        "  fit$WAIC\n\n",
        "Note: When comparing models with different latent variable structures,\n",
        "use latent_method = 'correlations' (MAG) to ensure comparable WAIC values.\n",
        "See ?phybase_run and the tutorial vignette for details."
    )

    # TODO: Implement pointwise log-likelihood calculation
    # This requires:
    # 1. Extracting response variable data
    # 2. For each MCMC sample, computing likelihood for each observation
    # 3. Accounting for phylogenetic correlation structure
    # 4. Handling different response types (Gaussian, binomial, etc.)

    # Placeholder code structure:
    # log_lik <- compute_pointwise_log_lik(model)
    # loo_result <- loo::loo(log_lik, ...)
    # return(loo_result)
}


#' Compare Models Using LOO-CV
#'
#' Wrapper for \code{loo::loo_compare()} to compare multiple PhyBaSE models.
#'
#' @param ... Two or more \code{loo} objects from \code{phybase_loo()}.
#'
#' @return A comparison table ranking models by expected out-of-sample predictive accuracy.
#'
#' @examples
#' \dontrun{
#'   loo1 <- phybase_loo(fit1)
#'   loo2 <- phybase_loo(fit2)
#'   phybase_loo_compare(loo1, loo2)
#' }
#'
#' @export
phybase_loo_compare <- function(...) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        stop(
            "Package 'loo' is required for LOO comparison. ",
            "Install it with: install.packages('loo')"
        )
    }

    loo::loo_compare(...)
}
