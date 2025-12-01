#' Calculate WAIC for a PhyBaSE Model
#'
#' Calculates the Widely Applicable Information Criterion (WAIC) for a fitted PhyBaSE model.
#'
#' @param model A fitted model object of class \code{"phybase"} returned by \code{\link{phybase_run}}.
#' @param n.iter Number of iterations for WAIC calculation. Default is 2000.
#' @param n.burnin Number of burn-in iterations. Default is 500.
#' @param n.thin Thinning interval. Default is 10.
#'
#' @return A named vector containing:
#' \item{waic}{The calculated WAIC value.}
#' \item{p_waic}{The effective number of parameters (penalty term).}
#'
#' @details
#' This function uses the \code{dic} module in JAGS to monitor \code{WAIC} and \code{deviance}.
#' The WAIC is calculated as: \code{WAIC = mean(deviance) + p_waic}.
#'
#' @examples
#' \dontrun{
#'   fit <- phybase_run(data, tree, equations)
#'   waic_res <- phybase_waic(fit, n.iter = 1000)
#'   print(waic_res)
#' }
#'
#' @export
#' @importFrom rjags jags.samples load.module
phybase_waic <- function(model, n.iter = 2000, n.burnin = 500, n.thin = 10) {
    if (!inherits(model, "phybase")) {
        stop("Input must be a 'phybase' model object.")
    }

    # Ensure dic module is loaded
    rjags::load.module("dic")

    samples <- rjags::jags.samples(
        model$model,
        c("WAIC", "deviance"),
        type = "mean",
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin
    )

    if (is.null(samples$WAIC)) {
        stop("JAGS did not return 'WAIC' monitor. Is 'dic' module loaded?")
    }

    samples$p_waic <- samples$WAIC
    samples$waic <- samples$deviance + samples$p_waic
    tmp <- sapply(samples, sum)
    waic <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 1)

    return(waic)
}
