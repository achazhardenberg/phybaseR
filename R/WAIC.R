#' Compute WAIC for a JAGS model
#'
#' This function computes the WAIC (Widely Applicable Information Criterion) from a fitted JAGS model
#' by extracting the deviance and WAIC-related components from posterior samples. The function recompiles
#' the model internally to allow sampling.
#'
#' @param model An object returned by \code{\link[R2jags]{jags}} (from \pkg{R2jags} or \code{phybase_run()}).
#' @param n.iter Number of MCMC iterations for WAIC computation (default: 2000).
#' @param n.burnin Number of burn-in iterations to discard (default: 500).
#' @param n.thin Thinning interval (default: 10).
#'
#' @return A named numeric vector with WAIC and effective number of parameters (\code{p_waic}).
#' @export
#'
#' @examples
#' # Assuming you have a model fitted with phybase_run()
#' # WAIC(CB2.pb)
WAIC <- function(model, n.iter = 2000, n.burnin = 500, n.thin = 10) {
  if (!"model" %in% names(model)) {
    stop("The input must be the output of R2jags::jags or phybase_run().")
  }

  # Recompile the JAGS model
  recompiled_model <- rjags::jags.model(
    file = model$model,
    data = model$model$data,
    n.chains = model$nchain,
    n.adapt = 0,
    quiet = TRUE
  )

  # Burn-in
  update(recompiled_model, n.iter = n.burnin)

  # Sample WAIC and deviance
  samples <- rjags::jags.samples(recompiled_model,
                                 variable.names = c("WAIC", "deviance"),
                                 n.iter = n.iter,
                                 thin = n.thin,
                                 type = "mean")

  if (!all(c("WAIC", "deviance") %in% names(samples))) {
    stop("WAIC and/or deviance not found in model output. Ensure they are computed inside the JAGS model.")
  }

  p_waic <- sum(samples$WAIC)
  deviance <- sum(samples$deviance)
  waic <- deviance + p_waic

  return(round(c(waic = waic, p_waic = p_waic), 1))
}
