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

    jags_model <- model$model

    # Check if model object is valid (it might be a broken pointer from parallel execution)
    is_valid_model <- FALSE
    if (!is.null(jags_model) && inherits(jags_model, "jags")) {
        # Try to access the pointer to see if it's valid
        tryCatch(
            {
                # Accessing state() or similar method usually triggers error if pointer is nil
                # But rjags objects are external pointers.
                # Best check is to try a simple operation or check if it's NULL
                # Actually, for parallel runs, model$model might be NULL or a broken pointer
                # If it's from parallel execution without recompilation, it might be invalid.
                # Let's try to use it, and if it fails, recompile.
                is_valid_model <- TRUE
            },
            error = function(e) {
                is_valid_model <- FALSE
            }
        )
    }

    # Function to run jags.samples with error handling
    run_jags_samples <- function(mod) {
        rjags::jags.samples(
            mod,
            c("WAIC", "deviance"),
            type = "mean",
            n.iter = n.iter,
            n.burnin = n.burnin,
            n.thin = n.thin
        )
    }

    # Try to run with existing model
    samples <- NULL
    if (is_valid_model) {
        tryCatch(
            {
                samples <- run_jags_samples(jags_model)
            },
            error = function(e) {
                # If error indicates recompilation needed, set flag to recompile
                if (grepl("recompiled", e$message, ignore.case = TRUE)) {
                    is_valid_model <<- FALSE
                } else {
                    stop(e)
                }
            }
        )
    }

    # If model is missing, invalid, or failed to run, recompile it
    if (!is_valid_model || is.null(samples)) {
        if (is.null(model$model_code) || is.null(model$data)) {
            stop(
                "Model object is invalid and cannot be recompiled (missing code or data)."
            )
        }

        message("Recompiling model for WAIC calculation...")

        # Create temporary model file
        model_file <- tempfile(fileext = ".jg")
        writeLines(model$model_code, model_file)

        # Recompile with 2 chains (required for WAIC/DIC)
        jags_model <- rjags::jags.model(
            model_file,
            data = model$data, # We need to ensure data is stored in the object!
            n.chains = 2,
            n.adapt = 100, # Short adapt
            quiet = TRUE
        )

        # Short burn-in
        update(jags_model, n.iter = 100)

        # Try running again with recompiled model
        samples <- run_jags_samples(jags_model)
    }

    if (is.null(samples$WAIC)) {
        stop("JAGS did not return 'WAIC' monitor. Is 'dic' module loaded?")
    }

    samples$p_waic <- samples$WAIC
    samples$waic <- samples$deviance + samples$p_waic
    tmp <- sapply(samples, sum)
    waic <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 1)

    return(waic)
}
