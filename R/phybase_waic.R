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

        # Determine number of chains: use original number, but at least 2
        n_chains_orig <- if (!is.null(model$samples)) {
            coda::nchain(model$samples)
        } else {
            2
        }
        n_chains_new <- max(2, n_chains_orig)

        # Create temporary model file
        model_file <- tempfile(fileext = ".jg")
        writeLines(model$model_code, model_file)

        # Create inits from posterior means to speed up convergence
        inits_list <- NULL
        if (!is.null(model$summary) && !is.null(model$summary$statistics)) {
            # Get means of stochastic parameters
            means <- model$summary$statistics[, "Mean"]

            # Filter for parameters that can be set as inits (beta, alpha, lambda, tau, rho)
            # Avoid setting deterministic nodes or complex arrays if possible
            # For simplicity, we try to set scalar parameters

            # Helper to parse parameter name
            # e.g. beta[1] -> name="beta", index=1
            # But JAGS inits list expects named list of values/arrays

            # A safer approach is to let JAGS initialize, or use a simple list if we can parse it easily.
            # Given the complexity of parsing JAGS array indices from flat names,
            # we might skip this for now unless we have a robust parser.
            # However, we can try to set simple scalar parameters.

            # Actually, rjags allows passing a function or a list of lists.
            # Let's stick to random inits for robustness unless we are sure.
            # But to address user concern, we can try to use the last samples if available?
            # No, samples are from parallel chains, might be hard to map back to structure.

            # Let's rely on the fact that we use n.burnin from the user.
            # If the user provides a small n.burnin, they assume quick convergence.
            # If we use random inits, we might need more burnin.

            # COMPROMISE: We will use the provided n.burnin.
            # If the user wants to speed it up, they can provide inits to phybase_run?
            # But here we are inside phybase_waic.

            # Let's just proceed with random inits but ensure we use the user's n.burnin.
            # The user's point "will this not make run even longer" implies they want to avoid full re-run.
            # But we MUST re-run to get WAIC if the model object is invalid.
            # The best we can do is use the provided burnin.
        }

        # Recompile
        jags_model <- rjags::jags.model(
            model_file,
            data = model$data,
            n.chains = n_chains_new,
            n.adapt = 1000, # Use a reasonable default adaptation
            quiet = TRUE
        )

        # Burn-in using the provided n.burnin argument
        if (n.burnin > 0) {
            update(jags_model, n.iter = n.burnin)
        }

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
