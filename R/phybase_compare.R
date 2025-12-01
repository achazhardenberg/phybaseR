#' Compare Multiple PhyBaSE Models in Parallel
#'
#' Fits multiple PhyBaSE models in parallel and returns a comparison table
#' based on WAIC and DIC.
#'
#' @param model_specs A named list where each element is a list containing the
#'   \code{equations} argument for \code{\link{phybase_run}}.
#'   Example: \code{list(m1 = list(equations = list(Y ~ X)), m2 = list(equations = list(Y ~ X + Z)))}
#' @param data A list containing the data for the models.
#' @param tree A phylo object or list of phylo objects.
#' @param n.cores Number of cores to use for parallel execution. Default is 1.
#' @param cl Optional cluster object created by \code{parallel::makeCluster()}.
#'   If \code{NULL} and \code{n.cores > 1}, a cluster is created and stopped automatically.
#' @param ... Additional arguments passed to \code{\link{phybase_run}} (e.g., \code{n.iter}, \code{n.burnin}).
#'   Note: \code{parallel} will be forced to \code{FALSE} for individual runs to avoid nested parallelism.
#'
#' @return A list containing:
#'   \item{results}{A list of fitted \code{phybase} objects.}
#'   \item{comparison}{A data frame comparing models by WAIC and DIC (if computed).}
#'
#' @examples
#' \dontrun{
#'   # Define models
#'   models <- list(
#'     m1 = list(equations = list(Y ~ X)),
#'     m2 = list(equations = list(Y ~ X + Z))
#'   )
#'
#'   # Run comparison
#'   comp <- phybase_compare(models, data, tree, n.cores = 2, n.iter = 1000)
#'   print(comp$comparison)
#' }
#'
#' @export
#' @importFrom parallel makeCluster stopCluster parLapply clusterExport clusterEvalQ
phybase_compare <- function(
    model_specs,
    data,
    tree,
    n.cores = 1,
    cl = NULL,
    ...
) {
    # Input validation
    if (!is.list(model_specs) || is.null(names(model_specs))) {
        stop("model_specs must be a named list of model specifications.")
    }

    # Capture dot arguments
    dot_args <- list(...)

    # Function to run a single model
    run_model <- function(name, spec, data, tree, dot_args) {
        # Prepare arguments
        args <- list(
            data = data,
            tree = tree,
            equations = spec$equations,
            parallel = FALSE # Force sequential for inner runs
        )

        # Add dot arguments (override if present)
        args <- c(args, dot_args)

        # Ensure parallel is FALSE even if passed in dots
        args$parallel <- FALSE

        # Run model
        do.call(phybaseR::phybase_run, args)
    }

    # Execution
    if (n.cores > 1) {
        message(sprintf(
            "Running %d models in parallel on %d cores...",
            length(model_specs),
            n.cores
        ))

        # Setup cluster
        if (is.null(cl)) {
            cl <- parallel::makeCluster(n.cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
        }

        # Export necessary objects
        parallel::clusterExport(
            cl,
            c("data", "tree", "dot_args", "run_model"),
            envir = environment()
        )

        # Load package on workers
        parallel::clusterEvalQ(cl, {
            library(phybaseR)
        })

        # Run in parallel
        results <- parallel::parLapply(cl, names(model_specs), function(name) {
            run_model(name, model_specs[[name]], data, tree, dot_args)
        })
        names(results) <- names(model_specs)
    } else {
        message("Running models sequentially...")
        results <- lapply(names(model_specs), function(name) {
            run_model(name, model_specs[[name]], data, tree, dot_args)
        })
        names(results) <- names(model_specs)
    }

    # Compile comparison table
    waic_vals <- sapply(results, function(x) {
        if (!is.null(x$WAIC)) x$WAIC["waic"] else NA
    })

    dic_vals <- sapply(results, function(x) {
        if (!is.null(x$DIC)) sum(x$DIC$deviance) + sum(x$DIC$penalty) else NA
    })

    comparison <- data.frame(
        Model = names(results),
        WAIC = waic_vals,
        DIC = dic_vals,
        stringsAsFactors = FALSE
    )

    # Calculate Delta and Weights if values exist
    if (!all(is.na(comparison$WAIC))) {
        min_waic <- min(comparison$WAIC, na.rm = TRUE)
        comparison$Delta_WAIC <- comparison$WAIC - min_waic

        # Akaike weights for WAIC
        exp_delta <- exp(-0.5 * comparison$Delta_WAIC)
        comparison$Weight_WAIC <- exp_delta / sum(exp_delta, na.rm = TRUE)
    }

    # Sort by WAIC if available, else DIC
    if (!all(is.na(comparison$WAIC))) {
        comparison <- comparison[order(comparison$WAIC), ]
    } else if (!all(is.na(comparison$DIC))) {
        comparison <- comparison[order(comparison$DIC), ]
    }

    list(results = results, comparison = comparison)
}
