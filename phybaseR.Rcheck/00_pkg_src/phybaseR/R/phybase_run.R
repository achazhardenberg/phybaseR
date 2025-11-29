#' Run a Phylogenetic Bayesian Structural Equation model (PhyBaSE)
#'
#' @param data A named list of data for JAGS.
#' @param tree A single phylogenetic tree of class \code{"phylo"} or a list of trees of class \code{"multiPhylo"} (i.e., a list of \code{"phylo"} objects). If multiple trees are provided, phylogenetic uncertainty is incorporated by sampling among them during model fitting.
#' @param equations A list of model formulas describing the structural equation model.
#' @param monitor Optional character vector of parameter names to monitor. If \code{NULL}, parameters will be selected automatically based on model structure.
#' @param n.chains Number of MCMC chains (default = 3).
#' @param n.iter Total number of MCMC iterations (default = 2000).
#' @param n.burnin Number of burn-in iterations (default = n.iter / 2).
#' @param n.thin Thinning rate (default = max(1, floor((n.iter - n.burnin) / 1000))).
#' @param DIC Logical; whether to compute DIC using \code{dic.samples()} (default = TRUE).
#' @param WAIC Logical; whether to sample values for WAIC and deviance (default = FALSE).
#' @param n.adapt Number of adaptation iterations (default = 100).
#' @param quiet Logical; suppress JAGS output (default = FALSE).
#' @param dsep Logical; if \code{TRUE}, monitor only the first beta in each structural equation (used for d-separation testing).
#'
#' @return A list of class \code{"phybase"} with model output and diagnostics.
#' @export
#' @importFrom ape vcv.phylo branching.times
phybase_run <- function(data, tree, equations,
                        monitor = NULL,
                        n.chains = 3,
                        n.iter = 2000,
                        n.burnin = floor(n.iter / 2),
                        n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                        DIC = TRUE,
                        WAIC = FALSE,
                        n.adapt = 100,
                        quiet = FALSE,
                        dsep = FALSE) {

  # Handle tree(s)
  is_multiple <- inherits(tree, "multiPhylo") || (is.list(tree) && all(sapply(tree, inherits, "phylo")))

  if (is_multiple) {
    for (i in seq_along(tree)) {
      tree[[i]]$edge.length <- tree[[i]]$edge.length / max(ape::branching.times(tree[[i]]))
    }
    multiVCV <- sapply(tree, ape::vcv.phylo, simplify = "array")
    ID <- diag(nrow(multiVCV[, , 1]))
    VCV <- NULL
    N <- dim(multiVCV)[1]
    data$multiVCV <- multiVCV
    data$Ntree <- dim(multiVCV)[3]
  } else {
    tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
    VCV <- ape::vcv.phylo(tree)
    ID <- diag(nrow(VCV))
    N <- nrow(VCV)
    data$VCV <- VCV
  }

  data$ID <- ID
  data$N <- N

  # JAGS model code
  model_string <- phybase_model(equations, multi.tree = is_multiple)
  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # Monitor parameters
  if (is.null(monitor)) {
    lines <- unlist(strsplit(model_string, "\n"))
    if (dsep) {
      mu_lines <- grep("^\\s*mu\\w+\\[i\\] <-", lines, value = TRUE)
      first_betas <- sapply(mu_lines, function(line) {
        betas <- regmatches(line, gregexpr("beta\\w+", line))[[1]]
        if (length(betas) > 0) betas[1] else NA
      })
      monitor <- unique(na.omit(first_betas))
    } else {
      extract_names <- function(pattern) {
        out <- grep(pattern, lines, value = TRUE)
        out <- grep("(<-|~)", out, value = TRUE)
        gsub("^\\s*(\\w+)\\s*(<-|~).*", "\\1", out)
      }
      monitor <- unique(c(
        extract_names("^\\s*beta"),
        extract_names("^\\s*alpha"),
        extract_names("^\\s*lambda"),
        extract_names("^\\s*tau")
      ))
    }
  }

  # Add response variables
  matches <- regmatches(model_string, gregexpr("\\b(\\w+)\\[1:N\\]\\s*~", model_string))[[1]]
  response_vars <- unique(gsub("\\[1:N\\]\\s*~", "", matches))
  for (v in response_vars) {
    if (!v %in% names(data)) {
      base <- sub("[0-9]+$", "", v)
      if (base %in% names(data)) data[[v]] <- data[[base]]
    }
  }

  # Compile model
  model <- rjags::jags.model(model_file, data = data, n.chains = n.chains, n.adapt = n.adapt, quiet = quiet)
  update(model, n.iter = n.burnin)

  # Sample posterior
  samples <- rjags::coda.samples(model, variable.names = monitor, n.iter = n.iter - n.burnin, thin = n.thin)

  # Initialize result object
  result <- list(
    model = model,
    samples = samples,
    monitor = monitor,
    modfile = model_file
  )

  # Add DIC
  if (DIC) {
    result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
  }

  # Optionally sample for WAIC
  if (WAIC) {
    waic_samples <- rjags::jags.samples(model,
                                        c("WAIC", "deviance"),
                                        type = "mean",
                                        n.iter = n.iter - n.burnin,
                                        thin = n.thin)
    waic_samples$p_waic <- waic_samples$WAIC
    waic_samples$waic <- waic_samples$deviance + waic_samples$p_waic
    tmp <- sapply(waic_samples, sum)
    result$WAIC <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]), 1)
  }

  class(result) <- "phybase"
  return(result)
}
