#' Run a Phylogenetic Bayesian Structural Equation model (PhyBaSE)
#'
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
#' @param DIC Logical; whether to compute DIC (default = TRUE).
#' @param pD Logical; whether to compute effective number of parameters (default = FALSE).
#' @param n.iter.pd Number of iterations for pD estimation (optional).
#' @param n.adapt Number of adaptation iterations (default = 100).
#' @param quiet Logical; suppress JAGS output (default = FALSE).
#' @param dsep Logical; if \code{TRUE}, monitor only the first beta in each structural equation (used for d-separation testing).
#'
#' @return Output from \code{R2jags::jags()}, including posterior samples and diagnostics.
#' @export
#' @importFrom ape vcv.phylo branching.times

phybase_run <- function(data, tree, equations,
                        monitor = NULL,
                        n.chains = 3,
                        n.iter = 2000,
                        n.burnin = floor(n.iter / 2),
                        n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                        DIC = TRUE, pD = FALSE, n.iter.pd = NULL,
                        n.adapt = 100,
                        quiet = FALSE,
                        dsep = FALSE) {

  # Check if multiple trees provided
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

  # Generate JAGS model string
  if (is_multiple) {
    jags_model_string <- phybase_model(equations, multi.tree = TRUE)
  } else {
    jags_model_string <- phybase_model(equations, multi.tree = FALSE)
  }

  # Write model to temp file
  model_file <- tempfile(fileext = ".jg")
  writeLines(jags_model_string, model_file)

  # Determine parameters to monitor
  if (is.null(monitor)) {
    parameter_lines <- unlist(strsplit(jags_model_string, "\n"))
    if (dsep) {
      struct_lines <- grep("^\\s*mu\\w+\\[i\\] <-", parameter_lines, value = TRUE)
      first_betas <- sapply(struct_lines, function(line) {
        betas <- regmatches(line, gregexpr("beta\\w+", line))[[1]]
        if (length(betas) > 0) betas[1] else NA
      })
      monitor <- unique(na.omit(first_betas))
    } else {
      beta_params   <- grep("^\\s*beta\\w+", parameter_lines, value = TRUE)
      alpha_params  <- grep("^\\s*alpha\\w+", parameter_lines, value = TRUE)
      lambda_params <- grep("^\\s*lambda\\w+", parameter_lines, value = TRUE)
      tau_params    <- grep("^\\s*tau\\w+", parameter_lines, value = TRUE)

      extract_param_names <- function(lines) {
        lines <- grep("(<-|~)", lines, value = TRUE)
        gsub("^\\s*(\\w+)\\s*(<-|~).*", "\\1", lines)
      }

      monitor <- unique(c(
        extract_param_names(beta_params),
        extract_param_names(alpha_params),
        extract_param_names(lambda_params),
        extract_param_names(tau_params)
      ))
    }
  }

  # Add any missing response variables
  matches <- regmatches(jags_model_string,
                        gregexpr("\\b(\\w+)\\[1:N\\]\\s*~", jags_model_string))[[1]]
  response_vars <- unique(gsub("\\[1:N\\]\\s*~", "", matches))

  for (v in response_vars) {
    if (!v %in% names(data)) {
      base <- sub("[0-9]+$", "", v)
      if (base %in% names(data)) {
        data[[v]] <- data[[base]]
      }
    }
  }

  # Run JAGS
  jags_model <- R2jags::jags(
    model.file = model_file,
    data = data,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    n.adapt = n.adapt,
    quiet = quiet,
    DIC = DIC,
    parameters.to.save = monitor
  )

  return(jags_model)
}
