#' Run a Phylogenetic Bayesian Structural Equation model (PhyBaSE)
#'
#' @param data A named list of data for JAGS
#' @param tree A phylogenetic tree of class \code{"phylo"}
#' @param equations A list of model formulas
#' @param n.chains Number of MCMC chains
#' @param n.iter Total number of MCMC iterations
#' @param n.burnin Number of burn-in iterations
#' @param n.thin Thinning rate
#' @param DIC Whether to compute DIC (default TRUE)
#' @param pD Whether to compute effective number of parameters (default FALSE)
#' @param n.iter.pd Number of iterations for pD estimation (if used)
#' @param n.adapt Number of adaptation iterations
#' @param quiet Suppress JAGS output
#' @param dsep If TRUE, monitor only the first beta in each structural equation (for d-separation testing)
#'
#' @return JAGS model output as returned by \code{R2jags}
#' @export

phybase_run <- function(data, tree, equations,
                        n.chains = 3,
                        n.iter = 2000,
                        n.burnin = floor(n.iter / 2),
                        n.thin = max(1, floor((n.iter - n.burnin) / 1000)),
                        DIC = TRUE, pD = FALSE, n.iter.pd = NULL,
                        n.adapt = 100,
                        quiet = FALSE,
                        dsep = FALSE) {

  # Standardize tree height and build VCV and ID
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  VCV <- ape::vcv.phylo(tree)
  ID <- diag(nrow(VCV))
  N <- nrow(VCV)

  # Generate JAGS model string
  jags_model_string <- phybase_model(equations)

  # Write the JAGS model to a temporary file
  model_file <- tempfile(fileext = ".jg")
  writeLines(jags_model_string, model_file)

  # Extract all lines from model
  parameter_lines <- unlist(strsplit(jags_model_string, "\n"))

  # ---- Choose parameters to monitor ----
  if (dsep) {
    # Only the first beta in each muX[i] <- ... line
    struct_lines <- grep("^\\s*mu\\w+\\[i\\] <-", parameter_lines, value = TRUE)
    first_betas <- sapply(struct_lines, function(line) {
      betas <- regmatches(line, gregexpr("beta\\w+", line))[[1]]
      if (length(betas) > 0) betas[1] else NA
    })
    parameters.to.save <- unique(na.omit(first_betas))
  } else {
    # Monitor only beta, alpha, lambda, and tau parameters
    beta_params   <- grep("^\\s*beta\\w+", parameter_lines, value = TRUE)
    alpha_params  <- grep("^\\s*alpha\\w+", parameter_lines, value = TRUE)
    lambda_params <- grep("^\\s*lambda\\w+", parameter_lines, value = TRUE)
    tau_params    <- grep("^\\s*tau\\w+", parameter_lines, value = TRUE)

    extract_param_names <- function(lines) {
      lines <- grep("(<-|~)", lines, value = TRUE)
      gsub("^\\s*(\\w+)\\s*(<-|~).*", "\\1", lines)
    }

    parameters.to.save <- unique(c(
      extract_param_names(beta_params),
      extract_param_names(alpha_params),
      extract_param_names(lambda_params),
      extract_param_names(tau_params)
    ))
  }

  # ---- Ensure all response variables used in model are in the data ----
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

  # Add VCV, ID and N to data list
  data <- append(data, list(VCV = VCV, ID = ID, N = N))

  # ---- Run JAGS ----
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
    parameters.to.save = parameters.to.save
  )

  return(jags_model)
}
