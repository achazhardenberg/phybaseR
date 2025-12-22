#' @title Family Definition Generics and Default Methods
#' @description S3 generics and methods for distribution families in because.
#' This enables custom distributions to be added by defining S3 methods.
#' @name family_definitions
NULL

#' Generate JAGS likelihood code for a distribution family
#'
#' This generic allows S3 dispatch to generate the appropriate JAGS
#' likelihood code for any distribution family.
#'
#' @param family A family object (created by \code{\link{get_family_object}} or custom constructor)
#' @param response Character string; name of the response variable
#' @param predictors Character vector; names of predictor variables (may be NULL)
#' @param suffix Character string; suffix for variable names (e.g., "1" for multiple responses)
#' @param has_structure Logical; whether the model includes a structure (e.g., phylogenetic)
#' @param link Character string; link function ("identity", "log", "logit")
#' @param ... Additional arguments passed to methods
#'
#' @return A list with:
#' \describe{
#'   \item{likelihood_code}{Character vector of JAGS likelihood statements}
#'   \item{prior_code}{Character vector of JAGS prior statements (or NULL)}
#'   \item{data_requirements}{Character vector of required data elements (or NULL)}
#' }
#'
#' @export
jags_family_likelihood <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    UseMethod("jags_family_likelihood")
}

#' @export
jags_family_likelihood.default <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    # Get family name
    fam_name <- if (is.character(family)) family else family$family
    stop(paste0(
        "No jags_family_likelihood method defined for family '",
        fam_name,
        "'.\n",
        "Use because_family() to create custom family methods, or check if ",
        "a module is required (e.g., because.occupancy)."
    ))
}

#' @export
jags_family_likelihood.because_family_gaussian <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    # Standard normal likelihood
    mu_var <- paste0("mu_", response, suffix)
    tau_var <- paste0("tau_e_", response, suffix)

    likelihood_code <- paste0(
        "    ",
        response,
        "[i] ~ dnorm(",
        mu_var,
        "[i], ",
        tau_var,
        ")"
    )

    prior_code <- paste0("  ", tau_var, " ~ dgamma(1, 1)")

    list(
        likelihood_code = likelihood_code,
        prior_code = prior_code,
        data_requirements = NULL
    )
}

#' @export
jags_family_likelihood.because_family_binomial <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "logit",
    ...
) {
    # Bernoulli likelihood with logit link
    p_var <- paste0("p_", response, suffix)
    mu_var <- paste0("mu_", response, suffix)

    likelihood_code <- c(
        paste0("    ", p_var, "[i] <- ilogit(", mu_var, "[i])"),
        paste0("    ", response, "[i] ~ dbern(", p_var, "[i])")
    )

    list(
        likelihood_code = likelihood_code,
        prior_code = NULL,
        data_requirements = NULL
    )
}

#' @export
jags_family_likelihood.because_family_poisson <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "log",
    ...
) {
    # Poisson likelihood with log link
    lambda_var <- paste0("lambda_", response, suffix)
    mu_var <- paste0("mu_", response, suffix)

    likelihood_code <- c(
        paste0("    ", lambda_var, "[i] <- exp(", mu_var, "[i])"),
        paste0("    ", response, "[i] ~ dpois(", lambda_var, "[i])")
    )

    list(
        likelihood_code = likelihood_code,
        prior_code = NULL,
        data_requirements = NULL
    )
}

#' Create a family object for a given distribution
#'
#' Returns a family object that can be dispatched on via S3.
#'
#' @param name Character string; family name (e.g., "gaussian", "binomial")
#' @return A family object with appropriate S3 class
#' @export
get_family <- function(name) {
    # Normalize name
    name <- tolower(name)

    # Map aliases
    name <- switch(name, "normal" = "gaussian", "bernoulli" = "binomial", name)

    # Create family object with S3 class
    structure(
        list(family = name),
        class = c(paste0("because_family_", name), "because_family")
    )
}
