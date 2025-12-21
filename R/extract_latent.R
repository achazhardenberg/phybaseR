#' Extract Latent Parameters from a because Model
#'
#' This function extracts posterior distributions of latent variables such as
#' occupancy probability (psi), detection probability (p), and latent occupancy state (z).
#' It is designed to handle multispecies models by mapping stacked results back to species and sites.
#'
#' @param object A \code{because} model object.
#' @param type Character; the type of latent variable to extract. Currently supports "occupancy".
#' @param variables Optional character vector specifying which latent variables to extract
#'   (e.g., \code{c("psi", "p", "z")}). If \code{NULL}, extracts all available for the given type.
#'
#' @return A tidy data frame (tibble) containing:
#' \itemize{
#'   \item \code{Variable}: The base name of the latent parameter (e.g., "psi").
#'   \item \code{Response}: The response variable name in the model (e.g., "Y").
#'   \item \code{SpeciesID}: The species identifier (for multispecies models).
#'   \item \code{SiteID}: The site identifier.
#'   \item \code{Mean, SD, Q2.5, Q50, Q97.5}: Posterior summary statistics.
#' }
#'
#' @details
#' For occupancy models, the function looks for parameters starting with \code{psi_}, \code{p_}, and \code{z_}.
#' It uses the internal \code{results$data} or \code{results$species_order} to map indices back to
#' meaningful identifiers.
#'
#' @export
extract_latent <- function(object, type = "occupancy", variables = NULL) {
    if (!inherits(object, "because")) {
        stop("Object must be of class 'because'.")
    }

    if (type != "occupancy") {
        stop("Currently only 'type = \"occupancy\"' is supported.")
    }

    # Identification of occupancy variables
    dist <- object$input$distribution
    occ_vars <- names(dist)[dist == "occupancy"]

    if (length(occ_vars) == 0) {
        warning("No occupancy variables found in the model.")
        return(data.frame())
    }

    # Get summary stats
    sum_stats <- object$summary
    if (is.null(sum_stats)) {
        stop(
            "Model summary not found. Ensure the model was run with monitoring enabled."
        )
    }

    results_list <- list()

    # Get data info for mapping
    # For stacked data, we check both `data` (if preserved as DF) or `stacked_data` (newly added)
    stacked_df <- if (!is.null(object$stacked_data)) {
        object$stacked_data
    } else if (is.data.frame(object$data) && !is.null(object$data$SpeciesID)) {
        object$data
    } else {
        NULL
    }

    is_stacked <- !is.null(stacked_df)

    for (v in occ_vars) {
        target_prefixes <- if (is.null(variables)) {
            c("psi", "p", "z")
        } else {
            variables
        }

        for (prefix in target_prefixes) {
            param_base <- paste0(prefix, "_", v)

            # Find all matching parameters in summary
            # Format: param_base[i]
            all_params <- rownames(sum_stats$statistics)
            pattern <- paste0("^", param_base, "\\[(\\d+)\\]$")
            matches <- grep(pattern, all_params, value = TRUE)

            if (length(matches) == 0) {
                next
            }

            # Extract indices
            indices <- as.integer(gsub(pattern, "\\1", matches))

            for (idx in seq_along(indices)) {
                i <- indices[idx]
                p_name <- matches[idx]

                stats <- sum_stats$statistics[p_name, ]
                quants <- sum_stats$quantiles[p_name, ]

                row_data <- data.frame(
                    Variable = prefix,
                    Response = v,
                    Index = i,
                    Mean = stats["Mean"],
                    SD = stats["SD"],
                    Q2.5 = quants["2.5%"],
                    Q50 = quants["50%"],
                    Q97.5 = quants["97.5%"],
                    stringsAsFactors = FALSE
                )

                # Add Species and Site info if available
                if (is_stacked && i <= nrow(stacked_df)) {
                    row_data$SpeciesID <- as.character(stacked_df$SpeciesID[i])
                    row_data$SiteID <- as.character(stacked_df$SiteID[i])
                } else if (
                    !is.null(object$species_order) &&
                        i <= length(object$species_order)
                ) {
                    row_data$SpeciesID <- object$species_order[i]
                }

                results_list[[length(results_list) + 1]] <- row_data
            }
        }
    }

    if (length(results_list) == 0) {
        warning("No latent occupancy parameters found in the summary.")
        return(data.frame())
    }

    res <- do.call(rbind, results_list)
    rownames(res) <- NULL

    # Reorder columns for readability
    cols <- c(
        "Variable",
        "Response",
        "SpeciesID",
        "SiteID",
        "Mean",
        "SD",
        "Q2.5",
        "Q50",
        "Q97.5"
    )
    existing_cols <- intersect(cols, names(res))
    res[, existing_cols]
}
