#' Plot DAG from Equations or Fitted Model
#'
#' Visualizes the Directed Acyclic Graph (DAG) implied by a set of equations
#' or a fitted `because` model. If a fitted model is provided, it can display
#' path coefficients (standardized if available) on the edges.
#'
#' @param x A list of formulas (equations), a `because` model object, or a list of these.
#' @param layout The layout algorithm to use (default "nicely"). See \code{\link[ggdag]{ggdag}}.
#' @param node_size Size of the nodes (default 10).
#' @param text_size Size of the labels (default 5).
#' @param edge_width_range Vector of length 2 defining the range of arrow widths (min, max) based on effect size.
#' @param show_coefficients Logical; whether to print coefficient values on edges (only for fitted models).
#' @param ... Additional arguments passed to \code{\link[ggdag]{expand_plot}}.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom stats terms formula
plot_dag <- function(
    x,
    layout = "nicely",
    node_size = 12,
    text_size = 4,
    edge_width_range = c(0.5, 2),
    show_coefficients = TRUE,
    ...
) {
    # Check dependencies
    if (
        !requireNamespace("dagitty", quietly = TRUE) ||
            !requireNamespace("ggdag", quietly = TRUE) ||
            !requireNamespace("ggraph", quietly = TRUE) || # Explicitly check ggraph
            !requireNamespace("ggplot2", quietly = TRUE) ||
            !requireNamespace("dplyr", quietly = TRUE)
    ) {
        stop(
            "Packages 'dagitty', 'ggdag', 'ggraph', 'ggplot2', and 'dplyr' are required for plot_dag."
        )
    }

    # Normalize input to a list of objects
    if (
        inherits(x, "because") ||
            inherits(x, "list") && all(sapply(x, inherits, "formula"))
    ) {
        x <- list(Model = x)
    }

    # Build Tidy DAG Data Frame
    combined_dag_data <- NULL

    for (i in seq_along(x)) {
        obj <- x[[i]]
        label <- names(x)[i]
        if (is.null(label)) {
            label <- paste("Model", i)
        }

        # 1. Extract Equations
        if (inherits(obj, "because")) {
            eqs <- obj$parameter_map$equations %||% obj$input$equations # Fallback if map missing
            if (is.null(eqs)) {
                stop("Could not find equations in 'because' object.")
            }
        } else {
            # List of formulas
            eqs <- obj
        }

        # 2. Convert to dagitty syntax
        dag_str <- equations_to_dag_string(eqs)
        dag_obj <- dagitty::dagitty(dag_str)

        # 3. Tidy it up using ggdag
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = layout)

        # Extract data frame for plotting (nodes + edges flattened)
        dag_data <- dplyr::as_tibble(tidy_dag)

        # 4. If fitted model, extract coefficients
        if (inherits(obj, "because") && !is.null(obj$summary)) {
            stats <- obj$summary$statistics
            edges <- dagitty::edges(dag_obj) # v, w

            edges$beta <- NA
            edges$label <- NA
            edges$val <- NA # Initialize

            for (r in seq_len(nrow(edges))) {
                predictor <- edges$v[r]
                response <- edges$w[r]

                # Try both beta_Y_X and other patterns?
                # Standard in because: beta_Response_Predictor
                param_name <- paste0("beta_", response, "_", predictor)

                if (param_name %in% rownames(stats)) {
                    val <- stats[param_name, "Mean"]
                    edges$beta[r] <- abs(val)
                    edges$val[r] <- val
                    edges$label[r] <- round(val, 2)
                }
            }

            # Join to dag_data
            # dag_data has 'name' (start) and 'to' (end)
            if ("to" %in% names(dag_data)) {
                dag_data$beta_abs <- 0.5
                dag_data$edge_label <- NA
                dag_data$val <- NA

                edge_rows <- which(!is.na(dag_data$to))
                for (idx in edge_rows) {
                    v <- dag_data$name[idx]
                    w <- dag_data$to[idx]
                    match <- which(edges$v == v & edges$w == w)
                    if (length(match) > 0) {
                        dag_data$beta_abs[idx] <- edges$beta[match]
                        dag_data$edge_label[idx] <- edges$label[match]
                        dag_data$val[idx] <- edges$val[match]
                    }
                }
                # Replace NA betas/alphas for plotting
                dag_data$beta_abs[
                    is.na(dag_data$beta_abs) & !is.na(dag_data$to)
                ] <- 0.2
            }
        }

        dag_data$model_label <- label

        if (is.null(combined_dag_data)) {
            combined_dag_data <- dag_data
        } else {
            combined_dag_data <- rbind(combined_dag_data, dag_data)
        }
    }

    # Build Plot
    p <- ggplot2::ggplot(
        combined_dag_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    ) +
        ggdag::geom_dag_node(size = node_size, color = "grey85") +
        ggdag::geom_dag_text(size = text_size) +
        ggdag::theme_dag()

    if (length(unique(combined_dag_data$model_label)) > 1) {
        p <- p + ggplot2::facet_wrap(~model_label)
    }

    if (show_coefficients && "val" %in% names(combined_dag_data)) {
        p <- p +
            ggdag::geom_dag_edges(
                ggplot2::aes(
                    edge_width = beta_abs,
                    edge_alpha = beta_abs,
                    label = edge_label
                )
            ) +
            # Note: geom_dag_edges typically handles arrows automatically.
            # If label positioning is poor, consider geom_dag_edge_text separately,
            # but direct label mapping is supported in ggraph/ggdag for edges.
            ggraph::scale_edge_width_continuous(
                range = edge_width_range,
                guide = "none"
            ) +
            ggraph::scale_edge_alpha_continuous(
                range = c(0.3, 1),
                guide = "none"
            )
    } else {
        p <- p + ggdag::geom_dag_edges()
    }

    return(p)
}


#' Convert Equations List to DAGitty String
#' @noRd
equations_to_dag_string <- function(equations) {
    edges <- c()
    for (eq in equations) {
        resp <- all.vars(eq)[1]
        predictors <- all.vars(eq)[-1]

        # Check for random effects or special terms and stripping them might be needed
        # but all.vars() usually gets basic variables.
        # Note: If intercept only (Y~1), predictors is empty.

        for (pred in predictors) {
            # Check against pure numbers (intercepts)
            # all.vars shouldn't capture "1".
            edges <- c(edges, paste(resp, "<-", pred))
        }
    }
    return(paste("dag {", paste(edges, collapse = "; "), "}"))
}
