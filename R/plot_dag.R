#' Plot DAG from Equations or Fitted Model
#'
#' Visualizes the Directed Acyclic Graph (DAG) implied by a set of equations
#' or a fitted `because` model. If a fitted model is provided:
#' * Path coefficients are displayed on the edges.
#' * Edges are colored black if the parameter's 95% Credible Interval excludes zero, and grey otherwise.
#' * Edge thickness scales with the absolute effect size.
#'
#' @param x A list of formulas (equations), a `because` model object, or a list of these.
#' @param layout The layout algorithm to use (default "nicely"). See \code{\link[ggdag]{ggdag}}.
#' @param latent Character vector of latent variable names. Overrides the model's latent variables if provided.
#' @param node_size Size of the nodes (default 14).
#' @param node_color Color of the node border (default "black").
#' @param node_fill Color of the node interior (default "white").
#' @param text_size Size of the labels (default 4).
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
    latent = NULL,
    node_size = 14,
    node_color = "black",
    node_fill = "white",
    text_size = 4,
    edge_width_range = c(0.5, 2),
    show_coefficients = TRUE,
    ...
) {
    # Check dependencies
    if (
        !requireNamespace("dagitty", quietly = TRUE) ||
            !requireNamespace("ggdag", quietly = TRUE) ||
            !requireNamespace("ggraph", quietly = TRUE) ||
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

        # 1. Extract Equations and Latent Info
        current_latent <- latent

        if (inherits(obj, "because")) {
            eqs <- obj$parameter_map$equations %||%
                obj$input$equations %||%
                obj$parameter_map$equations # Fallback
            if (is.null(eqs)) {
                stop(
                    "Could not find equations in 'because' object. Please refit the model or manually pass equations."
                )
            }
            # Use object's latent vars if not overridden
            if (is.null(current_latent)) {
                current_latent <- obj$input$latent
            }
        } else {
            # List of formulas
            eqs <- obj
        }

        # 2. Convert to dagitty syntax
        induced_cors <- NULL
        if (inherits(obj, "because")) {
            induced_cors <- obj$induced_correlations %||%
                obj$input$induced_correlations
        }

        dag_str <- equations_to_dag_string(eqs, induced_cors)
        dag_obj <- dagitty::dagitty(dag_str)

        # 3. Tidy it up using ggdag
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = layout)

        # Extract data frame for plotting (nodes + edges flattened)
        dag_data <- dplyr::as_tibble(tidy_dag)

        # Mark node types (Observed vs Latent)
        dag_data$type <- "Observed"
        if (!is.null(current_latent)) {
            dag_data$type[dag_data$name %in% current_latent] <- "Latent"
        }

        # Initialize columns for ALL edges
        dag_data$edge_type <- NA # "->", "<->"
        dag_data$weight_abs <- 0.5
        dag_data$edge_label <- NA
        dag_data$val <- NA

        # Get edges metadata from dagitty
        edges_meta <- dagitty::edges(dag_obj) # v, w, e

        # 4. If fitted model, extract coefficients/correlations
        stats <- NULL
        quantiles <- NULL
        if (inherits(obj, "because") && !is.null(obj$summary)) {
            stats <- obj$summary$statistics
            quantiles <- obj$summary$quantiles
        }

        # Process edges to find parameters
        if ("to" %in% names(dag_data) && nrow(edges_meta) > 0) {
            edge_rows <- which(!is.na(dag_data$to))

            for (idx in edge_rows) {
                v <- dag_data$name[idx]
                w <- dag_data$to[idx]

                # Find corresponding edge via match
                meta_match <- which(
                    (edges_meta$v == v & edges_meta$w == w) |
                        (edges_meta$e == "<->" &
                            edges_meta$v == w &
                            edges_meta$w == v)
                )

                if (length(meta_match) > 0) {
                    m_idx <- meta_match[1]
                    e_type <- edges_meta$e[m_idx]
                    dag_data$edge_type[idx] <- e_type

                    if (!is.null(stats)) {
                        val <- NA
                        sig <- FALSE
                        pname <- NULL

                        if (e_type == "->") {
                            # Beta: beta_w_v
                            try_pname <- paste0("beta_", w, "_", v)
                            if (try_pname %in% rownames(stats)) {
                                pname <- try_pname
                            }
                        } else if (e_type == "<->") {
                            # Rho: rho_v_w or rho_w_v
                            pname1 <- paste0("rho_", v, "_", w)
                            pname2 <- paste0("rho_", w, "_", v)
                            if (pname1 %in% rownames(stats)) {
                                pname <- pname1
                            } else if (pname2 %in% rownames(stats)) {
                                pname <- pname2
                            }
                        }

                        if (!is.null(pname)) {
                            val <- stats[pname, "Mean"]

                            # Check significance if quantiles available
                            if (
                                !is.null(quantiles) &&
                                    pname %in% rownames(quantiles)
                            ) {
                                lower <- quantiles[pname, "2.5%"]
                                upper <- quantiles[pname, "97.5%"]
                                if (sign(lower) == sign(upper)) {
                                    sig <- TRUE
                                }
                            } else {
                                # Fallback if no quantiles (unlikely for MCMC)
                                sig <- TRUE # Default to black if unsure? Or Grey? Let's say black.
                            }
                        }

                        if (!is.na(val)) {
                            dag_data$val[idx] <- val
                            dag_data$weight_abs[idx] <- abs(val)
                            dag_data$edge_label[idx] <- round(val, 2)
                            # Store significance for coloring
                            dag_data$significant[idx] <- sig
                        }
                    }
                }
            }
        }

        # Fill defaults for plotting if missing
        dag_data$weight_abs[
            is.na(dag_data$weight_abs) & !is.na(dag_data$to)
        ] <- 0.2
        # Default edge type to -> if not found (robustness)
        dag_data$edge_type[
            is.na(dag_data$edge_type) & !is.na(dag_data$to)
        ] <- "->"
        # Default significance to TRUE (black) for structural edges without parameters
        if (!"significant" %in% names(dag_data)) {
            dag_data$significant <- TRUE
        }
        dag_data$significant[is.na(dag_data$significant)] <- TRUE
        dag_data$significant <- as.character(dag_data$significant)

        dag_data$model_label <- label

        if (is.null(combined_dag_data)) {
            combined_dag_data <- dag_data
        } else {
            combined_dag_data <- rbind(combined_dag_data, dag_data)
        }
    }

    # Ensure type is a factor
    combined_dag_data$type <- factor(
        combined_dag_data$type,
        levels = c("Observed", "Latent")
    )

    # Build Plot
    p <- ggplot2::ggplot(
        combined_dag_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    ) +
        # Nodes using standard geom_point to ensure single border
        ggplot2::geom_point(
            ggplot2::aes(shape = type),
            size = node_size,
            color = node_color,
            fill = node_fill,
            stroke = 1.5,
            show.legend = FALSE
        ) +
        ggdag::geom_dag_text(size = text_size, color = node_color) +
        ggdag::theme_dag() +
        # Map shapes: Square (22) for Observed, Circle (21) for Latent
        ggplot2::scale_shape_manual(
            values = c(Observed = 22, Latent = 21),
            guide = "none" # Remove legend
        )

    if (length(unique(combined_dag_data$model_label)) > 1) {
        p <- p + ggplot2::facet_wrap(~model_label)
    }

    # Edges Layer: Split into Directed and Bidirected
    if ("edge_type" %in% names(combined_dag_data)) {
        # 1. Directed Edges (Solid)
        p <- p +
            ggdag::geom_dag_edges_link(
                data = function(x) dplyr::filter(x, edge_type == "->"),
                mapping = ggplot2::aes(
                    edge_width = weight_abs,
                    edge_alpha = weight_abs,
                    edge_colour = significant, # Map color to significance
                    label = edge_label
                ),
                angle_calc = "along",
                label_dodge = ggplot2::unit(3, "mm")
            )

        # 2. Bidirected Edges (Dashed, Curved, Double Arrow)
        p <- p +
            ggdag::geom_dag_edges_arc(
                data = function(x) dplyr::filter(x, edge_type == "<->"),
                mapping = ggplot2::aes(
                    edge_width = weight_abs,
                    edge_alpha = weight_abs, # Or keep constant? User asked for rho on top.
                    edge_colour = significant, # Map color to significance
                    label = edge_label
                ),
                curvature = 0.3,
                edge_linetype = "dashed",
                angle_calc = "along",
                label_dodge = ggplot2::unit(3, "mm"),
                arrow = ggplot2::arrow(
                    length = ggplot2::unit(2.5, "mm"),
                    type = "closed",
                    ends = "both"
                )
            )

        # Scales
        p <- p +
            ggraph::scale_edge_width_continuous(
                range = edge_width_range,
                guide = "none"
            ) +
            ggraph::scale_edge_alpha_continuous(
                range = c(0.3, 1),
                guide = "none"
            ) +
            ggraph::scale_edge_colour_manual(
                values = c("TRUE" = "black", "FALSE" = "grey70"),
                guide = "none"
            )
    } else {
        p <- p + ggdag::geom_dag_edges()
    }

    return(p)
}


#' Convert Equations List to DAGitty String
#' @param equations List of formulas
#' @param induced_cors List of character vectors (pairs) for bidirected edges
#' @noRd
equations_to_dag_string <- function(equations, induced_cors = NULL) {
    edges <- c()
    for (eq in equations) {
        resp <- all.vars(eq)[1]
        predictors <- all.vars(eq)[-1]
        for (pred in predictors) {
            edges <- c(edges, paste(resp, "<-", pred))
        }
    }

    if (!is.null(induced_cors)) {
        for (pair in induced_cors) {
            if (length(pair) == 2) {
                # Add bidirected edge A <-> B
                edges <- c(edges, paste(pair[1], "<->", pair[2]))
            }
        }
    }

    return(paste("dag {", paste(edges, collapse = "; "), "}"))
}
