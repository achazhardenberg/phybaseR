#' Plot DAG from Equations or Fitted Model
#'
#' Visualizes the Directed Acyclic Graph (DAG) implied by a set of equations
#' or a fitted `because` model. If a fitted model is provided:
#' * Path coefficients are displayed on the edges.
#' * Edge colors are determined by the `edge_color_scheme` argument (e.g., Red/Blue for directional, Black/Grey for binary).
#' * Structural edges (without fitted data) are colored black.
#' * Edge thickness scales with the absolute effect size.
#'
#' @param x A list of formulas (equations), a `because` model object, or a list of these.
#' @param layout The layout algorithm to use (default "nicely"). See \code{\link[ggdag]{ggdag}}.
#' @param latent Character vector of latent variable names. Overrides the model's latent variables if provided.
#' @param node_size Size of the nodes (default 14).
#' @param node_color Color of the node border (default "black").
#' @param node_fill Color of the node interior (default "white").
#' @param node_stroke Thickness of the node border (default 1.5).
#' @param text_size Size of the labels (default 4).
#' @param edge_width_range Vector of length 2 defining the range of arrow widths (min, max) based on effect size.
#' @param edge_color_scheme Character; one of "directional" (default), "binary", or "monochrome".
#' "directional" colors edges red/blue/grey based on effect direction and whether the 95\% CI excludes zero.
#' "binary" colors edges black/grey based on whether the 95\% CI excludes zero (black) or includes it (grey).
#' "monochrome" colors all edges black.
#' @param show_coefficients Logical; whether to print coefficient values on edges (only for fitted models).
#' @param coords Optional named list of coordinates for the nodes, e.g. `list(A = c(1, 1), B = c(2, 2))`.
#' If provided, these will override the `layout` algorithm.
#'
#' @return A `ggplot` object that can be further customized with standard ggplot2 functions (e.g., `+ theme_...()`, `+ ggtitle(...)`).
#' @export
#' @importFrom stats terms formula
#'
#' @examples
#' \dontrun{
#' # Basic plotting
#' eq <- list(y ~ x + z, x ~ z)
#' plot_dag(eq)
#'
#' # Custom Layout
#' my_coords <- list(
#'   y = c(1, 1),
#'   x = c(0, 0),
#'   z = c(2, 0)
#' )
#' plot_dag(eq, coords = my_coords)
#' }
plot_dag <- function(
    x,
    layout = "nicely",
    latent = NULL,
    node_size = 14,
    node_color = "black",
    node_fill = "white",
    node_stroke = 1.5,
    text_size = 4,
    edge_width_range = c(0.5, 2),
    edge_color_scheme = c("directional", "binary", "monochrome"),
    show_coefficients = TRUE,
    coords = NULL
) {
    edge_color_scheme <- match.arg(edge_color_scheme)

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

        # Apply coordinates if provided
        if (!is.null(coords)) {
            if (!is.list(coords)) {
                stop(
                    "coords must be a named list of numeric vectors, e.g. list(node = c(x, y))."
                )
            }
            # Extract x and y vectors with names
            x_coords <- sapply(coords, function(v) v[1])
            y_coords <- sapply(coords, function(v) v[2])
            names(x_coords) <- names(coords)
            names(y_coords) <- names(coords)

            dagitty::coordinates(dag_obj) <- list(x = x_coords, y = y_coords)
        }

        # 3. Tidy it up using ggdag
        # If coords are provided, we should probably trust them, but ggdag generally behaves well if dagitty has coords.
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = layout)

        # Extract data frame for plotting (nodes + edges flattened)
        dag_data <- as.data.frame(dplyr::as_tibble(tidy_dag))

        # Label Processing: Wrap text (replace _ with \n)
        dag_data$label_display <- gsub("_", "\n", dag_data$name)

        # Determine Node Size based on longest label (Heuristic)
        # Assuming circular node, diameter needs to cover the rectangular text block.
        # We find the max characters in any single line of the wrapped labels.
        max_chars <- max(nchar(unlist(strsplit(dag_data$label_display, "\n"))))
        # Heuristic: base size + scaling factor * chars * text_size
        # Typically node_size=14 covers ~3 chars at text_size=4.
        # We'll use a conservative multiplier.
        calc_size <- 10 + (max_chars * 2.5 * (text_size / 4))

        # Use simple logic: manual node_size is a baseline, but we ensure it fits.
        # However, user asked to "make the boxes size the same... according to the longest".
        # So we use the calculated max size for ALL nodes stringently.
        # We will take the maximum of the default/user input and the calculated requirement.
        current_node_size <- max(node_size, calc_size)

        # Initialize columns for ALL edges
        dag_data$edge_type <- NA_character_
        dag_data$weight_abs <- 1.0 # Default to 1 (solid/opaque) for structural plots
        dag_data$val <- NA_real_
        dag_data$edge_label <- NA_character_ # Initialize edge_label column
        dag_data$significant <- NA # Logical placeholder

        # Mark node types (Observed vs Latent)
        dag_data$type <- "Observed"
        if (!is.null(current_latent)) {
            dag_data$type[dag_data$name %in% current_latent] <- "Latent"
        }

        # 3. Add coefficients if available
        # Get edges metadata from dagitty for finding parameters
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
                        sig_cat <- "default"
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

                            # Determine significance category based on scheme
                            if (
                                !is.null(quantiles) &&
                                    pname %in% rownames(quantiles)
                            ) {
                                if (edge_color_scheme == "monochrome") {
                                    sig_cat <- "default"
                                } else {
                                    lower <- quantiles[pname, "2.5%"]
                                    upper <- quantiles[pname, "97.5%"]

                                    if (sign(lower) == sign(upper)) {
                                        # Significant
                                        if (
                                            edge_color_scheme == "directional"
                                        ) {
                                            if (val > 0) {
                                                sig_cat <- "pos"
                                            } else {
                                                sig_cat <- "neg"
                                            }
                                        } else {
                                            # Binary (sig vs ns)
                                            sig_cat <- "sig" # Maps to black
                                        }
                                    } else {
                                        # Non-significant
                                        sig_cat <- "ns"
                                    }
                                }
                            } else {
                                sig_cat <- "default"
                            }
                        }

                        if (!is.na(val)) {
                            dag_data$val[idx] <- val
                            dag_data$weight_abs[idx] <- abs(val)
                            dag_data$edge_label[idx] <- round(val, 2)
                            # Store category
                            dag_data$significant[idx] <- sig_cat
                        }
                    }
                }
            }
        }

        # Fill defaults for plotting if missing
        dag_data$weight_abs[
            is.na(dag_data$weight_abs) & !is.na(dag_data$to)
        ] <- 1.0
        # Default edge type to -> if not found (robustness)
        dag_data$edge_type[
            is.na(dag_data$edge_type) & !is.na(dag_data$to)
        ] <- "->"
        # Default significance to "default" (black)
        if (!"significant" %in% names(dag_data)) {
            dag_data$significant <- "default"
        }
        dag_data$significant[is.na(dag_data$significant)] <- "default"

        # Ensure factor levels
        dag_data$significant <- factor(
            dag_data$significant,
            levels = c("pos", "neg", "sig", "ns", "default")
        )

        dag_data$final_node_size <- current_node_size
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
            stroke = node_stroke,
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
        # Deduplicate bidirected edges to prevent double plotting
        # Create a sorted ID for bidirected edges
        combined_dag_data <- combined_dag_data %>%
            dplyr::mutate(
                edge_id = dplyr::case_when(
                    edge_type == "<->" ~ paste(
                        pmin(name, to),
                        pmax(name, to),
                        sep = "_"
                    ),
                    TRUE ~ paste(name, to, sep = "_")
                )
            ) %>%
            dplyr::distinct(edge_id, edge_type, .keep_all = TRUE) %>%
            dplyr::select(-edge_id)

        # 1. Directed Edges (Solid)
        p <- p +
            ggdag::geom_dag_edges_link(
                data = function(x) dplyr::filter(x, edge_type == "->"),
                mapping = ggplot2::aes(
                    edge_width = weight_abs,
                    edge_colour = significant, # Map color to significance category
                    label = edge_label
                ),
                angle_calc = "along",
                label_dodge = ggplot2::unit(3, "mm")
            )

        # 2. Bidirected Edges (Dotted, Curved, Double Arrow)
        p <- p +
            ggdag::geom_dag_edges_arc(
                data = function(x) dplyr::filter(x, edge_type == "<->"),
                mapping = ggplot2::aes(
                    edge_width = weight_abs,
                    edge_colour = significant, # Map color to significance category
                    label = edge_label
                ),
                curvature = 0.3,
                edge_linetype = "dotted", # Changed to dotted for clarity
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
            ggraph::scale_edge_colour_manual(
                values = c(
                    "pos" = "firebrick", # Red for positive sig
                    "neg" = "dodgerblue", # Blue for negative sig
                    "sig" = "black", # Black for significant (no direction)
                    "ns" = "grey70", # Grey for non-sig
                    "default" = "black" # Black for structural/none
                ),
                guide = "none"
            )
    } else {
        p <- p + ggdag::geom_dag_edges()
    }

    # Calculate uniform node size
    # Handle case where combined_dag_data might be missing final_node_size or empty
    if (is.null(combined_dag_data$final_node_size)) {
        uniform_node_size <- node_size
    } else {
        uniform_node_size <- max(
            combined_dag_data$final_node_size,
            na.rm = TRUE
        )
        # Fallback if max is -Inf
        if (!is.finite(uniform_node_size)) uniform_node_size <- node_size
    }

    # Nodes Layer
    p <- p +
        ggdag::geom_dag_node(
            size = uniform_node_size,
            shape = 21, # Circle with border
            color = node_color,
            fill = node_fill,
            stroke = node_stroke
        ) +
        ggdag::geom_dag_text(
            ggplot2::aes(label = label_display),
            size = text_size,
            colour = "black" # Text color
        )

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
