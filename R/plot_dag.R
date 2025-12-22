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
#' @param layout The layout algorithm to use (default "sugiyama"). See \code{\link[ggdag]{ggdag}}.
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
#'
plot_dag <- function(
    x,
    layout = "sugiyama",
    latent = NULL,
    node_size = 14,
    node_color = "black",
    node_fill = "white",
    node_stroke = 1.5,
    text_size = 4,
    edge_width_range = c(0.5, 2),
    edge_color_scheme = c("directional", "binary", "monochrome"),
    show_coefficients = TRUE,
    coords = NULL,
    family = NULL
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
        current_family <- family

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
            # Also get family if not provided
            if (is.null(current_family)) {
                current_family <- obj$input$family
            }
        } else {
            # List of formulas
            eqs <- obj
        }

        # Handle occupancy expansion for visualization
        occ_vars <- c()
        if (!is.null(current_family)) {
            occ_vars <- names(current_family)[current_family == "occupancy"]
            if (length(occ_vars) > 0) {
                # Add p_ and psi_ to latent set so they become circles
                latent_to_add <- unlist(lapply(occ_vars, function(v) {
                    c(paste0("p_", v), paste0("psi_", v))
                }))
                current_latent <- unique(c(current_latent, latent_to_add))
            }
        }

        # 2. Convert to dagitty syntax
        induced_cors <- NULL
        if (inherits(obj, "because")) {
            induced_cors <- obj$induced_correlations %||%
                obj$input$induced_correlations
        }

        dag_str <- equations_to_dag_string(
            eqs,
            induced_cors,
            family = current_family
        )
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

        # Identify occupancy nodes and assign to groups for box drawing
        dag_data$occ_species <- NA_character_
        if (length(occ_vars) > 0) {
            for (v in occ_vars) {
                # Pattern match p_Species and psi_Species
                dag_data$occ_species[grepl(
                    paste0("^(p_|psi_)", v, "$"),
                    dag_data$name
                )] <- v
                # Also include the observation node if it exists
                dag_data$occ_species[dag_data$name == v] <- v
            }
        }

        # Label Processing: Wrap text (replace _ with \n)
        dag_data$label_display <- gsub("_", "\n", dag_data$name)
        # Override for occupancy latent nodes to just 'p' and 'psi'
        if (length(occ_vars) > 0) {
            dag_data$label_display <- ifelse(
                grepl("^p_", dag_data$name) & !is.na(dag_data$occ_species),
                "p",
                dag_data$label_display
            )
            dag_data$label_display <- ifelse(
                grepl("^psi_", dag_data$name) & !is.na(dag_data$occ_species),
                "psi",
                dag_data$label_display
            )
        }

        # Determine Node Size based on longest label (Heuristic)
        # Assuming circular node, diameter needs to cover the rectangular text block.
        # We find the max characters in any single line of the wrapped labels.
        max_chars <- max(nchar(unlist(strsplit(dag_data$label_display, "\n"))))
        # Heuristic: base size + scaling factor * chars * text_size
        # Typically node_size=14 covers ~3 chars at text_size=4.
        # Reduced multiplier to prevent excessive size in compact/Rmd plots.
        calc_size <- 8 + (max_chars * 1.8 * (text_size / 4))

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

    # Calculate Bounding Boxes for Occupancy Groups
    occupancy_boxes <- NULL
    if (
        !is.null(combined_dag_data) &&
            any(!is.na(combined_dag_data$occ_species))
    ) {
        # Check if we have coordinates
        if (any(!is.na(combined_dag_data$x))) {
            occupancy_boxes <- combined_dag_data |>
                dplyr::filter(
                    !is.na(occ_species) & grepl("^(p_|psi_)", name)
                ) |>
                dplyr::group_by(model_label, occ_species) |>
                dplyr::summarize(
                    xmin = min(x) - 0.25,
                    xmax = max(x) + 0.25,
                    ymin = min(y) - 0.25,
                    ymax = max(y) + 0.25,
                    .groups = "drop"
                )
        }
    }

    # Ensure type is a factor
    combined_dag_data$type <- factor(
        combined_dag_data$type,
        levels = c("Observed", "Latent")
    )

    # Calculate uniform node size for plotting and caps
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

    # Define Caps based on node size
    # Reduced multiplier from 1.5 to 1.15 to decrease white space/margins.
    # This prevents edges from disappearing in compact plots.
    cap_size <- ggraph::circle(uniform_node_size * 1.15, "pt")

    p <- ggplot2::ggplot(
        combined_dag_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    )

    # 0. Shadowed Boxes for Occupancy (Background Layer)
    if (!is.null(occupancy_boxes)) {
        p <- p +
            ggplot2::geom_rect(
                data = occupancy_boxes,
                ggplot2::aes(
                    xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax
                ),
                inherit.aes = FALSE,
                fill = "grey93",
                color = "grey80",
                alpha = 0.5,
                linewidth = 0.4,
                lty = "dashed"
            ) +
            # Species name label for the box
            ggplot2::geom_text(
                data = occupancy_boxes,
                ggplot2::aes(
                    x = (xmin + xmax) / 2,
                    y = ymax + 0.05,
                    label = occ_species
                ),
                inherit.aes = FALSE,
                size = text_size * 1.1,
                fontface = "bold",
                vjust = 0
            )
    }

    p <- p +
        ggdag::theme_dag() +
        # Add margins back (ggdag removes them)
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10, "mm")) +
        # prevent clipping of large nodes
        ggplot2::coord_cartesian(clip = "off") +
        # Expand axes to give space for nodes
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2))

    if (length(unique(combined_dag_data$model_label)) > 1) {
        p <- p + ggplot2::facet_wrap(~model_label)
    }

    # Edges Layer: Split into Directed and Bidirected
    if ("edge_type" %in% names(combined_dag_data)) {
        # Normalize bidirected edges to curve AWAY from the graph center.
        # Heuristic:
        # 1. Calculate Graph Centroid.
        # 2. For each edge, check if the "Right" side (curvature > 0) points towards or away from centroid.
        # 3. Swap direction if it points towards center.

        is_bidirected <- combined_dag_data$edge_type == "<->" &
            !is.na(combined_dag_data$edge_type)
        if (any(is_bidirected)) {
            # Calculate Centroid of the graph layout
            # Use unique node positions to avoid weighting by edge count
            unique_nodes <- unique(combined_dag_data[, c("name", "x", "y")])
            centroid_x <- mean(unique_nodes$x, na.rm = TRUE)
            centroid_y <- mean(unique_nodes$y, na.rm = TRUE)

            # Vectorized check
            # Current Vector P1 -> P2
            v_x <- combined_dag_data$xend - combined_dag_data$x
            v_y <- combined_dag_data$yend - combined_dag_data$y

            # Midpoint
            m_x <- (combined_dag_data$x + combined_dag_data$xend) / 2
            m_y <- (combined_dag_data$y + combined_dag_data$yend) / 2

            # Normal Vector pointing "Right" (Curvature direction)
            # R = (dy, -dx)
            r_x <- v_y
            r_y <- -v_x

            # Test Point: Midpoint + Normal
            t_x <- m_x + r_x
            t_y <- m_y + r_y

            # Distances to Centroid
            # dist_current: distance from (Midpoint + Right) to Centroid
            dist_current <- (t_x - centroid_x)^2 + (t_y - centroid_y)^2

            # dist_swapped: distance from (Midpoint - Right) to Centroid
            # (Midpoint - Right) coincides with the curve apex if we swapped direction
            t_swapped_x <- m_x - r_x
            t_swapped_y <- m_y - r_y
            dist_swapped <- (t_swapped_x - centroid_x)^2 +
                (t_swapped_y - centroid_y)^2

            # If swapped distance is GREATER (farther from center), we should swap
            # to make the curve point that way.
            needs_swap <- is_bidirected & (dist_swapped > dist_current)

            if (any(needs_swap, na.rm = TRUE)) {
                swap_idx <- which(needs_swap)

                # Temporary storage to enable swapping
                tmp_name <- combined_dag_data$name[swap_idx]
                tmp_x <- combined_dag_data$x[swap_idx]
                tmp_y <- combined_dag_data$y[swap_idx]

                combined_dag_data$name[swap_idx] <- combined_dag_data$to[
                    swap_idx
                ]
                combined_dag_data$x[swap_idx] <- combined_dag_data$xend[
                    swap_idx
                ]
                combined_dag_data$y[swap_idx] <- combined_dag_data$yend[
                    swap_idx
                ]

                combined_dag_data$to[swap_idx] <- tmp_name
                combined_dag_data$xend[swap_idx] <- tmp_x
                combined_dag_data$yend[swap_idx] <- tmp_y
            }
        }

        # Deduplicate bidirected edges to prevent double plotting
        combined_dag_data <- dplyr::mutate(
            combined_dag_data,
            edge_id = dplyr::case_when(
                edge_type == "<->" ~ paste(
                    pmin(name, to), # Use sorting for ID just to be robust/safe
                    pmax(name, to),
                    sep = "_"
                ),
                TRUE ~ paste(name, to, sep = "_")
            )
        )
        combined_dag_data <- dplyr::distinct(
            combined_dag_data,
            edge_id,
            edge_type,
            .keep_all = TRUE
        )
        combined_dag_data <- dplyr::select(combined_dag_data, -edge_id)

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
                label_dodge = ggplot2::unit(3, "mm"),
                start_cap = cap_size,
                end_cap = cap_size
            )

        # 2. Bidirected Edges (Solid Grey, Curved, Double Arrow)
        p <- p +
            ggdag::geom_dag_edges_arc(
                data = function(x) dplyr::filter(x, edge_type == "<->"),
                mapping = ggplot2::aes(
                    label = edge_label # Keep label
                ),
                edge_width = 0.3, # Constant thin width
                curvature = 0.5, # Sharper curve to avoid corners/overlap
                edge_colour = "grey60", # Fixed grey color
                edge_linetype = "solid", # Continuous line
                angle_calc = "along",
                label_dodge = ggplot2::unit(3, "mm"),
                arrow = ggplot2::arrow(
                    length = ggplot2::unit(2.5, "mm"),
                    type = "closed",
                    ends = "both"
                ),
                start_cap = cap_size,
                end_cap = cap_size
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
        p <- p +
            ggdag::geom_dag_edges(
                start_cap = cap_size,
                end_cap = cap_size
            )
    }

    # Nodes Layer (Final, Single Layer)
    p <- p +
        ggplot2::geom_point(
            ggplot2::aes(shape = type), # Map shape to observed/latent
            size = uniform_node_size,
            color = node_color,
            fill = node_fill,
            stroke = node_stroke
        ) +
        ggdag::geom_dag_text(
            ggplot2::aes(label = label_display),
            size = text_size,
            colour = "black" # Text color
        ) +
        # Map shapes: Square (22) for Observed, Circle (21) for Latent
        ggplot2::scale_shape_manual(
            values = c(Observed = 22, Latent = 21),
            guide = "none" # Remove legend
        )

    return(p)
}


#' Convert Equations List to DAGitty String
#' @param equations List of formulas
#' @param induced_cors List of character vectors (pairs) for bidirected edges
#' @noRd
equations_to_dag_string <- function(
    equations,
    induced_cors = NULL,
    family = NULL
) {
    edges <- c()
    occ_vars <- c()
    if (!is.null(family)) {
        occ_vars <- names(family)[family == "occupancy"]
    }

    for (eq in equations) {
        resp <- all.vars(eq)[1]
        predictors <- all.vars(eq)[-1]

        # If response is occupancy, it's actually psi_response
        # unless it's already p_response
        actual_resp <- resp
        if (resp %in% occ_vars) {
            actual_resp <- paste0("psi_", resp)
        }

        for (pred in predictors) {
            actual_pred <- pred
            # Handle p_ and psi_ aliases
            if (pred %in% occ_vars) {
                actual_pred <- paste0("psi_", pred)
            }
            edges <- c(edges, paste(actual_resp, "<-", actual_pred))
        }
    }

    # Ensure p_Species exists for every occupancy variable
    # even if it's just an intercept-only 1
    for (v in occ_vars) {
        p_name <- paste0("p_", v)
        # We don't strictly need internal edges if they aren't in equations,
        # but for the visual "box", we need both nodes to exist.
        # Adding a self-loop or just an edge between p and psi if they are both in dag?
        # Actually, adding any edge involving p_name ensures it's in the DAG.
    }

    if (!is.null(induced_cors)) {
        for (pair in induced_cors) {
            if (length(pair) == 2) {
                p1 <- pair[1]
                p2 <- pair[2]
                if (p1 %in% occ_vars) {
                    p1 <- paste0("psi_", p1)
                }
                if (p2 %in% occ_vars) {
                    p2 <- paste0("psi_", p2)
                }
                edges <- c(edges, paste(p1, "<->", p2))
            }
        }
    }

    return(paste("dag {", paste(unique(edges), collapse = "; "), "}"))
}
