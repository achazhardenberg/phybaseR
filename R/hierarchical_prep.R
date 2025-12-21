#' Prepare Hierarchical Data for JAGS
#'
#' Transforms hierarchical dataframes into a flat list of vectors and ID indices suitable for JAGS
#'
#' @param hierarchical_info List containing 'data', 'levels', 'hierarchy', 'link_vars'
#' @param vars_needed Character vector of all variable names needed for the model
#' @return List with 'data_list' (for jags) and 'n_vec' (named vector of sample sizes)
#' @keywords internal
prepare_hierarchical_jags_data <- function(hierarchical_info, vars_needed) {
    data_list <- list()
    n_vec <- list()

    # 1. Extract variables from their native levels
    # Iterate through each level
    for (lvl_name in names(hierarchical_info$data)) {
        df <- hierarchical_info$data[[lvl_name]]
        n_vec[[paste0("N_", lvl_name)]] <- nrow(df)

        # Identify variables that belong to this level
        lvl_vars <- hierarchical_info$levels[[lvl_name]]

        # Filter for only those needed in the model
        vars_in_model <- intersect(lvl_vars, vars_needed)

        for (v in vars_in_model) {
            if (!v %in% names(df)) {
                stop(sprintf(
                    "Variable '%s' expected in level '%s' but not found in data.",
                    v,
                    lvl_name
                ))
            }
            data_list[[v]] <- df[[v]]
        }
    }

    # 2. Generate Indexing Vectors (Linking Child -> Parent)
    # Parse hierarchy (e.g., "Site > Species")
    # "Site" is child (finer), "Species" is parent (coarser) ??
    # Wait, usually hierarchy is written "Coarser > Finer" in my package?
    # Let's check: "Traits > Site_covs > Y"
    # User said: "Traits (50 rows) ... Site_covs (200 rows)"
    # So "Traits" is Top Level (Parent). "Site_covs" is Middle. "Y" is Bottom (Child).
    # Hierarchy string: "Traits > Site_covs > Y"
    # This implies Parent > Child.

    hierarchy_str <- hierarchical_info$hierarchy
    levels_ordered <- trimws(strsplit(hierarchy_str, ">")[[1]])

    # Iterate strictly through adjacent pairs in the hierarchy
    # We need to link Level[k+1] (Child) to Level[k] (Parent)
    if (length(levels_ordered) > 1) {
        for (k in 1:(length(levels_ordered) - 1)) {
            parent_lvl <- levels_ordered[k]
            child_lvl <- levels_ordered[k + 1]

            parent_df <- hierarchical_info$data[[parent_lvl]]
            child_df <- hierarchical_info$data[[child_lvl]]

            # Determine link variable
            # We assume the link variable is the ID column of the PARENT
            # It must exist in both df
            common_cols <- intersect(names(parent_df), names(child_df))

            # Prefer 'link_vars' if specified
            if (!is.null(hierarchical_info$link_vars)) {
                link_col <- intersect(common_cols, hierarchical_info$link_vars)
                if (length(link_col) > 0) {
                    link_col <- link_col[1]
                } else {
                    link_col <- NULL
                }
            } else {
                link_col <- if (length(common_cols) > 0) {
                    common_cols[1]
                } else {
                    NULL
                }
            }

            if (is.null(link_col)) {
                stop(sprintf(
                    "Cannot link levels '%s' and '%s'. No common column found.",
                    parent_lvl,
                    child_lvl
                ))
            }

            # Create Integer Index
            # Ensure factors match
            parent_vals <- as.character(parent_df[[link_col]])
            child_vals <- as.character(child_df[[link_col]])

            # Verify integrity
            missing_links <- setdiff(child_vals, parent_vals)
            if (length(missing_links) > 0) {
                stop(sprintf(
                    "Level '%s' contains values in '%s' not found in parent level '%s'.",
                    child_lvl,
                    link_col,
                    parent_lvl
                ))
            }

            idx_vec <- match(child_vals, parent_vals)

            # Naming convention: {Parent}_ID_in_{Child} ?? or just {Parent}_ID?
            # If we are in loop of Child, we access Parent[Parent_ID[i]]
            # But Child might not be the 'main' loop?
            # Actually, simple name: "{Parent}_ID" (vector of length Child)
            # Wait, if there are 3 levels A > B > C.
            # B needs A_ID (length N_B).
            # C needs B_ID (length N_C).
            # We name it specifically for the level it resides in?
            # No, JAGS models usually assume implicit context.
            # But here we have explicit loops.
            # Loop C uses: B_var[ B_index_in_C[i] ]
            # Loop B uses: A_var[ A_index_in_B[i] ]

            idx_name <- paste0(parent_lvl, "_idx_", child_lvl) # Explicit: Parent index vector residing in Child

            data_list[[idx_name]] <- idx_vec

            # Store mapping for model generator
            # We need to know: To access Parent from Child, use vector 'idx_name'
        }
    }

    return(list(data_list = data_list, n_vec = n_vec))
}
