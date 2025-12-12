#' Validate Hierarchical Data Structure
#'
#' @param data List of data.frames at different hierarchical levels
#' @param levels List mapping variable names to level names
#' @param hierarchy Character string specifying nesting (e.g., "site_year > individual")
#' @param link_vars Character vector of variables that link levels
#' @keywords internal
validate_hierarchical_data <- function(data, levels, hierarchy, link_vars) {
    # Check data is a named list of data.frames
    if (!is.list(data) || is.data.frame(data)) {
        stop(
            "For hierarchical data, 'data' must be a named list of data.frames"
        )
    }

    if (is.null(names(data)) || any(names(data) == "")) {
        stop("All elements of hierarchical 'data' must be named")
    }

    if (!all(sapply(data, is.data.frame))) {
        stop("All elements of hierarchical 'data' must be data.frames")
    }

    # Check levels is provided and valid
    if (is.null(levels)) {
        stop("'levels' argument required when using hierarchical data")
    }

    if (!is.list(levels) || is.null(names(levels))) {
        stop("'levels' must be a named list")
    }

    # Check that level names in 'levels' match dataset names in 'data'
    if (!all(names(levels) %in% names(data))) {
        missing <- setdiff(names(levels), names(data))
        stop("Level names not found in data: ", paste(missing, collapse = ", "))
    }

    # Check that all variables in levels exist in corresponding datasets
    for (level_name in names(levels)) {
        vars <- levels[[level_name]]
        dataset <- data[[level_name]]
        missing_vars <- setdiff(vars, colnames(dataset))

        if (length(missing_vars) > 0) {
            stop(
                "Variables not found in ",
                level_name,
                " dataset: ",
                paste(missing_vars, collapse = ", ")
            )
        }
    }

    # Check for variable overlap (variables should not appear in multiple levels)
    all_vars <- unlist(levels)
    if (any(duplicated(all_vars))) {
        dups <- all_vars[duplicated(all_vars)]
        stop(
            "Variables appear in multiple levels (not allowed): ",
            paste(unique(dups), collapse = ", ")
        )
    }

    # Validate link_vars if provided
    if (!is.null(link_vars)) {
        for (dataset in data) {
            missing_links <- setdiff(link_vars, colnames(dataset))
            if (length(missing_links) > 0) {
                stop(
                    "Link variables not found in all datasets: ",
                    paste(missing_links, collapse = ", ")
                )
            }
        }
    }

    # Validate hierarchy string if provided
    if (!is.null(hierarchy)) {
        if (!is.character(hierarchy) || length(hierarchy) != 1) {
            stop("'hierarchy' must be a single character string")
        }

        # Parse hierarchy (e.g., "site_year > individual")
        hierarchy_levels <- strsplit(hierarchy, "\\s*>\\s*")[[1]]

        # Check all hierarchy levels exist in data
        missing_h <- setdiff(hierarchy_levels, names(data))
        if (length(missing_h) > 0) {
            stop(
                "Hierarchy levels not found in data: ",
                paste(missing_h, collapse = ", ")
            )
        }
    }

    invisible(TRUE)
}


#' Infer Variable Level
#'
#' Determine which hierarchical level a variable belongs to
#'
#' @param var Character, variable name
#' @param levels List mapping variable names to level names
#' @return Character, level name
#' @keywords internal
infer_variable_level <- function(var, levels) {
    for (level_name in names(levels)) {
        if (var %in% levels[[level_name]]) {
            return(level_name)
        }
    }

    # Variable not found in any level
    stop("Variable '", var, "' not found in any hierarchical level")
}


#' Get Data for Variables
#'
#' Determine the finest grain level needed for a set of variables
#' and return the appropriate dataset (with joining if needed)
#'
#' @param variables Character vector of variable names
#' @param data List of data.frames at different levels
#' @param levels List mapping variable names to level names
#' @param hierarchy Character string specifying nesting
#' @param link_vars Character vector of linking variables
#' @return data.frame
#' @keywords internal
get_data_for_variables <- function(
    variables,
    data,
    levels,
    hierarchy,
    link_vars
) {
    # Determine which level each variable belongs to
    var_levels <- sapply(variables, function(v) {
        infer_variable_level(v, levels)
    })

    # Get unique levels needed
    needed_levels <- unique(var_levels)

    # If only one level, return that dataset
    if (length(needed_levels) == 1) {
        return(data[[needed_levels]])
    }

    # Multiple levels - need to determine which is finest grain and join
    # Parse hierarchy to get ordering
    hierarchy_order <- strsplit(hierarchy, "\\s*>\\s*")[[1]]

    # Find the finest grain level among those needed
    # (furthest right in hierarchy string)
    finest_idx <- max(match(needed_levels, hierarchy_order))
    finest_level <- hierarchy_order[finest_idx]

    # Start with finest grain dataset
    result <- data[[finest_level]]

    # Join data from coarser levels
    coarser_levels <- setdiff(needed_levels, finest_level)

    for (coarser_level in coarser_levels) {
        # Get variables from this level that we need
        vars_from_this_level <- names(var_levels)[var_levels == coarser_level]

        # Select those variables plus link vars from coarser dataset
        coarser_data <- data[[coarser_level]][,
            c(link_vars, vars_from_this_level),
            drop = FALSE
        ]

        # Helper: Prune colliding variables from result that are claimed by coarser level
        # This ensures 'WINt' from enviro overrides 'WINt' in individual (if present but unmapped)
        # and prevents them from becoming unintended join keys.
        vars_to_add <- vars_from_this_level
        potential_collisions <- intersect(vars_to_add, names(result))

        # Don't prune if it's explicitly a link var (we want to use it for joining)
        if (!is.null(link_vars)) {
            potential_collisions <- setdiff(potential_collisions, link_vars)
        }

        if (length(potential_collisions) > 0) {
            # Drop colliding columns from result to prefer the coarser (source) version
            result[potential_collisions] <- NULL
        }

        # Join to result
        # Handle by=NULL case explicitly to avoid Cross Join mangling (.x, .y)
        join_by <- link_vars
        if (is.null(join_by)) {
            join_by <- intersect(names(result), names(coarser_data))
        }

        result <- merge(result, coarser_data, by = join_by, all.x = TRUE)
    }

    return(result)
}


#' Parse Hierarchy from Random Effects
#'
#' Extract hierarchical nesting structure from random effects formula
#'
#' @param random Formula specifying random effects
#' @param data List of data.frames at different hierarchical levels (optional, currently unused)
#' @return Character string hierarchy (e.g., "site_year > individual") or NULL
#' @keywords internal
parse_hierarchy_from_random <- function(random, data = NULL) {
    if (is.null(random)) {
        return(NULL)
    }

    # Convert to character
    random_str <- deparse(random)

    # Look for explicit nesting syntax (1|A/B)
    # After expansion, this will be (1|A) + (1|A:B)
    # We want to extract A > A:B pattern
    nesting_pattern <- "\\(1\\s*\\|\\s*([^/)]+)/([^)]+)\\)"

    if (grepl(nesting_pattern, random_str)) {
        matches <- regmatches(
            random_str,
            regexec(nesting_pattern, random_str)
        )[[1]]

        if (length(matches) >= 3) {
            # matches[2] is outer level, matches[3] is inner level
            outer <- trimws(matches[2])
            inner <- trimws(matches[3])

            # Clean up any interaction terms (site:year -> site_year)
            outer <- gsub(":", "_", outer)
            inner <- gsub(":", "_", inner)

            return(paste(outer, ">", inner))
        }
    }

    return(NULL)
}
