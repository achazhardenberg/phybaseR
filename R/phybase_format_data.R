#' Format Data for PhyBaSE Analysis
#'
#' Converts data from long format (one row per observation) to the list format
#' required by \code{\link{phybase_run}}.
#'
#' @param data A data.frame in long format with one row per observation.
#' @param species_col Name of the column containing species identifiers (default: "SP").
#' @param tree A phylogenetic tree (class \code{phylo}). Required to determine species order.
#'
#' @return A named list where each element is either:
#'   \itemize{
#'     \item A numeric vector (if all species have exactly 1 observation)
#'     \item A numeric matrix with species in rows and replicates in columns
#'   }
#'   Species are ordered to match \code{tree$tip.label}.
#'
#' @details
#' This function handles:
#' \itemize{
#'   \item Different numbers of replicates per species (creates rectangular matrix with NA padding)
#'   \item Missing values (NA)
#'   \item Automatic alignment with phylogenetic tree tip labels
#' }
#'
#' When species have different numbers of replicates, the function creates a matrix
#' with dimensions (number of species) x (maximum number of replicates).
#' Species with fewer replicates are padded with NA values.
#'
#' Species in the tree but not in the data will have all NA values.
#' Species in the data but not in the tree will be excluded with a warning.
#'
#' @examples
#' \dontrun{
#' # Example data in long format
#' data_long <- data.frame(
#'   SP = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3"),
#'   BM = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8),
#'   NL = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)
#' )
#'
#' tree <- ape::read.tree(text = "(sp1:1,sp2:1,sp3:1);")
#' data_list <- phybase_format_data(data_long, species_col = "SP", tree = tree)
#'
#' # Use with phybase_run
#' fit <- phybase_run(data = data_list, tree = tree, equations = list(NL ~ BM))
#' }
#'
#' @export
phybase_format_data <- function(data, species_col = "SP", tree) {
    # Validate inputs
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame")
    }

    if (!species_col %in% names(data)) {
        stop(sprintf("Species column '%s' not found in data", species_col))
    }

    if (missing(tree) || is.null(tree)) {
        stop("'tree' is required to determine species order")
    }

    if (!inherits(tree, "phylo")) {
        stop("'tree' must be a phylo object")
    }

    # Get trait columns (everything except species column)
    trait_cols <- setdiff(names(data), species_col)

    if (length(trait_cols) == 0) {
        stop("No trait columns found in data")
    }

    # Detect and expand categorical variables (factor or character)
    categorical_vars <- list()
    for (col in trait_cols) {
        if (is.factor(data[[col]]) || is.character(data[[col]])) {
            # Get unique levels (excluding NA)
            levels <- sort(unique(data[[col]][!is.na(data[[col]])]))

            if (length(levels) < 2) {
                warning(sprintf(
                    "Categorical variable '%s' has fewer than 2 levels, skipping dummy coding",
                    col
                ))
                next
            }

            # Store categorical info
            categorical_vars[[col]] <- list(
                levels = levels,
                reference = levels[1],
                dummies = paste0(col, "_", levels[-1])
            )

            # Create dummy variables (K-1 for K levels)
            for (i in 2:length(levels)) {
                dummy_name <- paste0(col, "_", levels[i])
                data[[dummy_name]] <- as.numeric(data[[col]] == levels[i])
            }

            message(sprintf(
                "Categorical variable '%s' expanded to %d dummy variable(s) | Reference: '%s'",
                col,
                length(levels) - 1,
                levels[1]
            ))
        }
    }

    # Update trait_cols: remove categoricals, add dummies
    if (length(categorical_vars) > 0) {
        original_categoricals <- names(categorical_vars)
        all_dummies <- unlist(lapply(categorical_vars, function(x) x$dummies))
        trait_cols <- c(setdiff(trait_cols, original_categoricals), all_dummies)
    }

    # Check for species alignment
    data_species <- unique(data[[species_col]])
    tree_species <- tree$tip.label
    n_species <- length(tree_species)

    missing_in_tree <- setdiff(data_species, tree_species)
    missing_in_data <- setdiff(tree_species, data_species)

    if (length(missing_in_tree) > 0) {
        warning(sprintf(
            "Species in data but not in tree (will be excluded): %s",
            paste(missing_in_tree, collapse = ", ")
        ))
        # Filter out species not in tree
        data <- data[data[[species_col]] %in% tree_species, ]
    }

    if (length(missing_in_data) > 0) {
        message(sprintf(
            "Species in tree but not in data (will have NA values): %s",
            paste(missing_in_data, collapse = ", ")
        ))
    }

    # Count observations per species
    obs_counts <- table(data[[species_col]])
    max_reps <- max(obs_counts)

    # Initialize output list
    data_list <- list()

    # Process each trait
    for (trait in trait_cols) {
        # Create matrix: rows = species (in tree order), columns = replicates
        trait_matrix <- matrix(NA, nrow = n_species, ncol = max_reps)
        rownames(trait_matrix) <- tree_species

        # Fill matrix with observations
        for (i in seq_along(tree_species)) {
            sp <- tree_species[i]

            if (sp %in% names(obs_counts)) {
                # Get observations for this species
                sp_values <- data[data[[species_col]] == sp, trait, drop = TRUE]
                n_obs <- length(sp_values)

                # Fill in the observations (rest remain NA)
                trait_matrix[i, 1:n_obs] <- sp_values
            }
            # else: species not in data, row remains all NA
        }

        # If all species have exactly 1 observation, convert to vector
        if (max_reps == 1) {
            trait_matrix <- as.vector(trait_matrix)
            names(trait_matrix) <- tree_species
        }

        data_list[[trait]] <- trait_matrix
    }

    # Store categorical variable info as attribute for reference
    if (length(categorical_vars) > 0) {
        attr(data_list, "categorical_vars") <- categorical_vars
    }

    return(data_list)
}
