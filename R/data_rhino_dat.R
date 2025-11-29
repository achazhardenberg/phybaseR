#' Rhinograd life-history data
#'
#' A dataset containing simulated life-history trait data for 100 Rhinograd species.
#' This data is used to demonstrate phylogenetic Bayesian structural equation
#' modeling with the phybaseR package.
#'
#' @format A data frame with 100 rows (one per species) and 6 columns:
#' \describe{
#'   \item{SP}{Character. Species names matching the tip labels in rhino.tree}
#'   \item{BM}{Numeric. Body mass (log-transformed)}
#'   \item{LS}{Numeric. Litter size (log-transformed)}
#'   \item{NL}{Numeric. Number of litters per year (log-transformed)}
#'   \item{DD}{Numeric. Development duration (log-transformed)}
#'   \item{RS}{Numeric. Reproductive seasonality (binary: 0 or 1)}
#' }
#'
#' @details
#' This simulated dataset represents life-history traits for hypothetical
#' Rhinograd species. The data can be used to test causal hypotheses about
#' life-history evolution using phylogenetic structural equation models.
#'
#' Example causal model (Model 8 from the README):
#' \itemize{
#'   \item Body mass (BM) affects litter size (LS)
#'   \item Body mass (BM) and reproductive seasonality (RS) affect number of litters (NL)
#'   \item Number of litters (NL) affects development duration (DD)
#' }
#'
#' @source Simulated data for von Hardenberg and Gonzalez-Voyer (2025)
#'
#' @references
#' von Hardenberg, A. and Gonzalez-Voyer, A. (2025).
#' PhyBaSE: A Bayesian approach to Phylogenetic Structural Equation Models.
#' Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.70044
#'
#' @examples
#' data(rhino.dat)
#' head(rhino.dat)
#' summary(rhino.dat)
#'
"rhino.dat"
