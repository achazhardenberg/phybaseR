#' Rhinograd life-history data
#'
#' A dataset containing simulated life-history trait data for 100 Rhinograd species.
#' This data is used to demonstrate phylogenetic Bayesian structural equation
#' modeling with the becauseR package.
#'
#' @format A data frame with 100 rows (one per species) and 6 columns:
#' \describe{
#'   \item{SP}{Character. Species names matching the tip labels in rhino.tree}
#'   \item{BM}{Numeric. Body mass (standardised)}
#'   \item{LS}{Numeric. Litter size (standardised)}
#'   \item{NL}{Numeric. Nose length (standardised)}
#'   \item{DD}{Numeric. Dispersal distance (standardised)}
#'   \item{RS}{Numeric. Range size (standardised)}
#' }
#'
#' @details
#' This simulated dataset represents life-history traits for hypothetical
#' Rhinograd species. The data can be used to test causal hypotheses about
#' life-history evolution using phylogenetic structural equation models.
#'
#' Example causal model (Model 8 from Gonzalez-Voyer and von Hardenberg (2014)):
#' \itemize{
#'   \item Body mass (BM) affects litter size (LS)
#'   \item Body mass (BM) and range size (RS) affect nose length (NL)
#'   \item Nose length (NL) affects dispersal distance (DD)
#' }
#'
#' @source Simulated data for Gonzalez-Voyer & von Hardenberg (2014)
#'
#' @references
#' Gonzalez-Voyer, A., & von Hardenberg, A. (2014). An introduction to
#' phylogenetic path analysis. In: Modern phylogenetic comparative methods
#' and their application in evolutionary biology (pp. 201–229). Springer.
#'
#' @examples
#' data(rhino.dat)
#' head(rhino.dat)
#' summary(rhino.dat)
#'
"rhino.dat"
