#' Rhinograd Simulated Repeated Measures Data
#'
#' A simulated dataset containing repeated measurements of life-history traits for 100 Rhinograd species.
#' This dataset is in "long" format, with 10 observations per species, suitable for testing
#' measurement error models using `variability = "reps"`.
#'
#' @format A data frame with 1000 rows (10 replicates x 100 species) and 6 columns:
#' \describe{
#'   \item{SP}{Character. Species names matching the tip labels in `rhino.tree`. Each species appears 10 times.}
#'   \item{BM}{Numeric. Body mass (repeated measurements)}
#'   \item{LS}{Numeric. Litter size (repeated measurements)}
#'   \item{NL}{Numeric. Nose length (repeated measurements)}
#'   \item{DD}{Numeric. Dispersal distance (repeated measurements)}
#'   \item{RS}{Numeric. Range size (repeated measurements)}
#' }
#'
#' @details
#' This dataset was simulated by adding white noise to the original `rhino.dat` species values.
#' It represents a scenario where each species has 10 independent observations for each trait.
#' This allows using the `variability = "reps"` option in `because()` to explicitly model measurement error.
#'
#' The underlying causal model is the same as for `rhino.dat`.
#'
#' @seealso \code{\link{rhino.dat}}, \code{\link{RhinoMulti_se.dat}}
#'
#' @source Simulated based on Gonzalez-Voyer & von Hardenberg (2014) data.
"RhinoMulti.dat"

#' Rhinograd Means and Standard Errors
#'
#' A simulated dataset containing the mean and standard error of varying life-history traits for 100 Rhinograd species.
#' This dataset is in "wide" format (one row per species), suitable for measurement error models
#' using `variability = "se"`.
#'
#' @format A data frame with 100 rows and 11 columns:
#' \describe{
#'   \item{SP}{Character. Species names matching the tip labels in `rhino.tree`.}
#'   \item{BM}{Numeric. Mean Body mass}
#'   \item{BM_se}{Numeric. Standard error of Body mass (SD / sqrt(N))}
#'   \item{LS}{Numeric. Mean Litter size}
#'   \item{LS_se}{Numeric. Standard error of Litter size}
#'   \item{NL}{Numeric. Mean Nose length}
#'   \item{NL_se}{Numeric. Standard error of Nose length}
#'   \item{DD}{Numeric. Mean Dispersal distance}
#'   \item{DD_se}{Numeric. Standard error of Dispersal distance}
#'   \item{RS}{Numeric. Mean Range size}
#'   \item{RS_se}{Numeric. Standard error of Range size}
#' }
#'
#' @details
#' This dataset is derived from `RhinoMulti.dat` by calculating the mean and standard error
#' for each trait per species (N=10 replicates).
#' It is designed for use with the `variability = "se"` option in `because()`, where precision is modeled as `1/SE^2`.
#'
#' @seealso \code{\link{rhino.dat}}, \code{\link{RhinoMulti.dat}}
#'
#' @source Calculated from `RhinoMulti.dat`.
"RhinoMulti_se.dat"
