#' Storks and Babies Data
#'
#' Data form Matthews (2000) used to demonstrate spurious correlation.
#'
#' @format A data frame with 17 rows and 5 variables:
#' \describe{
#'   \item{Country}{Name of the country}
#'   \item{Area}{Area of the country in km^2}
#'   \item{Storks}{Number of stork pairs}
#'   \item{Humans}{Human population size in millions}
#'   \item{Birth}{Number of births in thousands per year}
#' }
#'
#' @source Matthews, R. (2000). Storks deliver babies (p= 0.008).
#' Teaching Statistics, 22(2), 36-38.
#'
#' @examples
#' data(storks.dat)
#' plot(Birth ~ Storks, data = storks.dat)
"storks.dat"
