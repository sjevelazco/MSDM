#' Virtual species occurrences.
#'
#' A dataset containing occurrences of 5 virtual species
#'
#' @format A data frame with 250 rows and 3 variables:
#' \describe{
#'   \item{sp}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   ...
#' }
"occurrences"

#' Virtual species absences
#'
#' A dataset containing absences of 5 virtual species
#'
#' @format A data frame with 250 rows and 3 variables:
#' \describe{
#'   \item{sp}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   ...
#' }
"absences"

#' A raster layer used to as expample in MSDM_Priori methods
#'
#' A dataset containing absences of 5 virtual species
#'
#' @format A raster of Americas continent from southern USA to southern South America:
#' @examples
#'
#' require(raster)
#' data("rlayer")
#' plot(rlayer)
"rlayer"



