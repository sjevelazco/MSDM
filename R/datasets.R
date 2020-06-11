#' Virtual species occurrences.
#'
#' A dataset containing occurrences of five virtual species
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
#' A dataset containing absences of five virtual species
#'
#' @format A data frame with 250 rows and 3 variables:
#' \describe{
#'   \item{sp}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   ...
#' }
"absences"

#' A dataset containing absences of five virtual species
#'
#' @format A raster of Americas from the southern United States to southern South America:
#' @examples
#'
#' require(raster)
#' data("sp_sdm")
#' plot(sp_sdm)
"sp_sdm"
