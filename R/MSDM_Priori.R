#' Create spatial predictor variables to reduce overprediction of species distribution models
#'
#' @description This function creates rasters that, together with environmental variables, can be used to construct constrained species distribution models spatial
#' @param records data.frame. A database with geographical coordinates of species presences used to create species distribution models.
#' @param x character. Column name with longitude values.
#' @param y character. Column name with latitude values.
#' @param sp character. Column name with species names.  It would be desirable that the species names are as simple as possible and with no space between the genus and the specific epithet (e.g., Alchornea_glandulosa).
#' Do not use author names, symbols, or accents. For example, substitute names like Senna chacoensis (L.Bravo) H.S.Irwin & Barneby or Erythrina crista-galli L., for Senna_chacoensis and Erythrina_cristagalli. The species names and the raster must be the same.
#' @param method character. A character string indicating which MSDM the method that must be used. The next methods are available: XY, MIN, CML, and KER. Usage method = 'CML'
#' @param rasterlayer raster object. A raster, stack, or brick object used to construct species distribution models. This object will be used to create MSDM variables with the same resolution, extent, and pattern of empty cells that the environmental variables. It is advisable to use a raster of an environmental layer that will be used in the future to create the species distribution models in order not to have problems (e.g., resolution, extent, cells with NA) between environmental and constraining raster.
#' @param dirsave character. A character string indicating the directory where the result must be saved.
#' @return This function saves raster files (with geotiff format) in a folder named with the MSDM method. Such raster/s have to be used together with environmental variables at the moment to construct species distribution models. The XY approach creates a single pair of raster layers that can be used for all species that share the same study region. Otherwise, CML, MIN, and KER create a species-specific raster layer.
#'
#' @details
#' XY (Latlong method). It assumes that spatial structure can partially explain species distribution (Bahn & Mcgill, 2007). Two raster layers will be created, containing the latitude and longitude of pixels in decimal degrees, respectively. These raster layers should be included as covariates with the environmental layers to construct species distribution models.
#'
#'
#' MIN (Nearest neighbor distance). Compiled and adapted from Allouche et al. (2008), this method calculates for each cell the Euclidian geographic distance to the nearest presence point.
#'
#'
#' CML (Cumulative distance method). Compiled and adapted from Allouche et al. (2008), it assumes that pixels closer to presences are likely included in species distributions. A raster layer will be created containing the sum of Euclidian geographic distances from each pixel to all occurrences of a species. Obtained values are normalized to vary from zero to one. This raster layer should be included as covariates with the environmental layers to construct species distribution models.
#'
#'
#' KER (Kernel method). Also compiled and adapted from Allouche et al. (2008), this method, like CML, assumes that pixels located in areas with a higher density of occurrences are likely included in the actual species distribution. Therefore, a raster layer will be created containing the Gaussian values based on the density of occurrences of a species. The standard deviation of the Gaussian distribution was the maximum value in a vector of minimum distances between pairs of occurrences of a species. Gaussian values are normalized to vary from zero to one. This raster layer should be included as covariates with the environmental layers to construct species distribution models.
#'
#'
#' Further methodological and performance information of these methods see Mendes et al. (2020).
#'
#' @references
#' \itemize{
#' \item Mendes, P.; Velazco S.J.E.; Andrade, A.F.A.; De Marco, P. (2020) Dealing with overprediction in species distribution models: how adding distance constraints can improve model accuracy, Ecological Modelling, in press. https://doi.org/10.1016/j.ecolmodel.2020.109180
#' \item Allouche, O.; Steinitz, O.; Rotem, D.; Rosenfeld, A.; Kadmon, R. (2008). Incorporating distance constraints into species distribution models. Journal of Applied Ecology, 45(2), 599-609. doi:10.1111/j.1365-2664.2007.01445.x
#' \item Bahn, V.; Mcgill, B. J. (2007). Can niche-based distribution models outperform spatial interpolation? Global Ecology and Biogeography, 16(6), 733-742. doi:10.1111/j.1466-8238.2007.00331.x
#' }
#'
#'
#' @examples
#' \dontrun{
#' library(MSDM)
#' library(raster)
#'
#' # Raster data and a data.frame with occurrences will be loaded
#' data("sp_sdm")
#' data("occurrences")
#' head(occurrences)
#'
#' # Some changes will be done on sp_sdm object to be used as example
#' plot(sp_sdm)
#' class(sp_sdm)
#' sp_sdm <- sp_sdm[[1]] # a layer of this RasterBrick will be selected
#' class(sp_sdm)
#' plot(sp_sdm)
#'
#'
#' tmdir <- tempdir()
#' tmdir # temporal directory where will be saves raster layers
#'
#' # # XY method----
#' MSDM_Priori(
#'   records = occurrences,
#'   x = "x", y = "y", sp = "sp", method = "XY",
#'   rasterlayer = sp_sdm, dirsave = tmdir
#' )
#'
#' # open directory were raster were saved
#' rdir <- paste(tmdir, "MSDM_XY", sep = "/")
#' rdir
#'
#' # plot results
#' new_var <- list.files(rdir, pattern = ".tif", full.names = TRUE)
#' new_var <- stack(new_var)
#' plot(new_var)
#'
#'
#' # CML method----
#' MSDM_Priori(
#'   records = occurrences,
#'   x = "x", y = "y", sp = "sp", method = "CML",
#'   rasterlayer = sp_sdm, dirsave = tmdir
#' )
#'
#' # open directory were raster were saved
#' rdir <- paste(tmdir, "MSDM_CML", sep = "/")
#' rdir
#' # shell.exec(rdir)
#'
#' # plot results
#' new_var <- list.files(rdir, pattern = ".tif", full.names = TRUE)
#' new_var <- stack(new_var)
#' plot(new_var)
#' # Note that a raster is created for each species
#' }
#'
#' @import raster
#' @import rgdal
#' @importClassesFrom raster RasterStack RasterLayer
#' @importFrom methods as
#' @importFrom flexclust dist2
#' @importFrom sp gridded<-
#' @seealso \code{\link{MSDM_Posteriori}}
#' @export
MSDM_Priori <- function(records,
                        x = NA,
                        y = NA,
                        sp = NA,
                        method = c("XY", "MIN", "CML", "KER"),
                        rasterlayer = NULL,
                        dirsave = NULL) {
  if (any(is.na(c(x, y, sp))) & method != "XY") {
    stop("Complete x, y and sp arguments to use MIN, CML or KER methods")
  }
  if (is.null(rasterlayer)) {
    stop("Complete rasterlayer argument")
  }
  if (is.null(dirsave)) {
    stop("Complete dirsave argument")
  }

  Species <- split(records[, c(x, y)], records[, sp])

  # Method 1: XY----
  if (method == "XY") {
    dir_pri <- "MSDM_XY"
    dir.create(file.path(dirsave, dir_pri))
    dir_pri <- file.path(dirsave, dir_pri)

    # Extract coordinates of raster layer
    rsd <- raster::coordinates(rasterlayer)
    long <- data.frame(rsd, val = NA)
    lat <- data.frame(rsd, val = NA)
    lins <- which(!is.na(rasterlayer[]))

    # Create raster with Longitude
    long[lins, 3] <- long[lins, 1]
    sp::gridded(long) <- ~ x + y
    long <- raster::raster(long)

    # Create raster with Latitude
    lat[lins, 3] <- lat[lins, 2]
    sp::gridded(lat) <- ~ x + y
    lat <- raster::raster(lat)

    # Create LatLong Stack
    result <- raster::stack(long, lat)
    names(result) <- c("Long", "Lat")
    rm(lat, long)

    raster::writeRaster(
      result,
      file.path(dir_pri, names(result)),
      format = "GTiff",
      bylayer = TRUE,
      overwrite = TRUE
    )
  }

  # Method 2- Minimum distance ----
  if (method == "MIN") {
    dir_pri <- "MSDM_MIN"
    dir.create(file.path(dirsave, dir_pri))
    dir_pri <- file.path(dirsave, dir_pri)

    spi <- as(rasterlayer, "SpatialPixels")@coords
    r <- lapply(Species, function(x) {
      x <- raster::rasterize(x, rasterlayer, field = 1)
      x <- as(x, "SpatialPixels")@coords
      return(x)
    })
    distr <- lapply(r, function(x) {
      xx <- flexclust::dist2(spi, x, method = "euclidean", p = 2)
      xx <- apply(xx, 1, min)
      xx <- (xx - min(xx)) / (max(xx) - min(xx))
      return(xx)
    })
    result <- list()
    for (b in 1:length(Species)) {
      msk <- rasterlayer
      msk[!is.na(msk[])] <- distr[[b]]
      result[[b]] <- msk
    }
    result <- raster::stack(result)
    names(result) <- names(Species)
    rm(msk)
    if (raster::nlayers(result) > 1) {
      raster::writeRaster(
        result,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        bylayer = TRUE,
        overwrite = TRUE
      )
    } else {
      raster::writeRaster(envM,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }

  # Method 3: Cummulative distance----
  if (method == "CML") {
    dir_pri <- "MSDM_CML"
    dir.create(file.path(dirsave, dir_pri))
    dir_pri <- file.path(dirsave, dir_pri)

    spi <- as(rasterlayer, "SpatialPixels")@coords
    r <- lapply(Species, function(x) {
      x <- raster::rasterize(x, rasterlayer, field = 1)
      x <- as(x, "SpatialPixels")@coords
      return(x)
    })

    distr <-
      lapply(r, function(x) {
        x <- flexclust::dist2(spi, x, method = "euclidean", p = 2)
        x <- x + 1
        x <- 1 / (1 / x^2)
        x <- apply(x, 1, sum)
        x <- (x - min(x)) / (max(x) - min(x))
        return(x)
      })
    envM <- list()
    for (b in 1:length(Species)) {
      spdist <- rasterlayer
      spdist[!is.na(spdist[])] <- distr[[b]]
      envM[[b]] <- spdist
    }
    envM <- raster::stack(envM)
    names(envM) <- names(Species)
    rm(spdist)
    if (nlayers(envM) > 1) {
      raster::writeRaster(
        envM,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        bylayer = TRUE,
        overwrite = TRUE
      )
    } else {
      raster::writeRaster(envM,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }

  # Method 4: Gaussian Kernel----
  if (method == "KER") {
    dir_pri <- "MSDM_KER"
    dir.create(file.path(dirsave, dir_pri))
    dir_pri <- file.path(dirsave, dir_pri)

    spi <- as(rasterlayer, "SpatialPixels")@coords
    r <- lapply(Species, function(x) {
      x <- raster::rasterize(x, rasterlayer, field = 1)
      x <- as(x, "SpatialPixels")@coords
      return(x)
    })

    distp <-
      lapply(r, function(x) {
        flexclust::dist2(x, x, method = "euclidean", p = 2)
      })
    distp1 <- lapply(distp, function(x) {
      matrix(0, nrow(x), 1)
    })
    result <- list()
    for (b in 1:length(Species)) {
      distsp <- distp[[b]]
      for (c in 1:nrow(distsp)) {
        vec <- distsp[c, ]
        distp1[[b]][c] <- min(vec[vec != min(vec)])
      }
      sd_graus <- max(distp1[[b]])
      distr <- flexclust::dist2(spi, r[[b]], method = "euclidean", p = 2)
      distr2 <- distr
      distr2 <-
        (1 / sqrt(2 * pi * sd_graus) * exp(-1 * (distr / (2 * sd_graus^2))))
      distr2 <- apply(distr2, 1, sum)
      spdist <- rasterlayer
      spdist[!is.na(spdist[])] <- distr2
      result[[b]] <- spdist
    }
    result <- lapply(result, function(x) {
      (x - raster::cellStats(x, min)) / (raster::cellStats(x, max) - raster::cellStats(x, min))
    })
    result <- raster::stack(result)
    names(result) <- names(Species)
    rm(spdist)
    if (nlayers(result) > 1) {
      raster::writeRaster(
        result,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        bylayer = TRUE,
        overwrite = TRUE
      )
    } else {
      raster::writeRaster(
        result,
        file.path(dir_pri, names(Species)),
        format = "GTiff",
        overwrite = TRUE
      )
    }
  }
  cat("results are in: \n", dir_pri, "\n")
}
