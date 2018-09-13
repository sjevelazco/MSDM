#' Methods to correct overprediction of species distribution models based on occurrences and suitability patterns of species.
#'
#' @param records data.frame. A database with geographical coordinates of species presences used to create species distribution models.
#' @param absences data.frame. A database with geographical coordinates of species absences used to create species distribution models.
#' @param x character. Column name with longitude values. This name must be the same for presences and absences database.
#' @param y character. Column name with latitude values. This name must be the same for presences and absences database.
#' @param sp character. Column name with species names. Species names must be the same for presences, absences, and raster layer (i.e. species distribution) databases.
#' @param method character. A character string indicating which MSDM method must be used create.
#' @param dirraster character. A character string indicating the directory where are species raster files, i.e. species distribution models.
#' Raster layer must be in geotiff format
#' @param threshold character. Select type of threshold (kappa, spec_sens, no_omission, prevalence, equal_sens_spec, sensitivty)
#' to get binary models (see \code{\link[dismo]{threshold}} help of dismo package for further information about differt thresholds). Default threshold value is "equal_sens_spec", it is the threshold at which sensitivity and specificity are equal.
#' @param buffer numeric. Buffer width in km to be used in BMCP approach.
#' @param dirsave character. A character string indicating the directory where result must be saved.
#' @return This function save raster files (with geotiff format) with continuous and binary species raster files separated in CONr and BINr folders respectively.

#' @details
#'
#'
#'Abbreviation list
#'
#'\itemize{
#'\item SDM: species distribution model
#'\item l: suitability patches that intercept species occurrences
#'\item k: suitability patches that do not intercept species occurrences
#'\item T: threshold distances used to select suitability patches
#'}
#'
#'
#'All the methods used by this function use raster layer with
#' the suitability values of a species distribution models (SDMs). After applied a threshold to binarize SDMs
#'  (see 'threshold' arguments), function create a binarized raster layer with suitable and unsuitable patches.
#'
#'
#' OBR(Occurrences based restriction)-
#' The method assumes that suitable patches intercepting species occurrences (l)
#' are likely a part of species distributions than suitable patches that do not
#' intercept any occurrence (k). Distance from all patches k species occurrences to the closest l
#' patch is calculated, later it is removed k patches that trespass a species-specific
#' distance threshold from SDMs models. This threshold (T) is calculated as the
#' maximum distance in a vector of minimal pairwise distances between occurrences.
#' Whenever a suitable pixel is within a k patch distant from the closest l in less than T,
#' the suitability of the pixel was reduced to zero. We assumed that this simple threshold
#' is a surrogate of the species-specific dispersal ability. If T is low, either the species
#' has been sampled throughout its distribution, or the species is geographically restricted,
#' justifying a narrow inclusion of k patches.
#'
#' PRES (Only occurrences based restriction). This is a more restrictive variant of the OBR method. It only retains those pixels in suitability patches intercepting occurrences (k)
#'
#' LQ (Lower Quantile). This method is similar to the OBR method, except by the
#' procedure to define a distance threshold to withdrawn k patches, which is the
#' lower quartile distance between k patches to the closest l patch. Whenever a suitable
#' pixel is within a k patch, i.e. not within this lower quartile, the suitability of the
#' pixel is reduced to zero. This means that 75% of k patches were withdrawn from the model.
#'
#' MCP (Minimum Convex Polygon). Compiled and adapted from
#' Kremen et al. (2008), this method excludes from SDMs climate suitable
#' pixels that do not intercept a minimum convex polygon,
#' with interior angles smaller than 180, enclosing all occurrences of a species.
#'
#' BMCP (Buffered Minimum Convex Polygon). Compiled and adapted
#' from Kremen et al. (2008), it is alike the MCP except by the inclusion of a
#' buffer zone surrounding minimum convex polygons. When used this method
#' function will ask for a value in km to be used as the buffer with.
#'
#'
#'@references
#'\itemize{
#'\item Kremen, C., Cameron, A., Moilanen, A., Phillips, S. J., Thomas, C. D.,
#'Beentje, H., . Zjhra, M. L. (2008). Aligning Conservation Priorities Across
#'Taxa in Madagascar with High-Resolution Planning Tools. Science, 320(5873),
#'222-226. doi:10.1126/science.1155193
#'}
#'
#'
#' @examples
#'
#' require(MSDM)
#' require(raster)
#' data("occurrences")
#' data("absences")
#'
#' dir_raster <- system.file("extdata", package="MSDM")
#' # List of raster layer with species suitability
#' list.files(dir_raster, full.names = TRUE)
#'
#' # Create a temporary MSDM folder
#' tmdir <- tempdir()
#' tmdir
#' dir.create(file.path(tmdir,"MSDM"))
#' tmdir <- file.path( tmdir,"MSDM")
#' tmdir
#'
#' # MCP method----
#' MSDM_Posteriori(records=occurrences, absences=absences,
#'                 x="x", y="y", sp="sp", method="MCP",
#'                 dirraster = dir_raster, threshold = "spec_sens",
#'                 dirsave = tmdir)
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#' # Categorical models corrected by MCP methods
#' cat_mcp <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_mcp)
#' # Continuous models corrected by MCP methods
#' con_mcp <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_mcp)
#'
#'
#' # BMCP method----
#' MSDM_Posteriori(records=occurrences, absences=absences,
#'                 x="x", y="y", sp="sp", method="BMCP",
#'                 dirsave = tmdir)
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#' # Categorical models corrected by MCP methods
#' cat_bmcp <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_bmcp)
#' # Continuous models corrected by MCP methods
#' con_bmcp <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_bmcp)
#'
#' # LQ method----
#' MSDM_Posteriori(records=occurrences, absences=absences,
#'                 x="x", y="y", sp="sp", method="LQ",
#'                 dirraster = dir_raster, threshold = "spec_sens",
#'                 dirsave = tmdir)
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#'
#' # Categorical models corrected by LQ methods
#' cat_lq <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_lq)
#' # Continuous models corrected by LQ methods
#' con_lq <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_lq)
#'
#'
#' # OBR method----
#' MSDM_Posteriori(records=occurrences, absences=absences,
#'                 x="x", y="y", sp="sp", method="OBR",
#'                 dirraster = dir_raster, threshold = "spec_sens",
#'                 dirsave = tmdir)
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#'
#' # Categorical models corrected by OBR methods
#' cat_obr <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_obr)
#' # Continuous models corrected by OBR methods
#' con_obr <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_obr)
#'
#'
#' # PRES method----
#' MSDM_Posteriori(records=occurrences, absences=absences,
#'                 x="x", y="y", sp="sp", method="PRES",
#'                 dirraster = dir_raster, threshold = "spec_sens",
#'                 dirsave = tmdir)
#'
#' d <- list.dirs(tmdir, recursive = FALSE)
#'
#' # Categorical models corrected by PRES methods
#' cat_pres <- stack(list.files(d[1], full.names = TRUE))
#' plot(cat_pres)
#' # Continuous models corrected by PRES methods
#' con_pres <- stack(list.files(d[2], full.names = TRUE))
#' plot(con_pres)
#'
#'
#' @import raster
#' @import rgdal
#' @importClassesFrom raster RasterStack RasterLayer
#' @importClassesFrom dismo ModelEvaluation ConvexHull CircleHull
#' @importFrom dismo threshold evaluate convHull circles
#' @importFrom flexclust dist2
#' @importFrom sp coordinates gridded<- "coordinates<-"
#'
#' @seealso \code{\link{MSDM_Priori}}
#' @export
MSDM_Posteriori <- function(records,
                            absences,
                            x = NA,
                            y = NA,
                            sp = NA,
                            method = c('OBR', 'PRES', 'LQ', 'MCP', 'BMCP'),
                            dirraster = NULL,
                            threshold = c('kappa',
                                          'spec_sens',
                                          'no_omission',
                                          'prevalence',
                                          'equal_sens_spec',
                                          'sensitivty'),
                            dirsave = NULL) {
  if (any(is.na(c(x, y, sp)))) {
    stop("Complete 'x', 'y' or 'sp' arguments")
  }
  if (is.null(dirraster)) {
    stop("Complete 'dirraster' argument")
  }
  if (is.null(dirsave)) {
    stop("Complete 'dirsave' argument")
  }
  # if((bynarymodels==FALSE & is.null(threshold))){
  #   stop("Complete 'bynarymodels' argument")
  # }
  if (any(
    threshold == c(
      'kappa',
      'spec_sens',
      'no_omission',
      'prevalence',
      'equal_sens_spec',
      'sensitivty'
    )
  ) == FALSE) {
    stop(
      "'threshold' argument can has on of the next values:
      'kappa', 'spec_sens', 'no_omission',
      'prevalence', 'equal_sens_spec' or 'sensitivty'"
    )
  }


  #Create Binary folder
  foldCat <- paste(dirsave, c("BINr"), sep = "/")
  foldCon <- paste(dirsave, c("CONr"), sep = "/")
  dir.create(foldCat)
  dir.create(foldCon)

  # creation of a data.frame wiht presences and absences
  SpData <-
    rbind(records[, c(sp, x, y)], absences[, c(sp, x, y)])
  SpData$pres_abse <-
    c(rep(1, nrow(records)), rep(0, nrow(absences)))

  # Data.frame wiht two columns 1-names of the species
  # 2-the directory of raster of each species
  if (is.null(dirraster) == TRUE) {
    stop("Give a directory in the dirraster argument")
  } else{
    RasterList <- list.files(dirraster, pattern = '.tif$')
    sps <- gsub('.tif$', '', RasterList)
    RasterList <-
      list.files(dirraster, pattern = '.tif$', full.names = TRUE)
    RasterList <- data.frame(sps, RasterList, stringsAsFactors = FALSE)
    colnames(RasterList) <- c("sp", 'RasterList')
  }

  # Vector with species names, and proving if records and raster layer
  # have the same specie names
  SpNames <- as.character(unique(records[, sp]))
  SpNamesR <- RasterList[, "sp"]


  if (any(!SpNames %in% SpNamesR)) {
    message(sum(!SpNames %in% SpNamesR),
            ' species names differ between records and raster files')
    message("Next names were not found in records database: ",
            paste0(SpNamesR[!SpNames %in% SpNamesR], " "))
    message("Next names were not found in raster layers: ",
            paste0(SpNames[!SpNames %in% SpNamesR], " "))
    stop("species names must be the same in records, absences and raster layers")
  }

  #### threshold for BMCP method
  if (method == "BMCP" & is.null(buffer)) {
    stop("If BMCP approach is used a numeric value must by supplied to 'buffer' argument.")
  } else if (method == "BMCP" & is.numeric(buffer)) {
    buffer <- buffer * 1000
  }

  # loop to process each species
  for (s in 1:length(SpNames)) {
    print(paste(s, "from", length(SpNames), ":", SpNames[s]))
    # Read the raster of the species
    Adeq <-
      raster(RasterList[RasterList[, "sp"] == SpNames[s], 'RasterList'])
    # if (is.na(crs(Adeq))) {
    #   crs(Adeq) <-
    #     "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    # }

    # Extract values for one species and calculate the threshold
    singleSpData <- SpData[SpData$sp == SpNames[s], ]
    PredPoint <- extract(Adeq, singleSpData[, c(x, y)])
    PredPoint <-
      data.frame(pres_abse = singleSpData[, 'pres_abse'], PredPoint)
    Eval <- evaluate(PredPoint[PredPoint$pres_abse == 1, 2],
                     PredPoint[PredPoint$pres_abse == 0, 2])
    Thr <- unlist(c(threshold(Eval))[threshold])

    # MCP method----
    if (method == "MCP") {
      hull <-
        convHull(singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- predict(Adeq, hull, mask = TRUE)
      Adeq[(hull[] == 0)] <- 0
      writeRaster(
        Adeq,
        paste(foldCon, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
      Mask <- Adeq > Thr
      writeRaster(
        Mask,
        paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
    }

    # BMCP method-----
    if (method == "BMCP") {
      hull <-
        convHull(singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")], lonlat = TRUE)
      hull <- predict(Adeq, hull, mask = TRUE)

      pts1 <-
        singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")]
      spraster <- rasterize(pts1, Adeq, field = 1)
      sps <- as(spraster, 'SpatialPixels')@coords
      dist <- dist2(sps, sps, method = 'euclidean', p = 2)
      dist[dist == 0] <- NA
      distmin <- apply(dist, 1, function(x)
        min(x, na.rm = TRUE))

      hull2 <- hull
      hull2[hull2[] == 0] <- NA
      hull2 <- boundaries(hull2)
      hull2[hull2[] == 0] <- NA
      df <- rasterToPoints(hull2)
      df <- df[df[, 3] == 1,-3]
      buf <- circles(df, lonlat = TRUE, d = buffer)
      buf <- predict(Adeq, buf,  mask = TRUE)
      buf[(hull[] == 1)] <- 1
      buf[(!is.na(Adeq[]) & is.na(buf[]))] <- 0
      Adeq[which(buf[] != 1)] <- 0

      writeRaster(
        Adeq,
        paste(foldCon, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
      Mask <- (Adeq > Thr)
      writeRaster(
        Mask,
        paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
        format = "GTiff",
        overwrite = TRUE
      )
    }

    if (method %in% c("OBR", "LQ", "PRES")) {
      # Transform coordinate in a SpatialPoints object
      pts1 <-
        singleSpData[singleSpData[, "pres_abse"] == 1, c("x", "y")]
      coordinates(pts1) <- ~ x + y
      crs(pts1) <- crs(Adeq)

      # Raster with areas equal or grater than the threshold
      AdeqBin <- Adeq >= as.numeric(Thr)
      AdeqBin[AdeqBin[] == 0] <- NA
      # A "SpatialPolygonsDataFrame" which each adequability patch is a feature
      AdeqBin2 <-
        rasterToPolygons(
          AdeqBin,
          fun = NULL,
          n = 8,
          na.rm = TRUE,
          digits = 12,
          dissolve = TRUE
        )
      AdeqBin2 <- disaggregate(AdeqBin2)
      AdeqBin2$layer <- NULL
      # Individualize each patch with a number
      AdeqBin2$ID <- 1:length(AdeqBin2)
      # create a data.frame wiht coordinate and patch number
      AdeqPoints <- rasterToPoints(AdeqBin)[, 1:2]
      AdeqPoints <-
        cbind(AdeqPoints, ID = extract(AdeqBin2, AdeqPoints)[, 'ID'])
      # Find the patches that contain presences records
      polypoint <- intersect(AdeqBin2, pts1)

      # PRES methods------
      if (method == "PRES") {
        Adeq2 <- Adeq
        Msk <- rasterize(polypoint, Adeq, background = 0)
        Msk[is.na(Adeq[])] <- NA
        Adeq2[Msk != 1] <- 0

        writeRaster(
          Adeq2,
          paste(foldCon, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
        Mask <- Adeq2 >= Thr
        writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )

      } else{
        # Create a vector wich contain the number (e.i. ID) of the patches
        # with presences
        filter1 <- unique(polypoint$ID)
        # In this step are created two data.frame one with the patches coordinates
        # that contain presences and another with patches coordinates without presences
        CoordPathP <-
          as.data.frame(AdeqPoints[AdeqPoints[, 3] %in% filter1,])
        CoordPathNP <-
          as.data.frame(AdeqPoints[!AdeqPoints[, 3] %in% filter1,])
        # Here is created a matrix wiht all combination between ID of patches
        # with and without presences

        if (ncol(CoordPathP) == 1) {
          CoordPathP <- data.frame(t(CoordPathNP))
          rownames(CoordPathP) <- NULL
        }

        if (ncol(CoordPathNP) == 1) {
          CoordPathNP <- data.frame(t(CoordPathNP))
          rownames(CoordPathNP) <- NULL
        }
        npatch1 <- unique(CoordPathP[, 3])
        npatch2 <- unique(CoordPathNP[, 3])

        DistBetweenPoly0 <- expand.grid(npatch1, npatch2)
        DistBetweenPoly0$Distance <- NA
        DistBetweenPoly0 <- as.matrix(DistBetweenPoly0)
        # Euclidean Distance between patches wiht and without presences
        for (i in 1:nrow(DistBetweenPoly0)) {
          comb <- (DistBetweenPoly0[i, 1:2])
          A <- CoordPathP[CoordPathP[, 3] == comb[1], 1:2]
          B <- CoordPathNP[CoordPathNP[, 3] == comb[2], 1:2]

          if (nrow(A) >= 40) {
            SEQ <- round(seq(0, nrow(A), by = (nrow(A)) / 20))
            dist <- rep(NA, length(SEQ))
            for (j in 2:length(SEQ)) {
              SEQ2 <- (SEQ[(j - 1)] + 1):SEQ[j]
              dist[j] <-
                min(dist2(A[SEQ2, ], B, method = 'euclidean', p = 2), na.rm = TRUE)
            }
            eucdist <- min(dist[2:length(SEQ)], na.rm = TRUE)
          } else{
            eucdist <- min(dist2(A, B, method = 'euclidean', p = 2))
          }
          DistBetweenPoly0[i, 3] <- eucdist
        }

        DistBetweenPoly0 <-
          DistBetweenPoly0[order(DistBetweenPoly0[, 2]),]
        # Minimum Euclidean Distance between patches wiht and without presences
        DistBetweenPoly <-
          tapply(X = DistBetweenPoly0[, 3], DistBetweenPoly0[, 2], min)
        # Adding value of distance patches to cells
        AdeqBin2$Eucldist <- 0
        AdeqBin2$Eucldist[!AdeqBin2$ID %in% filter1] <-
          round(DistBetweenPoly, 4)

        # OBR method------
        if (method == 'OBR') {
          # method based on the maximum value of the minimum distance
          spraster <- rasterize(pts1, Adeq, field = 1)
          sps <- as(spraster, 'SpatialPixels')@coords
          dist <- dist2(sps, sps, method = 'euclidean', p = 2)
          dist[dist == 0] <- NA
          distmin <- apply(dist, 1, function(x)
            min(x, na.rm = TRUE))#
          CUT <- max(distmin)
        }
        # LQ method------
        if (method == "LQ") {
          # method based the lower quartile distance
          CUT <- c(summary(DistBetweenPoly0[, 3]))[2]
        }

        AdeqPoints <- rasterToPoints(AdeqBin)[, 1:2]
        AdeqPoints <- extract(AdeqBin2, AdeqPoints)[, 'Eucldist']
        fist <- AdeqBin
        fist[fist[] == 1] <- AdeqPoints
        # Threshold based on Maximum value of minimum ditance between ocurrences
        final <- fist <= CUT
        final[final == 0] <- NA
        Adeq2 <- Adeq
        Mask <- Adeq >= as.numeric(Thr)
        Mask[Mask == 1] <- 0
        Mask[!is.na(final[])] <- 1
        Adeq2[Mask != 1] <- 0
        # Save results as raster object
        writeRaster(
          Adeq2,
          paste(foldCon, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
        writeRaster(
          Mask,
          paste(foldCat, paste(SpNames[s], '.tif', sep = ""), sep = "/"),
          format = "GTiff",
          overwrite = TRUE
        )
      }
    }
  }
  }
