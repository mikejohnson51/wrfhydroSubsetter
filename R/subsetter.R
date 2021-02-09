#' @title WRF-Hydro Subsetter
#' @description Create a cutout of the NWM domain files
#' @param AOI sf POLYGON. The area to cut out
#' @param domain_files character. A path to the local domain files
#' @param outDir Where to put the cutout data should be a named directory!
#' @return files on disk
#' @export
#' @importFrom sf st_bbox st_transform st_as_sf as_Spatial
#' @importFrom raster rowColFromCell cellFromXY cellFromRowCol xyFromCell
#' @importFrom dplyr filter mutate left_join group_by summarize n

subset_wrf_hydro_domain = function(AOI, domain_files, outDir){

  # create the direcotry if does not exist.
  myPath = paste0(outDir, "/")
  dir.create(myPath, showWarnings = FALSE, recursive = TRUE)

#******  Specify the path to the ORIGINAL (full extent) domain files:
  fullHydFile      <- paste0(domain_files, "/Fulldom_hires_netcdf_1km.nc")
  fullGeoFile      <- paste0(domain_files, "/geo_em.d01_1km.nc")
  fullWrfFile      <- paste0(domain_files, "/wrfinput_d01_1km.nc")
  fullRtlinkFile   <- paste0(domain_files, "/RouteLink_NHDPLUS.nc")
  fullSpwtFile     <- paste0(domain_files, "/spatialweights_1km_all_basins.nc")
  fullGwbuckFile   <- paste0(domain_files, "/GWBUCKPARM_CONUS.nc")
  fullSoilparmFile <- paste0(domain_files, "/soil_veg_properties_LongRange.nc")
  fullHydro2dFile  <- paste0(domain_files, "/HYDRO_TBL_2D.nc")
  geoSpatialFile   <- paste0(domain_files, "/WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc")

  coordProj = geo_grid_proj(fullGeoFile)

  # Multiplier between routing grid and LSM grid
  # (e.g., 1-km LSM and 250-m routing means a value of 4)
  dxy <- 1

  # Specify the NEW (subset extent) domain files:
  # Routing domain file
  subHydFile        <- paste0(myPath, "/Fulldom_hires.nc")
  subGeoFile        <- paste0(myPath, "/geo_em.d0x.nc")
  subWrfFile        <- paste0(myPath, "/wrfinput_d0x.nc")
  subRtlinkFile     <- paste0(myPath, "/Route_Link.nc")
  subSpwtFile       <- paste0(myPath, "/spatialweights.nc")
  subGwbuckFile     <- paste0(myPath, "/GWBUCKPARM.nc")
  subSoilparmFile   <- paste0(myPath, "/soil_properties.nc")
  subHydro2dFile    <- paste0(myPath, "/hydro2dtbl.nc")
  subGeoSpatialFile <- paste0(myPath, "/GEOGRID_LDASOUT_Spatial_Metadata.nc")


  if(
    all(file.exists(c(subHydFile, subGeoFile,
                  subWrfFile, subRtlinkFile,
                  subSpwtFile, subGwbuckFile,
                  subSoilparmFile, subHydro2dFile,
                  subGeoSpatialFile)))) {
   return(myPath)
  }

  # Setup coordinates df
  bb     <- st_bbox(st_transform(AOI, coordProj))
  coords <- data.frame(id=c(1,2,3,4),
                       lat=c(bb$ymin, bb$ymax, bb$ymax, bb$ymin),
                       lon=c( bb$xmin,  bb$xmin, bb$xmax, bb$xmax))

  # Generate spatial coords
  bb_pts  = st_as_sf(coords, coords = c("lon", 'lat'), crs = coordProj)

  r = make_empty_geogrid_raster(fullGeoFile)

  geoindex <- as.data.frame(raster::rowColFromCell(r, raster::cellFromXY(r, as_Spatial(bb_pts))))

  geoindex$we <- geoindex$col
  geoindex$sn <- nrow(r) - geoindex$row + 1
  geoindex$id <- coords[,"id"]

  # Get subsetting dimensions
  geo_w <- min(geoindex[,"we"])
  geo_e <- max(geoindex[,"we"])
  geo_s <- min(geoindex[,"sn"])
  geo_n <- max(geoindex[,"sn"])

  hyd_w <- (geo_w-1)*dxy+1
  hyd_e <- geo_e*dxy
  hyd_s <- (geo_s-1)*dxy+1
  hyd_n <- geo_n*dxy

  hyd_min <- (min(geoindex$row)-1)*dxy+1
  hyd_max <- max(geoindex$row)*dxy
  geo_min <- min(geoindex$row)
  geo_max <- max(geoindex$row)

  sp_new_buff_nad83 <- st_as_sf(xyFromCell(r,
                                           cellFromRowCol(r,
                                                          c(geo_max, geo_min, geo_min, geo_max),
                                                          c(geo_w, geo_w, geo_e, geo_e)), spatial=TRUE)) %>%
    st_transform("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs") %>%
    st_coordinates()

################# SUBSET DOMAIN Files
  # ROUTING GRID
  if (!is.null(fullHydFile)) {
    if  (!file.exists(fullHydFile)) stop(paste0("The fullHydFile : ", fullHydFile, " does not exits"))
    system(paste0("ncks -O -d x,", hyd_w-1, ",", hyd_e-1, " -d y,", hyd_min-1, ",", hyd_max-1, " ", fullHydFile, " ", subHydFile))
  }

  # Geo Spatial File
  if (!is.null(geoSpatialFile)) {
    if (!file.exists(geoSpatialFile)) stop(paste0("The geoSpatialFile :", geoSpatialFile, " does not exits"))
    system(paste0("ncks -O -d x,", geo_w-1, ",",
                  geo_e-1, " -d y,", geo_s-1, ",",
                  geo_n-1, " ", geoSpatialFile, " ",
                  subGeoSpatialFile))
    }


  # GEO GRID
  if (!is.null(fullGeoFile)) {
    if (!file.exists(fullGeoFile)) stop(paste0("The fullGeoFile :", fullGeoFile, " does not exits"))
    system(paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1,
                  " -d west_east_stag,", geo_w-1, ",", geo_e, " -d south_north_stag,",geo_s-1, ",", geo_n, " ",
                  fullGeoFile, " ", subGeoFile))
  }

  corner_lats <- c()

  for (latName in c('XLAT_M', 'XLAT_U', 'XLAT_V', 'XLAT_C')) {
    if (latName %in% var_names(subGeoFile)) {
      a = RNetCDF::var.get.nc(RNetCDF::open.nc(subGeoFile), latName)
      corners = c(a[1,1], a[1, ncol(a)], a[nrow(a), ncol(a)], a[nrow(a), 1])
    }else{
      corners = c(0,0,0,0)
    }
    corner_lats = c(corner_lats, corners)
  }

  corner_lons <- c()
  for (lngName in c('XLONG_M', 'XLONG_U', 'XLONG_V', 'XLONG_C')) {
    if (lngName %in% var_names(subGeoFile)) {
      a = RNetCDF::var.get.nc(RNetCDF::open.nc(subGeoFile), lngName)
      corners = c(a[1,1], a[1, ncol(a)], a[nrow(a), ncol(a)], a[nrow(a), 1])
    }else{
      corners = c(0,0,0,0)
    }
    corner_lons = c(corner_lons, corners)
  }

  # Attribute updates
  system(paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subGeoFile))
  system(paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subGeoFile))
  system(paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subGeoFile))
  system(paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subGeoFile))
  system(paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subGeoFile))
  system(paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subGeoFile))
  system(paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subGeoFile))
  system(paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subGeoFile))
  system(paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subGeoFile))
  system(paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subGeoFile))
  system(paste0("ncatted -O -a corner_lons,global,o,f,", paste(corner_lons, collapse  = ","), " ", subGeoFile))
  system(paste0("ncatted -O -a corner_lats,global,o,f,", paste(corner_lats, collapse  = ","), " ", subGeoFile))

  #HYDRO_TBL_2D GRID
  if (!is.null(fullHydro2dFile)) {
    if (!file.exists(fullHydro2dFile)) stop(paste0("The fullHydro2dFile : ", fullHydro2dFile, " does not exits"))
    system(paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullHydro2dFile, " ", subHydro2dFile))
  }

  # WRFINPUT GRID

  if (!is.null(fullWrfFile)) {
    if (!file.exists(fullWrfFile)) stop(paste0("The fullWrfFile : ", fullWrfFile, " does not exits"))
    system(paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullWrfFile, " ", subWrfFile))
    system(paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subWrfFile))
    system(paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subWrfFile))
  }

  ################# SUBSET PARAMS
if (!file.exists(fullSpwtFile)) stop(paste0("The fullSpwtFile : ", fullSpwtFile, " does not exits"))
if (!file.exists(fullRtlinkFile)) stop(paste0("The fullRtlinkFile : ", fullRtlinkFile, " does not exits"))

    SpatWts <- read_wt_file(fullSpwtFile)

    keepIdsPoly <-  dplyr::filter(SpatWts$data,
                         i_index >= hyd_w &
                         i_index <= hyd_e &
                         j_index >= hyd_s &
                         j_index <= hyd_n) %>%
      dplyr::group_by(IDmask) %>%
      dplyr::summarise(sumBas = sum(weight)) %>%
      dplyr::filter(sumBas > .999) %>%
      dplyr::pull(IDmask)

    rt_link = nc_to_df(fullRtlinkFile)

    keepIdsLink = rt_link %>%
      dplyr::filter(
          lon >= min(sp_new_buff_nad83[,1]) &
          lon <= max(sp_new_buff_nad83[,1]) &
          lat >= min(sp_new_buff_nad83[,2]) &
          lat <= max(sp_new_buff_nad83[,2])) %>%
      dplyr::distinct(link) %>%
      dplyr::pull(link)

    keepIds <- unique(c(keepIdsPoly, keepIdsLink))

    # SPATIAL WEIGHT
    subWts = subset_weights(SpatWts, keepIdsPoly, hyd_w, hyd_e, hyd_s, hyd_n)

    fs::file_copy(fullSpwtFile, subSpwtFile, overwrite = TRUE)

    system(paste0("ncks -O -d polyid,1,", nrow(subWts[[2]]),
                       " -d data,1,", nrow(subWts[[1]]), " ", subSpwtFile, " ",
                     subSpwtFile))

    update_nc(subSpwtFile, subWts[[1]])
    update_nc(subSpwtFile, subWts[[2]])

    # ROUTE LINK
    subRtlink <- dplyr::filter(rt_link, link %in% keepIds) %>%
      dplyr::mutate(to = ifelse(to %in% unique(link), to, 0))

    # reorder the ascendingIndex if ascendingIndex exists in the variables
    if ("ascendingIndex" %in% names(subRtlink)){
      subRtlink = dplyr::mutate(subRtlink, ascendingIndex = (rank(ascendingIndex) - 1))
    }

    #rwrfhydro::UpdateLinkFile(subRtlinkFile, subRtlink)
    fs::file_copy(fullRtlinkFile, subRtlinkFile, overwrite = TRUE)

    system(paste0("ncks -O -d feature_id,1,", nrow(subRtlink), " ",
                  subRtlinkFile, " ", subRtlinkFile))

    update_nc(path = subRtlinkFile, df = subRtlink)

  # GWBUCK PARAMETER

if (!is.null(fullGwbuckFile)) {

    if (is.null(fullSpwtFile) | is.null(fullRtlinkFile)) {
      stop("To subset the fullSpwtFile, you need fullSpwtFile and fullRtlinkFile")
    }
    if (!file.exists(fullGwbuckFile)) stop(paste0("the fullGwbuckFile : ",
                                                  fullGwbuckFile, " does not exits"))

   subGwbuck = nc_to_df(fullGwbuckFile) %>%
      dplyr::filter(ComID %in% keepIdsPoly) %>%
      dplyr::mutate(Basin = 1:n())

    file.copy(fullGwbuckFile, subGwbuckFile, overwrite = TRUE)

    system(paste0("ncks -O -d BasinDim,1,", nrow(subGwbuck),  " ",
                  subGwbuckFile, " ", subGwbuckFile))

    update_nc(path = subGwbuckFile, df = subGwbuck)

  }

  # SOIL PARAMETER

  if (!is.null(fullSoilparmFile)) {
    if (!file.exists(fullSoilparmFile)) stop(paste0("the fullSoilparmFile : ", fullSoilparmFile, " does not exits"))
    system(paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullSoilparmFile, " ", subSoilparmFile))
  }

  message("Finshied! All files in: ", myPath)

}

