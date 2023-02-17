#' @title Subset sequence wrapper
#' @description Initiates WrfhydroSubsetter workflow
#' @param comid comid of outlet reach
#' @param siteID siteID for NWIS
#' @param FULLDOMAINDIR the value representing "NO DATA", default is NA
#' @param outDir output directory (e.g. /home/subsetDOMAINS/)
#' @param nlcdDir directory for nlcd dataset
#' @param methodList resampling algorithm to use for resampling
#' @importFrom sf st_transform st_buffer
#' @importFrom AOI aoi_get bbox_get
#' @importFrom dataRetrieval findNLDI
#' @importFrom spex qm_rasterToPolygons
#' @importFrom raster crop
#' @importFrom ncdf4 nc_open ncvar_put nc_close
#' @export

subset_sequence = function(comid, siteID, FULLDOMAINDIR, outDir, nlcdDir = NA, methodList = NA){
  
  # Defining the area * This could go out and just stay in the loop.
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = comid), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, comid, siteID, sep = "_")))
  basin = findNLDI(comid = comid, find = c("basin"))$basin
  
  # Subsetting NetCDF DOMAIN files
  subset_files = paste0(outDir, "/", name , "/")
  subset_wrf_hydro_domain(AOI = basin,  domain_files = FULLDOMAINDIR,  outDir = subset_files, config = 'LongRange')
  
  
  # Resampling sequence
  if (exists("methodList")==T & exists("nlcdDir")==T) {
    
    # Reading in look-up table. This is somewhat problematic at this point.
    data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
      dplyr::select(nlcd = Class, 
                    description = Description,
                    nwm = NWM)
    
    # Cropping nlcd to geogrid resolution
    geo  = list.files(subset_files, "geo_em.nc", full.names = TRUE)
    output_geo = spex::qm_rasterToPolygons(wrfhydroSubsetter::make_empty_geogrid_raster(geo))
    o2 = output_geo %>%
      AOI::bbox_get() %>%
      sf::st_transform(5070) %>%
      sf::st_buffer(40)
    nlcdObj = raster(nlcdDir)
    methodList = methodList
    nlcdCrop = raster::crop(nlcdObj, o2)
    
    # Reclassifying nlcd LULC class into nwm class *DK: is it necessary to do this at this phase tho?
    nlcd_nwm = raster::reclassify(nlcdCrop, rcl = dplyr::select(data,from = nlcd, to = nwm ))
    
    # Resampling sequence.
    for(j in 1:length(methodList)){
      new_lu = resampleDataMod(input = nlcd_nwm, output = output_geo, method = methodList[j])
      out_file = paste0(subset_files, "geo_", methodList[j], ".nc")
      
      file.copy(geo, out_file, overwrite = TRUE)
      nc = ncdf4::nc_open(out_file, write = TRUE)
      ncdf4::ncvar_put(nc, "LU_INDEX",
                       vals = as.vector(t(apply(as.matrix(new_lu), 2, rev))))
      ncdf4::nc_close(nc)
    }
  }
}




#' @title mask_geogrid_byBasin
#' @description Mask GEOGRID by Watershed Boundary. This can be used to perturb LULC within WS boundary or to calculate LULC statistics within WS.
#' @param RasterLayer GEOGRID as RasterLayer stack
#' @param basin Watershed boundary as sf
#' @importFrom raster rasterize mask stack
#' @importFrom magrittr %>%
#' @return a dataframe


mask_geogrid_byWS = function(RasterLayer, basin){
  require(raster)
  require(dplyr)
  
  basin_coord = st_transform(basin, st_crs(RasterLayer))

  # Regular "raster::mask" function only selects  cells those have their centroid falls inside the polygon (watershed boundary).
  # We want to choose EVERY cell that intersect with watershed boundary (polygon).
  basin_coord_ras <- raster::rasterize(basin_coord, RasterLayer, getCover=TRUE) # getCover is the parameter that needs to be set True
  basin_coord_ras[basin_coord_ras==0] = NA
  RasterLayer_mask <- stack(raster::mask(RasterLayer, basin_coord_ras))
  RasterLayer_mask_df <- stack(raster::mask(RasterLayer, basin_coord_ras)) %>% as.data.frame(xy = TRUE)
  # Resources:
  # https://stackoverflow.com/questions/68781519/errors-when-clipping-raster-stacks-with-sfst-crop-and-rastercrop
  # https://gis.stackexchange.com/questions/255025/r-raster-masking-a-raster-by-polygon-also-remove-cells-partially-covered

  return(RasterLayer_mask_df)  
}