#' Find Latest Version of NWM on NCEP
#' @return character
#' @export

latest_nwm_version = function(){
  ncep = 'https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/'
  ver = grep("nwm", readLines(ncep), value = TRUE)
  gsub('^.*href=\"\\s*|\\s*/.*$', '', ver)
}

#' @title Download National Water Model Domain Files
#' @description Download the domain files (CONUS) used in the most
#' current operational NWM. It will download
#' ~17GB of data so plan accordingly!
#' @details This function will download the following national (CONUS) domain
#'  files for the latest version of the National Water Model on NCEP.
#' \itemize{
##'  \item{"Fulldom_hires_netcdf_1km.nc"}
##'  \item{"geo_em.d01_1km.nc"}
##'  \item{"wrfinput_d01_1km.nc"}
##'  \item{"RouteLink_NHDPLUS.nc"}
##'  \item{"spatialweights_1km_all_basins.nc"}
##'  \item{"GWBUCKPARM_CONUS.nc"}
##'  \item{"soil_veg_properties_LongRange.nc"}
##'  \item{"HYDRO_TBL_2D.nc"}
##'  \item{"WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc"}
##' }
##'
##' This function will NOT download:
##'\itemize{
##'  \item{"Fulldom_hires_netcdf_250m.nc"}("~13GB)
##'  \item{"GWBUCKPARM_CONUS_LongRange.nc"}("92B")
##'  \item{"HYDRO_TBL_2D_LongRange.nc"}{"338MB"}
##'  \item{"nudgingParams.nc"}{"2.4MB"}
##'  \item{"soil_veg_properties_ASM.nc"}{"3GB}
##'  \item{"soil_veg_properties_LongRange.nc"}{"3GB"}
##'  \item{"spatialweights_250m_all_basins.nc"}{"4.5GB"}
##' }
#' @param outDir where to write the files. Again ~17GB worth!
#' @importFrom httr GET write_disk progress
#' @return path to domain file downloads
#' @export

download_conus_nwm = function(outDir = NULL){

  ver = latest_nwm_version()

  nwmDir = paste0('https://www.nco.ncep.noaa.gov/pmb/codes/nwprod/',
                  ver,
                  '/parm/domain')

  needed_files = c("Fulldom_CONUS_FullRouting.nc",
                   "Fulldom_CONUS_LongRange.nc",
                   "GEOGRID_LDASOUT_Spatial_Metadata_CONUS.nc",
                   "GWBUCKPARM_CONUS_FullRouting.nc",
                   "GWBUCKPARM_CONUS_LongRange.nc",
                   "RouteLink_CONUS.nc",
                   "WRF_Hydro_NWM_geospatial_data_template_land_GIS.nc",
                   "geo_em.d01_1km.nc",
                   "geo_em_CONUS.nc",
                   "hydro2dtbl_CONUS_FullRouting.nc",
                   "hydro2dtbl_CONUS_LongRange.nc",
                   "nudgingParams_CONUS.nc",
                   "reservoir_index_AnA.nc",
                   "reservoir_index_Extended_AnA.nc",
                   "reservoir_index_Medium_Range.nc",
                   "reservoir_index_Short_Range.nc",
                   "reservoir_index_Standard_AnA.nc",
                   "soilproperties_CONUS_FullRouting.nc",
                   "soilproperties_CONUS_LongRange.nc",
                   "spatialweights_CONUS_FullRouting.nc",
                   "spatialweights_CONUS_LongRange.nc",
                   "wrfinput_CONUS.nc")
				   # This list of files were updated by Donny Kim to match the new filenames in repo.
  
  

  localDir = paste0(outDir, "/nwmCONUS-", gsub("nwm", "", gsub("[.]", "",ver)))

  dir.create(localDir, showWarnings = FALSE, recursive = TRUE)

  urls  = paste(nwmDir, needed_files, sep = "/")
  local = paste(localDir, needed_files, sep = "/")

  for( i in 1:length(urls)){
    if(!file.exists(local[i])){
      message('Downloading: ', basename(local[i]))
      httr::GET(urls[i], httr::write_disk(local[i]), httr::progress())
    } else {
      message(basename(local[i]),   ' already exists')
    }
  }

  message("\n\nAll files located: ", localDir)
  localDir
}

#' Get GeoGrid Projection
#' @description From a geogrid file, extract the proj4string.
#' @details Extracts projection data from global attributes
#' @param path a path to the local geogrid.nc file
#' @importFrom RNetCDF open.nc att.get.nc
#' @return a PROJ4STRING
#' @export

geo_grid_proj = function(path)
{
  suppressWarnings({
  nc = RNetCDF::open.nc(path)
  map_proj <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", "MAP_PROJ")
  cen_lat  <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", "CEN_LAT")
  cen_lon  <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", "STAND_LON")
  truelat1 <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", "TRUELAT1")
  truelat2 <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", "TRUELAT2")
  if (map_proj == 1) {
    crs.here <- paste0("+proj=lcc +lat_1=", truelat1,
                       " +lat_2=", truelat2,
                       " +lat_0=", cen_lat,
                       " +lon_0=", cen_lon,
                       " +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs")
  }
  return(crs.here)
  })
}


#' Return GeoGrid Raster Structure
#' @description Given a NetCDF path, a raster object is returned with the
#' appropriate structure, valeus are NA, unless a 'var' is defined. S
#' ee var_names for help indetifing availiable varible names
#' @param path a path to a NetCDF file
#' @param var an optional variable to extract
#' @return a rasterLayer
#' @export
#' @importFrom RNetCDF open.nc dim.inq.nc att.get.nc var.get.nc
#' @importFrom sf st_point st_sfc st_transform st_coordinates
#' @importFrom raster raster values
#' @importFrom magrittr %>%

make_empty_geogrid_raster = function(path, var = NULL)
{
  crs.here = geo_grid_proj(path)
  nc = RNetCDF::open.nc(path)
  x  = RNetCDF::dim.inq.nc(nc, 'west_east')$length
  y  = RNetCDF::dim.inq.nc(nc, 'south_north')$length
  dx = RNetCDF::att.get.nc(nc, "NC_GLOBAL", "DX")
  dy = RNetCDF::att.get.nc(nc, "NC_GLOBAL", "DY")
  lat_max = RNetCDF::var.get.nc(nc, "XLAT_M", start = c(1,y,1), count = c(1,1,1))
  lng_max = RNetCDF::var.get.nc(nc, "XLONG_M", start = c(1,y,1), count = c(1,1,1))

  pts = st_sfc(st_point(c(lng_max,lat_max)), crs = 4326) %>%
    st_transform(crs.here) %>%
    st_coordinates()

  xmin = pts[1] - dx/2
  ymax = pts[2] + dy/2
  xmax = xmin + x*dx
  ymin = ymax - y*dy

  suppressWarnings({
    r = raster(res = c(dx,dy), xmn = xmin, xmx = xmax,
               ymn = ymin, ymx = ymax,
               crs = crs.here)
  })

  if(!is.null(var)){
    var = RNetCDF::var.get.nc(nc, var)
    raster::values(r) = apply(t(var),2,rev)
  }

  r
}


#' @title NetCDF Variable Names from a path
#' @description Given a path to a NetCDF file, list the variables
#' @param path a path to a NetCDF file
#' @return a vectore of variable names
#' @export
#' @importFrom RNetCDF open.nc file.inq.nc var.inq.nc

var_names = function(path){
  nc = RNetCDF::open.nc(path)
  nvar <- RNetCDF::file.inq.nc(nc)$nvar
  varnames <- character(nvar)
  for(i in seq_len(nvar)) {
    varnames[i] <- RNetCDF::var.inq.nc(nc, i-1)$name
  }
  varnames
}

#' @title NetCDF Dimension Names from a path
#' @description Given a path to a NetCDF file, list the diminsions
#' @param path a path to a NetCDF file
#' @return a vector of dimension names
#' @export
#' @importFrom RNetCDF open.nc file.inq.nc var.inq.nc

dim_names = function(path){
  nc = RNetCDF::open.nc(path)
  ndim <- RNetCDF::file.inq.nc(nc)$ndim
  dimnames <- character(ndim)
  for(i in seq_len(ndim)) {
    dimnames[i] <- RNetCDF::dim.inq.nc(nc, i-1)$name
  }
  dimnames
}


#' @title Re-write of rwfhydro ReadWtFile
#' @description  reads in a spatial weight file and generates
#' a list of two dataframes: 1 for the gridded info (data) and 2 for
#' the poly info (poly)
#' @param path path to the weight file
#' @return list of 2 data.frames
#' @export
#' @importFrom RNetCDF open.nc file.inq.nc var.inq.nc
#
read_wt_file = function (wtFile)
{
  nc <- RNetCDF::open.nc(wtFile)

   list(data = data.frame(
    i_index = RNetCDF::var.get.nc(nc, "i_index"),
    j_index = RNetCDF::var.get.nc(nc, "j_index"),
    IDmask  = RNetCDF::var.get.nc(nc, "IDmask"),
    weight  = RNetCDF::var.get.nc(nc, "weight"),
    regridweight = RNetCDF::var.get.nc(nc, "regridweight")),

    polys = data.frame(polyid = RNetCDF::var.get.nc(nc, "polyid"),
                      overlaps = RNetCDF::var.get.nc(nc, "overlaps"))
  )

}



#' @title Re-write of rwrfhydro::SubsetWts
#' @description SubsetWts takes a weight file object (as read in from ReadWtFile)
#' and subsets the data and poly dataframes based on a user-defined set of link
#' IDs and user-defined start and end indices. Returns a list of 2 dataframes
#' for the subsetted data and poly objects.
#' @param wts Weight file list object
#' @param rlids Vector of link IDs to keep in the subsetted output
#' @param istart Start index for the grid in the i (x) direction
#' @param iend End index for the grid in the i (x) direction
#' @param jstart Start index for the grid in the j (y) direction
#' @param jend End index for the grid in the j (y) direction
#' @return List of 2 dataframes
#' @importFrom magrittr %>%
#' @export

subset_weights = function (wts, rlids, istart, iend, jstart, jend)
{

  wtssubrect = dplyr::filter(wts$data,
                         IDmask %in% rlids) %>%
    dplyr::mutate(i_index = i_index - (istart - 1),
           j_index =j_index - (jstart - 1)) %>%
    dplyr::filter(
      i_index > 0 &
      i_index <= (iend - istart + 1) &
      j_index > 0 &
      j_index <= (jend - jstart + 1)
    )

  polywts2 = wtssubrect %>%
    dplyr::group_by(IDmask) %>%
    dplyr::summarise(mult = sum(weight))

  wtssubrect = dplyr::left_join(
    wtssubrect, polywts2, by = "IDmask"
  ) %>%
    dplyr::mutate(weight = weight / mult,
                  mult = NULL)

  polysub = dplyr::filter(wts$polys, polyid %in% rlids)
  names(polysub) <- c("polyid", "kill")

  tmp = dplyr::count(wtssubrect, IDmask)
  names(tmp) <- c("polyid", "overlaps")

  polysub <- dplyr::left_join(polysub, tmp, by = "polyid") %>%
    dplyr::mutate(overlaps = ifelse(is.na(overlaps), 0, overlaps),
                  kill = NULL)

  return(list(wtssubrect, polysub))
}




#' @title Update the attribute of NetCDF file
#' @description Update the NetCDF file in path with target dataframe.
#' @param path Directory of NetCDF file to be updated
#' @param df R dataframe  of link IDs to keep in the subsetted output
#' @return Nothing
#' @export
#' 
update_nc = function(path, df){

  nc <- RNetCDF::open.nc(path, write = TRUE)
  for (i in names(df)) {
    if (i %in% var_names(path)) {
      if(RNetCDF::var.inq.nc(nc,i)$ndim > 0){
        print(i)
        RNetCDF::var.put.nc(nc, var = i, df[[i]])
      }
    }
  }
  RNetCDF::close.nc(nc)

}



#' @title Read a NetCDF file into a data.frame
#' @param path a path to a NetCDF
#' @return a data.frame
#' @importFrom dplyr bind_cols
#' @importFrom RNetCDF open.nc var.get.nc
#' @export

nc_to_df = function(path){
  nc = RNetCDF::open.nc(path)
  vars = var_names(path)
  out = lapply(1:length(vars), function(i){
    RNetCDF::var.get.nc(nc, vars[i])
  })

  out = suppressMessages(dplyr::bind_cols(out))
  names(out) = vars
  out
}

