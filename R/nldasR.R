#' @title Download CONUS NLDAS Forcing Files
#' @description Download CONUS wide NLDAS forcing files for local use. Hourly files are downloaded for each day between the requested start and end date. One year of data requires ~13 GB of disk space. Each file is ~1.8 MB.
#' @details This function will download the following national or subset NLDAS forcing files for the latest version of the National Water Model on NCEP.
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
#'
#' @param startDate
#' @param endDate
#' @param outDir
#' @param netrcFile
#'
#' @return
#' @export
#'
#' @examples
#'

#type = c("FORA", "FORB", "NOAH", "MOS", "VIC")
timestep = c("M", "H")
startDate = "1991-01-01"
AOI = AOI::aoi_get(state = "OR")

download_nldas = function(AOI = NULL,
                          startDate = NULL, endDate = NULL,
                          type = "FORA", timestep = "H",
                          outDir,
                          netrcFile = getNetrcPath()){

  if(is.null(endDate)){endDate = startDate}
  startDate = paste(startDate, "00:00:00")
  endDate = paste(endDate, "23:00:00")

  df = data.frame(date = seq.POSIXt(as.POSIXct(startDate),
                                    as.POSIXct(endDate), by = 'h')) %>%
  mutate(jul   = format(date, "%j"),
         year  = format(date, "%Y"),
         day   = format(date, "%d"),
         month = format(date, "%m"),
         hour  = format(date, "%H"))

  mod = paste0('NLDAS_', type,'0125_',timestep)

  if(is.null(AOI)){
    df = mutate(df, urls  = paste0('https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/',
                       mod,    '.002/',
                       year,
                       "/",
                       jul,
                       '/',mod,'.A',
                       year,
                       month,
                       day,
                       ".",
                       hour,
                       '00.002.grb'))
  } else {
    bb = sf::st_bbox(AOI)

    df = mutate(df, urls  = paste0('https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FNLDAS%2F',
                          mod,'.002%2F',
                          year,'%2F',jul,'%2F',
                          mod,'.A',year,month,day,".",hour,"00.002.grb&FORMAT=bmM0Lw&",
                          'BBOX=',bb[2],'%2C',bb[1],'%2C',bb[4],'%2C',bb[3],
                          '&LABEL=',mod,'.A',year,month,day,'.',hour,
                          '00.002.grb.SUB.nc4&SHORTNAME=',mod,
                          '&SERVICE=L34RS_LDAS&VERSION=1.02&DATASET_VERSION=002'))
  }

message(nrow(df), " files requested")

for(i in 1:nrow(df)){

  dir = paste(outDir, df$year[i], df$month[i], df$day[i], sep = "/")

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  here = paste0(dir, "/NLDAS2_",df$year[i], df$month[i], df$day[i], df$hour[i], ".grb")

  if(!file.exists(here)){
    httr::GET(df$urls[i], httr::write_disk(here, overwrite = TRUE),
              httr::config(netrc = TRUE, netrc_file = netrcFile),
              httr::set_cookies("LC" = "cookies"))
  }

  message(i)
}

}


#s= stars::read_stars('/Users/mikejohnson/Downloads/force/1991/01/01/NLDAS2_1991010123.grb')

#plot(s)

#raster::plot(raster::raster('/Users/mikejohnson/Downloads/force/1991/01/01/NLDAS2_1991010123.nc', var.name = "TMP"))
