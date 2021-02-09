library(sf)
library(lubridate)
library(RNetCDF)
library(dplyr)
library(ggplot2)


## Base Material
#=======================================================
subPath = '/Volumes/Transcend/land-cover-tests/berkeley_west-virginia_5894384_01616500/geo_em.d0x.nc'

geogrid.raster = wrfhydroSubsetter::make_empty_geogrid_raster(subPath)
basin = dataRetrieval::findNLDI(nwis = '01616500', find = c("basin"))$basin

maskBas = fasterize::fasterize(st_transform(basin, st_crs(geogrid.raster)), geogrid.raster)

#These are the variable names, descriptions and units we need
df = data.frame(variableNames = c("ACCECAN", "ACCEDIR", "ACCETRAN", "ACCPRCP",
                                  "CANICE", "CANLIQ", "SFCRUNOFF", "UDRUNOFF", "SNEQV"),

                description = c("Accumulated canopy evaporation",
                                "Accumulated direct soil evaporation",
                                "Accumulated transpiration",
                                "Accumulated precipitation",
                                "Canopy ice water content",
                                "Canopy liquid water content",
                                "Accumulated surface runoff",
                                "Accumulated underground runoff",
                                "Snow water equivalent"),

                units = c("mm", "mm", "mm", "mm", "mm", "mm", "mm", "mm", "kgm2"))

#=======================================================
## FUNCTION WILL START HERE:
## INPUTS ARE FILELIST (dir), soil depths (vector, mm) and maskBas (raster):
#=======================================================

# USER PROVIDES DIRECTORY, TO OUTPUT
dir  = '/Volumes/Transcend/land-cover-tests/berkley-spinup'
# USER PROVIDES SOIL DEPTHS
soil_depths_mm = c(100, 300, 600, 1000)

# Read all files in directory, and pull out date and water year
fileList = data.frame(
  files = list.files(dir, pattern = "RESTART", full.names = TRUE)) %>%
  dplyr::mutate(date  = lubridate::ymd_h(gsub("_DOMAIN1", "", gsub("^.*\\.","", files)), tz = "UTC"),
                wy = ifelse(lubridate::month(date) >= 10,
                            lubridate::year(date) + 1,
                            lubridate::year(date)))

# The current water budget formulation is the first minus the last timestep.
# So, lets only process those...
fileList = fileList[c(1,nrow(fileList)),]

## Mask of 1, NA depending on if cell is in basin
maskMatrix = raster::as.matrix(maskBas)

#=======================================================

#Empty list to populate
out = list()

for(i in 1:nrow(fileList)){
  # open the current netcdf here
  tmp.nc = RNetCDF::open.nc(fileList$files[i])

  # Begin surface variable extraction
  surface_extract = lapply(seq_along(df$variableNames), function(i){
    # Read in variable, transpose, multiply by basin mask, and average
    # Result is average *unit* of *variable* in the basin at time of file
    mean(apply(t(RNetCDF::var.get.nc(tmp.nc, df$variableNames[i])), 2, rev) *
           maskMatrix, na.rm = TRUE)
  })

  # Open 3D SMC variable
  soil = RNetCDF::var.get.nc(tmp.nc, "SMC")

  # Begin surface variable extraction
  soil_extract = lapply(1:length(soil_depths_mm), function(i){
    # Read in soil layer, transpose, multiply by basin mask, average,
    # multiply by soil depth of layer
    # Result is average mm in each layer in the basin at time of file
    mean(apply(t(soil[,i,]), 2, rev) *
           maskMatrix, na.rm = TRUE) *
      soil_depths_mm[i]
  })

  # 1 row data.frame, 1 date, all variables
  # add these to the 'out' list
  out[[i]] = data.frame(t(c(unlist(surface_extract),
                            unlist(soil_extract)))) %>%
    setNames(c(variableNames,
               paste0('SOIL_', 1:length(soil_depths_mm))))
  RNetCDF::close.nc(tmp.nc)
}

# bind all rows in out list, column bind to fileList
o2          = cbind(fileList, dplyr::bind_rows(out))

# Not used in water balance calc
# # Build delta columns with a lag of 1
# for(i in seq_along(df$variableNames)){
#   var = df$variableNames[i]
#   o2[[paste0("DEL_", var)]] = dplyr::lag(o2[[var]], n = 1)
# }

# Here we are summing all soil columns
TOT_SOIL = rowSums(dplyr::select(o2, dplyr::starts_with("SOIL")))

# For all states subtract the inital from the final
# For canopy, ice and liquid water content are combined

wbDf = data.frame(
  # RAINFALL
  LSM_PRCP      = o2$ACCPRCP[nrow(o2)]  - o2$ACCPRCP[1],
  # Canopy evaporation
  LSM_ECAN      = o2$ACCECAN[nrow(o2)]  - o2$ACCECAN[1],
  # Canopy transporations
  LSM_ETRAN     = o2$ACCETRAN[nrow(o2)] - o2$ACCETRAN[1],
  # Soil evaporation
  LSM_EDIR      = o2$ACCEDIR[nrow(o2)]  - o2$ACCEDIR[1],
  # Snow water
  LSM_DELSWE    = o2$SNEQV[nrow(o2)]    - o2$SNEQV[1],
  # Canopy water/ice
  LSM_DELCANWAT = (o2$CANICE[nrow(o2)]  + o2$CANLIQ[nrow(o2)]) - (o2$CANICE[1] + o2$CANLIQ[1]),
  # Surface runoff
  LSM_SFCRNOFF  = o2$SFCRUNOFF[nrow(o2)] -  o2$SFCRUNOFF[1],
  # Underground runoff
  LSM_UGDRNOFF  = o2$UDRUNOFF[nrow(o2)]  -  o2$UDRUNOFF[1],
  # Soil moisture
  LSM_DELSOILM = TOT_SOIL[length(TOT_SOIL)] - TOT_SOIL[1]
) %>% dplyr::mutate(
  # Since we are only using the LSM we set the WB to the LSM
    WB_SFCRNOFF =  LSM_SFCRNOFF,
    WB_GWOUT    =  LSM_UGDRNOFF,
    # I dont get this but its not used so who cares
    #WB_DELGWSTOR = (o2$UDRUNOFF[nrow(o2)] - o2$UDRUNOFF[1]) - WB_GWOUT,
    ERROR = LSM_PRCP -
            LSM_ECAN - LSM_ETRAN - LSM_EDIR -
            WB_SFCRNOFF - WB_GWOUT - LSM_DELSOILM - LSM_DELSWE - LSM_DELCANWAT,
    RUN_FRAC = (WB_SFCRNOFF + WB_GWOUT)/LSM_PRCP,
    EVAP_FRAC = (LSM_ECAN + LSM_ETRAN + LSM_EDIR)/LSM_PRCP,
    STOR_FRAC = (LSM_DELSOILM + LSM_DELSWE + LSM_DELCANWAT)/LSM_PRCP)


#=======================================================
## To REPLICATE WB BARPLOT:
#=======================================================



to_plot  = data.frame(
  class = c(
    "Canopy Evap",
    "Transpiration",
    "Surface Evap",
    "Surface Runoff",
    "Groundwater Outflow",
    "Change in Storage"),
  pcts = with(
    wbDf,
    c(LSM_ECAN / LSM_PRCP * 100,
      LSM_ETRAN / LSM_PRCP * 100,
      LSM_EDIR / LSM_PRCP * 100,
      WB_SFCRNOFF / LSM_PRCP * 100,
      WB_GWOUT / LSM_PRCP * 100,
      (LSM_DELSOILM + LSM_DELSWE + LSM_DELCANWAT) / LSM_PRCP * 100))
) %>%
  dplyr::mutate(labels = paste0(class, "\n", round(pcts, 1), "%"), basin = "1") %>%
  dplyr::arrange(pcts)

ggplot(data = to_plot, aes(y = pcts, x = basin, fill = labels)) +
  geom_col(position="fill") +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(fill = "", x = "", y = "% of Water Budget")


