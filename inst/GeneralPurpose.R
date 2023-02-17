{
  library(dataRetrieval); library(nhdplusTools);
  library(raster); library(sf); library(fasterize);
  library(wrfhydroSubsetter); library(AOI);
  library(tidyverse); library(dplyr);
  #library(resample);
}

# This script is prone to crash during debugging. In such cases, just restart R session. It solves most (if not all) of the problems.
.rs.restartR()

### STEP 0======================================================================
## For package development purpose:
#library(devtools)
#devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
library(wrfhydroSubsetter)
#detach("package:wrfhydroSubsetter", unload=TRUE); remove.packages("wrfhydroSubsetter");
### STEP 0 END======================================================================


### STEP 1======================================================================
# Setting up locations
locs = data.frame(comids = c(191739, 1631587, 5894384, 5781369, 19389766, 23762661),
                  siteID = c('6709000', '08173000', '01616500', '08159000' ,'03118500' , '14190500'))

#test_loc = data.frame(comids = c(5894384), siteID = c('01616500')) # In case using a single watershed

namelist = c()
for (i in 1:nrow(locs)){
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  namelist = c(namelist, name)
  rm(name)
  }

# Setting directories
outDir = "/mnt/d/subsetDOMAINS/"
FULLDOMAINDIR = "/mnt/d/FULLDOMAIN/nwmCONUS-v216/"
### STEP 1 END======================================================================




### STEP 2======================================================================
# Current version of NWM is likely to use NLCD 2016 LC, according to (https://www.weather.gov/media/notification/pdf2/scn20-119nwm_v2_1aad.pdf)
nlcdDir <- "/mnt/d/nlcd_2016_land_cover_l48_20210604/nlcd_2016_land_cover_l48_20210604.img"
nlcdObj <- raster(nlcdDir)
methodList <- c("nn", "maj","rawarea","AP")

# Running the following wrapper do tasks as follow:
# 1. Subset full NWM domain files based on  given COMIDs & siteIDs
# 2. If nlcdDir and methodList given, it resamples nlcd into NWM LULC grid using selected methods
for(i in 1:nrow(locs)){
  subset_sequence(locs$comids[i], locs$siteID[i], FULLDOMAINDIR, outDir, nlcdDir, methodList)
}
### STEP 2 END======================================================================
