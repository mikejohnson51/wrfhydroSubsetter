### This example is specific to research article submitted to WRR
### Title: Untangling the impacts of land cover representation and resampling in distributed hydrological model predictions
### Authors: Dong-Hyun (Donny) Kim, J. Michael Johnson, Keith C. Clark, and Hilary K. McMillan
### For any questions about the workflow in this example script, contact Donny Kim via dkim8398@sdsu.edu
# PERT = SHUF-NWM, SHUF = SHUF-AP in the manuscript.

{
  library(dataRetrieval); library(nhdplusTools);
  library(raster); library(sf); library(fasterize);
  library(wrfhydroSubsetter); library(AOI);
  library(tidyverse); library(dplyr);
  #library(resample);
}

# This script is prone to crash during debugging. In such cases, just restart R session. It solves most (if not all) of the problems.
.rs.restartR()


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
## For package development purpose:
#library(devtools)
#detach("package:wrfhydroSubsetter", unload=TRUE); remove.packages("wrfhydroSubsetter");
#devtools::install("/mnt/d/GitHub/wrfhydroSubsetter", repos = NULL, type="source")
#library(wrfhydroSubsetter)
### STEP 2 END======================================================================




### STEP 3======================================================================
# Current version of NWM is likely to use NLCD 2016 LC, according to (https://www.weather.gov/media/notification/pdf2/scn20-119nwm_v2_1aad.pdf)
nlcdDir <- "/mnt/d/nlcd_2016_land_cover_l48_20210604/nlcd_2016_land_cover_l48_20210604.img"
nlcdObj <- raster(nlcdDir)
methodList <- c("nn", "maj","rawarea","AP")

# Running the following wrapper do tasks as follow:
# 1. Subset NWM domain files based on  given COMIDs & siteIDs
# 2. If nlcdDir and methodList given, it resamples nlcd into nwm LULC using selected methods
for(i in 1:nrow(locs)){
  subset_sequence(locs$comids[i], locs$siteID[i], FULLDOMAINDIR, outDir, nlcdDir, methodList)
}
### STEP 3 END======================================================================




### From step 4, it is project-specific code: not for typical use of wrfhydrosubsetter.
### STEP 4======================================================================
# Watershed Shuffle Loop
for(i in 1:nrow(locs)){
  # Specifying WS area
  state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin

  # Using Area Preserve sampled GEOGRID to shuffle LC within watershed boundary
  area_resample_Dir = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo_AP", full = TRUE)
  area_resample = wrfhydroSubsetter::make_empty_geogrid_raster(area_resample_Dir, "LU_INDEX")
  area_resample_shuffle = area_resample # RasterLayer
  area_resample_mask_df = mask_geogrid_byWS(area_resample, basin) #mask_geogrid_byWS is in subset_sequence.R

  # Shuffle LULC inside WS. Cells that are not in WS boundary has value of NA, and omitting them.
  # Then, perturb/shuffle cell values using "sample" function.
  shuffled = transform(na.omit(area_resample_mask_df), layer = sample(layer))
  values(area_resample_shuffle)[as.numeric(rownames(shuffled))] <- shuffled$layer # Overwriting shuffled value over RasterLayer.
  # This overwriting relies on "rownames" which is technically a row index.

  # Setting output file directory
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geo_SHUF.nc")
  file.copy(area_resample_Dir, out_file, overwrite = TRUE)

  #using ncdf4 package to write separate GEOGRID with shuffled LU_INDEX
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX", vals = as.vector(t(apply(as.matrix(area_resample_shuffle), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": Done"))
}


# Watershed NWM PERT Loop
for(i in 1:nrow(locs)){
  # Specifying WS area
  #state = st_transform(AOI::aoi_get(state = "all", county = "all"), 4269)[st_transform(findNLDI(comid = locs$comids[i]), crs = 4269),]
  #name = gsub(" ", "-", tolower(paste(state$name, state$state_name, locs$comids[i], locs$siteID[i], sep = "_")))
  name = namelist[i]
  basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin

  # Using Area Preserve sampled GEOGRID to shuffle LC within watershed boundary
  NWM_geo_Dir = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo_em", full = TRUE)
  NWM_geo = wrfhydroSubsetter::make_empty_geogrid_raster(NWM_geo_Dir, "LU_INDEX")
  NWM_PERT = NWM_geo
  NWM_geo_mask_df = mask_geogrid_byWS(NWM_geo, basin)

  # Shuffle LULC inside WS. Cells that are not in WS boundary has value of NA, and omitting them.
  # Then, perturb/shuffle cell values using "sample" function.
  perturbed = transform(na.omit(NWM_geo_mask_df), layer = sample(layer))
  values(NWM_PERT)[as.numeric(rownames(perturbed))] <- perturbed$layer # Overwriting shuffled value over RasterLayer.
  # This overwriting relies on "rownames" which is technically a row index.

  # Setting output file directory
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geo_PERT.nc")
  file.copy(NWM_geo_Dir, out_file, overwrite = TRUE)

  #using ncdf4 package to write separate GEOGRID with shuffled LU_INDEX
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX", vals = as.vector(t(apply(as.matrix(NWM_PERT), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": Done"))
  rm(NWM_geo_Dir, NWM_geo_mask_df, NWM_geo)
}

### STEP 4 END======================================================================




### STEP 4-1 ======================================================================
# HUC 12 Single LU Loop
for(i in 1:nrow(locs)){
  name = namelist[i]; basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin;

  # Load every geogrid netcdf files
  files = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo", full = TRUE);
  ss = stack()
  for(j in files){
    ss=addLayer(ss, wrfhydroSubsetter::make_empty_geogrid_raster(j, "LU_INDEX"))
  }
  names(ss) = basename(files)

  # crop them to geogrid box
  geogrid_box = spex::qm_rasterToPolygons(ss[[1]]) %>% AOI::bbox_get() %>% sf::st_transform(5070)# %>% sf::st_buffer(40)

  ### Resampling scheme that assigns single LULC values to  each HUC12 ws within study basins
  basin_geo_coord = st_transform(basin, st_crs(ss[[2]])) #ss[[2]] = area preserving geogrid.
  huc12s = get_huc12(AOI = geogrid_box) %>% st_transform(st_crs(basin_geo_coord))
  huc12s = huc12s %>% st_intersection(basin_geo_coord) #only when you want to "clip" huc12 polygons by basin-polygon

  out_geo = ss[[2]]; out_geo_df = as.data.frame(ss[[2]], xy=TRUE);

  # Application of step 4 & 5.
  for (n in 1:nrow(huc12s)){
    #huc12s_df =  mask_geogrid_byWS(out_geo, huc12s[n,]) %>% na.omit() #n stands for each HUC12 ws
    huc12s_df =  raster::mask(out_geo, huc12s[n,]) %>% as.data.frame(xy = TRUE) %>% na.omit()
    huc_dom_LU = huc12s_df %>%  group_by_at(3) %>%  tally() %>% filter (n==max(n)) %>% select_at(1)

    if (nrow(huc_dom_LU)>1) {
      huc12s_df =  mask_geogrid_byWS(out_geo, huc12s[n,]) %>% na.omit() #n stands for each HUC12 ws
      huc_dom_LU = huc12s_df %>%  group_by_at(3) %>%  tally() %>% filter (n==max(n)) %>% select_at(1) %>% as.integer()
    } else {
      huc_dom_LU = huc_dom_LU %>% as.integer() #finding the most dominant LULC
    }

    raster::values(out_geo)[as.numeric(rownames(huc12s_df))] = huc_dom_LU # Updating the raster cell values

  }


  # Setting output file directory
  out_file = paste0('/mnt/d/subsetDOMAINS/', name, "/geo_HUC12L.nc")

  file.copy(paste0('/mnt/d/subsetDOMAINS/', name, "/geo_AP.nc"),
            out_file,
            overwrite = TRUE)

  #using ncdf4 package to write separate GEOGRID with HUC12 ws having single LU representation
  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX", vals = as.vector(t(apply(as.matrix(out_geo), 2, rev))))
  ncdf4::nc_close(nc)
  print(paste0(name, ": HUC12 uniform LC sampling done"))
}


###### Plot
out_geo_mapping = as.data.frame(out_geo, xy=TRUE) %>% rename(LU = 3)
ggplot()+
  #geom_raster(as.data.frame(ss[[2]], xy=T),mapping = aes(x=x, y=y, fill=as.factor(geogrid_area_resample.nc)))+
  geom_raster(out_geo_mapping,mapping = aes(x=x, y=y, fill=as.factor(LU)))+
  #geom_sf(basin_geo_coord, mapping = aes(fill=NULL), color = 'red', alpha = 0)+
  geom_sf(huc12s, mapping = aes(fill=NULL), alpha = 0)#
###### Plot end
### STEP 4-1 END======================================================================




### STEP 5======================================================================
# WS LULC statistics
lookup_table = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class, description = Description, nwm = NWM) %>% dplyr::select(nwm, nlcd, description) %>% arrange(nwm, nlcd)

for(i in 1:nrow(locs)){
  name = namelist[i]; basin = findNLDI(comid = locs$comids[i], find = c("basin"))$basin;

  # Load every geogrid netcdf files
  files = list.files(paste0('/mnt/d/subsetDOMAINS/', name), "geo", full = TRUE)
  ss = stack()
  for(j in files){
    ss=addLayer(ss, wrfhydroSubsetter::make_empty_geogrid_raster(j, "LU_INDEX"))
  }
  names(ss) = basename(files)

  # crop them to geogrid box
  geogrid_box = spex::qm_rasterToPolygons(ss[[1]]) %>%
    AOI::bbox_get() %>%
    sf::st_transform(5070) %>%
    sf::st_buffer(40)
  nlcd_crop = raster::crop(nlcdObj, geogrid_box)

  # Masking (clipping) cropped nlcd into watershed boundary
  # I am not using mask_geogrid_byWS function, as it leads to error.
  basin_coord = st_transform(basin, st_crs(nlcd_crop))
  nlcd_mask_df = raster::mask(nlcd_crop, fasterize::fasterize(basin_coord, nlcd_crop)) %>% getValues() %>% as.data.frame(xy = TRUE)
  nlcd_basin_stat = nlcd_mask_df %>% na.omit() %>% group_by_at(1) %>% tally() %>% rename(nlcd = 1) %>% rename(nlcd_2016 = 2)
  nlcd_basin_stat[,2] = nlcd_basin_stat[,2]/(sum(nlcd_basin_stat[,2])) * 100 # percent

  # I didn't want to spend time to write code that uses map or apply
  mask_LU_stat1 = mask_geogrid_byWS(ss[[1]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[1] := 2)
  mask_LU_stat2 = mask_geogrid_byWS(ss[[2]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[2] := 2)
  mask_LU_stat3 = mask_geogrid_byWS(ss[[3]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[3] := 2)
  mask_LU_stat4 = mask_geogrid_byWS(ss[[4]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[4] := 2)
  mask_LU_stat5 = mask_geogrid_byWS(ss[[5]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[5] := 2)
  mask_LU_stat6 = mask_geogrid_byWS(ss[[6]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[6] := 2)
  mask_LU_stat7 = mask_geogrid_byWS(ss[[7]], basin) %>% na.omit() %>%  group_by_at(3) %>% tally() %>% rename(nwm = 1) %>%  rename(!!basename(files)[7] := 2)

  #merge all data frames in list
  df_list <- list(mask_LU_stat1, mask_LU_stat2,mask_LU_stat3, mask_LU_stat4, mask_LU_stat5, mask_LU_stat6, mask_LU_stat7) %>% reduce(full_join, by='nwm') %>% as.data.frame()
  rm(mask_LU_stat1, mask_LU_stat2,mask_LU_stat3, mask_LU_stat4, mask_LU_stat5, mask_LU_stat6, mask_LU_stat7)

  # make it percent
  for (k in 2:ncol(df_list)){
    df_list[,k] <- df_list[,k]/sum(df_list[,k], na.rm=TRUE) * 100
  }

  # join nlcd and different geogrid LULC stats
  LU_stat_merged = lookup_table %>% right_join(nlcd_basin_stat) %>% right_join(df_list)
  LU_stat_merged[,4:ncol(LU_stat_merged)] = round(LU_stat_merged[,4:ncol(LU_stat_merged)], digits = 1)

  # Let's also look into pearson spatial correlation between different resampling methods
  x = layerStats(ss, stat = 'pearson')$pearson %>% as.data.frame()

  # Save it as csv files
  write_csv(LU_stat_merged, file = paste0(outDir, "/", name, "/LU_WS_stats.csv"))
  write.csv(x, file = paste0(outDir, "/", name, "/Pearson.csv"), row.names= TRUE)
  print(paste0(name, ": Basin stats and Pearson spatial correlation saved as CSV"))
}
### STEP 5 END======================================================================




### STEP 6======================================================================
# Let's use create_wrfinput and create_soilproperties.
for(x in namelist){
  targetDir = paste0(outDir, x)
  files = list.files('/mnt/d/create_all')
  file.copy(list.files('/mnt/d/create_all', full.names = T), targetDir)
  setwd(targetDir)
  system("./create_wrfinput_AP.R --geogrid='geo_AP.nc' --filltyp=3 --laimo=8;
  ./create_wrfinput_SHUF.R --geogrid='geo_SHUF.nc' --filltyp=3 --laimo=8;
  ./create_wrfinput_HUC12L.R --geogrid='geo_HUC12L.nc' --filltyp=3 --laimo=8;
  ./create_wrfinput_PERT.R --geogrid='geo_PERT.nc' --filltyp=3 --laimo=8;
  ./create_wrfinput_LR.R --geogrid='geo_em.nc' --filltyp=3 --laimo=8;
         ")
  file.remove(files);  setwd('/mnt/d/')
}

./create_soilproperties_AP.R; ./create_soilproperties_SHUF.R; ./create_soilproperties_HUC12L.R; ./create_soilproperties_PERT.R; ./create_soilproperties_LR.R

  ### STEP 6 END======================================================================

# Plotting
## Simply checking out.
## https://datacarpentry.org/r-raster-vector-geospatial/02-raster-plot/
#ggplot() +
#  geom_raster(data = area_resample_mask_df, aes(x = x, y = y, fill = as.factor(layer)))+
#  geom_sf(data = basin_coord, colour = "black", fill = "NA")+
#  scale_fill_manual(name = "grp",values = myColors) +
#  coord_sf()


# Checking out if it is shuffled nicely.
{
  s = stack()
  s = addLayer(area_resample, area_resample_shuffle)
  names(s) = c("OG", "Shuffle")

  test = as(basin_coord, 'Spatial') %>% fortify()

  # https://stackoverflow.com/questions/68865682/overlay-polygon-layer-on-a-raster-stack-qplot
  rasterVis::gplot(s) +
    geom_tile(aes(fill = as.factor(value))) +
    geom_path(data=test, aes(long, lat, group=group), color = 'black')+
    facet_wrap(~ variable, ncol = 2) +
    scale_fill_manual(name = "grp",values = myColors) +
    coord_equal()
  #theme(legend.position = "none")
}




### Plotting========================================================================
library(rasterVis)
library(ggplot2)
data = readxl::read_excel('/mnt/d/GitHub/wrfhydroSubsetter/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class,
                description = Description,
                nwm = NWM)

basic_plot_prep = function(x, nlcdObj){
  files = list.files(paste0('/mnt/d/subsetDOMAINS/',namelist[x]), "geo", full = TRUE)
  my_crs = geo_grid_proj(files[1])
  s = stack()

  for(i in files){
    print(i)
    s=addLayer(s, wrfhydroSubsetter::make_empty_geogrid_raster(i, "LU_INDEX"))
  }

  #names(s) = basename(files);   s = s[[1:4]];
  s1 = s[[c(2, 1, 3, 5)]]; names(s1) = c("NWM", "AP", "HUC12L", "SHUF_AP");
  s2 = s[[c(2, 1, 6, 7)]]; names(s2) = c("NWM", "AP", "Maj", "NN");

  geogrid_box = spex::qm_rasterToPolygons(s[[1]]) %>% AOI::bbox_get() %>% sf::st_transform(5070) %>% sf::st_buffer(40)
  #geogrid_box = spex::qm_rasterToPolygons(s[[1]]) %>% AOI::bbox_get() %>% sf::st_buffer(40)

  nlcd_tmp = raster::crop(nlcdObj, geogrid_box) #%>% st_transform(st_crs(s[[1]]))
  #nlcd_tmp = projectRaster(nlcd_tmp, crs= 5070, method = "ngb")
  nlcd_tmp = projectRaster(nlcd_tmp, crs= my_crs, method = "ngb")

  geogrid_box = spex::qm_rasterToPolygons(s[[1]]) %>% AOI::bbox_get() %>% sf::st_transform(my_crs) %>% sf::st_buffer(40)
  nlcd_tmp = raster::crop(nlcd_tmp, geogrid_box)


  ##s=  projectRaster(s, crs= 5070, method = "ngb")
  ##extent(nlcd_tmp) = extent(s[[1]])
  ##s2 = stack(nlcd_tmp, s)


  #tmp = values(s) %>% data.frame() %>% mutate(cell = 1:n()) %>%  tidyr::pivot_longer(-cell)
  #tmp_basin = findNLDI(comid = locs$comids[x], find = c("basin"))$basin %>% st_transform(5070) # 5070 when working with nlcd
  #tmp_HUC12 = get_huc12(AOI = geogrid_box) %>% st_transform(st_crs(5070)) %>% st_intersection(tmp_basin)

  tmp_basin = findNLDI(comid = locs$comids[x], find = c("basin"))$basin %>% st_transform(my_crs) # 5070 when working with nlcd
  tmp_HUC12 = get_huc12(AOI = geogrid_box) %>% st_transform(my_crs) %>% st_intersection(tmp_basin)


  tmp_basin = tmp_basin %>% as('Spatial') %>% fortify();
  tmp_HUC12 = tmp_HUC12 %>% as('Spatial') %>% fortify();


  return(list(s1, tmp_basin, tmp_HUC12, nlcd_tmp, s2))
  #return(list(s, tmp_basin, tmp_HUC12))
}


plot_out = basic_plot_prep(3, nlcdObj);
plot_out[[4]] =  raster::reclassify(plot_out[[4]], rcl = dplyr::select(data, from = nlcd, to = nwm)); #update NLCD classification to NWM

col_lu <- data.frame(
  nlcd.code = c(0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
  nwm.code  = c(0, 16, 23, NA, NA,  1, NA, 19, 11, 14, 15, 22, 8,  7,  20, NA, NA,  5,  3,  18, 17),
  #nwm.code  = c(-.1, 0, 16, 23,  7, 1,  1,  1,  19, 11, 14, 15, 22, 8,  7,  20, NA, NA,  5,  3,  18, 17),

  color = c("#000000",
            "#476BA0", "#D1DDF9",
            "#DDC9C9", "#D89382", "#ff0000", "#AA0000",
            "#B2ADA3",
            "#68AA63", "#1C6330", "#B5C98E",
            "#A58C30", "#CCBA7C",
            "#E2E2C1", "#C9C977", "#99C147", "#77AD93",
            "#DBD83D", "#AA7028",
            "#BAD8EA", "#70A3BA") ,

  name = c(NA,
           "Open Water", "Ice/Snow",
           "Developed (Open)", "Developed (Low)", 'Developed (Medium)', 'Developed (High)',
           "Barren",
           "Deciduous Forest", "Evergreen Forest", "Mixed Forest",
           "Dwarf Scrub", "Shurb",
           "Grassland", "Sedge",'Lichens', "Moss",
           "Pasture", "Culitivated Crops",
           "Woody Wetlands", "Herbaceous Wetlands"),

  stringsAsFactors = FALSE)

tt = dplyr::left_join(data, col_lu, by = c("nlcd" = "nlcd.code")) %>% arrange(nwm) %>% na.omit()
#tt = dplyr::left_join(data, col_lu, by = c("nwm" = "nwm.code")) %>% arrange(nwm) %>% na.omit()
myColors <- tt$color
names(myColors) <- levels(tt$nwm)
colScale <- scale_colour_manual(name = "grp",values = myColors, breaks=tt$nwm)

plot_description = tt$name
plot_description[1] = "Developed"


# p1 is resampled LULC grids
p1 = rasterVis::gplot(plot_out[[1]]) +
  geom_tile(aes(fill = as.factor(value))) +
  facet_wrap(~ variable, ncol = 2) +
  geom_path(data=plot_out[[3]], aes(long, lat, group=group), color = 'black')+
  geom_path(data=plot_out[[2]], aes(long, lat, group=group), color = 'blue')+
  scale_fill_manual(name = "LULC",values = myColors, breaks=tt$nwm, na.value=NA) +
  colScale+
  coord_equal()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x=element_blank()
    ,axis.text.y=element_blank()
    ,legend.position = "none"
    ,axis.title.x=element_blank()
    ,axis.title.y=element_blank()
    ,axis.ticks = element_blank()
    ,plot.margin = margin(2, 2, 2, 2)
    )

# p2 = NLC
p2 = rasterVis::gplot(plot_out[[4]])+
  geom_tile(aes(fill = as.factor(value))) +
  geom_path(data=plot_out[[3]], aes(long, lat, group=group), color = 'black')+
  geom_path(data=plot_out[[2]], aes(long, lat, group=group), color = 'blue')+
  scale_fill_manual(name = "Noah-MP\nLULC",values = myColors, breaks=tt$nwm,
                    labels = plot_description) +
  colScale+
  coord_equal()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x=element_blank()
    ,axis.text.y=element_blank()
    ,axis.title.x=element_blank()
    ,axis.title.y=element_blank()
    ,axis.ticks = element_blank()
    ,legend.position = "none"
    , legend.spacing.x = unit(10, 'pt')
    ,plot.margin = margin(2, 2, 2, 2)
    )

p_grid = plot_grid(p2, p1, ncol=2, greedy=T, rel_widths = c(.45, .55))
legend_b <- get_legend(p2 +
                         guides(color = guide_legend(nrow = 1), label.hjust=0) +
                         theme(legend.position = "bottom")
                       )
p_grid2 = plot_grid(p_grid, legend_b, rel_heights = c(1.2, .3), ncol=1, greedy=T)
p_grid2

#gridExtra::grid.arrange(p2, p1, ncol=2)

p3 = rasterVis::gplot(plot_out[[1]]) +
  geom_tile(aes(fill = as.factor(value))) +
  facet_wrap(~ variable, ncol = 4) +
  geom_path(data=plot_out[[3]], aes(long, lat, group=group), color = 'black')+
  geom_path(data=plot_out[[2]], aes(long, lat, group=group), color = 'blue')+
  scale_fill_manual(name = "LULC",values = myColors, breaks=tt$nwm, na.value=NA) +
  colScale+
  coord_equal()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x=element_blank()
    ,axis.text.y=element_blank()
    ,legend.position = "none"
    ,axis.title.x=element_blank()
    ,axis.title.y=element_blank()
    ,axis.ticks = element_blank()
    ,plot.margin = margin(2, 2, 2, 2)
  )

p_grid_fig4 = plot_grid(p3, legend_b, ncol=1, rel_heights = c(.9, .3))
p_grid_fig4

ggplot(data = tmp) +
  geom_histogram(aes(x = as.factor(value), fill= name), stat = "count", position='dodge') +
  #ggpubr::fill_palette('aaas') +
  theme_linedraw() +
  theme(legend.position = "bottom")

hist(s, breaks = 20)
### test plot end========================================================================

