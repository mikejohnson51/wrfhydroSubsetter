library(dataRetrieval)
library(raster)
library(sf)
library(resample)
library(wrfhydroSubsetter)
library(dplyr)

locs = data.frame(comids = c(23762661, 1631587, 5894384, 19389766, 191739, 5781369),
           siteID = c('14190500', '14190500', '01616500', '03118500' ,'6709000' , '08159000'))

for(i in 3:nrow(locs)){
  create_lu_subs(comid = locs$comids[i],
                 siteID = locs$siteID[i],
                 outDir = "/Volumes/Transcend/land-cover-tests/",
                 nlcd = raster('/Volumes/Backup/NLCD/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img'),
                 method = c("area", "rawarea")
                 )
  message(i, " of ", nrow(locs))
}



create_lu_subs = function(comid, siteID, outDir, nlcd, method){

state = AOI::aoi_get(state = "all", county = "all")[st_transform(findNLDI(comid = comid), 4269),]
name = gsub(" ", "-", tolower(paste(state$name, state$state_name, comid, siteID, sep = "_")))
subset_files = paste0(outDir, "/", name , "/")

basin = findNLDI(comid = comid, find = c("basin"))$basin
##mapview::mapview(basin)
subset_wrf_hydro_domain(AOI = basin,  domain_files = download_conus_nwm("/Volumes/Transcend/"),  outDir = subset_files)

geo  = list.files(subset_files, "geo_em.d0x", full.names = TRUE)

output = spex::qm_rasterToPolygons(wrfhydroSubsetter::make_empty_geogrid_raster(geo))

data = readxl::read_excel('/Users/mikejohnson/Documents/nwm_land_use_mappings.xlsx') %>%
  dplyr::select(nlcd = Class,
                name = Name,
                description = Descrip,
                nwm = NWMv2.1)

o2 = output %>%
  AOI::bbox_get() %>%
  sf::st_transform(5070) %>%
  sf::st_buffer(40)

nlcdCrop = raster::crop(nlcd, o2)
nlcd_nwm = raster::reclassify(nlcdCrop, rcl = dplyr::select(data,from = nlcd, to = nwm ))

for(j in 1:length(method)){
  new_lu = resample::resampleData(input = nlcd_nwm, output = output, method = method[j])
  out_file = paste0(subset_files, "geogrid_", method[j], "_resample.nc")
  file.copy(geo, out_file, overwrite = TRUE)

  nc = ncdf4::nc_open(out_file, write = TRUE)
  ncdf4::ncvar_put(nc, "LU_INDEX",
                   vals = as.vector(t(apply(as.matrix(new_lu), 2, rev))))
  ncdf4::nc_close(nc)
}

}


files = list.files('/Volumes/Transcend/land-cover-tests/travis_texas_5781369_08159000', "geo", full = TRUE)
s = stack()
for(i in files){
  print(i)
  s=addLayer(s, wrfhydroSubsetter::make_empty_geogrid_raster(i, "LU_INDEX"))
}
names(s) = basename(files)
plot(s)

tmp = values(s) %>% data.frame() %>% mutate(cell = 1:n()) %>%  tidyr::pivot_longer(-cell)


gplot(s) +
  geom_tile(aes(fill = as.factor(value))) +
  facet_wrap(~ variable, ncol = 1) +
  scale_fill_manual(name = "grp",values = myColors) +
  coord_equal()+
  theme(legend.position = "none") +
ggplot(data = tmp) +
  geom_histogram(aes(x = as.factor(value), fill= name), stat = "count", position='dodge') +
  ggpubr::fill_palette('aaas') +
  theme_linedraw() +
  theme(legend.position = "bottom")

hist(s)


myColors <- tt$color
names(myColors) <- levels(tt$nwm)
colScale <- scale_colour_manual(name = "grp",values = myColors)

tt = left_join(data, col_lu, by = c("nlcd" = "nlcd.code"))

mapview::mapview(s)

col_lu <- data.frame(
  nlcd.code = c(-.1, 0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 72, 73, 74, 81, 82, 90, 95),
  nwm.code  = c(-.1, 0, 16, 23, NA, 1,  NA, NA, 19, 11, 14, 15, 22, 8,  7,  20, NA, NA, 2,  3,  18, 17),

  color = c("#000000",
            "#476BA0", "#D1DDF9",
            "#DDC9C9", "#D89382", "#ED0000", "#AA0000",
            "#B2ADA3",
            "#68AA63", "#1C6330", "#B5C98E",
            "#A58C30", "#CCBA7C",
            "#E2E2C1", "#C9C977", "#99C147", "#77AD93",
            "#DBD83D", "#AA7028",
            "#BAD8EA", "#70A3BA", NA) ,

  name = c(NA, "EMPTY", "Open Water", "Ice/Snow", "Developed (Open)", "Developed (Low)", 'Developed (Medium)', 'Developed (High)', "Barren",
           "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Dwarf Scrub", "Shurb", "Grassland", "Sedge", 'Lichens', "Moss",
           "Pasture", "Culitivated Crops", "Woody Wetlands", "Herbaceous Wetlands"),

  stringsAsFactors = FALSE)

