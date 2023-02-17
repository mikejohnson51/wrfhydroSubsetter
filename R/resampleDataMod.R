# These two functions originate from "resample" R package (https://github.com/mikejohnson51/resample).
# Minor edits were made by Donny Kim to accomodate the need of wrfhydroSubsetter.



#' @title Build output grid
#' @description Builds output grid
#' @param input Polygon or WRF-Hydro geogrid
#' @param cellsize resolution
#' @importFrom magrittr %>%
#' @importFrom sf st_bbox st_make_grid st_sf

output_grid = function(input, cellsize){
  bb     = st_bbox(input)
  cols   = seq(bb$xmin,  bb$xmax, cellsize)
  rows   = seq(bb$ymax,  bb$ymin, -cellsize)
  cls    = round((bb$xmax - bb$xmin) / cellsize)
  rws    = round((bb$ymax - bb$ymin) / cellsize)
  
  output = st_make_grid(input,
                        n        = c(cls, rws),
                        offset   = c(min(cols), min(rows)),
                        cellsize = c(cellsize, cellsize),
                        crs      = st_crs(input),
                        what     = "polygons") %>%
    st_sf()
  
  output$row_id = rep(c(1:rws), cls)
  output$col_id = rep(c(1:cls), each = rws)
  output
}




#' @title Summarize Output Grid Cells
#' @description Summarizies the percentage of each output grid cell that is covered by each disiticnt class in the input grid
#' @param input an input raster to resample
#' @param output the output grid to resample to. `sf` object
#' @param no_data the value representing "NO DATA", default is NA
#' @importFrom dplyr group_by summarize mutate select n rename ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom tabularaster cellnumbers as_tibble
#' @importFrom magrittr %>%
#' @return a tibble

summarize.cells = function(input, output, no_data = NA){
  require(resample)
  require(tabularaster)
  
  # Extract Raster as tibble, 1 column for cellvalues, 1 for cellindex
  values = as_tibble(input)
  
  # Determine the percentage of each category in the input
  # Identify the corresponding number of cells in the output
  count = nrow(output) * (table(values$cellvalue) / nrow(values))
  
  # Flood the needed output values
  output_count = floor(count)
  
  # Identify how many cells can be added to the floored values and add to those with the
  # largest decimal value
  index <- order((count - output_count), decreasing = TRUE)[1:(nrow(output) - sum(output_count))]
  output_count[index] = output_count[index] + 1
  
  # Sort and remove classes with 0 needed cells
  output_count = sort(output_count)
  output_count = output_count[output_count != 0]
  
  # Reprioritize  no_data value if provided
  output_count = if(is.na(no_data)){
    output_count
  } else {
    c(output_count[names(output_count) == no_data], output_count[names(output_count) != no_data] )
  }
  
  # Determine which input cells underlay each output object
  cn <- tabularaster::cellnumbers(raster(input), output$geometry) #DK: this is the point where it doesn't work!
  
  # Add new column matching the cell values form the input
  cn$cellvalue_ = values$cellvalue[match(cn$cell_, values$cellindex)]
  
  # Group by Output object and cell value, count, and determine percentage
  # Pivot into wide table and rename object to id
  areal_per = cn %>%
    dplyr::group_by(object_, cellvalue_) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::mutate(v.pct = count / sum(count)) %>%
    tidyr::pivot_wider(id_cols = object_, values_from = v.pct, names_from = cellvalue_) %>%
    dplyr::rename(id = object_ ) %>%
    dplyr::ungroup()

  
  #Clean up and return
  areal_per = areal_per[names(areal_per) %in% names(output_count)]
  
  return( list(areal_per = areal_per, output_count = output_count) )
}




#' @title Resample_mod
#' @description resample an input raster to an output grid using three different methods
#' @param input a raster of fine resolution
#' @param output_geo *a coarse resolution of target raster. Customized for wrf-hydro subsetter
#' @param cellsize cell resolution of output target grid. In units of input CRS
#' @param method the method for resampling, nearest neighbor (nn), majority rule (maj), Area-based (area), or raw area (rawarea)
#' @param no_data a value to be treated as NO_DATA. Default to NA
#' @param seed seed number for random class assignment in rawarea
#' @importFrom raster writeRaster raster extent crs
#' @importFrom sf st_bbox st_centroid st_crs st_make_grid st_coordinates st_sf
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @return a resampled raster
#' @export

resampleDataMod = function(input, output_geo=NULL, cellsize = 1000,  method = "area", no_data = NA, seed = 10291991){
  
  require(resample)
  
 # Added by DK
 if (!is.null(output_geo)) {
   message("Using geogrid as target raster")
   bb = extent(output_geo)
   extent(input) = bb
   #crs(input) <- paste0(st_crs(output_geo))[1]
   output = output_grid(output_geo, cellsize = cellsize)
 } else {
   output = output_grid(input, cellsize = cellsize)
 }
  
  # Below is original
  if (method == "nn") {
    
    message("Nearest Neighbor...")
    rr = call_gdal_mod(input, "near", cellsize = cellsize)
    
  } else if(method == "rawarea"){
    
    message("Raw Areal Proportions...")
    #output = output_grid(input, cellsize = cellsize)
    out    = summarize.cells(input, output, no_data = no_data)
    vals   = unlist(lapply(1:length(out$output_count),
                           function(x){ as.numeric(rep(names(out$output_count)[x], out$output_count[x])) }
    ))
    set.seed(seed)
    output$cat = sample(vals)
    output = output%>% 
      sf::st_transform(st_crs(output_geo))
    rr = resample::raster.from.vector(output)

    
  } else if(method == "AP"){
    
    message("Areal Proportions...")
    #output     = output_grid(input, cellsize = cellsize)
    out        = summarize.cells(input,
                                 output = output,
                                 no_data = no_data)
    output$cat = resample::reclass(areal_per = out$areal_per,
                         output_count = out$output_count,
                         no_data = no_data)
    #crs.here = paste0(sf::st_crs(output_geo)[1])
    #output = output%>% 
    #  sf::st_transform(crs.here)
    rr         = resample::raster.from.vector(output)

  } else if (method == 'maj') {
    
    message("Majority Rule...")
    rr = call_gdal_mod(input, "mode", cellsize)
    
  } else{
    
    message("RNA...")
    #output = output_grid(input, cellsize = cellsize)
    output$cat = resample::rna(input, output)
    rr = resample::raster.from.vector(output)
  }
  
  #rr
  crs(rr) = paste0(st_crs(output_geo))[1]
  rr
}




#' @title Call GDAL
#' @param input input grid (raster)
#' @param method resampling method
#' @param cellsize output cell diminsion
#' @keywords internal
#' @return RasterLayer
#' @importFrom raster writeRaster raster
#' @importFrom sf gdal_utils

call_gdal_mod = function(input, method, cellsize){
  #require(resample)
  
  tmp = tempfile(pattern = "input", fileext = '.tif')
  nn  = tempfile(pattern = "nn",    fileext = '.tif')
  
  raster::writeRaster(input, filename = tmp, overwrite = TRUE)
  
  sf::gdal_utils("translate",
                 source      = tmp,
                 destination = nn,
                 options = c("-r", method, "-tr", c(cellsize,cellsize)))
  
  rr = raster(nn)
  
  file.remove(tmp)
  
  rr
  
}