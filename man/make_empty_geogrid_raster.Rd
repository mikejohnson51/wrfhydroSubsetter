% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_empty_geogrid_raster}
\alias{make_empty_geogrid_raster}
\title{Return GeoGrid Raster Structure}
\usage{
make_empty_geogrid_raster(path, var = NULL)
}
\arguments{
\item{path}{a path to a NetCDF file}

\item{var}{an optional variable to extract}
}
\value{
a rasterLayer
}
\description{
Given a NetCDF path, a raster object is returned with the
appropriate structure, valeus are NA, unless a 'var' is defined. S
ee var_names for help indetifing availiable varible names
}
