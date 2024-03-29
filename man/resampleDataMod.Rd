% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resampleDataMod.R
\name{resampleDataMod}
\alias{resampleDataMod}
\title{Resample_mod}
\usage{
resampleDataMod(
  input,
  output_geo = NULL,
  cellsize = 1000,
  method = "area",
  no_data = NA,
  seed = 10291991
)
}
\arguments{
\item{input}{a raster of fine resolution}

\item{output_geo}{*a coarse resolution of target raster. Customized for wrf-hydro subsetter}

\item{cellsize}{cell resolution of output target grid. In units of input CRS}

\item{method}{the method for resampling, nearest neighbor (nn), majority rule (maj), Area-based (area), or raw area (rawarea)}

\item{no_data}{a value to be treated as NO_DATA. Default to NA}

\item{seed}{seed number for random class assignment in rawarea}
}
\value{
a resampled raster
}
\description{
resample an input raster to an output grid using three different methods
}
