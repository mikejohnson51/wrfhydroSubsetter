% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetter.R
\name{subset_wrf_hydro_domain}
\alias{subset_wrf_hydro_domain}
\title{WRF-Hydro Subsetter}
\usage{
subset_wrf_hydro_domain(AOI, domain_files, outDir, config)
}
\arguments{
\item{AOI}{sf POLYGON. The area to cut out}

\item{domain_files}{character. A path to the local domain files}

\item{outDir}{Where to put the cutout data should be a named directory!}
}
\value{
files on disk
}
\description{
Create a cutout of the NWM domain files
}
