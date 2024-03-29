% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{outraster}
\alias{outraster}
\title{Create raster objects from simulations}
\usage{
outraster(grdobj, grids = FALSE, emaps = TRUE)
}
\arguments{
\item{grdobj, }{string,  name of list, output of function `grdfile()`}

\item{grids, }{if grids = T  saves simulated maps in 'native' raster package format .grd}

\item{emaps, }{if emaps = T (default), saves e-type and uncertainty maps in format .grd}
}
\value{
\item{simulations}{a rasterstack (package 'raster') where each layer is a simulation}
\item{etype}{a rasterlayer representing median-etype map for disease risk}
\item{uncertainty}{a rasterlayer representing risk uncertainty map for disease risk}
}
\description{
Function read simulation files (.out) returned by `ssdpars()`
and returns a list with simulated maps (rasterstack object),
e-type and uncertainty maps (rasterlayers).
}
\details{
All .grd files are geographic (spatial) data in 'raster' format, and are stored in input folder.
}
