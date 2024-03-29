% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{ssdpars}
\alias{ssdpars}
\title{Creates a parameters file and generates the simulated maps}
\usage{
ssdpars(
  grdobj,
  maskobj,
  dfobj,
  varmobj,
  simulations = 1,
  nrbias = 20,
  biascor = c(1, 1),
  ndMin = 1,
  ndMax = 32,
  nodMax = 12,
  radius1,
  radius2,
  radius3 = 1,
  ktype = 1
)
}
\arguments{
\item{grdobj, }{string, name of list, output of function `grdfile()`}

\item{maskobj, }{string, name of list, output of function `maskfile()`}

\item{dfobj, }{string, name of list, output of function `irates()`}

\item{varmobj, }{string, name of list, output of function `varmodel()`}

\item{simulations, }{numeric, number of simulations}

\item{nrbias, }{numeric, nr simulations for bias correction}

\item{biascor, }{num vector, flag for (mean, variance) correction (yes = 1, no = 0)}

\item{ndMin, }{numeric, min number of neighbour observations used in  kriging}

\item{ndMax, }{numeric, max number of neighbour observations used in  kriging}

\item{nodMax, }{numeric, max number of previously simulated nodes used in  kriging}

\item{radius1, }{numeric, search radii in the major horizontal axe}

\item{radius2, }{numeric, search radii in the axe orthogonal (horizontal) to radius1}

\item{radius3, }{numeric, search radii in the vertical axe}

\item{ktype, }{numeric, the kriging type to be used (available are: 0 = simple, 1 = ordinary)}
}
\value{
Function returns the following list of objects:
\item{structures}{the number of spatial strutures (excluding nugget-effect)}
\item{parameters}{a dataframe with model information to be passed to ssdpars()}
\item{fittedval}{a numeric vector of fitted values as set by the semi-variogram model}
}
\description{
Creates a parameters file (.par) and invokes dss.c.64.exe to run block simulation program and
generate simulated map files (.out). As input you should provide name of objects
and parameter values required for the simulation process.
}
\details{
Both parameters file (.par) and simulations files (.out) are stored in input folder. Note that the simulation process may take a while, depending mostly on the number of simulation nodes and number of simulations specified.
}
