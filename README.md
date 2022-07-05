# EpiGeostats

EpiGeostats R package allows to fit a geostatistical model based on block direct sequential simulation, to obtain disease risk and assess associated risk uncertainty maps with high spatial resolution, and also to generate a single-map summarizing effectively these two key elements of disease risk mapping. 

It is designed for disease risk mapping, but may also be used for similar problems such as mortality risk, or applied to other fields like ecology or criminology. Details on the geostatistical algorithm used in EpiGeostats are described in Azevedo et al (https://doi.org/10.1186/s12942-020-00221-5).

# R wrapper

EpiGeostats R package is a wrapper for running dss.c.64.exe - a software tool for running geostatistical algorithms - and pixelate, an R package for creating geostatistical maps with varying pixel sizes to represent disease mapping and spatial uncertainty in a single-map.

The package provides convenient wrappers to set parameters required for geostatistical simulations and write some files in dss.c.64.exe readable format, after which block direct sequential simulation algorithm can be executed. Results generate a set of simulated disease risk maps, a median disease risk map, spatial uncertainty risk map, and a pixelated map version representing both disease risks and spatial uncertainty risks of disease in a single-map, that can be extracted and/or plotted in R 

# Getting started

To run the EpiGeostats R package for the first time, you must run three steps with the following code:

1. if not installed, installs `Rtools` software and `devtools` package:  
```r
# install Rtools from CRAN
if (!require("installr")){
    install.packages("installr")
    }

installr::install.Rtools()

# install devtools from CRAN
if (!require("devtools")){
    install.packages("devtools")
    }
``` 

2. downloads dss.c.64.exe from GitHub:
```r
# create folder input
inp <- "./input"
if(!file.exists(inp)) dir.create(inp, recursive = TRUE)

# download dss.c.64.exe.zip from github
gitURL <- "https://github.com/maluicr/dss/raw/main/DSS.C.64.exe.zip"
utils::download.file(url = gitURL, destfile = file.path(inp, "DSS.C.64.exe.zip"))

# unzip file
unzip(file.path(inp, "DSS.C.64.exe.zip"), exdir = inp)
``` 

3. installs required packages from GitHub (please accept any suggested package updates): 

```r
  # install packages from GitHub
  devtools::install_github("maluicr/EpiGeostats", build_vignettes = TRUE, dependencies = TRUE, upgrade = "always")
  devtools::install_github("aimeertaylor/pixelate", build_vignettes = TRUE, dependencies = TRUE)
```

# Example

As an example, EpiGeostats is used to map COVID-19 incidence in Portugal on 15 January 2021. After packages installation and dss.c.64.exe download, run the code below to load `EpiGeostats` package and to follow the example, as presented in vignette:

```r
# load EpiGeostats 
library(EpiGeostats)

# see vignette document
help(package = "EpiGeostats")
```
