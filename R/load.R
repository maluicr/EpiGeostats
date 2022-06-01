#' Creates an disease incidence rates file to be read by dss.64.c.exe
#'
#' Function computes disease incidence rates and variance-error terms by region,
#' and writes the result into a text file (.out) that is stored in `./input` folder.
#' As input you should provide a dataframe with disease cases, population size
#' and cartesian coordinates by region. Columns should include the following data:
#' id of region, x, y and z cartesian coordinates at regions mass center,
#' number of disease cases and population size by region.
#
#' @param dfobj string, dataframe name with disease data
#' @param oid character, fieldname for region id
#' @param xx character, fieldname for x-coordinates
#' @param yy character, fieldname for y-coordinates
#' @param zz character, fieldname for z-coordinates
#' @param cases character, fieldname for number of cases
#' @param pop character, fieldname for population size
#' @param casesNA, numeric, an integer used to replace rows with cases = NA
#' @param day character, string indicating the date (format 'yyyymmdd') of disease cases
#' @param perhab numeric, an integer indicating how much the rate is multiplied for, to express the disease rate (e.g 10000 or 100000)

#' @return The following list of objects:
#' \item{rates}{dataframe; results containing id region, x,y,z coordinates, incidence rate, variance-error term and population by region}
#' \item{mrisk}{numeric; estimated global risk (/`perhab`)}
#' \item{file}{list characters; indicating the day of data collection, filename and root folder where text file is stored }
#' \item{ssdpars}{list numeric; list of values to be passed to `ssdpars()`}
#'
#' @importFrom stats dist quantile sd
#' @importFrom utils head install.packages read.table write.table
#'
#'
#' @details The function also writes and stores a text file (.out). This file contains x, y, z and incidence rate by row (region) using GeoEAS file format.

#' @export
irates = function(dfobj = NA, oid = NA, xx = NA, yy = NA, zz = NA,
                  cases = NA, pop = NA, casesNA = 2, day = NA, perhab = 10^5) {

  # rate per phab habitants
  phab = perhab

  # index variables id, x, y, z, cases, risk pop
  ioid = grep(oid, colnames(dfobj))
  ix = grep(xx, colnames(dfobj))
  iy = grep(yy, colnames(dfobj))
  iz = grep(zz, colnames(dfobj))
  ic = grep(cases, colnames(dfobj))
  ip = grep(pop, colnames(dfobj))

  # compute overall mean rate (exclude rows w/ cases = NA)
  dfobjnas = base::subset(dfobj, dfobj[, ic]!= "NA")
  n = ncol(dfobjnas)
  dfobjnas$rate = phab * dfobjnas[, ic] / dfobjnas[, ip]
  poptnas = base::sum(dfobjnas[, ip])
  m = base::sum(dfobjnas[, "rate"] * dfobjnas[, ip]) / poptnas

  # error variance term (m/n_i)
  error = m / dfobj[, ip]

  # NA cases set to casesNA
  dfobj[, ic] = base::ifelse(is.na(dfobj[, ic]), casesNA, dfobj[, ic])

  # recalculate crude rates
  rate = phab * dfobj[, ic] / dfobj[, ip]

  tab = data.frame (dfobj[, ioid], x = dfobj[, ix], y = dfobj[, iy], z = dfobj[, iz], rate, error)

  # cases file for dss

  # set folder
  foldin = "input"

  # create folder for inputs
  if(!file.exists(foldin)) dir.create(foldin, recursive = FALSE)

  # store string with path for input files
  wkin = paste0(getwd(), "/", foldin)

  # prepare data to write file
  tabnotf = tab[, c("x", "y", "z", "rate")]

  # store nr of variables
  nvars = base::ncol(tabnotf)

  # store nr of observations
  nobs = base::nrow(tabnotf)

  # store vector variable names
  namevars = names(tabnotf)

  # create file path
  not_name = "rate"
  not_nameO = paste0(not_name, ".out")
  fnot = paste0(wkin, "/", day, "_", not_nameO)

  if (file.exists(fnot)){
    file.remove(fnot)
  }

  # create notification file
  file.create(fnot)

  # write file header
  cat(not_name, file = fnot, sep="\n")
  cat(nvars, file = fnot, sep="\n", append = TRUE)
  cat(namevars, file = fnot, sep="\n", append = TRUE)

  # write notification data
  utils::write.table(format(tabnotf, digits = NULL, justify = "right"),
              file = fnot, quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
  if(file.exists(fnot)) {
    cat(paste0(fnot, " created."))
  }

  # pars for ssdr.par
  xcol = grep("x", colnames(tabnotf))
  ycol = grep("y", colnames(tabnotf))
  zcol = grep("z", colnames(tabnotf))
  varcol = grep("rate", colnames(tabnotf))
  minval = min(tabnotf$rate)
  maxval = max(tabnotf$rate)

  # return data.frame object for semi- calcs
  tabvgm = data.frame (id = dfobj[, ioid], x = dfobj[, ix], y = dfobj[, iy],
                       z = dfobj[, iz], rate, err = error, pop = dfobj[, ip])
  listf = list(day = day, name = paste0(day, "_",not_nameO), folder = wkin)
  listpars = list(nvars = nvars, xcolumn = xcol, ycolumn = ycol, zcolumn = zcol,
                  varcol = varcol, minval = minval, maxval = maxval)
  return(list(rates = tabvgm, mrisk = m, file = listf, ssdrpars = listpars))
}

#' Function creates a grid file to be read by dss.64.c.exe
#'
#' Converts grid nodes (spatial grid) into a text file (.out) that is stored in `input` folder.
#' As input you should provide the spatial grid (SpatialPixelsDataFrame), with id region values at simulation grid nodes,
#' and a list returned by funtion irates().
#'
#' @param rateobj, string, name of list, output of function irates()
#' @param gridimage, string,  name of block file, a SpatialPixelsDataFrame object
#' @param NAval, numeric, integer with grid value for "No data"
#'
#' @return
#' \item{gridpars}{list parameters from block data: number rows and columns, resolution, (x,y), coordinates of the origin, NA value}
#' \item{outgrid}{list vector node values, vector id of regions, and number of regions   }
#' \item{file}{list characters; indicating the day of data collection, filename and root folder where file is stored}
#' \item{ingrid}{the SpatialPixelsDataFrame object used}
#'
#' @details The gridded values should refer to the region id's at simulation locations (nodes). Every regions recorded in disease data file should be assigned to 1 or more grid node.
#'
#'
#' @importFrom raster res extent as.matrix
#'
#' @export
grdfile = function (rateobj, gridimage, NAval = -999){

  # strings w/ path to store files
  day = as.character(rateobj[["file"]]["day"])
  folder = as.character(rateobj[["file"]]["folder"])

  grd = raster::raster(gridimage)

  # grid parameters
  ny = grd@nrows
  nx = grd@ncols
  resx = raster::res(grd)[1]
  resy = raster::res(grd)[2]
  ox = raster::extent(grd)[1]
  oy = raster::extent(grd)[3]

  # create matrix to store block coordinates
  matA = raster::as.matrix( grd, nrow = ny, ncol = nx )
  matA2 = matrix(NAval, nrow = ny, ncol = nx)

  # fill matrix with block coordinates
  for (i in 1:ny){
    matA2[i,] <- c(matA[ny-i+1,] )
  }

  # convert to vector str
  matA3 = as.data.frame(t(matA2))
  stac3 = raster::stack(matA3)
  stacf = stac3$values

  # set NAs values
  nas = NAval
  stacf[is.na(stacf)] = nas

  # create array for grdfile
  grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = FALSE)
  grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, by = resx), y = seq(oy, oy + (ny-1) * resy, resy))
  grx = grxy[,1]
  gry = grxy[,2]
  gridout = array (c (grdata, grx, gry), dim =c(ny, nx , 3))
  gridout = round(gridout, 4)

  # write blocksfile for dss

  # store vector grid id
  massid = raster::unique(grd)

  # store length vector
  massn = length(massid)

  # create file path
  blk_name = "grid"
  blk_nameO = paste0(blk_name,".out")
  fblk = paste0(folder, "/", day, "_", blk_nameO)

  if (file.exists(fblk)){
    file.remove(fblk)
  }

  # create notification file
  file.create(fblk)

  # write file header
  cat(blk_name, file = fblk, sep="\n")
  cat(massn, file = fblk, sep="\n", append = TRUE)

  # write block id, rate, error and x, y, z
  for (i in 1 : massn){
    id = massid[i]
    # store nr blocks
    nc = sum(gridout[,,1] == id)

    # get vectors from irates function
    oid = rateobj[["rates"]][, "id"]
    rate = rateobj[["rates"]][, "rate"]
    err = rateobj[["rates"]][, "err"]
    dftab = data.frame(oid, rate, err)

    # store rate and error
    dft = dftab[dftab$oid == id, c("rate","err")]
    rt = as.numeric(dft[1])
    er = as.numeric(dft[2])

    # write block id
    cat(paste0("Block#", id), file = fblk, sep="\n", append = TRUE)

    # write rate, error & nr blocks
    cat(rt, file = fblk, sep="\n", append = TRUE)
    cat(er, file = fblk, sep="\n", append = TRUE)
    cat(nc, file = fblk, sep="\n", append = TRUE)

    # data.table::fwrite
    blk_list = list()
    for (k in 1:ny) {
      for(l in 1:nx){
        if(gridout[k, l, 1] == id) {
          d = paste(gridout[k, l, 2], "\t", gridout[k, l, 3], "\t",  0)
          blk_list[[length(blk_list) + 1]] = list (d)
        }
      }
    }
    data.table::fwrite(blk_list, file = fblk, append = TRUE, sep="\n")
  }
  if(file.exists(fblk)) {
    cat(paste0(fblk, " created."))
  }

  listgrid = grd
  listgridpars = list(nodes = c(nx, ny), resolution = c(resx, resy), origin = c(ox, oy), NAs = NAval)
  listfile = list(day = day, name = paste0(day, "_", blk_nameO), folder = folder)
  listgridout = list(values = stacf, idblock = massid, nblock = massn)
  return(list(gridpars = listgridpars, outgrid = listgridout , file = listfile, ingrid = listgrid))
}

#' Function creates a mask file to be read by dss.64.c.exe
#'
#' Creates a mask file and writes the result into a text file (.out) that is stored in `input` folder.
#' The text file (.out) stores values {-1,0} where -1 are assigned to nodata locations and 0 are assigned to nodes with values (id region).
#' As input you should provide the name of list returned by `grdfile()`.
#'
#' @param grdobj, string, name of list, output of function `grdfile()`
#'
#' @return `maskfile()` also returns the following list of objects:
#' \item{file}{list characters; indicating filename and root folder where text file is stored}
#' \item{zones}{list list of values to be passed to ssdpars()}
#'
#' @export
maskfile = function(grdobj){

  obj = unlist(grdobj[["outgrid"]]["values"], use.names = FALSE)
  na = unlist(grdobj[["gridpars"]]["NAs"], use.names = FALSE)
  day = as.character(grdobj[["file"]]["day"])
  folder = as.character(grdobj[["file"]]["folder"])

  # write maskfile for dss

  # create mask vector
  mask = ifelse(obj == na, -1, 0)
  val = unique(mask)
  mask_zones = length(unique(mask))

  # prepare data to write file
  nvars = 1
  namevars = "values"
  nval = length(mask)

  # create file path
  msk_name = "mask"
  msk_nameO = paste0(msk_name, ".out")
  fmsk = paste0(folder, "/", day,"_", msk_nameO)

  if (file.exists(fmsk)){
    file.remove(fmsk)
  }

  # create mask file
  file.create(fmsk)

  # write file header
  cat(msk_name, file = fmsk, sep="\n")
  cat(nvars, file = fmsk, sep="\n", append = TRUE)
  cat(namevars, file = fmsk, sep="\n", append = TRUE)

  # write mask data
  write.table(mask, file = fmsk, append = TRUE, row.names = FALSE, col.names = FALSE)
  if(file.exists(fmsk)) {
    cat(paste0(fmsk, " created."))
    }
  listf = list(day = day, name = paste0(day, "_", msk_nameO), folder = folder)
  listz = list(nzones = mask_zones, zoneval = val)

  return(list(file = listf, zones = listz))
}

#' Calculates experimental population-weighted semi-variogram
#'
#' Calculates experimental population-weighted  semi-variograms from disease incidence rates.
#' For now is only implemented in omnidirectional case.
#'
#' @param dfobj, string, name of list, output of function irates()
#' @param lag, numeric, the lag distance used for semi-variogram estimates
#' @param nlags, numeric, the number of lags to calculate semi-variogram
#'
#' @return The function returns a list with the variance and semi-variogram estimates weighted by population size at nlags:
#' \item{weightsvar}{is the value of the weighted population variance}
#' \item{semivar}{a dataframe with a vector of distances, a vector of experimental semivariogram values and number of pairs}
#'
#' @export
varexp = function(dfobj, lag, nlags){

  # store nr of observations
  nobs = nrow(dfobj[["rates"]])

  # experimental semi-variogram : get data
  rate = dfobj[["rates"]][,"rate"]
  ratexy = cbind(dfobj[["rates"]][c("x","y")])
  pop = dfobj[["rates"]][,"pop"]
  # store pop vector as double
  pop = as.double(pop)
  m = dfobj[["mrisk"]]


  # weighted sample variance for sill estimation
  # no ref about this sill estimation
  # looked in goovaerts & cressie books, journel.
  # at end, used an adapted fortran code from mjp.

  # create integer num & den to calc weighted sample var
  nwsv = 0
  dwsv = 0

  # calc num and denominator
  for (i in 1:nobs){
    nwsv = nwsv + pop[i]^2 / (2 * pop[i]) * (rate[i] - m)^2
    dwsv = dwsv + pop[i]^2 / (2 * pop[i])
  }

  # weighted sample variance
  wsvar = nwsv / dwsv

  # experimental semi-variogram : calc distances
  # lag distance
  lagd = lag
  # cut distance
  lagend = lag * nlags
  # nr lags
  lagn = nlags + 1
  # vector of lags
  lags = c()
  for (i in 1:lagn) {
    lags[i] = lagd * i - lagd
  }

  # compute distance matrix
  matd = as.matrix(dist(ratexy))

  # experimental semi-variogram: create vectors to store results
  # n pairs dist(h)
  nh = vector( mode = "integer", length = (lagn-1))
  # total dist(h)
  dh = vector( mode = "numeric", length = (lagn-1))
  # mean dist(h)
  mh = vector( mode = "numeric", length = (lagn-1))
  # numerator gamma(h)
  num_gh = vector( mode = "double", length = (lagn-1))
  # denominator gamma(h)
  den_gh = vector( mode = "double", length = (lagn-1))
  # gamma(h)
  gammah = vector( mode = "double", length = (lagn-1))

  # experimental semi-variogram: compute
  for (k in 1:(lagn-1)) {
    for (i in 1:nobs) {
      for (j in 1:nobs) {
        if(matd[i,j] > lags[k] & matd[i,j] <= lags[k+1]){
          nh[k] = nh[k] + 1
          dh[k] = dh[k] + matd[i,j]
          num_gh[k] = num_gh[k] + ((pop[i] * pop[j]) * (rate[i]-rate[j])^2 - m) / (pop[i] + pop[j])
          den_gh[k] = den_gh[k] + pop[i] * pop[j] / (pop[i] + pop[j])
        }
      }
    }
    gammah[k] = num_gh[k]/(2*den_gh[k])
    mh[k] = dh[k] / nh[k]
  }
  v = data.frame(dist = mh,  semivariance = gammah, npair = nh)
  list(weightsvar = wsvar, semivar = v)
}

#' Fit a theoretical semi-variogram to data
#'
#' Function fits (manually) a theoretical semi-variogram.
#' For now, model types available are spherical ('sph') and exponential ('exp').
#' Users should provide the experimental semi-variogram, fit a model type and set the semi-variogram parameters (nugget, range and sill).
#' With base function plot(), results can be visualized.

#' @param varexp, string, name of object, output of function `varexp()`
#' @param mod, character, the semi-variogram model type ('sph' or 'exp')
#' @param nug, numeric, nugget-effect value of the semi-variogram
#' @param ran, numeric, range value of the semi-variogram
#' @param sill, numeric, sill (or partial sill) value of the semi-variogram
#'
#' @return Function returns the following list of objects:
#' \item{structures}{the number of spatial strutures (excluding nugget-effect)}
#' \item{parameters}{a dataframe with model information to be passed to ssdpars()}
#' \item{fittedval}{a numeric vector of fitted values as set by the semi-variogram model}
#'
#' @export
varmodel = function (varexp, mod = c("exp","sph"), nug, ran , sill) {

  x = seq(1, ran , 1)
  xmax = seq(1, max(varexp[["semivar"]][,"dist"]), 1)
  if (mod=="sph") {
    vtype = 1
    model = nug + sill * (1.5 * (x / ran) - 0.5 * (x / ran)^3)
    cutoff = max(xmax) - ran
    msill = rep(nug + sill, cutoff )
    model = c(0, model, msill)
  }
  if (mod=="exp") {
    vtype = 2
    model = nug + sill * (1 - exp(- 3 * xmax / ran))
    model = c(0, model)
  }

  pars = data.frame (model = mod, modeltype = vtype, nugget = nug, range = ran, psill = sill)
  nstruc = 1 # for now only 1 (sph or exp)
  return(list(structures = nstruc, parameters = pars, fittedval = model))
}

#' Creates a parameters file and generates the simulated maps
#'
#' Creates a parameters file (.par) and invokes dss.c.64.exe to run block simulation program and
#' generate simulated map files (.out). As input you should provide name of objects
#' and parameter values required for the simulation process.

#' @param grdobj, string, name of list, output of function `grdfile()`
#' @param maskobj, string, name of list, output of function `maskfile()`
#' @param dfobj, string, name of list, output of function `irates()`
#' @param varmobj, string, name of list, output of function `varmodel()`
#' @param simulations, numeric, number of simulations
#' @param nrbias, numeric, nr simulations for bias correction
#' @param biascor, num vector, flag for (mean, variance) correction (yes = 1, no = 0)
#' @param ndMin, numeric, min number of neighbour observations used in  kriging
#' @param ndMax, numeric, max number of neighbour observations used in  kriging
#' @param nodMax, numeric, max number of previously simulated nodes used in  kriging
#' @param radius1, numeric, search radii in the major horizontal axe
#' @param radius2, numeric, search radii in the axe orthogonal (horizontal) to radius1
#' @param radius3, numeric, search radii in the vertical axe
#' @param ktype, numeric, the kriging type to be used (available are: 0 = simple, 1 = ordinary)
#'
#' @return Function returns the following list of objects:
#' \item{structures}{the number of spatial strutures (excluding nugget-effect)}
#' \item{parameters}{a dataframe with model information to be passed to ssdpars()}
#' \item{fittedval}{a numeric vector of fitted values as set by the semi-variogram model}
#'
#' @details Both parameters file (.par) and simulations files (.out) are stored in input folder. Note that the simulation process may take a while, depending mostly on the number of simulation nodes and number of simulations specified.
#'
#' @export
ssdpars = function (grdobj, maskobj, dfobj, varmobj, simulations = 1, nrbias = 20, biascor = c(1,1),
                    ndMin = 1, ndMax = 32, nodMax = 12, radius1, radius2, radius3 = 1, ktype = 1) {

  day = as.character(grdobj[["file"]]["day"])
  folder = as.character(grdobj[["file"]]["folder"])

  # filenames
  rf = as.character(dfobj[["file"]]["name"])
  bf = as.character(grdobj[["file"]]["name"])
  mf = as.character(maskobj[["file"]]["name"])

  # hard data pars
  nvars = as.numeric(dfobj[["ssdrpars"]]["nvars"])
  xcolumn = as.numeric(dfobj[["ssdrpars"]]["xcolumn"])
  ycolumn = as.numeric(dfobj[["ssdrpars"]]["ycolumn"])
  zcolumn = as.numeric(dfobj[["ssdrpars"]]["zcolumn"])
  varcol = as.numeric(dfobj[["ssdrpars"]]["varcol"])
  minval = as.numeric(dfobj[["ssdrpars"]]["minval"])
  maxval = as.numeric(dfobj[["ssdrpars"]]["maxval"])

  # mask grid pars
  nzones = as.numeric(maskobj[["zones"]]["nzones"])

  # output grid pars
  nx = as.numeric(unlist(grdobj[["gridpars"]]["nodes"]))[1]
  ny = as.numeric(unlist(grdobj[["gridpars"]]["nodes"]))[2]
  ox = as.numeric(unlist(grdobj[["gridpars"]]["origin"]))[1]
  oy = as.numeric(unlist(grdobj[["gridpars"]]["origin"]))[2]
  rx = as.numeric(unlist(grdobj[["gridpars"]]["resolution"]))[1]
  ry = as.numeric(unlist(grdobj[["gridpars"]]["resolution"]))[2]

  # general pars
  nas = as.numeric(grdobj[["gridpars"]]["NAs"])

  # varigram pars
  nstruct = as.numeric(varmobj[["structures"]])
  nugget = as.numeric(varmobj[["parameters"]]["nugget"])
  range = as.numeric(varmobj[["parameters"]]["range"])
  psill = as.numeric(varmobj[["parameters"]]["psill"])
  nuggetp = nugget/(nugget + psill)
  psillp = psill/(nugget + psill)
  mtype = as.numeric(varmobj[["parameters"]]["modeltype"])

  # block kriging pars
  maxblocks = as.numeric(grdobj[["outgrid"]]["nblock"])


  # ------ write ssdr.par for dss ------

  # store number of simulations
  nsim = simulations

  # store nr simulations for bias correction
  nbias = nrbias

  # store flag for (mean, var) correction (yes = 1, no = 0)
  biascor = biascor

  # draw pseudo-random value
  pseudon = sample(10^8,1)

  ## Run external DSS exectuable and return realization
  #
  # [input]
  #   simType (string) - SGEMS or GEOEAS type
  #   avgCorr (0/1) - correct mean
  #   varCorr (0/1) - correct variance
  #   usebihist (0/1) - use joint probability distributions
  #   inputPath (string) - path for folder with input data and DSS executable
  #   varIn (string) - name of the file with input experimental data
  #   noSim - number of realizations,
  #   outputFilePath (string) - full path for output file name without extension
  #   krigType (0/1/2/3/4/5) - kriging type
  #   secVar (string) - fullpath for secondary variable
  #   bihistFile (string) - fullpath for joint distribution file
  #   bounds (noZones x 2) - min and max of variable to be simulated
  #   XX (1 x 3) - Number of cells, origin, size in XX
  #   YY (1 x 3) - Number of cells, origin, size in YY
  #	  ZZ (1 x 3) - Number of cells, origin, size in ZZ
  #   localCorr (string) - fullpath for collocated correlation coefficient file
  #   globalCorr - global correlation coefficient between 2 variables
  #   auxVar string) - fullpath for sec variable for joint distribution
  #   krigANG (1 x 3) - azimuth, rake and dip
  #   krigRANGE (1 x 3) - major, minor, vertical
  #   varANG (1 x 3) - azimuth, rake and dip
  #   varRANGE (noZones x 3) - major, minor, vertical
  #   varType(noZones x 1) - variogram type
  #   varNugget (noZones x 1) - nugget effect
  #   zoneFileName(string) - fullpath for zone file
  #   noZones - number of zones in zonefile
  #   noClasses - number of classes for joint distribution
  #   rescale (0/1) - rescale secondary variable for co-simulation
  #   lvmFile (string) - fullpath for local varying means
  #   usePseudoHard (0/1) - for point distributions
  #   pseudoHardFile (string) - fullpath for local pdfs
  #   pseudoCorr 0/1) - correct local pdfs
  #   covTab (1 x 3) - size of the covariance matrix to be stored in mem
  #
  # [output]
  #   var_out - realization
  #
  #
  # Leonardo Azevedo - 2019
  # CERENA/Instituto Superior Tecnico (Portugal)
  #

  # create file path
  ssd_name = "ssdr"
  ssd_nameO = paste0(ssd_name, ".par")
  fssd = paste0(folder, "/", day, "_", ssd_nameO)

  if (file.exists(fssd)){
    file.remove(fssd)
  }

  # create notification file
  file.create(fssd)

  cat("#*************************************************************************************#", file = fssd, sep = "\n", append = TRUE)
  cat("#                                                                                     #", file = fssd, sep = "\n", append = TRUE)
  cat("#             PARALLEL DIRECT SEQUENCIAL SIMULATION PARAMETER FILE                    #", file = fssd, sep = "\n", append = TRUE)
  cat("#                                                                                     #", file = fssd, sep = "\n", append = TRUE)                                                                                     #');
  cat("#*************************************************************************************#", file = fssd, sep = "\n", append = TRUE)
  cat("#                                                                                     #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#         Remember, # - Comment; [GROUP] - parameter group; CAPS - parameter          #", file = fssd, sep = "\n", append = TRUE)
  cat("#         Also, no space allowed in paths/filenames                                   #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#                     here we define the hardata parameters                           #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("\n", file = fssd, append = TRUE)
  cat("[ZONES]", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("ZONESFILE = ", mf), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NZONES = ", nzones), file = fssd, sep = "\n", append = TRUE)
  cat("", file = fssd, sep = "\n", append = TRUE)

  for (i in 1:nzones) {
    cat(paste0("[HARDDATA", i, "]" ), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("DATAFILE = ", rf), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("COLUMNS = ", nvars), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("XCOLUMN = ", xcolumn), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("YCOLUMN = ", ycolumn), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("ZCOLUMN = ", zcolumn), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("VARCOLUMN = ", varcol), file = fssd, sep = "\n", append = TRUE)
    cat("WTCOLUMN = 0", file = fssd, sep = "\n", append = TRUE)
    cat(paste0("MINVAL = ", minval), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("MAXVAL = ", maxval), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("USETRANS = 1"), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("TRANSFILE =  Cluster.trn"), file = fssd, sep = "\n", append = TRUE)
  }

  cat("\n", file = fssd, append = TRUE )

  cat("[HARDDATA]", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("ZMIN = ", minval), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("ZMAX = ", maxval), file = fssd, sep = "\n", append = TRUE)
  cat("LTAIL = 1", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("LTPAR = ", minval), file = fssd, sep = "\n", append = TRUE)
  cat("UTAIL = 1", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("UTPAR = ", maxval), file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#                here we define parameters for the simulation                         #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[SIMULATION]", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("OUTFILE = ", day, "_","sim"), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NSIMS = ", nsim), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NTRY = ", nbias), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("AVGCORR = ", biascor[1]), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("VARCORR = ", biascor[2]), file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define the output grid (and secondary info grid)                     #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[GRID]", file = fssd, sep = "\n", append = TRUE)
  cat("# NX, NY and NZ are the number of blocks per direction", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NX = ", nx), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NY = ", ny), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NZ = ", 1), file = fssd, sep = "\n", append = TRUE)

  cat("# ORIGX, ORIGY and ORIGZ are the start coordinate for each direction", file = fssd, sep = "\n", append = TRUE)

  cat(paste0("ORIGX = ", ox), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("ORIGY = ", oy), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("ORIGZ = ", 0), file = fssd, sep = "\n", append = TRUE)

  cat("# SIZEX, SIZEY and SIZEZ is the size of blocks in each direction ", file = fssd, sep = "\n", append = TRUE)

  cat(paste0("SIZEX = ", rx ), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("SIZEY = ", ry ), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("SIZEZ = ", 1), file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#                 here we define some general parameters                              #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[GENERAL]", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("NULLVAL = ", nas), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("SEED = ", pseudon), file = fssd, sep = "\n", append = TRUE)
  cat("USEHEADERS = 1", file = fssd, sep = "\n", append = TRUE)
  cat("FILETYPE = GEOEAS", file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#                 here we define the parameters for search                            #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[SEARCH]", file = fssd, sep = "\n", append = TRUE)

  # ------- ssdr.par: search parameters -------

  # hard-coded values
  sstrat = 1
  mults = 0
  nmults = 1
  noct = 0
  sang1 = 0
  sang2 = 0
  sang3 = 0

  # min nr of observed samples
  cat(paste0("NDMIN  = ", ndMin), file = fssd, sep = "\n", append = TRUE)
  # max nr of observed samples
  cat(paste0("NDMAX  = ", ndMax), file = fssd, sep = "\n", append = TRUE)
  # max nr of previouly simulated nodes
  cat(paste0("NODMAX = ", nodMax), file = fssd, sep = "\n", append = TRUE)
  # Two-part search / data nodes flag
  cat(paste0("SSTRAT = ", sstrat), file = fssd, sep = "\n", append = TRUE)
  # Multiple grid simulation flag
  cat(paste0("MULTS  = ", mults), file = fssd, sep = "\n", append = TRUE)
  # Nr of multiple grid refinements
  cat(paste0("NMULTS = ", nmults), file = fssd, sep = "\n", append = TRUE)
  # Nr of original data per octant
  cat(paste0("NOCT = ", noct), file = fssd, sep = "\n", append = TRUE)
  # Search radii in the major horizontal axe
  cat(paste0("RADIUS1 = ", radius1), file = fssd, sep = "\n", append = TRUE)
  # Search radii in the ortogonal horizontal axe (to major)
  cat(paste0("RADIUS2 = ", radius2), file = fssd, sep = "\n", append = TRUE)
  # Search radii in the vertical axe
  cat(paste0("RADIUS3 = ", radius3), file = fssd, sep = "\n", append = TRUE)
  # Orientation angle parameter of direction I (degrees)
  cat(paste0("SANG1 = ", sang1), file = fssd, sep = "\n", append = TRUE)
  # Orientation angle parameter of direction II (degrees)
  cat(paste0("SANG2 = ", sang2), file = fssd, sep = "\n", append = TRUE)
  # Orientation angle parameter of direction III (degrees)
  cat(paste0("SANG3 = ", sang3), file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#   here we define the kriging information, and secondary info when applicable        #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  # ------- ssdr.par: krige info -------

  # hard-coded values
  colorcorr = 0
  softfile = "no file"
  lvmfile = "no file"
  nvaril = 1
  icollvm = 1
  ccfile = "no file"
  rescale = 0

  cat("[KRIGING]", file = fssd, sep = "\n", append = TRUE)
  # Kriging type: 0 = simple, 1 = ordinary, 2 = simple with locally varying mean
  # 3 = external drift, 4 = collo-cokrig global CC, 5 = local CC
  cat(paste0("KTYPE = ", ktype), file = fssd, sep = "\n", append = TRUE)
  # Global coef correlation (ktype = 4)
  cat(paste0("COLOCORR = ", colorcorr), file = fssd, sep = "\n", append = TRUE)
  # Filename of the soft data (ktype = 2)
  cat(paste0("SOFTFILE = ", softfile), file = fssd, sep = "\n", append = TRUE)
  # For ktype = 2
  cat(paste0("LVMFILE  = ", lvmfile), file = fssd, sep = "\n", append = TRUE)
  # Number of columns in the secundary data file
  cat(paste0("NVARIL  = ", nvaril), file = fssd, sep = "\n", append = TRUE)
  # Column number of secundary variable
  cat(paste0("ICOLLVM  = ", icollvm), file = fssd, sep = "\n", append = TRUE)
  # Filename of correlation file for local correlations (ktype = 5)
  cat(paste0("CCFILE  = ", ccfile), file = fssd, sep = "\n", append = TRUE)
  # Rescale secondary variable (for ktype >= 4)
  cat(paste0("RESCALE  = ", rescale), file = fssd, sep = "\n", append = TRUE)

  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define the variogram to use. if more than 1, use [VARIOGRAM2]        #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  # ------- ssdr.par: semi-variogram models -------

  for (j in 1 : nzones) {
    cat(paste0("[VARIOGRAMZ", j, "]"), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("NSTRUCT = ", nstruct), file = fssd, sep = "\n", append = TRUE)
    cat(paste0("NUGGET = ", nuggetp), file = fssd, sep = "\n", append = TRUE)
    for (i in 1 : nstruct){
      cat(paste0("[VARIOGRAMZ", j, "S", i, "]"), file = fssd, sep = "\n", append = TRUE)
      # store struture type ; 1 = spherical, 2 = exponential
      cat(paste0("TYPE = ", mtype), file = fssd, sep = "\n", append = TRUE)
      # C parameter "COV + NUGGET = 1.0" (CC(i))
      cat("COV = 1", file = fssd, sep = "\n", append = TRUE)
      # Geometric anisotropy angle I (ANG1(i))
      cat("ANG1 = 0", file = fssd, sep = "\n", append = TRUE)
      # Geometric anisotropy angle II (ANG2(i))
      cat("ANG2 = 0", file = fssd, sep = "\n", append = TRUE)
      # Geometric anisotropy angle III (ANG3(i))
      cat("ANG3 = 0", file = fssd, sep = "\n", append = TRUE)
      # Maximum horizontal range (AA(i))
      cat(paste0("AA = ", range), file = fssd, sep = "\n", append = TRUE)
      # Minimum horizontal range (AA1)
      cat(paste0("AA1 = ", range), file = fssd, sep = "\n", append = TRUE)
      # Vertical range (AA2)
      cat("AA2 = 1", file = fssd, sep = "\n", append = TRUE)
    }
  }

  cat("#-------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define parameters for joint DSS        #", file = fssd, sep = "\n", append = TRUE)
  cat("#-------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  # hard-coded values
  usebihist = 0
  bihistfile = "no file"
  nclasses = 0
  auxfile = "no file"

  for (j in 1:nzones){
    cat(paste0("[BIHIST", j, "]"), file = fssd, sep = "\n", append = TRUE)
    #Use Bihist? 1-yes 0-no'
    cat(paste0("USEBIHIST = ", usebihist), file = fssd, sep = "\n", append = TRUE)
    # bihistogram file
    cat(paste0("BIHISTFILE = ", bihistfile), file = fssd, sep = "\n", append = TRUE)
    # number of classes to use
    cat(paste0("NCLASSES  = ", nclasses), file = fssd, sep = "\n", append = TRUE)
    # auxiliary image
    cat(paste0("AUXILIARYFILE  = ", auxfile), file = fssd, sep = "\n", append = TRUE)
  }

  cat("#---------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define the debug parameters - probably you wont need this        #", file = fssd, sep = "\n", append = TRUE)
  cat("#---------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[DEBUG]", file = fssd, sep = "\n", append = TRUE)
  # 1 to 3, use higher than 1 only if REALLY needed
  cat("DBGLEVEL = 2", file = fssd, sep = "\n", append = TRUE)
  # File to write debug
  cat("DBGFILE   = debug.dbg ", file = fssd, sep = "\n", append = TRUE)

  cat("#---------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define parameters for COVARIANCE TABLE - reduce if memory is a problem       #", file = fssd, sep = "\n", append = TRUE)
  cat("#---------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[COVTAB]", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("MAXCTX = ", range), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("MAXCTY = ", range), file = fssd, sep = "\n", append = TRUE)
  cat("MAXCTZ = 1", file = fssd, sep = "\n", append = TRUE)

  cat("#--------------------------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)
  cat("#        here we define parameters for BLOCK KRIGING - if you are not block kriging useblocks should be 0      #", file = fssd, sep = "\n", append = TRUE)
  cat("#--------------------------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = TRUE)

  cat("[BLOCKS]", file = fssd, sep = "\n", append = TRUE)
  cat("USEBLOCKS = 1", file = fssd, sep = "\n", append = TRUE)
  cat(paste0("BLOCKSFILE = ", bf), file = fssd, sep = "\n", append = TRUE)
  cat(paste0("MAXBLOCKS  = ", maxblocks), file = fssd, sep = "\n", append = TRUE)
  cat("[PSEUDOHARD]", file = fssd, sep = "\n", append = TRUE)
  # 1 use, 0 no  pseudo hard data is point distributions that are simulated before all other nodes
  cat("USEPSEUDO = 0", file = fssd, sep = "\n", append = TRUE)
  # file
  cat(paste0("PSEUDOFILE = ", "no file"), file = fssd, sep = "\n", append = TRUE)
  # correct simulated value with point
  cat("PSEUDOCORR = 0", file = fssd, sep = "\n", append = TRUE)

  # ------ run dss ------

  wd = getwd()
  parfile = paste0(day, "_", ssd_nameO)
  setwd("./input")
  if(file.exists(parfile)) {
    cat(paste0(parfile, " created.\n"))
  }

  message("Running dss.c.64.exe. This process may take a while.")
  packageStartupMessage("initializing ...", appendLF = FALSE)
  system(paste0("./DSS.C.64.exe ", parfile))
  packageStartupMessage(" done")

  setwd(wd)

}


#' Create raster objects from simulations
#'
#' Function read simulation files (.out) returned by `ssdpars()`
#' and returns a list with simulated maps (rasterstack object),
#' e-type and uncertainty maps (rasterlayers).

#' @param grdobj, string,  name of list, output of function `grdfile()`
#' @param grids, if grids = T  saves simulated maps in 'native' raster package format .grd
#' @param emaps, if emaps = T (default), saves e-type and uncertainty maps in format .grd
#'
#' @return
#' \item{simulations}{a rasterstack (package 'raster') where each layer is a simulation}
#' \item{etype}{a rasterlayer representing median-etype map for disease risk}
#' \item{uncertainty}{a rasterlayer representing risk uncertainty map for disease risk}
#'
#' @details All .grd files are geographic (spatial) data in 'raster' format, and are stored in input folder.
#'
#' @importFrom sp CRS
#' @importFrom raster raster
#'
#'
#' @export

 outraster = function (grdobj, grids = FALSE, emaps = TRUE) {

   day = grdobj[["file"]]["day"]
   folder = grdobj[["file"]]["folder"]

   # create prefix .out filenames
   simout = paste0(day, "_", "sim")

   # store list .out filenames
   lf = list.files(paste0(folder,"/"), pattern ="\\.out$")
   dss_list = list(simnames = Filter(function(x) grepl(simout, x), lf))

   # store .out NA value
   bNA <- as.integer(grdobj[["gridpars"]]["NAs"])


   # store number of simulations
   nsims = length(dss_list[["simnames"]])
   ssims = raster::stack()
   namesims = c(rep(0, nsims))

   # loop each simulation
   for (k in 1:nsims){
     print(dss_list[["simnames"]][k])
     out = read.table(file = paste0(folder, "/", dss_list[["simnames"]][k]), sep = " ", skip=3)
     out01 = as.data.frame(out)
     out01[out01 == bNA] <- NA

     # set grid size
     # mc : nr columns, mr : nr rows
     mc = grdobj[["gridpars"]][["nodes"]][1]
     mr = grdobj[["gridpars"]][["nodes"]][2]

     xmin = grdobj[["ingrid"]]@extent[1]
     xmax = grdobj[["ingrid"]]@extent[2]
     ymin = grdobj[["ingrid"]]@extent[3]
     ymax = grdobj[["ingrid"]]@extent[4]

     crsname = as.character(grdobj[["ingrid"]]@crs)

     # create matrix
     out02 = matrix(0, nrow = mr, ncol = mc)

     count = 1
     for (i in 1:mr){
       for (j in 1:mc){
         out02[i, j] = out01[count, 1]
         count = count + 1
         }
     }

     # create matrix
     out03 = matrix(NA, nrow = mr, ncol = mc)

     # populate matrix
     for (i in 1:mr){
       out03[i, ]<- c(out02[mr - i + 1,] )
     }
     r = raster::raster(out03)
     raster::extent(r)=c(xmin, xmax, ymin, ymax)
     raster::projection(r) = sp::CRS(crsname)

     if(grids == TRUE){
       gridname <- paste0(folder, "/", day, "_", "sim", k)
       raster::writeRaster(r, filename = gridname, overwrite = TRUE)
       if(file.exists(paste0(gridname, ".gri"))){
         cat(gridname, ".gri (and .grd) created.\n")
       }
     }

     # add layer to stack
     ssims = raster::stack(ssims, r)
     # change layer name
     names(ssims[[k]]) = paste0("sim", k)
   }
   etype = raster::calc(ssims, fun = function(x) {quantile(x, probs = .5,na.rm=TRUE)})
   uncer = raster::calc(ssims, fun = sd, na.rm = TRUE)
   if (emaps == TRUE){
     etypeName <- paste0(folder, "/", day, "_", "medn")
     raster::writeRaster(etype, filename = etypeName, overwrite = TRUE)
     if(file.exists(paste0(etypeName, ".gri"))){
       cat(etypeName, ".gri (and .grd) created.\n")
     }
     uncrName <- paste0(folder, "/", day, "_", "uncr")
     raster::writeRaster(uncer, filename = uncrName, overwrite = TRUE)
     if(file.exists(paste0(uncrName, ".gri"))){
       cat(uncrName, ".gri (and .grd) created.")
     }
   }
   listmaps = list(simulations = ssims, etype = etype, uncertainty = uncer )
   return(listmaps)
   }

 #' Creates a pixelated map
 #'
 #' This is a wrapper function calling R package `pixelate`^[https://github.com/aimeertaylor/pixelate]
 #' developed by Aimee Taylor and colleagues ^[https://doi.org/10.48550/arXiv.2005.11993] providing tools
 #' to plot disease risk mapping with a visual representation of spatial uncertainty,
 #' as a function of pixel size.
 #'
 #' As input you should provide the result returned by `outraster()` (a list of `RasterLayer` objects)
 #' and additional arguments required for pixelation.
 #
 #' @param mapobj string, name of list, output of function `outraster()`
 #' @param nbigk vector with length 2, specifies the minimum number of large pixels in the x and y directions, respectively.
 #' @param bigk integer, specifies the number of average quantile intervals (i.e. number of different pixel sizes).
 #' @param scaleft integer, specifies a factor (in units of observations) that features in either iterative multiplication or iterative exponentiation (see `scale_factor` arg in `pixelate::pixelate() for more details).
 #' @param legname character, legend title.
 #'
 #' @importFrom ggplot2 ggplot
 #'
 #' @return `pxmap()` is used to map disease risk with varying pixel sizes
 #' showing spatial uncertainty generated by
 #' block direct sequential simulation algorithm.
 #' @export
 pxmap = function(mapobj = maps, nbigk = c(17, 17), bigk = 4, scaleft = 1, legname = "Median\nCum. risk"){

   if (!require(devtools)) {
     utils::install.packages("devtools", type = "source")
   }

   if(!require(pixelate, quietly = TRUE)){
     devtools::install_github(repo = "aimeertaylor/pixelate", build_vignettes = TRUE, force = TRUE)
   }

   if (!require(gpclib)) install.packages("gpclib", type="source")
   maptools::gpclibPermit()

   # prepare EpiGeostats() maps for pixelate()
   px_u <- raster::as.data.frame(mapobj[["uncertainty"]], row.names = NULL, optional = FALSE, xy = TRUE, na.rm = FALSE)
   px_m <- raster::as.data.frame(mapobj[["etype"]] , row.names = NULL, optional = FALSE, xy = TRUE, na.rm = FALSE)
   px_in <- cbind(px_m, px_u[,3])
   names(px_in) <- c("x", "y", "z", "u")

   # Pixelate using default parameters
   px_def <- pixelate::pixelate(px_in, num_bigk_pix = nbigk, bigk = bigk, scale_factor = scaleft )

   # Inspect list returned by pixelate
   # str(px_def)

   # Inspect a sample of uncertain pixelated predictions
   uncertain_ind <- which(px_def$pix_df$u > 0)
   head(px_def$pix_df[uncertain_ind, ])

   # Plot pixelated map with ggplot
   ggplot2::ggplot(px_def$pix_df) +
     # Add raster surface
     ggplot2::geom_raster(mapping = ggplot2::aes(x = x, y = y, fill = pix_z)) +

     # Add gradient
     ggplot2::scale_fill_viridis_c(name = legname, na.value = 'white') +

     # Add axis labels
     ggplot2::ylab('Y') +
     ggplot2::xlab('X') +

     # Ensure the plotting space is not expanded
     ggplot2::coord_fixed(expand = FALSE) +

     # Modify the legend and add a plot border:
     ggplot2::theme(legend.justification = c(0, 0),
                    legend.position = c(1, 0),
                    legend.background = ggplot2::element_rect(fill = NA),
                    legend.title = ggplot2::element_text(size = 10),
                    legend.text = ggplot2::element_text(size = 10))
 }

 #' Creates elegant simulations, median e-type/uncertainty maps
 #'
 #' The function wraps functions from ggplot2 to provide a more elegant map.
 #'
 #' As input you should provide a `RasterLayer` object returned by `outraster()`.
 #
 #' @param mapobj string, name of list, output of function `outraster()`
 #' @param mapvar character, the map to be plotted, one of the following c("etype", "unc", "simulations")
 #' @param simid numeric, for mapvar = "simulations", is the simulation number to be plotted.
 #' @param legname character, legend title.
 #'
 #' @return `spmap()` is used to map spatially continuous maps representing disease risk simulations, median e-type disease risk
 #' or spatial uncertainty (inter-quartile difference) generated by the block direct sequential simulation algorithm.

 #' @export
 spmap = function(mapobj = maps, mapvar = "etype", simid = 1, legname = ""){

   m <- mapvar

   if(m == "simulations"){
     id <- simid
     if(dim(mapobj[["simulations"]])[3] < id){stop("simid > number of simulations")}
     # prepare map for ggplot()
     px_in <- raster::as.data.frame(mapobj[[m]][[id]],
                                    row.names = NULL,
                                    optional = FALSE,
                                    xy = TRUE,
                                    na.rm = FALSE)
   } else {
     # prepare map for ggplot()
     px_in <- raster::as.data.frame(mapobj[[m]],
                                    row.names = NULL,
                                    optional = FALSE,
                                    xy = TRUE,
                                    na.rm = FALSE)
   }
   # assign new names to data.frame columns
   names(px_in) <- c("x", "y", "z")

   # Plot map with ggplot
   p <- ggplot2::ggplot(px_in) +
     # Add raster surface
     ggplot2::geom_raster(mapping = ggplot2::aes(x = x, y = y, fill = z))

   # Add gradient
   if (m == "uncertainty"){
     px <- p +
       ggplot2::scale_fill_viridis_c(option = 'inferno', name = legname, na.value = 'white')
   } else {
     px <- p +
       ggplot2::scale_fill_viridis_c(option = 'viridis', name = legname, na.value = 'white')
   }

   # Add axis labels
   px +
     ggplot2::ylab('Y') +
     ggplot2::xlab('X') +
     # Ensure the plotting space is not expanded
     ggplot2::coord_fixed(expand = FALSE) +
     # Modify the legend and add a plot border:
     ggplot2::theme(legend.justification = c(0, 0),
                    legend.position = c(1,0),
                    legend.background = ggplot2::element_rect(fill = NA),
                    legend.title = ggplot2::element_text(size = 10),
                    legend.text = ggplot2::element_text(size = 10))
 }
