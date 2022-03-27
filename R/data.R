#' @title ptdata dataset
#' @description The ptdata dataset refers to a data.frame with the
#' the 14-day cumulative incidence of covid-19 in Portugal
#' on the 15th January 2021, by municipality. It contains geographic coordinates of
#' municipalities (population-weighted centroids) and population at risk.
#' This dataset can be used to follow the step-by-step guide
#' to obtain block direct simulation risk maps of COVID-19,
#' median risk map and spatial risk uncertainty map.
#' @format A data.frame with 278 rows and 6 columns
#' \describe{
#'   \item{x}{x-coordinates, in meters}
#'   \item{y}{y-coordinates, in meters}
#'   \item{t}{z-coordinates, in meters}
#'   \item{ncases}{number of new cases}
#'   \item{oid_}{number for region id (municipality id code)}
#'   \item{pop19}{population at risk}
#' }
#' @docType data
#' @usage data(ptdata)
#' @keywords datasets
#' @references Azevedo, L. et al. Geostatistical COVID-19 infection risk maps for Portugal. Int J Health Geogr 19, 25 (2020).
#' @source Portuguese Directorate-General for Health (DGS). Data preparation by Manuel Ribeiro & Leonardo Azevedo (leonardo.azevedo@tecnico.ulisboa.pt).
#' @examples
#' data(ptdata)
#' str(ptdata)
"ptdata"

#' @title ptgrid dataset
#' @description The ptgrid dataset refers to a data.frame representing
#' a rectangular grid with 2 km x 2 km spacing that covers Portugal mainland.
#' Grid nodes with real values refer to the municipality id (block data).
#' Grid nodes located on sea or Spain country are set to `NA`.
#'
#' @format A data.frame with 40608 rows and 3 columns
#' \describe{
#'   \item{x}{x-coordinates, in meters}
#'   \item{y}{y-coordinates, in meters}
#'   \item{oid_}{id, municipality's id number)}
#' }
#' @docType data
#' @usage data(ptgrid)
#' @keywords datasets
#' @references Azevedo, L. et al. Geostatistical COVID-19 infection risk maps for Portugal. Int J Health Geogr 19, 25 (2020).
#' @source Data preparation by Manuel Ribeiro & Leonardo Azevedo (leonardo.azevedo@tecnico.ulisboa.pt).
#' @examples
#' data(ptgrid)
#' str(ptgrid)
"ptgrid"
