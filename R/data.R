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
#'   \item{oid_}{name for region id (municipality)}
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
