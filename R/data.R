#' @title COVID-19 notification rate in Portugal mainland, 2022-01-15
#' @description Refers to the 14-day notification rate per
#' 100000 population of new COVID-19 cases in Portugal on 2022-01-15.
#'
#' @format A \code{"data.frame"} with 278 rows and 6 columns,
#' refers to the 14-day notification rate per
#' 100000 population of new COVID-19 cases in Portugal on 2022-01-15.
#' \describe{
#'   \item{x}{x-coordinates, in meters}
#'   \item{y}{y-coordinates, in meters}
#'   \item{t}{z-coordinates, in meters}
#'   \item{ncases}{number of cases}
#'   \item{oid_}{name for region id (municipality)}
#'   \item{pop19}{population size}
#' }
#'
#'
#' @docType data
#'
#' @usage data(ptdata)
#'
#'
#' @keywords datasets
#'
#' @references Azevedo, L. et al. Geostatistical COVID-19 infection risk maps for Portugal. Int J Health Geogr 19, 25 (2020).
#' (\href{https://doi.org/10.1186/s12942-020-00221-5})
#'
#' @source \href{https://cerena.pt/projects/scope-spatial-data-sciences-covid-19-pandemic}{CERENA-IST/UL}
#'
#' @examples
#' data(ptdata)
#' str(ptdata)
"ptdata"
