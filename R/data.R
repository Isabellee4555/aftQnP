#' Simulated Partly Interval-censored Data
#' @description A simulated partly interval-censored dataset.
#' @format This dataset includes 100 observations.
#' \describe{
#' \item{\code{X1 X2 X3}}{incidence covariates}
#' \item{\code{Z1 Z2 Z3}}{latency covariates, where \code{Z1} is an intercept term}
#' \item{\code{delta deltaL deltaL deltaR}}{indicators for event, left, right, and interval-censoring}
#' \item{\code{y_L y_R}}{event and censoring times. For events, \code{y_L = y_R}. \code{y_L = NA} if the observation is left-censored and \code{y_R = NA} if the observation is right-censored.}
#' }
#' @docType data
#'
#' @usage data(ptces)
#'
#'
#'
#'
"ptces"

