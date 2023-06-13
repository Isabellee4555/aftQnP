#' Simulated Partly Interval-censored Data
#' @description A simulated partly interval-censored dataset.
#' @format \code{simu_data} is a tibble with 500 rows of simulated observations and 11 columns.
#' \describe{
#' \item{\code{X1, X2}}{latency covariates}
#' \item{\code{Z1, Z2, Z3}}{incidence covariates, where \code{Z1} is an intercept term}
#' \item{\code{delta, deltaL, deltaL, deltaR}}{indicators for event, left, right, and interval-censoring}
#' \item{\code{y_L, y_R}}{event and censoring times. For events, \code{y_L = y_R}. If the observation is left-censored, \code{y_L = NA}. If the observation is right-censored, \code{y_R = NA}.}
#' }
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data("simu_data")
#'
#' @details
#' The latency covariates \code{X1} and \code{X2} are random samples generated from
#' \code{rbinom(n, 1, 0.5)} and \code{rnorm(n, 0, 0.5)} respectively.
#' The incidence covariates \code{Z2} and \code{Z3} are random samples generated from
#' \code{runif(n, -0.5, 0.5)}. The true parameters for the latency are set as
#' \code{beta = c(0.2, -0.5)}, and \code{gamma = c(1.5, -1.2, -1.5)} for the incidence.
#' The conditional event proportion (conditional on the non-cured observations) is set to 0.8.
#'
#'
#'
#' @examples
#' data("simu_data")
#' mydata <- simu_data
#'
#'
"simu_data"
