% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.aftsur.R
\name{predict.aftsur}
\alias{predict.aftsur}
\title{Predicted Survival or Hazard Rates}
\usage{
\method{predict}{aftsur}(fit, x = NULL, type, time = NULL, interval = FALSE)
}
\arguments{
\item{fit}{an object of class "aftsur", a model fitted by the \code{aftsur} function.}

\item{x}{a dataframe which contains latency variables with the same data structure as the training data. If \code{x = NULL}, the baseline survival or hazard will be predicted.}

\item{type}{a character string to indicate type of prediction. Either "survival" or "hazard" should be supplied.}

\item{time}{a vector of time points at which the survival and hazard functions should be evaluated. A limit is set for \code{time} to ensure that the estimates are obtained for the susceptible group; see details below.}

\item{interval}{a switch indicating if 95\% confidence intervals are required.
If \code{time = NULL} (default), a sequence of estimates will be provided from timepoint 0 to the upper limit.}
}
\value{
If \code{interval = FALSE}, a data frame that contain timepoints evaluated and estimates will be return.
If \code{interval = TRUE}, two additional columns, representing the lower and upper bounds of the 95\% confidence intervals, will also be returned.
}
\description{
a function that predicts survival or hazard from a semi-parametric AFT mixture cure model.
}
\details{
The upper limit for \code{time} is set as the 99\% empirical quantile of follow-up times among non-right-censored observations.
}
\examples{
require(survival)
# load data
data("simu_data")
# create a formula using a Surv
formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
# fit a model
aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
# make prediction for the baseline survival at time 0.5, 1 and 1.5
predict(aft_fit, x = NULL, type = "survival", time = c(0.5, 1, 1.5), interval = TRUE)
# make prediction using the first 5 observations from the training data at time 0.5
predict(aft_fit, x = simu_data[1:5, c("X1", "X2")], type = "survival", time = c(0.5), interval = TRUE)
}
