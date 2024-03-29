% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.aftsur.R
\name{print.aftsur}
\alias{print.aftsur}
\title{Print Fit Output for Semi-parametric Accelerated Failure Time Mixture Cure Model}
\usage{
\method{print}{aftsur}(fit)
}
\arguments{
\item{fit}{an object of class "aftsur", a model fitted by the \code{aftsur} function.}
}
\description{
a function that prints a model result.
}
\examples{
require(survival)
# load data
data("simu_data")
# create a formula using a Surv
formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
# fit a model
aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
print(aft_fit)
}
