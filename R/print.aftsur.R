#' Print Fit Output for Semi-parametric Accelerated Failure Time Mixture Cure Model
#' @description a function that prints a model result.
#' @param fit an object of class "aftsur", a model fitted by the \code{aftsur} function.
#' @examples
#' require(survival)
#' # load data
#' data("simu_data")
#' # create a formula using a Surv
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
#' # fit a model
#' aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
#' print(aft_fit)
#' @export

print.aftsur <- function(fit){

  aft_beta <- t(t(fit$beta))

  logit_gamma <- t(t(fit$gamma))

  colnames(aft_beta) <-  colnames(logit_gamma) <- "Estimate"

  rownames(aft_beta) <- colnames(fit$X)

  rownames(logit_gamma) <- colnames(fit$Z)

  cat("Call:\n")
  cat(paste0(as.character(fit$call$latency)[c(2,1,3)], collapse = " "))
  cat("\n")
  cat(paste("E", paste0(fit$call$incidence, collapse = " ")))
  cat("\n")
  cat("---------------------------------------------------------------------------")
  cat("\nSemi-parametric Accelerated Failure Time Mixture Cure Model Using MPL\n")
  cat("=====")
  cat("\n")
  cat("\nAccelerated Failure Time Model:\n")
  print(aft_beta)
  cat("\n")
  cat("=====")
  cat("\nLogistic Model:\n")
  print(logit_gamma)
  cat("---------------------------------------------------------------------------")
  cat("\n")
}
