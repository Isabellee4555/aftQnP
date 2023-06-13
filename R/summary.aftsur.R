#' Summarising a Model Fit for Semi-parametric Accelerated Failure Time Mixture Cure Model
#' @description \code{summary} method for class \code{"aftsur"}.
#' @param fit an object of class "aftsur", a model fitted by the \code{aftsur} function.
#' @examples
#' require(survival)
#' # load data
#' data("simu_data")
#' # create a formula using a Surv
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
#' # fit a model
#' aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
#' summary(aft_fit)
#' @export

summary.aftsur <- function(fit) {
  aft_beta <- t(t(fit$beta))

  logit_gamma <- t(t(fit$gamma))

  sd_ <- sqrt(diag(compute_covariance_matrix(fit, return_F = T)))

  betasd <-  t(t(sd_[1:length(fit$beta)]))

  gammasd <-
    t(t(sd_[(length(fit$beta) + 1):(length(fit$beta) + length(fit$gamma))]))

  betaz <- fit$beta / betasd

  gammaz <- fit$gamma / gammasd

  beta.pfit <- pnorm(abs(betaz), lower.tail = F) * 2

  gamma.pfit <- pnorm(abs(gammaz), lower.tail = F) * 2

  beta_out <- cbind(aft_beta, betasd, betaz, beta.pfit)

  gamma_out <- cbind(logit_gamma, gammasd, gammaz, gamma.pfit)

  colnames(beta_out) <-
    c("Estimate", "Std.Error", "Z Value", "Pr(>|Z|)")

  colnames(gamma_out) <-
    c("Estimate", "Std.Error", "Z Value", "Pr(>|Z|)")

  fit$pzllh <- penalised_likelihood(fit)

  # cat(
  #   "---------------------------------------------------------------------------"
  # )
  # cat("\nSemi-parametric Accelerated Failure Time Mixture Cure Model Using MPL\n")
  # cat(paste("\nPenalised log-likelihood:", round(fit$pzllh, 3), "\n"))
  # cat(paste(
  #   "\nEstimated smoothing parameter:",
  #   round(fit$lambda, 3),
  #   "\n"
  # ))
  # cat("=====")
  # cat("\n")
  # cat("\nAccelerated Failure Time Model:\n")
  # print(beta_out)
  # cat("\n")
  # cat("=====")
  # cat("\nLogistic Model:\n")
  # print(gamma_out)
  # cat(
  #   "---------------------------------------------------------------------------"
  # )
  # cat("\n")

  out <-
    list(
      beta = beta_out,
      gamma = gamma_out,
      pzllh = fit$pzllh,
      lambda = fit$lambda
    )

  class(out) <- "summary.aftsur"
  out
}

#' Printing Summary Output
#' @description A function that prints the summary output for a model result.
#' @param x an object of class "summary.aftsur".
#' @export
print.summary.aftsur <- function(x){
  cat(
    "---------------------------------------------------------------------------"
  )
  cat("\nSemi-parametric Accelerated Failure Time Mixture Cure Model Using MPL\n")
  cat(paste("\nPenalised log-likelihood:", round(x$pzllh, 3), "\n"))
  cat(paste(
    "\nEstimated smoothing parameter:",
    round(x$lambda, 3),
    "\n"
  ))
  cat("=====")
  cat("\n")
  cat("\nAccelerated Failure Time Model:\n")
  print(x$beta)
  cat("\n")
  cat("=====")
  cat("\nLogistic Model:\n")
  print(x$gamma)
  cat(
    "---------------------------------------------------------------------------"
  )
  cat("\n")
}
