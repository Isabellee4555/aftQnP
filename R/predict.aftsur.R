#' Predicted Survival or Hazard Rates
#' @description a function that predicts survival or hazard from a semi-parametric AFT mixture cure model.
#' @param fit an object of class "aftsur", a model fitted by the \code{aftsur} function.
#' @param type a character string to indicate type of prediction. Either "survival" or "hazard" should be supplied.
#' @param time a vector of time points at which the survival and hazard functions should be evaluated. A limit is set for \code{time} to ensure that the estimates are obtained for the susceptible group; see details below.
#' @param interval a switch indicating if 95% confidence intervals are required.
#' If \code{time = NULL} (default), a sequence of estimates will be provided from timepoint 0 to the upper limit.
#' @param x a dataframe which contains latency variables with the same data structure as the training data. If \code{x = NULL}, the baseline survival or hazard will be predicted.
#' @details
#' The upper limit for \code{time} is set as the 99% empirical quantile of follow-up times among non-right-censored observations.
#' @examples
#' require(survival)
#' # load data
#' data("simu_data")
#' # create a formula using a Surv
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
#' # fit a model
#' aft_fit <- aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
#' # make prediction for the baseline survival at time 0.5, 1 and 1.5
#' predict(aft_fit, x = NULL, type = "survival", time = c(0.5, 1, 1.5), interval = TRUE)
#' # make prediction using the first 5 observations from the training data at time 0.5
#' predict(aft_fit, x = simu_data[1:5, c("X1", "X2")], type = "survival", time = c(0.5), interval = TRUE)
#' @returns
#' If \code{interval = FALSE}, a data frame that contain timepoints evaluated and estimates will be return.
#' If \code{interval = TRUE}, two additional columns, representing the lower and upper bounds of the 95% confidence intervals, will also be returned.
#' @export

predict.aftsur <-
  function(fit,
           x = NULL,
           type,
           time = NULL,
           interval = FALSE) {
    if (!any(type %in% c("survival", "hazard"))) {
      stop("type needs to be either 'survival' or 'hazard'")
    }

    X_base <-  if (is.null(x)) {
      matrix(0, ncol = length(fit$beta))
    } else{

      if (!is.data.frame(x)){
        stop("x needs to be supplied as a dataframe that contains latency variable with the same structure as the training data.")
      }
      model.matrix(fit$call$latency[-2], x)[,-1]
    }

    if(!is.matrix(X_base)){
      X_base <- matrix(X_base, nrow = 1)
    }

    estimates_list <- list()

    for(i in 1:nrow(X_base)){
      X_base_i <- X_base[i,]

      max_event_t <-
        quantile(fit$k[fit$deltaR == 0] * c(exp(X_base_i %*% fit$beta)), 0.99)

      const <- exp(-X_base_i %*% fit$beta)

      if (!is.null(time)) {
        if (max(time) > max_event_t) {
          stop(paste0(
            "time evaluated is outside the limit ",
            round(max_event_t, 2)
          ))
        }
        t_grid <-  time
      } else{
        t_grid <- seq(0, max_event_t, 0.01)
      }

      if (type == "survival") {
        estimates <-  exp(-H0(t_grid %*% const, fit$theta, fit$basis_params))

      } else{
        estimates <-
          h0(t_grid %*% const, fit$theta, fit$basis_params) * c(const)
      }

      if (interval == TRUE) {
        if (type == "survival") {
          variance_H <-
            diag(
              compute_cumu_hazard_covariance(
                fit,
                compute_covariance_matrix(fit, return_F = T),
                t_grid %*% const
              )
            )
          se <- sqrt(estimates ^ 2 * variance_H)
        } else{
          se <-
            sqrt(diag(
              compute_hazard_covariance(
                fit,
                compute_covariance_matrix(fit, return_F = T),
                t_grid %*% const
              )
            ))
        }

        estimates <- data.frame(
          times = t_grid,
          estimates = estimates,
          "lower" = estimates - 1.959964 * se,
          "upper" = estimates + 1.959964 * se
        )

      } else{
        estimates <- data.frame(times = t_grid,
                                estimates = estimates)
      }
      estimates[estimates < 0] <- 0

      estimates_list[[i]] <- estimates
    }


    return(do.call(rbind, estimates_list))

  }
