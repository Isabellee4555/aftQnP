#' Semi-parametric Accelerated Failure Time Mixture Cure Model
#' @description A method to fit semi-parametric AFT mixture cure model using the maximum penalised likelihood estimation, where left, right, interval-censored data are allowed.
#'
#' @param data a data frame which includes survival times and variables for the AFT model and the logistic regression.
#' @param formula_aft a formula for AFT model ~ X1 + X2.
#' @param formula_logit a formula for AFT model ~ Z1 + Z2.
#' @param y_L a vector of left censored times.
#' @param y_R a vector of right censored times.
#' @param lambda initial value for smoothing parameter. By default, lambda = 1e-5.
#' @param knots the number of Gaussian basis function. If leave it as NULL, knots will be automatically chosen based on the sample size.
#' @param cens_status a formula for censored status in the data with sequence: event, left, right, interval censoring. ~ delta + deltaL + deltaR + deltaI.
#' @return fit parameter estimates, standard deviations and
#' @export
#' @importFrom stats dnorm lm pnorm quantile qweibull rbinom runif rweibull model.frame na.omit
#' @importFrom dplyr filter mutate
#'



aftsur <- function(formula_aft, formula_logit, cens_status, y_L, y_R, data, lambda = 1e-5, knots = NULL){
  deltaI <- NULL
  num_knots <- ifelse(!is.null(knots), knots, ifelse(dim(data)[1]<500, 4, ifelse(dim(data)[1]<1000, 5, 6)))
  X <- model.frame(formula_aft, data)
  colnames(X) <- paste0("X", 1:dim(X)[2])
  Z <- model.frame(formula_logit, data)
  colnames(Z) <- paste0("Z", 1:dim(Z)[2])
  cens_status <- model.frame(cens_status, data)
  colnames(cens_status) <- c("delta", "deltaL",  "deltaR", "deltaI")


  data <- tibble(X, Z, cens_status, y_L, y_R)
  data_tmp <- data %>% mutate(y = ifelse(deltaL!=1 , y_L, y_R))
  data_tmp_append <- data_tmp %>% filter(deltaI==1) %>% mutate(y = y_R, deltaI = 2)
  data <- rbind(data_tmp, data_tmp_append)

  X <- as.matrix(data[, colnames(X)])
  Z <- as.matrix(data[, colnames(Z)])
  y <- data$y
  delta <- data$delta
  deltaL <- data$deltaL
  deltaR <- data$deltaR
  deltaI_L <- as.numeric(data$deltaI == 1)
  deltaI_R <- as.numeric(data$deltaI == 2)
  MAX_CTR <- 10000

  gamma_update <- FALSE
  update_basis_params <- TRUE
  max_lambda_update <- 5
  df_prev <- -10
  beta <- rep(0, dim(X)[2])
  gamma <- rep(0, dim(Z)[2])
  theta <- rep(1,  num_knots)
  val <- initialise_values_list(X, Z, y,
                                delta, deltaL, deltaR, deltaI_L, deltaI_R,
                                beta, gamma, theta, lambda, num_knots)
  ctr <- j <- 1
  while(j < max_lambda_update)
  {

    while(ctr <= MAX_CTR)
    {

      val_previous <- val

      # beta update
      step <- 1.
      penalised_likelihood_ <- penalised_likelihood(val)
      beta_gradient_ <- beta_gradient(val)
      beta_hessian_ <- beta_hessian(val, TRUE)
      while(step >= 1e-2)
      {
        beta_new <- val$beta - step * c(solve(beta_hessian_) %*% beta_gradient_)
        val_proposal <- update_values_list(val, beta_new=beta_new, update_basis_params=update_basis_params)
        if(penalised_likelihood(val_proposal) > penalised_likelihood_)
        {
          val <- val_proposal
          break
        }
        step <- 0.5 * step
      }

      # gamma update
      if(sum(abs(val$beta - val_previous$beta))<1e-3){
        gamma_update <- TRUE
      }
      if(gamma_update == TRUE){
        step <- 1.
        penalised_likelihood_ <- penalised_likelihood(val)
        gamma_gradient_ <- gamma_gradient(val)
        gamma_hessian_ <- gamma_hessian(val, TRUE)
        while(step >= 1e-2)
        {
          gamma_new <- val$gamma - step * c(solve(gamma_hessian_) %*% gamma_gradient_)
          val_proposal <- update_values_list(val, gamma_new=gamma_new)
          if(penalised_likelihood(val_proposal) > penalised_likelihood_)
          {
            val <- val_proposal
            break
          }
          step <- 0.5 * step
        }
      }
      # theta update
      step <- 1
      penalised_likelihood_ <- penalised_likelihood(val)
      theta_gradient_ <- theta_gradient(val)
      theta_hessian_ <- theta_hessian(val, TRUE)

      while(step >= 1e-2)
      {
        theta_new <- val$theta - step * c(theta_hessian_ %*% theta_gradient_)
        val_proposal <- update_values_list(val, theta_new=theta_new)
        if(penalised_likelihood(val_proposal) > penalised_likelihood_)
        {
          val <- val_proposal
          break
        }
        step <- 0.5 * step
      }

      beta_delta <- sum(abs(val$beta - val_previous$beta))
      gamma_delta <- sum(abs(val$gamma - val_previous$gamma))
      theta_delta <- sum(abs(val$theta - val_previous$theta))

      if (update_basis_params & (beta_delta < 1e-4))
      {
        update_basis_params <- FALSE
        fix_basis_ctr <- ctr
      }

      if (max(beta_delta, gamma_delta) < 1e-7 & theta_delta< 1e-6)
      {
        break
      }

      ctr <- ctr + 1

    } # inner loop ends
    if(ctr > MAX_CTR)
    {
      print("Max Iteration Reached")
      break
    }
    # break
    df_lambda_list <- compute_lambda(val)
    if(df_lambda_list$lambda >= 0){
      if((j == 1 | (abs(df_lambda_list$df - df_prev) > 1)) & j!= max_lambda_update - 1){
        df_prev <- df_lambda_list$df
        val <- update_values_list(val,lambda_new=df_lambda_list$lambda)
      } else {
        break
      }
    }
    j <- j + 1
  } #end out outer loop



  class(val) <- c("aftsur")
  aft_beta <- t(t(val$beta))
  logit_gamma <- t(t(val$gamma))

  sd_ <- sqrt(diag(compute_covariance_matrix(val,return_F = T)))
  betasd <-  t(t(sd_[1:length(val$beta)]))
  gammasd <-  t(t(sd_[(length(val$beta)+1) :(length(val$beta)+length(val$gamma))]))
  betaz <- val$beta/betasd
  gammaz <- val$gamma/gammasd
  beta.pval <- pnorm(abs(betaz),lower.tail = F)*2
  gamma.pval <- pnorm(abs(gammaz),lower.tail = F)*2

  beta_out <- cbind(aft_beta, betasd, betaz, beta.pval)
  gamma_out <- cbind(logit_gamma, gammasd, gammaz, gamma.pval)
  colnames(beta_out) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
  colnames(gamma_out) <- c("Estimate","Std.Error","Z value","Pr(>|Z|)")
  val$pzllh <- penalised_likelihood(val)

  cat("---------------------------------------------------------------------------")
  cat("\nSemi-parametric Accelerated Failure Time Mixture Cured Model Using MPL\n")
  cat(paste("\nPenalised log-likelihood:", round(val$pzllh,3),"\n"))
  cat(paste("\nEstimated smoothing parameter:", round(val$lambda,3),"\n"))
  cat("=====")
  cat("\n")
  cat("\nAccelerated Failure Time Model:\n")
  print(beta_out)
  cat("\n")
  cat("=====")
  cat("\nLogistic Model:\n")
  print(gamma_out)
  cat("---------------------------------------------------------------------------")
  predicth <- list()
  predicth[["x"]] <- seq(0, max(val$k[val$delta==1]), 0.01)
  predicth[["h"]] <- h0(predicth[["x"]], val$theta, val$basis_params)
  val$predicth <- predicth
  val
  invisible(val)

}
