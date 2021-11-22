#' Semi-parametric Accelerated Failure Time Mixture Cure Model
#' @description A method to fit a semi-parametric AFT mixture cure model using maximum penalised likelihood estimation, where left, right and interval-censored data are allowed. Logistic regression is used to model the incidence.
#'
#' @param formula a formula object which specifies the covariates for the latency part of the mixture cure model. The response must be a \code{Surv} object which is from \code{survival} package.
#' @param cure_var a formula object which specifies the covariates for the incidence part of the mixture cure model.
#' @param offset if offset is \code{FALSE}, an intercept term will be added into the incidence covariates. By default, \code{offset = FALSE}.
#' @param lambda an initial value for the smoothing parameter. By default, lambda = $10^{-5}$.
#' @param knots an integer which specifies the number of basis functions used to estimate the baseline hazard function. By default, \code{knots = 4} if the sample size is less than 500; \code{knots = 5} if the sample size is greater than 500 but less than 1000 and \code{knots = 6} if the sample size is greater or equal to 1000.
#' @param data a data frame which includes survival times, covariates, censoring status.
#' @return \code{aftsur} returns an object of class \code{"aftsur"}.
#' @examples
#' require(survival)
#' # load data
#' data("ptces")
#' # create Surv object
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 + X3 - 1
#' # fit model
#' aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = ptces)
#' @export
#' @importFrom stats dnorm lm pnorm quantile qweibull rbinom runif rweibull model.frame na.omit model.extract model.matrix
#' @importFrom dplyr filter mutate
#' @importFrom survival Surv
#'



aftsur <- function(formula, cure_var, offset = FALSE, lambda = 1e-5, knots = NULL, data){
  deltaI <- y_tmp <- NULL
  n <- dim(data)[1]
  num_knots <- ifelse(!is.null(knots), knots, ifelse(n < 500, 4, ifelse(n < 1000, 5, 6)))

  m_f <- model.frame(formula, data)
  m_rsp <- model.extract(m_f,"response")
  time_tab <- m_rsp[,1:2]
  colnames(time_tab) <- c("y", "y_tmp")
  #Z X
  cure_name <- all.vars(cure_var)
  Z <- as.matrix(cbind(rep(1,n),data[,cure_name]))
  colnames(Z) <- c("intercept.z", cure_name)
  if(offset == TRUE){Z <- as.matrix(data[,cure_name])}
  X <- model.matrix(attr(m_f,"terms"), m_f)
  # cured
  fac_status <- factor(m_rsp[,3], levels = c(1, 2, 0, 3))
  ces_status <- model.matrix(~ fac_status - 1)
  colnames(ces_status) <- c("delta", "deltaL",  "deltaR", "deltaI")

  data_sur <-  tibble(data.frame(cbind(X, Z, time_tab, ces_status)))
  data_append <- tibble(data_sur %>% filter(y_tmp!=1) %>% mutate(y = y_tmp, deltaI = 2))
  data_ <- rbind(data_sur, data_append) %>% select(-y_tmp)

  X <- as.matrix(data_[, colnames(data.frame(X))])
  if("X.Intercept." %in% colnames(X)){colnames(X)[1] <- c("(intercept)")}
  Z <- as.matrix(data_[, colnames(Z)])
  if("intercept.z" %in% colnames(Z)){colnames(Z)[1] <- c("(intercept)")}
  y <- data_$y
  delta <- data_$delta
  deltaL <- data_$deltaL
  deltaR <- data_$deltaR
  deltaI_L <- as.numeric(data_$deltaI == 1)
  deltaI_R <- as.numeric(data_$deltaI == 2)
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
