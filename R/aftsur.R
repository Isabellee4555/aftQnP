#' Semi-parametric Accelerated Failure Time Mixture Cure Model
#' @description A method to fit  a semi-parametric accelerated failure time (AFT) mixture cure model using maximum penalised likelihood estimation. This method allows for the inclusion of left, right, and interval-censored data. Logistic regression is employed to model the incidence.
#'
#' @param formula a formula object to specify the covariates for the latency part of the mixture cure model. The response variable must be a \code{Surv} object obtained from the \code{survival} package.
#' @param cure_var a formula object in which to specify the covariates for the incidence part of the mixture cure model.
#' @param offset if offset is \code{FALSE}, an intercept term will be added into the incidence covariates. By default, \code{offset = FALSE}.
#' @param lambda an initial value for the smoothing parameter. By default, \code{lambda = 1e-5}.
#' @param knots an integer which specifies the number of basis functions used to estimate the baseline hazard function. By default, \code{knots = 4} if the sample size is less than 500; \code{knots = 5} if the sample size is greater than 500 but less than 1000 and \code{knots = 6} if the sample size is greater or equal to 1000.
#' @param data a data frame which includes survival times, covariates, censoring statuses. Covairates can be numerical or categorical. If categorical variables are supplied, these variables should be a class of \code{Factor}. Censoring statuses need to be either \code{0} or \code{1}.
#' @param basis_range A range of quantile basis functions should be chosen from the accelerated time.
#' @return \code{aftsur} returns an object of class \code{"aftsur"}.
#' @examples
#' require(survival)
#' # load data
#' data("simu_data")
#' # create a formula using a Surv
#' formula_aft <- Surv(y_L, y_R, type = "interval2") ~ X1 + X2 - 1
#' # fit a model
#' aftsur(formula = formula_aft, cure_var = ~ Z1 + Z2 + Z3, offset = TRUE, data = simu_data)
#' @export
#' @importFrom stats dnorm lm pnorm quantile qweibull rbinom runif rweibull model.frame na.omit model.extract model.matrix
#' @importFrom dplyr filter mutate case_when select
#' @importFrom tibble tibble
#' @importFrom survival Surv
#'



aftsur <- function(formula, cure_var, offset = FALSE, lambda = 1e-5, knots = NULL, data, basis_range = c(0.05, 0.9)){
  deltaI <- y_tmp <- NULL
  n <- dim(data)[1]

  num_knots <- ifelse(!is.null(knots), knots, dplyr::case_when(
    n < 500 ~ 4,
    n >= 500 & n < 1000 ~ 5,
    n >= 1000 ~ 6
  ))

  m_f <- model.frame(formula, data)
  m_rsp <- model.extract(m_f,"response")
  time_tab <- m_rsp[,1:2]
  colnames(time_tab) <- c("y", "y_tmp")

  #Z X
  Z <- model.matrix(cure_var, data = data)

  fac_ind <- data %>%
    select(attr(terms(cure_var), "term.labels")) %>%
    sapply(., is.factor)

  cure_name <- all.vars(cure_var)

  if(offset == TRUE){
    Z <- Z[, -1]
  }else{
    colnames(Z) <- c("intercept.z", colnames(Z)[-1])
  }

  # if(offset == TRUE | any(fac_ind)){
  #   Z <- Z[, -1]
  # }else{
  #   colnames(Z) <- c("intercept.z", colnames(Z)[-1])
  # }

  X <- model.matrix(attr(m_f,"terms"), m_f)

  fac_ind <- data %>%
    select(attr(terms(formula), "term.labels")) %>%
    sapply(., is.factor)

  if(any(fac_ind) & (!"(Intercept)" %in% colnames(X))){
    X <- X[,-1]
  }

  # cured
  fac_status <- factor(m_rsp[,3], levels = c(1, 2, 0, 3))
  ces_status <- model.matrix(~ fac_status - 1)
  colnames(ces_status) <- c("delta", "deltaL",  "deltaR", "deltaI")

  data_sur <-  tibble(data.frame(cbind(X, Z, time_tab, ces_status)))
  data_append <- tibble(data_sur %>% filter(y_tmp!=1) %>% mutate(y = y_tmp, deltaI = 2))
  data_ <- rbind(data_sur, data_append) %>% select(-y_tmp)

  X <- as.matrix(data_[, colnames(data.frame(X))])
  colnames(X)[colnames(X) == "X.Intercept."] <- "(intercept)"
  Z <- as.matrix(data_[, colnames(Z)])
  colnames(Z)[colnames(Z) == "intercept.z"] <- "(intercept)"
  y <- data_$y
  delta <- data_$delta
  deltaL <- data_$deltaL
  deltaR <- data_$deltaR
  deltaI_L <- as.numeric(data_$deltaI == 1)
  deltaI_R <- as.numeric(data_$deltaI == 2)
  MAX_CTR <- 10000

  gamma_update <- FALSE
  if(sum(deltaR) == 0){cured_exist <- FALSE; stop("no right-censoring is not implemented")}else{cured_exist <- TRUE}
  update_basis_params <- TRUE
  max_lambda_update <- 5
  df_prev <- -10
  beta <- rep(0, dim(X)[2])
  gamma <- rep(0, dim(Z)[2])
  theta <- rep(1,  num_knots)
  val <- initialise_values_list(X, Z, y,
                                delta, deltaL, deltaR, deltaI_L, deltaI_R,
                                beta, gamma, theta, lambda, num_knots, range = basis_range)
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
      if(sum(abs(val$beta - val_previous$beta))<1e-3 & cured_exist){
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
      # print("Max Iteration Reached")
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
  val$call <- list(latency = formula, incidence = cure_var)

  predicth <- list()
  predicth[["x"]] <- seq(0, max(val$k[val$deltaR==0]), 0.01)
  predicth[["h"]] <- h0(predicth[["x"]], val$theta, val$basis_params)
  val$predicth <- predicth
  val$num_iterations <- ctr
  val
  invisible(val)

}
