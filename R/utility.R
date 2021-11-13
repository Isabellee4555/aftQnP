
phi <- function(y) dnorm(y)

Phi <- function(y) pnorm(y)

psi <- function(y, basis_params)
{
  if(length(y)!=0){
    tmp <- sapply(1:length(basis_params$mu), function(i) sqrt(2 * pi) * phi((y - basis_params$mu[i]) / basis_params$sigma[i]))
  }else{
    tmp <- matrix(0, nrow = 0, ncol = length(basis_params$mu))
  }
  return(tmp)
}

Psi <- function(y, basis_params)
{
  tmp <- sapply(1:length(basis_params$mu), function(i) sqrt(2 * pi) * basis_params$sigma[i] * (Phi((y - basis_params$mu[i]) / basis_params$sigma[i]) - Phi(- basis_params$mu[i] / basis_params$sigma[i])))
  return(tmp)
}

psi_prime <- function(y, basis_params)
{
  psi_ <- psi(y, basis_params)
  if(length(y)!=0){
    tmp <- psi_*sapply(1:length(basis_params$mu), function(i) (-(y - basis_params$mu[i]) / basis_params$sigma[i]^2))
  }else{
    tmp <- matrix(0, nrow = 0, ncol = length(basis_params$mu))
  }
  return(tmp)
}

psi_prime_2 <- function(y, basis_params)
{
  tmp <- psi(y, basis_params) * sapply(1:length(basis_params$mu), function(i) (((y - basis_params$mu[i]) / basis_params$sigma[i]^2)^2 - 1 / basis_params$sigma[i]^2))
  return(tmp)
}

logit <- function(x) exp(x) / (1 + exp(x))

cure_prob <- function(Z, gamma) c(logit(Z %*% gamma))

h0 <- function(y, theta, basis_params) c(psi(y, basis_params) %*% theta)

H0 <- function(y, theta, basis_params) c(Psi(y, basis_params) %*% theta)

h0_prime <- function(y, theta, basis_params) c(psi_prime(y, basis_params) %*% theta)

h0_prime_2 <- function(y, theta, basis_params) c(psi_prime_2(y, basis_params) %*% theta)


compute_basis_params <- function(val)
{
   range_k <- quantile(unique(val$k), c(0.05, 0.95))
  mu <- unname(quantile(range_k, seq(0, 1, length=val$num_knots)))
  sigma <- rep(diff(mu)[1]*0.7, val$num_knots)

  return(list(mu=mu, sigma=sigma))
}

compute_R_matrix <- function(val, integral_delta=0.01, length_out = 500)#, length_out = 500
{
  upper_bound <- max(val$k)
  x_grid <- seq(min(val$k), upper_bound, length.out = length_out)
  integral_delta <- diff(x_grid[1:2])
  psi_prime_2_mat <- psi_prime_2(x_grid, val$basis_params)
  return(t(psi_prime_2_mat) %*% psi_prime_2_mat * integral_delta)
}

initialise_values_list <- function(X, Z, y, delta, deltaL, deltaR, deltaI_L, deltaI_R, beta, gamma, theta, lambda, num_knots=5, basis_params=NULL)
{
  val <- list(X=X, Z=Z, y=y,
              delta=delta, deltaL=deltaL, deltaR=deltaR, deltaI_L=deltaI_L, deltaI_R=deltaI_R,
              beta=beta, gamma=gamma, theta=theta, lambda=lambda, num_knots=num_knots)

  val$k <- c(y * exp(-X %*% beta))
  if(is.null(basis_params)) basis_params <- compute_basis_params(val)
  val$basis_params <- basis_params
  val$R_matrix <- compute_R_matrix(val)
  val$pi <- cure_prob(Z, gamma)
  val$h <- h0(val$k, theta, basis_params)
  val$H <- H0(val$k, theta, basis_params)
  val$S <- exp(-val$H)
  val$h_prime <- h0_prime(val$k, theta, basis_params)
  return(val)
}

update_values_list <- function(val, beta_new=NULL, gamma_new=NULL, theta_new=NULL, lambda_new=NULL, update_basis_params=TRUE)
{
  if(!is.null(beta_new))
  {
    val$beta <- beta_new
    val$k <- c(val$y * exp(-val$X %*% val$beta))
    if(update_basis_params)
    {
      val$basis_params <- compute_basis_params(val)
      val$R_matrix <- compute_R_matrix(val)
    }
    val$h <- h0(val$k, val$theta, val$basis_params)
    val$H <- H0(val$k, val$theta, val$basis_params)
    val$S <- exp(-val$H)
    val$h_prime <- h0_prime(val$k, val$theta, val$basis_params)
  }

  if(!is.null(gamma_new))
  {
    val$gamma <- gamma_new
    val$pi <- cure_prob(val$Z, val$gamma)
  }

  if(!is.null(theta_new))
  {
    val$theta <- theta_new
    val$h <- h0(val$k, val$theta, val$basis_params)
    val$H <- H0(val$k, val$theta, val$basis_params)
    val$S <- exp(-val$H)
    val$h_prime <- h0_prime(val$k, val$theta, val$basis_params)
  }

  if(!is.null(lambda_new)) val$lambda <- lambda_new

  return(val)
}

penalised_likelihood <- function(val) return(likelihood(val) - penalty(val))

likelihood <- function(val)
{
  # event
  term1 <- -c(as.matrix(val$X[val$delta == 1, ]) %*% val$beta) + log(pmax(val$h[val$delta == 1], 1e-12)) -
    val$H[val$delta == 1] + log(val$pi[val$delta == 1])

  # left
  term2 <- log(1 - val$S[val$deltaL == 1]) + log(val$pi[val$deltaL == 1])

  # right
  pi_right <- val$pi[val$deltaR == 1]
  term3 <- log(val$S[val$deltaR == 1] * pi_right + 1 - pi_right)

  # interval
  term4 <- log(val$S[val$deltaI_L == 1] - val$S[val$deltaI_R == 1]) + log(val$pi[val$deltaI_L == 1])

  likelihood_ <- sum(c(term1, term2, term3, term4))
  return(likelihood_)
}

penalty <- function(val) return(c(val$lambda * t(val$theta) %*% val$R_matrix %*% val$theta))

penalty_gradient <- function(val, use_quasi=FALSE)
{
  penalty_gradient_full <- c(2 * val$lambda * val$R_matrix %*% val$theta)
  if(!use_quasi) return(penalty_gradient_full)
  return(pmax(penalty_gradient_full, 0))
}

penalty_hessian <- function(val, use_quasi=FALSE) return(2 * val$lambda * val$R_matrix)

beta_gradient <- function(val)
{
  # event
  h_event <- val$h[val$delta == 1]
  term1 <- -1 - val$k[val$delta == 1] * (val$h_prime[val$delta == 1] - h_event^2) / pmax(h_event, 1e-12)

  # left
  S_left <- val$S[val$deltaL == 1]
  term2 <- -S_left * val$h[val$deltaL == 1] * val$k[val$deltaL == 1] / (1 - S_left)

  # right
  pi_right <- val$pi[val$deltaR == 1]
  S_right<- val$S[val$deltaR == 1]
  term3_num <- pi_right *  S_right * val$h[val$deltaR == 1] * val$k[val$deltaR == 1]
  term3_den <- pi_right * S_right + 1 - pi_right
  term3 <- term3_num / term3_den

  # interval
  S_L_interval <- val$S[val$deltaI_L == 1]
  S_R_interval <- val$S[val$deltaI_R == 1]
  term4_num <- (S_L_interval * val$h[val$deltaI_L == 1] * val$k[val$deltaI_L == 1] - S_R_interval * val$h[val$deltaI_R == 1] * val$k[val$deltaI_R == 1])
  term4_den <- (S_L_interval - S_R_interval)
  term4 <- term4_num / term4_den

  X_sort <- rbind(as.matrix(val$X[val$delta == 1, ]),
                  as.matrix(val$X[val$deltaL == 1, ]),
                  as.matrix(val$X[val$deltaR == 1, ]),
                  as.matrix(val$X[val$deltaI_L == 1, ]))
  term_all <- c(term1, term2, term3, term4)

  return(c(t(X_sort) %*% term_all))
}

beta_hessian <- function(val, use_quasi=FALSE)
{
  # extract values for event
  h_event <- val$h[val$delta == 1]
  h_prime_event <- val$h_prime[val$delta == 1]
  k_event <- val$k[val$delta == 1]
  h_prime_2_event <- if(length(k_event)==0){k_event}else{
    h0_prime_2(k_event, val$theta, val$basis_params)
  }

  # extract values for left
  k_left <- val$k[val$deltaL == 1]
  S_left <- val$S[val$deltaL == 1]
  h_left <- val$h[val$deltaL == 1]

  # extract values for right
  k_right <- val$k[val$deltaR == 1]
  S_right <- val$S[val$deltaR == 1]
  h_right <- val$h[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]

  # extract values for interval
  k_L_interval <- val$k[val$deltaI_L== 1]
  k_R_interval <- val$k[val$deltaI_R== 1]
  S_L_interval <- val$S[val$deltaI_L== 1]
  S_R_interval <- val$S[val$deltaI_R== 1]
  h_L_interval <- val$h[val$deltaI_L== 1]
  h_R_interval <- val$h[val$deltaI_R== 1]

  # negative components

  # event
  term1_negative <- k_event * (h_event + (k_event * h_prime_event^2) / pmax(h_event^2, 1e-12))

  # left
  term2_negative <- k_left * S_left * (k_left * h_left^2 / (1 - S_left) + S_left * h_left^2 * k_left / (1 - S_left)^2)

  # right
  term3_negative <- k_right * S_right * pi_right * (h_right / (pi_right * S_right + 1 - pi_right) + pi_right * S_right * h_right^2 * k_right / (pi_right * S_right + 1 - pi_right)^2)

  # interval
  u_L <- k_L_interval * S_L_interval * h_L_interval
  u_R <- k_R_interval * S_R_interval * h_R_interval
  v <- S_L_interval - S_R_interval
  u_L_prime <- k_L_interval * S_L_interval * (k_L_interval * h_L_interval^2 - k_L_interval * val$h_prime[val$deltaI_L == 1] - h_L_interval)
  u_R_prime <- k_R_interval * S_R_interval * (k_R_interval * h_R_interval^2 - k_R_interval * val$h_prime[val$deltaI_R == 1] - h_R_interval)
  v_prime <- u_L - u_R
  term4_negative <- (u_L + S_R_interval * h_R_interval^2 * k_R_interval^2) / v + (u_L^2 + u_R^2) / v^2

  # return if use_quasi
  X_sort <- rbind(as.matrix(val$X[val$delta == 1, ]),
                  as.matrix(val$X[val$deltaL == 1, ]),
                  as.matrix(val$X[val$deltaR == 1, ]),
                  as.matrix(val$X[val$deltaI_L == 1, ]))
  term_negative_all <- c(term1_negative, term2_negative, term3_negative, term4_negative)
  if(use_quasi) return(-t(X_sort) %*% diag(term_negative_all, nrow = length(term_negative_all)) %*% X_sort)

  # all components

  # event
  term1_full <- k_event * (h_prime_event / h_event - 2 * k_event * h_prime_event + k_event * h_prime_2_event / h_event + k_event * h_prime_event) - term1_negative

  # left
  term2_full <- k_left * S_left * (k_left * val$h_prime[val$deltaL == 1] + h_left) / (1 - S_left) - term2_negative

  # right
  term3_full <- k_right * S_right * pi_right * (k_right * h_right^2 - k_right * val$h_prime[val$deltaR == 1]) / (pi_right * S_right + 1 - pi_right) - term3_negative

  # interval
  term4_full <- (u_L_prime - u_R_prime) / v - (u_L - u_R) * v_prime / v^2

  # return if use_quasi == FALSE
  term_full_all <- c(term1_full, term2_full, term3_full, term4_full)
  return(t(X_sort) %*% diag(term_full_all, nrow = length(term_full_all)) %*% X_sort)
}

gamma_gradient <- function(val)
{
  # event, left, interval
  term1 <- c(t(val$Z[(val$deltaR != 1)  & (val$deltaI_R != 1), ]) %*% (1 - val$pi[(val$deltaR != 1)  & (val$deltaI_R != 1)]))

  # right
  pi_right <- val$pi[val$deltaR == 1]
  S_right <- val$S[val$deltaR == 1]
  term2 <- -c(t(val$Z[val$deltaR == 1, ]) %*% (pi_right * (1 - pi_right) * (1 - S_right) / (S_right * pi_right + 1 - pi_right)))

  return(term1 + term2)
}

gamma_hessian <- function(val, use_quasi=FALSE)
{
  # negative components
  Z <- val$Z[val$deltaI_R != 1, ]
  pi_ <-  val$pi[val$deltaI_R != 1]
  term_negative_const <- pi_ * (1 - pi_)
  term_negative <- -t(Z) %*% diag(term_negative_const, nrow=length(term_negative_const)) %*% Z
  if(use_quasi) return(term_negative)

  # other components
  Z_right <- val$Z[val$deltaR == 1, ]
  pi_right <-  val$pi[val$deltaR == 1]
  S_right <- val$S[val$deltaR == 1]
  term_other_constant <- pi_right * (1 - pi_right) * S_right / (pi_right * S_right + 1 - pi_right)^2
  term_other <- t(Z_right) %*% diag(term_other_constant, nrow = length(term_other_constant)) %*% Z_right
  return(term_other + term_negative)
}

theta_gradient <- function(val, use_quasi=FALSE)
{
  Psi_mat <- Psi(val$k, val$basis_params)

  # extract values for event
  k_event <- val$k[val$delta == 1]
  Psi_mat_event <- matrix(Psi_mat[val$delta == 1, ], ncol = val$num_knots)
  psi_mat_event <- matrix(psi(k_event, val$basis_params), ncol = val$num_knots)
  # extract values for left
  S_left <- val$S[val$deltaL == 1]
  Psi_mat_left <- matrix(Psi_mat[val$deltaL == 1, ], ncol = val$num_knots)

  # extract values for right
  S_right <- val$S[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]
  Psi_mat_right <- matrix( Psi_mat[val$deltaR == 1, ], ncol = val$num_knots)

  # extract values for interval
  S_L_interval <- val$S[val$deltaI_L == 1]
  S_R_interval <- val$S[val$deltaI_R == 1]

  Psi_mat_L_interval <- matrix( Psi_mat[val$deltaI_L == 1, ], ncol = val$num_knots)
  Psi_mat_R_interval <- matrix(Psi_mat[val$deltaI_R == 1, ], ncol = val$num_knots)

  # negative components
  # event
  term1_negative <- colSums(Psi_mat_event)

  # right
  term3_negative <- colSums(Psi_mat_right * (pi_right * S_right / (pi_right * S_right + 1 - pi_right)))

  # interval
  den <- S_L_interval - S_R_interval
  term4_negative <- colSums(Psi_mat_L_interval * (S_L_interval / den))

  # return here if use_quasi
  if (use_quasi == TRUE) {
    theta_gradient_negative <- -(term1_negative + term3_negative + term4_negative + penalty_gradient(val, TRUE))
    return(theta_gradient_negative)
    # return(pmin(theta_gradient_negative, -1))
  }

  # all components
  # event
  term1_full <- colSums(psi_mat_event / val$h[val$delta == 1]) - term1_negative

  # left
  term2_full <- colSums(Psi_mat_left * (S_left / (1 - S_left)))

  # right
  term3_full <- -term3_negative

  # interval
  term4_full <- colSums(Psi_mat_R_interval * (S_R_interval / den)) - term4_negative

  # return here if use_quasi == FALSE
  return(term1_full + term2_full + term3_full + term4_full - penalty_gradient(val))
}

theta_hessian <- function(val, use_quasi=FALSE)
{
  # return here if use_quasi
  if (use_quasi == TRUE)
  {
    theta_gradient_negative <- theta_gradient(val, TRUE)
    tmp <- val$theta / (theta_gradient_negative - 1e-6)
    tmp <- ifelse(theta_gradient(val)>0, 1/(theta_gradient_negative - 1e-6), tmp)
    return(diag(tmp, nrow=length(tmp)))
  }

  Psi_mat <- Psi(val$k, val$basis_params)

  # event
  h_event <- val$h[val$delta == 1]
  term1 <- -1 / h_event^2

  # left
  k_left <- val$k[val$deltaL == 1]
  S_left <- val$S[val$deltaL == 1]
  term2 <- -S_left / (1 - S_left)^2

  # right
  k_right <- val$k[val$deltaR == 1]
  S_right <- val$S[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]
  term3 <- pi_right * (1 - pi_right) * S_right / (pi_right * S_right + 1 - pi_right)^2

  # interval
  k_L_interval <- val$k[val$deltaI_L == 1]
  k_R_interval <- val$k[val$deltaI_R == 1]
  S_L_interval <- val$S[val$deltaI_L == 1]
  S_R_interval <- val$S[val$deltaI_R == 1]

  Psi_mat_L_interval <- matrix( Psi_mat[val$deltaI_L == 1, ], ncol = val$num_knots)
  Psi_mat_R_interval <- matrix(Psi_mat[val$deltaI_R == 1, ], ncol = val$num_knots)

  u_R <- S_R_interval * Psi_mat_R_interval
  u_L <- S_L_interval * Psi_mat_L_interval
  v <- S_L_interval - S_R_interval
  v_prime <- u_L - u_R

  u_R_prime_over_v <- t(Psi_mat_R_interval) %*% diag(S_R_interval / v, nrow=length(v)) %*% Psi_mat_R_interval
  u_L_prime_over_v <- t(Psi_mat_L_interval) %*% diag(S_L_interval / v, nrow=length(v)) %*% Psi_mat_L_interval

  term4_1 <- u_L_prime_over_v - u_R_prime_over_v
  term4_2 <- t(u_L - u_R) %*% diag(1 / v^2, nrow=length(v)) %*% v_prime
  term_interval <- term4_1 - term4_2

  Psi_mat_sort_non_interval <- rbind(psi(val$k[val$delta == 1], val$basis_params), matrix( Psi_mat[val$deltaL == 1, ], ncol = val$num_knots), matrix( Psi_mat[val$deltaR == 1, ], ncol = val$num_knots))
  term_full_non_interval <- c(term1, term2, term3)

  return(t(Psi_mat_sort_non_interval) %*% diag(term_full_non_interval, nrow=length(term_full_non_interval)) %*% Psi_mat_sort_non_interval + term_interval)
}

beta_theta_hessian <- function(val)
{
  Psi_mat <- Psi(val$k, val$basis_params)
  psi_mat <- psi(val$k, val$basis_params)

  # event
  k_event <- val$k[val$delta == 1]
  h_event <- val$h[val$delta == 1]
  psi_mat_event <- matrix(psi_mat[val$delta == 1, ], ncol = val$num_knots)
  term1_c <- k_event
  term1_left <- psi_mat_event - psi_prime(k_event, val$basis_params) / h_event + psi_mat_event * (val$h_prime[val$delta == 1] / h_event^2)

  # left
  k_left <- val$k[val$deltaL == 1]
  h_left <- val$h[val$deltaL == 1]
  S_left <- val$S[val$deltaL == 1]
  psi_mat_left <- matrix(psi_mat[val$deltaL == 1, ], ncol = val$num_knots)
  term2_c <- k_left * S_left / (1 - S_left)^2
  term2_left <- matrix( Psi_mat[val$deltaL == 1, ], ncol = val$num_knots) * h_left + psi_mat_left * S_left - psi_mat_left

  # right
  k_right <- val$k[val$deltaR == 1]
  h_right <- val$h[val$deltaR == 1]
  S_right <- val$S[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]

  term3_c <- S_right * k_right * pi_right
  term3_left_1 <- matrix(psi_mat[val$deltaR == 1, ], ncol = val$num_knots) / (pi_right * S_right + 1 - pi_right)
  term3_left_2 <- -matrix( Psi_mat[val$deltaR == 1, ], ncol = val$num_knots) * (h_right * (1 - pi_right) / (pi_right * S_right + 1 - pi_right)^2)
  term3_left <- term3_left_1 + term3_left_2


  # interval
  k_L_interval <- val$k[val$deltaI_L == 1]
  k_R_interval <- val$k[val$deltaI_R == 1]
  S_L_interval <- val$S[val$deltaI_L == 1]
  S_R_interval <- val$S[val$deltaI_R == 1]
  h_L_interval <- val$h[val$deltaI_L == 1]
  h_R_interval <- val$h[val$deltaI_R == 1]
  Psi_mat_L_interval <- matrix( Psi_mat[val$deltaI_L == 1, ], ncol = val$num_knots)
  Psi_mat_R_interval <- matrix(Psi_mat[val$deltaI_R == 1, ], ncol = val$num_knots)

  u <- k_L_interval * S_L_interval * h_L_interval - k_R_interval * S_R_interval * h_R_interval
  u_prime <-  (k_L_interval * S_L_interval * (matrix( psi_mat[val$deltaI_L == 1, ], ncol = val$num_knots) - h_L_interval * Psi_mat_L_interval)
               - k_R_interval * S_R_interval * (matrix(psi_mat[val$deltaI_R == 1, ], ncol = val$num_knots) - h_R_interval * Psi_mat_R_interval))
  v <- S_L_interval - S_R_interval
  v_prime <- S_R_interval * Psi_mat_R_interval - S_L_interval * Psi_mat_L_interval

  term4_c <- 1 / v^2
  term4_left <- v * u_prime - u * v_prime

  X_sort <- rbind(as.matrix(val$X[val$delta == 1, ]),
                  as.matrix(val$X[val$deltaL == 1, ]),
                  as.matrix(val$X[val$deltaR == 1, ]),
                  as.matrix(val$X[val$deltaI_L == 1, ]))
  term_full_c <- c(term1_c, term2_c, term3_c, term4_c)
  term_full_left <- rbind(term1_left, term2_left, term3_left, term4_left)

  return(t(X_sort) %*% diag(term_full_c, nrow = length(term_full_c)) %*% term_full_left)
}

beta_gamma_hessian <- function(val)
{
  S_right <- val$S[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]

  term_num <- val$k[val$deltaR == 1] * pi_right * (1 - pi_right) * S_right * val$h[val$deltaR == 1]
  term_den <- (pi_right * S_right + 1 - pi_right)^2
  term <- term_num / term_den

  return(t(val$X[val$deltaR == 1, ]) %*% diag(term, nrow = length(term)) %*% val$Z[val$deltaR == 1, ])
}

gamma_theta_hessian <- function(val)
{
  S_right <- val$S[val$deltaR == 1]
  pi_right <- val$pi[val$deltaR == 1]

  term_num <- -pi_right * (1 - pi_right) * S_right
  term_den <- (pi_right * S_right + 1 - pi_right)^2
  term <- term_num / term_den
  return(t(val$Z[val$deltaR == 1, ]) %*% diag(term, nrow = length(term)) %*% Psi(val$k[val$deltaR == 1], val$basis_params))
}

full_hessian <- function(val)
{
  beta_hessian_ <- beta_hessian(val)
  gamma_hessian_ <- gamma_hessian(val)
  theta_hessian_ <- theta_hessian(val)

  beta_gamma_hessian_ <- beta_gamma_hessian(val)
  beta_theta_hessian_ <- beta_theta_hessian(val)
  gamma_theta_hessian_ <- gamma_theta_hessian(val)

  row_1 <- cbind(beta_hessian_, beta_gamma_hessian_, beta_theta_hessian_)
  row_2 <- cbind(t(beta_gamma_hessian_), gamma_hessian_, gamma_theta_hessian_)
  row_3 <- cbind(t(beta_theta_hessian_), t(gamma_theta_hessian_), theta_hessian_)

  return(rbind(row_1,row_2, row_3))
}

compute_Q_matrix <- function(val)
{

  length_beta_gamma <- length(val$beta) + length(val$gamma)
  length_full <- length_beta_gamma + length(val$theta)
  Q <- matrix(0, length_full, length_full)
  Q[(length_beta_gamma + 1): length_full, (length_beta_gamma + 1): length_full] <- penalty_hessian(val)
  return(Q)
}

compute_constrained_theta <- function(val)
{
  theta_gradient_ <- theta_gradient(val)
  constrained_theta <- which(val$theta < 0.1 & theta_gradient_< 0)
  return(constrained_theta)
}

compute_lambda <- function(val)
{
  G_matrix <- full_hessian(val)
  Q_matrix <- compute_Q_matrix(val)
  F_matrix <- -G_matrix + Q_matrix

  length_beta_gamma <- length(val$beta) + length(val$gamma)
  index_to_remove <- length_beta_gamma + compute_constrained_theta(val)

  if(length(index_to_remove) == 0)
  {
    F_matrix_removed <- F_matrix
    Q_matrix_removed <- Q_matrix
  } else
  {
    F_matrix_removed <- F_matrix[-index_to_remove, -index_to_remove]
    Q_matrix_removed <- Q_matrix[-index_to_remove, -index_to_remove]
  }

  df <- length(val$theta) - sum(diag(solve(F_matrix_removed + (diag(dim(F_matrix_removed)[1])* 1e-8)) %*% Q_matrix_removed))
  sigma2 <- c(t(val$theta) %*% val$R_matrix %*% val$theta / df)
  lambda <- 0.5 / sigma2
  return(list(df=df, lambda=lambda))
}

compute_covariance_matrix <- function(val, return_F = FALSE)
{
  constraints <- compute_constrained_theta(val)
  G_matrix <- full_hessian(val)
  Q_matrix <- compute_Q_matrix(val)
  F_matrix <- -G_matrix + Q_matrix

  length_beta_gamma <- length(val$beta) + length(val$gamma)
  index_to_remove <- length_beta_gamma + compute_constrained_theta(val)

  if(length(index_to_remove) == 0)
  {
    F_matrix_removed <- F_matrix
    G_matrix_removed <- G_matrix
  } else
  {
    F_matrix_removed <- F_matrix[-index_to_remove, -index_to_remove]
    G_matrix_removed <- G_matrix[-index_to_remove, -index_to_remove]
  }

  F_matrix_removed_inv <- solve(F_matrix_removed + (diag(dim(F_matrix_removed)[1])* 1e-10))
  if(return_F == TRUE)return(F_matrix_removed_inv)
  covariance_matrix <- -F_matrix_removed_inv %*% G_matrix_removed %*% F_matrix_removed_inv
  return(covariance_matrix)
}

compute_hazard_covariance <- function(val, cov_mat=NULL, k=NULL)
{
  if(is.null(cov_mat)) cov_mat <- compute_covariance_matrix(val)
  if(is.null(k)) k <- qweibull(c(0.25, 0.5, 0.75), 3, 1)

  constraints <- compute_constrained_theta(val)
  length_beta_gamma <- length(val$beta) + length(val$gamma)
  length_full <- nrow(cov_mat)

  psi_mat <- psi(k, val$basis_params)
  if(length(constraints) > 0)
  {
    psi_mat_removed <- psi_mat[, -constraints]
  } else {
    psi_mat_removed <- psi_mat
  }
  theta_cov_mat <- cov_mat[(length_beta_gamma + 1):length_full,  (length_beta_gamma + 1):length_full]
  hazard_cov_mat <- psi_mat_removed %*% as.matrix(theta_cov_mat) %*% t(psi_mat_removed)
  return(hazard_cov_mat)
}

compute_cumu_hazard_covariance <- function(val, cov_mat=NULL, k=NULL)
{
  if(is.null(cov_mat)) cov_mat <- compute_covariance_matrix(val)
  if(is.null(k)) k <- qweibull(c(0.25, 0.5, 0.75), 3, 1)

  constraints <- compute_constrained_theta(val)
  length_beta_gamma <- length(val$beta) + length(val$gamma)
  length_full <- nrow(cov_mat)

  Psi_mat <- Psi(k, val$basis_params)
  if(length(constraints) > 0)
  {
    Psi_mat_removed <- Psi_mat[, -constraints]
  } else {
    Psi_mat_removed <- Psi_mat
  }
  theta_cov_mat <- cov_mat[(length_beta_gamma + 1):length_full,  (length_beta_gamma + 1):length_full]
  cumu_hazard_cov_mat <- Psi_mat_removed %*% as.matrix(theta_cov_mat) %*% t(Psi_mat_removed)
  return(cumu_hazard_cov_mat)
}
